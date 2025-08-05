from perqueue.queue import PersistentQueue
from perqueue.selection import Selection
from pathlib import Path
from ase.db import connect
from typing import Optional, List
from perqueue.task_classes.util_classes import Resources
from perqueue.task_classes.entry import Entry
from ase.db.sqlite import SQLite3Database
from sys import argv

# Get the timed out and calculation failed jobs encapsulated in a function
def handle_errors(keys: Optional[List[int]], db_path = 'structures/hexag_perovs_wdiscards.db') -> None:

    calc_failed = "CalculationFailed"
    mem_string = "malloc(): unsorted double linked list corrupted"
    time_out = "DUE TO TIME LIMIT"
    error_string = "I REFUSE TO CONTINUE WITH THIS SICK JOB"
    matching_str = "EEEEEEE  R     R" # String from the bottom row of the ERROR message.
    common_errors = {
        "eddav": "Call to ZHEGV failed",
        "zbrent" : "ZBRENT: fatal error in bracketing",
        "bravais": "Inconsistent Bravais lattice" ,
        "fexcf": "ERROR FEXCF: supplied exchange-correlation",
        "ibrion0" : "Fatal error! IBRION=0",
        "lapack" : "LAPACK: Routine ZPOTRF failed",
        "subrot": "ERROR in subspace rotation",
        "wavecar" : "ERROR: while reading WAVECAR", 
    }
        
    wavecar_restarts = [
    "ERROR: while reading eigenvalues from WAVECAR",
    "ERROR: while reading plane wave coef. from WAVECAR"
    ]

    # Common errors and their solutions by modifying the INCAR file
    common_solutions = {
        "eddav":    {'algo': 'VeryFast'},
        "zbrent":   {'potim': 0.25,
                    'addgrid' : True},
        "bravais":  {'symprec': 1e-09,
                    'isym': 0,
                    'algo': 'VeryFast'},
        "fexcf":    {'potim': 0.25},
        "ibrion0":  {'ibrion' : 2},
        "wavecar":  {'algo': 'VeryFast',
                    'istart': 1},
        "lapack":   {'algo': 'VeryFast',
                    'potim': 0.25},
        "subrot":   {'algo': 'Normal',
                     'potim': 0.05}, # Testing phase to see if it improves the convergence
    }
    db: SQLite3Database = connect(db_path)
    
    if keys:
        s = Selection(ids=keys)
    else:
        s = Selection(states='tf')
    with PersistentQueue() as pq:
        entries = pq.get_entries(selection=s)

    #targets = s.filter(entries)
    codes = [pq.get_code(en.key) for en in entries]
    codes = list(set(codes))

    # PQ assigns the newest MQ id to its entry while still keeping the last state before the run was completed.
    errors_dict = dict()
    restart_jobs = []
    for en in entries:
        mq_id = en.mq_id
        pq_key = en.key
        status = en.state.serialize()
        data = en.data
        data_ids = ['initial_id', 'final_id', 'neb_id']
        db_id = next((data[data_id] for data_id in data_ids if data and data_id in data), None)
        basename = db.get(db_id).name if db_id else db.get(pq.get_args(pq_key).get('db_id')).name
        #if data is not None:
        #    db_id = [data[data_id] for data_id in data_ids if data_id in data] 
        #    basename = db.get(db_id[0]).name 
        #else:
        #    basename = db.get(pq.get_args(pq_key).get('db_id')).name

        name_parts = basename.split('_') if basename is not None else None
        job_name = '_'.join(name_parts[:-1])
        print(f"Checking pq key: {pq_key}, {job_name} with status {status}")
        
        # Use the mq_id to get the directory of the calculation.
        error_file = Path(f"perqueue.runner.{mq_id}.err")
        try:
            if error_file.is_file() and error_file.stat().st_size > 0:
                readlines = error_file.read_text().split('\n')
                print(f"Checking {error_file} for entries {en.key} with name {job_name} and status {status}")
                if [line for line in readlines if time_out in line].count(time_out) > 0:
                    restart_jobs.append(pq_key)
                
                if mem_string in readlines:
                    # Memory corruption error: Attempt at changing the partition
                    resources = str(en._task.resources)
                    if "epyc96" in resources:
                        resources = Resources.from_string("112:1:xeon56:50h") # This can be further automated
                        en._task.resources = resources
                        print(f"Memory corruption error found in {pq_key} changing resources to {resources}")
                        with PersistentQueue() as pq:
                            pq.save_resources(pq_key, resources, False)
                            restart_jobs.append(pq_key)
                
                elif (calc_failed_dir := [line.split(' ') for line in readlines if calc_failed in line][-1][-5]):
                    calc_failed_dir = Path(calc_failed_dir)
                    vaspout = calc_failed_dir / "vasp.out"
                    if vaspout.is_file() and vaspout.stat().st_size > 0:
                        vaspout_lines = vaspout.read_text().split('\n')
                        next_lines = [line.strip().strip('|').strip(' ') for idx, line in enumerate(vaspout_lines) if matching_str in line for line in vaspout_lines[idx+2:-5]]
                        error_msg = " ".join(next_lines)
                        # Save the error message and the directory of the calculation
                        errors_dict[pq_key] = error_msg , calc_failed_dir 
                
                else:
                    pass
        except IndexError:
            print(f"Error reading {error_file} for {en.key} with name {job_name} and status {status}")
            pass
                    
    # XXX: 
    # TODO:
    
    # Resubmit the jobs that timed out
    with PersistentQueue() as pq:
        if len(restart_jobs) > 0:
            print(f"Restarting {len(restart_jobs)} jobs")
        for job in restart_jobs:
            en = pq.get_entry(job)
            pq.resubmit(en)

    # Now match the error messages with the common errors and resubmit the jobs
        if len(errors_dict) > 0:
            print(f"Resubmitting {len(errors_dict)} jobs because of calculation failed")
        for pq_key, _ in errors_dict.items():
            en = pq.get_entry(pq_key)
            error_msg = errors_dict[pq_key][0]
            job_dir = errors_dict[pq_key][1]
            new_args = pq.get_args(pq_key) 
            # Check if the error message contains any of the common errors
            if any([value in error_msg for value in wavecar_restarts]):           
                # Remove the WAVECAR file and resubmit the job
                print(f"Removing the WAVECAR file for {pq_key}")
                wavecar = job_dir / "WAVECAR"
                wavecar.unlink()
                
            for key, value in common_errors.items():
                if value in error_msg:
                    #new_args['vasp'] = common_solutions[key]
                    args = pq.get_args(pq_key)
                    # Update the arguments
                    args['vasp'].update(common_solutions[key])
                    # Save the new arguments
                    pq.save_args(pq_key, args, False)
                    print(f"Error {key} found in {pq_key} the solution is {common_solutions[key]}")
                    print(f"Resubmitting {pq_key} with new args {args}")

            # Resubmit the job
            pq.resubmit(en)
            
if __name__ == "__main__":
    args = [int(split_arg) for arg in argv[1:] for split_arg in arg.split(',')]
    if len(args) >= 1:
        handle_errors(args)
    else:    
        handle_errors(None)