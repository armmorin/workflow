from ase.db import connect
from perqueue.queue import PersistentQueue
from sys import argv
from herculestools.dft import RunConfiguration

home = RunConfiguration.home

keys = argv[1].split(',')

db_path = home / 'structures' / 'hexag_perovs_wdiscards.db'
neb_path = home / 'NEB'
preneb_path = home / 'preNEB'

# Given a key, or keys, return the corresponding directory of that entry.
for key in keys:
    with PersistentQueue() as pq:
        entry = pq.get_entry(int(key))
        if entry._task.args['db_id'] is not None:
            db_id = entry._task.args['db_id']
        else:
            db_id = entry.data['neb_id']
        
        with connect(db_path) as db:
            row = db.get(id=db_id)
            split_name = row.name.split('_')
            joint_name = '_'.join(split_name[:-1])
            if 'neb' in split_name:
                job_dir = home / neb_path / joint_name
            else:
                job_dir = home / preneb_path / joint_name
            
            print(job_dir.resolve())