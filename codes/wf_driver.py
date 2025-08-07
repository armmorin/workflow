from sys import argv
from time import sleep
from random import uniform
from perqueue import PersistentQueue
from perqueue.task_classes.task import Task
from perqueue.task_classes.task_groups import Workflow, StaticWidthGroup
from pathlib import Path

j = int(argv[1])
k = int(argv[2])

project_dir = Path.cwd()
codes_dir = f"{project_dir}/codes"

WIDTH = 3
# Defining the tasks for a specific range of entries.
for i in range(j,k):
	t1 = Task(code=f"{codes_dir}/generation.py", args={'index': i}, resources="1:local:10m")
	
	for j in range(WIDTH):
		t2 = Task(code=f"{codes_dir}/relax.py", 
            	args={'index': j}, 
             	resources="40:xeon40:1:50h")
  
		t3 = Task(code=f"{codes_dir}/convex_hull.py", 
            	args={}, 
             	resources="1:local:1:5m")
  
		t4 = Task(code=f"{codes_dir}/band_gap.py", 
            	args={}, 
             	resources="1:local:1:5m")
  
		t5 = Task(code=f"{codes_dir}/preneb.py", 
            	args={'initial_vac': 30, 'final_vac': 31}, 
             	resources="112:xeon56:1:50h")
  
		t6 = Task(code=f"{codes_dir}/neb.py", 
				args={'N_images': 3, 'climb': True, 'parallel' : False}, 
    			resources="168:xeon56:1:50h")

		swf = Workflow({t2: [], t3: [t2], t4: [t3], t5: [t4], t6: [t5]}) 
		
		swg = StaticWidthGroup(swf, width=WIDTH)

		wf = Workflow({t1: [], swg: [t1]})
		
		with PersistentQueue() as pq:
			pq.submit(wf)
			sleep(uniform(.2, 2))
