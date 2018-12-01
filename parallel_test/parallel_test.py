import os
import subprocess
import time
import multiprocessing as mp


def firstproc5(i):
    print('Start: first 5  '+str(i))
    time.sleep(5)
    print('End: first 5  '+str(i))


def firstproc10(i):
    print('Start: first 10  '+str(i))
    time.sleep(10)
    print('End: first 10  '+str(i))


def secondproc5(i):
    print('Start: second 5  '+str(i))
    time.sleep(5)
    print('End: second 5  '+str(i))


# This works
pool = mp.Pool(3)
results = pool.map(firstproc5, [1,2,3,4,5,6])

pool2 = mp.Pool(3)
results2 = pool2.map(secondproc5, [1,2,3,4])


# # This doesn't work
# processes = set()
# max_processes = 3
# for i in range(0,5):
#     processes.add(subprocess(firstproc5(i)))
#     if len(processes) >= max_processes:
#         os.wait()
#         processes.difference_update([p for p in processes if p.poll() is not None])

