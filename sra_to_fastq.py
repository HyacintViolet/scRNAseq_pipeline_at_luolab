import os
# import re
import multiprocessing as mp
import subprocess


def work(cmd):
    return subprocess.call(cmd, shell=True)


def main():
    wd = '/media/luolab/ZA1BT1ER/SRA/sra/'
    os.chdir(wd)

    cmd_fastqdump = []
    for item in os.listdir(wd):
        if item.endswith('.sra'):
            this_cmd = 'fastq-dump --gzip --split-3 -A ' + item
            cmd_fastqdump.append(this_cmd)
        elif item.endswith('.fastq.gz'):
            continue

    # Parallel run by Pool
    pool = mp.Pool(16)
    pool.map(work, cmd_fastqdump)


if __name__ == '__main__':
    main()
