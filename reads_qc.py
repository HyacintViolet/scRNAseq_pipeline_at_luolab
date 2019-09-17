import os
import re
import multiprocessing as mp
import subprocess


# Function to recursively fetch file path
def fetch_qc_file_list(wd, file_list):
    for item in os.listdir(wd):
        s = os.path.join(wd, item)
        if os.path.isdir(s):
            fetch_qc_file_list(s, file_list)
        elif os.path.basename(s).endswith(('1.fq.gz', '2.fq.gz', '1.fastq.gz', '2.fastq.gz',
                                           '1.clean.fq.gz', '2.clean.fq.gz')):
            match = re.search('^([^_]*)_([^_]*)_(\d)(\..*)$', os.path.basename(s))
            if not any(fname.endswith('_fastqc.html') for fname in os.listdir(os.path.dirname(s))):
                file_list.append(s)
    return file_list


def work(cmd):
    return subprocess.call(cmd, shell=True)


def main():
    wd = '/media/luolab/ZA1BT1ER/scRNAseq/yanting_all/data/yanting/'
    file_list = []

    file_list = fetch_qc_file_list(wd, file_list)

    cmd_fastqc = []
    for file in file_list:
        cmd_fastqc.append('fastqc ' + file)

    # Parallel run by Pool
    pool = mp.Pool(16)
    pool.map(work, cmd_fastqc)

    # The build-in parallel functionality of fastQC runs very slowly, abandoned
    # num_threads = 24
    # cmd = 'fastqc -t ' + str(num_threads) + ' ' + ' '.join(file_list)
    # subprocess.call(cmd, shell=True)

    print('fastqc parallel: finished.')


if __name__ == '__main__':
    main()
