# Used jfs's answer in this thread: https://stackoverflow.com/questions/16450788/python-running-subprocess-in-parallel

import os
# import re
import multiprocessing as mp  # use threads
from subprocess import check_output


def md5sum(filename):
    try:
        return check_output(["md5sum", filename]), None
    except Exception as e:
        return None, e


# Function to recursively fetch file path
def fetch_file_list(wd, file_list=[]):
    for item in os.listdir(wd):
        s = os.path.join(wd, item)
        if os.path.isdir(s):
            fetch_file_list(s, file_list)
        elif os.path.basename(s).endswith(('1.fq.gz', '2.fq.gz',
                                           '1.fastq.gz', '2.fastq.gz',
                                           '1.clean.fq.gz', '2.clean.fq.gz')):
            file_list.append(s)
            # match = re.search('^([^_]*)_([^_]*)_(\d)(\..*)$', os.path.basename(s))
            # if not any(fname.endswith('_fastqc.html') for fname in os.listdir(os.path.dirname(s))):
            #     file_list.append(s)
    return file_list


if __name__ == "__main__":

    # Input
    parent_dir = '/media/luolab/ZA1BT1ER/scRNAseq/yanting_all/data/yanting_200115/'
    files = fetch_file_list(parent_dir)

    # Output
    out_file = os.path.join(parent_dir, "md5sum.txt")

    number_of_processes = 32

    p = mp.Pool(number_of_processes)  # specify number of concurrent processes
    with open(out_file, "wb") as logfile:
        for output, error in p.imap(md5sum, files):  # provide filenames
            if error is None:
                logfile.write(output)
