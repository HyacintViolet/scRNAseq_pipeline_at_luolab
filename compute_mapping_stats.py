import sys
import pandas as pd
import numpy as np

whitelist = pd.read_csv('whitelist.txt', sep="\t", header=None)

whitelist.columns = ["cell", "candidate", "cell_counts", "candidate_counts"]



def if __name__ == '__main__':
