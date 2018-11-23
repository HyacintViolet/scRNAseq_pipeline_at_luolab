import os, re
import pandas as pd


whitelist = pd.read_csv('whitelist.txt', sep="\t", header=None)

whitelist.columns = ["cell", "candidate", "cell_counts", "candidate_counts"]


# Not sure how to use this yet. Deal later.
# def if __name__ == '__main__':
