import time
import sys

def firstproc5(i):
    print('Start: first 5'+str(i))
    time.sleep(5)
    print('End: first 5'+str(i))


if __name__ == '__main__':
    firstproc5(sys.argv)
