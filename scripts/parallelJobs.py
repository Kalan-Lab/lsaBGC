import os
import sys
import multiprocessing
import argparse
import subprocess
import time

def create_parser():
    """ Parse arguments """
    parser = argparse.ArgumentParser(description="""
        Program to simplify running parallel jobs.
        """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--task_file', type=str, help="Task file. Each line represents an individual task/command.", required=True)
    parser.add_argument('-p', '--pool_size', type=int, help='Number of jobs to run in parallel. Default is 1.', required=False, default=1)
    args = parser.parse_args()
    return args

def mp_worker(input_cmd):
    input_cmd_split = input_cmd.split()
    print('Running\n%s\n' % input_cmd)
    proc = subprocess.Popen(' '.join(input_cmd_split), shell=True)
    proc.wait()

if __name__ == '__main__':
    myargs = create_parser()
    task_file = os.path.abspath(myargs.task_file)
    pool_size = myargs.pool_size

    tasks = []
    with open(task_file) as otf:
        for line in otf:
            line = line.strip()
            tasks.append(line)

    p = multiprocessing.Pool(pool_size)
    p.map(mp_worker, tasks)
