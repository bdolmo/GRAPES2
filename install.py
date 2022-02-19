#!/usr/bin/env python3

import os
import sys
import subprocess
# import argparse
import yaml

curr_dir = os.getcwd()

def execute_setup():
    '''
        Execute setup.py install via pip
    '''

    msg = " INFO: Installing python packages via pip"
    print(msg)

    cmd = "pip3 install ."
    p1 = subprocess.run(cmd, shell=True, stdout=sys.stdout,
    stderr=subprocess.PIPE)

def build_cpp():
    '''
        Execute setup.py install via pip
    '''
    cpp_dirs = {
        'htslib' : curr_dir + '/bin/htslib',
        'seqlib' : curr_dir + '/bin/SeqLib',
        'targetdepth': curr_dir + '/bin/TargetDepth'
    }

    for project in cpp_dirs:
      os.chdir(cpp_dirs[project])

      if project is 'seqlib':
        cmd = "make clean"
        p1 = subprocess.run(cmd, shell=True, stdout=sys.stdout,
        stderr=subprocess.PIPE)

        cmd = "./configure"
        p1 = subprocess.run(cmd, shell=True, stdout=sys.stdout,
        stderr=subprocess.PIPE)

        cmd ="make CXXFLAGS=\"-std=c++11\""
        p1 = subprocess.run(cmd, shell=True, stdout=sys.stdout,
        stderr=subprocess.PIPE)

        cmd = "make install"
        p1 = subprocess.run(cmd, shell=True, stdout=sys.stdout,
        stderr=subprocess.PIPE)
      else:
        cmd = "make clean"
        p1 = subprocess.run(cmd, shell=True, stdout=sys.stdout,
        stderr=subprocess.PIPE)

        cmd = "make install"
        p1 = subprocess.run(cmd, shell=True, stdout=sys.stdout,
        stderr=subprocess.PIPE)


if __name__ == "__main__":
    build_cpp()
    sys.exit()
    execute_setup()
