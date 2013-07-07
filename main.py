#!/bin/python

import time;
import sys;
import os;
import commands;
import re;
import threading;
import logging;
import string;
from matrix import *;
from pylab import *;

def main():
  logging.basicConfig(filename='result.log', level=logging.INFO, 
      format='%(asctime)s %(levelname)s: %(message)s')
  logging.info("============ start ============")

  vmNum = int(sys.argv[1])
  pdNum = int(sys.argv[2])
  topology = matrixUtil(vmNum, pdNum)

  logging.info("number of VMs = %d, number of physical disks = %d", vmNum, pdNum)
  logging.info("Generated topology:")
  logging.info(topology.data)

  # calculate steady-state distribution for each VM
  topology.steadyDistribute(1000, 0.01)

  # reArrange
  topology.reArrange()

  # calculate transient matrix for each VM
  topology.calTranMatrix()

  # simulate generated sequence
  topology.generateSeq()

  # evaluate transient based on cos
  topology.evalByCos()

if __name__ == '__main__':
  main()

