#cmp.py
# Simple configuration script

# import m5
# from m5.objects import *
# import os, optparse, sys
# AddToPath('../')

from __future__ import print_function
from __future__ import absolute_import

import optparse
import sys
import os

import m5
from m5.defines import buildEnv
from m5.objects import *
from m5.util import addToPath, fatal, warn

addToPath('../')
# import Simulation
# from Caches import *
import Mybench

# Get paths we might need.  It's expected this file is in m5/configs/example.
config_path = os.path.dirname(os.path.abspath(__file__))
print(config_path)
config_root = os.path.dirname(config_path)+"/configs"
print(config_root)
m5_root = os.path.dirname(config_root)
print(m5_root)

parser = optparse.OptionParser()

# Benchmark options

parser.add_option("-b", "--benchmark", default="",
                 help="The benchmark to be loaded.")

parser.add_option("-c", "--chkpt", default="",
                 help="The checkpoint to load.")

execfile(os.path.join(config_root, "configs", "Options.py"))

(options, args) = parser.parse_args()

if args:
    print("Error: script doesn't take any positional arguments")
    sys.exit(1)

if options.benchmark == 'perlbench':
   process = Mybench.perlbench
elif options.benchmark == 'bzip2':
   process = Mybench.bzip2
elif options.benchmark == 'gcc':
   process = Mybench.gcc
elif options.benchmark == 'bwaves':
   process = Mybench.bwaves
elif options.benchmark == 'gamess':
   process = Mybench.gamess
elif options.benchmark == 'mcf':
   process = Mybench.mcf
elif options.benchmark == 'milc':
   process = Mybench.milc
elif options.benchmark == 'zeusmp':
   process = Mybench.zeusmp
elif options.benchmark == 'gromacs':
   process = Mybench.gromacs
elif options.benchmark == 'cactusADM':
   process = Mybench.cactusADM
elif options.benchmark == 'leslie3d':
   process = Mybench.leslie3d
elif options.benchmark == 'namd':
   process = Mybench.namd
elif options.benchmark == 'gobmk':
   process = Mybench.gobmk
elif options.benchmark == 'dealII':
   process = Mybench.dealII
elif options.benchmark == 'soplex':
   process = Mybench.soplex
elif options.benchmark == 'povray':
   process = Mybench.povray
elif options.benchmark == 'calculix':
   process = Mybench.calculix
elif options.benchmark == 'hmmer':
   process = Mybench.hmmer
elif options.benchmark == 'sjeng':
   process = Mybench.sjeng
elif options.benchmark == 'GemsFDTD':
   process = Mybench.GemsFDTD
elif options.benchmark == 'libquantum':
   process = Mybench.libquantum
elif options.benchmark == 'h264ref':
   process = Mybench.h264ref
elif options.benchmark == 'tonto':
   process = Mybench.tonto
elif options.benchmark == 'lbm':
   process = Mybench.lbm
elif options.benchmark == 'omnetpp':
   process = Mybench.omnetpp
elif options.benchmark == 'astar':
   process = Mybench.astar
elif options.benchmark == 'wrf':
   process = Mybench.wrf
elif options.benchmark == 'sphinx3':
   process = Mybench.sphinx3
elif options.benchmark == 'xalancbmk':
   process = Mybench.xalancbmk
elif options.benchmark == 'specrand_i':
   process = Mybench.specrand_i
elif options.benchmark == 'specrand_f':
   process = Mybench.specrand_f

if options.chkpt != "":
   process.chkpt = options.chkpt

(CPUClass, test_mem_mode, FutureClass) = Simulation.setCPUClass(options)

CPUClass.clock = '1.0GHz'

#np = options.num_cpus 
np = 1

system = System(cpu = [CPUClass(cpu_id=i) for i in xrange(np)],
                physmem = PhysicalMemory(range=AddrRange("4096MB")),
                membus = Bus(), mem_mode = 'timing')

system.physmem.port = system.membus.port

for i in xrange(np):   
    if options.caches:
        system.cpu[i].addPrivateSplitL1Caches(L1Cache(size = '64kB'),
                                              L1Cache(size = '64kB'))
    if options.l2cache:
        system.l2 = L2Cache(size='2MB')
        system.tol2bus = Bus()
        system.l2.cpu_side = system.tol2bus.port
        system.l2.mem_side = system.membus.port
        system.cpu[i].connectMemPorts(system.tol2bus)
    else:
        system.cpu[i].connectMemPorts(system.membus)
    system.cpu[i].workload = process[i]


root = Root(system = system)

Simulation.run(options, root, system, FutureClass)