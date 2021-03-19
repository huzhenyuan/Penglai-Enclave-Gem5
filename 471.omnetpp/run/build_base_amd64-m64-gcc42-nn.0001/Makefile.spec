TUNE=base
EXT=amd64-m64-gcc42-nn
NUMBER=471
NAME=omnetpp
SOURCES= EtherAppCli.cc EtherAppCli_n.cc EtherAppSrv.cc EtherAppSrv_n.cc \
	 EtherApp_m.cc EtherBus.cc EtherBus_n.cc EtherCtrl_m.cc EtherEncap.cc \
	 EtherEncap_n.cc EtherFrame_m.cc EtherHost_n.cc EtherHub.cc EtherHub_n.cc \
	 EtherLLC.cc EtherLLC_n.cc EtherMAC.cc EtherMAC_n.cc EtherSwitch_n.cc \
	 LargeNet_n.cc MACAddress.cc MACAddress_m.cc MACRelayUnitBase.cc \
	 MACRelayUnitNP.cc MACRelayUnitNP_n.cc MACRelayUnitPP.cc \
	 MACRelayUnitPP_n.cc MACRelayUnit_n.cc Networks_n.cc eth-index_n.cc \
	 utils.cc libs/cmdenv/cmdenv.cc libs/cmdenv/enumstr.cc \
	 libs/cmdenv/heap.cc libs/envir/akoutvectormgr.cc libs/envir/args.cc \
	 libs/envir/cenvir.cc libs/envir/cinifile.cc libs/envir/filemgrs.cc \
	 libs/envir/main.cc libs/envir/omnetapp.cc libs/envir/patmatch.cc \
	 libs/envir/platdep.cc libs/envir/seeds.cc libs/envir/slaveapp.cc \
	 libs/envir/speedmtr.cc libs/sim/carray.cc libs/sim/cexception.cc \
	 libs/sim/cmessage.cc libs/sim/cpar.cc libs/sim/ctypes.cc \
	 libs/sim/task.cc libs/sim/cchannel.cc libs/sim/cfsm.cc \
	 libs/sim/cmodule.cc libs/sim/cpsquare.cc libs/sim/cvarhist.cc \
	 libs/sim/util.cc libs/sim/ccoroutine.cc libs/sim/cgate.cc \
	 libs/sim/cmsgheap.cc libs/sim/cqueue.cc libs/sim/cwatch.cc \
	 libs/sim/cdensity.cc libs/sim/chead.cc libs/sim/cnetmod.cc \
	 libs/sim/csimul.cc libs/sim/distrib.cc libs/sim/cdetect.cc \
	 libs/sim/chist.cc libs/sim/cobject.cc libs/sim/cstat.cc \
	 libs/sim/errmsg.cc libs/sim/cdispstr.cc libs/sim/cksplit.cc \
	 libs/sim/coutvect.cc libs/sim/cstruct.cc libs/sim/onstartup.cc \
	 libs/sim/cenum.cc libs/sim/cllist.cc libs/sim/cpacket.cc \
	 libs/sim/ctopo.cc libs/sim/random.cc libs/sim/std/netpack.cc \
	 libs/spec/spec_qsort.cc
EXEBASE=omnetpp
NEED_MATH=yes
BENCHLANG=CXX
ONESTEP=
CXXONESTEP=

BENCH_FLAGS      = -I. -Iomnet_include -Ilibs/envir
CC               = /usr/bin/gcc
COPTIMIZE        = -O2
CXX              = /usr/bin/g++
CXXOPTIMIZE      = -O2 
FC               = /usr/bin/gfortran
FOPTIMIZE        = -O2
FPBASE           = yes
OS               = unix
PORTABILITY      = -DSPEC_CPU_LP64
abstol           = 1e-06
action           = validate
allow_extension_override = 0
backup_config    = 1
baseexe          = omnetpp
basepeak         = 0
benchdir         = benchspec
benchmark        = 471.omnetpp
binary           = 
bindir           = exe
calctol          = 0
changedmd5       = 0
check_integrity  = 1
check_md5        = 1
check_version    = 1
command_add_redirect = 0
commanderrfile   = speccmds.err
commandexe       = omnetpp_base.amd64-m64-gcc42-nn
commandfile      = speccmds.cmd
commandoutfile   = speccmds.out
commandstdoutfile = speccmds.stdout
compareerrfile   = compare.err
comparefile      = compare.cmd
compareoutfile   = compare.out
comparestdoutfile = compare.stdout
compile_error    = 0
compwhite        = 
configdir        = config
configpath       = /mnt/spec-cpu2006/config/default.cfg
copies           = 1
datadir          = data
delay            = 0
deletebinaries   = 0
deletework       = 0
difflines        = 10
dirprot          = 511
endian           = 12345678
env_vars         = 0
exitvals         = spec_exit
expand_notes     = 0
expid            = 
ext              = amd64-m64-gcc42-nn
fake             = 0
feedback         = 1
flag_url_base    = http://www.spec.org/auto/cpu2006/flags/
floatcompare     = 
help             = 0
http_proxy       = 
http_timeout     = 30
hw_avail         = Dec-9999
hw_cpu_char      = 
hw_cpu_mhz       = 3000
hw_cpu_name      = AMD Opteron 256
hw_disk          = SATA
hw_fpu           = Integrated
hw_memory        = 2 GB (2 x 1GB DDR333 CL2.5)
hw_model         = Tyan Thunder KKQS Pro (S4882)
hw_nchips        = 1
hw_ncores        = 2
hw_ncoresperchip = 2
hw_ncpuorder     = 1 chip
hw_nthreadspercore = 1
hw_ocache        = None
hw_pcache        = 64 KB I + 64 KB D on chip per chip
hw_scache        = 1 MB I+D on chip per chip
hw_tcache        = None
hw_vendor        = Tyan
ignore_errors    = yes
ignore_sigint    = 0
ignorecase       = 
info_wrap_columns = 50
inputdir         = input
iteration        = -1
iterations       = 3
license_num      = 9999
line_width       = 0
locking          = 1
log              = CPU2006
log_line_width   = 0
logname          = /mnt/spec-cpu2006/result/CPU2006.003.log
lognum           = 003
mach             = default
mail_reports     = all
mailcompress     = 0
mailmethod       = smtp
mailport         = 25
mailserver       = 127.0.0.1
mailto           = 
make             = specmake
make_no_clobber  = 0
makeflags        = 
max_active_compares = 0
mean_anyway      = 0
min_report_runs  = 3
minimize_builddirs = 0
minimize_rundirs = 0
name             = omnetpp
nc               = 0
need_math        = yes
no_input_handler = close
no_monitor       = 
notes0100        =  C base flags: -O2
notes0110        =  C++ base flags: -O2 
notes0120        =  Fortran base flags: -O2
notes25          =  PORTABILITY=-DSPEC_CPU_LP64 is applied to all benchmarks in base.
notes_wrap_columns = 0
notes_wrap_indent =     
num              = 471
obiwan           = 
os_exe_ext       = 
output           = asc
output_format    = asc, pdf, Screen
output_root      = 
outputdir        = output
path             = /mnt/spec-cpu2006/benchspec/CPU2006/471.omnetpp
plain_train      = 1
prefix           = 
prepared_by      = 
rate             = 0
rawfile          = 
rawformat        = 0
realuser         = your name here
rebuild          = 0
reftime          = reftime
reltol           = 1e-05
reportable       = 0
resultdir        = result
review           = 0
run              = all
runspec          = ./bin/runspec --noreportable --size=ref int
safe_eval        = 1
section_specifier_fatal = 1
sendmail         = /usr/sbin/sendmail
setpgrp_enabled  = 1
setprocgroup     = 1
shrate           = 0
sigint           = 2
size             = ref
skipabstol       = 
skipobiwan       = 
skipreltol       = 
skiptol          = 
smarttune        = base
specdiff         = specdiff
specmake         = Makefile.YYYtArGeTYYYspec
specrun          = specinvoke
speed            = 0
srcalt           = 
srcdir           = src
stagger          = 10
strict_rundir_verify = 1
subworkdir       = work
sw_auto_parallel = No
sw_avail         = Dec-9999
sw_base_ptrsize  = 64-bit
sw_compiler      = gcc , g++ & gfortran 4.2.0325 (for AMD64)
sw_file          = ext3
sw_os            = SUSE SLES9 (for AMD64)
sw_other         = None
sw_peak_ptrsize  = Not Applicable
sw_state         = runlevel 3
sysinfo_program  = 
table            = 1
teeout           = yes
teerunout        = yes
test_date        = Dec-9999
test_sponsor     = Turbo Computers
tester           = 
top              = /mnt/spec-cpu2006
tune             = base
uid              = 0
unbuffer         = 1
update-flags     = 0
use_submit_for_speed = 0
username         = root
vendor           = anon
vendor_makefiles = 0
verbose          = 5
version          = 0
version_url      = http://www.spec.org/auto/cpu2006/current_version
workdir          = run
worklist         = list
