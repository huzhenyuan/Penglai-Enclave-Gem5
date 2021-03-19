TUNE=base
EXT=amd64-m64-gcc41-nn
NUMBER=401
NAME=bzip2
SOURCES= spec.c blocksort.c bzip2.c bzlib.c compress.c crctable.c \
	 decompress.c huffman.c randtable.c
EXEBASE=bzip2
NEED_MATH=
BENCHLANG=C
ONESTEP=
CONESTEP=

CC               = riscv64-unknown-linux-gnu-gcc
COPTIMIZE        = -O2 
CXX              = riscv64-unknown-linux-gnu-g++
CXXOPTIMIZE      = -O2 
FC               = riscv64-unknown-linux-gnu-gfortran
FOPTIMIZE        = -O2
FPBASE           = yes
OS               = unix
PORTABILITY      = -DSPEC_CPU_LP64
abstol           = 
action           = build
allow_extension_override = 0
backup_config    = 1
baseexe          = bzip2
basepeak         = 0
benchdir         = benchspec
benchmark        = 401.bzip2
binary           = 
bindir           = exe
calctol          = 0
changedmd5       = 0
check_integrity  = 1
check_md5        = 1
check_version    = 1
command_add_redirect = 0
commanderrfile   = speccmds.err
commandexe       = bzip2_base.amd64-m64-gcc41-nn
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
configpath       = /home/clive/workspace/luxu/spec-cpu2006_backup/config/riscv.cfg
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
ext              = amd64-m64-gcc41-nn
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
hw_ncores        = 1
hw_ncoresperchip = 1
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
logname          = /home/clive/workspace/luxu/spec-cpu2006_backup/result/CPU2006.051.log
lognum           = 051
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
name             = bzip2
nc               = 0
need_math        = 
no_input_handler = close
no_monitor       = 
notes0100        =  C base flags: -O2
notes0110        =  C++ base flags: -O2
notes0120        =  Fortran base flags: -O2
notes25          =  PORTABILITY=-DSPEC_CPU_LP64 is applied to all benchmarks in base.
notes_wrap_columns = 0
notes_wrap_indent =     
num              = 401
obiwan           = 
os_exe_ext       = 
output           = asc
output_format    = asc
output_root      = 
outputdir        = output
path             = /home/clive/workspace/luxu/spec-cpu2006_backup/benchspec/CPU2006/401.bzip2
plain_train      = 1
prefix           = 
prepared_by      = 
rate             = 0
rawfile          = 
rawformat        = 0
realuser         = your name here
rebuild          = 0
reftime          = reftime
reltol           = 
reportable       = 1
resultdir        = result
review           = 0
run              = all
runspec          = /home/clive/workspace/luxu/spec-cpu2006_backup/bin/runspec --action=build --config=riscv.cfg --tune=base bzip2
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
sw_compiler      = gcc , g++ & gfortran 4.1.0 (for AMD64)
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
top              = /home/clive/workspace/luxu/spec-cpu2006_backup
tune             = base
uid              = 1000
unbuffer         = 1
update-flags     = 0
use_submit_for_speed = 0
username         = clive
vendor           = anon
vendor_makefiles = 0
verbose          = 5
version          = 0
version_url      = http://www.spec.org/auto/cpu2006/current_version
workdir          = run
worklist         = list
