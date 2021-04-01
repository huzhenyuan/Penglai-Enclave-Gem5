build/X86/gem5.opt --debug-flags=MYTRACE -r --stdout-file=simout-GemsFDTD configs/example/se.py -n 1 --caches --cpu-type=TimingSimpleCPU --mem-size=2GB --mem-channels=1 -c ./459.GemsFDTD/exe/GemsFDTD_base.amd64-m64-gcc42-nn
build/RISCV/gem5.opt --debug-flag=MYTRACE -r --stdout-file=simout-gobmk configs/example/se.py -n 1 --caches --cpu-type=TimingSimpleCPU --mem-size=2GB --mem-channels=1 -c ./445.gobmk/exe/gobmk_base.amd64-m64-gcc41-nn < ./445.gobmk/data/test/input/capture.tst
build/RISCV/gem5.opt --debug-flag=MYTRACE -r --stdout-file=simout-lbm configs/example/se.py -n 1 --caches --cpu-type=TimingSimpleCPU --mem-size=2GB --mem-channels=1 -c ./470.lbm/exe/lbm_base.amd64-m64-gcc41-nn -o '20 reference.bat 0 1 ./470.lbm/data/test/input/100_100_130_cf_a.of'
build/RISCV/gem5.opt --debug-flag=MYTRACE -r --stdout-file=simout-mcf configs/example/se.py -n 1 --caches --cpu-type=TimingSimpleCPU --mem-size=2GB --mem-channels=1 -c ./429.mcf/exe/mcf_base.amd64-m64-gcc41-nn -o ./429.mcf/data/train/input/inp.in
build/RISCV/gem5.opt --debug-flag=MYTRACE -r --stdout-file=simout-milc configs/example/se.py -n 1 --caches --cpu-type=TimingSimpleCPU --mem-size=2GB --mem-channels=1 -c ./433.milc/exe/milc_base.amd64-m64-gcc41-nn < ./433.milc/data/train/input/su3imp.in
build/RISCV/gem5.opt --debug-flag=MYTRACE -r --stdout-file=simout-sjeng configs/example/se.py -n 1 --caches --cpu-type=TimingSimpleCPU --mem-size=2GB --mem-channels=1 -c ./433.milc/exe/milc_base.amd64-m64-gcc41-nn < ./433.milc/data/train/input/su3imp.in
build/X86/gem5.opt --debug-flags=MYTRACE -r --stdout-file=simout-zeusmp configs/example/se.py -n 1 --caches --cpu-type=TimingSimpleCPU --mem-size=2GB --mem-channels=1 -c ./434.zeusmp/exe/zeusmp_base.amd64-m64-gcc42-nn
cd m5out
python ./show_mmt_result.py
cd ../