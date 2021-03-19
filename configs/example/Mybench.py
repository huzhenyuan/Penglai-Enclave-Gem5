mcf = LiveProcess()
mcf.executable =  '/home/tiger/workspace/gem5/429.mcf/exe/'+'mcf_base.amd64-m64-gcc41-nn'
data='/home/tiger/workspace/gem5/429.mcf/data/test/input/inp.in'
mcf.cmd = [mcf.executable] + [data]
mcf.output = 'inp.out'