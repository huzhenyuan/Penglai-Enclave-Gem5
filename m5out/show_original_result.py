def cal_cycles(filename):
    f = open(filename)
    f.seek(0, 2)
    reverse_line = ''
    total_cycles = 0
    while True:
        f.seek(-2, 1)
        r_char = f.read(1)
        if r_char != '\n':
            reverse_line += r_char
        else:
            line = ''
            for char in reverse_line[::-1]:
                line += char
            reverse_line = ''
            line=line.strip('\n').split(' ')
            if ((len(line) == 10) and (line[len(line) - 3] == "write")):
                access_cycles = long(line[0][0:len(line[0])-1])
                extra_cycles = long(line[4])
                total_cycles = access_cycles + extra_cycles
                print(filename, total_cycles)
                break
    return total_cycles

def main():
    cal_cycles("./result/simout-GemsFDTD")
    cal_cycles("./result/simout-gobmk")
    cal_cycles("./result/simout-lbm")
    cal_cycles("./result/simout-mcf")
    cal_cycles("./result/simout-milc")
    cal_cycles("./result/simout-sjeng")
    cal_cycles("./result/simout-zeusmp")

main()
