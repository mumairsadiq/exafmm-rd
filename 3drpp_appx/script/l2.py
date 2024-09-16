import os
import sys
from rich import print
import matplotlib.pyplot as plt
import pyrtfmm # see issue-2

result_dir = "./result"


if __name__ == "__main__":
    begin = sys.argv[1].rfind('/')
    end = sys.argv[1].rfind('.')
    name = sys.argv[1][begin+1:end]
    
    log_files = open(sys.argv[1]).readlines()
    
    txt = ""
    for log_file in log_files:
        log_file = log_file.strip()
        if log_file != "":
            txt += open(log_file).read()
    
    for group in ["fvd", "fve", "dve"]:
        title = name + "_" + group
        data = pyrtfmm.parse(txt)[group]
        table = pyrtfmm.make_table(data,title)
        print(table)
        
        p = [d[0] for d in data]
        i = [d[1] for d in data]
        l2p = [d[2] for d in data]
        l2f = [d[3] for d in data]
        l2e = [d[4] for d in data]
        
        plt.figure()
        plt.plot(p, l2p, '-o', label='potential')
        plt.plot(p, l2f, '-o', label='force')
        plt.plot(p, l2e, '-o', label='energy')
        plt.yscale('log')
        plt.grid()
        plt.xlabel('P')
        plt.ylabel('$L_2 error$')
        plt.title(title)
        plt.legend()
        if not os.path.exists(result_dir):
            os.system(f"mkdir {result_dir}")
        plt.savefig(f"./{result_dir}/{title}.png",dpi=200)