import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import figure
import re


def get_exafmmt_data(step,n=20):
    #res = open("/home/jiamian/exafmm-t/exafmm-t/tests/slurm-147765.out") # -T 8
    #res = open("/home/jiamian/exafmm-t/exafmm-t/tests/slurm-147767.out") # -T 12
    res = open("/home/jiamian/exafmm-t/exafmm-t/tests/slurm-149351.out") # -T 8 ad4_1
    txt = res.read()
    data1 = re.findall(step + r'                  : (.*)\n',txt)[:n]
    value = [eval(e) for e in data1]
    fmmdata = re.findall(r'P                    : (.*)', txt)[:n]
    p = [eval(fmm) for fmm in fmmdata]
    return p, value

def get_rtfmm_data(step,n=20):
    #res = open("slurm-147766.out") # --th_num 8 sse
    #res = open("slurm-148810.out") # --th_num 12 avx
    #res = open("slurm-148813.out") # --th_num 8 avx
    #res = open("slurm-149348.out") # --th_num 8 avx t
    #res = open("slurm-149361.out") # --th_num 8 avx t 2
    res = open("slurm-149364.out") # --th_num 8 avx t 3
    txt = res.read()
    data1 = re.findall(step + r' time measured : (.*) seconds.\]',txt)[:n]
    value = [eval(e) for e in data1]
    fmmdata = re.findall(r'  P                    = (.*)', txt)[:n]
    p = [eval(fmm) for fmm in fmmdata]
    return p, value

fig = plt.figure()
fig.set_figheight(7)
fig.set_figwidth(15)
stepnames = ['P2P', 'P2M', 'M2M', 'M2L', 'L2L', 'L2P']
cnt = 1
for stepname in stepnames:
    p1, v1 = get_exafmmt_data(stepname)
    p2, v2 = get_rtfmm_data(stepname)

    '''figure()
    plt.plot(p1,v1,'-d',label='exafmmt')
    plt.plot(p2,v2,'-d',label='rtfmm')
    plt.xticks(p1)
    plt.xlabel('P')
    plt.ylabel('time(s)')
    plt.title(stepname)
    plt.legend()
    plt.savefig("compare_" + stepname + ".png", dpi=300)'''

    ax = fig.add_subplot(2, 3, cnt)
    ax.plot(p1,v1,'-d',label='exafmmt')
    ax.plot(p2,v2,'-d',label='rtfmm')
    ax.set_xlabel('P')
    ax.set_ylabel('time(s)')
    ax.set_xticks(p1)
    ax.legend()
    ax.set_title(stepname)
    cnt+=1

fig.tight_layout()
fig.savefig("compare_allsteps.png", dpi=300)