import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import figure
import re
import sys

def get_exafmmt_data(n=20):
    #res = open("/home/jiamian/exafmm-t/exafmm-t/tests/slurm-147765.out") # -T 8
    #res = open("/home/jiamian/exafmm-t/exafmm-t/tests/slurm-147767.out") # -T 12
    #res = open("/home/jiamian/exafmm-t/exafmm-t/tests/slurm-149351.out") # -T 8 ad4_1
    res = open("/home/jiamian/exafmm-t/exafmm-t/tests/slurm-152666.out") # -T 8 rtx6000-ada_1
    txt = res.read()
    data1 = re.findall(r'Potential Error L2   : (.*)\n',txt)[:n]
    data2 = re.findall(r'Gradient Error L2    : (.*)\n',txt)[:n]
    data3 = re.findall(r'Evaluation           : (.*)\n',txt)[:n]
    data4 = re.findall(r'Total                : (.*)\n',txt)[:n]
    L2p = [eval(e) for e in data1]
    L2f = [eval(e) for e in data2]
    L2kt = [eval(e) for e in data3]
    L2tt = [eval(e) for e in data4]
    fmmdata = re.findall(r'P                    : (.*)', txt)[:n]
    p = [eval(fmm) for fmm in fmmdata]
    return p, L2p, L2f, L2kt, L2tt

def get_rtfmm_data(n=20):
    #res = open("slurm-147766.out") # --th_num 8 sse
    #res = open("slurm-147776.out") # --th_num 12 sse
    #res = open("slurm-148813.out") # --th_num 8 avx
    #res = open("slurm-148810.out") # --th_num 12 avx
    #res = open("slurm-149348.out") # --th_num 8 avx t
    #res = open("slurm-149361.out") # --th_num 8 avx t 2
    #res = open("slurm-149364.out") # --th_num 8 avx t 3
    res = open("slurm-152673.out") # --th_num 8 avx t rtx6000-ada_1
    txt = res.read()
    data1 = re.findall(r'L2  \(p\)  : (.*)L2  \(f\)',txt)[:n]
    data2 = re.findall(r'L2  \(f\)  : (.*)L2  \(e\)',txt)[:n]
    data3 = re.findall(r'\[FMM_kernels time measured : (.*) seconds.\]\n',txt)[:n]
    data4 = re.findall(r'\[FMM time measured : (.*) seconds.\]\n',txt)[:n]
    L2p = [eval(e) for e in data1]
    L2f = [eval(e) for e in data2]
    L2kt = [eval(e) for e in data3]
    L2tt = [eval(e) for e in data4]
    fmmdata = re.findall(r'  P                    = (.*)', txt)[:n]
    p = [eval(fmm) for fmm in fmmdata]
    return p, L2p, L2f, L2kt, L2tt

n = 20
if len(sys.argv) > 1:
    n = int(sys.argv[1])
p, L2p, L2f, L2kt, L2tt = get_rtfmm_data(n)

figure()
plt.plot(p,L2p,'-d',label='p')
plt.plot(p,L2f,'-d',label='f')
plt.xticks(p)
plt.yscale('log')
plt.xlabel('P')
plt.ylabel('L2-err')
plt.legend()
plt.savefig("P_L2.png", dpi=300)

figure()
plt.plot(p,L2kt,'-d',label='kernel time')
plt.xticks(p)
plt.xlabel('P')
plt.ylabel('kernel time(s)')
plt.legend()
plt.savefig("P_kernel-time.png", dpi=300)

figure()
plt.plot(p,L2tt,'-d',label='total time')
plt.xticks(p)
plt.xlabel('P')
plt.ylabel('total time(s)')
plt.legend()
plt.savefig("P_total-time.png", dpi=300)