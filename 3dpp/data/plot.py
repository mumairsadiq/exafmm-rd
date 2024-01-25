import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import figure
import re

def get_exafmmt_data():
    #res = open("/home/jiamian/exafmm-t/exafmm-t/tests/slurm-147765.out") # -T 8
    res = open("/home/jiamian/exafmm-t/exafmm-t/tests/slurm-147767.out") # -T 12
    txt = res.read()
    data1 = re.findall(r'Potential Error L2   : (.*)\n',txt)
    data2 = re.findall(r'Gradient Error L2    : (.*)\n',txt)
    data3 = re.findall(r'Evaluation           : (.*)\n',txt)
    data4 = re.findall(r'Total                : (.*)\n',txt)
    L2p = [eval(e) for e in data1]
    L2f = [eval(e) for e in data2]
    L2kt = [eval(e) for e in data3]
    L2tt = [eval(e) for e in data4]
    fmmdata = re.findall(r'P                    : (.*)', txt)
    p = [eval(fmm) for fmm in fmmdata]
    return p, L2p, L2f, L2kt, L2tt

def get_rtfmm_data():
    #res = open("slurm-147766.out") # --th_num 8 sse
    #res = open("slurm-147776.out") # --th_num 12 sse
    #res = open("slurm-148813.out") # --th_num 8 avx
    res = open("slurm-148810.out") # --th_num 12 avx
    txt = res.read()
    data1 = re.findall(r'L2  \(p\)  : (.*)L2  \(f\)',txt)
    data2 = re.findall(r'L2  \(f\)  : (.*)L2  \(e\)',txt)
    data3 = re.findall(r'\[FMM_kernels time measured : (.*) seconds.\]\n',txt)
    data4 = re.findall(r'\[FMM time measured : (.*) seconds.\]\n',txt)
    L2p = [eval(e) for e in data1]
    L2f = [eval(e) for e in data2]
    L2kt = [eval(e) for e in data3]
    L2tt = [eval(e) for e in data4]
    fmmdata = re.findall(r'  P                    = (.*)', txt)
    p = [eval(fmm) for fmm in fmmdata]
    return p, L2p, L2f, L2kt, L2tt

p, L2p_exafmmt, L2f_exafmmt, L2kt_exafmmt, L2tt_exafmmt = get_exafmmt_data()
p2, L2p_rtfmmt, L2f_rtfmmt, L2kt_rtfmmt, L2tt_rtfmmt = get_rtfmm_data()

figure()
plt.plot(p,L2p_exafmmt,'-d',label='p-exafmmt')
plt.plot(p,L2f_exafmmt,'-d',label='f-exafmmt')
plt.plot(p2,L2p_rtfmmt,'-d',label='p-rtfmm')
plt.plot(p2,L2f_rtfmmt,'-d',label='f-rtfmm')
plt.xticks(p)
plt.yscale('log')
plt.xlabel('P')
plt.ylabel('L2-err')
plt.legend()
plt.savefig("P_L2.png", dpi=300)

figure()
plt.plot(p,L2kt_exafmmt,'-d',label='kernel time(exafmm-t)')
plt.plot(p2,L2kt_rtfmmt,'-d',label='kernel time(rtfmm)')
plt.xticks(p)
plt.xlabel('P')
plt.ylabel('kernel time(s)')
plt.legend()
plt.savefig("P_kernel-time.png", dpi=300)

figure()
plt.plot(p,L2tt_exafmmt,'-d',label='total time(exafmm-t)')
plt.plot(p2,L2tt_rtfmmt,'-d',label='total time(rtfmm)')
plt.xticks(p)
plt.xlabel('P')
plt.ylabel('total time(s)')
plt.legend()
plt.savefig("P_total-time.png", dpi=300)