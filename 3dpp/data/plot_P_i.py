import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import figure
import re
import sys

class fmmdata:
    def __init__(self,p,i,l2p,l2f,kt,tt):
        self.p = p
        self.i = i
        self.l2p = l2p
        self.l2f = l2f
        self.kt = kt
        self.tt = tt
    def __str__(self) -> str:
        return f'P={self.p}, i={self.i}, l2p={self.l2p}, l2f={self.l2f}, kernel_time={self.kt}, total_time={self.tt}'
    
    def __repr__(self):
        return str(self)

def get_rtfmm_data():
    res = open("slurm-152683.out") # test image
    txt = res.read()
    data1 = re.findall(r'L2  \(p\)  : (.*)L2  \(f\)',txt)
    data2 = re.findall(r'L2  \(f\)  : (.*)L2  \(e\)',txt)
    data3 = re.findall(r'\[FMM_kernels time measured : (.*) seconds.\]\n',txt)
    data4 = re.findall(r'\[FMM time measured : (.*) seconds.\]\n',txt)
    data5 = re.findall(r'images               = (.*)', txt)
    data6 = re.findall(r'P                    = (.*)', txt)
    L2ps = [eval(e) for e in data1]
    L2fs = [eval(e) for e in data2]
    kts = [eval(e) for e in data3]
    tts = [eval(e) for e in data4]
    imgs = [eval(fmm) for fmm in data5]
    ps = [eval(fmm) for fmm in data6]
    fmms = [fmmdata(p,i,l2p,l2f,kt,tt) for p,i,l2p,l2f,kt,tt in zip(ps,imgs,L2ps,L2fs,kts,tts)]
    return fmms

def collect_fmms_by_p(fmms):
    res = {}
    all_P = {fmm.p for fmm in fmms}
    for p in all_P:
        res[p] = [fmm for fmm in fmms if fmm.p == p]
    return res


fmms = get_rtfmm_data()
p_fmms = collect_fmms_by_p(fmms)

color_codes = ['#002c53','#ffa510','#0c84c6','#628449','#f74d4d','#2455a4','#41b7ac']
#color_codes = ['#3b6291','#943c39','#779043','#624c7c','#388498','#bf7334','#3f6899','#9c403d','#7d9847','#675083','#3b8ba1','#c97937']

plt.rcParams["font.family"] = "DejaVu Serif" 
fig = plt.figure()
fig.set_figheight(4)
fig.set_figwidth(11)
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)
images = []
cnt = 0
for p in p_fmms:
    print(p)
    print(*p_fmms[p],sep='\n')
    c = color_codes[cnt % len(color_codes)]
    p_fmm = p_fmms[p]
    i = [fmm.i for fmm in p_fmm]
    if len(i) > len(images):
        images = i
    l2f = [fmm.l2f for fmm in p_fmm]
    kt = [fmm.kt for fmm in p_fmm]
    ax1.plot(i,l2f,'-o',label=f'P={p}',color=c)
    ax2.plot(i,kt,'-o',label=f'P={p}',color=c)
    cnt+=1
    
ax1.set_xticks(images)
ax1.set_yscale('log')
ax1.set_xlabel('images')
ax1.set_ylabel('L2-err')
ax1.grid()
ax1.set_title('L2 force')
ax1.legend()
ax2.set_xticks(images)
ax2.set_xlabel('images')
ax2.set_ylabel('time(s)')
ax2.legend()
ax2.grid()
ax2.set_title('kernel time')
ax2.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

fig.tight_layout()
plt.savefig("p_i.png", dpi=300)
