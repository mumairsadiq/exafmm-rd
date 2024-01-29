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
        return f'FMM(P={self.p}, i={self.i}, l2p={self.l2p}, l2f={self.l2f}, kernel_time={self.kt}, total_time={self.tt})'
    
    def __repr__(self):
        return str(self)
    
    def same_setting(self,fmm):
        return self.p == fmm.p and self.i == fmm.i
    
    def get_setting(self):
        return (self.p, self.i)
    
    def __add__(self, other):
        if self.get_setting() == other.get_setting():
            return fmmdata(self.p,self.i,(self.l2p+other.l2p),(self.l2f+other.l2f),(self.kt+other.kt),(self.tt+other.tt))
        raise 'FMMs with different settings cannot be added'
    
    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)
    
    def __truediv__(self, n):
        return fmmdata(self.p,self.i,self.l2p/n,self.l2f/n,self.kt/n,self.tt/n)

    
def get_fmmdatas(filename : str):
    res = open(filename)
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
    return {'filename':filename,'data':fmms}

def collect_fmms_by_p(fmms):
    res = {}
    all_P = {fmm.p for fmm in fmms}
    for p in all_P:
        res[p] = [fmm for fmm in fmms if fmm.p == p]
    return res

def collect_fmms_by_i(fmms):
    res = {}
    all_i = {fmm.i for fmm in fmms}
    for i in all_i:
        res[i] = [fmm for fmm in fmms if fmm.i == i]
    return res

def average_fmms(fmms):
    res = {'filenames':[],'fmms':[]}
    fmm_grouped_by_setting = {}
    for fmm in fmms:
        datas = fmm['data']
        res['filenames'].append(fmm['filename'])
        for data in datas:
            setting = data.get_setting()
            if setting not in fmm_grouped_by_setting:
                fmm_grouped_by_setting[setting] = []
            fmm_grouped_by_setting[setting].append(data)
    for setting in fmm_grouped_by_setting:
        datas = fmm_grouped_by_setting[setting]
        num = len(datas)
        data_averaged = sum(datas) / num
        res['fmms'].append(data_averaged)
    return res

def plot_fmms(fmms,fmm_file_num,filename,keyname):
    color_codes = ['#002c53','#ffa510','#0c84c6','#628449','#f74d4d','#2455a4','#41b7ac']
    markers = ['-o','-*','-x','-d']
    plt.rcParams["font.family"] = "DejaVu Serif" 
    fig = plt.figure()
    fig.set_figheight(4)
    fig.set_figwidth(11)
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
    xs = []
    cnt = 0
    if keyname == 'P' or keyname == 'p':
        xname = 'images'
    else:
        xname = 'P'
    for key in fmms:
        c = color_codes[cnt % len(color_codes)]
        marker = markers[int(cnt / len(color_codes)) % len(markers)]
        fmm = fmms[key]
        if keyname == 'P' or keyname == 'p':
            x = [fmm.i for fmm in fmm]
        else:
            x = [fmm.p for fmm in fmm]
        if len(x) > len(xs):
            xs = x
        l2f = [fmm.l2f for fmm in fmm]
        kt = [fmm.kt for fmm in fmm]
        ax1.plot(x,l2f,marker,label=f'{keyname}={key}',color=c)
        ax2.plot(x,kt,marker,label=f'{keyname}={key}',color=c)
        cnt+=1
        
    ax1.set_xticks(xs)
    ax1.set_yscale('log')
    ax1.set_xlabel(xname)
    ax1.set_ylabel('L2-err')
    ax1.grid()
    ax1.set_title(f'L2 error of Force(avg_num={fmm_file_num})')
    ax1.legend()
    ax2.set_xticks(xs)
    ax2.set_xlabel(xname)
    ax2.set_ylabel('time(s)')
    ax2.legend()
    ax2.grid()
    ax2.set_title(f'kernel time(avg_num={fmm_file_num})')
    ax2.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

    fig.tight_layout()
    plt.savefig('exp(' + filename + ").png", dpi=250)

def get_fmms(filenames : list):
    return [get_fmmdatas(filename) for filename in filenames]

def get_averaged_fmms(filenames : list):
    return average_fmms(get_fmms(filenames))