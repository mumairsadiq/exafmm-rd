import matplotlib.pyplot as plt
import re
from typing import List, Dict

class FmmData:
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
            return FmmData(self.p,self.i,(self.l2p+other.l2p),(self.l2f+other.l2f),(self.kt+other.kt),(self.tt+other.tt))
        raise 'FMMs with different settings cannot be added'
    
    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)
    
    def __truediv__(self, n):
        return FmmData(self.p,self.i,self.l2p/n,self.l2f/n,self.kt/n,self.tt/n)
    
class ExpData:
    def __init__(self, filename : str, datas : List[FmmData], machine_id : str):
        self.filename = filename
        self.datas = datas
        self.machine_id = machine_id

    def __str__(self) -> str:
        return f'filename={self.filename}, # of datas={self.datas}, machine_id={self.machine_id}'
    
    def __repr__(self) -> str:
        return str(self)

    
def get_expdata(filename : str) -> ExpData:
    res = open(filename)
    txt = res.read()
    data0 = re.findall(r'machine-id=(.*)',txt)
    data1 = re.findall(r'L2  \(p\)  : (.*)L2  \(f\)',txt)
    data2 = re.findall(r'L2  \(f\)  : (.*)L2  \(e\)',txt)
    data3 = re.findall(r'\[FMM_kernels time measured : (.*) seconds.\]\n',txt)
    data4 = re.findall(r'\[FMM time measured : (.*) seconds.\]\n',txt)
    data5 = re.findall(r'images               = (.*)', txt)
    data6 = re.findall(r'P                    = (.*)', txt)
    machine_id = 'unknown' if len(data0) == 0 else data0[0]
    L2ps = [eval(e) for e in data1]
    L2fs = [eval(e) for e in data2]
    kts = [eval(e) for e in data3]
    tts = [eval(e) for e in data4]
    imgs = [eval(fmm) for fmm in data5]
    ps = [eval(fmm) for fmm in data6]
    datas = [FmmData(p,i,l2p,l2f,kt,tt) for p,i,l2p,l2f,kt,tt in zip(ps,imgs,L2ps,L2fs,kts,tts)]
    return ExpData(filename, datas, machine_id)

def check_data_validity(expdatas : List[ExpData], forced : int):
    # check if data lengths are same
    datalen = [len(expdata.datas) for expdata in expdatas]
    if len(set(datalen)) > 1:
        if not forced:
            raise BaseException(f"Experimental data of different lengths : {datalen} cannot be averaged!")
        else:
            print(f"WARNING : averaging data with different lengths : {datalen}, this may cause unreasonable result!")
    # check if machines are same
    ids = [expdata.machine_id for expdata in expdatas]
    if len(set(ids)) > 1:
        if not forced:
            raise "Experimental data from different machines cannot be averaged!"
        else:
            print("WARNING : averaging data from different machines, this may cause unreasonable result!")
    elif 'unknown' in ids:
        if not forced:
            raise "Experimental data from unknown machine cannot be averaged!"
        else:
            print("WARNING : averaging data from unknown machines, this may cause unreasonable result!")

def average_expdatas(expdatas : List[ExpData], forced : int) -> ExpData:
    check_data_validity(expdatas, forced)
    averaged_data = []
    fmm_grouped_by_setting = {}
    for expdata in expdatas:
        datas = expdata.datas
        for data in datas:
            setting = data.get_setting()
            if setting not in fmm_grouped_by_setting:
                fmm_grouped_by_setting[setting] = []
            fmm_grouped_by_setting[setting].append(data)
    for setting in fmm_grouped_by_setting:
        datas = fmm_grouped_by_setting[setting]
        num = len(datas)
        data_averaged = sum(datas) / num
        averaged_data.append(data_averaged)
    return ExpData("averaged", averaged_data, "virutal")

def collect_fmmdatas_by_p(expdata : ExpData) -> Dict[int,List[FmmData]]:
    res = {}
    all_P = {fmmdata.p for fmmdata in expdata.datas}
    for p in all_P:
        res[p] = [fmmdata for fmmdata in expdata.datas if fmmdata.p == p]
    return res

def collect_fmmdatas_by_i(expdata : ExpData) -> Dict[int,List[FmmData]]:
    res = {}
    all_i = {fmmdata.i for fmmdata in expdata.datas}
    for i in all_i:
        res[i] = [fmmdata for fmmdata in expdata.datas if fmmdata.i == i]
    return res

def plot_fmms(fmms : Dict[int, List[FmmData]],fmm_file_num : int, keyname: str, axs, a : float, legend : int):
    color_codes = ['#002c53','#ffa510','#0c84c6','#628449','#f74d4d','#2455a4','#41b7ac']
    markers = ['-o','-*','-x','-d']
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
        lab = '' if not legend else f'{keyname}={key}'
        axs[0].plot(x,l2f,marker,label=lab,color=c, alpha=a)
        axs[1].plot(x,kt, marker,label=lab,color=c, alpha=a)
        cnt+=1
        
    axs[0].set_xticks(xs)
    axs[0].set_yscale('log')
    axs[0].set_xlabel(xname)
    axs[0].set_ylabel('L2-err')
    axs[0].set_title(f'L2 error of Force(avg_num={fmm_file_num})')
    
    axs[1].set_xticks(xs)
    axs[1].set_xlabel(xname)
    axs[1].set_ylabel('time(s)')
    axs[1].set_title(f'kernel time(avg_num={fmm_file_num})')
    axs[1].ticklabel_format(axis='y', style='sci', scilimits=(0,0))

    if legend:
        axs[0].legend()
        axs[1].legend()
