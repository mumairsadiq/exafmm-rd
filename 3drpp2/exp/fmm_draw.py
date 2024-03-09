import re
from typing import List
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

plt.rcParams["font.family"] = "DejaVu Serif"
output_dir = 'fmm_draw_output'
   
class FmmDataParser:
    def __init__(self):
        self.setting_names = []
        self.setting_regxs = []
        self.inherence_names = []
        self.inherence_regxs = []
        self.result_names = []
        self.result_regxs = []

    def add_setting_rule(self, name : str, regx : str):
        self.setting_names.append(name)
        self.setting_regxs.append(regx)

    def add_inherence_rule(self, name : str, regx : str):
        self.inherence_names.append(name)
        self.inherence_regxs.append(regx)

    def add_result_rule(self, name : str, regx : str):
        self.result_names.append(name)
        self.result_regxs.append(regx)
        
    def parse(self, txt : str):
        inherences = [re.findall(regx,txt)[0] if re.findall(regx,txt) else "unknown" for regx in self.inherence_regxs]
        inherence = ''.join([f'{name}={value}' for name,value in zip(self.inherence_names, inherences)])
        try:
            setting_values = np.array([[float(d) for d in re.findall(regx,txt)] for regx in self.setting_regxs])
        except:
            raise 'cannot find setting!'
        try:
            result_values = np.array([[float(d) for d in re.findall(regx,txt)] for regx in self.result_regxs])
        except:
            raise 'cannot find result!'
        func = {tuple(setting) : result for setting,result in zip(np.transpose(setting_values),np.transpose(result_values))}
        return inherence, func, self.setting_names, self.result_names
        
class FmmData:
    '''essentially a vector function from settings to results, e.g. (P=4,i=2,r=0,m='asd') -> (l2p=1e-5,l2f=1e-6,kt=0.01,tt=0.02)'''
    def __init__(self, data, name, inherence, setting_names, result_names):
        self.data = data
        self.name = name
        self.inherence = inherence
        self.setting_names = setting_names
        self.result_names = result_names

    @classmethod
    def parsefromfilename(cls, filename : str, parser : FmmDataParser):
        res = open(filename)
        txt = res.read()
        inherence, func, setting_names, result_names = parser.parse(txt)
        return cls(func, filename, inherence, setting_names, result_names)
    
    def __str__(self) -> str:
        res = f'FmmData({self.name}) : \n'
        for key in self.data:
            res += str(key) + ' -> ' + str(self.data[key]) + '\n'
        return res
    
    def __repr__(self) -> str:
        return str(self)

def averaged_FmmData(fmmfiledatas : List[FmmData], is_forced):
    datalen = [len(fmmdatafile.data) for fmmdatafile in fmmfiledatas]
    if len(set(datalen)) > 1:
        if not is_forced:
            raise "Experimental data of different lengths cannot be averaged!"
        else:
            print(f"WARNING : averaging data with different lengths : {datalen}, this may cause unreasonable result!")
    ids = [fmmdatafile.inherence for fmmdatafile in fmmfiledatas]
    if len(set(ids)) > 1:
        if not is_forced:
            raise "Experimental data from different machines cannot be averaged!"
        else:
            print("WARNING : averaging data from different machines, this may cause unreasonable result!")
    elif 'unknown' in ids[0]:
        if not is_forced:
            raise "Experimental data from unknown machine cannot be averaged!"
        else:
            print("WARNING : averaging data from unknown machines, this may cause unreasonable result!")
    all_data = {}
    for fmmfiledata in fmmfiledatas:
        for setting in fmmfiledata.data:
            if setting not in all_data:
                all_data[setting] = []
            all_data[setting].append(fmmfiledata.data[setting])
    averaged_data = {}
    for setting in all_data:
        avergaed_fmmresult = sum(all_data[setting]) / len(all_data[setting])
        averaged_data[setting] = avergaed_fmmresult
    return FmmData(averaged_data, f'avg_num={len(fmmfiledatas)}', ids[0], fmmfiledatas[0].setting_names, fmmfiledatas[0].result_names)

# {marker : [x,y]}
def get_functions(fmmfiledata : FmmData, xname, yname):
    '''
    extract x->y functions from FmmData, 
    the remaining part of settings are used as function's marker
    '''
    data = fmmfiledata.data
    xidx = fmmfiledata.setting_names.index(xname)
    yidx = fmmfiledata.result_names.index(yname)
    marker_names = fmmfiledata.setting_names[:xidx]+fmmfiledata.setting_names[xidx+1:]
    res = {'marker_names': marker_names, 'xname' : xname, 'yname' : yname, 'func':{}}
    for setting in data:
        marker = setting[:xidx]+setting[xidx+1:]
        if marker not in res['func']:
            res['func'][marker] = [[],[]]
        res['func'][marker][0].append(setting[xidx])
        res['func'][marker][1].append(data[setting][yidx])
    return res

def draw_functions(funcs : dict, xlabel : str, ylabel : str, title : str, ylog = 0):
    plt.figure()
    xs = []
    for marker in funcs['func']:
        func = funcs['func'][marker]
        mark = ''.join([f'{n}={m},' for m,n in zip(marker,funcs['marker_names']) if n != ''])[:-1]
        plt.plot(func[0], func[1], '-o', label=mark)
        if len(func[0]) > len(xs):
            xs = func[0]
    plt.xticks(xs) 
    plt.legend()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid()
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    if ylog:
        plt.yscale('log')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    plt.savefig(f'{output_dir}/fmm_draw_{xlabel}-{title}.png',dpi=300)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--file', help='expdata file, or tabel of many expdata files')
    parser.add_argument('--forced', help='is average expdata forcily', default=0)
    parser.add_argument('--mass', help='is show all expdatas as semi-transparent', default=0)
    parser.add_argument('--truncate', help='truncation of files', default=sys.maxsize)
    args = vars(parser.parse_args())
    filename = args['file']
    forced_average = int(args['forced'])
    mass = int(args['mass'])
    truncate = int(args['truncate'])

    if filename:
        if filename[-4:] == '.out':
            filenames = [filename]
        elif filename[-4:] == '.txt':
            file = open(filename)
            filenames = file.read()
            filenames = filenames.splitlines()
            comments = [filename for filename in filenames if len(filename) > 0 and filename[0] == '#']
            filenames = [filename for filename in filenames if len(filename) > 0 and filename[0] != '#']
            print(*comments,sep='\n')
        else:
            raise 'file type not supported'
    else:
        raise 'No data file was given'
    filenames = filenames[:truncate]
    num = len(filenames)
    print(f"plot from {filenames}")

    fmm_data_parser = FmmDataParser()
    fmm_data_parser.add_inherence_rule('machine-id', r'machine-id=(.*)')
    fmm_data_parser.add_setting_rule('P', r'P                    = (.*)')
    fmm_data_parser.add_setting_rule('i', r'images               = (.*)')
    fmm_data_parser.add_setting_rule('c', r'cycle                = (.*)')
    fmm_data_parser.add_result_rule('l2p', r'L2  \(p\)  : (.*)L2  \(f\)')
    fmm_data_parser.add_result_rule('l2f', r'L2  \(f\)  : (.*)L2  \(e\)')
    fmm_data_parser.add_result_rule('kernel time', r'\[FMM_kernels time measured : (.*) seconds.\]\n')
    fmm_data_parser.add_result_rule('total time', r'\[FMM time measured : (.*) seconds.\]\n')
    fmm_data_parser.add_result_rule('hadamard product', r'\[hadamard_8x8 time measured : (.*) seconds.\]\n')
    fmm_data_parser.add_result_rule('P2P', r'\[P2P time measured : (.*) seconds.\]\n')
    fmm_data_parser.add_result_rule('M2L', r'\[M2L time measured : (.*) seconds.\]\n')

    fmms = [FmmData.parsefromfilename(filename,fmm_data_parser) for filename in filenames]
    fmmavg = averaged_FmmData(fmms, forced_average)

    p_tt = get_functions(fmmavg,'P','total time')
    p_ha = get_functions(fmmavg,'P','hadamard product')
    p_m2l = get_functions(fmmavg,'P','M2L')
    i_kt = get_functions(fmmavg,'i','kernel time')
    i_ha = get_functions(fmmavg,'i','hadamard product')
    i_l2f = get_functions(fmmavg,'i','l2f')
    i_kt = get_functions(fmmavg,'i','kernel time')
    i_m2l = get_functions(fmmavg,'i','M2L')
    draw_functions(funcs = i_kt, xlabel = 'images', ylabel = 'time(s)', title = f'kernel time(avg_num={num})')
    draw_functions(funcs = p_tt, xlabel = 'P', ylabel = 'time(s)', title = f'total time(avg_num={num})')
    draw_functions(funcs = p_ha, xlabel = 'P', ylabel = 'time(s)', title = f'hadamard_8x8 time(avg_num={num})')
    draw_functions(funcs = p_m2l, xlabel = 'P', ylabel = 'time(s)', title = f'M2L time(avg_num={num})')
    draw_functions(funcs = i_ha, xlabel = 'images', ylabel = 'time(s)', title = f'hadamard_8x8 time(avg_num={num})')
    draw_functions(funcs = i_l2f, xlabel = 'images', ylabel = 'L2 error', title = f'L2 error of force(avg_num={num})', ylog = 1)
    draw_functions(funcs = i_m2l, xlabel = 'images', ylabel = 'time(s)', title = f'M2L time(avg_num={num})')