import fmmdata
import sys

if __name__ == '__main__':
    filenames = ['slurm-152683.out']
    if len(sys.argv) > 1:
        argname = sys.argv[1]
        if argname[-4:] == '.out':
            filenames = [argname]
        elif argname[-4:] == '.txt':
            file = open(argname)
            filenames = file.read()
            filenames = filenames.splitlines()
            comments = [filename for filename in filenames if len(filename) > 0 and filename[0] == '#']
            filenames = [filename for filename in filenames if len(filename) > 0 and filename[0] != '#']
            print(*comments,sep='\n')
        else:
            raise 'file type not supported'
    print(f"plot from {filenames}")
    fmm_avg = fmmdata.get_averaged_fmms(filenames)
    num = len(fmm_avg['filenames'])
    p_fmms = fmmdata.collect_fmms_by_p(fmm_avg['fmms'])
    i_fmms = fmmdata.collect_fmms_by_i(fmm_avg['fmms'])
    fmmdata.plot_fmms(p_fmms, num, 'i-p', 'P')
    fmmdata.plot_fmms(i_fmms, num, 'p-i', 'i')