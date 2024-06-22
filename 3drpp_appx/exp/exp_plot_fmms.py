import fmmdata
import sys
import matplotlib.pyplot as plt

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
    print(f"plot from {filenames}")

    num = len(filenames)
    fmms = [fmmdata.get_expdata(filename) for filename in filenames]
    fmm_avg = fmmdata.average_expdatas(fmms, forced_average)
    p_fmms_avg = fmmdata.collect_fmmdatas_by_p(fmm_avg)
    i_fmms_avg = fmmdata.collect_fmmdatas_by_i(fmm_avg)

    plt.rcParams["font.family"] = "DejaVu Serif" 

    fig, axs = plt.subplots(1,2)
    fig.set_figheight(4)
    fig.set_figwidth(12)
    if mass:
        for fmm in fmms:
            p_fmms = fmmdata.collect_fmmdatas_by_p(fmm)
            fmmdata.plot_fmms(p_fmms, num, 'P', axs, 0.1, 0)
    fmmdata.plot_fmms(p_fmms_avg, num, 'P', axs, 1, 1)
    axs[0].grid()
    axs[1].grid()
    fig.tight_layout()
    plt.savefig('exp(i-p).png', dpi=250)

    fig, axs = plt.subplots(1,2)
    fig.set_figheight(4)
    fig.set_figwidth(12)
    if mass:
        for fmm in fmms:
            i_fmms = fmmdata.collect_fmmdatas_by_i(fmm)
            fmmdata.plot_fmms(i_fmms, num, 'i', axs, 0.1, 0)
    fmmdata.plot_fmms(i_fmms_avg, num, 'i', axs, 1, 1)
    axs[0].grid()
    axs[1].grid()
    fig.tight_layout()
    plt.savefig('exp(p-i).png', dpi=250)