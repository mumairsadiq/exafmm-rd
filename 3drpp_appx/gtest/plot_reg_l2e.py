import re
import matplotlib.pyplot as plt

fig = plt.figure()
fig.set_figheight(8)
fig.set_figwidth(8)

p=6

files = [
    f'exp/debug0629_p={p}_seed=5.log',
    f'exp/debug0629_p={p}_seed=123.log'
]
labels = [
    'seed=5',
    'seed=123'
]

ax1 = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(2,1,2)

for file, label in zip(files,labels):
    res = open(file)
    txt = res.read()
    regdata = re.findall(r'rega                 = (.*)\n',txt)
    l2fdata = re.findall(r'L2  \(f\)  : (.*)L2.*\n',txt)
    l2edata = re.findall(r'L2  \(e\)  : (.*)\n',txt)
    reg = [eval(d) for d in regdata]
    l2f = [eval(d) for d in l2fdata]
    l2e = [eval(d) for d in l2edata]
    ax1.set_xticks(reg)
    ax1.plot(reg,l2e,'-o',label=label)
    ax2.set_xticks(reg)
    ax2.plot(reg,l2f,'-o',label=label)

ax1.set_xlabel('reg')
ax1.set_ylabel('l2e')
ax1.legend()
ax1.set_yscale('log')
ax1.grid()
ax2.set_xlabel('reg')
ax2.set_ylabel('l2f')
ax2.legend()
ax2.set_yscale('log')
ax2.grid()

fig.tight_layout()
plt.savefig(f'reg_l2e_l2f_p={p}.png',dpi=200)