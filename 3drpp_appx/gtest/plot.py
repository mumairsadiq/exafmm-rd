import re
import matplotlib.pyplot as plt

res = open('/home/jiamian/rtfmm/3drpp_appx/build/gtest/debug0627.txt')
txt = res.read()
data = re.findall(r'\[(.*)\] (.*) (.*) (.*)\((.*)\) (.*) \((.*)\) (.*) \((.*)\) (.*) \((.*)\) (.*)\n',txt)
errpq = [eval(d[3]) * eval(d[5]) - eval(d[3]) * eval(d[7]) for d in data]
en1 = sum([eval(d[3]) * eval(d[5]) for d in data])
en2 = sum([eval(d[3]) * eval(d[7]) for d in data])
print(f"energy1 = {en1}, energy2 = {en2}, l2e = {((en1 - en2) ** 2 / en2 / en2) ** 0.5}")

fig,ax = plt.subplots()
fig.set_figheight(6)
fig.set_figwidth(10)
ax.plot(range(len(errpq)),errpq)
ax.set_xlabel('bidx')
ax.set_ylabel('errpq')
plt.savefig('errpq.png',dpi=200)