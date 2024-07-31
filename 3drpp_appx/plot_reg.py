import re
import matplotlib.pyplot as plt

res = open('test_noreg.log')
txt = res.read()
z0data = re.findall(r'\[x0,y0,z0\]           = \[(.*),(.*),(.*)\] for body\[(.*)\], check body\[(.*)\]\n',txt)
z0 = [eval(d[2]) for d in z0data]
data = re.findall(r'\[          \] idx=(.*),p=(.*),f=(.*)\n',txt)
potential = [eval(d[1]) for d in data]
force = [eval(d[2]) for d in data]
print(z0)
print(potential)
print(force)
# errpq = [eval(d[3]) * eval(d[5]) - eval(d[3]) * eval(d[7]) for d in data]
# en1 = sum([eval(d[3]) * eval(d[5]) for d in data])
# en2 = sum([eval(d[3]) * eval(d[7]) for d in data])
# print(f"energy1 = {en1}, energy2 = {en2}, l2e = {((en1 - en2) ** 2 / en2 / en2) ** 0.5}")

fig,ax = plt.subplots()
fig.set_figheight(6)
fig.set_figwidth(10)
ax.plot(z0,potential,'-o')
ax.set_xlabel('z0')
ax.set_ylabel('potential')
plt.savefig('z0_potential_noreg.png',dpi=200)