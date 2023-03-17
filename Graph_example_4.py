from matplotlib import pyplot as plt
import numpy as np


N = 20

L = 0.5
x = np.linspace(0, L, N)


#data1 = np.loadtxt(f"Nonlinear1_T_0.1_N_{N}.txt")
#data2 = np.loadtxt(f"Nonlinear1_T_0.105_N_{N}.txt")
data3 = np.loadtxt(f"Nonlinear1_T_0.111_N_{N}.txt")


y1 = np.zeros(N)
t1 = 0.1
for i in range(N):
    y1[i] = ((2*(0.5-x[i])**2)/(4*(0.1125-t1)))**(0.5)


y2 = np.zeros(N)
t2 = 0.105
for i in range(N):
    y2[i] = ((2*(0.5-x[i])**2)/(4*(0.1125-t2)))**(0.5)


y3 = np.zeros(N)
t3 = 0.111
for i in range(N):
    y3[i] = ((2*(0.5-x[i])**2)/(4*(0.1125-t3)))**(0.5)


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.set_xlim(0, 0.51)
ax.set_ylim(0, 10)
ax.set_xlabel(r'$x,$ м')
ax.set_ylabel(r'$u(x,t_{stop}),$ K')
ax.set_title(u'Нераспространяющийся фронт , N=50')



#ax.plot(x,y1, '-', color ='red', label =r'$t=0.1$')
#ax.plot(x, y2, '-', color='green', label=r'$t=0.105$')
ax.plot(x,y3, '-', color ='y', label =r'$t=0.111$')

#ax.plot(x, data1, '--', color='k')
#ax.plot(x, data2, '-.', color='k')
ax.plot(x, data3, ':', color='k')


ax.legend(loc='best')
fig.show()

#fig.savefig(f"Example_4_{N}.png", orientation='landscape', dpi=300)
