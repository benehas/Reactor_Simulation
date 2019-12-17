#%matplotlib qt
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import style
import pandas as pd
import random as rd
from Reaktor_Simulation import Stoffstrom,CSTRReaktion,PFRReaktion,Stoffstrom_aus_Teilströmen,RueckvermischtePFRReaktion
import scipy.optimize as spo
import copy
style.use('fivethirtyeight')

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)

data=pd.DataFrame()
ny={'A':-1,'B':-1,'C':1,'D':1}
n={'A':1,'B':1,'C':0,'D':0}

SA={'A':5.,'B':0,'C':0,'D':0}
SB={'A':0,'B':5.,'C':0,'D':0}

S2={'A':2.5,'B':2.5,'C':0,'D':0}
s2=Stoffstrom(S2,20)
CSTR=CSTRReaktion(s2,0.2,ny,n,10)
deltat=0.1
csoll=9

def opt_funkt(X):
    sA=Stoffstrom(SA,10)
    sB=Stoffstrom(SB,X)
    sMix=Stoffstrom_aus_Teilströmen(sA,sB)
    CSTRtest=CSTRReaktion(sMix,0.2,ny,n,10)
    f=open('cein.txt','r')
    CSTRtest.cein[0]=4+float(f.read())/50.
    erg=CSTRtest.berechne_erg_instat_Zeitraum(2,set_c0=True,c0=copy.deepcopy(CSTR.cplot))
    return (pow(erg[2]*CSTRtest.VReaktor/CSTRtest.tau-csoll,2))#
X0=[10]


def animate(i,cein):
    global data
    if i%5==0 :
        x=spo.fmin_l_bfgs_b(opt_funkt,X0,approx_grad=True,bounds=[(0.001,1000)],disp=0)
        sA=Stoffstrom(SA,10)
        sB=Stoffstrom(SB,x[0])
        sMix=Stoffstrom_aus_Teilströmen(sA,sB)
        CSTR.set_new_streaminput(sMix)
        print(x)
    f=open('cein.txt','r')
    data=data.append(pd.DataFrame({'t':[i],'c0':[CSTR.cplot[0]],'c1':[CSTR.cplot[1]],'c2':[CSTR.cplot[2]*CSTR.VReaktor/CSTR.tau],'c3':[CSTR.cplot[3]]}))
    CSTR.cein[0]=4+float(f.read())/50.
    c0=CSTR.berechne_erg_instat_Zeitraum(0.1,set_c0=True, c0=CSTR.cplot)
    CSTR.cplot=c0
    ys = data['c2']
    xs = data['t']
    ax1.clear()
    ax1.plot(xs, ys)
    
ani = animation.FuncAnimation(fig, animate, interval=100,fargs=[5])
plt.show()

