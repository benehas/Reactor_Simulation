# -*- coding: utf-8 -*-
"""
Created on Sun Dec  1 15:39:59 2019

@author: Bened
"""
from scipy.optimize import fsolve
import math
import numpy as np
import copy


class Stoffstrom():
    def __init__(self,Stoffe,Volumenstrom):
        self.stoffe=Stoffe# Dictionariy mit Namen : Konzentration
        self.Vp=Volumenstrom
        

class Stoffstrom_aus_Teilströmen(Stoffstrom):
    def __init__(self,Teilstrom1,Teilstrom2):
        Vpges=Teilstrom1.Vp+Teilstrom2.Vp
        Stoffe=copy.deepcopy(Teilstrom1.stoffe)
        for i in Stoffe.keys():
            Stoffe[i]=Stoffe[i]*Teilstrom1.Vp/Vpges
            if i in Teilstrom2.stoffe.keys():
                Stoffe[i]=Stoffe[i]+Teilstrom2.stoffe[i]*Teilstrom2.Vp/Vpges
        for i in Teilstrom2.stoffe.keys():
            if i not in Stoffe.keys():
                Stoffe[i]=Teilstrom2.stoffe[i]*Teilstrom2.Vp/Vpges
        Stoffstrom(Stoffe,Vpges)
        self.stoffe=Stoffe# Dictionariy mit Namen : Konzentration
        self.Vp=Vpges
        
class CSTRReaktion:
    def __init__(self,Stoffstromein,k,ny,n,VReaktor):
        self.SEin=dict()
        for i in Stoffstromein.stoffe.keys():
            #print(i)
            #self.SEin[i]=[Stoffstromein.stoffe[i],ny[i],n[i]]
            self.SEin.update({i:[Stoffstromein.stoffe[i],ny[i],n[i]]})
        self.k=k #Kinetishcer Faktor
        self.ny=np.array(list(self.SEin.values()))[:,1].tolist() #Stöchometriesche Koeffizienten
        self.n=np.array(list(self.SEin.values()))[:,2].tolist() #Kinetische Teilordnung
        self.cein=np.array(list(self.SEin.values()))[:,0].tolist()
        self.tau=VReaktor/Stoffstromein.Vp
        self.VReaktor=VReaktor
        self.cplot=np.zeros(len(self.cein))

    def set_new_streaminput(self,Stoffstromein):
        self.SEin=dict()
        for i in Stoffstromein.stoffe.keys():
            self.SEin.update({i:[Stoffstromein.stoffe[i]]})
        self.cein=np.array(list(self.SEin.values()))[:,0].tolist()
        self.tau=self.VReaktor/Stoffstromein.Vp

        
    #Berechnet stationäre Lösung
    def berechne_erg_stat(self):
        X0=np.zeros(len(self.cein))
        X=fsolve(self.eq,X0)
        #self.SAus=copy.deepcopy(self.SEin)
        self.caus=X
        return X
    
    #Berechnet einen instationären Zeitschritt explizit über Gauß Vorwärts
    def berechne_erg_instat(self,c,deltat=0.01):
        for index,i in enumerate(c):
            tempr=self.k*self.ny[index]
            for j in range(len(c)):
                tempr=tempr*math.pow(c[j],self.n[j])
            c[index]=i+deltat*((self.cein[index]-i)/self.tau +tempr)
        return c
    
    #Berechnet instationär über Zeitraum
    def berechne_erg_instat_Zeitraum(self,maxT,deltat=0.01,AlleWerte=False,set_c0=False,c0=0):
        if set_c0:
            pass
        else:
            c0=np.zeros(len(self.cein))
        tmp=self.berechne_erg_instat(c0)
        c=[]
        maxN=int(maxT/deltat)
        a=np.zeros((maxN,len(self.cein)))
        for i in range(maxN):
            tmp=self.berechne_erg_instat(tmp)
            c.append(tmp)
            a[i]=tmp
        if AlleWerte==True:
            return a
        else:
            return a[maxN-1]
        
    def eq(self,X):
        temp=list()
        for i in range(len(X)):
            tempr=self.k*self.tau*self.ny[i]
            for j in range(len(X)):
                tempr=tempr*math.pow(X[j],self.n[j])
            temp.append(X[i]-tempr-self.cein[i])
        return temp
    
    def eq_imp(self,X,deltat=0.01):
        tmp=list()
        tempr=self.k
        for j in range(len(X)):
            tempr=tempr*math.pow(X[j],self.n[j])        
        for j in range(len(X)):
            tmp.append(1/self.tau*(self.cein[j]-X[j])+self.ny[j]*tempr+1/deltat*(self.c[j]-X[j]))
        return tmp

    #Berechnet instationär über Zeitraum implizit
    def berechne_erg_instat_Zeitraum_imp(self,maxT,deltat=0.01,AlleWerte=False,set_c0=False,c0=0):
        if set_c0:
            pass
        else:
            c0=np.zeros(len(self.cein))
        self.c=c0
        X0=np.ones((len(self.cein)))
        maxN=int(maxT/deltat)
        a=np.zeros((maxN,len(self.cein)))
        for i in range(maxN):
            self.c=fsolve(self.eq_imp,X0,deltat)
            a[i]=self.c
        if AlleWerte==True:
            return a
        else:
            return a[maxN-1]
            


    

class PFRReaktion:
    def __init__(self,Stoffstromein,k,ny,n,VReaktor,N):
        self.Stoffstromein=Stoffstromein
        self.SEin=dict()
        for i in Stoffstromein.stoffe.keys():
            self.SEin[i]=[Stoffstromein.stoffe[i],ny[i],n[i]]
        self.k=k #Kinetishcer Faktor
        self.ny=np.array(list(self.SEin.values()))[:,1].tolist() #Stöchometriesche Koeffizienten
        self.n=np.array(list(self.SEin.values()))[:,2].tolist() #Kinetische Teilordnung
        self.cein=np.array(list(self.SEin.values()))[:,0].tolist()
        self.N=N #Anzahl der CSTRs
        self.DeltaV=VReaktor/N
        self.tau=self.DeltaV/Stoffstromein.Vp #Raumzeit der einzelnen CSTRs
        self.ny_stat=ny
        self.n_stat=n
        self.cplot=np.zeros(len(self.cein))

    def set_new_streaminput(self,Stoffstromein):
        self.SEin=dict()
        for i in Stoffstromein.stoffe.keys():
            self.SEin.update({i:[Stoffstromein.stoffe[i]]})
        self.cein=np.array(list(self.SEin.values()))[:,0].tolist()
        self.tau=self.VReaktor/Stoffstromein.Vp


    def berechne_erg_stat(self):
        c=self.cein
        for i in range(1,self.N):
            c=self.berechne_cstr(c)            
        return np.array(c)

    def berechne_erg_instat(self,c,deltat=0.01):
        cp1=copy.deepcopy(c)
        for i in range(len(self.cein)):
            tempr=self.k*self.ny[i]
            for j in range(len(self.cein)):
                tempr=tempr*math.pow(c[0][j],self.n[j])
            cp1[0][i]=deltat*(1/self.tau*(self.cein[i]-c[0][i])+tempr)+c[0][i]
        for g in range(1,len(c)):
            for i in range(len(self.cein)):
                tempr=self.k*self.ny[i]
                for j in range(len(self.cein)):
                    tempr=tempr*math.pow(c[g][j],self.n[j])
                cp1[g][i]=deltat*(1/self.tau*(c[g-1][i]-c[g][i])+tempr)+c[g][i]
        return cp1

    
       #Berechnet instationär über Zeitraum
    def berechne_erg_instat_Zeitraum(self,maxT,deltat=0.001,AlleWerte=False):
        k=list()
        returnvalue=list()
        k=np.zeros((self.N,len(self.cein)))
        k=k.tolist()
        maxN=int(maxT/deltat)
        for i in range(maxN):
            k=self.berechne_erg_instat(k,deltat)
            self.ci=np.array(k)
            returnvalue.append(k[self.N-1])
        if AlleWerte==True:
            return np.array(returnvalue)
        else:
            return np.array(returnvalue[maxN-1])
    
    def eq(self,X,cein):
        temp=list()
        for i in range(len(X)):
            tempr=self.k*self.tau*self.ny[i]
            for j in range(len(X)):
                tempr=tempr*math.pow(X[j],self.n[j])
            temp.append(X[i]-tempr-cein[i])
        return temp

    def berechne_cstr(self,cein):
        X0=np.zeros(len(self.cein))
        X=fsolve(self.eq,X0,cein)
        return X 

    def eq2(self,cip1,deltat):
        ci=copy.deepcopy(self.ci)
        ci=np.array(ci).reshape((self.N,len(self.cein)))
        deltat=deltat
        eqlist=list()
        tempr=self.k
        cip1=cip1.reshape((self.N,len(self.cein)))
        for j in range(len(self.cein)):
            tempr=tempr*math.pow(cip1[0,j],self.n[j])
        for j in range(len(self.cein)):
            eqlist.append((self.cein[j]-cip1[0,j])/self.tau
                          +self.ny[j]*tempr
                          -(cip1[0,j]-ci[0,j])/deltat)
        for i in range(1,self.N):
            tempr=self.k
            for j in range(len(self.cein)):
                tempr=tempr*math.pow(cip1[i,j],self.n[j])
            for j in range(len(self.cein)):
                eqlist.append((cip1[i-1,j]-cip1[i,j])/self.tau+self.ny[j]*tempr-(cip1[i,j]-ci[i,j])/deltat)
        return eqlist

    def berechne_erg_instat_Zeitraum_imp(self,maxT,deltat=0.001,AlleWerte=False):
        ret=list()
        self.ci=np.zeros((self.N*len(self.cein)))
        X0=np.ones((self.N*len(self.cein)))
        maxN=int(maxT/deltat)
        for i in range(maxN):
            self.ci=np.array(fsolve(self.eq2,X0,deltat))
            self.ci=self.ci.reshape((self.N,len(self.cein)))
            ret.append(self.ci[self.N-1])
        if AlleWerte==True:
            return np.array(ret)
        else:
            return np.array(ret[maxN-1])               

class RueckvermischtePFRReaktion:
    def __init__(self,Stoffstromein,k,ny,n,VReaktor,N,ruecklauf):
        self.Stoffstromein=Stoffstromein
        self.SEin=dict()
        for i in Stoffstromein.stoffe.keys():
            self.SEin[i]=[Stoffstromein.stoffe[i],ny[i],n[i]]
        self.k=k #Kinetishcer Faktor
        self.ny=np.array(list(self.SEin.values()))[:,1].tolist() #Stöchometriesche Koeffizienten
        self.n=np.array(list(self.SEin.values()))[:,2].tolist() #Kinetische Teilordnung
        self.cein=np.array(list(self.SEin.values()))[:,0].tolist()
        self.N=N #Anzahl der CSTRs
        self.DeltaV=VReaktor/N
        self.tau=self.DeltaV/(Stoffstromein.Vp*(ruecklauf+1)) #Raumzeit der einzelnen CSTRs
        self.ny_stat=ny
        self.n_stat=n
        self.r=ruecklauf

    def set_new_streaminput(self,Stoffstromein):
        self.SEin=dict()
        for i in Stoffstromein.stoffe.keys():
            self.SEin.update({i:[Stoffstromein.stoffe[i]]})
        self.cein=np.array(list(self.SEin.values()))[:,0].tolist()
        self.tau=self.VReaktor/Stoffstromein.Vp


    def berechne_erg_stat(self):
        c=self.cein
        caus=np.zeros(len(c))
        causm1=np.zeros(len(c))
        for j in range(1000):
            c=(self.cein+caus*self.r)/(self.r+1)
            for i in range(1,self.N):
                c=self.berechne_cstr(c)
            caus=c
            if abs(causm1.sum()-caus.sum())<0.0001:
                break
            causm1=caus
        return np.array(c)

    def berechne_erg_instat(self,c,deltat=0.001):
        cp1=copy.deepcopy(c)
        for i in range(len(self.cein)):
            tempr=self.k
            for j in range(len(self.cein)):
                tempr=tempr*math.pow(c[0][j],self.n[j])
            cp1[0][i]=deltat*(1/self.tau*((self.cein[i]+self.r*c[-1][i])/(self.r+1)-c[0][i])+self.ny[i]*tempr)+c[0][i]
        for g in range(1,len(c)):
            for i in range(len(self.cein)):
                tempr=self.k
                for j in range(len(self.cein)):
                    tempr=tempr*math.pow(c[g][j],self.n[j])
                cp1[g][i]=deltat*(1/self.tau*(c[g-1][i]-c[g][i])+self.ny[i]*tempr)+c[g][i]        
        return cp1
    
       #Berechnet instationär über Zeitraum
       #instabil für zu große N und R
    def berechne_erg_instat_Zeitraum(self,maxT,deltat=0.001,AlleWerte=False):
        k=list()
        ret=list()
        k=np.zeros((self.N,len(self.cein)))
        k=k.tolist()
        maxN=int(maxT/deltat)
        for i in range(maxN):
            k=self.berechne_erg_instat(k,deltat=deltat)
            ret.append(k[self.N-1])
        if AlleWerte==True:
            return ret
        else:
            return np.array(ret[maxN-1])

    def eq(self,X,cein):
        temp=list()
        for i in range(len(X)):
            tempr=self.k*self.tau*self.ny[i]
            for j in range(len(X)):
                tempr=tempr*math.pow(X[j],self.n[j])
            temp.append(X[i]-tempr-cein[i])
        return temp

    def berechne_cstr(self,cein):
        X0=np.zeros(len(self.cein))
        X=fsolve(self.eq,X0,cein)
        return X                
    
    def eq2(self,cip1,deltat):
        ci=copy.deepcopy(self.ci)
        ci=np.array(ci).reshape((self.N,len(self.cein)))
        deltat=deltat
        eqlist=list()
        tempr=self.k
        cip1=cip1.reshape((self.N,len(self.cein)))
        for j in range(len(self.cein)):
            tempr=tempr*math.pow(cip1[0,j],self.n[j])
        for j in range(len(self.cein)):
            eqlist.append(((self.cein[j]+self.r*cip1[-1,j])/(1+self.r)-cip1[0,j])/self.tau
                          +self.ny[j]*tempr
                          -(cip1[0,j]-ci[0,j])/deltat)
        for i in range(1,self.N):
            tempr=self.k
            for j in range(len(self.cein)):
                tempr=tempr*math.pow(cip1[i,j],self.n[j])
            for j in range(len(self.cein)):
                eqlist.append((cip1[i-1,j]-cip1[i,j])/self.tau+self.ny[j]*tempr-(cip1[i,j]-ci[i,j])/deltat)
        return eqlist

    def berechne_erg_instat_Zeitraum_imp(self,maxT,deltat=0.001,AlleWerte=False):
        ret=list()
        self.ci=np.zeros((self.N*len(self.cein)))
        X0=np.ones((self.N*len(self.cein)))
        maxN=int(maxT/deltat)
        for i in range(maxN):
            self.ci=np.array(fsolve(self.eq2,X0,deltat))
            self.ci=self.ci.reshape((self.N,len(self.cein)))
            ret.append(self.ci[self.N-1])
        if AlleWerte==True:
            return np.array(ret)
        else:
            return np.array(ret[maxN-1])