# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 12:15:03 2019

@author: Bened
"""
from Reaktor_Simulation import Stoffstrom,CSTRReaktion,PFRReaktion,Stoffstrom_aus_Teilströmen,RueckvermischtePFRReaktion
import numpy as np
import matplotlib.pyplot as plt
S1={'A':5.,'B':0,'C':0,'D':0}
S2={'A':5.,'B':5.,'C':0,'D':0}
S3={'A':5.,'B':5.,'C':0,'D':0}
ny={'A':-1,'B':-1,'C':1,'D':1}
n={'A':1,'B':1,'C':0,'D':0}
s1=Stoffstrom(S1,10)
s2=Stoffstrom(S2,20)
s3=Stoffstrom_aus_Teilströmen(s1,s2)

CSTR=CSTRReaktion(s2,2,ny,n,10)
print(CSTR.cein)
print(CSTR.VReaktor/CSTR.tau)
CSTR.set_new_streaminput(s3)
print(CSTR.cein)
print(CSTR.VReaktor/CSTR.tau)
PFR=PFRReaktion(s2,2,ny,n,10,10)
R2=RueckvermischtePFRReaktion(s2,2,ny,n,10,10,0)
