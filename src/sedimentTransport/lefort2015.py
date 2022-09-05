from src.sedimentTransport.sedimentTransportLaw import SedimentTransportLaw
import math

class Lefort2015(SedimentTransportLaw):

    def __init__(self):
        return

    def compute_Qs_formula(self, b, granulometry, Q, I):
        dm = granulometry.dm
        Gr = granulometry.Gr
        dms=dm*(9.81*1.65/0.000001**2)**(1./3.)
        Cdms=0.0444*(1+15/(1+dms)-1.5*math.exp(-dms/75))
        qs=(Q/b)/(9.81*I*dm**3)**0.5
        if qs<200 :
            rkskr=0.75*(qs/200)**0.23
        else :
            rkskr=0.75
        n=1.6+0.06*math.log10(I)    
        m=1.8+0.08*math.log10(I)
        q0=(9.81*(1.65*dm)**3)**0.5*Cdms*(dm/b)**(1./3.)*rkskr**-0.5*I**-n
        if (dms<14 and rkskr<0.63) :
            cor=1-1.4*math.exp(-0.9*rkskr**2*(Q/b/q0)**0.5)
        else :
            cor=1.
        M=(qs+2.5)/200.
        Z=1+0.38/dms**0.45*(Q/b/(9.81*dm**3)**0.5)**0.192

        if (Q/b)<q0:
            F=0.06*M*(Q/b)/q0
        else:
            F=(6.1*(1-0.938*(q0*b/Q)**0.284)**1.66)**Z
        
        Cp=1700000*I**m*2.65/1.65**1.65*Gr**0.2*cor*F
        Qs=Cp/1000.*Q/2650*(2650./2000.)
        return Qs

    def __str__(self):
        return "Lefort 2015"