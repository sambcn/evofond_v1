from src.sedimentTransport.sedimentTransportLaw import SedimentTransportLaw
from src.utils import G

class PitonRecking2017(SedimentTransportLaw):

    def __init__(self):
        return

    def compute_Qs_formula(self, b, granulometry, Q, I):
        Q=max(Q,0.001)
        I=max(I,0.001)
        d84tb = granulometry.d84tb
        d84bs = granulometry.d84bs
        q=Q/b
        qe=q/(G*I*d84bs**3)**0.5
        if qe<100:
            p=0.24
        else:
            p=0.31
        taue=0.015*q**(2*p)*d84bs**(1-3*p)*I**(1-p)/(p**2.5*G**p*1.65*d84tb)
        taume=1.5*I**0.75
        phi=14*taue**2.5/(1+(taume/taue)**4)
        qsv=phi*(G*1.65*d84tb**3)**0.5
        Qsapp=qsv*b*2650/2000
        return Qsapp

    def __str__(self):
        return "Piton & Recking 2017"