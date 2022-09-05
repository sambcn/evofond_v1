from src.sedimentTransport.sedimentTransportLaw import SedimentTransportLaw
from src.utils import G

class Piton2016(SedimentTransportLaw):

    def __init__(self):
        return

    def compute_Qs_formula(self, b, granulometry, Q, I):
        Q=max(Q,0.001)
        I=max(I,0.001) 
        d84bs = granulometry.d84bs
        q=Q/b
        qe=q/(G*I*d84bs**3)**0.5
        if qe<100:
            p=0.24
        else:
            p=0.31
        Qsv=0.00058*b*d84bs**(1.5-7.5*p)*q**(5*p)*I**(2.5*(1-p))/p**6.25/G**(2.5*p)
        Qsapp=Qsv*2650/2000
        return Qsapp

    def __str__(self):
        return "Piton 2016"