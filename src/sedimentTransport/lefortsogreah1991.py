from src.sedimentTransport.sedimentTransportLaw import SedimentTransportLaw
import math

class LefortSogreah1991(SedimentTransportLaw):

    def __init__(self):
        return

    def compute_Qs_formula(self, b, granulometry, Q, I):
        dm = granulometry.dm
        d90 = granulometry.d90
        d30 = granulometry.d30

        I=min(0.83,max(I,0.001))
        Qcrit=0.0776*(9.81*dm**5)**0.5*1.65**(8./3.)/I**(13./6.)*(1-1.2*I)**(8./3.)
        Qs=4.45*Q*math.pow(I,1.5)/1.65*math.pow(d90/d30,0.2)*(1-math.pow((Qcrit/Q),0.375))
        return max(Qs, 0.01) # return Qs a la base


    def __str__(self):
        return "Lefort-Sogreah 1991"