from src.sedimentTransport.sedimentTransportLaw import SedimentTransportLaw

class Rickenmann1991(SedimentTransportLaw):

    def __init__(self):
        return

    def compute_Qs_formula(self, b, granulometry, Q, I):
        d50 = granulometry.d50
        q=Q/b
        qcr=0.065*1.65**1.67*9.81**0.5*d50**1.5*I**(-1.12)
        qs=1.5*(q-qcr)*I**1.5
        Qs=qs*b/0.75
        Qs=max(Qs,0.)
        return Qs

    def __str__(self):
        return "Rickenmann 1991"