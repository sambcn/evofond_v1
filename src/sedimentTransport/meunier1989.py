from src.sedimentTransport.sedimentTransportLaw import SedimentTransportLaw

class Meunier1989(SedimentTransportLaw):

    def __init__(self):
        return

    def compute_Qs_formula(self, b, granulometry, Q, I):
        return 8.2*Q*I**2

    def __str__(self):
        return "Meunier 1989"