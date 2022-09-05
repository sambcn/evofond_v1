from src.sedimentTransport.sedimentTransportLaw import SedimentTransportLaw
from src.sedimentTransport.rickenmann1991 import Rickenmann1991
from src.utils import G

class MeyerPeter1948(SedimentTransportLaw):

    def __init__(self):
        return

    def compute_Qs_formula(self, b, granulometry, Q, I):
        return Rickenmann1991().compute_Qs_formula(b, granulometry, Q, I)

    def compute_Qs(self, section, Q, y, y_down):
        Q=max(Q,0.001)
        if section.is_downstream():
            I = section.get_Sf(Q, y)
        else:
            I = (section.get_H(Q, y) - section.get_down_section().get_H(Q, y_down)) / (section.get_down_section().get_x() - section.get_x())
        I=max(I,0.001)
        # I=max(section.get_S0(),0.001)
        inter_section = section.interp_as_up_section(section.get_down_section())
        d50  = section.get_granulometry().d50
        d90  = section.get_granulometry().d90
        R = inter_section.get_R(Q, 0.5*(y+y_down))
        n_bis = (d90/26)**(1/6)
        n = (R**(2/3)*I**0.5)/inter_section.get_V(Q, 0.5*(y+y_down))
        correction_coef = (n_bis/n)**(3/2)
        correction_coef = min(1, correction_coef)

        tau_star = R*I/(1.65*d50)
        phi = 8*(correction_coef*tau_star-0.047)**(3/2)
        qsv = phi * (G*1.65*d50**3)**0.5
        Qs = qsv * inter_section.get_b(0.5*(y+y_down)) / 0.75 # debit apparent
        Qs = max(0, Qs)
        return Qs

    def __str__(self):
        return "Meyer-Peter & Muller 1948"