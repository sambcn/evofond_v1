from abc import ABC, abstractmethod
from src.perf import Performance

class SedimentTransportLaw(ABC):
    """
    This class is like an interface of Object-Oriented programming.
    It is only used to hide all the transport law under an unique class.
    It allows to avoid the 
    It could be done better with the right modules, but in the first place it will work well.
    """
     
    # computation functions

    @Performance.measure_perf
    def compute_Qs(self, section, Q, y, y_down): 
        Q=max(Q,0.001)
        if section.is_downstream():
            I = section.get_Sf(Q, y)
        else:
            I = (section.get_H(Q, y) - section.get_down_section().get_H(Q, y_down)) / (section.get_down_section().get_x() - section.get_x())
        I=max(I,0.001)
        # I=max(section.get_S0(),0.001)
        b = 0.5*(section.get_b(y) + section.get_down_section().get_b(y_down))
        return self.compute_Qs_formula(b, section.get_granulometry(), Q, I)

    @abstractmethod
    def compute_Qs_formula(self, b, granulometry, Q, I): 
        pass

    @abstractmethod
    def compute_Qs_formula(self, b, granulometry, Q, I): 
        pass

    # operators overloading

    def __str__(self):
        return "unknown SedimentTransportLaw"
