from src.utils import G
from src.irregularSection import IrregularSection
import numpy as np

class TrapezoidalSection(IrregularSection):

    def __init__(self, x, z, b, s, z_min=None, y_max=None , up_section=None, down_section=None, granulometry=None, manning=0.013, K_over_tauc=None, tauc_over_rho=None):
        self.__b = b
        self.__s = s
        self.__y_max = 1000 if y_max==None or y_max <= 0 else y_max # MAX_INT
        points = [(0, self.__y_max), (s*self.__y_max,0), (s*self.__y_max+b, 0), (2*s*self.__y_max+b, self.__y_max)]
        super().__init__(points, x, z, z_min=z_min, up_section=up_section, down_section=down_section, granulometry=granulometry, manning=manning, K_over_tauc=K_over_tauc, tauc_over_rho=tauc_over_rho)
        
    def interp_as_up_section(self, other_section, x=None):
        interpolated_section = super().interp(self, other_section, x=x)
        if x==None:
            x = 0.5*(self.get_x()+other_section.get_x())
        interpolated_section.__b = np.interp(x, [self.get_x(), other_section.get_x()], [self.get_b(), other_section.get_b()])
        interpolated_section.__s = np.interp(x, [self.get_x(), other_section.get_x()], [self.get_s(), other_section.get_s()])
        return interpolated_section

    def interp_as_down_section(self, other_section, x=None):
        interpolated_section = super().interp(other_section, self, x=x)
        if x==None:
            x = 0.5*(self.get_x()+other_section.get_x())
        interpolated_section.__b = np.interp(x, [other_section.get_x(), self.get_x()], [other_section.get_b(), self.get_b()])       
        interpolated_section.__s = np.interp(x, [other_section.get_x(), self.get_x()], [other_section.get_s(), self.get_s()])       
        return interpolated_section

    def copy(self):
        """return a safe copy of this section"""
        return TrapezoidalSection(self.get_x(), self.get_z(), self.get_b(), self.get_s(), z_min=self.get_z_min(), y_max=self.get_y_max(), up_section=self.get_up_section(), down_section=self.get_down_section(), granulometry=self.get_granulometry(), manning=self.get_manning())

    def update_bottom(self, Q, y, y_next, Qs_in, dt, law):
        raise(NotImplementedError("solid transport is not well defined for trapezoidal sections."))
        #TODO

    def get_stored_volume(self):
        raise(NotImplementedError("solid transport is not well defined for trapezoidal sections."))
        #TODO

    # getters and setters

    def get_wet_section(self, y):
        points = self.get_points()    
        if y >= self.get_y_max():
            print("WARNING : Water depth has gone upper than the maximum one")
            return points
        else:
            return [(self.__s*self.__y_max - self.__s*y, y)] + points[1:-1] + [(self.__s*self.__y_max+self.__b+y*self.__s, y)]

    def get_b(self, y=0, wet_points=None):
        return self.__b + 2*self.__s*y

    def set_b(self, b):
        if b <= 0:
            raise ValueError("width b can not be lower than 0")
        self.__b = b

    def get_S(self, y, wet_points=None):
        return (self.__b+self.__s*y)*y
    
    def get_P(self, y, wet_points=None):
        return self.__b + 2*y*np.sqrt(1+self.__s**2)

    def get_R(self, y, wet_points=None):
        return self.get_S(y) / self.get_P(y)

    def get_V(self, Q, y, wet_points=None):
        return Q / self.get_S(y)

    def get_H(self, Q, y, wet_points=None):
        return self.get_z() + y + (self.get_V(Q, y)**2)/(2*G)

    def get_Hs(self, Q, y, wet_points=None):
        return y + (self.get_V(Q, y)**2)/(2*G)

    def get_Sf(self, Q, y, wet_points=None):
        return (self.get_manning()*self.get_V(Q, y)/(self.get_R(y)**(2/3)))**2

    def get_Fr(self, Q, y, wet_points=None):
        return self.get_V(Q, y)/((G * y)**0.5)
    
    def get_dP(self, y, wet_points=None):
        return 2*np.sqrt(A+self.__s**2)

    def __get_A_for_coussot(self, y):
        return 1.93 - 0.6 * np.arctan((0.4*y/self.get_b(0))**20)
    
    def get_s(self):
        return self.__s

    # Operators overloading

    def __str__(self):
        return f'Trapezoidal section : x={self.get_x()}, z={self.get_z()}, b={self.__b}, s={self.__s}'