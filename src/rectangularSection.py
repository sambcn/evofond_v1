from src.perf import Performance
from src.utils import G, real_roots_cubic_function
from src.irregularSection import IrregularSection
import numpy as np

class RectangularSection(IrregularSection):

    def __init__(self, x, z, b, z_min=None, y_max=None , up_section=None, down_section=None, granulometry=None, manning=0.013, K_over_tauc=None, tauc_over_rho=None):
        self.__b = b
        self.__y_max = 1000 if y_max==None or y_max <= 0 else y_max # MAX_INT
        points = [(0, self.__y_max), (0,0), (b, 0), (b, self.__y_max)]
        super().__init__(points, x, z, z_min=z_min, up_section=up_section, down_section=down_section, granulometry=granulometry, manning=manning, K_over_tauc=K_over_tauc, tauc_over_rho=tauc_over_rho)
        
    def interp_as_up_section(self, other_section, x=None):
        interpolated_section = super().interp(self, other_section, x=x)
        if x==None:
            x = 0.5*(self.get_x()+other_section.get_x())
        interpolated_section.__b = np.interp(x, [self.get_x(), other_section.get_x()], [self.get_b(), other_section.get_b()])
        return interpolated_section

    def interp_as_down_section(self, other_section, x=None):
        interpolated_section = super().interp(other_section, self, x=x)
        if x==None:
            x = 0.5*(self.get_x()+other_section.get_x())
        interpolated_section.__b = np.interp(x, [other_section.get_x(), self.get_x()], [other_section.get_b(), self.get_b()])       
        return interpolated_section

    def copy(self):
        """return a safe copy of this section"""
        return RectangularSection(self.get_x(), self.get_z(), self.get_b(), z_min=self.get_z_min(), up_section=self.get_up_section(), down_section=self.get_down_section(), granulometry=self.get_granulometry(), manning=self.get_manning())

    def update_bottom(self, Q, y_up, y, y_down, Qs_in, dt, law):
        Qs_out = law.compute_Qs(self, Q, y, y_down)
        
        # computation of maximal sediment discharge Qs_out_max
        up_section = self.get_up_section()
        down_section = self.get_down_section()
        b = self.get_b()
        b_up = up_section.get_b()
        b_down = down_section.get_b()
        dx_up = abs(up_section.get_x()-self.get_x())
        dx_down = abs(down_section.get_x()-self.get_x())
        # Numerical convention for upstream and downstream section, not very important.
        if dx_up==0:
            dx_up = dx_down
        if dx_down==0:
            dx_down = dx_up
        ######
        S =0.5*((0.75*b+0.25*b_up)*dx_up+(0.75*b+0.25*b_down)*dx_down)
        Qs_out_max = (S * (self.get_z()-self.get_z_min()) /dt) #*0.3 # z can not be less than zmin
        Qs_out = min(Qs_out, Qs_out_max)
        
        if self.is_downstream():
            Qs_out = Qs_in # max(Qs_in, Qs_out) # Sediments can not stay on the last section. 
        
        z_new = (self.get_z() + (Qs_in-Qs_out)*dt/S)

        # # check CFL
        # if not(self.is_downstream()):
        #     delta_z = abs(z_new-self.get_z())
        #     inter_S = (0.5*y+delta_z) * 0.5 * (b_up+b_down)
        #     v = Qs_out/inter_S
        #     if v*dt > dx_down:
        #         pass # LOG
        #         # print(f"WARNING : CFL CHECK FAILED IN SEDIMENT TRANSPORT (x={self.get_x()}, dt={dt}>{dx_down/v:.3f})")

        self.set_z(z_new)
        return Qs_out

    def get_stored_volume(self):
        up_section = self.get_up_section()
        down_section = self.get_down_section()
        dx_up = abs(up_section.get_x() - self.get_x())
        dx_down = abs(self.get_x() - down_section.get_x())
        if dx_up == 0:
            dx_up = dx_down
        if dx_down == 0:
            dx_down = dx_up
        dz = self.get_z() - self.get_z_min()
        b = self.get_b()
        b_up = self.get_up_section().get_b()
        b_down = self.get_down_section().get_b()
        s = 0.5*((0.75*b+0.25*b_up)*dx_up+(0.75*b+0.25*b_down)*dx_down)
        return dz*s

    # getters and setters

    @Performance.measure_perf
    def get_wet_section(self, y):
        points = self.get_points()    
        if y >= self.get_y_max():
            print("WARNING : Water depth has gone upper than the maximum one")
            return points
        else:
            return [(0, y)] + points[1:-1] + [(self.__b, y)]

    def get_b(self, y=0, wet_points=None):
        """
        Constant width.
        """
        return self.__b

    def set_b(self, b):
        if b <= 0:
            raise ValueError("width b can not be lower than 0")
        self.__b = b

    def get_S(self, y, wet_points=None):
        return y*self.__b
    
    def get_P(self, y, wet_points=None):
        return self.__b + 2*y

    def get_R(self, y, wet_points=None):
        return self.get_S(y) / self.get_P(y)

    def get_V(self, Q, y, wet_points=None):
        return Q / self.get_S(y)

    def get_H(self, Q, y, wet_points=None):
        return self.get_z() + y + (self.get_V(Q, y)**2)/(2*G)

    def get_Hs(self, Q, y, wet_points=None):
        return y + (self.get_V(Q, y)**2)/(2*G)

    # def get_Sf(self, Q, y, wet_points=None, law="Manning-Strickler"):
    #     if law=="Manning-Strickler":
    #         return (self.get_manning()*self.get_V(Q, y)/(self.get_R(y)**(2/3)))**2
    #     else:
    #         return super().get_Sf(Q, y, wet_points, law)

    def get_Fr(self, Q, y, wet_points=None):
        return ((Q**2 * self.__b)/(G * self.get_S(y)**3))**0.5
    
    def get_dP(self, y, wet_points=None):
        return 2

    def get_A_for_coussot(self, y):
        ratio = y/self.get_b(y)
        # if ratio >= 1:
        #     print("y/L >= 1...")
        return 1.93 - 0.43 * np.arctan((10*ratio)**20)
    
    def get_Fs(self, Q, y, wet_points=None):
        area = self.get_S(y)
        return 0.5*y*area + Q**2 / (area*G)

    def get_yc(self, Q):
        return (Q**2 / (G*self.__b**2))**(1/3)
    
    # @Performance.measure_perf
    # def get_y_from_Hs(self, Q, Hs, supercritical, yc=None):
    #     positive_real_roots = real_roots_cubic_function(1, -Hs, 0, Q**2 / (2*G*self.get_b()**2))
    #     if supercritical:
    #         return min(positive_real_roots)
    #     else:
    #         return max(positive_real_roots)

    # operators overloading

    def __str__(self):
        return f'Rectangular section : x={self.get_x()}, z={self.get_z()}, b={self.__b}'

    # static methods

    def is_rectangular(self):
        return True