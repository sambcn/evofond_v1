import matplotlib.pyplot as plt
import numpy as np

from copy import copy
from scipy.misc import derivative
from scipy.optimize import brentq, newton
from src.perf import Performance
from src.utils import G, Y_MIN, get_centroid

class IrregularSection:

    def __init__(self, points, x, z, z_min=None, up_section=None, down_section=None, granulometry=None, manning=None, tauc_over_rho=None, K_over_tauc=None):
        """
        Constructor and initializations
        Args :
            points (list of tuples) : points which describe the cross section (for example [(0, 10), (1, 0), (10, 0), (10, 10)])
            z (float) : total height of the lowest point of the section
            z_min (float) : height minimal of the lowest point of the section
        """

        # initializations
        if len(points) < 3:
            raise(ValueError("Error : you need at least 3 points to describe a section."))
        self.__initialPoints = copy(points)                             # points of original section are saved in this hidden variable
        self.__points = copy(points)                                    # import of points of cross section
        self.setup_points()                                             # initialize self.__x_list, self.__y_list and setup self.__points to verify some conditions (more details in the methods)
        self.__index_min = self.__y_list.index(min(self.__y_list))      # index of the lowest point of the section
        self.__x_min = self.__x_list[self.__index_min]
        self.__y_min = self.__y_list[self.__index_min]
        self.__x = x                                                    # abscissa of the section
        self.__z = z
        if z_min==None or z_min>z:
            self.__z_min = z
        else:
            self.__z_min = z_min
        
        if up_section==None:
            self.__up_section = self
            self.__is_upstream = True
        else:
            self.__up_section = up_section
            self.__is_upstream = False

        if down_section==None:
            self.__down_section = self
            self.__is_downstream = True
        else:
            self.__down_section = down_section
            self.__is_downstream = False

        self.__granulometry = granulometry
        if manning != None:
            self.__manning = manning 
        elif granulometry.d84bs != None:
            self.__manning = (granulometry.d84bs)**(1/6)/26
        else:
            self.__manning = None
            print("WARNING : no value for manning coefficient : you did neither precise the value of it, nor d84 which is used for the default assignment. This might produce an error if you want to use Manning-Strickler friction law")

        self.__tauc_over_rho = tauc_over_rho
        self.__K_over_tauc = K_over_tauc


    def setup_points(self):
        """
        Change the attribute self.__points in order to make it satisfy the following conditions :
        - sorted by x increasing (1)
        - x[0] = 0 and y_min = 0 (normalization) (2)
        - the first and the last points have the same y (3)
        - there is no y greater than the first and the last one (4)

        it also updates the attributes self.__x_list and self.__y_list
        """

        self.__points.sort(key = lambda p : p[0]) # (1)
        # (2)
        x = [p[0] for p in self.__points]
        x = list(np.array(x)-x[0])
        y = [p[1] for p in self.__points]
        y = list(np.array(y)-min(y))

        # (3)
        if y[0] > y[-1]:
            y.append(y[0])
            x.append(x[-1])
        elif y[-1] > y[0]:
            y.insert(0, y[-1])
            x.insert(0, x[0])

        # (4)
        if y[0] < max(y):
            y.insert(0, max(y))
            x.insert(0, x[0])
        if y[-1] < max(y):
            y.append(max(y))
            x.append(x[-1])
        
        self.__x_list = x
        self.__y_list = y
        self.__y_max = max(y)
        return

    def interp_as_up_section(self, other_section, x=None):
        interpolated_section = IrregularSection.interp(self, other_section, x=x)
        return interpolated_section

    def interp_as_down_section(self, other_section, x=None):
        interpolated_section = IrregularSection.interp(other_section, self, x=x)
        return interpolated_section

    def copy(self):
        """return a safe copy of this section"""
        return IrregularSection(self.__points[:], self.get_x(), self.get_z(), z_min=self.get_z_min(), up_section=self.get_up_section(), down_section=self.get_down_section(), granulometry=self.get_granulometry(), manning=self.get_manning(), K_over_tauc=self.get_K_over_tauc(), tauc_over_rho=self.get_tauc_over_rho())

    def get_stored_volume(self):  
        print("get_stored_volume not defined yet for an irregular section. Return 0.")
        #TODO
        return 0

    def get_wet_section(self, y):
        """
        Args :
            y (float) : depth of water (m)
        Returns :
            wet_points (list of tuples) : list of points that describe wet section
        """
        points = self.__points    
        if y >= self.get_y_max():
            print(f"WARNING : Water depth has gone upper than the maximum one ({y}>{self.__y_max}m)")
            return points        
        # Look for the left point of intersection
        index=self.__index_min
        left=True
        while left:
            index+=-1
            if self.__y_list[index] >= y:
                left = False
                x1,y1 = points[index+1]
                x2,y2 = points[index]
                x_left = (y - y1) * (x2 - x1) / (y2 - y1) + x1 # Thales used to know abscissa of water surface left limit
                index_left=index
                
        # Look for the right point of intersection
        index=self.__index_min
        right=True
        while right:
            index+=1
            if self.__y_list[index] > y:
                right = False
                x1,y1 = points[index-1]
                x2,y2 = points[index]
                x_right = (y - y1) * (x2 - x1) / (y2 - y1) + x1
                index_right=index
            
        wet_points=[(x_left, y)] + self.__points[index_left+1:index_right] + [(x_right, y)]        
        return wet_points

    # getter and setter which require computations and which can be simplified in children classes

    def get_b(self, y, wet_points=None):
        """width"""
        wet_points = self.get_wet_section(y) if wet_points==None else wet_points
        return abs(wet_points[-1][0]-wet_points[0][0])

    def get_S(self, y, wet_points=None):
        """wet surface"""
        wet_points = self.get_wet_section(y) if wet_points==None else wet_points
        s = 0
        for i in range(len(wet_points)-1):
            p_current = wet_points[i]
            p_next = wet_points[i+1]
            dx = abs(p_next[0]-p_current[0])
            dy1 = abs(y-p_current[1])
            dy2 = abs(y-p_next[1])
            s += (dy1 + dy2)*dx*0.5
        return s    

    def get_P(self, y, wet_points=None):
        """wet perimeter"""
        wet_points = self.get_wet_section(y) if wet_points==None else wet_points
        p = 0
        for i in range(len(wet_points)-1):
            p_current = wet_points[i]
            p_next = wet_points[i+1]
            dx = p_next[0]-p_current[0]
            dy= p_next[1]-p_current[1]
            p += np.sqrt(dx**2 + dy**2)
        return p
    
    def get_R(self, y, wet_points=None):
        """hydraulic radius"""
        wet_points = self.get_wet_section(y) if wet_points==None else wet_points
        return self.get_S(y, wet_points) / self.get_P(y, wet_points)

    def get_V(self, Q, y, wet_points=None):
        """mean flow velocity"""
        return Q / self.get_S(y, wet_points)

    def get_H(self, Q, y, wet_points=None):
        """total head"""
        return self.get_z() + y + (self.get_V(Q, y, wet_points)**2)/(2*G)

    def get_Hs(self, Q, y, wet_points=None):
        """specific head"""
        return y + (self.get_V(Q, y, wet_points)**2)/(2*G)

    def get_Sf(self, Q, y, wet_points=None, friction_law="Ferguson"):
        """slope of energy grade line, computed with a friction law"""
        wet_points = self.get_wet_section(y) if wet_points==None else wet_points
        if friction_law=="Coussot":
            # check that K/tauc, taux/ro are defined... 
            sinus = (((self.get_V(Q, y, wet_points=wet_points)/y)**(3/10))*self.get_A_for_coussot(y)*(self.get_K_over_tauc()**(9/10))+1)*self.get_tauc_over_rho()/(G*self.get_R(y, wet_points))
            if sinus < 1 and sinus > -1:
                return np.sqrt(sinus**2/(1-sinus**2))
            else:
                return self.get_Sf(Q, y, wet_points)
            # print(f"x={self.get_x()} : {sinus}")
            # i =  np.tan(np.arcsin(sinus))
            # i =  np.sqrt(sinus**2/(1-sinus**2))
        else:
            return (self.get_Cf(Q, y, wet_points, friction_law=friction_law) * self.get_V(Q, y, wet_points)**2) / (G*self.get_R(y, wet_points)) #(self.get_manning()*self.get_V(Q, y, wet_points)/(self.get_R(y, wet_points)**(2/3)))**2
    
    def get_Cf(self, Q, y, wet_points=None, friction_law="Ferguson"):
        """dimensionless friction coefficient"""
        wet_points = self.get_wet_section(y) if wet_points==None else wet_points
        if friction_law=="Manning-Strickler":
            # check that manning coef is defined...
            return (G*self.get_manning()**2) / (self.get_R(y, wet_points)**(1/3)) # manning
        elif friction_law=="Ferguson":
            # check that granulometry is complete....
            coef = self.get_R(y, wet_points)/self.get_granulometry().d84bs
            return  (1+0.15*coef**(5/3)) / ((2.5*coef)**2) # Ferguson
        else:
            raise ValueError(f"unknown friction law : {friction_law}")

    def get_A_for_coussot(self, y):
        print(self)
        raise ValueError("you cannot compute geometrical coefficient of coussot formula for an irregular section")

    def get_Fr(self, Q, y, wet_points=None):
        """Froude number"""
        wet_points = self.get_wet_section(y) if wet_points==None else wet_points
        return ((Q**2 * self.get_b(y, wet_points))/(G * self.get_S(y, wet_points)**3))**0.5

    def get_Fs(self, Q, y, wet_points=None):
        """specific force"""      
        wet_points = self.get_wet_section(y) if wet_points==None else wet_points
        centroid = get_centroid(wet_points)
        centroid_depth = y-centroid[1]
        area = self.get_S(y, wet_points=wet_points)
        return centroid_depth*area + Q**2 / (area*G)

    def get_yc(self, Q):
        """critical water depth"""
        def equation_function(y): # Froude number is one for critical depth
            wet_points = self.get_wet_section(y)
            Fr = self.get_Fr(Q, y, wet_points=wet_points)
            return 1-Fr
        try:
            yc = brentq(equation_function, Y_MIN, 0.999*self.get_y_max())
        except ValueError:  # f(a) and f(b) must have different signs
            try:
                yc = newton(equation_function, Y_MIN, tol=0.001, maxiter=1000) # Utilisation de la méthode de la sécante
                print("Warning : brentq method failed to converge for yc, used Newton")
            except RuntimeError as e:  # Failed to converge
                print(f"WARNING : Newton method failed to converge.")
                raise(e)
        yc = max(yc, 0)
        return yc

    def get_yn(self, Q, friction_law="Ferguson"):
        """normal water depth"""
        b = self.is_downstream()
        s0 = self.get_S0(up_direction=b)
        if s0 < 0:
            return 0
        elif s0 == 0:
            return self.__y_max
        def equation_function(y): # Friction law for a uniform regim leads to this equation
            wet_points = self.get_wet_section(y)
            sf = self.get_Sf(Q, y, wet_points=wet_points, friction_law=friction_law)
            s0 = self.get_S0(up_direction=b)
            return sf-s0
        try:
            yn = brentq(equation_function, Y_MIN, 0.999*self.get_y_max())
        except ValueError:  # f(a) and f(b) must have different signs
            try:
                yn = newton(equation_function, Y_MIN, tol=0.001, maxiter=1000)
                print("Warning : brentq method failed to converge for yn, used Newton")
            except RuntimeError as e:  # Failed to converge
                print(f"WARNING : Newton method failed to converge.")
                raise(e)
        yn = max(yn, 0)
        return yn

    @Performance.measure_perf
    def get_y_from_Hs(self, Q, Hs, supercritical, yc=None):
        """
        Takes a given water discharge Q and a specific head Hs and compute the associate water depth y.
        It returns the supercritical solution if supercritical is True, subcritical else.
        yc is given to know the range of solution ([0,yc] for supercritical and [yc,+inf] for subcritical).
        """
        yc = yc if yc != None else self.get_yc(Q)
        def objective_function(y):
            return Hs - self.get_Hs(Q, y)
        try:
            if supercritical:
                return brentq(objective_function, Y_MIN, yc)
            else:
                return brentq(objective_function, yc, 0.999*self.get_y_max())
        except ValueError as e:
            y_solution = "Y_MIN" if supercritical else "Y_MAX"
            print(f"BIG WARNING : in computing y from Hs, could not find a fine solution, returned {y_solution} (x={self.get_x()})")
            return Y_MIN if supercritical else 0.999*self.get_y_max()

    # getters and setters that will not change for children classes

    def get_x(self):
        return self.__x

    def get_z(self):
        return self.__z
        
    def get_z_min(self):
        return self.__z_min

    def set_x(self, x):
        if x < self.__up_section.__x:
            raise ValueError("the abscissa of a section can not be lower than its up_section")
        if x > self.__down_section.__x:
            raise ValueError("the abscissa of a section can not be greater than its down_section")
        self.__x = x

    def set_z(self, z):
        if z < self.__z_min:
            raise ValueError("z can not be lower than z_min.")
        self.__z = z
    
    def set_z_min(self, z_min):
        if self.__z < z_min:
            raise ValueError("z_min can not be greater than z")
        self.__z_min = z_min

    def get_y_max(self):
        return max(self.__y_list)-min(self.__y_list)

    def get_manning(self):
        return self.__manning

    def get_K_over_tauc(self):
        return self.__K_over_tauc

    def get_tauc_over_rho(self):
        return self.__tauc_over_rho

    def is_upstream(self):
        return self.__is_upstream

    def is_downstream(self):
        return self.__is_downstream

    def get_up_section(self):
        return self.__up_section

    def get_down_section(self):
        return self.__down_section

    def set_down_section(self, down_section):
        if down_section.get_x() < self.__x:
            raise(ValueError("Error while setting downSection, it must be a greater or equal abscissa"))
        self.__down_section = down_section
    
    def set_up_section(self, up_section):
        if up_section.get_x() > self.__x:
            raise(ValueError("Error while setting upSection, it must be a lower or equal abscissa"))
        self.__up_section = up_section
    
    def set_is_downstream(self, boolean):
        self.__is_downstream = boolean
        if boolean:
            self.set_down_section(self)

    def set_is_upstream(self, boolean):
        self.__is_upstream = boolean
        if boolean:
            self.set_up_section(self)

    def get_s(self):
        raise(TypeError("get_s on a non trapezoidal section."))

    def get_granulometry(self):
        return self.__granulometry

    def get_points(self):
        return copy(self.__points)

    def get_S0(self, up_direction=False):
        """
        slope of the bottom of the channel. Can be computed toward the up_section if specified (toward the down_section by default)
        """
        if up_direction:
            if self.__is_upstream:
                dz = self.__z*1000 - self.__down_section.__z*1000
                dx = self.__down_section.__x - self.__x
            else:
                dz = self.__up_section.__z*1000 - self.__z*1000
                dx = self.__x - self.__up_section.__x
            return (dz/dx)/1000
        else:
            if self.__is_downstream:
                dz = self.__up_section.__z*1000 - self.__z*1000
                dx = self.__x - self.__up_section.__x
            else:
                dz = self.__z*1000 - self.__down_section.__z*1000
                dx = self.__down_section.__x - self.__x
            return (dz/dx)/1000

    # operators overloading

    def __str__(self):
        return f'IrregularSection : x={self.__x}, z={self.__z}'

    # plot stuff

    def plot(self, y=None):
        fig, ax = plt.subplots()
        ax.plot(self.__x_list, np.array(self.__y_list) + self.__z)
        if y != None:
            wet_points = self.get_wet_section(y)
            wpx = [p[0] for p in wet_points]
            wpy = [p[1] for p in wet_points]
            ax.plot(wpx, np.array([y for _ in range(len(wpx))]) + self.__z, "b--", label="water depth")
            ax.fill_between(wpx, np.array(wpy) + self.__z, np.array(wpy)+y+self.__z, color="cyan")
        plt.title(f"irregular section at x = {self.__x}")
        return fig

    # static methods

    @staticmethod
    def is_rectangular():
        return False

    @staticmethod
    def interp(up_section, down_section, x=None):
        """
        Args : 
            upSection (IrregularSection) : upper section
            downSection (IrregularSection) : lower section
            x (float) : abscissa of the interpolation, by default at the middle between the two section
        Returns :
            interpolatedSection (same section type than up_section) : section with interpolated attributes. 
            
            However be careful : this section will have the same point array, granulometry than the upper section
        """
        if x==None:
            x = 0.5*(up_section.__x+down_section.__x)
        interpolated_section = up_section.copy()
        interpolated_section.__x = x
        interpolated_section.__z = np.interp(x, [up_section.__x, down_section.__x], [up_section.__z, down_section.__z])
        interpolated_section.__z_min = np.interp(x, [up_section.__x, down_section.__x], [up_section.__z_min, down_section.__z_min])
        interpolated_section.__manning = up_section.get_manning()
        interpolated_section.__K_over_tauc = up_section.get_K_over_tauc()
        interpolated_section.__tauc_over_rho = up_section.get_tauc_over_rho()
        interpolated_section.__points = copy(up_section.__points)
        interpolated_section.x_list = copy(up_section.__x_list)
        interpolated_section.y_list = copy(up_section.__y_list)
        interpolated_section.__up_section = up_section
        interpolated_section.__down_section = down_section
        interpolated_section.__is_downstream = False
        interpolated_section.__is_upstream = False
        interpolated_section.__granulometry = up_section.__granulometry # interp for granulometry ?

        return interpolated_section

