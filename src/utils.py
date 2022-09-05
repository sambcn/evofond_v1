import numpy as np

# CONSTANT VALUES

G = 9.80665
Y_MIN = 0.001 # To avoid 0 on water depth, we assume that there is at least one millimeter

# Functions

def hydrogrammeLavabre(Qmax,tm,alpha,Qbase,t):
    """
    returns water discharge array matching t such that there is a maximum Qmax reached at tm by a alpha degree curve.
    It tends to Qbase as t tends to 0 or +inf.
    """
    Qbase = max(Qbase, 0.01)
    if alpha < 0:
        print(f"ERROR : while building lavabre hydrogram, alpha cannot be negative ({alpha})")
        raise ValueError()
    if not(tm >= t[0] and tm <= t[-1]): 
        print(f"ERROR : while building lavabre hydrogram, tm ({tm}) must be between t_begin ({t[0]}) and t_end ({t[-1]})")
        raise ValueError()
    if Qmax < Qbase:
        print(f"ERROR : while building lavabre hydrogram, Qm ({Qmax}) must be greater or equal to Qb ({Qbase})")
        raise ValueError()
    Q=np.array([(Qmax-Qbase)*2*np.power(ti/tm,alpha)/(1+np.power(ti/tm,2*alpha))+Qbase for ti in t])
    return np.array(Q)

def get_centroid(points):
    """returns gravity center of points"""
    x_coords = [p[0] for p in points]
    y_coords = [p[1] for p in points]
    _len = len(points)
    centroid_x = sum(x_coords)/_len
    centroid_y = sum(y_coords)/_len
    return (centroid_x, centroid_y)


def check_answer(answer, list_of_accepted_answers):
    while not(answer in list_of_accepted_answers):
        print("invalid answer, please be careful to uppercase letters.")
        answer = str(input(f"please your answer must be one of these ones {list_of_accepted_answers} : "))
    return answer

def input_float(question, positive=True):
    while True:
        try:
            answer = float(input(question))
            if positive and answer < 0:
                raise Exception()
                rais
            break
        except Exception:
            print(f"Write a {'positive' if positive else ''} float please.")
    return answer


def input_int(question, positive=True):
    while True:
        try:
            answer = int(input(question))
            if positive and answer < 0:
                raise Exception()
                rais
            break
        except Exception:
            print(f"Write a {'positive' if positive else ''} integer please.")
    return answer


def parse_datafile(filename):
    """
    read a txt file with columns x1, x2, ..., xn and return this data : [x1_list, x2_list, ..., xn_list]
    """
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
            if len(lines)==0:
                return []
            data_list = [[] for _ in range(len(lines[0].split()))]
            for k, line in enumerate(lines):
                data = line.split()
                if len(data) != len(data_list):
                    print(f"ERROR : line {k} has not the same amount of column than the first line ({filename})")
                    exit()
                try:
                    for i, x in enumerate(data):
                        data_list[i].append(float(x))
                except ValueError:
                    if k==0:
                        continue
                    print(f"WARNING : line {k+1} has some non-float data ({filename})")
            return data_list
    except FileNotFoundError:
        raise FileNotFoundError()

def write_datafile(filename, column_name_list, data):
    """
    write a txt file with as mush column as element in column_name_list (which contains columns' name).
    """
    if len(column_name_list) != len(data):
        print(f"ERROR : could not write {filename} because column_name_list has not the same length ({len(column_name_list)} than data {len(data)})")
        return
    if len(data)==0:
        print(f"WARNING : nothing to write in {filename}")
        return
    with open(filename, 'w') as f:
        f.writelines([d + " " for d in column_name_list]+["\n"])
        for i in range(len(data[0])):
            line = ""
            for c in range(len(data)):
                try:
                    line += f"{data[c][i]} "
                except IndexError:
                    print(f"ERROR : while writing {filename}, all the data must be of the same length.")
                    return
            line += "\n"
            f.writelines(line)
    return
        

def inter_xy(list_xy):
    """
    takes as input a list of plots (x1, y1), (x2, y2)... and returns same data with the common abscissa only :
    (x, y1bis), (x, y2bis),...
    """
    current_x = list_xy[0][0][:]
    for i in range(1, len(list_xy)):
        x = list_xy[i][0]
        for xi in current_x:
            if not(xi in x):
                current_x.pop(current_x.index(xi))
    result = []
    for i in range(len(list_xy)):
        x, y = list_xy[i]
        new_y = []
        for xi in current_x:
            index = x.index(xi)
            new_y.append(y[index])
        result.append((current_x, new_y))
    return result


def get_matrix_max(matrix):
    """return the maximum of a matrix. Used to set windows size in animation plots."""
    try:
        maxi = max(matrix[0])
        for col in matrix:
            if max(col) > maxi:
                maxi = max(col)
        return maxi
    except IndexError as e:
        raise IndexError(f"Matrix must have at least one element. {e}")

def time_to_string(t):
    """return a string to print a given time t in seconds with the format Ah Bm Cs"""
    hours = int(t) // 3600
    minutes = (int(t) - 3600*hours) // 60
    seconds = t%60
    if hours == 0:
        if minutes == 0:
            return f"{seconds:.2f}s"
        else:
            return f"{minutes}min {seconds:.2f}s"
    else:
        return f"{hours}h {minutes}min {seconds:.2f}s"

def read_hecras_data(filename):
    """
    returns two lists x, y which can be plot, by reading the filename which must be like :
    x0 y0
    x1 y1
    .. ..
    xn yn
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
        x = []
        y = []
        for line in lines:
            data = line.split("\t")
            x.append(float(data[0]))
            y.append(float(data[1]))
        return x, y

def reverse_data(x, y):
    """
    returns x' y' such that the curve y(x) is mirrored.
    """
    new_x = []
    new_y = []
    for i in range(len(x)-1, -1, -1):
        new_x.append(x[-1] - x[i])
        new_y.append(y[i])
    return new_x, new_y

def real_roots_cubic_function(a, b, c, d):
    """
    Computes the root(s) of a cubic function described by the argument coeff.
    If coeff=[a, b, c, d], the function is a*x**3 + b*x**2 + c*x + d.
    It returns the list of positive real roots.
    """
    p = (3*a*c - b**2) / (3*a**2)
    q = (2*b**3 - 9*a*b*c + 27*a**2 *d) / (27*a**3)
    delta1 = q**2 + (4*p**3) / (27)
    roots = []
    if delta1 < 0:
        j = complex(0, 1)
        x1 = ((-q-j*(-delta1)**0.5)/2)**(1/3) + ((-q+j*(-delta1)**0.5)/2)**(1/3) - (b/(3*a))
        if x1.imag < 1e-5:
            x1 = x1.real
            if x1 > 0:
                roots.append(x1)
        else:
            print("complex roots !")
    else:
        x1 = np.sign(-q-delta1**0.5)*(abs(-q-delta1**0.5)/2)**(1/3) + np.sign(-q+delta1**0.5)*(abs(-q+delta1**0.5)/2)**(1/3) - (b/(3*a))

        if x1 > 0:
            roots.append(x1)
    a2 = a
    b2 = b + a*x1
    c2 = c + b2*x1
    delta2 = b2**2 - 4*a2*c2
    if delta2 < -1e-8:
        return roots
    else:
        delta2 = max(0, delta2)
        x2 = (-b2-delta2**0.5)/(2*a2)
        if x2 > 0:
            roots.append(x2)
        x3 = (-b2+delta2**0.5)/(2*a2)
        if x3 > 0:
            roots.append(x3)
    roots.sort()
    return roots