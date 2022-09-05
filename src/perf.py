import time as t
import numpy as np
from prettytable import PrettyTable
from src.utils import time_to_string

class Performance:
    """
    class used for storing times spent by each function using the decorator @measure_perf
    """
    dict_of_perf = dict()
    time_start = None
    time_end = None

    @staticmethod
    def get_perf_table():
        """
        return the table of the stored performances
        """
        if Performance.time_end == None:
            print("automatic stop of perf measure.")
            Performance.time_end = t.time()
        te = Performance.time_end
        ts = Performance.time_start
        table = PrettyTable(['method name', '% of the execution time', 'total time spent (s)', 'number of call', 'mean time spent by call (s)'])
        for key, value in Performance.dict_of_perf.items():
            table.add_row([key, f"{100*(value[0]/(te-ts)):.2f}%", float(f"{value[0]:.8f}"), f"{value[1]}", f"{value[0]/value[1]:.8f}"])
        return table

    @staticmethod
    def print_perf():
        ts = Performance.time_start
        te = Performance.time_end
        if ts == None:
            print("nothing measured")
        else:
            if te == None:
                print("automatic stop of perf measure.")
                Performance.time_end = t.time()
                te = Performance.time_end
            print(f"total time : {time_to_string(te-ts)} \n")
            table = Performance.get_perf_table()
            # for key, value in Performance.dict_of_perf.items():
            #     table.add_row([key, f"{100*(sum(value)/(te-ts)):.2f}%", f"{sum(value):.2f}s", f"{len(value)}", f"{np.mean(value)}s"])
            print(table.get_string(sortby="total time spent (s)", reversesort=True))

    @staticmethod
    def save_perf(filename):
        table = Performance.get_perf_table()
        with open(filename, 'w') as f:
            f.writelines(f"total time = {Performance.time_end - Performance.time_start}\n")
            f.writelines(table.get_string(sortby="total time spent (s)", reversesort=True))

    @staticmethod
    def start():
        Performance.time_start = t.time()

    @staticmethod
    def stop():
        Performance.time_end = t.time()

    @staticmethod
    def measure_perf(func):
        """
        decorator home-made in order to store time spent in the decorated method during the full execution
        """
        def wrapper(*args, **kargs):
            if Performance.time_start == None:
                return func(*args, **kargs)
            ts = t.time()
            result = func(*args, **kargs)
            te = t.time()
            if func.__name__ in Performance.dict_of_perf:
                Performance.dict_of_perf[func.__name__] = Performance.dict_of_perf[func.__name__][0]+(te-ts), Performance.dict_of_perf[func.__name__][1]+1
            else:
                Performance.dict_of_perf[func.__name__] = (te-ts), 1
            return result
        return wrapper
