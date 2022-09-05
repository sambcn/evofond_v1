#!/usr/bin/env python3

import datetime as dt
import argparse as argp
import logging as lg
import os
import shutil
import datetime
import json
from src.run import run_back
from src.utils import parse_datafile, input_float, input_int, check_answer


def main():

    parser = argp.ArgumentParser(description="Evofond executable. Project in progress...")
    parser.add_argument("-v", "--verbose", action="store_true", help="increase the amount of information printed during the execution and create a log file of the execution")
    parser.add_argument("-q", "--quickstart", action="store_true", help="help to start a new project")
    parser.add_argument("-l", "--list", action="store_true", help="list of existing projects")
    parser.add_argument("-d", "--delete", nargs='?', help="delete a given project")
    parser.add_argument("-c", "--copy", nargs="?", help="copy a given project")
    parser.add_argument("-m", "--modify", nargs='?', help="modify a given project")
    parser.add_argument("-r", "--run", nargs='?', help="run a given project")
    parser.add_argument("--clear", action="store_true", help="remove all the existing log files")
    parser.add_argument("--hydrau", nargs='?', help="hydraulic computation for a given water discharge")
    args = parser.parse_args()

    ### CLEAR LOG ###
    if args.clear:
        clear()

    ### LOG INITIALIZATION ###
    if args.verbose:
        e = dt.datetime.now()
        try:
            lg.basicConfig(filename=f"./log/{e.year}_{e.month}_{e.day}_{e.hour}_{e.minute}_{e.second}.log", encoding='utf-8', level=lg.DEBUG)
        except FileNotFoundError:
            os.mkdir("log")
            lg.basicConfig(filename=f"./log/{e.year}_{e.month}_{e.day}_{e.hour}_{e.minute}_{e.second}.log", encoding='utf-8', level=lg.DEBUG)

    ### DELETE ###
    if args.delete:
        delete(args.delete)

    ### LIST ###
    if args.list:
        ls()

    ### QUICKSTART ###
    if args.quickstart:
        quickstart() 

    ### MODIFY ###
    if args.modify:
        modify(args.modify) 

    ### COPY ###
    if args.copy:
        copy(args.copy) 
    
    ### RUN ###
    if args.run:
        run(args.run)

    ### HYDRAU ###
    if args.hydrau:
        run(args.hydrau, hydrau=True)

    return

def ls():
    """
    print project list
    """
    try:
        print("project list : ")
        print(os.listdir("./projects"))
    except FileNotFoundError:
        print("no folder called 'projects', no projects found")

def copy(project_name):
    """
    copy a given project called project_name
    """
    try:
        os.chdir("./projects")
    except FileNotFoundError:
        print("there is no project to copy.")
        return
    if os.path.isdir(f"./{project_name}"):
        shutil.copytree(f"./{project_name}", f"./{project_name}_copy")
    else:
        print(f"there is no project called {project_name}. The available project list is : {os.listdir('./')}")
    os.chdir("..")
    return


def delete(project_name):
    """
    delete a given project called prokect_name
    """
    try:
        os.chdir("./projects")
    except FileNotFoundError:
        print("there is no project to delete.")
        return
    if os.path.isdir(f"./{project_name}"):
        shutil.rmtree(f"./{project_name}")
    else:
        print(f"there is no project called {project_name}. The available project list is : {os.listdir('./')}")
    os.chdir("..")
    return

def clear():
    """
    clear log files (in the log folder)
    """
    files_in_directory = os.listdir("./log")
    filtered_files = [file for file in files_in_directory if file.endswith(".log")]
    for file in filtered_files:
        path_to_file = os.path.join("./log", file)
        os.remove(path_to_file)

def modify(project_name):
    """
    help users to make some classical modifications on
    """
    try:
        os.chdir("./projects")
    except FileNotFoundError:
        print("no folder called 'projects', no projects found")
        return
    try:
        os.chdir(project_name)
    except FileNotFoundError:
        print(f"no project called {project_name}")
        return
    print(f"[modify] profile chosen : {json.load(open(f'{project_name}_conf.json', 'r'))['PROFILE_PATH']}")
    print("[modify] the available modifications are the following : ")
    print(" 1 - [adjust elevation] you can add or substract a given height to all sections' elevation")
    print(" 2 - [adjust width] you can add or substract a given length to all sections' width")
    print(" 3 - [adjust minimal heigt] you can add or substract a given height to all sections' minimal height")
    print(" 4 - [set minimal heigt] initialize minimal height data (default value of zmin will be z)")
    print(" 5 - [set granulometry to section] initialize granulometry association data (default granulometry will be granulometry_1)")
    print(f"[modify] for any other modification, you can check the {project_name}_conf.json file and change configuration on your own (be careful..) ")
    answer = input("[modify] write the number of the target modification : ")
    answer = int(check_answer(answer, [f"{i}" for i in range(1, 6)]))
    if answer == 1:
        adjust_data(json.load(open(f"{project_name}_conf.json", 'r'))["PROFILE_PATH"], 'z')
    elif answer == 2:
        adjust_data(json.load(open(f"{project_name}_conf.json", 'r'))["PROFILE_PATH"], 'b')
    elif answer == 3:
        adjust_data(json.load(open(f"{project_name}_conf.json", 'r'))["PROFILE_PATH"], 'zmin')
    elif answer == 4:
        set_data(json.load(open(f"{project_name}_conf.json", 'r'))["PROFILE_PATH"], 'zmin', 'z')
    elif answer == 5:
        set_data(json.load(open(f"{project_name}_conf.json", 'r'))["PROFILE_PATH"], 'granulometry', 1)

def set_data(profile_path, data_name, default_value):
    """
    used by 'modify' function to create a new data column in the profile
    """
    data = parse_datafile(profile_path)
    try:
        first_line = open(profile_path, 'r').readlines()[0].split()
    except IndexError:
        return
    try:
        data_index = first_line.index(data_name)
    except ValueError:
        data_index = len(first_line)
        data.append([])
        first_line.append(data_name)
    if type(default_value)==str:
        try:
            new_data = data[first_line.index(default_value)] 
        except ValueError:
            print(f"ERROR : while setting new data {data_name}, the default value {default_value} not found in other columns")
    elif type(default_value)==float or type(default_value)==int:
        new_data = [default_value for _ in range(len(data[0]))]
    data[data_index] = new_data
    nb_col = len(data)
    nb_line = len(data[0])
    with open(profile_path, 'w') as f:
        f.writelines([d+" " for d in first_line]+["\n"])
        for i in range(nb_line):
            line = ""
            for j in range(nb_col):
                line += f"{data[j][i]} "
            f.writelines(line+"\n")
    return

def adjust_data(profile_path, data_name):
    """
    used by 'modify' function to modify a given data column in the profile
    """
    delta = input_float(f"[modify] please choose your delta {data_name} (m) : ", positive=False)
    data = parse_datafile(profile_path)
    try:
        first_line = open(profile_path, 'r').readlines()[0]
    except IndexError:
        return
    if len(data)==0:
        return
    try:
        data_index = first_line.split().index(data_name)
    except ValueError:
        print(f"ERROR : no '{data_name}' column found when trying to adjust this data")
        return

    nb_col = len(data)
    nb_line = len(data[0])
    with open(profile_path, 'w') as f:
        f.writelines(first_line)
        for i in range(nb_line):
            line = ""
            for j in range(nb_col):
                line += f"{data[j][i]+delta*int(j==data_index)} "
            f.writelines(line+"\n")
    return


def quickstart():
    """
    function used by users to initialize a new project
    """
    try:
        os.chdir("projects")
    except FileNotFoundError:
        os.mkdir("projects")
        lg.info("directory 'projects' created successfully")
        os.chdir("projects")
    project_name = str(input("[QS] choose a project name : "))
    while True: 
        try:
            os.mkdir(project_name)
            os.chdir(project_name)
            lg.info(f"directory ./projects/{project_name} created successfully")
            break
        except FileExistsError:
            print("[QS] this project already exist : ")
            os.listdir()
            project_name = str(input("[QS] please choose another one : "))
    print(f"[QS] starting initialization of {project_name}_conf.json")
    data_files = []
    conf_dict = {"NAME":project_name}
    section = str(input("[QS] [profile] what kind of section ? [rectangular/trapezoidal] "))
    section = check_answer(section, ["rectangular", "trapezoidal"])
    conf_dict["SECTION"]=section
    f = open("./profile.txt", 'w')
    conf_dict["PROFILE_PATH"]="./profile.txt"
    column_title = "x z b"
    if section == "trapezoidal":
        column_title += " s"
    # column_title += " zmin granulo"
    f.writelines(column_title)
    f.close()
    answer = str(input("[QS] [hydrogram] do you want to initialize a Lavabre hydrogram ? [yes/no] "))
    answer = check_answer(answer, ["yes", "no"])
    conf_dict["LAVABRE"] = (answer=='yes')
    if answer == "yes":
        answer = input_float("[QS] [hydrogram] [Lavabre]  how long last the event ? [float answer expected (seconds)] ")
        conf_dict["DURATION"]=answer
        answer = input_float("[QS] [hydrogram] [Lavabre]  when is the maximum water discharge reached ? [float answer expected (seconds)] ")
        conf_dict["TM"]=answer
        answer = input_float("[QS] [hydrogram] [Lavabre]  what is the maximum water discharge ? [float answer expected (m3/s)] ")
        conf_dict["QM"]=answer
        answer = input_float("[QS] [hydrogram] [Lavabre]  what is the minimum water discharge ? [float answer expected (m3/s)] ")
        conf_dict["QB"]=answer
        answer = input_int("[QS] [hydrogram] [Lavabre]  what is the alpha exponent (which define the curve variation) ? [int answer expected] ")
        conf_dict["ALPHA"]=answer
        answer = input_float("[QS] [hydrogram] [Lavabre]  choose your time step ? [float answer expected] ")
        conf_dict["DT"]=answer
        conf_dict["HYDROGRAM_PATH"]=None
    else:
        conf_dict["DURATION"] = None
        conf_dict["TM"] = None
        conf_dict["QM"] = None
        conf_dict["QB"] = None
        conf_dict["ALPHA"] = None
        conf_dict["DT"] = None
        print("[QS] [hydrogram] \t \t manual hydrogram chosen, please fill hydrogram.txt with (t (s), Q (m3/s)) values")
        open("./hydrogram.txt", 'w')
        data_files.append("hydrogram.txt")
        conf_dict["HYDROGRAM_PATH"]="./hydrogram.txt"

    answer = str(input("[QS] [sediment transport] which sediment transport law do you want ? [Rickenmann1990/Rickenmann1991/Lefort2015/LefortSogreah1991/Meunier1989/MeyerPeter1948/PitonRecking2017/Piton2016] "))
    answer = check_answer(answer, ["Rickenmann1990", "Rickenmann1991", "Lefort2015", "LefortSogreah1991", "Meunier1989", "MeyerPeter1948", "PitonRecking2017", "Piton2016"])
    conf_dict["TRANSPORT_LAW"]=answer
    answer = input_float("[QS] [sediment transport] chose an upstream slope for sediment transport : [float expected (%)] ")
    conf_dict["UPSTREAM_SLOPE"]=answer    
    answer = input_float("[QS] [sediment transport] chose an upstream width for sediment transport : [float expected (m)] ")
    conf_dict["UPSTREAM_WIDTH"]=answer

    granulometry_files = []
    new_granulo = True
    nb_granulo = 0
    while new_granulo:
        nb_granulo += 1
        granulo_name = f"granulometry_{nb_granulo}.json"
        granulo_dict = dict()
        granulo_file = open(granulo_name, 'w')
        granulometry_files.append(granulo_name)
        print(f"[QS] [sediment transport] [granulometry {nb_granulo}] If you want to fill it later you can write any answer and change it later in 'granulometry_{nb_granulo}.json'")
        answer = input_float("[QS] [sediment transport] [granulometry_1.json] give dm value (mean diameter) : [float expected]")
        granulo_dict["dm"] = answer
        answer = input_float("[QS] [sediment transport] [granulometry_1.json]  give d30 value : [float expected]")
        granulo_dict["d30"] = answer
        answer = input_float("[QS] [sediment transport] [granulometry_1.json]  give d50 value : [float expected]")
        granulo_dict["d50"] = answer
        answer = input_float("[QS] [sediment transport] [granulometry_1.json]  give d90 value : [float expected]")
        granulo_dict["d90"] = answer
        answer = input_float("[QS] [sediment transport] [granulometry_1.json]  give d84tb value (travelling bedload) : [float expected]")
        granulo_dict["d84tb"] = answer
        answer = input_float("[QS] [sediment transport] [granulometry_1.json]  give d84bs value (bed surface) : [float expected]")
        granulo_dict["d84bs"] = answer
        answer = input_float("[QS] [sediment transport] [granulometry_1.json]  give Gr value (Gradation) : [float expected]")
        granulo_dict["Gr"] = answer
        json.dump(granulo_dict, granulo_file, indent=6)
        granulo_file.close()
        answer = str(input("[QS] [sediment transport] [granulometry] do you want to initialize a new granulometry ? [yes/no]"))
        answer = check_answer(answer, ["yes", "no"])
        new_granulo = (answer == "yes")
    if len(granulometry_files) > 1:
        print("[QS] [sediment transport] [granulometry] the first granulometry will be by default the granulometry of all sections, use -m option to change this")
    conf_dict["GRANULOMETRY_FILES"] = granulometry_files
    answer = str(input("[QS] [interpolation] do you want to interpolate your profile ? [yes/no]"))
    answer = check_answer(answer, ["yes", "no"])
    conf_dict["INTERPOLATION"] = answer=="yes"
    if answer == "yes":
        dx = input_float("[QS] [interpolation] choose your dx : ")
        conf_dict["DX"] = dx
    else:
        conf_dict["DX"] = None
    
    answer = str(input("[QS] [model] do you want to use critical depth hypothesis ? [yes/no] "))
    answer = check_answer(answer, ["yes", "no"])
    conf_dict["CRITICAL"]=answer=="yes"
    if answer == "no":
        answer = str(input("[QS] [model] then choose your friction law ? [Manning-Strickler/Ferguson] "))
        answer = check_answer(answer, ["Manning-Strickler", "Ferguson"])
        conf_dict["FRICTION_LAW"]=answer
    else:
        conf_dict["FRICTION_LAW"]=None

    if conf_dict["CRITICAL"]:
        conf_dict["UPSTREAM_CONDITION"]="critical_depth"
        conf_dict["DOWNSTREAM_CONDITION"]="critical_depth"
    else:
        answer = str(input("[QS] [Boundaries] choose your upstream boundary condition : [critical_depth/normal_depth] "))
        answer = check_answer(answer, ["critical_depth", "normal_depth"])
        conf_dict["UPSTREAM_CONDITION"]=answer
        answer = str(input("[QS] [Boundaries] choose your downstream boundary condition : [critical_depth/normal_depth] "))
        answer = check_answer(answer, ["critical_depth", "normal_depth"])
        conf_dict["DOWNSTREAM_CONDITION"]=answer

    answer = input_float("[QS] [SPEED_COEF] choose a SPEED_COEF value (please check documentation to check what this is, choose 1 in case you have no idea) : [float expected] ")
    conf_dict["SPEED_COEF"]=answer

    answer = input_float("[QS] [RESULTS BACKUP] choose your time step for the results : [float answer expected (seconds)] ")
    conf_dict["BACKUP_TIME_STEP"]=answer

    answer = str(input("[QS] [Perf] do you want to know the execution performance (it may make it last a little longer) ? : [yes/no] "))
    answer = check_answer(answer, ["yes", "no"])
    conf_dict["PERF"]=(answer=="yes")


    outfile = open(f"{project_name}_conf.json", 'w')
    json.dump(conf_dict, outfile, indent=6)
    outfile.close()
    print(f"[QS] {project_name}_conf.json intialized")
    print(f"[QS] please now fill the following data files : {data_files} \n note that you can use -m option to change quickly some configuration data.")
    return

def run(project_name, hydrau=False):
    """
    run a given project called project_name. If hydrau==True, it will only compute water depth and not the entire simulation.
    """
    try:
        os.chdir("./projects")
    except FileNotFoundError:
        print("there is no project to run.")
        return

    if os.path.isdir(f"./{project_name}"):
        os.chdir(f"./{project_name}")
    else:
        print(f"there is no project called {project_name}. The available project list is : {os.listdir('./')}")
        os.chdir("..")
        return

    try:
        f = open(f"{project_name}_conf.json", 'r')
        args_dict = json.load(f)
        f.close()
    except FileNotFoundError:
        print(f"ERROR : no conf file found. Please be sure that your conf file name is '{project_name}_conf.json'")
        return
    water_discharge = input_float("choose a water discharge (m3/s) [float expected] : ") if hydrau else None

    run_back(args_dict, hydrau=water_discharge)
    
    return



if __name__=='__main__':
    main()