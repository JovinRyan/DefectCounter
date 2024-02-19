import pandas as pd 
import numpy as np 
import os 

IDEAL_BOND_LENGTH = 2.36 #Angstrom

#MIN_DIST = IDEAL_BOND_LENGTH/3
#MAX_DIST = 2*IDEAL_BOND_LENGTH


def ExcitationCount(InputFile, IdealBondLength = IDEAL_BOND_LENGTH, Path = ""):

    MIN_DIST = IdealBondLength/2
    OutputFile = InputFile.split(".")[0] + ".csv"

    datf = open(os.path.join(Path, InputFile), "r").read().split("\n")
    if datf[-1] == "": datf.pop(-1)                                                  #There's a "" character at the end of the file for some reason
    ATOM_COUNT = int(datf[0])
    COL_NUM = 4                                                                      #4 is the default for .xyz files 

    indexes = [i for i, x in enumerate(datf) if x == datf[0]]

    indexes.reverse()
    TIMESTEP_LIST = []

    for i in indexes:
        TIMESTEP_LIST.append(datf[i+1])                                                                #Removes "Atom Count" and "TIme Step" lines from .xyz file  
        datf.pop(i+1)
        datf.pop(i)

    TIMESTEP_LIST.reverse()

    for i in range(len(TIMESTEP_LIST)):
        TIMESTEP_LIST[i] = int(TIMESTEP_LIST[i].split(" ")[-1])

    print("Save File: " + OutputFile)

    TypeList = []
    XList = []
    YList = []
    ZList = []

    for i in range(len(datf)):                                                       #Converting List to Dataframe, **try and optimize this better**
        datf[i] = datf[i].split(" ")
        datf[i][0] = int(datf[i][0])
        TypeList.append(datf[i][0])
        datf[i][1] = float(datf[i][1])
        XList.append(datf[i][1])
        datf[i][2] = float(datf[i][2])
        YList.append(datf[i][2])
        datf[i][3] = float(datf[i][3])
        ZList.append(datf[i][3])

    FRAMES = int(len(datf)/ATOM_COUNT)                                                    #Gets number of frames from cleaned list of data

    DATA_DF = pd.DataFrame({"Type" : TypeList, "X" : XList, "Y" : YList, "Z" : ZList})

    INIT_STRUCTURE = DATA_DF[:ATOM_COUNT]


    VACANCY_DEFECT = []
    #INTERSTITIAL_DEFECT = []
    #SUBSTITUTION_DEFECT = []

    for i in range(1, FRAMES):
        Temp_Struct = DATA_DF[ATOM_COUNT*i: ATOM_COUNT*(i+1)].reset_index(drop=True)    #Accesses "Frames". Number of atoms must be constant
        VACANCY_DEFECT.append(0) 
        for j in range(len(INIT_STRUCTURE)):
            if np.sqrt((Temp_Struct.loc[j, "X"] - INIT_STRUCTURE.loc[j, "X"])**2 + (Temp_Struct.loc[j, "Y"] - INIT_STRUCTURE.loc[j, "Y"])**2 + (Temp_Struct.loc[j, "Z"] - INIT_STRUCTURE.loc[j, "Z"])**2) > MIN_DIST :
                VACANCY_DEFECT[i-1] += 1
        print("Frame " + str(i))

    VACANCY_DEFECT.insert(0, 0)

    (pd.DataFrame({"TimeStep" : TIMESTEP_LIST, "Excitations" : VACANCY_DEFECT})).to_csv(os.path.join(Path, OutputFile), index=False)

Directory = r"D:\ZirconiumCarbide2\ProperOutputs\Zr"

Files = os.listdir(Directory)

# for i in range(len(Files)):
#     if Files[i].split(".")[-1] != "xyz":
#         Files.pop(i)

# for File in Files:
#     ExcitationCount(File, Path=Directory)

ExcitationCount()
