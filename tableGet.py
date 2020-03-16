
import pickle
from astropy.table import QTable
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

def tableGet(args):
    try:
        with open(getTableName(args), 'rb') as table_file:
            table = pickle.load(table_file)
            return table
    except:
        return tableCreate(args)

def tableCreate(args):
    with open('pulsars.table', 'rb') as pulsarDoc:
        pulsars = pickle.load(pulsarDoc)


    if "hasDM" in args:
        row = 0
        rows = []
        for i in pulsars.iterrows():
            if type(i[9]) == np.ma.core.MaskedConstant:
                rows.append(row)
            row +=1
        pulsars.remove_rows(rows)

    for i in args:

        if "posErrorASec" in i:
            row = 0
            rows = []
            error = int(i[-3:-1])
            for j in pulsars.iterrows():
                error_RA = j[4]
                error_DEC = j[7] #j is pulsar 

                if error_RA >= error or error_DEC >= error:
                    rows.append(row)
                row += 1
            pulsars.remove_rows(rows)
        if "posRange" in i:
            ranges = i[9:-1].split(',')
            for j in range(len(ranges)):
                range_part = list((int(k) for k in ranges[j].split(".")[0].split(":")))
                try:
                    range_part.append(float("0." + ranges[j].split(".")[1]))
                except:
                    pass
                if j == 2 or j == 3:
                    range_part[0] += 90
                ranges[j] = range_part

            rows = []
            row = -1
            for j in pulsars.iterrows():
                row += 1

                pos_RA = list((int(k) for k in str(j[3]).split(".")[0].split(":"))) # change to using skycoord objects
                try:
                    pos_RA.append(int(j[0].split(".")[1]))
                except:
                    pass
                pos_DEC = list((int(k) for k in str(j[6]).split(".")[0][1:].split(":")))
                if j[2][0] == "+":
                    pos_DEC[0] += 90
                else:
                    pos_DEC[0] = 90 - pos_DEC[0]

                try:
                    pos_DEC.append(int(j[2].split(".")[1]))
                except:
                    pass
                for k in range(len(pos_RA)):

                    if int(ranges[0][k]) > pos_RA[k]:
                        rows.append(row)
                        break
                    if int(ranges[0][k] != pos_RA[k]):
                        break

                if row in rows:
                    continue

                for k in range(len(pos_RA)):

                    if int(ranges[1][k]) < pos_RA[k]:

                        rows.append(row)
                        break
                    if int(ranges[1][k] != pos_RA[k]):
                        break

                if row in rows:
                    continue

                for k in range(len(pos_DEC)):
                    if int(ranges[2][k]) > pos_DEC[k]:
                        rows.append(row)
                        break
                    elif int(ranges[2][k] != pos_DEC[k]):
                        break
                if row in rows:
                    continue

                for k in range(len(pos_DEC)):
                    if int(ranges[3][k]) < pos_DEC[k]:
                        rows.append(row)
                        break
                    elif int(ranges[3][k]) != pos_DEC[k]:
                        break

            pulsars.remove_rows(rows)

    with open(getTableName(args), 'wb') as pulsarDoc:
        pickle.dump(pulsars, pulsarDoc)
    return pulsars

def tableRenew(args):
    pass

def getTableName(args):
    name ="pulsars.table."
    for i in args:
        name += i + ","
    return name



def fits_crop(np_array, sky_coords, wanted_cords, size_pic, size_area_wanted, coords_top_left = True, coords_center = False):
    if coords_top_left:
        pos_dif = [wanted_cords[0]-sky_coords[0], sky_coords[1]-wanted_cords[1]]
    pos_dif[0] /= size_pic[0] / len(np_array[0])
    pos_dif[1] /= size_pic[1] / len(np_array)
    pos_dif[0] = int(pos_dif[0])
    pos_dif[1] = int(pos_dif[1])

    max = [pos_dif[0] + (size_area_wanted / size_pic[0]) * len(np_array[0]),  pos_dif[1] - (size_area_wanted / size_pic[1]) * len(np_array)] # make sure to test and shit
    if max[0] > len(np_array[0]):
        max[0] = len(np_array)
    if max[1] < 0:
        max[1] = 0

    min =  [pos_dif[0] - (size_area_wanted / size_pic[0]) * len(np_array[0]),  pos_dif[1] + (size_area_wanted / size_pic[1]) * len(np_array)]

    if min[1] > len(np_array[0]):
        min[1] = len(np_array)
    if min[0] < 0:
        min[0] = 0

    max[0] = int(max[0])
    max[1] = int(max[1])

    min[0] = int(min[0])
    min[1] = int(min[1])


    reduced_array = np.array([np_array[max[1]][min[0]:max[0]]])

    print(min[1]-max[1])

    for i in range(min[1]-max[1]):
        if i == 0:
            continue
        reduced_array = np.append(reduced_array, [(np_array[i + max[1]][min[0]:max[0]])], 0)

    return reduced_array


def range_check(pos1, pos2, maxD):
    pos1[0] = list((int(i) for i in pos1[0].split(".")[0].split(":")))
    pos1_holder = pos1[1]
    pos1[1] = list((int(i) for i in pos1[1].split(".")[0][1:].split(":")))

    if pos1_holder[1] == "+":
        pos1[1][1] += 90
    else:
        pos1[1][1] = 90 - pos1[1][1]



    pos2[0] = list((int(i) for i in pos2[0].split(".")[0].split(":")))
    pos2_holder = pos2[1]
    pos2[1] = list((int(i) for i in pos2[1].split(".")[0][1:].split(":")))

    if pos2_holder[1] == "+":
        pos2[1][1] += 90
    else:
        pos2[1][1] = 90 - pos2[1][1]


    if abs(pos1[0][0] - pos2[0][0]) > 2:
        return False
    if abs(pos1[1][0] - pos2[1][0]) > 2:
        return False
    if abs((pos1[0][0] * 60 * 60) + pos1[0][1] * 60 + pos1[0][2] - pos2[0][0] * 60  * 60 - pos2[0][1] * 60 + pos2[0][2]) > maxD:
        return False
    if abs((pos1[1][0] * 4 * 60) + pos1[1][1] * 60 + pos1[1][2] - pos2[1][0] * 4  * 60 - pos2[1][1] * 60 + pos2[1][2]) > maxD:
        return False
    return True


#table = tableGet(["V[test4]","hasDM", "posErrorASec[01]", "posRange[1:0:0.0,2:0:0.0,-90:0:0.0,90:0:0.0]"])

