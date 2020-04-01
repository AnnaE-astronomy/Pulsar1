from astropy.coordinates import SkyCoord
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
        rows = []

        for i in range(len(pulsars)):
            if type(pulsars["DM"][i]) == np.ma.core.MaskedConstant:
                rows.append(i)
        pulsars.remove_rows(rows)

    for i in args:

        if "posErrorASec" in i:
            rows = []
            error = int(i[-3:-1])
            for j in range(len(pulsars)):
                error_RA = pulsars["RAJ_ERR"][j]
                error_DEC = pulsars["DECJ_ERR"][j]

                if error_RA >= error or error_DEC >= error:
                    rows.append(j)
            pulsars.remove_rows(rows)
        if "posRange" in i:
            ranges = i[9:-1].split(',')
            print(ranges)
            coord_min = SkyCoord(ranges[0], ranges[2], unit=("hour","deg"))
            coord_max = SkyCoord(ranges[1], ranges[3], unit=("hour","deg"))


            rows = []
            coords = SkyCoord(pulsars['RAJ'], pulsars['DECJ'], unit=('hour', 'deg'))
            for j in range(len(pulsars)):

                if coords[j].ra.deg < coord_min.ra.deg or coords[j].ra.deg > coord_max.ra.deg or \
                        coords[j].dec.deg < coord_min.dec.deg or coords[j].dec.deg > coord_max.dec.deg:
                    rows.append(j)

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

print("STARTING")
table = tableGet(["V[test23]","hasDM", "posErrorASec[01]", "posRange[1:0:0.0,23:59:59.59,-90:0:0.0,90:0:0.0]"])
print(table[0])
for i in table["DM"]:
    #print(type(i))
    pass
