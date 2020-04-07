from astropy.coordinates import SkyCoord
import pickle
from astropy.table import QTable
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

# Checks if table of specific set of pulsars exists and loads if does, otherwise creates it, saves it for future use, and loads.
def tableGet(args):
    try:
        with open("saved_reductions/" + getTableName(args), 'rb') as table_file:
            table = pickle.load(table_file)
            return table
    except:
        return tableCreate(args)


# Takes existing table of all pulsars and cuts specific pulsars out based on specs from position to whether or not it has a DM on file

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

    with open("saved_reductions/" + getTableName(args), 'wb') as pulsarDoc:
        pickle.dump(pulsars, pulsarDoc)
    return pulsars


# The table name is just all of the inputs added together into a long string.
def getTableName(args):
    name ="pulsars.table."
    for i in args:
        name += i + ","
    return name





