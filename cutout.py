import sys
from astropy.table import Table, vstack, Column, QTable
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Polygon
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.visualization import ZScaleInterval, ImageNormalize
from PIL import Image
import scipy.stats
import numpy as np

import math


def cutOut(fits_file, radius, pulsarDataASKAP, pulsarDataStokesV, pulsarDataStokesI, stokesVersion ="", stokesType=""):
    secondMarkers = False
    if stokesType == "both":
        pulsarDataStokesVislands = pulsarDataStokesV[1]
        pulsarDataStokesV = pulsarDataStokesV[0]
        pulsarDataStokesIislands = pulsarDataStokesI[1]
        pulsarDataStokesI = pulsarDataStokesI[0]
        secondMarkers = True

    # Make a SkyCoord object for the source position
    ra_unit = u.hourangle
    position = SkyCoord(ra=pulsarDataASKAP["RAJ"], dec=pulsarDataASKAP["DECJ"], unit=(ra_unit, u.deg))

    # Open the FITS image and make a World Coordinate System object
    try:
        with fits.open(fits_file) as hdul:
            data, header = hdul[0].data, hdul[0].header
            wcs = WCS(header, naxis=2)
    except Exception as e:
        print(e, fits_file)
        return
    # Make a postage stamp cutout at the area of interest
    # and change the data units from Jy to mJy
    cutout = Cutout2D(data, position, radius * u.deg, wcs=wcs)
    data = cutout.data * 1000
    wcs = cutout.wcs

    # Make a figure using the cutout WCS
    fig = plt.figure()
    ax = fig.add_subplot(111, projection=wcs)

    # Stretch the data for better contrast and plot the image
    norm = ImageNormalize(data, interval=ZScaleInterval(contrast=0.2))
    im = ax.imshow(data, cmap=plt.cm.gray_r, norm=norm)
    fig.colorbar(im, label=r"Flux Density (mJy beam$^{-1}$)")
    ax.set_ylabel("Dec (J2000)")
    ax.set_xlabel("RA (J2000)")
    ax.set_title("Graph of Flux Density vs Position for: " + pulsarDataASKAP["NAME"] + " (Stokes" +stokesVersion.upper() + ")")


    #Saves the image into its own folder to be used along with the other images


    xmin, xmax, ymin, ymax = plt.axis()
    xAverage = (xmin + xmax)/2
    yAverage = (ymin + ymax)/2


    if stokesType == "components" or stokesType == "both":
        size =  math.sqrt( (((pulsarDataStokesV["ra_err"]  * math.sin(pulsarDataStokesV["dec_err"] )) ** 2
                         + pulsarDataStokesV["dec_err"]  **2 ) ** 1/2)) * 500
    else:
        size = 20
    stokesV = plt.plot(pulsarDataStokesV["skycoord"].ra.deg, pulsarDataStokesV["skycoord"].dec.deg,
                             color="magenta",marker="x", markersize=size, transform=ax.get_transform("fk5"))

    ATNF = plt.axvline(x=xAverage, color='red', ymin=0, ymax=1, linewidth = 1)
    plt.axhline(y=yAverage, color='red', xmin=0, xmax=1, linewidth = 1)

    if stokesType == "components" or stokesType == "both":
        stokesIError = [pulsarDataStokesI["ra_err"] if pulsarDataStokesI["ra_err"] != 0 else 1,
                    pulsarDataStokesI["dec_err"] if pulsarDataStokesI["dec_err"] != 0 else 1]
    else:
        stokesIError = [1,1]
    stokesI = Ellipse((pulsarDataStokesI["skycoord"].ra.deg, pulsarDataStokesI["skycoord"].dec.deg),
                          (math.sqrt(stokesIError[0]) / 3600 * 500 * math.sqrt(radius)),
                          (math.sqrt(stokesIError[1]) / 3600 * 500 * math.sqrt(radius)), linewidth=2.5, edgecolor="orange", facecolor="none", transform=ax.get_transform("fk5"))
    ax.add_artist(stokesI)

    if secondMarkers:

        stokesIislands = Ellipse((pulsarDataStokesI["skycoord"].ra.deg, pulsarDataStokesI["skycoord"].dec.deg),
                          500 / 3600 * math.sqrt(radius),
                          500 / 3600 * math.sqrt(radius), linewidth=2.5,
                          edgecolor="green", facecolor="none", transform=ax.get_transform("fk5"))
        stokesVislands = plt.plot(pulsarDataStokesV["skycoord"].ra.deg, pulsarDataStokesV["skycoord"].dec.deg,
                           color="blue", marker="x", markersize=size, transform=ax.get_transform("fk5"))
        ax.add_artist(stokesIislands)

        ax.legend((ATNF, stokesV[0], stokesI ,stokesVislands[0], stokesIislands),
                  ("ATNF", "StokesVcomponents", "StokesIcomponents", "StokesIislands", "StokesIislands"))

    else:
        ax.legend((ATNF, stokesV[0], stokesI),("ATNF", "StokesV", "StokesI"))


    plt.savefig('pulsars/' + pulsarDataASKAP["NAME"] + "/cutout" + str(radius) + stokesVersion +"_" + stokesType+ ".png")

    plt.show()
    plt.close()

def create_cutout(fits_filePath, title, radius, markers, racsType, position=None, legend = False, stokestype="", pulsarFilePath="pulsars"):
    # Creates a cutout of a specific racs mosiac based on positions of a pulsar in ATNF
    # and adds position markers for its position in different catalogs

    # Opens the FITS image and make a World Coordinate System object
    try:
        with fits.open(fits_filePath) as hdul:
            data, header = hdul[0].data, hdul[0].header
            wcs = WCS(header, naxis=2)
            if position == None:
                position = SkyCoord(ra=header["CRVAL1"], dec=header["CRVAL2"], unit=(u.deg, u.deg))

    except Exception as e:
        print(e, fits_filePath)
        return
    # Make a postage stamp cutout at the area of interest
    # and change the data units from Jy to mJy
    cutout = Cutout2D(data, position, radius * u.deg, wcs=wcs)
    data = cutout.data * 1000
    wcs = cutout.wcs

    # Make a figure using the cutout WCS
    fig = plt.figure()
    ax = fig.add_subplot(111, projection=wcs)

    # Stretch the data for better contrast and plot the image
    norm = ImageNormalize(data, interval=ZScaleInterval(contrast=0.2))
    im = ax.imshow(data, cmap=plt.cm.gray_r, norm=norm)
    fig.colorbar(im, label=r"Flux Density (mJy beam$^{-1}$)")
    ax.set_ylabel("Dec (J2000)")
    ax.set_xlabel("RA (J2000)")
    ax.set_title(title)
    ylim = plt.ylim()
    xlim = plt.xlim()


    # creates a marker for each position and if legend is set creates a legend from them`
    markers_names = []
    markers_order = []
    name = ""
    first = {"II": True, "IC": True, "VI": True, "VC": True, "A": True}
    for i in markers:
        if i["type"] == "RACS_StokesI":
            artist = stokesIMarker(i, radius, ax)
            ax.add_artist(artist)
            test = testMarker(i, radius, ax)
            try:
                ax.add_artist(test)
            except Exception as e:
                print(e)
            if legend == "total":
                markers_order.append(artist)
                markers_names.append(i["type"])
            elif legend == "first":
                if i["stokesType"] == "components":
                    if first["IC"]:
                        first["IC"] = False
                        markers_names.append(i["stokesType"] + "-" + i["type"])
                        markers_order.append(artist)
                else:
                    if first["II"]:
                        first["II"] = False
                        markers_names.append(i["stokesType"] + "-" + i["type"])
                        markers_order.append(artist)


        elif i["type"] == "RACS_StokesV":
            artist = stokesVMarker(i, radius, ax)
            ax.add_artist(artist)
            if legend == "total":
                markers_order.append(artist)
                markers_names.append(i["type"])
            elif legend == "first":
                if i["stokesType"] == "components":
                    if first["VC"]:
                        first["VC"] = False
                        markers_names.append(i["stokesType"] + "-" + i["type"])
                        markers_order.append(artist)
                else:
                    if first["VI"]:
                        first["VI"] = False
                        markers_names.append(i["stokesType"] + "-" + i["type"])
                        markers_order.append(artist)

        elif i["type"] == "ATNF":
            artist = atnfMarker(i, ax)
            ax.add_artist(artist)
            if legend == "total":
                markers_order.append(artist)
                markers_names.append(i["type"])
            elif legend == "first":
                if first["A"]:
                    first["A"] = False
                    markers_names.append(i["stokesType"] + "-" + i["type"])
                    markers_order.append(artist)

            name = i["NAME"]

        else:
            print(i["type"])



    if len(markers_names) != 0:
        ax.legend(markers_order, markers_names, loc=4)
    print(stokestype, "STOKES TYPE")
    plt.ylim(ylim[0],ylim[1])
    plt.xlim(xlim[0],ylim[1])
    plt.savefig(pulsarFilePath + "/" + name + '/cutout' + str(radius) + "," + racsType + "," + stokestype + ".png")
    plt.close()



def stokesIMarker(marker, radius, ax):
    # Creates an Ellipse marker and sizes based on zoom and error in position if component
    if marker["stokesType"] == "components":
        stokesIError = [marker["ra_err"] if marker["ra_err"] != 0 else 0,
                        marker["dec_err"] if marker["dec_err"] != 0 else 0]

    else:
        stokesIError = [1, 1]
    if marker["stokesType"] == "components":
        color = "orange"
    else:
        color = "green"
    #stokesI = Ellipse((marker["skycoord"].ra.deg, marker["skycoord"].dec.deg),
     #                 (math.sqrt(stokesIError[0]) / 3600 * 2000 * math.sqrt(radius)),
      #                (math.sqrt(stokesIError[1]) / 3600 * 2000 * math.sqrt(radius)) , linewidth=2.5, edgecolor=color,
       #               facecolor="none", transform=ax.get_transform("fk5"), alpha=0.5)
    stokesI = Ellipse((marker["skycoord"].ra.deg, marker["skycoord"].dec.deg), (stokesIError[0] * u.arcsec).to(u.deg).value,
                      (stokesIError[1] * u.arcsec).to(u.deg).value,
                linewidth=2.5, edgecolor=color, facecolor="none", transform=ax.get_transform("fk5"), alpha=0.5)

    return stokesI


def testMarker(marker, radius, ax):
    # Creates an Ellipse marker and sizes based on zoom and error in position if component
    if marker["stokesType"] == "components":
        stokesIError = [marker["ra_err"] if marker["ra_err"] != 0 else 0,
                        marker["dec_err"] if marker["dec_err"] != 0 else 0]

    color = "blue"

    #stokesI = Ellipse((marker["skycoord"].ra.deg, marker["skycoord"].dec.deg), stokesIError[0], stokesIError[1],
     #                 linewidth=2.5, edgecolor=color, facecolor="none", transform=ax.get_transform("fk5"), alpha=0.5)
    try:
        stokesI = Polygon(polygon_coords([marker["skycoord"].ra,marker["skycoord"].dec], stokesIError[0] * u.arcsec, stokesIError[1] * u.arcsec, 10),
                          color=color, transform=ax.get_transform("fk5"))
    except Exception as e:
        print(e)

    return stokesI

def stokesVMarker(marker, radius, ax):
    # Maybe use a hexagon here instead? Can make it close to an ellipse but different enough that it is visually distinct
    if marker["stokesType"] == "components":
        size = math.sqrt((((marker["ra_err"] * math.sin(marker["dec_err"])) ** 2
                                   + marker["dec_err"] ** 2) ** 1 / 2) * radius) * 800
        color = "magenta"
        markerType = "x"
    else:
        size = 15
        color = "blue"
        markerType = "o"

    stokesV = plt.plot(marker["skycoord"].ra.deg, marker["skycoord"].dec.deg, color=color, marker=markerType, markersize=size, transform=ax.get_transform("fk5"), alpha=0.5)[0]
    return stokesV

def atnfMarker(marker, ax):

    ATNF = plt.plot(marker["skycoord"].ra.deg, marker["skycoord"].dec.deg, color="red", marker="+", markersize=25, transform=ax.get_transform("fk5"), alpha=0.5)[0]

    return ATNF


def get_rms(pulsarData, mosaicFilePaths, pulsarFilePaths, stokesType, size=20, creating_upper_limits = True):
    # loads image and creates sub-image
    # takes sub-image and finds standard deviation
    # it then looks at the standard deviation from the inter-quartile range
    sigma_stds = np.zeros(len(pulsarData["ASKAP"]))
    sigma_iqrs = np.zeros(len(pulsarData["ASKAP"]))

    for i in range(len(pulsarData["ASKAP"])):
        try:
            position = pulsarData["ASKAP"]["skycoord"][i]
            if stokesType == "racsi":
                file_name = mosaicFilePaths[1] + "RACS_" + pulsarData["racsi"]["fileName"][i] + ".EPOCH00.I.fits"
            elif stokesType == "racsv":
                file_name = mosaicFilePaths[0] + "RACS_test4_1.05_" + pulsarData["racsi"]["fileName"][i] + ".fits"
            with fits.open(file_name) as hdul:

                data, header = hdul[0].data, hdul[0].header
                wcs = WCS(header, naxis=2)
            cutout = Cutout2D(data, position, size, wcs=wcs)
            data = cutout.data * 1000
            std = data.std()
            plt.clf()
            data_array = np.zeros(size ** 2)
            for j in range(len(data)):
                for k in range(len(data[j])):
                    data_array[size * j + k] = data[j][k]
            plt.hist(data_array, bins=40, density=True, stacked=True)
            x_vals = np.arange(np.amin(data_array), np.amax(data_array), 0.00001)
            y_vals = scipy.stats.norm.pdf(x_vals,np.mean(data_array), std)
            plt.plot(x_vals, y_vals)
            #print(x_vals)
            #print(y_vals)
            plt.xlabel(r'Flux Density (mJy beam$^{-1}$)')
            plt.ylabel('Portion of Pixels Around x Flux Density')
            plt.title("Histogram of Flux Density Around the Pulsar {} (StokesV)".format(pulsarData["ASKAP"]["NAME"][i]))
            for j in [-2, -1, 1, 2]:
                plt.vlines(np.mean(data) + std * j, 0, scipy.stats.norm.pdf([np.mean(data) + std * j], np.mean(data_array), std)[0], color="black")

            plt.savefig(pulsarFilePaths + "/" + pulsarData["ASKAP"]["NAME"][i] + "/std_" + stokesType + ".png")
            print("did STD calculation: " + str(i) + "/" + str(len(pulsarData["ASKAP"])))
            plt.close()
            sigma_iqr = scipy.stats.iqr(data) / 1.349
            print(std, sigma_iqr)
            sigma_stds[i] = std
            sigma_iqrs[i] = sigma_iqr

            try:
                file = open(pulsarFilePaths + "/" + pulsarData["ASKAP"]["NAME"][i] + "/std_info_" + stokesType, "x")
            except Exception as e:
                file = open(str(pulsarFilePaths) + "/" + str(pulsarData["ASKAP"]["NAME"][i]) + "/std_info_" + stokesType, "w")

            info = """dataValues: {}
            sigma: {}
            iqr_sigma: {}""".format(np.array2string(data), str(std), str(sigma_iqr))

            file.write(info)
            file.close()
        except Exception as e:
            print(e, i)
    pulsarData[stokesType].add_column(Column(sigma_stds, name="sigmaOfMosaic"))
    pulsarData[stokesType].add_column(Column(sigma_iqrs, name="sigmaOfMosaicIQR"))

    if creating_upper_limits == True:
        RACS_has_upper_limit = np.array(np.zeros(len(pulsarData["ASKAP"])), dtype=bool)
        for i in range(len(pulsarData[stokesType + "_match"])):
            RACS_has_upper_limit[i] = False if pulsarData[stokesType + "_match"][i] else True
            if RACS_has_upper_limit[i]:
                pulsarData[stokesType]["flux_peak"][i] = 3 * sigma_iqrs[i]
                pulsarData[stokesType]["flux_int"][i] = 3 * sigma_iqrs[i]
        pulsarData[stokesType].add_column(Column(RACS_has_upper_limit, name="hasUpperLimit"))




    #pulsarData[stokesType].add_column(Column(sigma_iqrs, name="detectionSigmas"))


    return pulsarData









def createImage(images, dimensions, imageSize):
    #combines the images into a 2x2 grid adding one at a time.
    #return addImage(1,1,addImage(0,1,addImage(1, 0, images[0], images[1]), images[2]), images[3])

    imageMain = addWhiteSpace(dimensions, imageSize)
    for i in range(dimensions[0]):
        for j in range(dimensions[1]):
            imageMain = addImagePart(i, j, imageMain, images[i + j * dimensions[0]], imageSize)
    return imageMain



def addImagePart(x, y, imageMain, imageToAdd, imageSize):
    imageMain.paste(imageToAdd, (int(imageSize[0] * x), int(imageSize[1] * y)))
    return imageMain


def addWhiteSpace(dimensions, imageSize):
    return Image.new('RGB', (imageSize[0] * dimensions[0], imageSize[1] * dimensions[1]))




def addImage(x, y, image1, image2):
    #adds images together horizontally, vertically, or both
    if x == 1 and y == 0:
        dst = Image.new('RGB', (image1.width + image2.width, image1.height))
        dst.paste(image1, (0, 0))
        dst.paste(image2, (image1.width, 0))
    elif y == 1 and x == 0:
        dst = Image.new('RGB', (image1.width, image1.height + image2.height))
        dst.paste(image1, (0, 0))
        dst.paste(image2, (0, image2.height))
    else:
        dst = Image.new('RGB', (image1.width, image1.height))
        dst.paste(image1, (0, 0))
        dst.paste(image2, (int(image1.width/2), int(image1.height/2)))
    return dst



def polygon_coords(center, error_RA, error_DEC, sides):
    angles = list(i * 2 * math.pi / sides for i in range(sides))
    coords = np.zeros((sides, 2)) * u.deg
    for i in range(sides):
        coords[i][0] = center[0] + np.cos(angles[i]) * error_RA
        coords[i][1] = center[1] + np.sin(angles[i]) * error_DEC

    return coords
