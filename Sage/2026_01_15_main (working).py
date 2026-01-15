import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import spiceypy
import pyoorb as oorb
import pyarrow
import pandas as pd
import astropy
import astropy.units as u
import sbpy
import pymongo
import time
import re
from astropy.time import Time
from sbpy.data import Orbit, Ephem
import warnings
warnings.filterwarnings('ignore')

# w = argument of periapsis (I think)
# Omega = node
def orbitFromDict(targetname, a, e, i, w, omega, M, epoch, H, G):
    return Orbit.from_dict({''
    'targetname': targetname,
    'a': a * u.au,
    'e': e,
    'i': i * u.deg,
    'w': w * u.deg,
    'Omega': omega * u.deg,
    'M': M * u.deg,
    'epoch': epoch,
    'H': H * u.mag,
    'G': G});

# ----------- LSST_MAIN(): -----------
# [orbit] is an sbpy Orbit type (NOTE: This function currently only takes one orbit at a time!)
# [epoch_start] & [epoch_end] are astropy Time objects (can be in any scale/format, i.e. utc, mjd)
# [step] is the time between measurements, in days
# [printALot] will print more information if set to True
# [output_type] will return a Pandas dataframe if "false", a .csv file if set to "csv", and a .parquet file if set to "parquet"
# [output_name] is the name of the file that is outputted. If output_type="false" you can set this to "" (or whatever you want really)
# This function will return 0 if any known errors occur
def LSST_MAIN(orbit, epoch_start, epoch_end, step, printALot, output_name, output_type, d, R):
    epoch_start_mjd = epoch_start.mjd   # Float version of epoch_start
    epoch_end_mjd = epoch_end.mjd       # Float version of epoch_end
    object_coords = pd.DataFrame(columns=["Name", "RA", "DEC", "Epoch"])
    name = orbit['targetname'][0] if len(orbit) > 0 else None
    print(f"epochStart: {epoch_start_mjd}       epochEnd: {epoch_end_mjd}")
    print(f"R: {R}      d: {d}      step: {step}")

    if(epoch_end_mjd < epoch_start_mjd):
        print(f"getCoordsBetweenEpochs(): Error! The starting epoch ({epoch_start_mjd}) is further in time than the ending epoch ({epoch_end_mjd}). Try flipping them")
        return 0
    
    start = time.time()
    epoch_current_mjd = epoch_start_mjd
    while(epoch_current_mjd <= epoch_end_mjd):
        epoch_current = Time(epoch_current_mjd, format='mjd')
        coords = Ephem.from_oo(orbit, epoch_current)
        RA = coords['RA'][0].value
        DEC = coords['DEC'][0].value
        object_coords.loc[len(object_coords)] = {
            "Name": name,
            "RA": RA,
            "DEC": DEC,
            "Epoch": epoch_current.mjd
        }
        epoch_current_mjd = epoch_current_mjd + step
    end = time.time()
    # print(f"Coordinates: {object_coords}")
    print(f"LSST_MAIN(): Coordinates have been generated for the object! Took {end-start} seconds")

    # PART 2 -----------------------------------------------------------------------------------------------------------------------------------------------------------
    currentRA = 0
    currentDEC = 0
    currentRA_LSST = 0
    currentDEC_LSST = 0
    start = time.time()
    output = pd.DataFrame(columns=["Name", "Epoch_OBJ", "Epoch_LSST", "RA_OBJ", "DEC_OBJ", "RA_LSST", "DEC_LSST", "P", "R"])

    for i in range(len(object_coords)):
        currentRA = object_coords['RA'][i]
        currentDEC = object_coords['DEC'][i]
        currentEpoch = object_coords['Epoch'][i]
        upperBound = currentEpoch + d/2
        lowerBound = currentEpoch - d/2
        lsst_obs = pd.read_parquet('lsst_snapshots.parquet',dtype_backend="pyarrow",columns=["fieldRA","fieldDec","observationStartMJD"],filters=[('observationStartMJD','>=',lowerBound),('observationStartMJD','<=',upperBound)])

        for ii in range(len(lsst_obs)):
            currentEpoch_LSST = lsst_obs['observationStartMJD'][ii]
            if(upperBound > currentEpoch_LSST > lowerBound):
                currentRA_LSST = lsst_obs['fieldRA'][ii]
                currentDEC_LSST = lsst_obs['fieldDec'][ii]

                currentRA_rad = np.radians(currentRA)
                currentDEC_rad = np.radians(currentDEC)
                currentRA_LSST_rad = np.radians(currentRA_LSST)
                currentDEC_LSST_rad = np.radians(currentDEC_LSST)

                P = np.arccos(np.sin(currentDEC_rad)*np.sin(currentDEC_LSST_rad) + np.cos(currentDEC_rad)*np.cos(currentDEC_LSST_rad)*np.cos(currentRA_rad - currentRA_LSST_rad))
                P = np.degrees(P)
                if P <= R:
                    output.loc[len(output)] = {
                        "Name": name,
                        "Epoch_OBJ": currentEpoch,
                        "Epoch_LSST": currentEpoch_LSST,
                        "RA_OBJ": currentRA,
                        "DEC_OBJ": currentDEC,
                        "RA_LSST": currentRA_LSST,
                        "DEC_LSST": currentDEC_LSST,
                        "P": P,
                        "R": R
                    }
                    # print(f"Asteroid will be in view of the LSST at MJD of {currentEpoch_LSST}")
        # print(f"currentRA = {currentRA}     |      currentDEC = {currentDEC}")
    print(output)
    # output.to_csv("testing_output.csv")
    end = time.time()
    print("DONE! Time:",(end-start),"seconds")
# Potential next steps:
# - Make a GUI
# - Optimize the program
# - Add more graphing capabilities
# - Start thinking about re-implementing it into MONSOON (do we even need MONSOON...?)
# - GET IT WORKING ON MY LAPTOP
# - Get the gang all set up locally
# - Reinstate our MONSOON access
# - TALK TO TRILLING

epoch_init = Time('2026-05-05T00:00:00', scale='utc')
epoch_final = Time('2028-06-05T00:00:00', scale='utc')
ceres = orbitFromDict('Pallas', 2.76992582511479, 0.2306429787781384, 34.92832687077855, 310.9333840114307, 172.8885963367437, 211.5297778033731, Time('2025-11-21T00:00:00', scale='utc'), 4.11, 0.11) # THIS IS PALLAS
LSST_MAIN(ceres, epoch_init, epoch_final, 1, True, "test2", "parquet", 1, 1.75)




























# Should we do objectS or just one object at a time?
# [orbit] is the orbit of the object (AstroPy Orbit object)
# [csvFile] is the name of the CSV file that contains the generated coordinate (use generateCoordsBetweenEpochs() to generate new ones) (String)
# [startEpoch] and [endEpoch] are the start and end dates for the search (sbpy Time object)
# [R] is how close the asteroid needs to be to the center of the LSST observation to "count" (float, in degrees). Default is 1.75
# [d] is how close in time the LSST observation and the asteroid's coords need to be. Recommended: 
# [step] is the step count in days (float)

# def listWhenObjectIsInLSST(orbit, file, startEpoch, endEpoch, R, d, printALot, output_name, output_type):
#     startEpoch_mjd = startEpoch.mjd
#     endEpoch_mjd = endEpoch.mjd
#     parquet_pattern = re.compile(r'.*\.parquet$', re.IGNORECASE)
#     csv_pattern = re.compile(r'.*\.csv$', re.IGNORECASE)
#     output = pd.DataFrame(columns=["Name", "RA", "DEC", "Epoch"])

#     # Setting up the dataframes
#     lsst_obs = pd.read_parquet('lsst_snapshots.parquet',dtype_backend="pyarrow",columns=["fieldRA","fieldDec","observationStartMJD"],filters=[('observationStartMJD','>=',startEpoch_mjd),('observationStartMJD','<=',endEpoch_mjd)])
#     if(parquet_pattern.search(file)):
#         # if(printALot): print(f"listWhenObjectIsInLSST(): Input file loaded as a .parquet ({file})")
#         object = pd.read_parquet(file,columns=["Name","RA","DEC","Epoch"])
#     elif(csv_pattern.search(file)):
#         # if(printALot): print(f"listWhenObjectIsInLSST(): Input file loaded as a .csv ({file})")
#         object = pd.read_csv(file)
#     else:
#         print("listWhenObjectIsInLSST(): Error! Inputted file ({file}) is not a .csv or a .parquet")
#         return 0

#     currentRA = 0
#     currentDEC = 0
#     currentRA_LSST = 0
#     currentDEC_LSST = 0
#     start = time.time()
#     for i in range(len(object)):
#         currentRA = object['RA'][i]
#         currentDEC = object['DEC'][i]
#         currentEpoch = object['Epoch'][i]
#         upperBound = currentEpoch + d/2
#         lowerBound = currentEpoch - d/2
#         # print(f"upperBound ({upperBound})    lowerBound ({lowerBound})    currentEpoch ({currentEpoch})    currentRA ({currentRA})    currentDEC ({currentDEC})")
#         # 
#         lsst_obs = pd.read_parquet('lsst_snapshots.parquet',dtype_backend="pyarrow",columns=["fieldRA","fieldDec","observationStartMJD"],filters=[('observationStartMJD','>=',lowerBound),('observationStartMJD','<=',upperBound)])
#         # ^ THIS WILL NEED TO BE UPDATED TO WORK WITH CSVS AS WELL
#         for ii in range(len(lsst_obs)):
#             currentEpoch_LSST = lsst_obs['observationStartMJD'][ii]
#             if(upperBound < currentEpoch_LSST < lowerBound):
#                 currentRA_LSST = lsst_obs['fieldRA'][ii]
#                 currentDEC_LSST = lsst_obs['fieldDec'][ii]
#                 P = math.sqrt(((currentRA - currentRA_LSST)* np.cos(currentDEC-currentDEC_LSST))**2 + (currentDEC-currentDEC_LSST)**2)
#                 if P <= R:
#                     print(f"Asteroid will be in view of the LSST at MJD of {currentEpoch_LSST}")
#         # print(f"currentRA = {currentRA}     |      currentDEC = {currentDEC}")
#     end = time.time()
#     print("DONE! Time:",(end-start),"seconds")
#     return
# listWhenObjectIsInLSST(ceres, "test.parquet", epoch_init, epoch_final, 1.75, 1, True, "pallas.csv", True)







# This function currently only works with asteroids that are in the "generated_coordinates.parquet" file.
#   It does not have the capability to generate new orbit data:
#       targetname, a, e, i, Node, Peri, M, G, H, epoch_og
#   A new function must be created to generate orbit data

# def listWhenAstIsInLSST(astname, afterEpoch, d):
#     returnedValues = pd.DataFrame(columns=["ssnamenr", "RA", "DEC", "Epoch", "RA_LSST", "DEC_LSST", "Epoch_LSST"])
#     start = time.time()
#     afterEpoch = afterEpoch.mjd
#     R = 1.75 # Degrees
#     ast_df = pd.read_parquet('generated_coordinates.parquet',dtype_backend="pyarrow",columns=["ssnamenr","RA","DEC","observationMJD"],filters=[('ssnamenr', '=', astname)])
#     end = time.time()
#     print("Initialization took ",(end-start),"seconds")

#     if(len(ast_df) == 0):
#         print("ERROR LOADING ASTEROID'S DATA! SSNAMENR OF",astname,"NOT IN DATABASE")
#     else:
#         in_or_not = np.zeros(len(ast_df), dtype=bool)
#         currentRA = 0
#         currentDEC = 0
#         currentRA_LSST = 0
#         currentDEC_LSST = 0

#         # Repeats over the length of the asteroid's RA/DEC list (800 right now)
#         for i in range(len(ast_df)): # len(ast_df)
#             currentEpoch = ast_df['observationMJD'][i]
#             if(currentEpoch > afterEpoch): # Makes sure the observation MJD is not in the past
#                 currentRA = ast_df['RA'][i]
#                 currentDEC = ast_df['DEC'][i]
#                 upperBound = currentEpoch + d/2
#                 lowerBound = currentEpoch - d/2
#                 LSST_snapshots = pd.read_parquet('lsst_snapshots.parquet',dtype_backend="pyarrow",columns=["fieldRA","fieldDec","observationStartMJD"],filters=[('observationStartMJD','>=',lowerBound),('observationStartMJD','<=',upperBound)]) # Is this really the fastest way to do this? Calling the parquet file every time?
#                 for ii in range(len(LSST_snapshots)):
#                     currentRA_LSST = LSST_snapshots['fieldRA'][ii]
#                     currentDEC_LSST = LSST_snapshots['fieldDec'][ii]
#                     P = math.sqrt(((currentRA - currentRA_LSST)* np.cos(currentDEC-currentDEC_LSST))**2 + (currentDEC-currentDEC_LSST)**2)
#                     if P <= R:
#                         currentEpoch_LSST = LSST_snapshots['observationStartMJD'][ii]
#                         returnedValues.loc[len(returnedValues)] = [str(astname), currentRA, currentDEC, currentEpoch, currentRA_LSST, currentDEC_LSST, currentEpoch_LSST]
#                         print("Asteroid #",astname,"will be in view of the LSST at MJD of",currentEpoch,"(asteroid time)")
#             else:
#                 print(currentEpoch,"is in the past: skipping")
#     end = time.time()
#     print("DONE! Time:",(end-start),"seconds")
#     return returnedValues



# def generatePrettyPlotOfObject(parquetFile, plotName):
#     object_coords = pd.read_parquet(parquetFile)
#     RA_plot = np.radians(object_coords['RA'])
#     DEC_plot = np.radians(object_coords['DEC'])
#     RA_plot = np.remainder(RA_plot + np.pi, 2 * np.pi) - np.pi
#     print(RA_plot)
#     name = object_coords['Name'][0]

#     plt.figure(figsize=(10, 6))
#     plt.subplot(projection="aitoff")
#     plt.plot(RA_plot, DEC_plot, 'o', color='black', markersize=2, alpha=0.2)
#     plt.xlabel('RA')
#     plt.ylabel('DEC')
#     plt.title(name)
#     plt.grid(axis='x')
#     plt.grid(axis='y')
#     plt.savefig(plotName)
#     plt.grid(True)
#     plt.show()
# generatePrettyPlotOfObject("test2.parquet", "Hiiii")
