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

# ----------- generateCoordsBetweenEpochs(): -----------
# [orbit] is an sbpy Orbit type (NOTE: This function currently only takes one orbit at a time!)
# [epoch_start] & [epoch_end] are astropy Time objects (can be in any scale/format, i.e. utc, mjd)
# [step] is the time between measurements, in days
# [printALot] will print more information if set to True
# [output_type] will return a Pandas dataframe if "false", a .csv file if set to "csv", and a .parquet file if set to "parquet"
# [output_name] is the name of the file that is outputted. If output_type="false" you can set this to "" (or whatever you want really)
# This function will return 0 if any known errors occur
def generateCoordsBetweenEpochs(orbit, epoch_start, epoch_end, step, printALot, output_name, output_type):
    epoch_start_mjd = epoch_start.mjd   # Float version of epoch_start
    epoch_end_mjd = epoch_end.mjd       # Float version of epoch_end
    output = pd.DataFrame(columns=["Name", "RA", "DEC", "Epoch"])
    name = orbit['targetname'][0] if len(orbit) > 0 else None

    if(epoch_end_mjd < epoch_start_mjd):
        print(f"getCoordsBetweenEpochs(): Error! The starting epoch ({epoch_start_mjd}) is further in time than the ending epoch ({epoch_end_mjd}). Try flipping them")
        return 0
    
    start = time.time()
    epoch_current_mjd = epoch_start_mjd
    while(epoch_current_mjd < epoch_end_mjd):
        epoch_current = Time(epoch_current_mjd, format='mjd')
        coords = Ephem.from_oo(orbit, epoch_current)
        RA = coords['RA'][0].value
        DEC = coords['DEC'][0].value
        output.loc[len(output)] = {
            "Name": name,
            "RA": RA,
            "DEC": DEC,
            "Epoch": epoch_current.mjd
        }
        # if(printALot): print(f"RA: {RA}   DEC: {DEC}     Epoch: {epoch_current}")
        epoch_current_mjd = epoch_current_mjd + step
    end = time.time()

    if(printALot): print(output)
    print(f"generateCoordsBetweenEpochs(): DONE! Took {end-start} seconds")

    if(output_type == "csv"):
        output_name = output_name + ".csv"
        output.to_csv(output_name)
        if(printALot): print(f"Outputting to {output_name}")
    elif(output_type == "parquet"):
        output_name = output_name + ".parquet"
        output.to_parquet(output_name)
        if(printALot): print(f"Outputting to {output_name}")
    elif(output_type == "false"):
        return output
    else:
        if(printALot): print(f"generateCoordsBetweenEpochs(): Error! The input variable (output_type) is an unknown value: {output_type}. Known values are 'csv', 'parquet', or 'false'")
        return 0

# This is an example of how to run the function
epoch_init = Time('2026-05-05T00:00:00', scale='utc')
epoch_final = Time('2031-05-05T00:00:00', scale='utc')
ceres = orbitFromDict('Pallas', 2.7701937, 0.2305404, 34.92402, 310.91037, 172.8953, 168.79869, Time(60808, format='mjd'), 4.11, 0.15) # THIS IS PALLAS
# generateCoordsBetweenEpochs(ceres, epoch_init, epoch_final, 0.1, True, "test2", "parquet")

# Should I do objectS or just one object at a time?
# [orbit] is the orbit of the object (AstroPy Orbit object)
# [csvFile] is the name of the CSV file that contains the generated coordinate (use generateCoordsBetweenEpochs() to generate new ones) (String)
# [startEpoch] and [endEpoch] are the start and end dates for the search (sbpy Time object)
# [R] is how close the asteroid needs to be to the center of the LSST observation to "count" (float, in degrees). Default is 1.75
# [d] is how close in time the LSST observation and the asteroid's coords need to be. Recommended: 
# [step] is the step count in days (float)
def listWhenObjectIsInLSST(orbit, file, startEpoch, endEpoch, R, d, printALot, output_name, output_type):
    startEpoch_mjd = startEpoch.mjd
    endEpoch_mjd = endEpoch.mjd
    parquet_pattern = re.compile(r'.*\.parquet$', re.IGNORECASE)
    csv_pattern = re.compile(r'.*\.csv$', re.IGNORECASE)
    output = pd.DataFrame(columns=["Name", "RA", "DEC", "Epoch"])

    # Setting up the dataframes
    lsst_obs = pd.read_parquet('lsst_snapshots.parquet',dtype_backend="pyarrow",columns=["fieldRA","fieldDec","observationStartMJD"],filters=[('observationStartMJD','>=',startEpoch_mjd),('observationStartMJD','<=',endEpoch_mjd)])
    if(parquet_pattern.search(file)):
        # if(printALot): print(f"listWhenObjectIsInLSST(): Input file loaded as a .parquet ({file})")
        object = pd.read_parquet(file,columns=["Name","RA","DEC","Epoch"])
    elif(csv_pattern.search(file)):
        # if(printALot): print(f"listWhenObjectIsInLSST(): Input file loaded as a .csv ({file})")
        object = pd.read_csv(file)
    else:
        print("listWhenObjectIsInLSST(): Error! Inputted file ({file}) is not a .csv or a .parquet")
        return 0

    currentRA = 0
    currentDEC = 0
    currentRA_LSST = 0
    currentDEC_LSST = 0
    start = time.time()
    for i in range(len(object)):
        currentRA = object['RA'][i]
        currentDEC = object['DEC'][i]
        currentEpoch = object['Epoch'][i]
        upperBound = currentEpoch + d/2
        lowerBound = currentEpoch - d/2
        # print(f"upperBound ({upperBound})    lowerBound ({lowerBound})    currentEpoch ({currentEpoch})    currentRA ({currentRA})    currentDEC ({currentDEC})")
        # 
        lsst_obs = pd.read_parquet('lsst_snapshots.parquet',dtype_backend="pyarrow",columns=["fieldRA","fieldDec","observationStartMJD"],filters=[('observationStartMJD','>=',lowerBound),('observationStartMJD','<=',upperBound)])
        # ^ THIS WILL NEED TO BE UPDATED TO WORK WITH CSVS AS WELL
        for ii in range(len(lsst_obs)):
            currentEpoch_LSST = lsst_obs['observationStartMJD'][ii]
            if(upperBound < currentEpoch_LSST < lowerBound):
                currentRA_LSST = lsst_obs['fieldRA'][ii]
                currentDEC_LSST = lsst_obs['fieldDec'][ii]
                P = math.sqrt(((currentRA - currentRA_LSST)* np.cos(currentDEC-currentDEC_LSST))**2 + (currentDEC-currentDEC_LSST)**2)
                if P <= R:
                    print(f"Asteroid will be in view of the LSST at MJD of {currentEpoch_LSST}")
        # print(f"currentRA = {currentRA}     |      currentDEC = {currentDEC}")
    end = time.time()
    print("DONE! Time:",(end-start),"seconds")
    return
listWhenObjectIsInLSST(ceres, "test.parquet", epoch_init, epoch_final, 1.75, 1, True, "pallas.csv", True)


# This function currently only works with asteroids that are in the "generated_coordinates.parquet" file.
#   It does not have the capability to generate new orbit data:
#       targetname, a, e, i, Node, Peri, M, G, H, epoch_og
#   A new function must be created to generate orbit data
def listWhenAstIsInLSST(astname, afterEpoch, d):
    returnedValues = pd.DataFrame(columns=["ssnamenr", "RA", "DEC", "Epoch", "RA_LSST", "DEC_LSST", "Epoch_LSST"])
    start = time.time()
    afterEpoch = afterEpoch.mjd
    R = 1.75 # Degrees
    ast_df = pd.read_parquet('generated_coordinates.parquet',dtype_backend="pyarrow",columns=["ssnamenr","RA","DEC","observationMJD"],filters=[('ssnamenr', '=', astname)])
    end = time.time()
    print("Initialization took ",(end-start),"seconds")

    if(len(ast_df) == 0):
        print("ERROR LOADING ASTEROID'S DATA! SSNAMENR OF",astname,"NOT IN DATABASE")
    else:
        in_or_not = np.zeros(len(ast_df), dtype=bool)
        currentRA = 0
        currentDEC = 0
        currentRA_LSST = 0
        currentDEC_LSST = 0

        # Repeats over the length of the asteroid's RA/DEC list (800 right now)
        for i in range(len(ast_df)): # len(ast_df)
            currentEpoch = ast_df['observationMJD'][i]
            if(currentEpoch > afterEpoch): # Makes sure the observation MJD is not in the past
                currentRA = ast_df['RA'][i]
                currentDEC = ast_df['DEC'][i]
                upperBound = currentEpoch + d/2
                lowerBound = currentEpoch - d/2
                LSST_snapshots = pd.read_parquet('lsst_snapshots.parquet',dtype_backend="pyarrow",columns=["fieldRA","fieldDec","observationStartMJD"],filters=[('observationStartMJD','>=',lowerBound),('observationStartMJD','<=',upperBound)]) # Is this really the fastest way to do this? Calling the parquet file every time?
                for ii in range(len(LSST_snapshots)):
                    currentRA_LSST = LSST_snapshots['fieldRA'][ii]
                    currentDEC_LSST = LSST_snapshots['fieldDec'][ii]
                    P = math.sqrt(((currentRA - currentRA_LSST)* np.cos(currentDEC-currentDEC_LSST))**2 + (currentDEC-currentDEC_LSST)**2)
                    if P <= R:
                        currentEpoch_LSST = LSST_snapshots['observationStartMJD'][ii]
                        returnedValues.loc[len(returnedValues)] = [str(astname), currentRA, currentDEC, currentEpoch, currentRA_LSST, currentDEC_LSST, currentEpoch_LSST]
                        print("Asteroid #",astname,"will be in view of the LSST at MJD of",currentEpoch,"(asteroid time)")
            else:
                print(currentEpoch,"is in the past: skipping")
    end = time.time()
    print("DONE! Time:",(end-start),"seconds")
    return returnedValues

# POTENTIAL PROBLEM: THE ORBITAL TYPE IS 'MB II' FOR PALLAS
# ceres = orbitFromDict('Ceres', 'KEP', 2.77, 0.0786, 10.6, 73.6, 80.3, 188.7, epoch_calc, 3.34, 0.15) REAL CERES

# def getCoordsAtEpoch(orbit, epoch):
#     epoch = epoch + np.arange(1) * u.day
#     coords = Ephem.from_oo(orbit, epoch)
#     # print(coords)
#     o_RA = coords['RA'].value
#     o_DEC = coords['DEC'].value
#     return {"RA":o_RA, "DEC":o_DEC};