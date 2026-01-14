import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import spiceypy
import pyoorb as oorb
import pandas as pd
import astropy
import astropy.units as u
import sbpy
import pymongo
import time
from astropy.time import Time
from sbpy.data import Orbit, Ephem

import warnings
warnings.filterwarnings('ignore')

# from concurrent.features import ProcessPoolExecutor

# TO RUN ME:
# cd .conda/envs/find_orbits
# conda activate find_orbits
# python main.py

# This is the (pre-processed)
df = pd.read_csv('mpc_df.csv')

# w = argument of periapsis (I think)
# Omega = node
def orbitFromDict(targetname, a, e, i, w, omega, M, epoch_og, H, G):
    return Orbit.from_dict({''
    'targetname': targetname,
    'a': a * u.au,
    'e': e,
    'i': i * u.deg,
    'w': w * u.deg,
    'Omega': omega * u.deg,
    'M': M * u.deg,
    'epoch': epoch_og,
    'H': H * u.mag,
    'G': G});

def getCoordsAtEpoch(orbit, epoch):
    epoch = epoch + np.arange(1) * u.day
    coords = Ephem.from_oo(orbit, epoch)
    # print(coords)
    o_RA = coords['RA'].value
    o_DEC = coords['DEC'].value
    return {"RA":o_RA, "DEC":o_DEC};

epoch_init = Time('2026-05-05T00:00:00', scale='utc')
epoch_final = Time('2026-06-05T00:00:00', scale='utc')
ceres = orbitFromDict('Pallas', 2.7701937, 0.2305404, 34.92402, 310.91037, 172.8953, 168.79869, Time(60808, format='mjd'), 4.11, 0.15) # THIS IS PALLAS

# ----------- generateCoordsBetweenEpochs(): -----------
# [orbit] is an sbpy Orbit type
# [epoch_start] & [epoch_end] are astropy Time objects (can be in any scale/format, i.e. utc, mjd)
# [step] is the time between measurements, in days
# [printALot] will print more information if set to True
# [outputCSV] will return a Pandas dataframe if False, and output a CSV of [CSVname] if True (with no dataframe returned in this case)
# [CSVname] is the name of the CSV file that is outputted. If outputCSV=False you can set this to "" (or whatever you want really)
# This function will return 0 if any errors occur (but currently the only error that is checked for is if epoch_end < epoch_start)
def generateCoordsBetweenEpochs(orbit, epoch_start, epoch_end, step, printALot, CSVname, outputCSV):
    epoch_start_mjd = epoch_start.mjd   # Float version of epoch_start
    epoch_end_mjd = epoch_end.mjd       # Float version of epoch_end
    output = pd.DataFrame(columns=["Name", "RA", "DEC", "Epoch"])
    name = orbit['targetname'][0] if len(orbit) > 0 else None

    if(epoch_end_mjd < epoch_start_mjd):
        if(printALot): print(f"getCoordsBetweenEpochs(): Error! The starting epoch ({epoch_start_mjd}) is further in time than the ending epoch ({epoch_end_mjd}). Try flipping them")
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
            "Epoch": epoch_current
        }
        # if(printALot): print(f"RA: {RA}   DEC: {DEC}     Epoch: {epoch_current}")
        epoch_current_mjd = epoch_current_mjd + step
    
    end = time.time()
    if(printALot): print(output)
    print(f"generateCoordsBetweenEpochs(): DONE! Took {end-start} seconds")
    if(outputCSV):
        if(printALot): print(f"Outputting to {CSVname}")
        output.to_csv(CSVname)
    else:
        return output
generateCoordsBetweenEpochs(ceres, epoch_init, epoch_final, 2, True, "test.csv", False)

epoch_calc = Time('2026-05-05T00:00:00', scale='utc')
epoch_TEST = Time('2026-10-15T00:00:00', scale='utc')

# POTENTIAL PROBLEM: THE ORBITAL TYPE IS 'MB II' FOR PALLAS
# ceres = orbitFromDict('Ceres', 'KEP', 2.77, 0.0786, 10.6, 73.6, 80.3, 188.7, epoch_calc, 3.34, 0.15) REAL CERES