from scipy.io import loadmat
import numpy as np


#Recorded data extraction from matlab file
#'Parameters' file has the following format:
#Name = recdata(datapath)[0]
#Unit = recdata(datapath)[1]
#Values = recdata(datapath)[2]

#List of names and units:
# 0 | name:  Time | unit: sec
# 1 | name:  Angle of attack | unit: deg
# 2 | name:  Deflection of elevator trim | unit: deg
# 3 | name:  Force on elevator control wheel | unit: N
# 4 | name:  Engine 1: Fuel mass flow | unit: lbs/hr
# 5 | name:  Engine 2: Fuel mass flow | unit: lbs/hr
# 6 | name:  Engine 1: Inter Turbine Temperature (ITT) | unit: deg C
# 7 | name:  Engine 2: Inter turbine temperature (ITT) | unit: deg C
# 8 | name:  Engine 1: Oil pressure | unit: psi
# 9 | name:  Engine 2: Oil pressure | unit: psi
# 10 | name:  Deflection of the control column (Se or DCOC) | unit: deg
# 11 | name:  Engine 1: Fan speed (N1) | unit: %
# 12 | name:  Engine 1: Turbine speed (N2) | unit: %
# 13 | name:  Engine 2: Fan speed (N1) | unit: %
# 14 | name:  Engine 2: Turbine speed (N2) | unit: %
# 15 | name:  c | unit: l
# 16 | name:  c | unit: l
# 17 | name:  Deflection of aileron (right wing?) | unit: deg
# 18 | name:  Deflection of elevator | unit: deg
# 19 | name:  Deflection of rudder | unit: deg
# 20 | name:  UTC Date DD:MM:YY | unit: ddmmyy
# 21 | name:  UTC Seconds | unit: sec
# 22 | name:  Roll Angle | unit: deg
# 23 | name:  Pitch Angle | unit: deg
# 24 | name:  <no description> | unit: <no units>
# 25 | name:  GNSS Latitude | unit: deg
# 26 | name:  GNSS Longitude | unit: deg
# 27 | name:  Body Roll Rate | unit: deg/s
# 28 | name:  Body Pitch Rate | unit: deg/s
# 29 | name:  Body Yaw Rate | unit: deg/s
# 30 | name:  Body Long Accel | unit: g
# 31 | name:  Body Lat Accel | unit: g
# 32 | name:  Body Norm Accel | unit: g
# 33 | name:  Along Heading Accel | unit: g
# 34 | name:  Cross Heading Accel | unit: g
# 35 | name:  Vertical Accel | unit: g
# 36 | name:  Static Air Temperature | unit: deg C
# 37 | name:  Total Air Temperature | unit: deg C
# 38 | name:  Pressure Altitude (1013.25 mB) | unit: ft
# 39 | name:  Baro Corrected Altitude #1 | unit: ft
# 40 | name:  <no description> | unit: <no units>
# 41 | name:  Mach | unit: mach
# 42 | name:  Computed Airspeed | unit: knots
# 43 | name:  True Airspeed | unit: knots
# 44 | name:  Altitude Rate | unit: ft/min
# 45 | name:  Measurement Running | unit: no unit
# 46 | name:  Number of Measurements Ready | unit: no unit
# 47 | name:  Status of graph | unit: no unit
# 48 | name:  Active Screen | unit: no unit



datapath = "/home/wieger/Documents/FlightDynamics/FTISxprt-20200310_flight3.mat"
def recdata(datapath):
    data = loadmat(datapath)
    data = data['flightdata'][0][0]
    parameters = []

    #Time
    name = data[-1][0][0][2][0]
    unit = data[-1][0][0][1][0]
    seconds = data[-1][0][0][0].transpose()
    parameters.append([name, unit, seconds])

    #Parameters
    for i in range(48):
        value = data[i][0][0][0]
        name = data[i][0][0][2][0][0][0]

        unit = data[i][0][0][1][0][0]
        if len(unit) > 0:
            unit = unit[0]
        else:
            unit = "no unit"

        parameters.append([name, unit, value])
    return parameters

parameters = recdata(datapath)
