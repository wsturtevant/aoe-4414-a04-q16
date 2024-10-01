# ecef_to_sez.py
#
# Usage: python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km
#  Converts ECEF (Earth-Centered, Earth-Fixed) coordinates to SEZ (South, East, Zenith) frame
# Parameters:
#  o_x_km: ECEF X-coordinate for origin of SEZ frame in km
#  o_y_km: ECEF Y-coordinate for origin of SEZ frame in km
#  o_z_km: ECEF Z-coordinate for origin of SEZ frame in km
#  x_km: ECEF X-coordinate of position in km
#  y_km: ECEF Y-coordinate of position in km
#  z_km: ECEF Z-coordinate of position in km
# Output:
#  Prints the SEZ s, e, and z coordinates in km
#
# Written by Wheat Sturtevant
# Other contributors: None
#
# This work is licensed under CC BY-SA 4.0

import math  # math module for trigonometric functions
import sys   # argv for command-line arguments

# Constants for LLH calculation
R_E_KM = 6378.137
E_E    = 0.081819221456

# Helper function for denominator in LLH calculation
def calc_denom(ecc, lat_rad):
    return math.sqrt(1.0-(ecc**2)*(math.sin(lat_rad)**2))

# Convert ECEF to LLH (Latitude, Longitude, HAE)
def ecef_to_llh(r_x_km, r_y_km, r_z_km):
    # Calculate longitude
    lon_rad = math.atan2(r_y_km, r_x_km)
    lon_deg = lon_rad*180.0/math.pi

    # Initialize latitude
    lat_rad = math.asin(r_z_km/math.sqrt(r_x_km**2+r_y_km**2+r_z_km**2))
    r_lon_km = math.sqrt(r_x_km**2+r_y_km**2)
    prev_lat_rad = float('nan')

    # Find latitude via iteration
    count = 0
    c_E = float('nan')
    while math.isnan(prev_lat_rad) or abs(lat_rad-prev_lat_rad)>1e-7:
        denom = calc_denom(E_E, lat_rad)
        c_E = R_E_KM/denom
        prev_lat_rad = lat_rad
        lat_rad = math.atan((r_z_km+c_E*(E_E**2)*math.sin(lat_rad))/r_lon_km)
        count+=1
    
    # Calculate HAE
    hae_km = r_lon_km/math.cos(lat_rad)-c_E
    
    # Out: latitude, longitude, and HAE
    return lat_rad, lon_rad, hae_km

# Convert from ECEF to SEZ using the updated lat/lon values from LLH
def ecef_to_sez(o_x_km, o_y_km, o_z_km, x_km, y_km, z_km):
    # Calculate observer's latitude and longitude using LLH conversion
    lat_o_rad, lon_o_rad, _ = ecef_to_llh(o_x_km, o_y_km, o_z_km)

    # Find difference in ECEF
    dx = x_km-o_x_km
    dy = y_km-o_y_km
    dz = z_km-o_z_km

    # Rotation matrix to convert ECEF to SEZ
    sin_lat_o = math.sin(lat_o_rad)
    cos_lat_o = math.cos(lat_o_rad)
    sin_lon_o = math.sin(lon_o_rad)
    cos_lon_o = math.cos(lon_o_rad)

    # Apply the rotation matrix to get SEZ coordinates
    s_km =-dz*cos_lat_o + dx*cos_lon_o*sin_lat_o + dy*sin_lat_o*sin_lon_o
    e_km = dy*cos_lon_o - dx*sin_lon_o
    z_km = dx*cos_lat_o*cos_lon_o + dz*sin_lat_o + dy*cos_lat_o*sin_lon_o

    return s_km, e_km, z_km

# Initialize script arguments
if len(sys.argv) == 7:
    o_x_km = float(sys.argv[1])
    o_y_km = float(sys.argv[2])
    o_z_km = float(sys.argv[3])
    x_km = float(sys.argv[4])
    y_km = float(sys.argv[5])
    z_km = float(sys.argv[6])
else:
    print('Usage: python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km')
    exit()

# Convert ECEF to SEZ
s_km, e_km, z_km = ecef_to_sez(o_x_km, o_y_km, o_z_km, x_km, y_km, z_km)

# Print SEZ coordinates in km
print('s: '+str(s_km)+' km')
print('s: '+str(e_km)+' km')
print('s: '+str(z_km)+' km')