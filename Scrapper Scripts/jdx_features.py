from jcamp import jcamp_readfile
from os import listdir
from scipy import interpolate
from rdkit import Chem
# from smarts import fg_list_original, fg_list_extended
from PIL import Image
import os
import numpy as np
import sys
import csv
import re
import pandas as pd

feature_vectors = []
cas_no = []

jdx = pd.read_csv("sdf_left.csv")
jdx = set(jdx['cas_no'])

def convert_x(x_in, unit_from, unit_to):
    """Convert between micrometer and wavenumber."""
    if unit_to == 'micrometers' and unit_from == 'MICROMETERS':
        x_out = x_in
        return x_out
    elif unit_to == 'cm-1' and unit_from in ['1/CM', 'cm-1', '1/cm', 'Wavenumbers (cm-1)']:
        x_out = x_in
        return x_out
    elif unit_to == 'micrometers' and unit_from in ['1/CM', 'cm-1', '1/cm', 'Wavenumbers (cm-1)']:
        x_out = np.array([10 ** 4 / i for i in x_in])
        return x_out
    elif unit_to == 'cm-1' and unit_from == 'MICROMETERS':
        x_out = np.array([10 ** 4 / i for i in x_in])
        return x_out

def convert_y(y_in, unit_from, unit_to):
    """Convert between absorbance and trasmittance."""
    if unit_to == 'transmittance' and unit_from in ['% Transmission', 'TRANSMITTANCE', 'Transmittance']:
        y_out = y_in
        return y_out
    elif unit_to == 'absorbance' and unit_from == 'ABSORBANCE':
        y_out = y_in
        return y_out
    elif unit_to == 'transmittance' and unit_from == 'ABSORBANCE':
        y_out = np.array([1 / 10 ** j for j in y_in])
        return y_out
    elif unit_to == 'absorbance' and unit_from in ['% Transmission', 'TRANSMITTANCE', 'Transmittance']:
        y_out = np.array([np.log10(1 / j) if j != 0 else 0 for j in y_in])
        return y_out

def get_unique(x_in, y_in):
    """Removes duplicates in x and takes smallest y value for each x value."""
    x_out = sorted(list(set(x_in)), reverse=True)
    y_out = []
    for i in x_out:
        y_temp = []
        for ii, j in zip(x_in, y_in):
            if i == ii:
                y_temp.append(j)
        y_out.append(min(y_temp))
    return x_out, y_out

# Specify the number of data points you want, in this case, 4000.
num_data_points = 4000

# here make list which will contain all the file paths
parent_directory = 'jdx35/'
file_list = [parent_directory + file for file in os.listdir(parent_directory) if file.endswith('_0.jdx')]
new_file_list = []
for i in jdx:
    if 'jdx35/' + str(i) + '_0.jdx' in file_list:
        new_file_list.append('jdx35/' + str(i) + '_0.jdx')

for i in new_file_list:
    try:
        jcamp_dict = jcamp_readfile(i)
        print(i.split("_0.jdx")[0].split("jdx35/")[1])
    except:
        print("error")
    if jcamp_dict['x'] is None or len(jcamp_dict['x']) == 0:
        print("error")
    if jcamp_dict['yunits'] is not None:
        if jcamp_dict['yunits'] in ['dispersion index', 'absorption index', '(micromol/mol)-1m-1 (base 10)']:
            print("error")
    elif jcamp_dict['ylabel'] is not None:
        if jcamp_dict['ylabel'] in ['dispersion index', 'absorption index', '(micromol/mol)-1m-1 (base 10)']:
            print("error")
    if 'xunits' in jcamp_dict:
        xunit = jcamp_dict['xunits']
    if 'yunits' in jcamp_dict:
        yunit = jcamp_dict['yunits']
    if 'xlabel' in jcamp_dict:
        xunit = jcamp_dict['xlabel']
    if 'ylabel' in jcamp_dict:
        yunit = jcamp_dict['ylabel']

    x = jcamp_dict['x']
    y = jcamp_dict['y']

    if 'xunits' in jcamp_dict:
        xunit = jcamp_dict['xunits']
    if 'yunits' in jcamp_dict:
        yunit = jcamp_dict['yunits']
    if 'xlabel' in jcamp_dict:
        xunit = jcamp_dict['xlabel']
    if 'ylabel' in jcamp_dict:
        yunit = jcamp_dict['ylabel']

    x = convert_x(x, xunit, 'cm-1')
    y = convert_y(y, yunit, 'absorbance')

    y_min = min(y)
    y_max = max(y)
    x_min = min(x)
    x_max = max(x)

    if y_max > 1:
        y = [1 if j > 1 else j for j in y]
    if y_min < 0:
        y = [0 if j < 0 else j for j in y]

    if x_max > 4000:
        idx = next(i for i, xx in enumerate(x) if xx >= 4000)
        x = x[:idx + 1]
        y = y[:idx + 1]
    if x_min < 400:
        idx = next(i for i, xx in enumerate(x) if xx >= 400)
        x = x[idx - 1:]
        y = y[idx - 1:]
    if x_max < 4000:
        x = np.append(x, np.linspace(x_max, 4000, num_data_points - len(x)))
        y = np.append(y, [y[-1]] * (num_data_points - len(y)))
    if x_min > 400:
        x = np.insert(x, 0, np.linspace(400, x_min, num_data_points - len(x)))
        y = np.insert(y, 0, [y[0]] * (num_data_points - len(y)))

    spectrum = get_unique(x, y)
    x = np.linspace(4000, 1, num_data_points)
    f = interpolate.interp1d(spectrum[0], spectrum[1], kind='slinear', fill_value='extrapolate')
    y = f(x)

    print(i.split("_0.jdx")[0].split("jdx35/")[1])
    feature_vectors.append(y)
    cas_no.append(i.split("_0.jdx")[0].split("jdx35/")[1])

df = pd.DataFrame(feature_vectors, columns=None)
df.insert(loc=0, column='cas_no', value=cas_no)
df.to_csv("jdx_4000.csv", index=False)
