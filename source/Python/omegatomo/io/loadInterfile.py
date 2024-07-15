# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 13:36:22 2024

@author: Mad Gigerdi
"""
import re
import numpy as np
import os

def loadInterfile(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    n_dim1 = 1
    n_dim2 = 1
    n_dim3 = 1
    n_dim4 = 1
    n_dim5 = 1
    type = 'uint16'
    machinefmt = 'l'

    for line in lines:
        if 'number format' in line:
            if 'float' in line or 'short float' in line:
                type = 'float32'
            elif 'double' in line or 'long float' in line:
                type = 'float64'
            elif 'unsigned integer' in line:
                for line2 in lines:
                    if 'number of bytes per pixel' in line2:
                        joku = re.search(r'\d+', line2).group()
                if joku == '2':
                    joku = '16'
                elif joku == '1':
                    joku = '8'
                elif joku == '4':
                    joku = '32'
                else:
                    joku = '64'
                type = 'uint' + joku
            elif 'signed integer' in line:
                for line2 in lines:
                    if 'number of bytes per pixel' in line2:
                        joku = re.search(r'\d+', line2).group()
                if joku == '2':
                    joku = '16'
                elif joku == '1':
                    joku = '8'
                elif joku == '4':
                    joku = '32'
                else:
                    joku = '64'
                type = 'int' + joku
        elif 'matrix size[1]' in line or 'matrix size [1]' in line:
            n_dim1 = int(re.findall(r'\d+', line)[1])
        elif 'matrix size[2]' in line or 'matrix size [2]' in line:
            n_dim2 = int(re.findall(r'\d+', line)[1])
        elif 'matrix size[3]' in line or 'matrix size [3]' in line:
            n_dim3 = int(re.findall(r'\d+', line)[1])
        elif 'matrix size[4]' in line or 'matrix size [4]' in line:
            n_dim4 = int(re.findall(r'\d+', line)[1])
        elif 'matrix size[5]' in line or 'matrix size [5]' in line:
            n_dim5 = int(re.findall(r'\d+', line)[1])
        elif 'imagedata byte order' in line or 'image data byte order' in line:
            if 'LITTLEENDIAN' in line:
                machinefmt = 'l'
            elif 'BIGENDIAN' in line:
                machinefmt = 'b'
        elif 'name of data file' in line:
            f_name = line.split('=')[1].strip()
            if '/' in f_name:
                f_name = f_name.split('/')[-1]

    if n_dim3 == 1:
        for line in lines:
            if 'number of slices' in line:
                n_dim3 = int(re.findall(r'\d+', line)[0])
                break
        if n_dim3 == 1 and n_dim4 == 1 and n_dim5 == 1:
            for line in lines:
                if 'total number of images' in line:
                    n_dim3 = int(re.findall(r'\d+', line)[0])
                    break

    if filename.endswith('.hdr') or filename.endswith('.h33'):
        filename = filename[:-3] + 'img'
    path = os.path.split(filename)[0] + '/'

    try:
        with open(path + f_name, 'rb') as f:
            output = np.fromfile(f, dtype=type, count=-1)
    except FileNotFoundError: 
        with open(filename, 'rb') as f:
            output = np.fromfile(f, dtype=type, count=-1)
    
    try:
        output = np.reshape(output, (n_dim1, n_dim2, n_dim3, n_dim4, n_dim5), order='F')
    except ValueError:
        output = np.reshape(output, (n_dim1, n_dim2, -1), order='F')
    
    return output