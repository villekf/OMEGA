# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 13:29:39 2024

@author: Ville-Veikko Wettenhovi
"""



def loadProjectionData(ftype, fpath = '', dims = None, binning = 1, headerBytes = 0, loadAll = True):
    import numpy as np
    import os
    import glob
    import tkinter as tk
    from tkinter.filedialog import askopenfilename
    import re
    if dims == None:
        dims = np.zeros(3, dtype=np.uint64)
    def atoi(text):
        return int(text) if text.isdigit() else text
    
    def natural_keys(text):
        '''
        alist.sort(key=natural_keys) sorts in human order
        http://nedbatchelder.com/blog/200712/human_sorting.html
        (See Toothy's implementation in the comments)
        '''
        return [ atoi(c) for c in re.split(r'(\d+)', text) ]
    if len(fpath) == 0:
        root = tk.Tk()
        root.withdraw()
        fpath = askopenfilename(title='Select projection image')
        if not fpath:
            raise ValueError('No file was selected')
    else:
        if not os.path.exists(fpath):
            print('Specified file was not found! Please select (first) projection image')
            root = tk.Tk()
            root.withdraw()
            fpath = askopenfilename(title='Select projection image')
            if not fpath:
                raise ValueError('No file was selected')
    filename, suffix = os.path.splitext(fpath)
    fpath, file = os.path.split(fpath)
    if len(suffix) == 0:
        flist = glob.glob(fpath + '/*')
    else:
        flist = glob.glob(fpath + '/*' + suffix)
    if headerBytes >=0:
        A = np.fromfile(fpath + "/" + file, dtype = ftype, offset=headerBytes)
    else:
        fsize = os.path.getsize(fpath + "/" + file)
        dt = np.dtype(ftype)
        nRead = (fsize - headerBytes) // dt.itemsize
        A = np.fromfile(fpath + "/" + file, dtype = ftype, count=nRead)
    if loadAll == True:
        nFiles = len(flist)
        if dims[0].item() > 0 and dims[1].item() > 0 and dims[2].item() > 0:
            projData = np.zeros((dims[0].item(),dims[1].item(),dims[2].item()), dtype=ftype)
        else:
            projData = np.zeros((np.size(A)*nFiles/binning**2, 1), dtype=ftype)
        flist.sort(key=natural_keys)
    else:
        nFiles = 1
        if binning > 1:
            if dims[0].item() > 0 and dims[1].item() > 0 and dims[2].item() > 0:
                projData = np.zeros((dims[0].item(),dims[1].item(),dims[2].item()), dtype=ftype)
            else:
                projData = np.zeros((np.size(A)*nFiles/binning**2, 1), dtype=ftype)
    for ll in range(nFiles):
        if nFiles > 1:
            file = flist[ll]
            if headerBytes >=0:
                A = np.fromfile(file, dtype = ftype, offset=headerBytes)
            else:
                fsize = os.path.getsize(fpath + "/" + file)
                dt = np.dtype(ftype)
                nRead = (fsize - headerBytes) // dt.itemsize
                A = np.fromfile(file, dtype = ftype, count=nRead)
        if dims[0].item() > 0 and dims[1].item() > 0:
            if nFiles == 1:
                A = np.reshape(A, (dims[0].item() * binning, dims[1].item() * binning, dims[2].item()))
            else:
                A = np.reshape(A, (dims[0].item() * binning, dims[1].item() * binning))
        if dims[0].item() > 0 and dims[1].item() > 0:
            if nFiles == 1:
                if binning > 1:
                    for kk in range(dims[2].item()):
                        B = np.sum(np.reshape(A[:,:,kk],(binning,-1)),axis=0)
                        B = np.squeeze(np.sum(np.reshape(B,(A.shape[0]//binning,binning,-1)),axis=1))
                        projData[:,:,kk] = B.astype(ftype)
                else:
                    projData = A
            else:
                if binning > 1:
                    B = np.sum(np.reshape(A,(binning,-1)),axis=0)
                    B = np.squeeze(np.sum(np.reshape(B,(A.shape[0]//binning,binning,-1)),axis=1))
                    projData[:,:,ll] = B.astype(ftype)
                else:
                    projData[:,:,ll] = A.astype(ftype)
        else:
            if ll == 0:
                projData = A.astype(ftype)
            else:
                projData = np.append(projData, A.astype(ftype))
    return projData