# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 14:01:49 2024

@author: Ville-Veikko Wettenhovi
"""

import numpy as np

    
def rampFilt(N, ww, c, sigma, use2D = False):
    filt = np.linspace(0, N, int(N) // 2 + 1, dtype=np.float32) / N
    w = 2 * np.pi * np.arange(filt.shape[0]) / N
    
    # Choose the window type
    if ww == 'hamming':
        w = (.54 + .46 * np.cos(w / c))
        filt[1:] *= w[1:]
    elif ww == 'hann':
        w = (.5 + .5 * np.cos(w / c))
        filt[1:] *= w[1:]
    elif ww == 'blackman':
        w = (.42659 - .49656 * np.cos(w / c) + .076849 * np.cos(2*w / c))
        filt[1:] *= w[1:]
    elif ww == 'nuttal':
        w = (.355768 - .487396 * np.cos(w / c) + .144232 * np.cos(2*w / c) - .012604 * np.cos(3*w / c))
        filt[1:] *= w[1:]
    elif ww == 'gaussian':
        sigma = 0.1  # Adjust the sigma value accordingly
        w = np.exp((-1/2) * (((np.arange(filt.shape[0] - 1, -1, -1) - N/2)/(sigma*(N/2))) ** 2))
        filt[1:] *= w[1:]
    elif ww == 'shepp-logan':
        w = np.sin(w / (2 * c)) / (w / (2 * c))
        filt[1:] *= w[1:]
    elif ww == 'cosine':
        filt[1:] *= np.cos(w[1:] / (2 * c))
    elif ww == 'parzen':
        L = N + 1
        w = np.arange(filt.shape[0])
        w[np.abs(w) <= L/4] = 1 - 6 * (w[np.abs(w) <= L/4] / (L/2))**2 * (1 - np.abs(w[np.abs(w) <= L/4]) / (L/2))
        w[np.logical_and(np.abs(w) > L/4, np.abs(w) <= L/2)] = 2 * (1 - np.abs(w[np.logical_and(np.abs(w) > L/4, np.abs(w) <= L/2)]) / (L/2))**3
        filt[1:] *= w[1:]
    
    # Truncate frequencies outside the range
    filt[w > np.pi * c] = 0
    
    # Mirror filter coefficients for even length
    if N % 2:
        filt = np.concatenate((filt, filt[-2::-1]))
    else:
        filt = np.concatenate((filt, filt[-2:0:-1]))
    
    if use2D:
        filt = np.tile(filt, (1, int(N)))
        filt = filt * filt.T
        filt = filt / np.max(filt)
    return filt