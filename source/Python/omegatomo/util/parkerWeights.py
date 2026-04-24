# -*- coding: utf-8 -*-
import numpy as np
import warnings

def _S(x):
    """Smooth S-function transition."""
    y = np.zeros_like(x)
    m1 = (x > -0.5) & (x < 0.5)
    m2 = (x >= 0.5)
    y[m1] = 0.5 * (1 + np.sin(np.pi * x[m1]))
    y[m2] = 1
    return y

def ParkerWeights(options):
    """
    Computes Parker weights as defined in DOI: 10.1118/1.1450132
    
    Parameters:
    options (struct): Struct containing all scan parameters and projection data.
    
    Returns:
    options: The modified options struct where options.SinM has been weighted.
    """
    # Required inputs
    betaIn = options.angles.flatten()
    DSD = options.sourceToDetector
    nU = int(options.nRowsD)
    du = options.dPitchX
    
    # Optional parameters with defaults
    if hasattr(options, 'ParkerWeight') and options.ParkerWeight > 0:
        q = options.ParkerWeight
    else:
        q = 0.25
        
    # Address the trailing space typo in the MATLAB code, checking both
    if hasattr(options, 'detOffsetRow') and np.sum(options.detOffsetRow) > 0:
        detOffset = options.detOffsetRow
    else:
        detOffset = 0
        
    # Input validation
    if not (np.isscalar(DSD) and np.isscalar(nU) and np.isscalar(du)):
        raise ValueError('sourceToDetector, nRowsD, and dPitchX must be scalars.')
        
    if q <= 0 or q > 1:
        raise ValueError('options.ParkerWeight must satisfy 0 < q <= 1.')
        
    if len(betaIn) < 2:
        raise ValueError('options.angles must contain at least two projection angles.')
        
    # Flat-panel detector coordinate u
    u = (np.arange(nU) - (nU - 1) / 2) * du + detOffset
    
    # Horizontal fan angle alpha for each detector pixel
    alpha = np.arctan2(u, DSD)
    
    # Use actual maximum fan half-angle from detector edges
    delta = np.max(np.abs(alpha))
    
    # Relative angles and direction checking
    betaRel = betaIn - betaIn[0]
    if betaRel[-1] < 0:
        betaRel = np.abs(betaRel)
        alpha = -alpha
        
    scanRange = betaRel[-1]
    
    # Short-scan validity
    minShortScan = np.pi + 2 * delta
    if scanRange + 1e-8 < minShortScan:
        raise ValueError(f'Angular range is too short. Need at least pi + 2*delta = {minShortScan:.6f} rad, got {scanRange:.6f} rad.')
        
    if scanRange > 2 * np.pi + 1e-8:
        warnings.warn('Angular range exceeds 2*pi. Weights are intended for short/over-scan up to 2*pi.')
        
    epsilon = max(scanRange - (np.pi + 2 * delta), 0)
    
    # Compute weights
    for iu in range(nU):
        g = alpha[iu]
        
        B = epsilon + 2 * delta - 2 * g
        b = q * B
        B2 = epsilon + 2 * delta + 2 * g
        b2 = q * B2
        
        if b <= 0 or b2 <= 0:
            raise ValueError('Encountered nonpositive transition length b. Check geometry.')
            
        x1 = betaRel / b - 0.5
        x2 = (betaRel - B) / b + 0.5
        x3 = (betaRel - np.pi + 2 * g) / b2 - 0.5
        x4 = (betaRel - np.pi - 2 * delta - epsilon) / b2 + 0.5
        
        w2 = 0.5 * (_S(x1) + _S(x2) - _S(x3) - _S(x4))
        w2 = np.clip(w2, 0, 1)
        options.SinM[iu, :, :] *= w2.reshape(1, -1)
            
    # return options