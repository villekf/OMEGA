This folder contains the example scripts that can be used to reproduce the results in the OMEGA V2 article (link and DOI will be included later). Files containing the word `PET` are PET related, 
`CBCT` ones are (CB)CT related, and `SPECT` SPECT related. Same examples are available for both MATLAB/Octave and Python. If you have extracted all the below example data to your working directory, 
you should be able to run the codes without changes, but if the example data is not in the working directory or in path, you'll need to specify the folder of the data at the beginning of each script.
There is a specific variable (`path`) that takes this path to the folder. Each script can be freely modified or adapted to other cases, but not all adjustable parameters are present. For a more complete
selection of examples, see the `main-files` for MATLAB/Octave and `/source/Python` for Python.

PET example data is available from: https://doi.org/10.5281/zenodo.17185907

SPECT example data is available from: https://doi.org/10.5281/zenodo.17315440

CBCT example data is available from: https://doi.org/10.5281/zenodo.12722386

- `PET_exampleFig2a.m/py` contains an example script that can be used to reproduce the results of Figure 2 (a) from the article. This reconstruction uses PET TOF sinogram data, reconstructed using the improved version
of the Siddon's algorithm (`projector_type = 1` in OMEGA), and using the OSEM algorithm with 1 iteration and 30 subsets.

- `PET_exampleFig2b.m/py` contains an example script that can be used to reproduce the results of Figure 2 (b) from the article. This reconstruction uses PET TOF sinogram data, reconstructed using a volume of intersection based projector
(`projector_type = 3` in OMEGA), and using the OSEM algorithm with 1 iteration and 30 subsets.

- `PET_exampleFig2c.m/py` contains an example script that can be used to reproduce the results of Figure 2 (c) from the article. This reconstruction uses PET TOF list-mode data, reconstructed using the improved version
of the Siddon's algorithm (`projector_type = 1` in OMEGA), and using the PKMA algorithm with 1 iteration and 30 subsets, enhanced with RDP regularization.

- `SPECT_exampleFig3a.m/py` contains an example script that can be used to reproduce the results of Figure 3 (a) from the article. This reconstruction uses SPECT projection data, reconstructed using the improved version
of the Siddon's algorithm (`projector_type = 1` in OMEGA), and using the OSEM algorithm with 5 iterations and 8 subsets.

- `SPECT_exampleFig3b.m/py` contains an example script that can be used to reproduce the results of Figure 3 (b) from the article. This reconstruction uses SPECT projection data, reconstructed using an orthogonal distance-based projector
(`projector_type = 2` in OMEGA), and using the OSEM algorithm with 5 iterations and 8 subsets.

- `SPECT_exampleFig3c.m/py` contains an example script that can be used to reproduce the results of Figure 3 (c) from the article. This reconstruction uses SPECT projection data, reconstructed using the rotation-based projector
(`projector_type = 6` in OMEGA), and using the OSEM algorithm with 5 iterations and 8 subsets.

- `CBCT_exampleFig4a.m/py` contains an example script that can be used to reproduce the results of Figure 4 (a) from the article. This reconstruction uses linearized CBCT projection data, reconstructed using the interpolation-based projector 
(`projector_type = 4` in OMEGA), and using the PDHG algorithm with 8 iterations and 20 subsets.

- `CBCT_exampleFig4b.m/py` contains an example script that can be used to reproduce the results of Figure 4 (b) from the article. This reconstruction uses linearized CBCT projection data, reconstructed using the interpolation-based projector 
(`projector_type = 4` in OMEGA), and using the CV algorithm with 8 iterations and 20 subsets, enhanced with NLRDP regularization.

- `CBCT_exampleFig4c.m/py` contains an example script that can be used to reproduce the results of Figure 4 (c) from the article. This reconstruction uses original CBCT projection data, reconstructed using the interpolation-based projector 
(`projector_type = 4` in OMEGA), and using the PKMA algorithm with 50 iterations and 20 subsets, enhanced with RDP regularization.