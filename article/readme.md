This folder contains the example scripts that can be used to reproduce the results in the second OMEGA article (link and DOI will be included later). Files containing the word `PET` are PET related, `CBCT` ones are (CB)CT related and `SPECT` SPECT related.

Same examples are available for both MATLAB/Octave and Python. There are two different types of examples. Files starting with `OMEGAv2_` are examples that allow to easily run all three example cases. You simply need to select the number (0 corresponds to the left image, 
1 to the middle and 2 to the right), specify the path to the example data (not required for CBCT data on MATLAB/Octave as long as the data is in path), and optionally the GPU number in case you have multiple GPU devices. For CBCT cases, you can also
select the "high-dimensional computing" option that allows the reconstructions to be performed even on GPUs with less than 12 GB of memory (if set to true). If you are interested in simply reproducing the results, the `OMEGAv2_` files are highly recommended.

Alternatively, you can run each example result as an individual script. These can also be useful as general examples since you can easily modify them to suit your needs. These are labeled as `_exampleX` files, where X can be 1, 2 or 3. 1 is the 
left image, 2 the middle and 3 the right image. For example, `PET_example2.m/py` can be used to reconstruct the middle PET example image. 

Note that `PET_TOF_example` is the function called by `OMEGAv2_PET` and the unnumbered `CBCT_example` is called by `OMEGAv2_CBCT`.

PET example data is available from: https://doi.org/10.5281/zenodo.17185907
SPECT example data is available from: 
CBCT example data is available from: https://doi.org/10.5281/zenodo.12722386