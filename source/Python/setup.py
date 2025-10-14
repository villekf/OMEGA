# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

setup(
    name='Open-source multi-dimensional tomographic reconstruction software',
    version='2.0.0',
    description='OMEGA, or open-source multi-dimensional tomographic reconstruction software, is an image reconstruction toolkit for tomographic imaging, optimized for positron emission tomography (PET), single photon emission computed tomography (SPECT) and computed tomography (CT). However, any ray-tracing based tomographic imaging is supported. ',
    author='Ville-Veikko Wettenhovi',
    author_email='villewe@uef.fi',
    url='https://github.com/villekf/OMEGA',  # Update with your repository URL
    packages=find_packages(),  # Automatically find packages in the directory
    install_requires=[
        'numpy',
        'pymatreader',
        'scikit-image',
        'arrayfire',
        'SimpleITK',
        'h5py',
        # Add other dependencies here
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GPL-3.0 license',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.8',  # Specify the minimum Python version required
)