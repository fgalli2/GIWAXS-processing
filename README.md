# GIWAXS-processing
A repository containing the python project for the analysis of GIWAXS edf (raw) files.


## Core features
The code heavily relies on the `pyFAI` module for processing, to get from the `edf`
files to the general plots in GIWAXS measurements


## Installment and Requirements
To install the application clone the repository with
```
git clone https://github.com/fgalli2/GIWAXS-processing.git
```

The library is entirely written and tested using `python 3.13.9`, but should
work with any older version of `python 3`. The running code relies on the
following standard libraries:
- `configparser`
- `argparse`
- `datetime`
- `time`
- `sys`
- `os`
- `math`
- `copy`

and on the following scientific libraries:
- `numpy`
- `PIL.Image`
- `matplotlib.pyplot`
- `matplotlib.colors`
- `fabio`
- `pyFAI`


## Usage
To use the library, the user can set all the input information (which will be
later specified throughout a `Tutorial_`) in a configuration file, i.e. the
file with the `.ini` extention. \
This regards the name of the main data folder and the sample file
name, the folder and file names of the `.poni` and the mask, all the
intetgration and plotting ranges and the plotting colormap.\
To run the desired file with the chosen configuration file, the user can use
the command in the command line
```
python3 Image_treatment.py --config=<configuration_file_name>
```
or
```
run Image_treatment.py --config=<configuration_file_name>
```
if they are in Spyder. `<confi_file_name>` is the name, without the
`.ini` extention, of the configuration file. This is done in order to keep
track of the input settings given to the program, but it is optional, with
the default value of `configuration_settings`. 


## Structure
The library is divided in 2 modules:\
[Image_treatment.py]()
consists of the worflow for data extreactiona and general plotting of the
edf file (raw data). It can perform the addition of repetition files, step
needed in order to compensate for the low acquisition time set to avoid
beam damage on organic samples. It converts the data from pixels maps to
ploar coordinates and further to q space. It plots the reciprocal 2D maps,
the chi-q (azimuthal) map and the out-of-plane (OOP), in-plane (IP) and
full azimuthal range integration, saving them both as `.dat` and `.jpg`
files.

[Functions_pyFAI.py]()
Core module with helper functions. It contains three classes: \
- `ReadingInputFunctions` all the functions to read the input information are stored here.
- `ConversionFunctions` functions needed for the data converion.
- `PlottingFunctions` functions needed for platting and saving the 2D maps and 1D cuts.
