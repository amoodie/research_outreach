This flooding risk module is intended to help educators demonstrate some basic physics of why rivers flood. A computer code with a fully interactive graphical user interface module (GUI) is interacted with alongside worksheets designed for students of varying ages.



# Installation and setup guides

Visit the section of the text below for more information on installing and executing the FloodingModule program on your computer. 

Note that if you have a full Matlab license, there is a section at the end for running the module from source. This method is recommended over using the MCR if possible.

## Windows (Windows 10)
* Visit the 'download' folder of this repository ([https://github.com/amoodie/research_outreach/tree/master/flooding_risk/download](https://github.com/amoodie/research_outreach/tree/master/flooding_risk/download)) and download the 64-bit Windows 10 FloodingModule installer.
* Run the installer by navigating to the download directory and executing the installer by double clicking. Note that you may have to run this as administrator, depending on where you intend to install the program and Matlab MCR.
* Follow the default options. If you must make any changes during the install process, be sure to note them (e.g., the install location), as you will may this information to run the module program.
* You can now run the program by navigating to the install location. 
  * note that the program may take up to a minute to launch the first time.

## Mac (OSX)


## Linux (Ubuntu)
* Visit the 'download' folder of this repository ([https://github.com/amoodie/research_outreach/tree/master/flooding_risk/download](https://github.com/amoodie/research_outreach/tree/master/flooding_risk/download)) and download the 64-bit Linux FloodingModule installer.
* Run the installer by `cd`ing into the download location and running `./FloodingModule_lnx64.install`. Note that you may have to run this as root (i.e., with `sudo`) depending on where you intend to install the program and Matlab MCR.
* Follow the default options. If you must make any changes during the install process, be sure to note them (e.g., the install location), as you will need this information to run the module program.
* You can now run the program with `<FloodingModuleRoot>/application/run_FloodingModule.sh <MCRpath>`.
  * or, if you followed the default options, `/usr/local/FloodingModule/application/run_FloodingModule.sh /usr/local/MATLAB/MATLAB_Runtime/v92/`.
  * note that the program may take up to a minute to launch the first time.

## Run from source
* The source Matlab code (`.m` file) is included in the repository in the 'source' folder. 
Two options for running the module
1. Save the file to a folder on your Matlab path and type 

# Module worksheet information

Worksheets are being written to accompany the GUI module. The aim of the worksheets is to help guide a discussion about flooding on deltas or rivers, as may be discussed in Earth Science courses. The modules currently available are targeted at sprecific age groups:

* High School Earth Science -- 9th to 12th graders




# Disclaimer

This module utilizes a Matlab Compiled Runtime (MCR) program for backend calculations and plotting of the model results.
The program accompanying this module is to be run as a standalone application, but relies on the Matlab MCR; the program can be be downloaded at \url{http://www.coastalsustainability.rice.edu/outreach/}.
For help troubleshooting the Matlab executable, please visit \url{https://www.mathworks.com/products/compiler/mcr.html}.
The source code can be viewed at \url{https://github.com/amoodie/research_outreach}.

The module was created by Andrew J. Moodie as part of an National Science Foundation funded research project assessing the sustainability of anthropogenically influenced deltas.
The research is supported by Grant No. 1427262 and an NSF Graduate Research Fellowship to A.J.M. under Grant No. 1450681.
Any opinion, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
For more information, please visit \url{http://www.coastalsustainability.rice.edu/outreach/}.
