# TRANSCOMPP

TRANSCOMPP is a tool that allows users to quantify phenotypic plasticity starting from single-cell or bulk measurements of phenotype. TRANSCOMPP uses Markov modeling and optimization to estimate best-fit values for stochastic transition rates, transition rate intervals, and phenotype-specific proliferation parameters.

## Input Data

Input data to TRANSCOMPP is either in the form of single-cell measurements of phenotype-related attributes, or in the form of summary statistics from bulk data indicating the relative abundance of each phenotypic state in the population. More details about the types and formatting requirements for input data can be found in the *Sample Data* folder

## Installation

### Users with MATLAB
If you have access to MATLAB and its optimization toolbox, please install TRANSCOMPP by navigating to the *MATLAB Package* folder and executing the *Transcompp_MATLAB.mlappinstall* file.

### Users without MATLAB - Windows

Navigate to the *Windows* folder. Run the *Transcompp.exe* file and follow the on-screen instructions to install. The program first installs the MATLAB Runtime environment (Free) that is required to execute MATLAB-based software. For more details about the runtime environment please click [here](https://www.mathworks.com/products/compiler/matlab-runtime.html).

### Users without MATLAB - Mac OS
Navigate to the *Mac OS* folder and unzip the the *Transcompp_MacOS.zip* file. Run the unpackaged *Transcompp.app* file to install. The program first installs the MATLAB Runtime environment (Free) that is required to execute MATLAB-based software. For more details about the runtime environment please click [here](https://www.mathworks.com/products/compiler/matlab-runtime.html).


## Usage
When using the MATLAB package, check for the installed Transcompp App in the Apps toolbar. 

When using the Windows or Mac OS versions, double click the executable files to start TRANSCOMPP. 

For more details about using TRANSCOMPP, please refer to the attached User guide.


## License
[GPLv3](https://choosealicense.com/licenses/gpl-3.0/)