# Reproducibility 

## Compiling 
Code for the simulation model as well as the Markov Chain Monte Carlo (MCMC)
is written in C++ language. The files contained at the project folder: 
*random.cpp*, *utils.cpp*, *utils.h*, *random.h* and *json.hpp*, 
contain useful functions and random number generators used in the
stochastic simulations of the learning model. The source file for the 
analysis (*ActCritMCMC.cpp*) can be found also in the project folder.
Executables get JSON files as parameter inputs. JSON files are
integrated into c++ using the library contained in *json.hpp* (for more info see:
https://github.com/nlohmann/json). This library is written with c++11 
standards. With all these files, executables can be compile for the simulation 
with standard c++ compilers. For example for compilation with g++ use command:

_g++ ActCrip.cpp random.cpp utils.cpp -o $PROGRAMNAME -std=c++11_

Likewise, executable file necessary to generate contour plots can be generated
by compiling the file *ActCrit.cpp*

## Generating parameter Values for the MCMC
Files containing parameter values can be generated using the 
R file: parTOjson.R. Comments inside the file indicate how to set the different 
parameters of the models and the MCMC. In that same, file following the comments 
it is possible to generate parameter values for generating the predictions
used to produce Figure 2. 

## Running MCMCs
The executable files compiled in section 1, get the parameter files as input. 
With the settings used as default, the executable for the MCMC
runs one chain for each parameter files that gets as input,
which produces one output file and places it in the directory specified 
in the parameter file. Likewise, the executable that produces the predictions
gets one parameter file for each parameter combination. Using the 
parameter generating file (*parTOjson.R*) this files can be systematically 
generated to produce the data corresponding to the contour plots 
(Figure 2, left hand side panels). 

## Visualization
Description and analysis of the resuts are embebbed in R chunks in the 
Rmd file *manuscript_1.0.Rmd*. These R chunks source other R files contained
in the repository. Bellow a short description of the files necessary for the 
analysis and their role.

## Description of relevan files
*ActCrit.cpp*: c++ code of the learning model. Used to generate thre predictions
  visualized as contour plots Figure 2, left panels.
*ActCritMCMC.cpp*: c++ code of the bayesian analysis, the code runs the 
  Monte Carlo Markov Chains. It also generates the predictions for each data 
  point of the data set.
*aesth_par*: Defines a set of aesthetic parameters to be used in 
  the visualizations.
*analysisMCMC*: Visualization of the marginal posterior distribution and MCMC
  diagnostics. It generates Fig. 1, Fig S. 2,3,4 and 5. 
*data2interp.R*: Definition of function used in the interpolation that 
  generates contour plots in Fig. 2. 
*json.hpp*: Header file necessary to use json files in the c++ code.
*loadData.R*: Define accesory functions used to load simulation output.
*loadFiledData.R*: Loads the two sources of data (field and experimental),
  extracts the preference for the visitor clients (according to the
  criteri, see Statistical analysis in the main text) and exports the data 
  file that is used by the executables. It also uses the raw data to produce 
  Fig S. 1.
*manFigures.R*: Loads data from the chains and predictions. Computes the 
  interpolation to generate the countor plots. It produces Fig. 2.
*parTOjson.R*: generate different types of json files with the parameter values.
*random.cpp*, *random.h*, *utils.cpp*, *utils.h*: Headers and cpp files 
  necessary to compile the executables. They define convenient functions and
  random number generators. 
*manuscript_1.0/*: files necessary to compile the R markdown file that 
  renders the manuscript. Among those *Cleanerlearning.bib* the reference list.