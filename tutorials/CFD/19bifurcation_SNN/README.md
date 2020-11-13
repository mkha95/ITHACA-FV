# TUTORIAL 19: BIFURCATION ANALYSIS FOR STEADY NS

This is the case folder for a bifurcation analysis on the coanda effect.
In order to effectively run the tutorial, after having compiled the ITHACA-FV
library with the instruction provided in https://github.com/mathLab/ITHACA-FV.
Make sure to source the file ITHACA-FV/etc/bash.
Go to the ITHACA_BIF folder and run wmake in order to compile the shared
library linked to the case test.
Then from the folder 19bifurcation_SNN run again wmake to compile the
bifurcationNSS.exe executable;

In order to run the executable after the compilation it is sufficient to
execute the command:


### 1. Prerequisites
Install **ITHACA-FV** using the instruction in  (http://mathlab.github.io/ITHACA-FV/).
Move to the main ITHACA-FV folder and source the bashrc file as follow:
```
source etc/bashrc
```

### 2.  Important files and folders:

- mesh.sh: this is a basic bash script that allows changing the current mesh,
it must be run with one of the available option: refined, coarse, structured.
The relevant information are extracted from the utilities folder.
Example of usage

```
./mesh.sh coarse
```


- system folder: contains all the dictionary necessary to run the case, in
particular:
           BIFURCATIONdict: needed to set all the user defined bifurcation
           parameters

           ITHACAdict: containing information about the POD phase
-ITHACAoutput: provided all the relevant ouputs of the program run, in
particular one can find the following subdirectoris:
           Offline: which will contain the fields computed for the FOM problem
           during the offline phase

           Online: which will contain the fields computed for the ROM problem
           during the online phase




