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
- Install **ITHACA-FV** using the instruction in  (http://mathlab.github.io/ITHACA-FV/).
Move to the main ITHACA-FV folder and source the bashrc file as follow:
```
source etc/bashrc
```


### 2. Running the test case
This tutorial uses a shared library which is not compiled during the
installation of ITHACA-FV. In order to compile it move to the ITHACA_BIF folder
and run:
```
wmake
```
This will create a shared libray object which will be stored in $(FOAM_USER_LIBBIN)

After having compiled the library, one needs to compile the executable for the
test case, this can be done again with:

```
wmake
```

After having executed the former command an executable called
bifurcationNSS.exe will be created in the current directory. In order to run
the case you run:
```
./Allrun
```

In order to clean up the directory and removing the output data you can run:

```
./Allclean
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
* [**BIFURCATION**]needed to set all the user defined bifurcation parameters
* [**ITHACAdict**] containing information about the POD phase
-ITHACAoutput: provided all the relevant ouputs of the program run, in
particular one can find the following subdirectoris:
* [**Offline**] which will contain the fields computed for the ROM problem during the online phase
* [**Online**] which will contain the fields computed for the ROM problem during the online phase
* [**POD**] containing information about the POD phase
-ITHACA_BIF: folder containing the files needed to compile the ITHACA_BIF
shared library which will be stored in $(FOAM_USER_LIBBIN)




