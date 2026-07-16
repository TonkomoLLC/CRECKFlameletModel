## Steady-state laminar flamelet model port for OpenFOAM 14

The libraries you are using here are developed by Alberto Cuoci and his team (CRECK Modeling group). For more information have a look at the official website http://creckmodeling.chem.polimi.it/

## Supported OpenFOAM versions
 - OpenFOAM v14 (this port)
 - OpenFOAM v13 (previous port)
 - OpenFOAM v10 (Tonkomo GitHub)
 - OpenFOAM v8
 - OpenFOAM v7
 - OpenFOAM 4.x
 - OpenFOAM 2.3.x
 - OpenFOAM 2.2.x
 - OpenFOAM 2.1.x

## Introduction
Tobias worked with that flamelet-tool over 8 month during his Master thesis in 2012 and analysed the whole thermo model and the flamelet extraction libraries. For that Tobias made some validations and compared the results with ANSYS-CFX and measurements. There are some "bugs" in the original version which are removed in this version. Furthermore, you can use the SIMPLEC algorithm for solving steady-state combustions and the PIMPLE algorithm for transient calculations.

The library only contains a steady-state solver. If you create the transient one, feel free to make a push request.


## Discussion
A discussion can be followed here: http://www.cfd-online.com/Forums/openfoam-programming-development/99645-libopensmoke.html

## Compiling

The legacy documents remain useful for model background, but their build instructions predate the modular solver framework. For OpenFOAM 14, load the OpenFOAM environment and use the supplied build script.

```bash
source /opt/openfoam14/etc/bashrc
cd CRECKFlameletModel-OpenFOAM-14
./Allwmake -j4

# Run the supplied OF14 tutorial
cd tutorials/flameletSimpleFoam/Sandia_COH2N2
./Allrun
```

## Changes | Features | Documentation
+ Have a look into the documentation folder

## Important | Validation

+ This OF14 port has been statically checked against the OpenFOAM 14 interfaces. The supplied SANDIA case should still be compared with the original validation data before relying on quantitative results.

## Notice | Warranty
+ No warranty is provided for numerical accuracy or fitness for a particular application.

## Older versions
+ the modified version of 7.x can be found inside this repository
+ the modified version of 4.x can be found inside this repository
+ the modified version of 2.3.x can be found inside this repository
+ the modified version of 2.2.x can be found inside this repository
+ the modified version of 2.1.x can be found inside this repository
+ the origin version of 2.1.x can be found here: https://github.com/wyldckat/libOpenSMOKE
+ the origin version of 1.7.x can be found here: https://github.com/wyldckat/libOpenSMOKE


