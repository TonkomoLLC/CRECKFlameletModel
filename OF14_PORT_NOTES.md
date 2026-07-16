# CRECK Flamelet Model: OpenFOAM 14 port notes

## Applied source changes

1. Replaced the removed `applications/modules/fluidSolver` include path with
   `applications/modules/basicFluidSolver`.
2. Replaced the OF13 module libraries `fluidSolver` and `isothermalFluid` with
   the OF14 libraries `basicFluidSolver` and `isothermalFluidSolver`.
3. Registered the flamelet thermo package with `makeFluidThermo` instead of
   `makeThermo`. In OF14, fluid thermos must be added to the `basicThermo`,
   `fluidThermo`, and derived thermo run-time selection tables.
4. Added OF14 thermophysical library dependencies, including
   `fluidThermophysicalModels` and `thermophysicalProperties`.
5. Migrated the supplied tutorial dictionaries:
   - `thermophysicalProperties` -> `physicalProperties`
   - `turbulenceProperties` -> `momentumTransport`
   - `RASModel` -> `model`
   - `Hf` -> `hf`
   - custom module loaded through `libs ("libSLFMFluid.so")`
6. Added top-level `Allwmake`/`Allwclean` and tutorial `Allrun`/`Allclean` scripts.

## Build

```bash
source /opt/openfoam14/etc/bashrc
cd CRECKFlameletModel-OpenFOAM-14
./Allwclean
./Allwmake -j4
```

Expected user libraries:

```text
$FOAM_USER_LIBBIN/libOpenSMOKE_common.so
$FOAM_USER_LIBBIN/libOpenSMOKE_turbulent_flamelets.so
$FOAM_USER_LIBBIN/libflameletThermophysicalModels.so
$FOAM_USER_LIBBIN/libSLFMFluid.so
```

## Tutorial

```bash
cd tutorials/flameletSimpleFoam/Sandia_COH2N2
./Allclean
./Allrun
```

The case is executed by `foamRun`; the `solver SLFMFluid` entry selects the
module, and the `libs` entry loads the user library.

## Validation status

This archive was statically checked against the OpenFOAM 14 source interfaces
and against the supplied compile logs. The current execution environment does
not contain an OpenFOAM 14 installation, so an actual `wmake` and tutorial run
could not be performed here. Any subsequent compiler error should be captured
from the first error onward; later errors are often cascades.
