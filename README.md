# preCICE-adapter for the CFD code OpenFOAM - Modified version for Ishaan

This works only for FF and does not require yaml-cpp. The config is hard-coded into Adapter.C, method `configFileRead()`.

There are two variants, the one commented out: `Fluid1` and `Fluid2`. When changing between the two, also change the name of the produced library in the file `Make/files` (variable `LIB`).
