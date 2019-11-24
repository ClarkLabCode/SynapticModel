# SynapticModel
Code for simulations of the minimal synaptic model. Details can be found in our [preprint](https://doi.org/10.1101/833970).


## Getting started

Before these simulations can be run, the `rootDataPath` and `sceneSourcePath` directories must be set in the `SetConfiguration` function, or provided as arguments to that function. A subset of these analyses require the natural scene database from Meyer _et al_. 2014, which may be downloaded [here (click for link)](https://pub.uni-bielefeld.de/data/2689637). The database is provided as a set of .rar archives, each containing a set of .mat files. Extract the .mat files from the archives using a tool such as [7-Zip](https://www.7-zip.org/), and copy the resulting .mat files into a folder named `imageData` within the `sceneSourcePath` root directory.
