# SynapticModel
Code for simulations of a minimal synaptic model of _Drosophila_ T4 direction-selective neurons, as described in our [paper](https://doi.org/10.1167/jov.20.2.2). If you find this useful, please consider citing: 

```
@article{zavatone-veth2020,
    author = {Zavatone-Veth, Jacob A. and Badwan, Bara A. and Clark, Damon A.},
    title = {A minimal synaptic model for direction selective neurons in Drosophila},
    journal = {Journal of Vision},
    volume = {20},
    number = {2},
    pages = {2-2},
    year = {2020},
    month = {02},
    issn = {1534-7362},
    doi = {10.1167/jov.20.2.2},
    url = {https://doi.org/10.1167/jov.20.2.2},
    eprint = {https://arvojournals.org/arvo/content\_public/journal/jov/938384/i0035-8711-217-1-07027.pdf},
}
```

## Getting started

Before these simulations can be run, the `rootDataPath` and `sceneSourcePath` directories must be set in the `SetConfiguration` function, or provided as arguments to that function. A subset of these analyses require the natural scene database from Meyer _et al_. 2014, which may be downloaded [here (click for link)](https://pub.uni-bielefeld.de/data/2689637). The database is provided as a set of .rar archives, each containing a set of .mat files. Extract the .mat files from the archives using a tool such as [7-Zip](https://www.7-zip.org/), and copy the resulting .mat files into a folder named `imageData` within the `sceneSourcePath` root directory.
