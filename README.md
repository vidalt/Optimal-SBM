# Optimal-SBM
This project contains the source code and data needed to reproduce all the experiments of the paper "Community Detection in the Stochastic Block Model by Mixed Integer Programming", authored by Breno Serrano and Thibaut Vidal.  
The code is provided as a Julia package.

## Test Environment

This code was tested on Ubuntu 20.04 using Julia Version 1.5.3.  
The `Project.toml` file describes all package dependencies.

## Folder Structure

The repository contains the following folders:

src<br>
instances<br>
experiments<br>

#### src:

* This folder contains the Julia implementation of the exact and heuristic algorithms.

#### instances:

* This folder contains two subfolders, one for each data group (S1 and S2) described in the paper.  
* Each data group subfolder contains two folders, named `in` and `w`. 
The subfolder `in` contains the files that represent the network instances of each data group (as an edge list).
The subfolder `w` contains the ground-truth affinity matrix used in the generation of each network instance.

#### experiments:

* This folder contains Julia scripts which can be used to run the algorithms on all network instances.
The results of the scripts are stored in this same folder by default. 


## Team

Contributors to this code:
* <a href="https://github.com/BSAraujo" target="_blank">`Breno Serrano`</a>
* <a href="https://github.com/vidalt" target="_blank">`Thibaut Vidal`</a>

## License

[![License](http://img.shields.io/:license-mit-blue.svg?style=flat-square)](http://badges.mit-license.org)

- **[MIT license](http://opensource.org/licenses/mit-license.php)**
- Copyright 2021 © Breno Serrano and Thibaut Vidal
