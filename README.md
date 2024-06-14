# pC-SAC Installation 

This guide provides instructions for installing the pC-SAC. 
pC-SAC is a tool that uses adaptive importance sampling technique with sequential Monte Carlo to generate a set of three-dimensional chromatin chains to enhance low-resolution Hi-C data. 

## Prerequisites

Before installing pC-SAC, make sure the following dependencies are installed on your system:

- CMake
- Boost
- Eigen

These can typically be loaded on an HPC environment using module commands, like so:

    module load cmake
    module load boost
    module load eigen

## Installation

Follow these steps to install pC-SAC:

### 1. Setup

First, set the main directory where pC-SAC will be installed. Replace `/path/to/pCSAC` with the actual path on your system.

    export MAIN_DIR="/path/to/pCSAC"

### 2. Download and Install Eigen

Eigen is a dependency for pC-SAC. Download and install it using the following commands:

    wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    tar xvfz eigen-3.4.0.tar.gz
    rm eigen-3.4.0.tar.gz
    cd eigen-3.4.0
    mkdir build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=../ -DINCLUDE_INSTALL_DIR=../include/
    make install
    cd $MAIN_DIR

### 3. Compile pC-SAC

Edit the CMakeLists.txt in the pC-SAC directory to include the path to the Eigen include directory:

    nano ${MAIN_DIR}/CMakeLists.txt
    # Modify the include path to /path/to/pCSAC/eigen-3.4.0/include/

Now, compile pC-SAC:

    mkdir build
    cd build
    cmake ..
    make

## Testing

To verify that pC-SAC has been installed correctly, you can run a test using provided example files:

### Files Description

Within the `test_toy` directory, you will find several example input files required for running pC-SAC. Please explore these files to understand each required input:
Additionally, `input_generation` directory contains Python build functions to generate these files.

- `int_mat_seg_len.txt`: Length (number of rows) of this file represents the length of any chain for a given reconstruction. Each row represents the amount of DNA in each node (i.e. the radius) of a chromatin chain, either on amstrongs or base-pairs.
- `interaction_matrix.txt`: Initial low-resolution probabilities.
- `test.conf`: Configuration file with pC-SAC parameters.

Navigate to the test data directory:

    cd ${MAIN_DIR}/script/test_toy

Set the LD_LIBRARY_PATH to include the Boost library path:

    export LD_LIBRARY_PATH="/nfs/sw/boost/boost-1.72.0/lib"

Load the Boost module:

    module load boost/1.72.0

Run the test:

    ${MAIN_DIR}/build/bin/chromatin.sis.coarse -conf ${MAIN_DIR}/script/test.conf -prefix test_pCSAC
