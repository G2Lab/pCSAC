# pC-SAC 


This guide provides instructions for installing pC-SAC. 
pC-SAC is a tool that uses adaptive importance sampling technique with sequential Monte Carlo to generate a set of three-dimensional chromatin chains to enhance low-resolution Hi-C data. 

## Installation via Conda

pC-SAC is available on conda by executing the following command.

    conda install pcsac

## Build locally   
### Prerequisites

Before installing pC-SAC, ensure that the following dependencies and software are installed on your system:

- GCC compiler (12.2)
- Python (> 3.10)

- CMake (>= 2.6)
- Boost (>= 1.81)
- Eigen (= 3.4)


These can typically be loaded on an HPC environment using module commands, like so:

    module load cmake
    module load boost
    module load eigen

### Installation

Follow these steps to install pC-SAC:

#### 1. Setup

Clone this repository:

    git clone https://github.com/G2Lab/pCSAC.git

    MAIN_DIR="/path/to/pCSAC"

    cd ${MAIN_DIR}

Set the main directory where pC-SAC will be installed. Replace `/path/to/pCSAC` with the actual path on your system.

    export MAIN_DIR="/path/to/pCSAC"

#### 2. Download and Install Eigen

Eigen is a dependency for pC-SAC. Download and install it using the following commands:

    wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    tar xvfz eigen-3.4.0.tar.gz
    rm eigen-3.4.0.tar.gz
    cd eigen-3.4.0
    mkdir build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=../ -DINCLUDE_INSTALL_DIR=../include/
    make install
    cd ${MAIN_DIR}

#### 3. Compile pC-SAC

Edit the CMakeLists.txt in the pC-SAC directory to include the path to the Eigen include directory:

    nano ${MAIN_DIR}/CMakeLists.txt
    
    # Modify the include path to /path/to/pCSAC/eigen-3.4.0/include/

Now, compile pC-SAC:

    cd ${MAIN_DIR}
    mkdir build
    cd build
    cmake ..
    make

## Toy Testing

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

## Configuration File 

In order to execute pC-SAC the user must provide a configuration file with the required information to generate a chromatin ensemble.

- out_path: Path where chromatin chains will be stored

- start_end_file: This file defines the number of nodes and the radius of each node for the reconstruction of a genomic locus.
The file specifies the total number of nodes \(N\) and the radius \(r\) of each node. The number of nodes \(N\) depends on the total length \(L\) of the genomic locus being reconstructed and the desired resolution \(R\) for each node. The radius \(r\) of each node corresponds to the physical size of the node.

    Example:

    1. **Number of Nodes (N):**  
        **N = L / R**

    - \(N\): Number of nodes (monomers)  
    - \(L\): Total length of the genomic locus (e.g., 2 Mb = 2,000 kb)  
    - \(R\): Resolution or size of each node (e.g., 5 kb)  

    2. **Radius of Each Node (r):**  

        **r = R * k**

    - \(r\): Radius of each node (e.g., in Angstroms)  
    - \(R\): Resolution or size of each node (e.g., 5 kb)  
    - \(k\): Conversion factor (e.g., 1 kb = 10 Angstroms)

    For a target genomic locus of **2 Mb** with a desired resolution of **5 kb** for each node, the calculations are as follows:

    1. **Number of Nodes (N):**  
        **N = (2,000 kb) / (5 kb) = 400**
    
        Therefore, there will be **400 nodes** per chain in total.

    2. **Radius of Each Node (R):**  
    If 1 kb = 10 Angstroms, then:  
    **r = 5 kb × 10 A = 50 A**

    Therefore, the radius of each node is **50 A**.

    This will generate a start_end_file with 400 lines, where each line represents a node, and each node is assigned a value of 50.


- pval_file: Path to the file containing all projected restrictions (see **generate_input** module for more details). This file maps low-resolution interactions onto high-resolution bin coordinates. 

    When an interaction occurs between bins \(j\) and \(k\) in a low-resolution matrix, it is projected onto the high-resolution matrix using the following approach:  
    - Each low-resolution bin \(j\) is expanded into a range of high-resolution bins from:  

        (j * r + 1) to  (j * r + r)

    - Similarly, bin \(k\) is expanded from:  

        (k * r + 1) to  (k * r + r)

    Here, \(r\) represents the ratio between the low-resolution bin size and the high-resolution bin size.  


    If the low-resolution bin size is **20 kb** and the high-resolution bin size is **5 kb**, then:  

        r = 20/5 = 4

    This means each 20 kb bin corresponds to four 5 kb bins.  

    An interaction between bins \(j\) and \(k\) in the 20 kb matrix will be projected to include all interactions (and corresponding probabilities) between the high-resolution bin ranges:  

        [j * 4 + 1, j * 4 + 4] and [k * 4 + 1, k * 4 + 4]


- collision_length: Defines the minimum allowable distance between two nodes to prevent overlap or collisions during chain generation. It is defined based on the previously defined radius.

- interaction_distance: Specifies the maximum distance within which two nodes can interact or be considered in proximity for calculating interactions.  This distance can be calculated as: **d=r*2** where r is the radius of the monomers.

- packing_density: Defines the proportion of the total available space within the nucleus that is occupied by the chain. A value of **0.1** indicates that 10% of the nucleus volume is occupied with the chain, reflecting its compactness within the nuclear environment.

- nucleus_sphere_diameter = Specifies the diameter of the spherical nucleus for a chain to be generated. A value of **6,000,000,000,000,000** represents the total distance across the nucleus where nodes can be positioned, providing a spatial constraint for chain generation.

- M_max: The desired number of chains to be generated in the ensemble. (see Supplementary Figure S11E for recommended values).

- number_sample_points: The lattice value (s) representing the number of possible positions (\(m_j = \{a_j, b_j, c_j\}\)) where the \(j^\text{th}\) node can be added. This ensures the selected constraints are satisfied with respect to the partially generated chain \(B_x^{t-1}\) (see Supplementary Figure S11C for recommended values). 

- monte_carlo_coefficient: The coefficient used to control the weight of random sampling in the Monte Carlo process. Since random sampling is used to generate chain configurations, this parameter ensures that the generated samples adhere to the desired probability distribution. A value of **1.0** indicates a standard weighting

- rho_1: Represents the weight assigned to the collision constraint term  f_1(x_t^{(l)}) in the minimizing function. Controls how strongly collisions are penalized during optimization. A higher value  (**1.0**) increases the importance of avoiding collisions.

- rho_2: Represents the weight assigned to the distance constraint term  f_2(x_t^{(l)}) in the minimizing function.Controls how strongly the desired interaction distances between nodes are enforced. A higher value (**1.0**) makes the system prioritize satisfying interaction distances.

- tau_t: A temperature-like parameter that influences the overall contribution of the weighted terms to the minimization function.
    - At higher values of  \tau, the function smoothens the penalties, allowing for more exploration and less strict adherence to constraints.

	- At lower values, the function sharpens, penalizing violations more strongly and favoring convergence to optimal solutions.

    (see Supplementary Figure S11C for suggested values)


## STDOUT 
For each generated chain, pC-SAC reports the following key metrics:

* Number of Corrections: The number of times the chain generation process required structural adjustments to maintain consistency with the expected distribution.

* Error: Represents the proportion of missed interactions in the generated chain relative to the total expected interactions in the trial distribution. This metric reflects the chain’s accuracy.

* LogWeight: Statistical score reflecting the likelihood or importance of the generated chain.

* Iterations per Node: The ratio of the total number of iterations performed to the number of nodes in the generated chain. 

