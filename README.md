# ParSOM
Parallel Self-organizing map

ParSOM is developed at Czech national supercomputing centre IT4Innovations. The goal of this application is to contribute to the area of Artificial Neural networks by speed up computing by using High Performance Computing and improvements for sparse data with higher dimensions. 

The code has been developed and optimized for the [IT4Innovations Salomon cluster](https://docs.it4i.cz/salomon/hardware-overview/).

# Dependence
1. GCC/6.3.0-2.27
2. OpenMPI/1.10.7


## Build instructions
1. Clone the repository
2. Run make

After successfull compilation the `run` directory should contain one executable 'som_linux'.

## Usage


Parameters:
-c <name> - configuration file
