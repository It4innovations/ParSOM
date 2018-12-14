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

**Usage**

mpirun -n 20 som_linux -c <file> -o <file> -i <input>

List of basic parameters:

 -c <name> - name of configuration file\
 -o <name> - name of output file\
 -i <name> - name of input file\
 -debug    - enable debug reports (Default false)\
 -info     - continuously writes status of calculation (Default false)\
 -t <number>      - number of thread (Default 1)\
 -somHybrid <range 0 - 100>\
 -hybridInc        -\
 -hybridDec        -\
 -hybridStep <number>\
 -cosin           - Use cosin metrik instead of euklid (Default false)\
 -version <number> - Version type of input file (Default 2)\
 -update <number>  - Version of update (Default 0)\
 -onOutput - Turn off any output file (Default false)


