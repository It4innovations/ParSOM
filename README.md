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

mpirun -n 20 som_linux -c <file> -o <file> -i <input>

**List of basic parameters:**
```
-c <name> - name of configuration file\
-o <name> - name of output file\
-i <name> - name of input file\
-debug    - enable debug reports (Default false)\
-info     - continuously writes status of calculation (Default false)\
-t <number>      - number of thread (Default 1)\
-TC              - Enabled hybrid learning (Default false) \
-somHybrid <range 0 - 100> \
-hybridInc        \
-hybridDec        \
-hybridStep <number>\
-cosin           - Use cosine similarity instead of euclidean distance (Default false)\
-version <number> - Version type of input file (Default 2)\
-update <number>  - Version of update (Default 0)\
-onOutput - Turn off any output file (Default false)
```
>Hybrid learning is decribe [here](https://www.scopus.com/record/display.uri?eid=2-s2.0-84959217069&origin=resultslist&sort=plf-f&src=s&sid=a17258611de1fb30fefc2aeb51f1fa5d&sot=autdocs&sdt=autdocs&sl=18&s=AU-ID%2853985639600%29&relpos=4&citeCnt=2&searchTerm=)

**Input file example**

Each record is on a separate line. 

`<position>:<value> <position>:<value> ...`
```
0:0.222222222 1:0.625 2:0.06779661 3:0.041666667
0:0.166666667 1:0.416666667 2:0.06779661 3:0.041666667
0:0.111111111 1:0.5 2:0.050847458 3:0.041666667
0:0.083333333 1:0.458333333 2:0.084745763 3:0.041666667
```

**Configuration file example**

`<x dimension>  <y dimension> <number of iteration> <number of input records> <default name of output files> <version>`

**Output file example**

Each record describes one neuron.

`0:<value> 1:<value> ... `
