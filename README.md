### DESCRIPTION
	- we present, the Spatial Contextualization for Closed Itemset Mining (SCIM) algorithm, a mining procedure that builds a metric space for  the targetdatabase in such a way that relevant closed itemsets can beretrieved regarding the relative spatial location of their items, proposed by Mantuan et al. 

### REQUIREMENTS
	- Eigen 3.2.0 or later;
  	- Boost 1.5.0 or later
	- GCC 4.9.1 / MSVC 2013 or later; 

### Building

The basic steps for configuring and building the library look like this:

```bash
$ git clone https://github.com/google/testeSCIM.git
$ mkdir build && cd build
$ cmake -G <generator> [options] ../benchmark
# Assuming a makefile generator was used
$ make
```

### USAGE: 
	SCIM.exe [INPUT_PATH] [DR] [OUTPUT_PATH]

In case of doubts about the code, send it to amantuan@id.uff.br / altobelli.bm@gmail.com

If you use this implementation, please refer to:
```
@inproceedings{MANTUAN2018,
 author = {Altobelli B. Mantuan and Leandro A.F. Fernandes },
 title = {Spatial Contextualization for Closed Itemset Mining},
 booktitle = {Proc. IEEE ICDM},
 year = {1999},
 pages = {507--516},
} 

```
