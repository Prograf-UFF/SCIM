# Spatial Contextualization for Closed Itemset Mining
The Spatial Contextualization for Closed Itemset Mining (SCIM) algorithm is a mining procedure that builds a metric space for the target database in such a way that relevant closed itemsets can be retrieved regarding the relative spatial location of their items.

The SCIM algorithm uses [Dual Scaling](https://www.taylorfrancis.com/books/9781317781943) to map the items of the database to a multidimensional metric space called Solution Space. The representation of the database in the Solution Space assists in the interpretation and definition of overlapping clusters of related items. The distances of the items to the centers of the clusters are used as criteria for generating itemsets. Therefore, instead of using the minimum support threshold, a distance threshold is defined concerning the reference and the maximum distances computed per cluster during the mapping procedure.

The approach was developed by Altobelli B. Mantuan and [Leandro A. F. Fernandes](http://www.ic.uff.br/~laffernandes). Check out the [project's website](http://www.ic.uff.br/~laffernandes/projects/sodm) for details.

This repository includes the C++ implementation of the SCIM algorithm, and a sample application using this implementation.

Please cite our IEEE ICDM 2018 paper if you use this code in your research:
```
@InProceedings{mantuan_fernandes-icdm-2018,
  author    = {Mantuan, Altobelli B. and Fernandes, Leandro A. F.},
  title     = {Spatial contextualization for closed itemset mining},
  booktitle = {Proceedings of the 2018 IEEE International Conference on Data Mining (ICDM)},
  year      = {2018},
  pages     = {to appear},
  url       = {http://www.ic.uff.br/~laffernandes/projects/sodm},
}
```

Do not exitate to contact Altobelli B. Mantuan ([amantuan@ic.uff.br](mailto:amantuan@ic.uff.br), [altobelli.bm@gmail.com](mailto:altobelli.bm@gmail.com)) if any problems are encountered.


## Licence
All code is released under the [GNU General Public License](https://www.gnu.org/licenses/), version 3, or (at your option) any later version.


## Platforms
We have compiled and tested the sample application on Linux and Windows using GCC 4.9.1 and Microsoft Visual C++ 2013.


## Requirements
Make sure that you have all the following tools and libraries installed and working before attempting to compile the SCIM implementation.

Required tools:
- [GCC](https://pt.wikipedia.org/wiki/GNU_Compiler_Collection) 4.9.1 or later (Linux or Windows)
- [Microsoft Visual C++](https://en.wikipedia.org/wiki/Microsoft_Visual_C%2B%2B) 2013 or later (Windows)
- [CMake](https://cmake.org)

Required C++ libraries:
- [Boost](https://www.boost.org) 1.5.0 or later (header-only libraries)
- [Eigen](https://eigen.tuxfamily.org) 3.2.0 or later


## Building, Compiling, and Running
Use the [git clone](https://git-scm.com/docs/git-clone) command to download the project:
```bash
$ git clone https://github.com/Prograf-UFF/SCIM.git SCIM
$ cd SCIM
```

Make sure you have the environment variables for Eigen (``EIGEN3_INCLUDE_DIR``) and Boost (``BOOST_ROOT``) defined in your system.

The basic steps for configuring and building the sample application look like this:
```bash
$ mkdir build
$ cd build

$ cmake [-G <generator>] [options] -DCMAKE_BUILD_TYPE=Release ..
```

Assuming a makefile generator was used:
```bash
$ make
```

To run the sample application, just call:
```bash
$ SCIM <database-file-path> <dr-threshold-value> <output-folder-path>
```

The database file must follow the ``.num`` format used by [The LUCS-KDD Discretised/normalised ARM and CARM Data Library](http://www.csc.liv.ac.uk/~frans/KDD/Software/LUCS_KDD_DN).

The distance ratio threshold (``dr-threshold-value``) must be in the \[0, 1\] range. We believe that 0 (zero) is an excellent initial guess value. The user may increase the parameter value slightly in an exploratory fashion in order to detect more closed itemsets.
