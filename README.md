# orthpol

A package of routines for generating orthogonal polynomials and Gauss-type quadrature rules

## Description

This is a modern Fortran frontend to the ORTHPOL package by W. Gautschi [1]. 

The original ORTHPOL code can be found [here](https://www.cs.purdue.edu/archives/2001/wxg/codes/ORTHPOL). A copy of the old package for verification purposes can be found in the folder [legacy](./legacy). The folder contains an updated Makefile to compile the original tests with gfortran.


## Building orthpol

1. Download the Fortran source code:
```
git clone https://github.com/ivan-pi/orthpol
```

2. Create a build folder and build the Makefiles using CMake:
```
cd orthpol
mkdir build && cd build
cmake ..
```

3. Run the make script to create the orthpol library:
```
make orthpol
```

## Documentation

Coming soon...

## References

1. Gautschi, W. (1994). Algorithm 726: ORTHPOL–a package of routines for generating orthogonal polynomials and Gauss-type quadrature rules. ACM Transactions on Mathematical Software (TOMS), 20(1), 21–62. doi:10.1145/174603.174605