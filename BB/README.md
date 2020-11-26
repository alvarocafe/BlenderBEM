# BEM base
---
A Julia language implementation of the a fast isogeometric boundary element method formulation (IGABEM) with Hierarchical Matrices partitioning and approximation using Lagrange polynomials and the Adapative-Cross Approximation (ACA). 
The goal of the BEM_base project is to provide a platform that can be used alongside FreeCad and Gmsh to quickly solve bidimensional boundary element problems.
Initially, this implementation will solve the Helmholtz and Laplace equations.

To solve a problem with a geometry file 'file.msh' and boundary conditions in each face described by an array 'BCFace', to use the BEM_base to solve the problem, simply run:  
    `include("BEM_base.jl")`

    T,qT = BEM_base("file.msh",BCFace, k, "heat") # For solving the Laplace equation

    phi,qphi = BEM_base("file.msh",BCFace, k, "wave") # For solving the Helmholtz equation

where k is the thermal condutivity in the Laplace equation and the wavenumber in the Helmholtz equation.
BEM_base will analyse the mesh and choose the solver accordingly.

## Installation
---
The project runs on Julia 0.6.4. To include the dependecies, simply clone the repository and include the files.  

    git clone https://github.com/alvarocafe/BEM_base
And run the code given above.

## Contribute
---
This project is still in a very early stage and contributions are welcome.

## NOTICE
---
The program is not yet usable in its current form, I'll keep updating it and I'll update this notice when the program is usable. You can still use the subroutines to build BEM models and run them, but that's harder than the intended usage of the package. 

## License
---
This project is currently licensed under GNU GPL v.3
