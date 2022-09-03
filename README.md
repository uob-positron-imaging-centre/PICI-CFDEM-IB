# UoB Positron Imaging Centre's Improved CFDEM Distribution

This is a modified version of CFDEM with an improved fully resolved CFD-DEM solver.

**Email**: <chehanqiao@outlook.com>

The solver is developed by combining the framework of CFDEM with the existing methodology proposed by C. Zhang (2018). C. Zhang's code can be found [here](https://github.com/ChenguangZhang/sdfibm). We have carefully considered his algorithm and have adopted what we believe to be a more correct solution.

## New submodels :
### **Force model**
* **ReactionIB**  
  The fluid-particle interaction force is calculated based on the momentum source term. An alternative to Model ShirgaonkarIB.
### **Voidfraction model**
* **IBVoidFractionKempe**
* **IBVoidFractionSDF**
## Solver:
* **cfdemSolverIBPICI**
## Installation:
Tested with OpenFOAM-5.x only. Details on installation are given on the <www.cfdem.com> website. Just replace the CFDEM folder with this one and build the solver by runing `cfdemCompCFDEM`.
## Examples:
* **Flow through the ordered packings**  
A validation case has been added in the folder ./tutorials/. The calculated fluid-particle interaction force can be validated using the existing analytical data by Zick&Homsy (1995). In the tutorial case, simple cubic (SC), body centred cubic (BCC) and face centred cubic (FCC) configurations over a range of overall solid packing fracitons were considered.
## License:

GPL-3.0

## References:
C. Zhang, C. Wu, K. Nandakumar, Effective Geometric Algorithms for Immersed Boundary Method Using Signed Distance Field, Journal of Fluids Engineering, 141 (2018). 

T. Kempe, S. Schwarz, J. Fröhlich, Modelling of spheroidal particles in viscous flows,  Proceedings of the Academy Colloquium Immersed Boundary Methods: Current Status and Future Research Directions (KNAW, Amsterdam, The Netherlands, 15–17 June 2009), 2009.

Zick, A., & Homsy, G.  Stokes flow through periodic arrays of spheres. Journal of Fluid Mechanics, 115(1982), 13-26. doi:10.1017/S0022112082000627.


