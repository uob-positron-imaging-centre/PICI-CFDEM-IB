# UoB Positron Imaging Centre's Improved CFDEM Distribution

This represents a modified iteration of CFDEM, featuring an enhanced fully resolved Computational Fluid Dynamics-Discrete Element Method (CFD-DEM) solver.

**Email**: <chehanqiao@outlook.com>

The solver has been meticulously crafted by integrating the CFDEM framework with the established methodology introduced by C. Zhang (2018). You can access C. Zhang's code [https://github.com/ChenguangZhang/sdfibm]. Our team has thoroughly examined his algorithm, and we have incorporated what we consider to be a more accurate and refined solution.

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
Compatibility verified with OpenFOAM-5.x exclusively. For installation instructions, please refer to the <www.cfdem.com> website. Simply substitute the existing 'CFDEMcoupling-PUBLIC/' with this version and initiate the solver build by executing the command 'cfdemCompCFDEM'
## Examples:
* **Flow through the ordered packings**  
A validation case has been added in the folder ./tutorials/. The calculated fluid-particle interaction force can be validated using the existing analytical data by Zick&Homsy (1995). In the tutorial case, simple cubic (SC), body centred cubic (BCC) and face centred cubic (FCC) configurations over a range of overall solid packing fracitons were considered. Contact me for more tutorial cases.
## License:

GPL-3.0

## References:
C. Zhang, C. Wu, K. Nandakumar, Effective Geometric Algorithms for Immersed Boundary Method Using Signed Distance Field, Journal of Fluids Engineering, 141 (2018). 

T. Kempe, S. Schwarz, J. Fröhlich, Modelling of spheroidal particles in viscous flows,  Proceedings of the Academy Colloquium Immersed Boundary Methods: Current Status and Future Research Directions (KNAW, Amsterdam, The Netherlands, 15–17 June 2009), 2009.

Zick, A., & Homsy, G.  Stokes flow through periodic arrays of spheres. Journal of Fluid Mechanics, 115(1982), 13-26. doi:10.1017/S0022112082000627.

## Citation

If you use PICI-CFDEM-IB in a publication, please cite the following article (the full text is freely accessible):

> Hanqiao Che et al._ "**A novel semi-resolved CFD-DEM method with two-grid mapping: methodology and validation.**" AIChE Journal, 2023.
> [https://doi.org/10.1002/aic.18321].

Also, we ask that you "star" :star: this repository to help us gauge the level of interest in this project. Your star not only signifies your support but also assists us in tracking user engagement, a crucial factor for securing future grant funding. Your contribution in this regard is highly appreciated.


