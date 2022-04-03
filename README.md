# Smooth Finite Element Modeling of Poration in Lipid Membranes

## Summery
Pore formation in lipid membranes is a natural phenomenon that can be triggered by various chemical or physical cues such as detergents, electric fields (electroporation), illumination, and mechanical loads. Under normal physiological conditions, transient pores can open in the boundaries of mammalian cells such as cardiac muscle cells, and they may lead to cell death. On the other hand, controlled pore formation in lipid membranes has important biomedical applications, including drug delivery, gene therapy, and cancer treatment. Therefore, it is crucial to understand the mechanisms that lead to pore nucleation and evolution in lipid membranes.
Equilibrium configurations of lipid membranes can be accurately approximated by using the so-called curvature models proposed by [Helfrich](https://www.researchgate.net/publication/281560666_Elastic_Properties_of_Lipid_Bilayers_Theory_and_Possible_Experiments). The free energy functional of an open membrane can be found [here](https://link.springer.com/article/10.1007/s00285-007-0118-2).


## The Problem & the Solution Method
The program provided here allows one to capture the shape evolution of closed membranes in the presence of Gaussian energy. Moreover, one can investigate the closure of a pore on the boundary of lipid membranes, which is controlled by the line tension energy of membranes, and see the effect of Gaussian curvature on the smooth transformation of an open membrane to a closed one. 

This program employs a C<sup>1</sup> conforming finite element to discretize the nonlinear equations and uses the standard Newton's method to solve the problem. The time-consuming step in Newton's method is solving the linear system J(x<sub>i</sub>) &Delta;x = R(x<sub>i</sub>), when the Jacobian matrix is not sparse, which is the case here due to the penalty and regularization terms in the energy functional of open lipid membranes. To overcome this issue, we use the [Sherman-Morrison formula](https://en.wikipedia.org/wiki/Sherman%E2%80%93Morrison_formula).



