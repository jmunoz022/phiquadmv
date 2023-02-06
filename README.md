phiquadmv
============

Matlab software for computing actions of $\varphi$-functions of Kronecker sum matrices
======================================================================================

Dependencies
------------

* Matlab >= r2021b
* Chebfun (available for free at https://www.chebfun.org/)

Installation
------------

Just make sure that chebfun and the contents of the src folder are in the Matlab path.

Usage
------
In the examples folder: 
* HochOster_2D_gauss.m --> Script to compute the Hochbruck-Osterman equation in 2D+time with three different Exponential Runge-Kutta methods (orders 1 to 3) and check the convergence of the error in time. 
* AllenCahn_2D_gauss.m --> Script to compute the 2D+time Allen-Cahn equation with a star-shaped initial condition with the Exponential Euler method.


References
----------

M. Croci & J. Muñoz-Matute, Exploiting Kronecker structure in exponential integrators: fast approximation of the action of $φ$-functions of matrices via quadrature, ArXiV preprint (https://arxiv.org/abs/2211.00696), 2022.

Acknowledgements
----------------

This material is based upon work supported by the Department of Energy, NNSA under Award Number DE-NA0003969 and also by the European Union’s Horizon 2020 research and innovation program under the Marie Sklodowska-Curie individual fellowship No. 101017984 (GEODPG).
