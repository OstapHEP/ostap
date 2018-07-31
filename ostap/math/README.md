# Math

* [ostap.math](README.md)

Collection of various math classes and utilities 
- `math_ve.py` 
   - `VE`: very powerful `Ostap::Math::ValueWithError` object, a heart of  many other `ostap` modules 
   - zillions of *standard* math-fcuntions that deals with uncertainties 
- `bernstein.py` : definition of *Bernstein* polymonials (or, more correctly, polynomials in *Bernstein form*) and their partners
   - operations with polynomials 
   - conversions to/from *Legendre*, *Chebyshev*, *monomial*, *Newton* and *Lagrange* forms 
   - derivatives and integrals 
   - root finding
   - polynomial (long) division, deflating, the greatest common divisor, the least common multiple, ...
   - *Sturm's* sequence
   - convex hulls 
- `derivative.py`     : numerical differentiation and related functions 
- `interpolation.py`  : powerful utilities for interpolation
- `integral.py`       : numerical integrations and related functions. Mainly as a replacement for several functions from  `scipy`, when/if `scipy` is not accessible. 
- `geometry.py`       : easy-to-use but powerful collection of functions for operations with 3D-points, vectors, lines and planes    
- `kinematics.py`     : easy-to-use but powerful collection of functions for operations with 4D-vectors 
- `linalg.py`         : easy-to-use but powerful collection of functions for operations with matrices and vectors  
- `rootfinder.py`     : root finder. Mainly as a replacement for `scipy.optimise.brentq` function from when/if `scipy` is not accessible. 
- `models.py`         : decorations for many classes/models from `Ostap::Math`-namespace 
- `param.py`          : utilities for parameterization of fnuction (or historgams) in terms of *Legendre*, *Chebyshev*, *Fourier* , *Fourier/Cosine* , *Bernstein/Bezier* sums 
- `base.py`           : few basis decorations, including transparent communication of `complex` and `std::complex<double>` types 

