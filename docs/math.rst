.. _math:

.. include:: ./references.rst

math
====

Collection of various math classes and utilities 


* ```ve_`` and ``math_ve`` 

  * `VE`: very powerful ``Ostap::Math::ValueWithError`` object, a heart of  many other `ostap` modules 
  * zillions of *standard* math-functions that deals with uncertainties 

* ``base``      : few basis decorations, including transparent communication of ``complex`` and ``std::complex<double>`` types 

* ``bernstein`` : definition of *Bernstein* polymonials (or, more correctly, polynomials in *Bernstein form*\ ) and their partners

  * operations with polynomials 
  * conversions to/from other polynomials forms
      + *Legendre*\ 
      + *Chebyshev*\ 
      + *monomial*\ 
      + *Newton* (interpolation) 
      + *Lagrange* (interpolation)  
  * derivatives and integrals 
  * root finding
  * polynomial algebra
      + (long) division
      + deflating
      + the greatest common divisor
      + the least common multiple, ...
  * *Sturm's* sequence
  * convex hulls 

* ``bspline``        : operations with *BSpline*
* ``dalitz``         : useful fnuctions for manipulations with `Dalitz plot variables`
* ``derivative``     : numerical differentiation and related functions 
* ``geometry``       : easy-to-use but powerful collection of functions for operations with 3D-points, 3D-vectors, 3D-lines and 3D-planes    
* ``integral``       : numerical integrations and related functions. Mainly as a replacement for several functions from  scipy_ , when/if scipy_ is not accessible. 
* ``interpolation``  : powerful utilities for interpolation
* ``kinematics``     : easy-to-use but powerful collection of functions for operations with 4D-vectors 
* ``linalg``         : easy-to-use but powerful collection of functions for operations with matrices and vectors  
* ``minimize``       : home-made replacement for ``scipy.minimize.minimize_scalar`` when/if scipy_ is not accessible
* ``models``         : collection of properly decorated C++ functions, models, pdfs, parameterizations from ``Ostap::Math``\-namespace 
* ``operations``     : collection of simple technical classes to wrap the basic math operations for callables
* ``other``          : collection of other, non-classified math-related classes & functions
* ``param``          : utilities for parameterization of fnuction (or historgams) in terms of *Legendre*\ , *Chebyshev*\ , *Fourier* , *Fourier/Cosine* , *Bernstein/Bezier* sums 
* ``primes``         : utilities for dealing with prime numbers 
* ``random_ext``     : extewnsionof ``python.random`` module 
* ``rootfinder``     : root finder. Mainly as a replacement for ``scipy.optimise.brentq`` function from when/if scipy_ is not accessible. 



.. toctree::
   :maxdepth: 2
   :caption: Modules

   modules/math_ve

   :caption: classes 

   classes/ValueWithError
   classes/OstapMath


