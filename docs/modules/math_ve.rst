.. _ve:

.. include:: ../references.rst

:mod:`ostap.math.ve` module
===========================

      
This module decorates the C++ class :class:`Ostap::Math::ValueWithError`, referred hereafter to as ``VE`` :

- adding nice printout 
- comparison functions, compare between ``VE`` objects and simple numerical types
- random number generation 

The main object
---------------

.. class:: VE

   Python alias for C++ class ``Ostap::Math::ValueWithError``, that correspons to a very 
   simple notion of the value and the associated uncertainty

   .. code-block::

      v = VE ( 1 , 0.1**1 )
      print ( 'VE           :', v ) 
      print ( 'Value        :', v.value ()  ) 
      print ( 'Error        :', v.error ()  ) 
      print ( 'Cov2/error^2 :', v.cov2  ()  ) 

   .. method:: b2s ()   

      Calculate the "effective" background-to-signal ratio from the value and the associated uncertainty 
      using the identity :math:`\frac{\sigma(S)}{S} = \frac{\sqrt{S}}{S}\sqrt{1+\frac{B}{S}}`. 
      From this identity one gets :math:`\left.\frac{B}{S}\right|_{\mathrm{eff}} \equiv \frac{\sigma^2(S)}{S}-1`:
      
      .. code-block::

         v = VE( ... )
         print('B/S=', v.b2s())

   .. method:: purity ()   

      Calculate the "effective purity" ratio using the 
      identity :math:`p_{\mathrm{eff}} = \frac{S}{S+B} = \frac{1}{1+\frac{B}{S}}`, 
      where  the effective "background-to-signal" ratio is estimated as 
      :math:`\left.\frac{B}{S}\right|_{\mathrm{eff}} = \frac{\sigma^2(S)}{S} -1`.
      Finally one gets :math:`p_{\mathrm{eff}} \equiv \frac{S}{\sigma^2(S)}`

      .. code-block::

         v = VE( ... )
         print('purity=', v.purity())
      
   .. method:: precision ()   

      Get the precision with "some" error estimate:

      .. code-block::
         
         v = VE( ... )
         print('precision=', v.precision())
         
   .. method:: minmax( n = 1 ) 

      Get an easy and coherent way to access ``min/max`` for the object: 
      :math:`v_{min}=v-n\sigma`, :math:`v_{max}=v+n\sigma`:
      
      .. code-block::

         v = VE ( ... ) 
         vmin  , vmax  = v.minmax()     #: :math:`\pm1\sigma`
         vmin3 , vmax3 = v.minmax( 3 )  #: :math:`\pm2\sigma`

   .. method:: gauss ( accept = lambda s : True , nmax = 1000 ) 
   
      Get the Gaussian random numbers.

      .. code-block::

         v = VE ( ... ) 
         for i in range(10) : print ( v.gauss() ) ##   get 10 random numbers 

      Get only positive Gaussian random:

      .. code-block::

         v = VE ( ... ) 
         for i in range(10) : print ( v.gauss ( accept = lambda s : s > 0 ) ) 
      

   .. method:: poisson ( fluctuate , accept = lambda s : True ) 

      Get Poisson random number accoring to parameters
      
      Use the value as the expected valur :math:`\mu` for Poisson distribution
      .. code-block:: 

         v = VE ( ... ) 
         for i in  range ( 10 ) : print ( v.poisson( fluctuate = False ) )
      
      Fluctuate the mean value witjng Gaussian uncerttainrty:
      .. code-block:: 

         v = VE ( ... ) 
         for i in  range ( 10 ) : print ( v.poisson( fluctuate = True  ) )
         
      Accept only even numbers:

      .. code-block:: 

         v = VE ( ... ) 
         even = lambda s : 0 == s%2
         for i in  range ( 10 ) : print ( v.poisson( fluctuate = True  , accept = even ) )


.. toctree::

   :maxdepth: 2

   :caption: Classes:

   ../classes/ValueWithError
   ../classes/OstapMath
