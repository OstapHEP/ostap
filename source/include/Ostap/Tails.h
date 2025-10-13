// ============================================================================
#ifndef OSTAP_TAILS_H 
#define OSTAP_TAILS_H 1
// ============================================================================
//  Include  files 
// ============================================================================
// STD & STL 
// ============================================================================
#include <cmath>
// ============================================================================
/** @file Ostap/Tails.h
 *  Tails for CrystalBall-like functions 
 *  @see Ostap::Math::CrystalBall
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2010-04-19
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /** @class Tail
     *  Representanton of the tail for CrystalBall-like functions 
     */
    class Tail
    {
      // =====================================================================
    public : 
      // =====================================================================
      /** \f$ n \rightarrow N \f$ transformation
       *
       *  It is done to ensure that the actual-power-law exponent 
       *  does not excees 1 
       *
       *  @param    n (input) n-paameter (external) 
       *  @return transformed N-parameter (internal)
       */
      static double N ( const double n ) ;
      // =====================================================================      
    public :
      // =====================================================================
      /** Tail parameters 
       *  @param alpha alpha-parameter
       *  @param n     n-parameter
       */
      Tail
      ( const double alpha = 2 ,
	const double n     = 1 ) ;
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// tail parameter alpha 
      inline double alpha  () const { return m_alpha   ; }
      /// tail parameter alpha squared 
      inline double alpha2 () const { return m_alpha * m_alpha  ; }
      /// tail n-parameter n     (external) 
      inline double n      () const { return m_n       ; }
      /// true tail N-parameter (internal ) 
      inline double N      () const { return N ( m_n ) ; } // internal N parameter      
      // ======================================================================
    public:  // setters 
      // ======================================================================
      /// set alpha-parameter
      bool setAlpha ( const double value ) ;
      /// set n-parameter
      bool setN     ( const double value ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private :
      // ======================================================================
      /// tail parameter alpha 
      double m_alpha { 2 } ; // tail parameter alpha 
      /// tail parameter n 
      double m_n     { 1 } ; // tail parameter n       
      // ======================================================================
    } ; //                                   The end of class Ostap::Math::Tail
    // ========================================================================
    /** @class LeftTail 
     *  The left tail of Crystal Ball-like function:
     *  \f$ f ( x ) \equiv F \left( 1 - \frac{F^\prime}{F}\frac { x - x_0 } { N } \right)^{-N} \f$
     *  such as  
     *  - \f$ f ( x_0 )       = F        \f$ 
     *  - \f$ f^\prime ( x_0 ) = F^\prime \f$ 
     *  - \f$ N = \sqrt { 1 + n^2 } \f$ 
     */
    class LeftTail : public Ostap::Math::Tail 
    {
      // ======================================================================
    public:
      // ======================================================================
      /** Tail parameters 
       *  @param alpha alpha-parameter
       *  @param n     n-parameter
       */
      LeftTail
      ( const double alpha = 2 ,
	const double n     = 1 ) ;
      // ======================================================================
      /// Tail parameters 
      LeftTail ( const Ostap::Math::Tail& tail ) ;
      // ======================================================================
    public:
      // ======================================================================
      /** evaluate the (left) tail function
       *  @param x    the x point 
       *  @param x0   normalization point \f$ f(x) \equiv 0 \f$ for \f$ x>x_0\f$ 
       *  @param F    function value \f$ f(x_0) \f$ at normalization point 
       *  @param dFoF value of log-dervative \f$ \frac{f^\prime(x_0)}{f(x_0)}\f$  at normalization point 
       */
      double evaluate
      ( const double x      ,
	const double x0     ,
	const double F      ,
	const double dFoF   ) const ;
      // ======================================================================
      /** evaluate the (left) tail function
       *  @param x    the x point 
       *  @param x0   normalization point \f$ f(x) \equiv 0 \f$ for \f$ x>x_0\f$ 
       *  @param F    function value \f$ f(x_0) \f$ at normalization point 
       *  @param dFoF value of log-dervative \f$ \frac{f^\prime(x_0)}{f(x_0)}\f$  at normalization point 
       */
      inline double operator() 
      ( const double x      ,
	const double x0     ,
	const double F      ,
	const double dFoF   ) const { return evaluate ( x , x0 , F , dFoF ) ; } 
      // ======================================================================
    public:
      // ======================================================================
      /** get the integral of power-law function
       *  @paral low  low  integral edge 
       *  @paral high high integral edge
       *  @param x0   normalization point \f$ f(x) \equiv 0 \f$ for \f$ x>x_0\f$ 
       *  @param F    function value \f$ f(x_0) \f$ at normalization point 
       *  @param dFoF value of log-dervative \f$ \frac{f^\prime(x_0)}{f(x_0)}\f$  at normalization point 
       */
      double integral
      ( const double low  ,
	const double high ,
	const double x0   ,
	const double F    ,
	const double dFoF ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    } ; //                               The end of class Ostap::Math::LeftTail
    // ========================================================================
    /** @class RightTail 
     *  The left tail of Crystal Ball-like function:
     *  \f$ f ( x ) \equiv F \left( 1 - \frac{F^\prime}{F}\frac { x - x_0 } { N } \right)^{-N} \f$
     *  such as  
     *  - \f$ f ( x_0 )       = F        \f$ 
     *  - \f$ f^\prime ( x_0 ) = F^\prime \f$ 
     *  - \f$ N = \sqrt { 1 + n^2 } \f$ 
     */
    class RightTail: public Ostap::Math::Tail 
    {
      // ======================================================================
    public:
      // ======================================================================
      /** Tail parameters 
       *  @param alpha alpha-parameter
       *  @param n     n-parameter
       */
      RightTail
      ( const double alpha = 2 ,
	const double n     = 1 ) ;
      // ======================================================================
      /// Tail parameters 
      RightTail ( const Ostap::Math::Tail& tail ) ;
      // ======================================================================
    public:
      // ======================================================================
      /** evaluate the (right) tail function
       *  @param x    the x point 
       *  @param x0   normalization point \f$ f(x) \equiv 0 \f$ for \f$ x<x_0\f$ 
       *  @param F    function value \f$ f(x_0) \f$ at normalization point 
       *  @param dFoF value of log-dervative \f$ \frac{f^\prime(x_0)}{f(x_0)}\f$  at normalization point 
       */
      double evaluate
      ( const double x      ,
	const double x0     ,
	const double F      ,
	const double dFoF   ) const ;
      // ======================================================================
      /** evaluate the (left) tail function
       *  @param x    the x point 
       *  @param x0   normalization point \f$ f(x) \equiv 0 \f$ for \f$ x<x_0\f$ 
       *  @param F    function value \f$ f(x_0) \f$ at normalization point 
       *  @param dFoF value of log-dervative \f$ \frac{f^\prime(x_0)}{f(x_0)}\f$  at normalization point 
       */
      inline double operator() 
      ( const double x      ,
	const double x0     ,
	const double F      ,
	const double dFoF   ) const { return evaluate ( x , x0 , F , dFoF ) ; } 
      // ======================================================================
    public:
      // ======================================================================
      /** get the integral of power-law function
       *  @paral low  low  integral edge 
       *  @paral high high integral edge
       *  @param x0   normalization point \f$ f(x) \equiv 0 \f$ for \f$ x<x_0\f$ 
       *  @param F    function value \f$ f(x_0) \f$ at normalization point 
       *  @param dFoF value of log-dervative \f$ \frac{f^\prime(x_0)}{f(x_0)}\f$  at normalization point 
       */
      double integral
      ( const double low  ,
	const double high ,
	const double x0   ,
	const double F    ,
	const double dFoF ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    } ; //                              The end of class Ostap::Math::RightTail
    // ========================================================================    
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_TAILS_H
// ============================================================================
