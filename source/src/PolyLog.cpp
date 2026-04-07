// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <cmath>
#include <array>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/PolyLog.h"
#include "Ostap/Choose.h"
#include "Ostap/MoreMath.h"
// ============================================================================
// Local 
// ============================================================================
#include "local_math.h"
#include "status_codes.h"
// ============================================================================
// Polylogaritm & friends    
// ============================================================================
/*  Imaginary part of polylogarithm function \f$ Im Li_n(x) f\$
 *  @see Ostap::Math::Li
 *  @see https://en.wikipedia.org/wiki/Polylogarithm
 *  @see David Wood, "The computation of polylogarithms",
 *       Technical Report 15-92*, University of Kent, Computing Laboratory, Canterbury, UK, June 1992.
 *  @see https://www.cs.kent.ac.uk/pubs/1992/110/
 *  @see Eq (3.1) in Wood
 *  @see Ostap::Math::Li 
 *  @param n parameter
 *  @param x argument
 *  @return Imaginary part of polylogarithm function \f$ Li_n(x) \f$
 */
// ============================================================================
double Ostap::Math::ImLi
( const short  n , 
  const double x )
{
  //
  if      ( n <= 0 || x <= 1 || s_equal ( x , 1 ) ) { return  0    ; } 
  else if ( 1 == n                                ) { return -s_pi ; }
  //
  const double u = std::log ( x ) ;
  return - s_pi * std::pow ( u , n - 1 ) * igamma ( n ) ;
}
// ============================================================================
/*  Imaginary part of polylogarithm function \f$ Im Li_s(x) f\$
 *  @see Ostap::Math::Li
 *  @see https://en.wikipedia.org/wiki/Polylogarithm
 *  @see David Wood, "The computation of polylogarithms",
 *       Technical Report 15-92*, University of Kent, Computing Laboratory, Canterbury, UK, June 1992.
 *  @see https://www.cs.kent.ac.uk/pubs/1992/110/
 *  @see Eq (3.1) in Wood
 *  @see Ostap::Math::Li 
 *  @param s parameter
 *  @param x argument
 *  @return Imaginary part of polylogarithm function \f$ Li_s(x) \f$
 */
// ============================================================================
double Ostap::Math::ImLi
( const double s , 
  const double x )
{
  if      ( x <= 1 || s_equal ( x , 1 ) ) { return 0 ; } 
  else if ( isshort ( s ) )
  {
    const short n = round ( s ) ;
    return ImLi ( n , x ) ;
  }
  //
  const long double u = std::log ( x * 1.0L ) ;
  return -s_pi * std::pow ( u , s - 1.0L ) * igamma ( s ) ;
}
#include <iostream> 
// ============================================================================
namespace Ostap
{
  namespace Math
  {
    inline double  eta ( const long double x )
    { return  eta ( static_cast<double> ( x ) ) ; }
    inline double zeta ( const long double x )
    { return zeta ( static_cast<double> ( x ) ) ; }
    //
    inline std::complex<long double> zeta
    ( const std::complex<long double>&  x )
    {
      const std::complex<double> xx { x } ;
      return std::complex<long double> ( zeta ( xx ) ) ;
    }
    inline std::complex<long double> gamma 
    ( const std::complex<long double>&  x )
    {
      const std::complex<double> xx { x } ;
      return std::complex<long double> ( gamma ( xx ) ) ;
    }
  }
}
namespace
{
  // ==========================================================================
  inline bool alternating ( const long double            x    ) { return x <= 0 ; }
  template <typename TYPE>
  inline bool alternating ( const std::complex<TYPE>& /* z */ ) { return false  ; }
  // ==========================================================================
  inline bool is_not      ( const long double            x    ) { return !x     ; }
  template <typename TYPE>
  inline bool is_not      ( const std::complex<TYPE>&    z    )
  { return !z.real() && !z.imag () ; } 
  // ==========================================================================  
  template <typename TYPES, 
	    typename TYPEX>
  TYPEX Li_power  
  ( const TYPES       s  , 
    const TYPEX       x  , 
    const std::size_t NN ) 
  {
    /// mandatory condition for convergency
    Ostap::Assert ( std::abs ( x ) < 1  , 
                    "Argument is not small enough for the power series!" , 
                    "::Li_power"        , 
                    INVALID_ARGUMENT , __FILE__ , __LINE__ ) ;

    //
    static const Ostap::Math::Equal_To<TYPEX> xequal {} ;
    static const Ostap::Math::Zero<TYPEX>     xzero  {} ; 
    //
    TYPEX          result      =       1 ; // NB: no x here 
    TYPEX          term        =       1 ; // NB: no x here 
    unsigned short nSmall      =       0 ;
    //
    // alternating series? 
    const bool     xalt        = alternating ( x ) ;
    //
    for ( std::size_t k = 2 ; NN ; ++ k )
    {
      term *= x ;
      TYPEX delta = term * std::pow ( TYPEX ( 1.0L / k ) , s ) ;
      //
      if   ( is_not ( delta ) || xzero ( 2.0L * delta ) || xequal ( result + 2.0L * delta , result ) ) { ++nSmall     ; }
      else                                                                                             {   nSmall = 0 ; }
      //
      result += delta ;
      //
      // Alternating series ?
      // - one "small" term is enough,
      // - otherwise requre several  consequitive small terms are requred 
      if ( ( xalt && nSmall ) || ( 2 <= nSmall ) ) 
      {
	std::cerr << "POWER # terms " << k << " " << delta << " " << x << "# " << nSmall <<std::endl ;
	break ;
      }
    }
    std::cerr << "POWER x=" << x << std::endl ;
    return x * result ; // NB: x here
  } 
  // ==========================================================================
  long double Li_power 
  ( const int    n , 
    const double x )
  {
    const std::size_t NN  = 1000 + std::max ( 0 , -2 * n ) ;
    return Li_power ( n , 1.0L * x , NN ) ;
  }
  // ==========================================================================
  long double Li_power 
  ( const double s , 
    const double x )
  {
    const std::size_t NN  = 100000 ;
    return Li_power ( 1.0L * s , 1.0L * x , NN ) ;
  }
  // ========================================================================== 
  std::complex<double> 
  Li_power 
  ( const int                   n , 
    const std::complex<double>& z )
  {
    const std::size_t NN  = 1000 + std::max ( 0 , -2 * n ) ;
    const std::complex<long double> zz { z } ;
    const std::complex<long double> result { Li_power ( n , zz , NN ) } ;    
    return std::complex<double> ( result ) ;
  }
  // ==========================================================================  
  std::complex<double> 
  Li_power 
  ( const double                 s , 
    const std::complex<double>& z )
  {
    const std::size_t NN  = 10000 ; ;
    const std::complex<long double> zz { z } ;
    const std::complex<long double> result { Li_power ( 1.0L * s , zz , NN ) } ;    
    return std::complex<double> ( result ) ;
  }
  // ==========================================================================
  std::complex<double> 
  Li_power 
  ( const std::complex<double>& s , 
    const std::complex<double>& z )
  {
    const std::size_t NN  = 10000 ; ;
    const std::complex<long double> zz { z } ;
    const std::complex<long double> ss { s } ;    
    const std::complex<long double> result { Li_power ( s , zz , NN ) } ;    
    return std::complex<double> ( result ) ;
  }
  // ===========================================================================
  template <typename TYPES,
	    typename TYPEX> 
  inline TYPEX Li_eq_9_2_
  ( const TYPES       s , 
    const TYPEX       x , 
    const TYPEX       w ,
    const std::size_t N ) 
  {
    Ostap::Assert ( std::abs( w / s_pi ) < 1 ,
		    "Invalid range for x or w" , 
		    "::Li_eq_9_2"    , 
		    INVALID_ARGUMENT , __FILE__ , __LINE__ );
    
    //
    static const Ostap::Math::Equal_To<TYPEX> xequal {} ;
    static const Ostap::Math::Zero<TYPEX>     xzero  {} ; 
    //
    
    TYPEX          result = -1 * Ostap::Math::eta ( s ) ;
    TYPEX          term   =  1 ;
    unsigned short nSmall =  0 ;
    ///
    for ( unsigned int k = 1  ; k < N ; ++k )
    {
      term *=  ( w / k  ) ;
      const TYPES nmk   = 1 * s - k ;
      const TYPEX eta_v = Ostap::Math::eta ( nmk ) ;
      if ( is_not ( eta_v ) )  { continue ; }
      //
      const TYPEX delta = term * eta_v ;
      // 
      if ( is_not ( delta ) || xzero ( 2.0 * delta ) ||  xequal ( result - 2.0 * delta , result ) ) { ++nSmall     ; }
      else                                                                                          {   nSmall = 0 ; }
      //
      result -= delta ; 
      // 
      if ( 2 <= nSmall ) 
      {
	std::cerr << "POWER-W (9.2) #terms " << k << " " << delta << " " << x << "# " << nSmall <<std::endl ;
	break ; // NB: x here 
      }
    }
    std::cerr << "POWER-W (9.2) x=" << x << std::endl ;
    return result ;   
  }
  // ===========================================================================
  inline long double Li_eq_9_2
  ( const int         n , 
    const long double x , 
    const long double w ) 
  {    
    Ostap::Assert ( x < 0 && std::abs( w / s_pi ) < 1 ,
		    "Invalid range for x or w" , 
		    "::Li_eq_9_2"    , 
		    INVALID_ARGUMENT , __FILE__ , __LINE__ );
    return Li_eq_9_2_ ( n , 1.0L * x , 1.0L * w , 10000 ) ;
  }
  // ===========================================================================
  inline long double Li_eq_9_2
  ( const long double s , 
    const long double x , 
    const long double w ) 
  {
    Ostap::Assert ( x < 0 && std::abs( w / s_pi ) < 1 ,
		    "Invalid range for x or w" , 
		    "::Li_eq_9_2"    , 
		    INVALID_ARGUMENT , __FILE__ , __LINE__ );
    return Li_eq_9_2_ ( 1.0L * s , 1.0L * x , 1.0L * w , 1000 ) ;
  }
  // ===========================================================================
  template <typename TYPES,
            typename TYPEX>
  inline TYPEX Li_eq_9_3_
  ( const TYPES       s  , 
    const TYPEX       x  , 
    const TYPEX       w  ,
    const std::size_t NN ) 
  {
    Ostap::Assert ( std::abs( w / s_2pi ) < 1 ,
		    "Invalid range for n,x,w" , 
		    "::Li_eq_9_3"  , 
		    INVALID_ARGUMENT , __FILE__ , __LINE__ ) ;    
    //
    static const Ostap::Math::Equal_To<TYPEX> xequal {} ;
    static const Ostap::Math::Zero<TYPEX>     xzero  {} ; 
    //    
    TYPEX result = Ostap::Math::zeta ( s ) ;
    TYPEX term   = 1 ;
    unsigned short nSmall = 0 ;  
    for ( unsigned int k = 1  ; k < NN ; ++k )
    {
      term *= ( w / ( 1.0L * k )  ) ;
      const TYPES nmk  = s - TYPES ( k ) ;
      const TYPEX zeta_v = Ostap::Math::zeta ( nmk ) ;
      if ( is_not ( zeta_v ) )  { continue ; } 
      const TYPEX delta = term * zeta_v ;
      // 
      if ( is_not ( delta ) || xzero ( 2.0L * delta ) ||  xequal ( result + 2.0L * delta , result ) ) { ++nSmall     ; }
      else                                                                                            {   nSmall = 0 ; }
      //
      result += delta ; 
      // 
      if ( 3 <= nSmall ) 
      {
	std::cerr << "POWER-W (9.3) #terms " << k << " " << delta << " " << x << "# " << nSmall <<std::endl ;
	break ; // NB: x here 
      }
    }
    std::cerr << "POWER-W (9.3) x=" << x << std::endl ;
    const TYPEX tt = Ostap::Math::gamma ( 1.0L - s ) ; 
    return result + tt * std::pow ( -w , s - 1.0L ) ; 
  }
  // =========================================================================================   
  inline double Li_eq_9_3
  ( const int         n , 
    const long double x , 
    const long double w ) 
  {
    Ostap::Assert ( ( n <= 0 ) && ( x > 0 ) && std::abs( w / s_2pi ) < 1 ,
		    "Invalid range for n,x,w" , 
		    "::Li_eq_9_3"  , 
		    INVALID_ARGUMENT , __FILE__ , __LINE__ );
    //
    return Li_eq_9_2_ ( n , 1.0L * x , 1.0L * w , 10000 ) ;
  } ; 
  // =========================================================================================
  inline double Li_eq_9_3
  ( const long double s , 
    const long double x , 
    const long double w ) 
  {
    Ostap::Assert ( ( x > 0 ) && std::abs( w / s_2pi ) < 1 ,
		    "Invalid range for n,x,w" , 
		    "::Li_eq_9_3"  , 
		    INVALID_ARGUMENT , __FILE__ , __LINE__ );
    //
    return Li_eq_9_3_ ( 1.0L * s  , 1.0L * x , 1.0L * w , 10000 ) ;
  }
  // =========================================================================
  inline std::complex<double> Li_eq_9_3
  ( const int                  n , 
    const std::complex<double> x , 
    const std::complex<double> w ) 
  {
    Ostap::Assert ( ( n <= 0 ) && std::abs( w ) / s_2pi < 1 ,
                  "Invalid range for n,x,w" , 
                  "::Li_eq_9_3"  , 
                  INVALID_ARGUMENT , __FILE__ , __LINE__ );
    
    const std::complex<long double> xx { x } ;
    const std::complex<long double> ww { w } ;
    const std::complex<long double> rr { Li_eq_9_3_ ( n , xx , ww , 10000 ) } ;
    return std::complex<double> ( rr ) ;
  }    
  // =========================================================================
  inline std::complex<double> Li_eq_9_3
  ( const long double          s , 
    const std::complex<double> x , 
    const std::complex<double> w ) 
  {
    Ostap::Assert ( std::abs( w ) / s_2pi < 1 ,
		    "Invalid range for x,w" , 
		    "::Li_eq_9_3"  , 
		    INVALID_ARGUMENT , __FILE__ , __LINE__ );
    
    const std::complex<long double> xx { x } ;
    const std::complex<long double> ww { w } ;
    const std::complex<long double> rr { Li_eq_9_3_ ( 1.0L * s , xx , ww  , 10000 ) } ;
    return std::complex<double> ( rr ) ;
  }      
  // =========================================================================
  inline std::complex<double> Li_eq_9_3
  ( const std::complex<double> s , 
    const std::complex<double> x , 
    const std::complex<double> w ) 
  {
    Ostap::Assert ( std::abs( w ) / s_2pi < 1 ,
                  "Invalid range for x,w" , 
                  "::Li_eq_9_3"  , 
                  INVALID_ARGUMENT , __FILE__ , __LINE__ );
    
    const std::complex<long double> xx { x } ;
    const std::complex<long double> ww { w } ;
    const std::complex<long double> ss { s } ;
    const std::complex<long double> rr { Li_eq_9_3_ ( ss , xx , ww , 10000 ) } ;
    return std::complex<double> ( rr ) ;
  }    
  // ==========================================================================
  template <typename TYPEX>  
  inline TYPEX Li_eq_9_5_ 
  ( const int         n , 
    const TYPEX       x , 
    const TYPEX       w ,
    const std::size_t N ) 
  {
    Ostap::Assert ( ( 1 <= n ) && std::abs( w / s_2pi ) < 1 ,
		    "Invalid range for n,x,w"  , 
		    "::Li_eq_9_5"              , 
		    INVALID_ARGUMENT           , __FILE__ , __LINE__ ) ;

    static const Ostap::Math::Equal_To<TYPEX> xequal {} ;
    static const Ostap::Math::Zero<TYPEX>     xzero  {} ; 
    //
    TYPEX result  = Ostap::Math::zeta ( n ) ;
    TYPEX term    = 1 ;
    unsigned int nSmall = 0 ;
    for ( unsigned int k = 1 ; k < N ; ++k ) 
    {
      term *= ( w / ( 1.0L * k ) ) ; 
      if ( 1 * n  == k + 1 ) { continue ; } // skip singularity 
      const int   nmk    = 1 * n - k    ;
      const TYPEX zeta_v = Ostap::Math::zeta ( nmk ) ;
      std::cerr << "EQ9.5 " << k << " " << nmk  << " " << zeta_v << std::endl ;
      if ( is_not ( zeta_v ) )  { continue ; }
      //
      const TYPEX delta = term * zeta_v ;
      //
      if ( is_not ( delta ) || xzero ( 2.0L * delta ) ||  xequal ( result + 2.0L * delta , result ) ) { ++nSmall     ; }
      else                                                                                            {   nSmall = 0 ; }
      //
      result += delta ; 
      // 
      if ( 2 <= nSmall ) 
      {
	std::cerr << "# EQ (9.5) terms " << k << " " << delta << " " << x << "# " << nSmall <<std::endl ;
	break ;  
      }
    }
    //
    // const unsigned short nm1 = n - 1 ; 
    // result += ( Ostap::Math::harmonic ( nm1 ) - std::log ( std::abs ( w ) ) ) * std::pow ( w , nm1 ) * Ostap::Math::igamma ( nm1 ) ;
    return result ;   
  }
  // ==========================================================================
  inline double Li_eq_9_5
  ( const int         n , 
    const long double x ,
    const long double w )
  {
    Ostap::Assert ( ( 1 <= n ) && std::abs( w ) / s_2pi < 1 ,
		    "Invalid range for n,w"    , 
		    "::Li_eq_9_5"              , 
		    INVALID_ARGUMENT           , __FILE__ , __LINE__ ) ;
    //
    auto sum = Li_eq_9_5_ ( n , 1.0L * x , 1.0L * w , 10000 ) ;
    //
    const unsigned short nm1 = n - 1 ;
    sum += ( Ostap::Math::harmonic ( nm1 ) - std::log ( std::abs ( w ) ) ) * std::pow ( w , nm1 ) * Ostap::Math::igamma ( nm1 ) ;
    return sum ;
  }
  // ==========================================================================
  inline std::complex<double> Li_eq_9_5
  ( const int                   n , 
    const std::complex<double>& x ,
    const std::complex<double>& w )
  {
    Ostap::Assert ( ( 1 <= n ) && std::abs( w ) / s_2pi < 1 ,
		    "Invalid range for n,w"    , 
		    "::Li_eq_9_5"              , 
		    INVALID_ARGUMENT           , __FILE__ , __LINE__ ) ;

    const std::complex<long double> xx  { x } ;
    const std::complex<long double> ww  { w } ;
    //
    auto sum { Li_eq_9_5_ ( n , xx , ww , 10000 ) } ;
    //
    const unsigned short nm1 = n - 1 ;
    const long double    hm  =  Ostap::Math::harmonic ( nm1 )  ;
    const long double    ig  =  Ostap::Math::igamma   ( nm1 )  ;    
    sum += ( hm - std::log ( -ww ) ) * std::pow ( ww , nm1  ) * ig ;
    //
    return std::complex<double> ( sum ) ;
  }
  // ==========================================================================
}
// ============================================================================
/*  Polylogarithm function  \f$ Li_n(x)  = \sum \frac{x^k}{k^n} f\$
 *  @see https://en.wikipedia.org/wiki/Polylogarithm
 *  @see https://www.cs.kent.ac.uk/pubs/1992/110/
 *  @param n parameter
 *  @param x argument
 *  @return the value of polylogarithm function \f$ \f$
 *  @attention for \$ 1 < x \f$ the real part is returned, for imaginary part see      
 *  @see  Ostap::Math::ImLi
 */
// ============================================================================
double Ostap::Math::Li
( const short  n ,
  const double x )
{
  if       ( !x || s_zero ( x ) ) { return 0          ; }
  else if  ( s_equal ( x ,  1 ) ) { return zeta ( n ) ; } 
  else if  ( s_equal ( x , -1 ) ) { return -eta ( n ) ; }
  //
  
  // explicit case:
  if      ( 1  == n ) { return - std::log ( 1.0L - std::abs ( x ) )         ; } // Eq (6.1)
   /// dilogarithm
  else if ( 2  == n )
  {

    static const short two = 2 ;
    
    // helper constant 
    static const long double s_c { s_pi2 / 6.0L } ;
    
    /// very special case
    static const long double s_Li2_2 = s_pi2 / 4 ;  
    if       ( s_equal (  2 , x ) ) { return                          s_Li2_2 ; }  // RETURN 
    else if  ( s_equal ( -2 , x ) ) { return Li ( two , 2 * 2 ) / 2 - s_Li2_2 ; }  // RETURN  

    if       ( -1 > x )
    {
      /// Eq. (7.3)
      const long double l1 = std::log ( - 1.0L * x ) ;
      const long double l2 = l1 * l1 ; 
      return -s_c - 0.5L * l2 - Li ( two , 1.0 / x ) ;
    }
    else if ( -0.5 > x )
    {
      /// Eq. (7.3)
      const long double l2  = std::log ( 1.0L - x ) * std::log ( ( 1.0L - x ) / ( x * x ) ) ;
      return -s_c +0.5L * l2  +  Li ( two , 1 / ( 1 - x ) ) ;
    }
    else if ( 1 < x )
    {
      /// Eq. (7.3)
      const long double l1 = std::log ( 1.0L * x ) ;
      const long double l2 = l1 * l1 - s_pi2       ;        // ATTENTION: pi^2 here fro m complex logarithm! 
      return -s_c - 0.5L * l2 - Li ( two , 1.0 / x ) ;
    }
    else if ( 0.5 < x )
    {
      /// Eq. (7.3)
      const long double LL = std::log ( 1.0L * x ) * std::log ( 1.0L - x ) ;
      return s_c - LL  - Li ( two , 1 - x ) ;
    }

    /// couple of very special cases 
    static const long double s_Li2_half = s_pi2/12.0L - 0.5L * s_ln2 * s_ln2 ;
    if      ( s_equal (  0.5 , x ) ) { return                              s_Li2_half ; } // RETURN 
    else if ( s_equal ( -0.5 , x ) ) { return Li ( two , 0.5 * 0.5 ) / 2 - s_Li2_half ; } // RETURN
    //

    /// ATTENTION! 
    if ( 0 > x  )
    {
      // Eq (14.1) 
      return Li ( two , x * x ) / 2 - Li ( two , -x ) ; 
    }
    
    std::cerr << "dilog case " << x << std::endl ;
  }
  //

  std::cerr << "Li(" << n << "," << x << ") / 1" << std::endl ;
  
  // RATIONAL CASES 
  //  
  if      (  0 == n ) { return x                /          ( 1.0L - x )     ; } // Eq.(6.2)
  else if ( -1 == n ) { return x                / std::pow ( 1.0L - x , 2 ) ; } // Eq.(6.2)
  else if ( -2 == n ) { return x * ( x + 1.0L ) / std::pow ( 1.0L - x , 3 ) ; } // Eq.(6.2)
  // Rational case 
  else if  ( -1 * N_EULERIAN_MAX <= n && n < 0 ) 
  {
    const unsigned short NN = -n ;
    Eulerian eu { NN } ;
    std::cerr << "Eulerian explicit " << x << std::endl ;
    return x * eu ( x ) / std::pow ( 1.0L - x , NN + 1 ) ;
  }
  //

  std::cerr << "Li(" << n << "," << x << ") / 2" << std::endl ;

  const long double w = std::log ( 1.0L * std::abs ( x ) ) ;

  // for small z use power series 
  
  const long double absw1pi = std::abs ( w * s_1_pi ) ; 
  const long double absw2pi = 0.5L* absw1pi           ;  
  const long double absx    = std::abs ( x          ) ;

  //  use simple power series 
  const bool xsmall_1 = ( 0 > x ) && ( absx <= absw1pi ) ; // for 
  const bool xsmall_2 = ( 0 < x ) && ( absx <= absw2pi ) ; // for Eq (9.3) & Eq.(9.5) 

  std::cerr << "Li(" << n << "," << x << ") / 3"
	    << " small1:" << xsmall_1 
	    << " small2:" << xsmall_2 
	    << std::endl ;
  
  /// (A) simple power serie  
  if ( std::abs ( x ) <= 0.25 || xsmall_1 || xsmall_2 ) { return Li_power  ( n , x ) ; }

  std::cerr << "Li(" << n << "," << x << ") / 4"
	    << " abs(w/pi):" << absw1pi  << std::endl ;

  /// (B) use Eq (9.2)
  if ( ( x < 0 )  && ( absw1pi <= 0.5 ) )
  { return Li_eq_9_2 ( n , 1.0L * x , 1.0L * w ) ; }
  
  std::cerr << "Li(" << n << "," << x << ") / 5"
	    << " abs(w/2pi):" << absw2pi 
	    << std::endl ;

  /// (C) use Eq (9.3)
  if ( ( n < 0 ) && ( 0 < x ) && ( absw2pi <= 0.512 ) )
  { return Li_eq_9_3 ( n , 1.0L * x , 1.0L * w ) ; }
  
  std::cerr << "Li(" << n << "," << x << ") / 6"
	    << " abs(w/2pi):" << absw2pi
	    << std::endl ;

  /// (D) use Eq (9.5)
  if ( ( n > 0 ) && ( 0 < x ) && ( absw2pi <= 0.512 ) )
  { return Li_eq_9_5 ( n , 1.0L * x , 1.0L * w ) ; }

  std::cerr << "Li(" << n << "," << x << ") / 7" << std::endl ;

  /// (E) Eq. (10.3)
  if ( 0 > n && 1 < std::abs ( x ) )
  {
    std::cerr << "Reciprocal 10.3 " << x << std::endl ;
    return ( 0 == n % 2 ? -1 : +1 ) * Li ( n , 1.0 / x ) ;
  }
  
  std::cerr << "Li(" << n << "," << x << ") / 8" << std::endl ;

  /// (F) Eq. (10.1)
  if ( 0 < n && 1 < x ) 
  {
    long double          rr = zeta ( 0 ) ;
    long double          tt = 1 ; 
    const long double    w2 = w * w ;
    for ( unsigned short k  = 1 ; 2 * k <= n ; ++k )
    {
      tt *= ( n + 2 - 2 * k ) * ( n + 1 - 2 * k ) / w2 ;
      rr += tt * zeta ( 2 * k ) ;      
    }
    const long double R = std::pow ( w , n ) * igamma ( n + 1 ) * rr ; 
    std::cerr << "Reciprocal 10.1 " << x << std::endl ;
    return R + ( 0 == n % 2 ? -1 : +1 ) * Li ( n , 1.0 / x ) ;
  }    

  std::cerr << "Li(" << n << "," << x << ") / 9" << std::endl ;
  
  /// (G) negative argument? (Eq 14.1)
  if ( x < 0  )
  {
    // Use the square formula Eq.(14.1) to convert to the positive argument
    std::cerr << "Square14.1 " << x << std::endl ;
    return std::pow ( 2 , 1 - n ) * Li ( n , x * x ) - Li ( n , -x ) ;
  }
  
  std::cerr << "Li(" << n << "," << x << ") / 10" << std::endl ;
  
  Ostap::Assert ( false ,
		  "Not yet implemented!" ,
		  "Ostap::Math::Li"      ,
		  NOT_IMPLEMENTED_YET    , __FILE__ , __LINE__ ) ;
  
  return std::numeric_limits<double>::quiet_NaN() ; 
}
// ============================================================================
/*  Polylogarithm function  \f$ Li_s(x)  = \sum \frac{x^k}{k^s} f\$
 *  @see https://en.wikipedia.org/wiki/Polylogarithm
 *  @see https://www.cs.kent.ac.uk/pubs/1992/110/
 *  @param s parameter
 *  @param x argument
 *  @return the value of polylogarithm function \f$ \f$
 *  @attention for \$ 1 < x \f$ the real part is returned, for imaginary part see      
 *  @see  Ostap::Math::ImLi
 */
// ============================================================================
double Ostap::Math::Li
( const double s ,
  const double x )
{
  if      ( !x || s_zero ( x ) ) { return 0  ; }
  else if ( isshort ( s ) )
  {
    const short n = round ( s ) ;
    return Li ( n , x ) ;
  }
  else if ( s_equal ( x ,  1 ) ) { return zeta ( s ) ; } 
  else if ( s_equal ( x , -1 ) ) { return -eta ( s ) ; }  
  //
  
  //
  // Here s is *not*  (short) integer
  //

  const long double w = std::log ( 1.0L * std::abs ( x ) ) ;
  
  ///
  const long double absw1pi = std::abs ( w * s_1_pi ) ; 
  const long double absw2pi = 0.5L* absw1pi           ;  
  const long double absx    = std::abs ( x          ) ;
  
  ///  use simple power series 
  const bool xsmall_1 = ( 0 > x ) && ( absx <= absw1pi ) ; // for 
  const bool xsmall_2 = ( 0 < x ) && ( absx <= absw2pi ) ; // for Eq (9.3) & Eq.(9.5) 

  /// (A) simple power serie  for small x 
  if ( std::abs ( x ) <= 0.25 || xsmall_1 || xsmall_2 ) { return Li_power  ( s , x ) ; }
  
  /// (B) use Eq (9.2)
  if ( ( x < 0 )  && ( absw1pi <= 0.5 ) )
  { return Li_eq_9_2 ( 1.0L * s , 1.0L * x , 1.0L * w ) ; }
  
  /// (C) use Eq (9.3)
  if ( ( 0 < x ) && ( absw2pi <= 0.512 ) )
  { return Li_eq_9_3 ( 1.0L * s , 1.0L * x , 1.0L * w ) ; }
  
  /// use the square formula to reduce the modulus 
  if ( 0 < x )
  {
    const double sqx = std::sqrt ( x ) ; 
    return ( Li ( s , +sqx ) + Li ( s , -sqx ) ) * std::pow ( 2.0 , s - 1 ) ;
  }
  else
  {
    /// attention: here we jump into complex variant 
    const double sqx = std::sqrt ( std::abs ( x ) ) ;
    const std::complex<double> l1 { Li ( s , +s_j * sqx ) } ;
    const std::complex<double> l2 { Li ( s , -s_j * sqx ) } ;
    return std::pow ( 2.0L , s - 1 ) * ( l1.real () + l2.real () ) ;
  }  
  //
  Ostap::Assert ( false ,
		  "Not yet implemented!" ,
		  "Ostap::Math::Li"      ,
		  NOT_IMPLEMENTED_YET , __FILE__ , __LINE__ ) ;
  //
  return std::numeric_limits<double>::quiet_NaN () ;   
}
// ===========================================================================
/** Polylogarithm function  \f$ Li_n(x)  = \sum \frac{x^k}{k^n} f\$
 *  @see https://en.wikipedia.org/wiki/Polylogarithm
 *  @see David Wood, "The computation of polylogarithms",
 *       Technical Report 15-92*, University of Kent, Computing Laboratory, Canterbury, UK, June 1992.
 *  @see https://www.cs.kent.ac.uk/pubs/1992/110/
 *  @param n parameter
 *  @param x argument
 *  @return value of polylogarithm function \f$ Li_n(x) \f$
 *  @see  Ostap::Math::ImLi
 *
 */
// ===========================================================================
std::complex<double>
Ostap::Math::Li
( const short                 n ,
  const std::complex<double>& z )
{
  const double x = z.real () ;
  const double y = z.imag () ;
  //
  if ( !y || s_zero ( y ) )
  {
    return 0 <= y ?
      std::complex<double> ( Li ( n , x ) , +ImLi ( n , x ) ) : 
      std::complex<double> ( Li ( n , x ) , -ImLi ( n , x ) ) ;
  }
  //
  const std::complex<long double> zz { z } ;
  
  typedef std::complex<double> RR ;
  // explicit case:
  if      (  1 == n ) { return RR ( - std::log ( 1.0L - zz ) ) ; }                       // Eq (6.1)
  else if (  0 == n ) { return RR ( zz                /          ( 1.0L - zz     ) ) ; } // Eq.(6.2)
  else if ( -1 == n ) { return RR ( zz                / std::pow ( 1.0L - zz , 2 ) ) ; } // Eq.(6.2)
  else if ( -2 == n ) { return RR ( zz * ( x + 1.0L ) / std::pow ( 1.0L - zz , 3 ) ) ; } // Eq.(6.2)
  // Rational case 
  else if  ( -1 * N_EULERIAN_MAX <= n && n < 0 ) 
  {
    const unsigned short NN = -n ;
    Eulerian eu { NN } ;
    std::cerr << "Eulerian explicit " << x << std::endl ;
    return eu ( z ) * RR ( zz / std::pow ( 1.0L - zz , NN + 1 ) ) ;
  }
  
  const std::complex<double> w { std::log ( z ) } ;
  
  ///
  const long double absw1pi = std::abs ( w ) * s_1_pi ; 
  const long double absw2pi = 0.5L* absw1pi           ;  
  const long double absz    = std::abs ( z ) ;
  
  ///  use simple power series 
  const bool xsmall_2 = absz <= absw2pi ;
  
  /// (A) simple power serie  for small x 
  if ( absz <= 0.25 || xsmall_2 ) { return Li_power  ( n , z ) ; }

  /// (B) Use Eq (9.3) 
  if ( ( n < 0 ) && ( absw2pi <= 0.512 ) )
  { return Li_eq_9_3 ( n , z , w  ) ; }

  /// (C) Use Eq (9.5) 
  if ( ( n > 0 ) && ( absw2pi <= 0.512 ) )
  { return Li_eq_9_5 ( n , z , w  ) ; }

  /// (D) Use the square formula top reduce modulus 
  if ( 1 < absz )
  {
    const std::complex<double> sqz { std::sqrt ( z ) } ;
    return ( Li ( n , +sqz ) + Li ( n  , -sqz ) ) * std::pow ( 2.0 , n - 1 ) ;
  }
  
  Ostap::Assert ( false ,
		  "Not yet implemented!" ,
		  "Ostap::Math::Li"      ,
		  NOT_IMPLEMENTED_YET , __FILE__ , __LINE__ ) ;
  //
  return std::numeric_limits<double>::quiet_NaN () ;     
}
// ===========================================================================
/** Polylogarithm function  \f$ Li_s(x)  = \sum \frac{x^k}{k^s} f\$
 *  @see https://en.wikipedia.org/wiki/Polylogarithm
 *  @see David Wood, "The computation of polylogarithms",
 *       Technical Report 15-92*, University of Kent, Computing Laboratory, Canterbury, UK, June 1992.
 *  @see https://www.cs.kent.ac.uk/pubs/1992/110/
 *  @param s parameter
 *  @param x argument
 *  @return value of polylogarithm function \f$ Li_s(x) \f$
 *  @see  Ostap::Math::ImLi
 *
 */
// ===========================================================================
std::complex<double>
Ostap::Math::Li
( const double                s ,
  const std::complex<double>& z )
{
  //
  if ( isshort ( s ) )
  {
    const short n = round ( s ) ;
    return Li ( n , z ) ;
  }
  //
  const double x = z.real () ;
  const double y = z.imag () ;
  //
  if ( !y || s_zero ( y ) )
  {
    return 0 <= y ?
      std::complex<double> ( Li ( s , x ) , +ImLi ( s , x ) ) :
      std::complex<double> ( Li ( s , x ) , -ImLi ( s , x ) ) ;
  }
  //
  // Here s is *NOT* (short) integer!  
  //
  const std::complex<double> w { std::log ( z ) } ;
  
  const long double absw1pi = std::abs ( w ) * s_1_pi ; 
  const long double absw2pi = 0.5L* absw1pi           ;  
  const long double absz    = std::abs ( z ) ;
  
  ///  use simple power series 
  const bool xsmall_2 = absz <= absw2pi ;
  
  /// (A) simple power serie  for small x 
  if ( absz <= 0.25 || xsmall_2 ) { return Li_power ( s , z ) ; }
  
  
  /// (B) Use Eq (9.3) 
  if (  absw2pi <= 0.512 ) { return Li_eq_9_3 ( 1.0L * s , z , w  ) ; }

  /// (D) Use the square formula top reduce modulus 
  if ( 1 < absz )
  {
    const std::complex<double> sqz { std::sqrt ( z ) } ;
    return ( Li ( s , +sqz ) + Li ( s  , -sqz ) ) * std::pow ( 2.0 , s - 1 )  ;
  }

  //
  Ostap::Assert ( false ,
		  "Not yet implemented!" ,
		  "Ostap::Math::Li"      ,
		  NOT_IMPLEMENTED_YET , __FILE__ , __LINE__ ) ;
  //
  return std::numeric_limits<double>::quiet_NaN () ;     
}
// ===========================================================================
/** Polylogarithm function  \f$ Li_s(x)  = \sum \frac{x^k}{k^s} f\$
 *  @see https://en.wikipedia.org/wiki/Polylogarithm
 *  @see David Wood, "The computation of polylogarithms",
 *       Technical Report 15-92*, University of Kent, Computing Laboratory, Canterbury, UK, June 1992.
 *  @see https://www.cs.kent.ac.uk/pubs/1992/110/
 *  @param s parameter
 *  @param x argument
 *  @return value of polylogarithm function \f$ Li_s(x) \f$
 *  @see  Ostap::Math::ImLi
 *
 */
// ===========================================================================
std::complex<double>
Ostap::Math::Li
( const std::complex<double>& s , 
  const std::complex<double>& z )
{
  //
  const double sx = s.real () ;
  const double sy = s.imag () ;
  //
  if ( !sy || s_zero ( sy ) ) { return Li ( sx , z ) ; } 
  //
  
  std::cerr << "(0) Li(C,C)<1" << z << std::endl ;  
  
  //
  // Here s is *NOT* (short) integer!  
  //
  const std::complex<double> w { std::log ( z ) } ;
  
  const long double absw1pi = std::abs ( w ) * s_1_pi ; 
  const long double absw2pi = 0.5L* absw1pi           ;  
  const long double absz    = std::abs ( z ) ;
  
  ///  use simple power series 
  const bool xsmall_2 = absz <= absw2pi ;
  
  /// (A) simple power serie  for small x 
  if ( absz <= 0.25 || xsmall_2 ) { return Li_power  ( s , z ) ; }
  
  /// (B) Use Eq (9.3) 
  if (  absw2pi <= 0.512 ) { return Li_eq_9_3 ( s , z , w  ) ; }
  
  /// (C) Use the square formula to reduce modulus 
  if ( 1 < absz )
  {
    const std::complex<double> sqz { std::sqrt ( z ) } ;
    std::cerr << "(1) Li(C,C)<1" << z << std::endl ;  
    return ( Li ( s , +sqz ) + Li ( s  , -sqz ) ) * std::pow ( 2.0 , s - 1.0 ) ;    
  }
  
  std::cerr << "(2) Li(C,C)<1" << z << std::endl ;  
  
  Ostap::Assert ( false ,
		  "Not yet implemented!" ,
		  "Ostap::Math::Li"      ,
		  NOT_IMPLEMENTED_YET , __FILE__ , __LINE__ ) ;
  //
  return std::numeric_limits<double>::quiet_NaN () ;     
}
// ===========================================================================
//                                                                     The END 
// ===========================================================================
