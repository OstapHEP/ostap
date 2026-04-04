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
/* Imaginary part of polylogarithm function \f$ Im Li_n(x) f\$
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
/* Imaginary part of polylogarithm function \f$ Im Li_s(x) f\$
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
/* polylogarithm function  \f$ Li_n(x)  = \sum \frac{x^k}{k^n} f\$
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
  
  ///
  if      ( 0  == n ) { return x /        ( 1.0L - x )              ; } // Eq.(6.2)
  else if ( 1  == n ) { return - std::log ( 1.0L - std::abs ( x ) ) ; } // Eq (6.1)
  /// dilogarithm
  else if ( 2  == n )
  {
    static const short two = 2u ;
    
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
  
  if       ( -1 == n  ) { return x                / std::pow ( 1.0L - x , 2 ) ; } // Eq.(6.2)
  else if  ( -2 == n  ) { return x * ( x + 1.0L ) / std::pow ( 1.0L - x , 3 ) ; } // Eq.(6.2)
  // Rational case 
  else if  ( N_EULERIAN_MAX <= -n && n < 0 )
  {
    const unsigned short NN = -n ;
    Eulerian eu { NN } ;
    return x * eu ( x ) / std::pow ( 1.0L - x , NN + 1 ) ;
  }
  //

  //
  const long double w = std::log ( 1.0L * std::abs ( x ) ) 

  // use Eq. (9.2)
  if ( 0 > x && std::abs ( w ) < s_pi_2 )
  {
    const long double w1p    = w * s_1_pi ; 
    long double       result = eta ( n )  ;
    long double       term   = 1 ;
    unsigned short nSmall = 0 ;  
    for ( insigned int k = 1  ; k < 100 ; ++k )
    {
      term /=  ( w1p / k  ) ;
      const long double delta = term * eta ( n - k ) ; 
      if ( !delta || s_zero ( delta ) ||  s_equal ( result + delta , result ) ) { ++nSmall     ; }
      else                                                                    ) {   nSmall = 0 ; }
      //
      result += delta ; 
      // 
      if ( 4 <= nSmall ) 
      {
	      std::cerr << "# terms " << k << " " << delta << " " << x << "# " << nSmall <<std::endl ;
	      return result ; // NB: x here 
      }
    return result ;   
    }
  }
  if ( 0 < x && n < 0 && atd::abs ( w / s_2pi ) < 0.512 i )
  {
  
  }

  // Simple serie for small argument Eq. (8.1)
  // 
  if ( 2 <= n && std::abs ( s_2pi * x ) <=  w ) 
    {
    long double    result      = 1 ; // NB: no x here 
    long double    term        = 1 ; // NB: no x here 
    unsigned short nSmall      = 0 ;
    const bool     alternating = x < 0 ;
    for ( unsigned short k = 2 ; k < 50 ; ++ k )
    {
      term *= x ;
      const long double delta = term * std::pow ( 1.0L / k , n ) ;
      //
      if ( !delta || s_zero ( delta ) || s_equal ( result + delta , result ) ) { ++nSmall     ; }
      else                                                                     {   nSmall = 0 ; }
      //
      result += delta ;
      //
      // Alternating series ?
      // - one "small" term is enough,
      // - otherwise requre several  consequitive small terms are requred 
      if ( ( 2 <= nSmall && alternating ) || ( 4 <= nSmall ) ) 
      {
	std::cerr << "# terms " << k << " " << delta << " " << x << "# " << nSmall <<std::endl ;
	return x * result ; // NB: x here 
      }
    }

    std::cerr << "slow convergency" << x << std::endl ;
    return x * result ; // NB: x here
  }
  // 
  
  //
  // Eq (10.3) 
  if ( 0 > n && 1 < std::abs ( x ) )
  { return ( ( n + 1 ) % 2 ? 1 : -1 ) * Li ( n , 1.0 /x ) ; }
  //

  // Eq (9.5) 
  if ( 0 < n && 0 < x ) 
  {
    const double w = std::log ( x ) ;
    /// actual convergency for abs(w)<2*pi
    /// Eq (9.7) 
    if ( std::abs ( w ) < 1.5 * s_pi )
    {
      long double result = zeta ( n ) ;
      long double term   = 1 ;
      for ( unsigned short k = 1 ; k + 2 <= n  ; ++k )
      {
	term *= ( w / k ) ;
	result += term ;
      }
      //
      term   *= w / ( n - 1 ) ;
      const unsigned short n1 = n -1 ;
      result +=  ( harmonic ( n1 ) - std::log ( std::abs ( w ) ) ) * term ;
      term   *= w /   n       ;      
      result -= 0.5 * term ;  
      //
      const long double   wpi = w / s_2pi ;
      long double result2 = zeta ( 2 ) * igamma ( n ) ;
      long double term2   = 1 ;      
      for ( unsigned short j   = 1 ; j < 1000 ; ++j )
      {
	/// NOTE YEt
      }      
    }
  }
  
  
  
  
  Ostap::Assert ( false ,
		  "Not yet implemented!" ,
		  "Ostap::Math::Li"      ,
		  NOT_IMPLEMENTED_YET , __FILE__ , __LINE__ ) ;
  
  // not yet ;

  if ( 1 == n ) { return std::numeric_limits<double>::quiet_NaN() ; }

  return std::numeric_limits<double>::quiet_NaN() ; 
  
}
// ============================================================================
/* polylogarithm function  \f$ Li_s(x)  = \sum \frac{x^k}{k^s} f\$
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
  if       ( !x || s_zero ( x ) ) { return 0             ; }
  //  
  Ostap::Assert ( false ,
		  "Not yet implemented!" ,
		  "Ostap::Math::Li"      ,
		  NOT_IMPLEMENTED_YET , __FILE__ , __LINE__ ) ;
  
  return std::numeric_limits<double>::quiet_NaN() ;   
}

// ===========================================================================
//                                                                     The END 
// ===========================================================================