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
  // RATIONAL CASES 
  //  
  if      (  0  == n ) { return x                /          ( 1.0L - x )     ; } // Eq.(6.2)
  else if ( -1 == n  ) { return x                / std::pow ( 1.0L - x , 2 ) ; } // Eq.(6.2)
  else if ( -2 == n  ) { return x * ( x + 1.0L ) / std::pow ( 1.0L - x , 3 ) ; } // Eq.(6.2)
  // Rational case 
  else if  (  0 >  n && N_EULERIAN_MAX >= std::abs ( n ) ) 
  {
    const unsigned short NN = -n ;
    Eulerian eu { NN } ;
    return x * eu ( x ) / std::pow ( 1.0L - x , NN + 1 ) ;
  }
  //

  const long double w = std::log ( 1.0L * std::abs ( x ) ) ;

  // for small z use power series 
  
  const long double absw1pi = std::abs ( w * s_1_pi ) ; 
  const long double absw2pi = 0.5L* absw1pi           ;  
  const long double absx    = std::abs ( x          ) ;

  //  use simple power series 
  const bool xsmall_1 = ( 0 > x ) && ( absx <= absw1pi ) ; // for 
  const bool xsmall_2 = ( 0 < x ) && ( absx <= absw2pi ) ; // for Eq (9.3) & Eq.(9.5) 

  /// (A) simple power serie  
  if ( std::abs ( x ) < 0.25 || xsmall_1 || xsmall_2 )
  {
    long double    result      =       1 ; // NB: no x here 
    long double    term        =       1 ; // NB: no x here 
    unsigned short nSmall      =       0 ;
    const bool     alternating = ( x < 0 ) ;
    for ( unsigned short k = 2 ; k < std::abs ( n ) + 50 ; ++ k )
    {
      term *= x ;
      const long double kterm = 0 <= n ? std::pow ( 1.0L/k , n ) : std::pow ( 1.0L * k , std::abs ( n ) ) ;
      const long double delta = term * kterm ;
      //
      if ( !delta || s_zero ( 2 * delta ) || s_equal ( result + 2 * delta , result ) ) { ++nSmall     ; }
      else                                                                             {   nSmall = 0 ; }
      //
      result += delta ;
      //
      // Alternating series ?
      // - one "small" term is enough,
      // - otherwise requre several  consequitive small terms are requred 
      if ( ( alternating && nSmall ) || ( 2 <= nSmall ) ) 
      {
	      std::cerr << "# terms " << k << " " << delta << " " << x << "# " << nSmall <<std::endl ;
	      return x * result ; // NB: x here 
      }
    }
    std::cerr << "slow convergency" << x << std::endl ;
    return x * result ; // NB: x here
  }
  
  //  (B)  Use Eq 9.2 
  if  ( ( 0 > x ) && ( absw1pi <= 0.5 ) ) 
  {
    long double       result = eta ( n ) ;
    long double       term   = 1 ;
    unsigned short nSmall = 0 ;  
    for ( unsigned int k = 1  ; k < std::abs ( n ) + 100 ; ++k )
    {
      term /=  ( w / k  ) ;
      const int         nmk   = 1 * n - k ;
      const long double eta_v = eta ( nmk ) ;
      if ( !eta_v )  { continue ; } 
      const long double delta = term * eta_v ;
      // 
      if ( !delta || s_zero ( 2 * delta ) ||  s_equal ( result + 2 * delta , result ) ) { ++nSmall     ; }
      else                                                                              {   nSmall = 0 ; }
      //
      result += delta ; 
      // 
      if ( 2 <= nSmall ) 
      {
	      std::cerr << "# terms " << k << " " << delta << " " << x << "# " << nSmall <<std::endl ;
	      return result ; // NB: x here 
      }
    return result ;   
    }
  }
  // 
  //  (c) use Eq ( 9.5 )  
  if ( 0 < x && 2 <= n && absw2pi < 0.512 )
  {
    long double result  = zeta ( n ) ;
    long double term    = 1 ;
    unsigned int nSmall = 0 ;
    for ( unsigned int k = 1 ; k < n + 100 ; ++k ) 
    {
      term /= ( w / k ) ; 
      if ( 1 * n  == k + 1 ) { continue ; } // skip singularity 
      const int         nmk    = 1 * n - k    ;
      const long double zeta_v = zeta ( nmk ) ;
      if ( !zeta_v )  { continue ; }
      const long double delta = term * zeta_v ;
      //
      if ( !delta || s_zero ( 2 * delta ) ||  s_equal ( result + 2 * delta , result ) ) { ++nSmall     ; }
      else                                                                              {   nSmall = 0 ; }
      //
      result += delta ; 
      // 
      if ( 2 <= nSmall ) 
      {
	      std::cerr << "# terms " << k << " " << delta << " " << x << "# " << nSmall <<std::endl ;
	      break ;  
      }
    }
    //
    const unsigned short nm1 = n - 1 ; 
    result += ( harmonic ( nm1 ) - std::log ( std::abs ( w ) ) ) *  std::pow ( w , n - 1 ) * igamma ( nm1 ) ;
    return result ;   
  }


  /// negative argument?
  if ( x <= -0.25 )
  {
    // Use square formula to convert to positive arguments
    return std::pow ( 2 , 1 - n ) * Li ( n , x * x ) - Li ( n , -x ) ;
  }


  // Eq. (10.3)
  if      ( 0 > n && 1 < std::abs ( x ) )
  { return ( ( n + 1 ) % 2 ? 1 : -1 ) * Li ( n , 1.0 / x ) ; }
  // Eq. (10.1)
  else if ( 0 < n && 1 <  x  ) 
  {
    long double       rr = zeta ( 0 ) ;
    long double       tt = 1 ; 
    const long double w2 = w * w ;
    for ( unsigned short k = 1 ; 2 * k <= n ; ++k )
    {
      tt *= ( n + 2 - 2 * k ) * ( n + 1 - 2 * k ) / w2 ;
      rr = tt * zeta ( 2 * k ) ; 
    }
    const double f1 = std::pow ( w , n - 1 ) * igamma ( n - 2 ) ;

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