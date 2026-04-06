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
namespace
{
  // ==========================================================================
  template <class TYPE1, 
            class TYPE2, 
            class ZERO , 
            class EQUAL>
  double Li_power  
  ( const TYPE1       s         , 
    const TYPE2       x         , 
    const ZERO        cmp_zero  ,
    const EQUAL       cmp_equal , 
    const std::size_t NN  = 100 ) 
  {
    /// mandatory condition for convergency
    Ostap::Assert ( std::abs ( x ) < 1  , 
                    "Argument is not small enough to use the power series!" , 
                    "::Li_power"  , 
                    INVALID_ARGUMENT , __FILE__ , __LINE__ ) ;

    TYPE2 result      =       1 ; // NB: no x here 
    TYPE2 term        =       1 ; // NB: no x here 
    unsigned short nSmall      =       0 ;
    const bool     alternating = ( x < 0 ) ;
    //
    for ( std::size_t k = 2 ; NN ; ++ k )
    {
      term *= x ;
      TYPE2 delta = term * std::pow ( 1.0L / k , s ) ;
      //
      if ( !delta || cmp_zero ( 2 * delta ) || cmp_equal ( result + 2 * delta , result ) ) { ++nSmall     ; }
      else                                                                                 {   nSmall = 0 ; }
      //
      result += delta ;
      //
      // Alternating series ?
      // - one "small" term is enough,
      // - otherwise requre several  consequitive small terms are requred 
      if ( ( alternating && nSmall ) || ( 2 <= nSmall ) ) 
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
  ( const short       n , 
    const long double x )
  {
    const std::size_t NN  = 100 + std::max ( 0 , -2 * n ) ;
    return Li_power ( n , 1.0L * x , s_zero , s_equal , NN ) ;
  }
  // ===========================================================================
  inline long double Li_eq_9_2
  ( const short       n , 
    const long double x , 
    const long double w ) 
  {
    Ostap::Assert ( x < 0 && std::abs( w / s_pi ) < 1 ,
                  "Invalid ranage for x or w" , 
                  "::Li_eq_9_2_"  , 
                  INVALID_ARGUMENT , __FILE__ , __LINE__ );

    long double    result = -1 * Ostap::Math::eta ( n ) ;
    long double    term   =  1 ;
    unsigned short nSmall =  0 ;  
    for ( unsigned int k = 1  ; k < 1000 ; ++k )
    {
      term /=  ( w / k  ) ;
      const int         nmk   = 1 * n - k ;
      const long double eta_v = Ostap::Math::eta ( nmk ) ;
      if ( !eta_v )  { continue ; } 
      const long double delta = term * eta_v ;
      // 
      if ( !delta || s_zero ( 2 * delta ) ||  s_equal ( result - 2 * delta , result ) ) { ++nSmall     ; }
      else                                                                              {   nSmall = 0 ; }
      //
      result -= delta ; 
      // 
      if ( 3 <= nSmall ) 
      {
	      std::cerr << "POWER-W (9.2) #terms " << k << " " << delta << " " << x << "# " << nSmall <<std::endl ;
	      break ; // NB: x here 
      }
    }
    std::cerr << "POWER-W (9.2) x=" << x << std::endl ;
    return result ;   
  }
  // ===========================================================================
  inline double Li_eq_9_3
  ( const short       n , 
    const long double x , 
    const long double w ) 
  {
    Ostap::Assert ( ( n <= 0 ) && ( x > 0 ) && std::abs( w / s_2pi ) < 1 ,
                  "Invalid ranage for n,x or w" , 
                  "::Li_eq_9_3_"  , 
                  INVALID_ARGUMENT , __FILE__ , __LINE__ );

    long double    result = Ostap::Math::zeta ( n ) ;
    long double    term   = 1 ;
    unsigned short nSmall = 0 ;  
    for ( unsigned int k = 1  ; k < 1000 ; ++k )
    {
      term /=  ( w / k  ) ;
      const int         nmk   = 1 * n - k ;
      const long double eta_v = Ostap::Math::zeta ( nmk ) ;
      if ( !eta_v )  { continue ; } 
      const long double delta = term * eta_v ;
      // 
      if ( !delta || s_zero ( 2 * delta ) ||  s_equal ( result + 2 * delta , result ) ) { ++nSmall     ; }
      else                                                                              {   nSmall = 0 ; }
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
    return result + Ostap::Math::gamma ( 1.0 + std::abs ( n ) ) * std::pow ( -w , 1 * n - 1 ) ; 
  }
  // =========================================================================
  inline double Li_eq_9_5
  ( const short       n , 
    const long double x , 
    const long double w ) 
  {
    Ostap::Assert ( ( 1 <= n ) && ( x > 0 ) && std::abs( w / s_2pi ) < 1 ,
                  "Invalid ranage for n,x or w" , 
                  "::Li_eq_9_5_"  , 
                  INVALID_ARGUMENT , __FILE__ , __LINE__ ) ;

    long double result  = Ostap::Math::zeta ( n ) ;
    long double term    = 1 ;
    unsigned int nSmall = 0 ;
    for ( unsigned int k = 1 ; k < n + 1000 ; ++k ) 
    {
      term /= ( w / k ) ; 
      if ( 1 * n  == k + 1 ) { continue ; } // skip singularity 
      const int         nmk    = 1 * n - k    ;
      const long double zeta_v = Ostap::Math::zeta ( nmk ) ;
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
	      std::cerr << "# EQ (9.5) terms " << k << " " << delta << " " << x << "# " << nSmall <<std::endl ;
	      break ;  
      }
    }
    //
    const unsigned short nm1 = n - 1 ; 
    // result += ( Ostap::Math::harmonic ( nm1 ) - std::log ( std::abs ( w ) ) ) * std::pow ( w , nm1 ) * Ostap::Math::igamma ( nm1 ) ;
    return result ;   
  }
 //const long double sB99 =  -94598037819122125295227433069493721872702841533066936133385696204311395415197247711.0L/33330;
}
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
  if      (  0 == n ) { return x                /          ( 1.0L - x )     ; } // Eq.(6.2)
  else if ( -1 == n ) { return x                / std::pow ( 1.0L - x , 2 ) ; } // Eq.(6.2)
  else if ( -2 == n ) { return x * ( x + 1.0L ) / std::pow ( 1.0L - x , 3 ) ; } // Eq.(6.2)
  // Rational case 
  else if  ( -1 * N_EULERIAN_MAX <= n && n < 0 ) 
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
  { return Li_power  ( n , 1.0 * x ) ; }

  // use Eq (9.2)
  if ( ( x < 0 )  && ( absw1pi <= 0.5 ) )
  { return Li_eq_9_2 ( n , 1.0L * x , 1.0L * w ) ; }

  
  // use Eq (9.3)
  if ( ( n < 0 ) && ( x < 0 )  && ( absw2pi <= 0.512 ) )
  { return Li_eq_9_3 ( n , 1.0L * x , 1.0L * w ) ; }


  // use Eq (9.5)
  if ( ( n > 0 ) && ( x < 0 )  && ( absw2pi <= 0.512 ) )
  { return Li_eq_9_5 ( n , 1.0L * x , 1.0L * w ) ; }



  //  (c) use Eq ( 9.5 )  
  if ( 0 < x && 2 <= n && absw2pi < 0.512 )
  {
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
// double Ostap::Math::Li
// ( const double s ,
//   const double x )
// {
//   if       ( !x || s_zero ( x ) ) { return 0  ; }
//  //  
//  Ostap::Assert ( false ,
//		  "Not yet implemented!" ,
//		  "Ostap::Math::Li"      ,
//		  NOT_IMPLEMENTED_YET , __FILE__ , __LINE__ ) ;
//  
//  return std::numeric_limits<double>::quiet_NaN() ;   
// }

// ===========================================================================
//                                                                     The END 
// ===========================================================================