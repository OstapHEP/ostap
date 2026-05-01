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
#include "Ostap/Power.h"
#include "Ostap/Differences.h"
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
// ============================================================================
namespace
{
  // ==========================================================================
  constexpr std::complex<long double>  s_J { 0.0L , 1.0L } ;
  // ==========================================================================
  constexpr long double s_Li_DELTA = 1e-3L ;
  // ==========================================================================
  static_assert ( 0 < s_Li_DELTA && s_Li_DELTA < 1 ,
		  "Invalid Li_DELTA parameter!"   ) ; 
  // ==========================================================================
  constexpr unsigned short s_borwain_JMAX = 7 ;
  static_assert ( 3 <= s_borwain_JMAX  && s_borwain_JMAX <= 10 ,
		  "Invalid s_borwain_JMAX parameter!"   ) ; 
  // ==========================================================================
  // Borwain' formulae
  // @see https://doi.org/https://doi.org/10.1016/j.jat.2014.10.004
  // page 119, Eqs (11-14)
  // ==========================================================================
  constexpr inline long double borwain_f 
  ( const unsigned int   k , 
    const unsigned short q )
  {
    if      ( !k && !q ) { return 1 ; }  // f(0,0) = 1 
    else if ( !k       ) { return 0 ; }  // f(0,q) = 0 
    else if ( !q       ) { return 1 ; }  // f(k,0) = 1 
    // 
    long double r = 0 ;
    for ( unsigned short h = 0 ; h <= q ; ++h  )
      { r += Ostap::Math::sign ( h ) * Ostap::Math::POW ( 1.0L / k , h ) * borwain_f ( k - 1 , q - h ) ; }
    return r ; 
  }
  // ==========================================================================
  template <typename  TYPE>
  inline TYPE borwain_b 
  ( const unsigned int   k ,
    const unsigned short j , 
    const TYPE&          L )
  {
    TYPE result       = 0 ;
    //
    std::vector<TYPE>   LP ( j + 1 ) ; 
    LP [ 0 ] = 1  ;
    TYPE t1  = 1 ;
    for ( unsigned short p = 1 ; p <= j ; ++p )
    {
      t1 *= ( L / ( 1.0L * p ) ) ; 
      LP [ p ] = t1 ;
    }
    //
    std::vector<long double> GP ( j + 1 ) ;
    GP [ 0 ] = Ostap::Math::dgamma_at_1 ( 0 ) ;
    double t2 = 1 ;
    for ( unsigned short t = 1 ; t <= j ; ++t )
    {
      t2 /= t ;
      GP [ t ] = Ostap::Math::dgamma_at_1 ( t ) * t2 * Ostap::Math::sign ( t );
    }
    //
    for ( unsigned short p = 0 ; p <= j ; ++p )
    {
      for ( unsigned short t = 0 ; p + t <= j ; ++t )
      {
        const unsigned short q = j - p - t ; 
        result += LP [ p ] * GP [ t ] * borwain_f ( k , q ) ;
      }
    }
    //
    return result ;
  } 
  // ==========================================================================
  template <typename  TYPE>
  inline TYPE borwain_c 
  ( const unsigned int   k ,
    const unsigned short j , 
    const TYPE&          L )
  {
    const long double t = Ostap::Math::stieltjes ( j ) * Ostap::Math::sign ( j ) * Ostap::Math::igamma ( 1 + j ) ;
    // return TYPE ( t ) - borwain_b ( k , j + 1 , L ) ;
    const auto r = TYPE ( t ) - borwain_b ( k , j + 1 , L ) ;
    //
    return r ;           
  }
  // ==========================================================================
  template <typename TYPE>
  inline TYPE borwain_Q 
  ( const unsigned int   k     ,
    const TYPE&          L     , 
    const long double    t     , 
    const unsigned short J = s_borwain_JMAX )
  {
    TYPE        result = 0 ;
    long double tau    = 1 ;
    for ( unsigned int j = 0 ; j <= J ; ++j )
    {
      result += borwain_c ( k , j , L ) * tau ;
      tau    *= t ;  
    }
    //
    return result ;
  }  
  // ==========================================================================
  inline bool alternating ( const long double            x    ) { return x <= 0 ; }
  template <typename TYPE>
  inline bool alternating ( const std::complex<TYPE>& /* z */ ) { return false  ; }
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
    TYPEX        result =  1 ; // NB: no x here 
    TYPEX        term   =  1 ; // NB: no x here 
    std::uint8_t nSmall =  0 ;
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
      if ( ( xalt && nSmall ) || ( 2 <= nSmall ) ) { break ; }
    }
    //
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
  /// use formula from Wikipedia
  template <typename TYPES ,
            typename TYPEX >
  inline TYPEX Li_zz
  ( const TYPES       s ,
    const TYPEX       z , 
    const std::size_t N ) 
  {
    const TYPEX zz = -z / ( 1.0L - z ) ;

    /// Condition for convergency
    Ostap::Assert ( std::abs ( zz ) < 0.5  , 
                    "Argument is not small enough for the z-power series!" , 
                    "::Li_zz"        , 
                    INVALID_ARGUMENT , __FILE__ , __LINE__ ) ;

    //    
    static const Ostap::Math::Equal_To<TYPEX> xequal {} ;
    static const Ostap::Math::Zero<TYPEX>     xzero  {} ; 
    //
    TYPEX        result = 0 ;
    TYPEX        term   = 1 ; 
    std::uint8_t nSmall = 0 ;
    //
    auto fun = [s]( const long double q ) -> long double 
    { return std::pow ( 1.0L / ( q + 1.0L ) , s ) ; } ;

    //
    for ( unsigned short k = 0 ; k < N ; ++k )
    {
      term *= zz ;
      const long  double diff  = Ostap::Math::Differences::forward_( fun , k ,  0 , 1 );
      const TYPEX        delta = term * ( diff * Ostap::Math::sign ( k + 1 ) ) ;
      //
      if   ( is_not ( delta ) || xzero ( 2.0L * delta ) || xequal ( result + 2.0L * delta , result ) ) { ++nSmall     ; }
      else                                                                                             {   nSmall = 0 ; }
      //
      result += delta ;
      //
      if ( 2 <= nSmall ) { return result ; } 
      }
    //
    return result ;
  }
  // ===========================================================================
  double Li_zz 
  ( const short  n , 
    const double x )
  {
    constexpr std::size_t N = 200 ; 
    return Li_zz ( n , 1.0L * x , N ) ;
  }
 // ===========================================================================
  double Li_zz 
  ( const double s , 
    const double x )
  {
    constexpr std::size_t N  = 200  ;
    return Li_zz ( 1.0L * s  , 1.0L * x , N ) ;
  }
  // ==========================================================================
  std::complex<double> Li_zz 
  ( const short                 n , 
    const std::complex<double>& z )
  {
    constexpr std::size_t N = 200 ; 
    const std::complex<long double> zz { z } ;
    auto result = Li_zz ( n , zz , N ) ;
    return std::complex<double> ( result ) ;
  }
  // ===========================================================================
  std::complex<double> Li_zz 
  ( const double                s , 
    const std::complex<double>& z )
  {
    constexpr std::size_t N = 200 ; 
    const std::complex<long double> zz { z } ;
    auto result = Li_zz ( s , zz , N ) ;
    return std::complex<double> ( result ) ;
  }
  // ============================================================================
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
    
    TYPEX        result = -1 * Ostap::Math::eta ( s ) ;
    TYPEX        term   =  1 ;
    std::uint8_t nSmall =  0 ;
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
      if ( 2 <= nSmall ) { break ; }
    }
    //
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
		    "Invalid range for w"     , 
		    "::Li_eq_9_3"             , 
		    INVALID_ARGUMENT          , __FILE__ , __LINE__ ) ;    
    //
    static const Ostap::Math::Equal_To<TYPEX> xequal {} ;
    static const Ostap::Math::Zero<TYPEX>     xzero  {} ; 
    //    
    TYPEX        result = Ostap::Math::zeta ( s ) ;
    TYPEX        term   = 1 ;
    std::uint8_t nSmall = 0 ;
    //   
    for ( unsigned int k = 1  ; k < NN ; ++k )
    {
      term *= ( w / ( 1.0L * k )  ) ;
      const TYPES nmk    = s - TYPES ( k ) ;
      const TYPEX zeta_v = Ostap::Math::zeta ( nmk ) ;
      //
      if ( is_not ( zeta_v ) )  { continue ; } 
      const TYPEX delta = term * zeta_v ;
      // 
      if ( is_not ( delta ) || xzero ( 2.0L * delta ) ||  xequal ( result + 2.0L * delta , result ) ) { ++nSmall     ; }
      else                                                                                            {   nSmall = 0 ; }
      //
      result += delta ; 
      // 
      if ( 3 <= nSmall ) { break ; }
    }
    //
    return result ;
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
    
    // series 
    const long double sum = Li_eq_9_3_ ( n , 1.0L * x , 1.0L * w , 10000 ) ;
    //
    const long double gg  = Ostap::Math::gamma ( 1.0L - n  ) ;
    //
    // get the Re ( pow ( -w , n - 1 ) )
    const long double pp  = w < 0 ?
      std::pow ( -w , n - 1 ) :
      std::pow (  w , n - 1 ) * Ostap::Math::sign ( n - 1 ) ;
    //
    return sum + gg * pp ; 
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
    /// series 
    const long double sum = Li_eq_9_3_ ( 1.0L * s  , 1.0L * x , 1.0L * w , 10000 ) ;
    //
    const long double gg  = Ostap::Math::gamma ( 1.0L - s  ) ;
    //    
    const long double pp  = w < 0 ?
      std::pow ( -w , s - 1.0L ) :
      std::pow (  w , s - 1.0L ) * std::cos ( s_pi * ( s - 1.0L ) ) ;
    //
    return sum + gg * pp ; 
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
    //
    const std::complex<long double> xx  { x } ;
    const std::complex<long double> ww  { w } ;
    const std::complex<long double> sum { Li_eq_9_3_ ( n , xx , ww , 10000 ) } ;
    //
    const long double               gg  { Ostap::Math::gamma ( 1.0L - n  ) } ;
    const std::complex<long double> pp  { std::pow ( -w , n - 1 )          } ;
    //
    return std::complex<double> (  sum + gg * pp ) ; 
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
    
    const std::complex<long double> xx  { x } ;
    const std::complex<long double> ww  { w } ;
    const std::complex<long double> sum { Li_eq_9_3_ ( 1.0L * s , xx , ww  , 10000 ) } ;
    //    
    const long double               gg  { Ostap::Math::gamma ( 1.0L - s  ) } ;
    const std::complex<long double> pp  { std::pow ( -w , s - 1          ) } ;
    //
    return std::complex<double> ( sum + gg * pp ) ;
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
    
    const std::complex<long double> xx  { x } ;
    const std::complex<long double> ww  { w } ;
    const std::complex<long double> ss  { s } ;
    const std::complex<long double> sum { Li_eq_9_3_ ( ss , xx , ww , 10000 ) } ;
    //
    const std::complex<long double> gg  { Ostap::Math::gamma ( 1.0L - ss ) } ;
    const std::complex<long double> pp  { std::pow ( -w , ss - 1.0L      ) } ;
    //
    return std::complex<double> ( sum + gg * pp ) ;
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
    TYPEX        result  = ( 1 < n ) ? Ostap::Math::zeta ( n ) : 0.0L ; 
    TYPEX        term    = 1 ;
    std::uint8_t nSmall = 0 ;
    //
    for ( unsigned int k = 1 ; k < N ; ++k ) 
    {
      term *= ( w / ( 1.0L * k ) ) ; 
      if ( 1 * n  == k + 1 ) { continue ; } // skip singularity 
      const int   nmk    = 1 * n - k    ;
      const TYPEX zeta_v = Ostap::Math::zeta ( nmk ) ;
      if ( is_not ( zeta_v ) )  { continue ; }
      //
      const TYPEX delta = term * zeta_v ;
      //
      if ( is_not ( delta ) || xzero ( 2.0L * delta ) ||  xequal ( result + 2.0L * delta , result ) ) { ++nSmall     ; }
      else                                                                                            {   nSmall = 0 ; }
      //
      result += delta ; 
      // 
      if ( 2 <= nSmall ) { break ; }
    }
    //
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
    const long double    sum = Li_eq_9_5_ ( n , 1.0L * x , 1.0L * w , 10000 ) ;
    //
    const unsigned short nm1  = n - 1 ;
    const long double    hm  =  Ostap::Math::harmonic ( nm1 )  ;
    const long double    ig  =  Ostap::Math::igamma   ( n   )  ;
    //
    const long double    t1 =  w < 0 ? 
      hm - std::log ( -w            )         :
      hm - std::log ( -w + 0.0L*s_J ).real()  ;
    //
    const long double    t2 = std::pow ( w , nm1 ) * ig ;
    //
    return sum + t1 * t2 ;
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
    const std::complex<long double> sum { Li_eq_9_5_ ( n , xx , ww , 10000 ) } ;
    //
    const unsigned short nm1 = n - 1 ;
    const long double    hm  =  Ostap::Math::harmonic ( nm1 )  ;
    const long double    ig  =  Ostap::Math::igamma   ( n   )  ;
    //
    const std::complex<long double> t1 = hm - std::log ( -ww ) ;
    const std::complex<long double> t2 = std::pow ( ww , nm1 ) * ig ; 
    //
    return std::complex<double> ( sum + t1 * t2 ) ;
  }
  // ==========================================================================
  /// Li ( k + 1 + tau , e^w ) 
  long double Li_eq_BB_11
  ( const unsigned int k   ,
    const long double  tau ,
    const long double  w   )
  {
    //
    Ostap::Assert ( std::abs( w ) / s_2pi < 1 ,
		    "Invalid range for w"     , 
		    "::Li_eq_BB_11"           , 
		    INVALID_ARGUMENT          , __FILE__ , __LINE__ ) ;
    //
    typedef long double TYPE ;
    static const Ostap::Math::Equal_To<TYPE> xequal {} ;
    static const Ostap::Math::Zero<TYPE>     xzero  {} ; 
    //
    TYPE         result = 0 ;
    TYPE         term   = 1 ;
    std::uint8_t nSmall = 0 ;
    //
    for ( unsigned int n  = 0 ; n <= k + 1000 ; ++n )
    {
      if ( n      ) { term *= ( w / ( 1.0L * n ) ) ; }
      if ( n == k ) { continue ; }                      // ATTENTION!
      const TYPE arg_zeta = 1.0L + k + tau - n ;
      const TYPE zeta_v   = Ostap::Math::zeta ( arg_zeta ) ;
      if ( is_not ( zeta_v ) ) { continue ; }
      const TYPE delta  = term * zeta_v   ;
      //      
      if ( is_not ( delta ) || xzero ( 2.0L * delta ) ||  xequal ( result - 2.0L * delta , result ) ) { ++nSmall     ; }
      else                                                                                            {   nSmall = 0 ; }
      //
      result += delta ;
      //
      if ( 2 <= nSmall ) { break ; }      
    }
    //
    const long double pp = Ostap::Math::POW    ( w , k ) ; 
    const long double gg = Ostap::Math::igamma ( k + 1 ) ;  
    //
    if ( 0 <= w )
    {
      const std::complex<long double> L { std::log ( -w + 0.0L * s_J ) }  ;
      const long double               T { tau } ;      
      const std::complex<long double> Q = borwain_Q ( k , L , T , s_borwain_JMAX ) ;
      return result + pp * gg * Q.real() ;
    }
    //
    const long double L = std::log ( w ) ;    
    const long double T = tau            ;    
    const long double Q = borwain_Q ( k , L , T , s_borwain_JMAX ) ;
    //
    return result + pp * gg * Q ;
  }
  // ==========================================================================
  /// Li ( k + 1 + tau , e^w ) 
  std::complex<double> Li_eq_BB_11
  ( const unsigned int          k    ,
    const long double           tau  ,
    const std::complex<double>& w_   )
  {
    //
    Ostap::Assert ( std::abs ( w_ ) / s_2pi < 1 ,
		    "Invalid range for w"     , 
		    "::Li_eq_BB_11"           , 
		    INVALID_ARGUMENT          , __FILE__ , __LINE__ ) ;
    //    
    typedef std::complex<long double> TYPE ;
    //
    static const Ostap::Math::Equal_To<TYPE> xequal {} ;
    static const Ostap::Math::Zero<TYPE>     xzero  {} ; 
    //
    const TYPE w   { w_   } ;
    //
    TYPE         result = 0 ;
    TYPE         term   = 1 ;
    std::uint8_t nSmall = 0 ;
    //
    for ( unsigned int n  = 0 ; n <= k + 1000 ; ++n )
    {
      //
      if ( n      ) { term *= ( w / ( 1.0L * n ) ) ; }
      if ( n == k ) { continue ; }                      // ATTENTION!
      const long double arg_zeta = 1.0L + k - n + tau ;
      const long double zeta_v   = Ostap::Math::zeta ( arg_zeta ) ;
      if ( is_not ( zeta_v ) ) { continue ; }
      const TYPE delta  = term * zeta_v   ;
      //
      if ( is_not ( delta ) || xzero ( 2.0L * delta ) ||  xequal ( result - 2.0L * delta , result ) ) { ++nSmall     ; }
      else                                                                                            {   nSmall = 0 ; }
      //
      result += delta ;
      //
      if ( 2 <= nSmall ) { break ; }      
    }
    //
    const TYPE        pp = Ostap::Math::POW    ( w , k ) ; 
    const long double gg = Ostap::Math::igamma ( k + 1 ) ;  
    //
    const TYPE        L = std::log  ( - w ) ;    
    const TYPE        Q = borwain_Q ( k , L , tau , s_borwain_JMAX ) ;
    //
    return std::complex<double> ( result + pp * gg * Q ) ;
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
  if      ( 1  == n ) // Eq (6.1)
  {
    return x < 1 ?
      - std::log (                            1.0L - x       )         :
      - std::log ( std::complex<long double>( 1.0L - x , 0 ) ).real () ;
  }
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
    //
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
  if ( std::abs ( x ) <= 0.25 || xsmall_1 || xsmall_2 ) { return Li_power  ( n , x ) ; }

  /// (B) use Eq (9.2)
  if ( ( x < 0 )  && ( absw1pi <= 0.5 ) )
  { return Li_eq_9_2 ( n , 1.0L * x , 1.0L * w ) ; }
  
  /// (C) use Eq (9.3)
  if ( ( n < 0 ) && ( 0 < x ) && ( absw2pi <= 0.512 ) )
  { return Li_eq_9_3 ( n , 1.0L * x , 1.0L * w ) ; }
  
  /// (D) use Eq (9.5)
  if ( ( n > 0 ) && ( 0 < x ) && ( absw2pi <= 0.512 ) )
  { return Li_eq_9_5 ( n , 1.0L * x , 1.0L * w ) ; }

  /// (E) Eq. (10.3)
  if ( 0 > n && 1 < std::abs ( x ) )
  { return ( 0 == n % 2 ? -1 : +1 ) * Li ( n , 1.0 / x ) ; }
  
  /// (F) Eq. (10.1)
  if ( 1 <= n && 1 < x ) 
  {
    long double          rr = zeta ( 0 ) ;
    long double          tt = 1 ; 
    const long double    w2 = w * w ;
    for ( unsigned short k  = 1 ; 2 * k <= n ; ++k )
    {
      tt *= ( n + 2 - 2 * k ) * ( n + 1 - 2 * k ) / w2 ;
      rr += tt * zeta ( 2 * k ) ;      
    }
    const long double R = 2 * std::pow ( w , n ) * igamma ( n + 1 ) * rr ; 
    return R + ( 0 == n % 2 ? -1 : +1 ) * Li ( n , 1.0 / x ) ;
  }

  /// (G) use formula from Wiki-page, it has some useful range 
  if ( x < 0.5 )
  {
    const long double xx = -x / ( 1.0L - x ) ;
    /// use 
    if ( std::abs ( xx ) < 0.5 ) { return Li_zz  ( n , x ) ;  }
  }


  /// (H) negative argument? (Eq 14.1)
  if ( x < 0  )
  {
    // Use the square formula Eq.(14.1) to convert to the positive argument
    return std::pow ( 2 , 1 - n ) * Li ( n , x * x ) - Li ( n , -x ) ;
  }
  
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
  //
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
  // Here s is *not* (short) integer
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
  
  /// (C) use Eq (9.3) or (9.4) 
  if ( ( 0 < x ) && ( absw2pi <= 0.512 ) )
  {
    /// Is s close to the positive integer ?
    const bool sposint =
      ( 1 -  s_Li_DELTA ) < s                                        &&
      ( s <= s_Li_DELTA + std::numeric_limits<unsigned int>::max() ) && 
      std::abs ( s - round ( s ) ) <= s_Li_DELTA    ;
    //
    if ( !sposint ) { return Li_eq_9_3 ( 1.0L * s , 1.0L * x , 1.0L * w ) ;}
    
    /// close to positive integer
    /// From here follow Bailey & Borwain, p 119 Eq(11)
    const unsigned int kp1 = round ( s )    ;
    const unsigned int k   = kp1 - 1        ;
    const long double  tau = 1.0L * s - kp1 ;
    //
    return Li_eq_BB_11 ( k , tau , w ) ; 
    //
  }
  //

  /// (G) use formula from Wiki-page, it has some useful range 
  if ( x < 0.5 )
  {
    const long double xx = -x / ( 1.0L - x ) ;
    /// use 
    if ( std::abs ( xx ) < 0.5 ) { return Li_zz  ( s , x ) ;  }
  }

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
/*  Polylogarithm function  \f$ Li_n(x)  = \sum \frac{x^k}{k^n} f\$
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
  // REAL case ?
  if ( !y || s_zero ( y ) || s_equal ( x + y , x ) )
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


  /// (G) use formula from Wiki-page, it has some useful range 
  if ( z.real () < 0.5 )
  {
    const std::complex<long double> zw = -zz / ( 1.0L - zz ) ;
    /// use 
    if ( std::abs ( zw ) < 0.5 ) { return Li_zz  ( n , z ) ;  }
  }


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
/*  Polylogarithm function  \f$ Li_s(x)  = \sum \frac{x^k}{k^s} f\$
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
  if ( !y || s_zero ( y ) || s_equal ( x + y , x ) )
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
  if (  absw2pi <= 0.512 )
  {
    /// close to positive integer ? 
    const bool sposint =
      ( 1 -  s_Li_DELTA ) < s                                        &&
      ( s <= s_Li_DELTA + std::numeric_limits<unsigned int>::max() ) && 
      std::abs ( s - round ( s ) ) <= s_Li_DELTA    ;
    //
    if ( !sposint ) { return Li_eq_9_3 ( 1.0L * s , z , w  ) ; }
    //
    
    /// close to positive integer
    /// From here follow Bailey & Borwain, p 119 Eq(11)
    const unsigned int    kp1 = round ( s ) ;
    const unsigned int    k   = kp1 - 1     ;
    const long     double tau =  s - kp1    ;
    //
    return Li_eq_BB_11 ( k , tau , w ) ; 
  }
  
  /// (G) use formula from Wiki-page, it has some useful range 
  if ( z.real () < 0.5 )
  {
    const std::complex<long double> zz { z } ; 
    const std::complex<long double> zw = - zz / ( 1.0L - zz ) ;
    /// use 
    if ( std::abs ( zw ) < 0.5 ) { return Li_zz  ( s , z ) ;  }
  }

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

// ===========================================================================
/*  (Real part of) inverse tangent integral
 *  \f[ Ti_s(z) = \frac{Li_s(iz) - Li_s(-iz)}{2i} \right)\f]
 *  @See https://en.wikipedia.org/wiki/Inverse_tangent_integral
 */ 
// ===========================================================================
double Ostap::Math::Ti
( const short  n ,
  const double x )
{
  //
  if      ( 0 == n ) { return x / ( 1 + x * x ) ; }
  else if ( 1 == n ) { return std::atan ( x )   ; }
  //
  const std::complex<double> z1 { 0.0 , +x } ;
  const std::complex<double> z2 { 0.0 , -x } ;
  //
  const std::complex<double> l1 { Li ( n , z1 ) } ;
  const std::complex<double> l2 { Li ( n , z2 ) } ;
  //
  return 0.5 * std::imag ( l1 - l2 ) ; 
}
// ===========================================================================
/*  (Real part of) inverse tangent integral
 *  \f[ Ti_s(z) = \frac{Li_s(iz) - Li_s(-iz)}{2i} \right)\f]
 *  @See https://en.wikipedia.org/wiki/Inverse_tangent_integral
 */ 
// ===========================================================================
double Ostap::Math::Ti
( const double s ,
  const double x )
{
  if ( isshort ( s ) )
  {
    const short n = round ( s ) ;
    return Ti ( n , x ) ;
  }
  //
  const std::complex<double> z1 { 0.0 , +x } ;
  const std::complex<double> z2 { 0.0 , -x } ;
  //
  const std::complex<double> l1 { Li ( s , z1 ) } ;
  const std::complex<double> l2 { Li ( s , z2 ) } ;
  //
  return 0.5 * std::imag ( l1 - l2 ) ; 
}
// ===========================================================================
/*  Inverse tangent integral
 *  \f[ Ti_s(z) = \frac{Li_s(iz) - Li_s(-iz)}{2i} \right)\f]
 *  @See https://en.wikipedia.org/wiki/Inverse_tangent_integral
 */ 
// ===========================================================================
std::complex<double>
Ostap::Math::Ti
( const short                 n ,
  const std::complex<double>& z )
{
  //
  if      ( 0 == n ) { return z / ( 1.0 + z * z ) ; }
  else if ( 1 == n ) { return std::atan ( z )   ; }
  //
  const std::complex<double> z1 { +s_j * z } ;
  const std::complex<double> z2 { -s_j * z } ;
  //
  const std::complex<double> l1 { Li ( n , z1 ) } ;
  const std::complex<double> l2 { Li ( n , z2 ) } ;
  //
  return 0.5 * ( l1 - l2 ) / s_j ; 
}
// ===========================================================================
/*  Inverse tangent integral
 *  \f[ Ti_s(z) = \frac{Li_s(iz) - Li_s(-iz)}{2i} \right)\f]
 *  @See https://en.wikipedia.org/wiki/Inverse_tangent_integral
 */ 
// ===========================================================================
std::complex<double>
Ostap::Math::Ti
( const double                s ,
  const std::complex<double>& z )
{
  if ( isshort ( s ) )
  {
    const short n = round ( s ) ;
    return Ti ( n , z ) ;
  }
  // 
  const std::complex<double> z1 { +s_j * z } ;
  const std::complex<double> z2 { -s_j * z } ;
  //
  const std::complex<double> l1 { Li ( s , z1 ) } ;
  const std::complex<double> l2 { Li ( s , z2 ) } ;
  //
  return 0.5 * ( l1 - l2 ) / s_j ; 
}
// ===========================================================================

// ===========================================================================
/* (Real part of) Legendre chi-function
 *  \f[ \chi_s(z) = \frac{Li_s(z) - Li_s(-z)}{2} \f] 
 *  https://en.wikipedia.org/wiki/Legendre_chi_function
 */
// ===========================================================================
double Ostap::Math::legendre_chi
( const short  n ,
  const double x )
{
  if ( 0 == n ) { return x / ( 1 - x * x ) ; }
  return 0.5 * ( Li ( n , +x ) - Li ( n , -x ) ) ;
}
// ===========================================================================
/*  (Real part of) Legendre chi-function
 *  \f[ \chi_s(z) = \frac{Li_s(z) - Li_s(-z)}{2} \f] 
 *  https://en.wikipedia.org/wiki/Legendre_chi_function
 */
// ===========================================================================
double Ostap::Math::legendre_chi
( const double s ,
  const double x ) 
{
  if ( isshort ( s ) )
  {
    const short n = round ( s ) ;
    return legendre_chi ( n , x ) ;
  }
  //
  return 0.5 * ( Li ( s , +x ) - Li ( s , -x ) ) ;
}
// ===========================================================================
/*  Legendre chi-function
 *  \f[ \chi_s(z) = \frac{Li_s(z) - Li_s(-z)}{2} \f] 
 *  https://en.wikipedia.org/wiki/Legendre_chi_function
 */
// ===========================================================================
std::complex<double>
Ostap::Math::legendre_chi
( const short                 n ,
  const std::complex<double>& z )
{
  if ( 0 == n ) { return z / ( 1.0 - z * z ) ; }
  return 0.5 * ( Li ( n , +z ) - Li ( n , -z ) ) ;  
}
// ===========================================================================
/*  Legendre chi-function
 *  \f[ \chi_s(z) = \frac{Li_s(z) - Li_s(-z)}{2} \f] 
 *  https://en.wikipedia.org/wiki/Legendre_chi_function
 */
// ===========================================================================
std::complex<double>
Ostap::Math::legendre_chi
( const double                s ,
  const std::complex<double>& z )
{
  if ( isshort ( s ) )
  {
    const short n = round ( s ) ;
    return legendre_chi ( n , z ) ;
  }
  //
  return 0.5 * ( Li ( s , +z ) - Li ( s , -z ) ) ;  
}
// ===========================================================================

// ===========================================================================
//                                                                     The END 
// ===========================================================================
