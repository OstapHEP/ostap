// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================
// GSL 
// ============================================================================
#include "gsl/gsl_sf_clausen.h"
// ============================================================================
// LHCbMath
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/MakeArray.h"
#include "Ostap/Clenshaw.h"
#include "Ostap/Bernulli.h"
#include "Ostap/Clausen.h"
#include "Ostap/PolyLog.h"
#include "Ostap/Gamma.h"
// ============================================================================
// Local
// ============================================================================
#include "GSL_sentry.h"
#include "local_math.h"
#include "local_gsl.h"
// ============================================================================
/** @file
 *  implementation file for Clausen functions 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-03-26
 */
// ============================================================================
// initialize the static arrays&objects 
// ============================================================================
namespace
{
  // ==========================================================================
  // precomputed scaled  bernulli's polynomials 
  const std::array<long double,2>  s_S1 { -0.50L      ,  s_pi / 2.0L  } ;
  const std::array<long double,3>  s_C2 {  0.25L      , -s_pi / 2.0L  ,  s_pi * s_pi / 6.0L  } ;
  const std::array<long double,4>  s_S3 {  1.0L /  12 , -s_pi / 4.0L  ,  s_pi * s_pi / 6.0L  , 0.0L } ;
  const std::array<long double,5>  s_C4 { -1.0L /  48 ,  s_pi / 12.0L , -s_pi * s_pi / 12.0L , 0.0L , std::pow ( s_pi , 4 ) / 90.0L } ;
  const std::array<long double,6>  s_S5 { -1.0L / 240 ,  s_pi / 48.0L , -s_pi * s_pi / 36.0L , 0.0L , std::pow ( s_pi , 4 ) / 90.0L , 0.0L } ;
  // ==========================================================================
} //                                            The END of anonymous namespace 
// ============================================================================
// S-functions 
// ============================================================================
double Ostap::Math::Clausen::S0 ( const double x )
{ return 0.5 * std::tan ( 0.5 * s_pi - x * 0.5L ) ; }
// ============================================================================
double Ostap::Math::Clausen::S1 ( const double x )
{  
  const long double y = Ostap::Math::reduce ( 1.0L * x , 0.0L , s_2pi ).first ;
  return Ostap::Math::Clenshaw::monomial_sum ( s_S1.begin () , s_S1.end () , y ) . first ;
}
// ============================================================================
double Ostap::Math::Clausen::S2 ( const double x )
{ return Ostap::Math::clausen ( x ) ; }
// ============================================================================
double Ostap::Math::Clausen::S3 ( const double x )
{  
  const long double y = Ostap::Math::reduce ( 1.0L * x , 0.0L , s_2pi ).first ;
  return Ostap::Math::Clenshaw::monomial_sum ( s_S3.begin () , s_S3.end () , y ) . first ;
}
// ============================================================================
double Ostap::Math::Clausen::S4 ( const double x )
{  
  const double y = Ostap::Math::reduce ( 1.0L * x , -1.0L * s_pi , +1.0L * s_pi ).first ;
  static const short n4 { 4 } ;
  return Ostap::Math::Li ( n4 , std::exp ( y * s_j ) ) .imag () ;
}
// ============================================================================
double Ostap::Math::Clausen::S5 ( const double x )
{  
  const long double y = Ostap::Math::reduce ( 1.0L * x , 0.0L , s_2pi ).first ;
  return Ostap::Math::Clenshaw::monomial_sum ( s_S5.begin () , s_S5.end () , y ) . first ;
}
// ============================================================================
double Ostap::Math::Clausen::S6 ( const double x )
{
  const double y = Ostap::Math::reduce ( 1.0L * x , -1.0L * s_pi , +1.0L * s_pi ).first ;
  static const short n6 { 6 } ;
  return Ostap::Math::Li ( n6 , std::exp ( y * s_j ) ) .imag () ;
}
// =============================================================================
double Ostap::Math::Clausen::S7 ( const double x )
{
  const unsigned int m = 4 ;
  static const Ostap::Math::Bernulli s_B { 7 } ;
  static const long double           s_C { std::pow ( 2.0L * s_pi , 7 ) / ( 2 * 5040 ) } ;
  //
  const long double y = Ostap::Math::reduce ( 1.0L * x , 0.0L , s_2pi ).first ;
  return s_C * s_B ( y * s_1_2pi ) ;
}
// ============================================================================
double Ostap::Math::Clausen::S8 ( const double x )
{  
  const double y = Ostap::Math::reduce ( 1.0L * x , -1.0L * s_pi , +1.0L * s_pi ).first ;
  static const short n8 { 8 } ;
  return Ostap::Math::Li ( n8 , std::exp ( y * s_j ) ) .imag () ;
}
// =============================================================================
double Ostap::Math::Clausen::S9 ( const double x )
{
  const unsigned int m = 5 ;
  static const Ostap::Math::Bernulli s_B { 9 } ;
  static const long double           s_C { std::pow ( 2.0L * s_pi , 9 ) / ( 2 * 362880L ) } ;
  //
  const long double y = Ostap::Math::reduce ( 1.0L * x , 0.0L , s_2pi ).first ;
  return - s_C * s_B ( y * s_1_2pi ) ;
}
// ============================================================================
double Ostap::Math::Clausen::S10 ( const double x )
{  
  const double y = Ostap::Math::reduce ( 1.0L * x , -1.0L * s_pi , +1.0L * s_pi ).first ;
  static const short n10 { 10 } ;
  return Ostap::Math::Li ( n10 , std::exp ( y * s_j ) ) .imag () ;
}
// ============================================================================
// C-functions 
// ============================================================================
double Ostap::Math::Clausen::C0 ( const double /* x */ ) { return - 0.5 ; }
// ============================================================================
double Ostap::Math::Clausen::C1 ( const double x )
{ return - std::log ( std::abs ( 2 * std::sin ( 1.0L * x / 2 ) ) ) ; }
// ============================================================================
double Ostap::Math::Clausen::C2 ( const double x )
{  
  const long double y = Ostap::Math::reduce ( 1.0L * x , 0.0L , 2.0L * s_pi ).first ;
  return Ostap::Math::Clenshaw::monomial_sum ( s_C2.begin () , s_C2.end () , y ) . first ;
}
// ============================================================================
double Ostap::Math::Clausen::C3 ( const double x )
{  
  const double y = Ostap::Math::reduce ( 1.0L * x , -1.0L * s_pi , +1.0L * s_pi ).first ;
  static const short n3 { 3 } ;
  return Ostap::Math::Li ( n3 , std::exp ( y * s_j ) ) .real () ;
  // 
  // const long double y = Ostap::Math::reduce ( 1.0L * x , -1.0L * s_pi , +1.0L * s_pi ).first ;
  // return s_c3 ( std::abs ( y ) ) ;
}
// ============================================================================
double Ostap::Math::Clausen::C4 ( const double x )
{  
  const long double y = Ostap::Math::reduce ( 1.0L * x , 0.0L , s_2pi ).first ;
  return Ostap::Math::Clenshaw::monomial_sum ( s_C4.begin () , s_C4.end () , y ) . first ;
}
// ============================================================================
double Ostap::Math::Clausen::C5 ( const double x )
{  
  const double y = Ostap::Math::reduce ( 1.0L * x , -1.0L * s_pi , +1.0L * s_pi ).first ;
  static const short n5 { 5 } ;
  return Ostap::Math::Li ( n5 , std::exp ( y * s_j ) ) .real () ;
  // 
  // const long double y = Ostap::Math::reduce ( 1.0L * x , -1.0L * s_pi , +1.0L * s_pi ).first ;
  // return s_c5 ( std::abs ( y ) ) ;
}
// ============================================================================
double Ostap::Math::Clausen::C6 ( const double x )
{
  const unsigned int m = 3 ;
  static const Ostap::Math::Bernulli s_B { 2 * m  } ;
  static const long double           s_C { std::pow ( 2.0L * s_pi , 6 ) / ( 2 * 720 ) } ;
  //
  const long double y = Ostap::Math::reduce ( 1.0L * x , 0.0L , 2.0L * s_pi ).first ;
  return s_C * s_B ( y / ( 2.0L * s_pi ) ) ;
}
// ============================================================================
double Ostap::Math::Clausen::C7 ( const double x )
{  
  const double y = Ostap::Math::reduce ( 1.0L * x , -1.0L * s_pi , +1.0L * s_pi ).first ;
  static const short n7 { 7 } ;
  return Ostap::Math::Li ( n7 , std::exp ( y * s_j ) ) .real () ;
}
// =============================================================================
double Ostap::Math::Clausen::C8 ( const double x )
{
  constexpr unsigned int m = 4 ;
  static const Ostap::Math::Bernulli s_B { 2 * m  } ;
  static const long double           s_C { std::pow ( s_2pi , 2 * m ) / ( 2 * 40320 ) } ;
  //
  const long double y = Ostap::Math::reduce ( 1.0L * x , 0.0L , s_2pi ).first ;
  return - s_C * s_B ( y * s_1_2pi ) ;
}
// ============================================================================
double Ostap::Math::Clausen::C9 ( const double x )
{  
  const double y = Ostap::Math::reduce ( 1.0L * x , -1.0L * s_pi , +1.0L * s_pi ).first ;
  static const short n9 { 9 } ;
  return Ostap::Math::Li ( n9 , std::exp ( y * s_j ) ) .real () ;
}
// ============================================================================
double Ostap::Math::Clausen::C10 ( const double x )
{
  constexpr unsigned int m = 5 ;
  static const Ostap::Math::Bernulli s_B { 2 * m } ;
  static const long double           s_C { std::pow ( s_2pi , 2 * m ) / ( 2 * 3628800L ) } ;
  //
  const long double y = Ostap::Math::reduce ( 1.0L * x , 0.0L , s_2pi ).first ;
  return s_C * s_B ( y * s_1_2pi ) ;
}
// ============================================================================

// ============================================================================
// \f$ \sum_{i=1} \frac{ \sin kx }{ k^n }\f$ 
// ============================================================================
double Ostap::Math::Clausen::S
( const unsigned int n ,
  const double       x )
{
  // 
  switch ( n )
    {
    case  0 : return S0  ( x ) ;
    case  1 : return S1  ( x ) ;
    case  2 : return S2  ( x ) ;
    case  3 : return S3  ( x ) ;
    case  4 : return S4  ( x ) ;
    case  5 : return S5  ( x ) ;
    case  6 : return S6  ( x ) ;
    case  7 : return S7  ( x ) ;
    case  8 : return S8  ( x ) ;
    case  9 : return S9  ( x ) ;
    case 10 : return S10 ( x ) ;
    default : break ; 
    }
  //
  // Bernulli's polynomial
  if ( n < 20 && 1 == n % 2 )
  {
    const unsigned short        nn = n    ; 
    const unsigned short        m  = ( nn + 1 ) / 2 ;
    const Ostap::Math::Bernulli b  { nn } ;
    const double                c = 0.5L * std::pow ( s_2pi , n ) / factorial ( n ) ;
    const double                y = Ostap::Math::reduce ( 1.0L * x , 0.0L , s_2pi ).first ;
    return sign ( m + 1 ) * c * b ( y * s_1_2pi ) ;    
  }
  //
  // use polylogarithm
  if ( 50 <= n ) { return Ostap::Math::Li ( 1.0 * n , std::exp ( x * s_j ) ).imag () ; }
  //
  // use straightforward summation of Fourier sine series  
  constexpr std::size_t N = 35 ;
  const std::array<long double,N> s_C { Ostap::Math::make_array
    ( [ n ] ( const std::size_t  k ) -> long double  
    { return std::pow ( 1.0L / ( k + 1 ) , n ) ; } ,
      std::make_index_sequence<N>() ) } ;
  //
  return Ostap::Math::Clenshaw::sine_sum ( s_C.begin() , s_C.end() , x ) ; 
}
// ============================================================================
// \f$ \sum_{i=1} \frac{ \cos kx }{ k^n }\f$ 
// ============================================================================
double Ostap::Math::Clausen::C
( const unsigned int n ,
  const double       x )
{
  switch ( n )
    {
    case  0 : return C0  ( x ) ;
    case  1 : return C1  ( x ) ;
    case  2 : return C2  ( x ) ;
    case  3 : return C3  ( x ) ;
    case  4 : return C4  ( x ) ;
    case  5 : return C5  ( x ) ;
    case  6 : return C6  ( x ) ;
    case  7 : return C7  ( x ) ;
    case  8 : return C8  ( x ) ;
    case  9 : return C9  ( x ) ;
    case 10 : return C10 ( x ) ;
    default : break ; 
    }
  //
  // Bernulli's polynomial
  if ( n < 20 && 0 == n % 2 )
  {
    const unsigned short        nn = n     ; 
    const unsigned short        m  =nn / 2 ;
    const Ostap::Math::Bernulli b { nn }   ;
    const double                c = 0.5L * std::pow ( s_2pi , n ) / factorial ( n ) ;
    const double                y = Ostap::Math::reduce ( 1.0L * x , 0.0L , s_2pi ).first ;
    return sign ( m + 1 ) * c * b ( y * s_1_2pi ) ;    
  }
  //
  // use polylogarithm
  if ( 50 <= n ) { return Ostap::Math::Li ( 1.0 * n , std::exp ( x * s_j ) ).real () ; }
  //  
  // use straightforward summation of Fourier cosine series  
  constexpr std::size_t N = 35 ;
  const std::array<long double,N> s_C { Ostap::Math::make_array
    ( [ n ] ( const std::size_t  k ) -> long double  
    { return 0 == k ? 0.0 : std::pow ( 1.0L / k , n ) ; } ,
      std::make_index_sequence<N>() ) } ;
  //
  return Ostap::Math::Clenshaw::cosine_sum ( s_C.begin() , s_C.end() , x ) ; 
}
// ============================================================================
/*  Generalized Clausen' function
 *  \f[ S_s(x) = \Im Li ( s , \mathrm{e}^{ix} ) = \sum \frac{ \sin kx}{k^s}\f]
 */
// ============================================================================
double Ostap::Math::Clausen::S
( const double s ,
  const double x )
{
  if ( isuint ( s ) )
  {
    const unsigned int n = round ( s ) ; 
    return S ( n , x ) ;
  }
  // 
  return Li ( s , std::exp ( x * s_j ) ) .imag () ;  
}
// ============================================================================
/*  Generalized Clausen' function
 *  \f[ C_s(x) = \Re Li ( s , \mathrm{e}^{ix} ) = \sum \frac{ \sin kx}{k^s}\f]
 */
// ============================================================================
double Ostap::Math::Clausen::C
( const double s ,
  const double x )
{
  if ( isuint ( s ) )
  {
    const unsigned int n = round ( s ) ;
    return C ( n , x ) ;
  }
  // 
  return Li ( s , std::exp ( x * s_j ) ) .real () ; 
}
// ============================================================================
/* standard Clausen functions
 *  \f[ \begin{array}{lcc}
 *      Cl_{2m+2} ( x ) & = \sum_k \frac{ \sin kx }{k^{2m+2}}& \\ 
 *      Cl_{2m+1} ( x ) & = \sum_k \frac{ \cos kx }{k^{2m+1}}& \\ 
 *      \end{array}   \f] 
 */
// ============================================================================
double Ostap::Math::Cl
( const unsigned int n , 
  const double       x )
{
  return 0 == n % 2 ?
    Ostap::Math::Clausen::S ( n , x ) :
    Ostap::Math::Clausen::C ( n , x ) ;     
}
// ============================================================================
/*  standard Clausen functions, aka Gleisher-Clausen fnuctions  
 *  \f[ \begin{array}{lcc}
 *      Sl_{2m+2} ( x ) & = \sum_k \frac{ \cos kx }{k^{2m+2}}& \\ 
 *      Sl_{2m+1} ( x ) & = \sum_k \frac{ \sin kx }{k^{2m+1}}& \\ 
 *      \end{array}   \f] 
 * The function are related to Bernulli polynomials 
 */
// ============================================================================
double Ostap::Math::Sl
( const unsigned int  n , 
  const double        x )
{
  return 0 == n % 2 ?
    Ostap::Math::Clausen::C ( n , x ) :
    Ostap::Math::Clausen::S ( n , x ) ;     
}
// ============================================================================

// ============================================================================
//                                                                      The END 
// ============================================================================
