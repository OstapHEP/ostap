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
#include "Ostap/ChebyshevApproximation.h"
// ============================================================================
// Local
// ============================================================================
#include "GSL_sentry.h"
#include "local_math.h"
#include "local_gsl.h"
// ============================================================================
/** @file
 *  implementation fiel for Clausen functions 
 *  @see LHCbMath/MoreFunctions.h
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-03-26
 */
// ============================================================================
// initialize the static arrays&objects 
// ============================================================================
namespace
{
  // ==========================================================================
  const std::array<long double,2>  s_S1 { -0.50L      ,  M_PI / 2.0L  } ;
  const std::array<long double,3>  s_C2 {  0.25L      , -M_PI / 2.0L  ,  M_PI * M_PI / 6.0L  } ;
  const std::array<long double,4>  s_S3 {  1.0L /  12 , -M_PI * M_PI / 4.0L , M_PI * M_PI / 6.0L , 0.0L } ;
  const std::array<long double,5>  s_C4 {  1.0L /  48 ,  M_PI / 12.0L , -M_PI * M_PI / 12.0L , 0.0L , std::pow ( M_PI , 4 ) / 90.0L } ;
  const std::array<long double,6>  s_S5 { -1.0L / 240 ,  M_PI / 48.0L , -M_PI * M_PI / 36.0L , 0.0L , std::pow ( M_PI , 4 ) / 90.0L , 0.0L } ;
  // ==========================================================================
  const Ostap::Math::ChebyshevSum  s_c3
    { Ostap::Math::ChebyshevApproximation
      ( Ostap::Math::Clausen::C_<3>  () , 0.0 , M_PI , 70 ).polynomial ( 1.e-10 ) } ;
  // ==========================================================================
  const Ostap::Math::ChebyshevSum  s_s4
    { Ostap::Math::ChebyshevApproximation
      ( Ostap::Math::Clausen::S_<4>  () , 0.0 , M_PI , 50 ).polynomial ( 1.e-10 ) } ;
  // ==========================================================================
  const Ostap::Math::ChebyshevSum  s_c5
    { Ostap::Math::ChebyshevApproximation
      ( Ostap::Math::Clausen::C_<5>  () , 0.0 , M_PI , 25 ).polynomial ( 1.e-10 ) } ;
  // ==========================================================================
  const Ostap::Math::ChebyshevSum  s_s6
    { Ostap::Math::ChebyshevApproximation
      ( Ostap::Math::Clausen::S_<6>  () , 0.0 , M_PI , 20 ).polynomial ( 1.e-10 ) } ;
  // ==========================================================================
  const Ostap::Math::ChebyshevSum  s_c7
    { Ostap::Math::ChebyshevApproximation
      ( Ostap::Math::Clausen::C_<7>  () , 0.0 , M_PI , 15 ).polynomial ( 1.e-10 ) } ;
  // ==========================================================================
  const Ostap::Math::ChebyshevSum  s_s8
    { Ostap::Math::ChebyshevApproximation
      ( Ostap::Math::Clausen::S_<9>  () , 0.0 , M_PI , 15 ).polynomial ( 1.e-10 ) } ;
  // ==========================================================================
  const Ostap::Math::ChebyshevSum  s_c9
    { Ostap::Math::ChebyshevApproximation
      ( Ostap::Math::Clausen::C_<9>  () , 0.0 , M_PI , 15 ).polynomial ( 1.e-10 ) } ;
  // ==========================================================================
  const Ostap::Math::ChebyshevSum  s_s10
    { Ostap::Math::ChebyshevApproximation
      ( Ostap::Math::Clausen::S_<10> () , 0.0 , M_PI , 15 ).polynomial ( 1.e-10 ) } ;
  // ==========================================================================
} //                                            The END of anonymous namespace 
// ============================================================================
// S-functions 
// ============================================================================
double Ostap::Math::Clausen::S0 ( const double x )
{ return 0.5 * std::tan ( 0.5 * M_PI - x * 0.5L ) ; }
// ============================================================================
double Ostap::Math::Clausen::S1 ( const double x )
{  
  const long double y = Ostap::Math::reduce ( 1.0L * x , 0.0L , 2.0L * M_PI ).first ;
  return Ostap::Math::Clenshaw::monomial_sum ( s_S1.begin () , s_S1.end () , y ) . first ;
}
// ============================================================================
double Ostap::Math::Clausen::S2 ( const double x )
{ return Ostap::Math::clausen ( x ) ; }
// ============================================================================
double Ostap::Math::Clausen::S3 ( const double x )
{  
  const long double y = Ostap::Math::reduce ( 1.0L * x , 0.0L , 2.0L * M_PI ).first ;
  return Ostap::Math::Clenshaw::monomial_sum ( s_S3.begin () , s_S3.end () , y ) . first ;
}
// ============================================================================
double Ostap::Math::Clausen::S4 ( const double x )
{  
  const long double y = Ostap::Math::reduce ( 1.0L * x , -1.0L * M_PI , +1.0L * M_PI ).first ;
  const double z = s_s4 ( std::abs ( y ) ) ;
  return 0 < y ? z : -z ;
}
// ============================================================================
double Ostap::Math::Clausen::S5 ( const double x )
{  
  const long double y = Ostap::Math::reduce ( 1.0L * x , 0.0L , 2.0L * M_PI ).first ;
  return Ostap::Math::Clenshaw::monomial_sum ( s_S5.begin () , s_S5.end () , y ) . first ;
}
// ============================================================================
double Ostap::Math::Clausen::S6 ( const double x )
{  
  const long double y = Ostap::Math::reduce ( 1.0L * x , -1.0L * M_PI , +1.0L * M_PI ).first ;
  const double z = s_s6 ( std::abs ( y ) ) ;
  return 0 < y ? z : -z ;
}
// =============================================================================
double Ostap::Math::Clausen::S7 ( const double x )
{
  const unsigned int m = 4 ;
  static const Ostap::Math::Bernulli s_B { 7 } ;
  static const long double           s_C { std::pow ( 2.0L * M_PI , 7 ) / ( 2 * 5040 ) } ;
  //
  const long double y = Ostap::Math::reduce ( 1.0L * x , 0.0L , 2.0L * M_PI ).first ;
  return s_C * s_B ( y / ( 2.0L * M_PI ) ) ;
}
// ============================================================================
double Ostap::Math::Clausen::S8 ( const double x )
{  
  const long double y = Ostap::Math::reduce ( 1.0L * x , -1.0L * M_PI , +1.0L * M_PI ).first ;
  const double z = s_s8 ( std::abs ( y ) ) ;
  return 0 < y ? z : -z ;
}
// =============================================================================
double Ostap::Math::Clausen::S9 ( const double x )
{
  const unsigned int m = 5 ;
  static const Ostap::Math::Bernulli s_B { 9 } ;
  static const long double           s_C { std::pow ( 2.0L * M_PI , 9 ) / ( 2 * 362880L ) } ;
  //
  const long double y = Ostap::Math::reduce ( 1.0L * x , 0.0L , 2.0L * M_PI ).first ;
  return - s_C * s_B ( y / ( 2.0L * M_PI ) ) ;
}
// ============================================================================
double Ostap::Math::Clausen::S10 ( const double x )
{  
  const long double y = Ostap::Math::reduce ( 1.0L * x , -1.0L * M_PI , +1.0L * M_PI ).first ;
  const double z = s_s10 ( std::abs ( y ) ) ;
  return 0 < y ? z : -z ;
}
// ============================================================================
// C-functions 
// ============================================================================
double Ostap::Math::Clausen::C0 ( const double /* x */ ) { return 0.5 ; }
// ============================================================================
double Ostap::Math::Clausen::C1 ( const double x )
{ return - std::log ( std::abs ( 2 * std::sin ( 1.0L * x / 2 ) ) ) ; }
// ============================================================================
double Ostap::Math::Clausen::C2 ( const double x )
{  
  const long double y = Ostap::Math::reduce ( 1.0L * x , 0.0L , 2.0L * M_PI ).first ;
  return Ostap::Math::Clenshaw::monomial_sum ( s_C2.begin () , s_C2.end () , y ) . first ;
}
// ============================================================================
double Ostap::Math::Clausen::C3 ( const double x )
{  
  const long double y = Ostap::Math::reduce ( 1.0L * x , -1.0L * M_PI , +1.0L * M_PI ).first ;
  return s_c3 ( std::abs ( y ) ) ;
}
// ============================================================================
double Ostap::Math::Clausen::C4 ( const double x )
{  
  const long double y = Ostap::Math::reduce ( 1.0L * x , 0.0L , 2.0L * M_PI ).first ;
  return Ostap::Math::Clenshaw::monomial_sum ( s_C4.begin () , s_C4.end () , y ) . first ;
}
// ============================================================================
double Ostap::Math::Clausen::C5 ( const double x )
{  
  const long double y = Ostap::Math::reduce ( 1.0L * x , -1.0L * M_PI , +1.0L * M_PI ).first ;
  return s_c5 ( std::abs ( y ) ) ;
}
// ============================================================================
double Ostap::Math::Clausen::C6 ( const double x )
{
  const unsigned int m = 3 ;
  static const Ostap::Math::Bernulli s_B { 2 * m  } ;
  static const long double           s_C { std::pow ( 2.0L * M_PI , 6 ) / ( 2 * 720 ) } ;
  //
  const long double y = Ostap::Math::reduce ( 1.0L * x , 0.0L , 2.0L * M_PI ).first ;
  return s_C * s_B ( y / ( 2.0L * M_PI ) ) ;
}
// ============================================================================
double Ostap::Math::Clausen::C7 ( const double x )
{  
  const long double y = Ostap::Math::reduce ( 1.0L * x , -1.0L * M_PI , +1.0L * M_PI ).first ;
  return s_c7 ( std::abs ( y ) ) ;
}
// =============================================================================
double Ostap::Math::Clausen::C8 ( const double x )
{
  const unsigned int m = 4 ;
  static const Ostap::Math::Bernulli s_B { 2 * m  } ;
  static const long double           s_C { std::pow ( 2.0L * M_PI , 8 ) / ( 2 * 40320 ) } ;
  //
  const long double y = Ostap::Math::reduce ( 1.0L * x , 0.0L , 2.0L * M_PI ).first ;
  return - s_C * s_B ( y / ( 2.0L * M_PI ) ) ;
}
// ============================================================================
double Ostap::Math::Clausen::C9 ( const double x )
{  
  const long double y = Ostap::Math::reduce ( 1.0L * x , -1.0L * M_PI , +1.0L * M_PI ).first ;
  return s_c9 ( std::abs ( y ) ) ;
}
// ============================================================================
double Ostap::Math::Clausen::C10 ( const double x )
{
  const unsigned int m = 5 ;
  static const Ostap::Math::Bernulli s_B { 2 * m  } ;
  static const long double           s_C { std::pow ( 2.0L * M_PI , 6 ) / ( 2 * 3628800L  ) } ;
  //
  const long double y = Ostap::Math::reduce ( 1.0L * x , 0.0L , 2.0L * M_PI ).first ;
  return s_C * s_B ( y / ( 2.0L * M_PI ) ) ;
}
// ============================================================================

// ============================================================================
// \f$ \sum_{i=1} \frac{ \sin kx }{ k^n }\f$ 
// ============================================================================
double Ostap::Math::Clausen::S
( const unsigned int n ,
  const double       x )
{
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
  // use straightforward summation of Fourier series  
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
  // use straightforward summation of Fourier series  
  constexpr std::size_t N = 35 ;
  const std::array<long double,N> s_C { Ostap::Math::make_array
    ( [ n ] ( const std::size_t  k ) -> long double  
    { return 0 == k ? 0.0 : std::pow ( 1.0L / k , n ) ; } ,
      std::make_index_sequence<N>() ) } ;
  //
  return Ostap::Math::Clenshaw::cosine_sum ( s_C.begin() , s_C.end() , x ) ; 
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
