// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <functional>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Hash.h"
#include "Ostap/StatusCode.h"
#include "Ostap/Kinematics.h"
#include "Ostap/Dalitz.h"
#include "Ostap/DalitzIntegrator.h"
// ============================================================================
// local
// ============================================================================
#include "local_math.h"
#include "local_hash.h"
// ============================================================================
/** @file
 *  Implementation file for class Ostap::Kinematics::Dalitz
 *  @see Ostap::Kinematics::Dalitz0
 *  @see Ostap::Kinematics::Dalitz
 *  @date 2019-07-15 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// ============================================================================
// Dalitz 
// ============================================================================
Ostap::Kinematics::Dalitz0::Dalitz0
( const double m1 , 
  const double m2 , 
  const double m3 ) 
  : m_m1 ( s_zero ( m1 ) || s_zero ( m1 * m1 ) ? 0.0 : std::abs ( m1 ) )
  , m_m2 ( s_zero ( m2 ) || s_zero ( m2 * m2 ) ? 0.0 : std::abs ( m2 ) )
  , m_m3 ( s_zero ( m3 ) || s_zero ( m3 * m3 ) ? 0.0 : std::abs ( m3 ) )
    // precalculated quantities: s1_min/max, s2_min/max , s3_min/max, sum(s_i) & m_i^2
  , m_cache {
  // s1_min, s2_min, s3_min
  ( m_m1 + m_m2 ) * ( m_m1 + m_m2 )         , // [0]
    ( m_m2 + m_m3 ) * ( m_m2 + m_m3 )       , // [1] 
    ( m_m3 + m_m1 ) * ( m_m3 + m_m1 )       , // [2] 
    // mass-squared 
    m_m1 * m_m1                             , // [3] 
    m_m2 * m_m2                             , // [4]
    m_m3 * m_m3                             , // [5] 
    m_m1 * m_m1 + m_m2 * m_m2 + m_m3 * m_m3 , // [6]
    // sum of masses
    m_m1 + m_m2 + m_m3                      , // [7]
    // sum of masses
    std::pow  ( m_m1 + m_m2 + m_m3 , 2 )      // [8]
    }
  , m_cacheb { 
    s_zero ( m_m1 ) ,
      s_zero ( m_m2 ) ,
      s_zero ( m_m3 )
      }
  , m_tag ( Ostap::Utils::hash_combiner ( m_m1  , m_m2 , m_m3 ) ) 
{}
// ============================================================================
/** Is point \f$ (s_1,s_2)\f$ "inside" the Dalizt plot?
 *  Get the sign of G-function 
 *  \f$ g(s_1,s_2) = G ( s_1, s_2 , s , m_2^2, m_1^2, m_3^2) \f$
 *  @see Ostap::Math::PhaseSpace2::G
 *  Physical region corresponds to \f$ g\le0 \f$  
 */
// ============================================================================
bool Ostap::Kinematics::Dalitz0::inside
( const double s  ,
  const double s1 ,
  const double s2 ) const 
{
  return
    s1_min () <= s1 && s1 <= s &&   
    s2_min () <= s2 && s2 <= s &&
    sqsumm () <= s             && 
    s1 + s2   <= s + summ2 ()  && 
    0 >= Ostap::Kinematics::G ( s1 , s2 , s , m2sq() , m1sq (), m3sq () ) ;
}
// ============================================================================
/*  get the measure of the distance from the point to the boundary of Dalitz plot. 
 *  Distance is defined as \f$ d \equiv = \lambda ( P_1^2, P_2^2, P_3^2) \f$ 
 */
// ============================================================================
double Ostap::Kinematics::Dalitz0::distance
( const double s  ,
  const double s1 ,
  const double s2 ) const 
{
  //
  const double s3 = s + summ2 () - s1 - s2 ;
  //
  const double p1 = Ostap::Kinematics::triangle ( s , m1sq () , s2 ) ;
  const double p2 = Ostap::Kinematics::triangle ( s , m2sq () , s3 ) ;
  const double p3 = Ostap::Kinematics::triangle ( s , m3sq () , s1 ) ;
  //
  const double scale = 0.25 / s ;
  //
  return Ostap::Kinematics::triangle ( scale * p1 , scale * p2 , scale * p3 );
}
// ============================================================================
/*  density of Dalitz plot 
 *  \f$ \frac{d^2}{ds_1 ds_2} R_3  = \frac{\pi^2}{4s} \f$ 
 */
// ============================================================================
double Ostap::Kinematics::Dalitz0::density 
( const double s  ,
  const double s1 , 
  const double s2 ) const 
{
  static const double s_dens = 0.25 * M_PI * M_PI ;
  return
    s  <= s_min  () ? 0.0 : 
    s1 <= s1_min () ? 0.0 : 
    s2 <= s2_min () ? 0.0 : 
    !inside ( s , s1 , s2 ) ? 0.0 : s_dens / s ;
}
// ============================================================================
/*  density of Dalitz plot as function of masses 
 *  \f$ \frac{d^2}{dm_{12} dm_{23}} R_3  = \frac{\pi^2}{4s} \f$ 
 */
// ============================================================================
double Ostap::Kinematics::Dalitz0::density_mass
( const double M   ,
  const double m12 , 
  const double m23 ) const 
{
  return 
    M <= 0 || m12 <= 0 || m23 <= 0 ? 0.0 :
    density ( M * M , m12 * m12 , m23 * m23 ) * 4 * m12 * m23 ;  
}
// ============================================================================
// maximal value of momenta of the first  particle for the given s 
// ============================================================================
double Ostap::Kinematics::Dalitz0::p1_max 
( const double s ) const
{
  return s <= s_min () ? 0.0 : 
    0.5  * std::sqrt ( Ostap::Kinematics::triangle ( s , m1sq () , s2_min () ) / s ) ;
}
// ============================================================================
// maximal value of momenta of the second particle for the given s 
// ============================================================================
double Ostap::Kinematics::Dalitz0::p2_max 
( const double s ) const 
{
  return s <= s_min () ? 0.0 : 
    0.5  * std::sqrt ( Ostap::Kinematics::triangle ( s , m2sq () , s3_min () ) / s ) ;
}
// ============================================================================
// maximal value of momenta of the third  particle for the given s 
// ============================================================================
double Ostap::Kinematics::Dalitz0::p3_max 
( const double s ) const 
{
  return s <= s_min () ? 0.0 : 
    0.5  * std::sqrt ( Ostap::Kinematics::triangle ( s , m3sq () , s1_min () ) / s ) ;
}
// ============================================================================
// Dalitz plot boundaries \f$ s_1^{min/max} ( s, s_2 ) \f$ 
// ============================================================================
std::pair<double,double>  
Ostap::Kinematics::Dalitz0::s1_minmax_for_s_s2 
( const double s  , 
  const double s2 ) const 
{
  //
  static const std::pair<double,double> s_BAD { s1_min () , -s1_min () } ;
  //
  // wrong arguments ?
  if ( s < sqsumm () || s < s2 || s2 < s2_min () ) { return s_BAD ; }
  //
  const bool m12_zero = m1_zero () && m2_zero () ;
  const bool m23_zero = m2_zero () && m3_zero () ;
  const bool m31_zero = m3_zero () && m1_zero () ;
  //
  const bool all_zero  = m12_zero && m23_zero        ;
  //
  const double sqs   = std::sqrt ( s ) ;
  //
  // simple case : all masses are zero 
  if      ( all_zero )
  {
    //
    const double s_min = 0 ;
    const double s_max = s ;
    //
    return
      s2 < s_min || s2 > s_max ? s_BAD :
      s_equal ( s2 , s_min ) ? std::make_pair ( s_min , s_max ) :
      s_equal ( s2 , s_max ) ? std::make_pair ( s_min , s_min ) : std::make_pair ( 0.0 , s_max - s2 ) ; 
  }
  // two masses are zero 
  else if ( m12_zero )
  {
    //
    const double s_min = m3sq () ;
    const double s_max = s       ;
    //
    return
      s2 < s_min || s2 > s_max ? s_BAD :
      s_equal ( s2 , s_min ) ? std::make_pair ( 0.0 , 0.0 ) :
      s_equal ( s2 , s_max ) ? std::make_pair ( 0.0 , 0.0 ) : std::make_pair ( 0.0 , ( s - s2 ) * ( s2 - m3sq () ) / s2 )  ;
  }
  // two masses are zero 
  else if ( m23_zero ) 
  {
    const double s_min = 0 ;
    const double s_max = s2_max ( sqs ) ;
    //
    if      ( s2 < s_min || s2 > s_max ) { return s_BAD ; }
    else if ( s_equal ( s2 , s_min )   ) { return std::make_pair ( m1 () , s ) ; }
    else if ( s_equal ( s2 , s_max )   ) { const double q = m1 () * sqs ; return std::make_pair ( q , q ) ; }
    //
    const double a  = 1 ;
    const double b  = s2 - s - m1sq () ;
    const double c  = s * m1sq () ;
    const double d  = std::sqrt ( b * b - 4  * a * c  )  ;
    const double sc =  -0.5 / a ;
    //
    return std::make_pair ( sc * ( -b - d ) , sc * ( -b + d ) ) ;
  }
  // two masses are zero 
  else if ( m31_zero ) 
  {
    //
    const double s_min = m2sq () ;
    const double s_max = s       ;
    //
    return 
      s2 < s_min || s2 > s_max  ? s_BAD :
      s_equal ( s2 , s_min )    ? std::make_pair ( s_max , s_max ) :
      s_equal ( s2 , s_max )    ? std::make_pair ( s_min , s_min ) :
      std::make_pair ( s * s_min /  s2 , s + s_min - s2 ) ; 
  }
  //
  // generic case 
  //
  const double s_a       = m1sq()  + m2sq()  ;
  //
  const double s_min     = s2_min () ;
  const double s_max     = s2_max ( sqs ) ;
  //
  if ( s2 < s_min || s2 > s_max ) { return s_BAD  ; }
  //
  const double f1 = Ostap::Kinematics::triangle ( s2 , s       , m1sq () ) ;
  const double f2 = Ostap::Kinematics::triangle ( s2 , m2sq () , m3sq () ) ;
  //
  if ( f1 < 0 || f2 < 0         ) { return s_BAD ; }
  //
  const double b  = ( s2 - s  + m1sq() ) * ( s2 + m2sq () - m3sq () ) ;
  const double c  = std::sqrt ( f1 * f2 ) ;
  //
  const double q1 = s_a - ( b - c ) / ( 2 * s2 ) ;
  const double q2 = s_a - ( b + c ) / ( 2 * s2 ) ;
  //
  return std::make_pair ( q2 , q1 ) ;
}



// ============================================================================
// Dalitz plot boundaries \f$ s_2^{min/max} ( s , s_1 ) \f$ 
// ============================================================================
std::pair<double,double>  
Ostap::Kinematics::Dalitz0::s2_minmax_for_s_s1 
( const double s  ,
  const double s1 )
  const 
{
  //
  static const std::pair<double,double> s_BAD { s2_min () , -s2_min () } ;
  //
  // wrong arguments ?
  if ( s < sqsumm () || s < s1 || s1 < s1_min () ) { return s_BAD ; }
  //
  const bool m12_zero  = m1_zero () && m2_zero() ;
  const bool m23_zero  = m2_zero () && m3_zero() ;
  const bool m31_zero  = m3_zero () && m1_zero() ;
  //
  const bool all_zero  = m12_zero        && m23_zero        ;
  //
  const double sqs   = std::sqrt ( s ) ;
  //
  if      ( all_zero )
  {
    //
    const double s_min = 0 ;
    const double s_max = s ;
    //
    return
      s1 < s_min || s1 > s_max ? s_BAD :
      s_equal ( s1 , s_min ) ? std::make_pair ( s_min , s_max ) :
      s_equal ( s1 , s_max ) ? std::make_pair ( s_min , s_min ) : std::make_pair ( 0.0 , s_max - s1 ) ; 
  }
  //
  else if  ( m23_zero )
  {
    //
    const double s_min = m1sq () ;
    const double s_max = s       ;
    //
    return
      s1 < s_min || s1 > s_max ? s_BAD :
      s_equal ( s1 , s_min ) ? std::make_pair ( 0.0 , 0.0 ) :
      s_equal ( s1 , s_max ) ? std::make_pair ( 0.0 , 0.0 ) : std::make_pair ( 0.0 , ( s - s1 ) * ( s1 - m1sq () ) / s1 )  ;
  }
  //
  else if ( m12_zero ) 
  {
    const double s_min = 0               ;
    const double s_max = s1_max ( sqs )  ;
    //
    if      ( s1 < s_min || s1 > s_max ) { return s_BAD ; }
    else if ( s_equal ( s1 , s_min )   ) { return std::make_pair ( m1 () , s ) ; }
    else if ( s_equal ( s1 , s_max )   ) { const double q = m3 () * sqs ; return std::make_pair ( q , q ) ; }
    //
    const double a  = 1 ;
    const double b  = s1 - s - m3sq () ;
    const double c  = s * m3sq () ;
    const double d  = std::sqrt ( b * b - 4  * a * c  )  ;
    const double sc =  -0.5 / a ;
    //
    return std::make_pair ( sc * ( -b - d ) , sc * ( -b + d ) ) ;
  }
  else if ( m31_zero ) 
  {
    //
    const double s_min = m2sq () ;
    const double s_max = s       ;
    //
    return 
      s1 < s_min || s1 > s_max  ? s_BAD :
      s_equal ( s1 , s_min )    ? std::make_pair ( s_max , s_max ) :
      s_equal ( s1 , s_max )    ? std::make_pair ( s_min , s_min ) :
      std::make_pair ( s * s_min /  s1 , s + s_min - s1 ) ; 
  }
  //
  const double s_a       = m3sq () + m2sq()  ;
  //
  const double s_min     = s1_min (     ) ;
  const double s_max     = s1_max ( sqs ) ;
  //
  if ( s1 < s_min || s1 > s_max ) { return std::make_pair ( 1.0 , -1.0 ) ; }
  //
  const double f1 = Ostap::Kinematics::triangle ( s1 , s       , m3sq () ) ;
  const double f2 = Ostap::Kinematics::triangle ( s1 , m2sq () , m1sq () ) ;
  //
  if ( f1 < 0 || f2 < 0 ) { return std::make_pair ( 1.0 , -1.0 ) ; }
  //
  const double b  = ( s1 - s + m3sq () ) * ( s1 + m2sq () - m1sq()  ) ;
  const double c  = std::sqrt ( f1 * f2 ) ;
  //
  const double q1 = s_a - ( b - c ) / ( 2 * s1 ) ;
  const double q2 = s_a - ( b + c ) / ( 2 * s1 ) ;
  //
  return std::make_pair ( q2 , q1 ) ; 
}
// ============================================================================
/** the first x-variable is just \f$ x_1 = \cos_{R23}(12) \f$ 
 *  - cosine on the angle between 1st and 2nd particles in the  (2,3) rest frame
 *  \f$ \cos \theta_{12}^{R(2,3)}
 */
// ============================================================================
double Ostap::Kinematics::Dalitz0::x1
( const double s  ,
  const double s1 ,
  const double s2 ) const
{
  //
  if  ( !inside ( s , s1  , s2 ) ) { return -1000 ; }
  //
  const double f1 = Ostap::Kinematics::triangle ( s  , s2      , m1sq () ) ;
  const double f2 = Ostap::Kinematics::triangle ( s2 , m2sq () , m3sq () ) ;
  //
  if ( 0 >= f1 || 0 >= f2       ) { return  -1000 ; }
  //
  const double f = ( s - s2 - m1sq () ) * ( s2 + m2sq () - m3sq () ) 
                            + 2 * s2 * ( m1sq () + m2sq () - s1 ) ;
  //
  return f / std::sqrt ( f1 * f2 ) ;
}
// ============================================================================
/* (inverse) variable transformation  
 *   \f[ \begin{array{l} 
 *        s_1 = f_1 ( x_1 ,x_2  )  \\ 
 *        s_2 = f_2 ( x_1 ,x_2  )
 *       \end{array}\f]
 *   where 
 *   \f[ \begin{array{l} 
 *        x_1 = \cos_{R23)(12)  \\ 
 *        x_2 =s_2 
 *       \end{array} \f]
 *  @code
 *  Dalitz d = ... ;
 *  const double x1 = 0.1  ;
 *  const double x2 = 14.5 ;
 *  double s1, s2 ;
 *   std::tie ( s1, s2 ) = d.x2s ( s , x1 , x2 ) ; 
 *  @endcode  
 */
// ============================================================================
std::pair<double,double> Ostap::Kinematics::Dalitz0::x2s
( const double s  ,
  const double x1 ,
  const double x2 ) const 
{
  if ( s  < sqsumm () ) { return std::make_pair ( -1 , -1  ) ; }
  // adjust to the allowed boundaries 
  const double s2 = std::min ( std::max ( x2 , s2_min () ) , s2_max ( std::sqrt ( s ) ) ) ;
  const double ct = std::min ( std::max ( x1 , -1.0      ) , 1.0       ) ;
  //
  const double f1 = Ostap::Kinematics::triangle ( s  , s2      , m1sq () ) ;
  const double f2 = Ostap::Kinematics::triangle ( s2 , m2sq () , m3sq () ) ;
  //
  const double f  = ( s - s2 - m1sq () ) * ( s2 + m2sq () - m3sq () ) 
    + 2 * s2 * ( m1sq () + m2sq ()  ) ;
  //
  const double s1 = ( f - ct * std::sqrt ( f1 * f2 ) ) / ( 2 * s2 ) ;
  //
  return std::make_pair ( s1 , s2  ) ;
}
// ============================================================================
/** (inverse) variable transformation  
 *   \f[ \begin{array{l} 
 *        s   = f_1 ( y_1 ,y_2  )  \\ 
 *        s_1 = f_2 ( y_1 ,y_2  )
 *       \end{array}\f]
 *   where 
 *   \f[ \begin{array{l} 
 *        y_1 = s \\  
 *        y_2 = \cos_{R23)(12) 
 *       \end{array} \f]
 *  @code
 *  Dalitz d = ... ;
 *  const double y1 = 0.1  ;
 *  const double y2 = 0.5 ;
 *  double s, s1 ;
 *   std::tie ( s , s1 ) = d.y2s ( s2 , y1 , y2 ) ; 
 *  @endcode  
 */
// ============================================================================
std::pair<double,double> Ostap::Kinematics::Dalitz0::y2s
( const double s2 ,
  const double y1 ,
  const double y2 ) const
{
  const double s = std::max  ( y1 , sqsumm () ) ;
  if ( s2 < s2_min() || s2 > s2_max  ( s ) ) { return std::make_pair ( -1 , -1  ) ; }
  //
  const double ct = std::min ( std::max ( y2 , -1.0 ) , 1.0 ) ;
  const double f1 = Ostap::Kinematics::triangle ( s  , s2      , m1sq () ) ;
  const double f2 = Ostap::Kinematics::triangle ( s2 , m2sq () , m3sq () ) ;
  //
  const double f  = ( s - s2 - m1sq () ) * ( s2 + m2sq () - m3sq () ) 
    + 2 * s2 * ( m1sq () + m2sq ()  ) ;
  //
  const double s1 = ( f - ct * std::sqrt ( f1 * f2 ) ) / ( 2 * s2 ) ;
  //
  return std::make_pair ( s , s1  ) ;
}
// ============================================================================
/**  absolute value of the jacobian  
 *   \f$ J(s, s_1,s_2) = \left| \frac{\partial(s_1,s_2) }{\partial(x_1,x_2)} \right| \f$ 
 */
// ============================================================================
double Ostap::Kinematics::Dalitz0::J
( const double s  ,
  const double s1 ,
  const double s2 ) const
{
  if  ( !inside ( s , s1 , s2 ) ) { return 0 ; }
  //
  const double f1 = Ostap::Kinematics::triangle ( s  , s2      , m1sq () ) ;
  const double f2 = Ostap::Kinematics::triangle ( s2 , m2sq () , m3sq () ) ;
  //
  return f1 <= 0 || f2 <= 0? 0.0 : std::sqrt ( f1 * f2 ) / ( 2 * s2 ) ;
}
// ============================================================================
/*  "transpose it", such that \f$ s_{i1} \f$ and \f$ s_{i2}\f$ 
 *  becomes  the main  variable
 */
// ============================================================================
Ostap::Kinematics::Dalitz0
Ostap::Kinematics::Dalitz0::transpose 
( const unsigned short i1 , 
  const unsigned short i2 ) const 
{
  Ostap::Assert ( 1 <= i1 && i1 <= 3 , "Invalid i1" , 
                  "Ostap::Kinematics::Dalitz0::transpose" ) ;
  Ostap::Assert ( 1 <= i2 && i2 <= 3 , "Invalid i2" , 
                  "Ostap::Kinematics::Dalitz0::transpose" ) ;
  Ostap::Assert ( i1 != i2           , "Invalid i1/i2" , 
                  "Ostap::Kinematics::Dalitz0::transpose" ) ;
  
  if      ( 1 == i1 && 2 == i2 ) { return Dalitz0 ( m1 () , m2 () , m3 () ) ; }
  else if ( 1 == i1 && 3 == i2 ) { return Dalitz0 ( m2 () , m1 () , m3 () ) ; }
  else if ( 2 == i1 && 1 == i2 ) { return Dalitz0 ( m3 () , m2 () , m1 () ) ; }
  else if ( 2 == i1 && 3 == i2 ) { return Dalitz0 ( m2 () , m3 () , m1 () ) ; }
  else if ( 3 == i1 && 1 == i2 ) { return Dalitz0 ( m3 () , m1 () , m2 () ) ; }
  else if ( 3 == i1 && 2 == i2 ) { return Dalitz0 ( m1 () , m3 () , m2 () ) ; }
  //
  return Dalitz0 ( m1 () , m2 () , m3 () ) ;
}
// ============================================================================
// momentum of the 1st particle 
// ============================================================================
double Ostap::Kinematics::Dalitz0::P1 
( const double    s      , 
  const double /* s1  */ , 
  const double    s2     ) const 
{
  const double v = Ostap::Kinematics::triangle ( s , m1sq () , s2 ) ;
  return 0 < v ? 0.5 * std::sqrt ( v / s ) :  0.0 ;
}
// ============================================================================
// momentum of the 2nd particle 
// ============================================================================
double Ostap::Kinematics::Dalitz0::P2 
( const double s  , 
  const double s1 , 
  const double s2 ) const 
{
  const double v = Ostap::Kinematics::triangle ( s , m2sq () , s3 ( s , s1 , s2 ) ) ;
  return 0 < v ? 0.5 * std::sqrt ( v / s ) : 0.0 ;
}
// =============================================================================
/// momentum of the 3rd particle 
// =============================================================================
double Ostap::Kinematics::Dalitz0::P3 
( const double    s     , 
  const double    s1    , 
  const double /* s2 */ ) const 
{
  const double v = Ostap::Kinematics::triangle ( s  , m3sq () , s1 ) ;
  return 0 < v ? 0.5 * std::sqrt ( v /s ) : 0.0 ;
}
// ============================================================================
/*   \f$ \theta^{*}_{12}\f$ is angle between 
 *  \f$ p_1\f$  and \f$ p_2 \f$  in the rest frame:
 *   \f$ \cos \theta^{*}_{12} = 
 *   \left\frac { p_1o_2}{P_1P_2}\right|_{\vec{P}=0}\f$
 */
// ============================================================================
double Ostap::Kinematics::Dalitz0::cos_12 
( const double s  , 
  const double s1 , 
  const double s2 ) const 
{
  const double s3_ = s3 ( s , s1 , s2 ) ;
  const double f1  = Ostap::Kinematics::triangle ( s , m1sq () , s2  ) ;
  const double f2  = Ostap::Kinematics::triangle ( s , m2sq () , s3_ ) ;
  //
  if ( 0 >= f1 || 0 >= f2 ) { return  -1 ; }
  //
  const double f = ( s  + m1sq ()  -  s2 ) * ( s + m2sq () - s3_ ) 
    + 2 * s * ( m1sq()  + m2sq () - s1 ) ;
  //
  return  f / std::sqrt ( f1 * f2 ) ;
}
// ============================================================================
/*   \f$ \theta^{*}_{23}\f$ is angle between 
 *  \f$ p_2\f$  and \f$ p_3 \f$  in the rest frame:
 *   \f$ \cos \theta^{*}_{23} = 
 *   \left\frac { p_2p_3}{P_2P_3}\right|_{\vec{P}=0}\f$
 */
// ============================================================================
double Ostap::Kinematics::Dalitz0::cos_23 
( const double s  , 
  const double s1 , 
  const double s2 ) const 
{
  const double s3_ = s3 ( s , s1 , s2 ) ;
  const double f1  = Ostap::Kinematics::triangle ( s , m2sq () , s3_ ) ;
  const double f2  = Ostap::Kinematics::triangle ( s , m3sq () , s1  ) ;
  //
  if ( 0 >= f1 || 0 >= f2 ) { return  -1 ; }
  //
  const double f = ( s + m2sq () -  s3_ ) * ( s + m3sq () - s1 ) 
    + 2 * s * ( m2sq () + m3sq ()  - s2 ) ;
  //
  return  f / std::sqrt ( f1 * f2 ) ;
}
// ============================================================================
/*   \f$ \theta^{*}_{31}\f$ is angle between 
 *  \f$ p_3\f$  and \f$ p_1 \f$  in the rest frame:
 *   \f$ \cos \theta^{*}_{31} = 
 *   \left\frac { p_3p_1}{P_3P_1}\right|_{\vec{P}=0}\f$
 */
// ============================================================================
double Ostap::Kinematics::Dalitz0::cos_31 
( const double s  , 
  const double s1 , 
  const double s2 ) const 
{
  const double s3_ = s3 ( s , s1 , s2 ) ;
  const double f1  = Ostap::Kinematics::triangle ( s , m3sq () , s1  ) ;
  const double f2  = Ostap::Kinematics::triangle ( s , m1sq () , s2  ) ;
  //
  if ( 0 >= f1 || 0 >= f2 ) { return  -1 ; }
  //
  const double f = ( s + m3sq () -  s1 ) * ( s + m1sq () - s2 ) 
    + 2 * s * ( m3sq ()  + m1sq ()  - s3_ ) ;
  //
  return f / std::sqrt ( f1 * f2 ) ;
}
// ============================================================================
/*  \f$ \theta^{*}_{12}\f$ is angle between 
 *  \f$ p_1\f$  and \f$ p_2 \f$  in the rest frame:
 *   \f$ \sin^2 \theta^{*}_{12} = 
 *   -4s \frac{ G(s_1,s_2,s, m_2^2, m_1^2,m_3^2) }
 *   {\lambda(s, m_1^2, s2) \lambda(s, m_2^2, s3)} \f$
 */
// ============================================================================
double Ostap::Kinematics::Dalitz0::sin2_12 
( const double s  , 
  const double s1 , 
  const double s2 ) const 
{
  const double s3_ = s3 ( s ,  s1 , s2 ) ;
  const double f1  = Ostap::Kinematics::triangle ( s , m1sq () , s2  ) ;
  const double f2  = Ostap::Kinematics::triangle ( s , m2sq () , s3_ ) ;
  const double g   = Ostap::Kinematics::G (  s1 , s2 , s , m1sq () , m2sq () , m3sq () ) ;
  //
  return -4 * s * g / ( f1 * f2 ) ;
}
// ============================================================================
/*  \f$ \theta^{*}_{23}\f$ is angle between 
 *  \f$ p_2\f$  and \f$ p_3 \f$  in the rest frame:
 *   \f$ \sin^2 \theta^{*}_{23} = 
 *   -4s \frac{ G(s_2,s_3,s, m_3^2, m_2^2,m_1^2) }
 *   {\lambda(s, m_2^2, s3 ) \lambda(s, m_3^2, s1 ) } \f$
 */
// ============================================================================
double Ostap::Kinematics::Dalitz0::sin2_23 
( const double s  , 
  const double s1 , 
  const double s2 ) const 
{
  const double s3_ = s3 ( s , s1 , s2 ) ;
  const double f1  = Ostap::Kinematics::triangle ( s , m2sq () , s3_  ) ;
  const double f2  = Ostap::Kinematics::triangle ( s , m3sq () , s1   ) ;
  const double g   = Ostap::Kinematics::G (  s2 , s3_ , s , m2sq () , m3sq () , m1sq () ) ;
  //
  return -4 * s * g / ( f1 * f2 ) ;
}
// ============================================================================
/*   \f$ \theta^{*}_{31}\f$ is angle between 
 *  \f$ p_3\f$  and \f$ p_1 \f$  in the rest frame:
 *   \f$ \sin^2 \theta^{*}_{31} = 
 *   -4s \frac{ G(s_3,s_1,s, m_1^2, m_3^2,m_2^2) }
 *   {\lambda(s, m_3^2, s1 ) \lambda(s, m_1^2, s2 ) } \f$
 */
// ============================================================================
double Ostap::Kinematics::Dalitz0::sin2_31 
( const double s  , 
  const double s1 , 
  const double s2 ) const 
{
  const double s3_ = s3 ( s , s1 , s2 ) ;
  const double f1  = Ostap::Kinematics::triangle ( s , m3sq ()  , s1  ) ;
  const double f2  = Ostap::Kinematics::triangle ( s , m1sq () , s2  ) ;
  const double g   = Ostap::Kinematics::G (  s3_ , s1 , s , m3sq () , m1sq () , m2sq () ) ;
  //
  return -4 * s * g / ( f1 * f2 ) ;
}
// =============================================================================
// total momentum in (1,2) restframe
// =============================================================================
double Ostap::Kinematics::Dalitz0::P_R12  
( const double    s     , 
  const double    s1    , 
  const double /* s2 */ ) const
{
  const double f1 = Ostap::Kinematics::triangle ( s , s1 , m3sq () ) ;
  return 0 < f1 ? 0.5 * std::sqrt ( f1 / s1 ) : 0.0 ;
}
// =============================================================================
// momentum of 1st particle in (1,2) restframe
// =============================================================================
double Ostap::Kinematics::Dalitz0::P1_R12 
( const double /* s  */ , 
  const double    s1    , 
  const double /* s2 */ ) const
{
  const double f1 = Ostap::Kinematics::triangle ( s1 , m1sq() , m2sq () ) ;
  return 0 < f1 ? 0.5 * std::sqrt ( f1 / s1 ) : 0.0 ;
}
// ============================================================================
/*  cosine on the angle between 3rd and 1st particles in the  (1,2) rest frame
 *  \f$ \cos \theta_{21}^{R(1,2)}
 */
// ============================================================================
double Ostap::Kinematics::Dalitz0::cos_31_R12 
( const double s  , 
  const double s1 , 
  const double s2 ) const
{
  //
  const double s3_ = s3 ( s , s1 , s2 ) ;
  const double f1 = Ostap::Kinematics::triangle ( s  , s1      , m3sq () ) ;
  const double f2 = Ostap::Kinematics::triangle ( s1 , m1sq () , m2sq () ) ;
  //
  if ( 0 >= f1 || 0 >= f2 ) { return  -1 ; }
  //
  const double f  = ( s - s1 - m3sq()  ) * ( s1 + m1sq() - m2sq ()  ) 
    + 2 * s1 * ( m3sq()  + m1sq()  - s3_ ) ;
  //
  return f / std::sqrt ( f1 * f2 ) ;
}
// ============================================================================
/*  sine squared  on the angle between 3rd and 1st particles in the  (1,2) rest frame
 *  \f$ \sin^2 \theta_{12}^{R(2,3)}
 */
// ============================================================================
double Ostap::Kinematics::Dalitz0::sin2_31_R12 
( const double s  , 
  const double s1 , 
  const double s2 ) const 
{
  const double s3_ = s3 ( s , s1 , s2 ) ;
  const double f1 = Ostap::Kinematics::triangle ( s  , s1      , m3sq () ) ;
  const double f2 = Ostap::Kinematics::triangle ( s1 , m1sq () , m2sq () ) ;  
  const double g  = Ostap::Kinematics::G ( s3_ , s1 , s , m1sq () , m3sq () , m2sq () ) ;
  //
  return -4 * s1 * g / ( f1 * f2 ) ;
}
// =============================================================================
// total momentum in (2,3) restframe
// =============================================================================
double Ostap::Kinematics::Dalitz0::P_R23  
( const double    s     , 
  const double /* s1 */ , 
  const double    s2 ) const
{
  const double f1 = Ostap::Kinematics::triangle ( s , s2 , m1sq () ) ;
  return 0 < f1 ? 0.5 * std::sqrt ( f1 / s2 ) : 0.0 ;
}
// =============================================================================
// momentum of 2nd particle in (2,3) restframe
// =============================================================================
double Ostap::Kinematics::Dalitz0::P2_R23  
( const double /* s  */ , 
  const double /* s1 */ , 
  const double    s2    ) const
{
  const double f1 = Ostap::Kinematics::triangle ( s2 , m2sq () , m3sq () ) ;
  return 0 < f1 ? 0.5 * std::sqrt ( f1 / s2 ) : 0.0 ;
}
// ============================================================================
/*  cosine on the angle between 1st and 2nd particles in the  (2,3) rest frame
 *  \f$ \cos \theta_{12}^{R(2,3)}
 */
// ============================================================================
double Ostap::Kinematics::Dalitz0::cos_12_R23 
( const double s  , 
  const double s1 , 
  const double s2 ) const
{
  //
  const double f1 = Ostap::Kinematics::triangle ( s  , s2      , m1sq () ) ;
  const double f2 = Ostap::Kinematics::triangle ( s2 , m2sq () , m3sq () ) ;
  //
  if ( 0 >= f1 || 0 >= f2 ) { return  -1 ; }
  //
  const double f  = ( s - s2 - m1sq () ) * ( s2 + m2sq () - m3sq () ) 
    + 2 * s2 * ( m1sq () + m2sq () - s1 ) ;
  //
  return f / std::sqrt ( f1 * f2 ) ;
}
// ============================================================================
/*  sine squared  on the angle between 1st and 2nd particles in the  (2,3) rest frame
 *  \f$ \sin^2 \theta_{12}^{R(2,3)}
 */
// ============================================================================
double Ostap::Kinematics::Dalitz0::sin2_12_R23 
( const double s  , 
  const double s1 , 
  const double s2 ) const 
{
  const double f1 = Ostap::Kinematics::triangle ( s  , s2      , m1sq ()  ) ;
  const double f2 = Ostap::Kinematics::triangle ( s2 , m2sq () , m3sq ()  ) ;  
  const double g  = Ostap::Kinematics::G ( s1 , s2 , s , m2sq () , m1sq () , m3sq () ) ;
  //
  return -4 * s2 * g / ( f1 * f2 ) ;
}
// =============================================================================
// total momentum in (3,1) restframe
// =============================================================================
double Ostap::Kinematics::Dalitz0::P_R31 
( const double s  , 
  const double s1 , 
  const double s2 ) const
{
  const double s3_ = s3 ( s , s1 , s2 ) ;
  const double f1 = Ostap::Kinematics::triangle ( s , s3_ , m2sq () ) ;
  return 0 < f1 ? 0.5 * std::sqrt ( f1 / s3_ ) : 0.0 ;
}
// =============================================================================
// momentum of 3rd particle in (3,1) restframe
// =============================================================================
double Ostap::Kinematics::Dalitz0::P3_R31  
(const double s  , 
 const double s1 , 
 const double s2 ) const
{
  const double s3_ = s3 ( s , s1 , s2 ) ;
  const double f1 = Ostap::Kinematics::triangle ( s3_ , m3sq () , m1sq () ) ;
  return 0 < f1 ? 0.5 * std::sqrt ( f1 / s3_ ) : 0.0 ;
}
// ============================================================================
/* cosine on the angle between 2nd and 3rd particles in the  (3,1) rest frame
 *  \f$ \cos \theta_{23}^{R(3,1)}
 */
// ============================================================================
double Ostap::Kinematics::Dalitz0::cos_23_R31 
( const double s  , 
  const double s1 , 
  const double s2 ) const 
{
  const double s3_ = s3 ( s , s1 , s2 ) ;
  //
  const double f1 = Ostap::Kinematics::triangle ( s   , s3_     , m2sq () ) ;
  const double f2 = Ostap::Kinematics::triangle ( s3_ , m3sq () , m1sq () ) ;
  //
  if ( 0 >= f1 || 0 >= f2 ) { return  -1 ; }
  //
  const double f  = ( s - s3_ - m2sq () ) * ( s3_ + m3sq () - m1sq () ) 
    + 2 * s3_ * ( m2sq () + m3sq() - s2 ) ;
  //
  return f / std::sqrt ( f1 * f2 ) ;
}
// ============================================================================
/*  sine squared  on the angle between 2nd and 3rd particles in the  (3,1) rest frame
 *  \f$ \sin^2 \theta_{23}^{R(3,1)}
 */
// ============================================================================
double Ostap::Kinematics::Dalitz0::sin2_23_R31 
( const double s  , 
  const double s1 , 
  const double s2 ) const 
{
  const double s3_ = s3 ( s , s1 , s2 ) ;
  //
  const double f1 = Ostap::Kinematics::triangle ( s    , s3_     , m2sq () ) ;
  const double f2 = Ostap::Kinematics::triangle ( s3_  , m3sq () , m1sq () ) ;  
  const double g  = Ostap::Kinematics::G ( s2 , s3_ , s , m3sq () , m2sq () , m1sq () ) ;
  //
  return -4 * s3_ * g / ( f1 * f2 ) ;
}
// ============================================================================
/*  \f$ \cos \zeta_{1(3)}^{1} \f$ 
 *   @see M.Mikhasenko et al, "Dalitz plot decomposiiton for three-body decays", 
 *                            Phys. Rev. D 101, 034033 (2020)
 *   @see https://doi.org/10.48550/arXiv.1910.04566
 *   @see https://arxiv.org/abs/1910.04566
 */
// ============================================================================
double Ostap::Kinematics::Dalitz0::cos_zeta131 
( const double s  , 
  const double s1 , 
  const double s2 ) const 
{
  //
  const double sig1 = sigma1 ( s , s1 , s2 ) ;
  const double sig2 = sigma2 ( s , s1 , s2 ) ;
  const double sig3 = sigma3 ( s , s1 , s2 ) ;
  //
  const double f1   = Ostap::Kinematics::triangle ( s    , m1sq () , sig1    ) ;
  const double f2   = Ostap::Kinematics::triangle ( sig3 , m1sq () , m2sq () ) ;
  //
  if ( 0 >= f1 || 0 >= f2 ) { return -1 ; }
  //
  return 
    ( 2 * m1sq() * ( sig2 - s - m2sq () ) + ( s + m1sq() - sig1 ) * ( sig3 - m1sq() - m2sq() ) )
    / std::sqrt ( f1 * f2 ) ;

}
// ============================================================================
/*  \f$ \cos \zeta_{2(1)}^{1} \f$ 
 *   @see M.Mikhasenko et al, "Dalitz plot decomposiiton for three-body decays", 
 *                            Phys. Rev. D 101, 034033 (2020)
 *   @see https://doi.org/10.48550/arXiv.1910.04566
 *   @see https://arxiv.org/abs/1910.04566
 */
// ============================================================================
double Ostap::Kinematics::Dalitz0::cos_zeta211
( const double s  , 
  const double s1 , 
  const double s2 ) const 
{
  //
  const double sig1 = sigma1 ( s , s1 , s2 ) ;
  const double sig2 = sigma2 ( s , s1 , s2 ) ;
  const double sig3 = sigma3 ( s , s1 , s2 ) ;
  //
  const double f1   = Ostap::Kinematics::triangle ( s    , m1sq () , sig1    ) ;
  const double f2   = Ostap::Kinematics::triangle ( sig2 , m1sq () , m3sq () ) ;
  //
  if ( 0 >= f1 || 0 >= f2 ) { return -1 ; }
  //
  return 
    ( 2 * m1sq() * ( sig3 - s - m3sq () ) + ( s + m1sq() - sig1 ) * ( sig2 - m1sq() - m3sq() ) )
    / std::sqrt ( f1 * f2 ) ;
}
// ============================================================================
/*  \f$ \cos \zeta_{2(1)}^{2} \f$ 
 *   @see M.Mikhasenko et al, "Dalitz plot decomposiiton for three-body decays", 
 *                            Phys. Rev. D 101, 034033 (2020)
 *   @see https://doi.org/10.48550/arXiv.1910.04566
 *   @see https://arxiv.org/abs/1910.04566
 */
// ============================================================================
double Ostap::Kinematics::Dalitz0::cos_zeta212
( const double s  , 
  const double s1 , 
  const double s2 ) const 
{
  //
  const double sig1 = sigma1 ( s , s1 , s2 ) ;
  const double sig2 = sigma2 ( s , s1 , s2 ) ;
  const double sig3 = sigma3 ( s , s1 , s2 ) ;
  //
  const double f1   = Ostap::Kinematics::triangle ( s    , m2sq () , sig2    ) ;
  const double f2   = Ostap::Kinematics::triangle ( sig1 , m2sq () , m3sq () ) ;
  //
  if ( 0 >= f1 || 0 >= f2 ) { return -1 ; }
  //
  return 
    ( 2 * m2sq() * ( sig3 - s - m3sq () ) + ( s + m2sq() - sig2 ) * ( sig1 - m2sq() - m3sq() ) )
    / std::sqrt ( f1 * f2 ) ;
}
// ============================================================================
/*  \f$ \cos \zeta_{3(2)}^{2} \f$ 
 *   @see M.Mikhasenko et al, "Dalitz plot decomposiiton for three-body decays", 
 *                            Phys. Rev. D 101, 034033 (2020)
 *   @see https://doi.org/10.48550/arXiv.1910.04566
 *   @see https://arxiv.org/abs/1910.04566
 */
// ============================================================================
double Ostap::Kinematics::Dalitz0::cos_zeta322
( const double s  , 
  const double s1 , 
  const double s2 ) const 
{
  //
  const double sig1 = sigma1 ( s , s1 , s2 ) ;
  const double sig2 = sigma2 ( s , s1 , s2 ) ;
  const double sig3 = sigma3 ( s , s1 , s2 ) ;
  //
  const double f1   = Ostap::Kinematics::triangle ( s    , m2sq () , sig2    ) ;
  const double f2   = Ostap::Kinematics::triangle ( sig3 , m2sq () , m1sq () ) ;
  //
  if ( 0 >= f1 || 0 >= f2 ) { return -1 ; }
  //
  return 
    ( 2 * m2sq() * ( sig1 - s - m1sq () ) + ( s + m2sq() - sig2 ) * ( sig3 - m2sq() - m1sq() ) )
    / std::sqrt ( f1 * f2 ) ;
}
// ============================================================================
/*  \f$ \cos \zeta_{3(2)}^{3} \f$ 
 *   @see M.Mikhasenko et al, "Dalitz plot decomposiiton for three-body decays", 
 *                            Phys. Rev. D 101, 034033 (2020)
 *   @see https://doi.org/10.48550/arXiv.1910.04566
 *   @see https://arxiv.org/abs/1910.04566
 */
// ============================================================================
double Ostap::Kinematics::Dalitz0::cos_zeta323
( const double s  , 
  const double s1 , 
  const double s2 ) const 
{
  //
  const double sig1 = sigma1 ( s , s1 , s2 ) ;
  const double sig2 = sigma2 ( s , s1 , s2 ) ;
  const double sig3 = sigma3 ( s , s1 , s2 ) ;
  //
  const double f1   = Ostap::Kinematics::triangle ( s    , m3sq () , sig3    ) ;
  const double f2   = Ostap::Kinematics::triangle ( sig2 , m3sq () , m1sq () ) ;
  //
  if ( 0 >= f1 || 0 >= f2 ) { return -1 ; }
  //
  return 
    ( 2 * m3sq() * ( sig1 - s - m1sq () ) + ( s + m3sq() - sig3 ) * ( sig2 - m3sq() - m1sq() ) )
    / std::sqrt ( f1 * f2 ) ;
}
// ============================================================================
/*  \f$ \cos \zeta_{1(3)}^{3} \f$ 
 *   @see M.Mikhasenko et al, "Dalitz plot decomposiiton for three-body decays", 
 *                            Phys. Rev. D 101, 034033 (2020)
 *   @see https://doi.org/10.48550/arXiv.1910.04566
 *   @see https://arxiv.org/abs/1910.04566
 */
// ============================================================================
double Ostap::Kinematics::Dalitz0::cos_zeta133
( const double s  , 
  const double s1 , 
  const double s2 ) const 
{
  //
  const double sig1 = sigma1 ( s , s1 , s2 ) ;
  const double sig2 = sigma2 ( s , s1 , s2 ) ;
  const double sig3 = sigma3 ( s , s1 , s2 ) ;
  //
  const double f1   = Ostap::Kinematics::triangle ( s    , m3sq () , sig3    ) ;
  const double f2   = Ostap::Kinematics::triangle ( sig1 , m3sq () , m2sq () ) ;
  //
  if ( 0 >= f1 || 0 >= f2 ) { return -1 ; }
  //
  return 
    ( 2 * m3sq() * ( sig2 - s - m2sq () ) + ( s + m3sq() - sig3 ) * ( sig1 - m3sq() - m2sq() ) )
    / std::sqrt ( f1 * f2 ) ;
}
// ============================================================================
/*  \f$ \cos \zeta_{2(3)}^{1} \f$ 
 *   @see M.Mikhasenko et al, "Dalitz plot decomposiiton for three-body decays", 
 *                            Phys. Rev. D 101, 034033 (2020)
 *   @see https://doi.org/10.48550/arXiv.1910.04566
 *   @see https://arxiv.org/abs/1910.04566
 */
// ============================================================================
double Ostap::Kinematics::Dalitz0::cos_zeta231
( const double s  , 
  const double s1 , 
  const double s2 ) const 
{
  //
  const double sig1 = sigma1 ( s , s1 , s2 ) ;
  const double sig2 = sigma2 ( s , s1 , s2 ) ;
  const double sig3 = sigma3 ( s , s1 , s2 ) ;
  //
  const double f1   = Ostap::Kinematics::triangle ( sig2 , m3sq () , m1sq () ) ;
  const double f2   = Ostap::Kinematics::triangle ( sig3 , m1sq () , m2sq () ) ;
  //
  if ( 0 >= f1 || 0 >= f2 ) { return -1 ; }
  //
  return 
    ( 2 * m1sq() * ( m2sq () + m3sq() - sig1 ) + ( sig2 - m1sq() - m3sq() ) * ( sig3 - m1sq() - m2sq() ) )
    / std::sqrt ( f1 * f2 ) ;
}
// ============================================================================
/*  \f$ \cos \zeta_{3(1)}^{2} \f$ 
 *   @see M.Mikhasenko et al, "Dalitz plot decomposiiton for three-body decays", 
 *                            Phys. Rev. D 101, 034033 (2020)
 *   @see https://doi.org/10.48550/arXiv.1910.04566
 *   @see https://arxiv.org/abs/1910.04566
 */
// ============================================================================
double Ostap::Kinematics::Dalitz0::cos_zeta312
( const double s  , 
  const double s1 , 
  const double s2 ) const 
{
  //
  const double sig1 = sigma1 ( s , s1 , s2 ) ;
  const double sig2 = sigma2 ( s , s1 , s2 ) ;
  const double sig3 = sigma3 ( s , s1 , s2 ) ;
  //
  const double f1   = Ostap::Kinematics::triangle ( sig3 , m1sq () , m2sq () ) ;
  const double f2   = Ostap::Kinematics::triangle ( sig1 , m2sq () , m3sq () ) ;
  //
  if ( 0 >= f1 || 0 >= f2 ) { return -1 ; }
  //
  return 
    ( 2 * m2sq() * ( m3sq () + m1sq() - sig2 ) + ( sig3 - m2sq() - m1sq() ) * ( sig1 - m2sq() - m3sq() ) )
    / std::sqrt ( f1 * f2 ) ;
}
// ============================================================================
/*  \f$ \cos \zeta_{1(2)}^{3} \f$ 
 *   @see M.Mikhasenko et al, "Dalitz plot decomposiiton for three-body decays", 
 *                            Phys. Rev. D 101, 034033 (2020)
 *   @see https://doi.org/10.48550/arXiv.1910.04566
 *   @see https://arxiv.org/abs/1910.04566
 */
// ============================================================================
double Ostap::Kinematics::Dalitz0::cos_zeta123
( const double s  , 
  const double s1 , 
  const double s2 ) const 
{
  //
  const double sig1 = sigma1 ( s , s1 , s2 ) ;
  const double sig2 = sigma2 ( s , s1 , s2 ) ;
  const double sig3 = sigma3 ( s , s1 , s2 ) ;
  //
  const double f1   = Ostap::Kinematics::triangle ( sig1 , m2sq () , m3sq () ) ;
  const double f2   = Ostap::Kinematics::triangle ( sig2 , m3sq () , m1sq () ) ;
  //
  if ( 0 >= f1 || 0 >= f2 ) { return -1 ; }
  //
  return 
    ( 2 * m3sq() * ( m1sq () + m2sq() - sig3 ) + ( sig1 - m3sq() - m2sq() ) * ( sig2 - m3sq() - m1sq() ) )
    / std::sqrt ( f1 * f2 ) ;
}
// ============================================================================



// ============================================================================
// Dalitz 
// ============================================================================
Ostap::Kinematics::Dalitz::Dalitz
( const double M  ,
  const Ostap::Kinematics::Dalitz::Dalitz0& b )   
  : Dalitz0 ( b )  
  , m_M     ( std::abs ( M  ) )
    // precalculated quantities: s1_min/max, s2_min/max , s3_min/max, sum(s_i) & m_i^2
  , m_cache2 { Dalitz0::s1_max ( m_M ) , // [0]
               Dalitz0::s2_max ( m_M ) , // [1]
               Dalitz0::s3_max ( m_M ) ,   // [2] 
    // sum of all invariants 
    m_M * m_M + summ2 () ,                                                          // [3] 
    // mass-squared               
    m_M  * m_M           ,                                                          // [4]
    // max e1 , e2 , e3 
    ( m_M * m_M + m1sq () - ( m2 () + m3 () ) * ( m2 () + m3 () ) ) / ( 2 * m_M ) , // [5] 
    ( m_M * m_M + m2sq () - ( m1 () + m3 () ) * ( m1 () + m3 () ) ) / ( 2 * m_M ) , // [6] 
    ( m_M * m_M + m3sq () - ( m1 () + m2 () ) * ( m1 () + m2 () ) ) / ( 2 * m_M )   // [7]
    }
                 , m_tag2 ( Ostap::Utils::hash_combiner ( Dalitz0::tag()  , m_M ) ) 
                 {
                   Ostap::Assert ( m_M > m1 ()  + m2 () + m3 ()  , 
                                   "Invalid masses for Dalitz" , 
                                   "Ostap::Kinematics::Dalitz" ) ;
}
// ============================================================================
/** Is point \f$ (s_1,s_2)\f$ "inside" the Dalizt plot?
 *  Get the sign of G-function 
 *  \f$ g(s_1,s_2) = G ( s_1, s_2 , s , m_2^2, m_1^2, m_3^2) \f$
 *  @see Ostap::Math::PhaseSpace2::G
 *  Physical region correspoinds to \f$ g\le0 \f$  
 */
// ============================================================================
bool Ostap::Kinematics::Dalitz::inside 
( const double s1 , const double s2 ) const 
{
  if ( s1  < s1_min () || s1  > s1_max () ) { return false ; }
  if ( s2  < s2_min () || s2  > s2_max () ) { return false ; }
  const double s3_ = s3 (  s1  , s2 ) ;
  if ( s3_ < s3_min () || s3_ > s3_max () ) { return false ; }
  //
  return 
    0 >= Ostap::Kinematics::G ( s1 , s2 , s () , m2sq() , m1sq (), m3sq () ) ;
}
// ============================================================================
/** Dalitz density in 1-dimension:
 *  \f$  \frac{d R_3}{d_s2} = 
 *  \frac{\pi^2}{2ss_2}
 *  \lambda^{1/2} ( s_2, s, m_1^2) 
 *  \lambda^{1/2} ( s_2, m_2^2, m_3^2) 
 *  \f$ 
 */
// ============================================================================
double Ostap::Kinematics::Dalitz::dRds2   ( const double s2 ) const 
{
  //
  if ( s2 < s2_min () || s2 > s2_max () ){ return 0 ; }
  //
  const double f1 = Ostap::Kinematics::triangle ( s2 , s ()    , m1sq () ) ;
  if ( f1 < 0 ) {  return 0 ; }
  const double f2 = Ostap::Kinematics::triangle ( s2 , m2sq () , m3sq () ) ;  
  if ( f1 < 0 ) {  return 0 ; }
  //
  static const double s_norm = 0.25 * M_PI * M_PI ;
  //
  return 0 < f1 && 0 < f2 ? s_norm * std::sqrt ( f1 * f2 ) / ( s () * s2 ) : 0.0 ; 
}
// ============================================================================
/*  Dalitz density in 1-dimension:
 *  \f$  \frac{\deriv R_3}{\deriv s3} = 
 *  \frac{\pi^2}{2ss_3}
 *  \lambda^{1/2} ( s_3, s, m_2^2) 
 *  \lambda^{1/2} ( s_3, m_3^2, m_1^2) 
 *  \f$ 
 */
// ============================================================================
double Ostap::Kinematics::Dalitz::dRds3   ( const double s3 ) const 
{
  //
  if ( s3 < s3_min () || s3 > s3_max () ){ return 0 ; }
  //
  const double f1 = Ostap::Kinematics::triangle ( s3 , s ()   , m2sq () ) ;
  if ( f1 < 0 ) {  return 0 ; }
  const double f2 = Ostap::Kinematics::triangle ( s3 , m3sq() , m1sq () ) ;  
  if ( f1 < 0 ) {  return 0 ; }
  //
  static const double s_norm = 0.25 * M_PI * M_PI ;
  //
  return 0 < f1 && 0 < f2 ? s_norm * std::sqrt ( f1 * f2 ) / ( s () * s3 ) : 0.0 ; 
}
// ============================================================================
/*  Dalitz density in 1-dimension:
 *  \f$  \frac{\deriv R_3}{\deriv s1} = 
 *  \frac{\pi^2}{2ss_1}
 *  \lambda^{1/2} ( s_1, s, m_3^2) 
 *  \lambda^{1/2} ( s_1, m_1^2, m_2^2) 
 *  \f$ 
 */
// ============================================================================
double Ostap::Kinematics::Dalitz::dRds1   ( const double s1 ) const 
{
  //
  if ( s1 < s1_min () || s1 > s1_max () ){ return 0 ; }
  //
  const double f1 = Ostap::Kinematics::triangle ( s1 , s ()   , m3sq () ) ;
  if ( f1 < 0 ) {  return 0 ; }
  const double f2 = Ostap::Kinematics::triangle ( s1 , m1sq() , m2sq () ) ;  
  if ( f1 < 0 ) {  return 0 ; }
  //
  static const double s_norm = 0.25 * M_PI * M_PI ;
  //
  return 0 < f1 && 0 < f2 ? s_norm * std::sqrt ( f1 * f2 ) / ( s () * s1 ) : 0.0 ; 
}
// ============================================================================
/*  density of Dalitz plot 
 *  \f$ \frac{d^2}{ds_1 ds_2} R_3  = \frac{\pi^2}{4s} \f$ 
 */
// ============================================================================
double Ostap::Kinematics::Dalitz::density 
( const double s1 , 
  const double s2 ) const 
{
  static const double s_dens = 0.25 * M_PI * M_PI ;
  return
    s1 <= s1_min ()     ? 0.0 : 
    s2 <= s2_min ()     ? 0.0 : 
    !inside ( s1 , s2 ) ? 0.0 : s_dens / s () ;
}
// ============================================================================
/*  density of Dalitz plot as function of masses 
 *  \f$ \frac{d^2}{dm_{12} dm_{23}} R_3  = \frac{\pi^2}{4s} \f$ 
 */
// ============================================================================
double Ostap::Kinematics::Dalitz::density_mass
( const double m12 , 
  const double m23 ) const 
{
  return 
    m12 <= 0 || m23 <= 0 ? 0.0 :
    density ( m12 * m12 , m23 * m23 ) * 4 * m12 * m23 ;  
}
// ============================================================================
//                                                                      The END 
// ============================================================================
