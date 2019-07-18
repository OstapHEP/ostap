// ============================================================================
// Include files
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Kinematics.h"
#include "Ostap/Dalitz.h"
// ============================================================================
// local
// ============================================================================
#include "Exception.h"
#include "local_math.h"
// ============================================================================
/** @file
 *  Implementation file for class Ostap::Kinematics::Dalitz
 *  @see Ostap::Kinematics::Dalitz
 *  @date 2019-07-15 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// ============================================================================
// Dalitz 
// ============================================================================
Ostap::Kinematics::Dalitz::Dalitz
( const double M  , 
  const double m1 , 
  const double m2 , 
  const double m3 ) 
  : m_M  ( std::abs ( M  ) )
  , m_m1 ( std::abs ( m1 ) )
  , m_m2 ( std::abs ( m2 ) )
  , m_m3 ( std::abs ( m3 ) )
    // precalculated quantities: s1_min/max, s2_min/max , s3_min/max, sum(s_i) & m_i^2
  , m_cache {{  
      ( m_m1 + m_m2 ) * ( m_m1 + m_m2 ) , // [0]
      ( m_M  - m_m3 ) * ( m_M  - m_m3 ) , // [1]
      ( m_m2 + m_m3 ) * ( m_m2 + m_m3 ) , // [2] 
      ( m_M  - m_m1 ) * ( m_M  - m_m1 ) , // [3]
      ( m_m3 + m_m1 ) * ( m_m3 + m_m1 ) , // [4] 
      ( m_M  - m_m2 ) * ( m_M  - m_m2 ) , // [5] 
      // sum of all invariants 
      m_M * m_M + m_m1 * m_m1 + m_m2 * m_m2 + m_m3 * m_m3  , //    [6] 
      // mass-squared 
      m_m1 * m_m1 ,  // [7] 
      m_m2 * m_m2 ,  // [8]
      m_m3 * m_m3 ,  // [9] 
      m_M  * m_M     // [10]
      }}
{
  Ostap::Assert ( m_M > m_m1 + m_m2 + m_m3 , 
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
/*  get the meaure of the distance from the point to the boundary of Dalitz plot. 
 *  Distance is defined as \f$ d \equiv = \lambda ( P_1^2, P_2^2, P_3^2) \f$ 
 */
// ============================================================================
double Ostap::Kinematics::Dalitz::distance ( const double s1 , const double s2 ) const 
{
  const double s3_ = s3 ( s1 , s2 ) ;
  //
  const double p1 = Ostap::Kinematics::triangle ( s () , m1sq () , s2  ) ;
  const double p2 = Ostap::Kinematics::triangle ( s () , m2sq () , s3_ ) ;
  const double p3 = Ostap::Kinematics::triangle ( s () , m3sq () , s1  ) ;
  //
  const double scale = 0.25 / s() ;
  //
  return Ostap::Kinematics::triangle ( scale * p1 , scale * p2 , scale * p3 );
} 
// ============================================================================
// momentum of the 1st particle 
// ============================================================================
double Ostap::Kinematics::Dalitz::P1 
( const double /* s1  */ , const double s2 ) const 
{
  const double v = Ostap::Kinematics::triangle ( s () , m1sq () , s2 ) ;
  return 0 <= v ? std::sqrt ( v ) / ( 2 * sqs () ) : -std::sqrt ( -v  ) / ( 2 * sqs () ) ;
}
// ============================================================================
// momentum of the 2nd particle 
// ============================================================================
double Ostap::Kinematics::Dalitz::P2 
( const double s1  , const double s2 ) const 
{
  const double v = Ostap::Kinematics::triangle ( s () , m2sq () , s3 ( s1 , s2 ) ) ;
  return 0 <= v ? std::sqrt ( v ) / ( 2 * sqs () ) : -std::sqrt ( -v ) / ( 2 * sqs () ) ;
}
// =============================================================================
/// momentum of the 3rd particle 
// =============================================================================
double Ostap::Kinematics::Dalitz::P3 
( const double s1  , const double /* s2 */ ) const 
{
  const double v = Ostap::Kinematics::triangle ( s () , m3sq () , s1 ) ;
  return 0 <= v ? std::sqrt ( v ) / ( 2 * sqs () ) : -std::sqrt ( -v ) / ( 2 * sqs () ) ;
} 
// ============================================================================
/*   \f$ \theta^{*}_{12}\f$ is angle between 
 *  \f$ p_1\f$  and \f$ p_2 \f$  in the rest frame:
 *   \f$ \cos \theta^{*}_{12} = 
 *   \left\frac { p_1o_2}{P_1P_2}\right|_{\vec{P}=0}\f$
 */
// ============================================================================
double Ostap::Kinematics::Dalitz::cos_12 ( const double s1 , const double s2 ) const 
{
  const double s3_ = s3 ( s1 , s2 ) ;
  //
  const double f1  = Ostap::Kinematics::triangle ( s() ,  m1sq () , s2  ) ;
  const double f2  = Ostap::Kinematics::triangle ( s() ,  m2sq () , s3_ ) ;
  //
  if ( 0 >= f1 || 0 >= f2 ) { return  -1 ; }
  //
  const double f = ( s () + m1sq ()  -  s2 ) * ( s () + m2sq ()  - s3_ ) 
    + 2 * s () * ( m1sq()  + m2sq () - s1 ) ;
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
double Ostap::Kinematics::Dalitz::cos_23 ( const double s1 , const double s2 ) const 
{
  const double s3_ = s3 ( s1 , s2 ) ;
  const double f1  = Ostap::Kinematics::triangle ( s() ,  m2sq () , s3_  ) ;
  const double f2  = Ostap::Kinematics::triangle ( s() ,  m3sq () , s1   ) ;
  //
  if ( 0 >= f1 || 0 >= f2 ) { return  -1 ; }
  //
  const double f = ( s () + m2sq () -  s3_ ) * ( s () + m3sq () - s1 ) 
    + 2 * s () * ( m2sq () + m3sq ()  - s2 ) ;
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
double Ostap::Kinematics::Dalitz::cos_31 ( const double s1 , const double s2 ) const 
{
  const double s3_ = s3 ( s1 , s2 ) ;
  const double f1  = Ostap::Kinematics::triangle ( s() ,  m3sq () , s1  ) ;
  const double f2  = Ostap::Kinematics::triangle ( s() ,  m1sq () , s2  ) ;
  //
  if ( 0 >= f1 || 0 >= f2 ) { return  -1 ; }
  //
  const double f = ( s () + m3sq () -  s1 ) * ( s () + m1sq () - s2 ) 
    + 2 * s () * ( m3sq ()  + m1sq ()  - s3_ ) ;
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
double Ostap::Kinematics::Dalitz::sin2_12 ( const double s1 , const double s2 ) const 
{
  const double s3_ = s3 ( s1 , s2 ) ;
  //
  const double f1 = Ostap::Kinematics::triangle ( s (), m1sq () , s2  ) ;
  const double f2 = Ostap::Kinematics::triangle ( s (), m2sq () , s3_ ) ;
  const double g  = Ostap::Kinematics::G (  s1 , s2 , s () , m1sq () , m2sq () , m3sq () ) ;
  //
  return -4 * s() * g / ( f1 * f2 ) ;
}
// ============================================================================
/*  \f$ \theta^{*}_{23}\f$ is angle between 
 *  \f$ p_2\f$  and \f$ p_3 \f$  in the rest frame:
 *   \f$ \sin^2 \theta^{*}_{23} = 
 *   -4s \frac{ G(s_2,s_3,s, m_3^2, m_2^2,m_1^2) }
 *   {\lambda(s, m_2^2, s3 ) \lambda(s, m_3^2, s1 ) } \f$
 */
// ============================================================================
double Ostap::Kinematics::Dalitz::sin2_23 ( const double s1 , const double s2 ) const 
{
  const double s3_ = s3 ( s1 , s2 ) ;
  //
  const double f1 = Ostap::Kinematics::triangle ( s (), m2sq () , s3_  ) ;
  const double f2 = Ostap::Kinematics::triangle ( s (), m3sq () , s1   ) ;
  const double g  = Ostap::Kinematics::G (  s2 , s3_ , s () , m2sq () , m3sq () , m1sq () ) ;
  //
  return -4 * s() * g / ( f1 * f2 ) ;
}
// ============================================================================
/*   \f$ \theta^{*}_{31}\f$ is angle between 
 *  \f$ p_3\f$  and \f$ p_1 \f$  in the rest frame:
 *   \f$ \sin^2 \theta^{*}_{31} = 
 *   -4s \frac{ G(s_3,s_1,s, m_1^2, m_3^2,m_2^2) }
 *   {\lambda(s, m_3^2, s1 ) \lambda(s, m_1^2, s2 ) } \f$
 */
// ============================================================================
double Ostap::Kinematics::Dalitz::sin2_31 ( const double s1 , const double s2 ) const 
{
  const double s3_ = s3 ( s1 , s2 ) ;
  //
  const double f1 = Ostap::Kinematics::triangle ( s (), m3sq ()  , s1  ) ;
  const double f2 = Ostap::Kinematics::triangle ( s (), m1sq () , s2  ) ;
  const double g  = Ostap::Kinematics::G (  s3_ , s1 , s () , m3sq () , m1sq () , m2sq () ) ;
  //
  return -4 * s() * g / ( f1 * f2 ) ;
}
// =============================================================================
// total momentum in (2,3) restframe
// =============================================================================
double Ostap::Kinematics::Dalitz::P_R23  ( const double /* s1 */ , const double s2 ) const
{
  const double f1 = Ostap::Kinematics::triangle ( s (), s2 , m1sq () ) ;
  return 0 < f1 ? 0.5 * std::sqrt ( f1 / s() ) : 0.0 ;
}
// =============================================================================
// momentum of 2nd particle in (2,3) restframe
// =============================================================================
double Ostap::Kinematics::Dalitz::P2_R23  ( const double /* s1 */ , const double s2 ) const
{
  const double f1 = Ostap::Kinematics::triangle ( s2 , m2sq () , m3sq () ) ;
  return 0 < f1 ? 0.5 * std::sqrt ( f1 / s2 ) : 0.0 ;
}
// ============================================================================
/*  cosine on the angle between 1st and 2nd particles in the  (2,3) rest frame
 *  \f$ \cos \theta_{12}^{R(2,3)}
 */
// ============================================================================
double Ostap::Kinematics::Dalitz::cos_12_R23 ( const double    s1    , const double s2 ) const
{
  //
  const double f1 = Ostap::Kinematics::triangle ( s () , s2      , m1sq () ) ;
  const double f2 = Ostap::Kinematics::triangle ( s2   , m2sq () , m3sq () ) ;
  //
  if ( 0 >= f1 || 0 >= f2 ) { return  -1 ; }
  //
  const double f  = ( s () - s2 - m1sq () ) * ( s2 + m2sq () - m3sq () ) 
    + 2 * s2 * ( m1sq () + m2sq () - s1 ) ;
  //
  return f / std::sqrt ( f1 * f2 ) ;
}
// ============================================================================
/*  sine squared  on the angle between 1st and 2nd particles in the  (2,3) rest frame
 *  \f$ \sin^2 \theta_{12}^{R(2,3)}
 */
// ============================================================================
double Ostap::Kinematics::Dalitz::sin2_12_R23 ( const double    s1    , const double s2 ) const 
{
  const double f1 = Ostap::Kinematics::triangle ( s () , s2      , m1sq ()  ) ;
  const double f2 = Ostap::Kinematics::triangle ( s2   , m2sq () , m3sq ()  ) ;  
  const double g  = Ostap::Kinematics::G ( s1 , s2 , s () , m2sq () , m1sq () , m3sq () ) ;
  //
  return -4 * s2 * g / ( f1 * f2 ) ;
}
// =============================================================================
// total momentum in (1,2) restframe
// =============================================================================
double Ostap::Kinematics::Dalitz::P_R12  ( const double s1 , const double /* s2 */ ) const
{
  const double f1 = Ostap::Kinematics::triangle ( s (), s1 , m3sq () ) ;
  return 0 < f1 ? 0.5 * std::sqrt ( f1 / s() ) : 0.0 ;
}
// =============================================================================
// momentum of 1st particle in (1,2) restframe
// =============================================================================
double Ostap::Kinematics::Dalitz::P1_R12  ( const double s1 , const double /* s2 */ ) const
{
  const double f1 = Ostap::Kinematics::triangle ( s1 , m1sq() , m2sq () ) ;
  return 0 < f1 ? 0.5 * std::sqrt ( f1 / s1 ) : 0.0 ;
}
// ============================================================================
/*  cosine on the angle between 3rd and 1st particles in the  (1,2) rest frame
 *  \f$ \cos \theta_{21}^{R(1,2)}
 */
// ============================================================================
double Ostap::Kinematics::Dalitz::cos_31_R12 ( const double s1 , const double s2 ) const
{
  //
  const double s3_ = s3 ( s1 , s2 ) ;
  //
  const double f1 = Ostap::Kinematics::triangle ( s () , s1      , m3sq () ) ;
  const double f2 = Ostap::Kinematics::triangle ( s1   , m1sq () , m2sq () ) ;
  //
  if ( 0 >= f1 || 0 >= f2 ) { return  -1 ; }
  //
  const double f  = ( s() - s1 - m3sq()  ) * ( s1 + m1sq() - m2sq ()  ) 
    + 2 * s1 * ( m3sq()  + m1sq()  - s3_ ) ;
  //
  return f / std::sqrt ( f1 * f2 ) ;
}
// ============================================================================
/*  sine squared  on the angle between 3rd and 1st particles in the  (1,2) rest frame
 *  \f$ \sin^2 \theta_{12}^{R(2,3)}
 */
// ============================================================================
double Ostap::Kinematics::Dalitz::sin2_31_R12 ( const double s1 , const double s2 ) const 
{
  const double s3_ = s3 ( s1 , s2 ) ;
  //
  const double f1 = Ostap::Kinematics::triangle ( s () , s1      , m3sq () ) ;
  const double f2 = Ostap::Kinematics::triangle ( s1   , m1sq () , m2sq () ) ;  
  const double g  = Ostap::Kinematics::G ( s3_ , s1 , s () , m1sq () , m3sq () , m2sq () ) ;
  //
  return -4 * s1 * g / ( f1 * f2 ) ;
}
// =============================================================================
// total momentum in (3,1) restframe
// =============================================================================
double Ostap::Kinematics::Dalitz::P_R31  ( const double s1 , const double s2 ) const
{
  const double s3_ = s3 ( s1 , s2 ) ;
  //
  const double f1 = Ostap::Kinematics::triangle ( s (), s3_ , m2sq () ) ;
  return 0 < f1 ? 0.5 * std::sqrt ( f1 / s() ) : 0.0 ;
}
// =============================================================================
// momentum of 3rd particle in (3,1) restframe
// =============================================================================
double Ostap::Kinematics::Dalitz::P3_R31  ( const double s1 , const double s2 ) const
{
  const double s3_ = s3 ( s1 , s2 ) ;
  //
  const double f1 = Ostap::Kinematics::triangle ( s3_ , m3sq () , m1sq () ) ;
  return 0 < f1 ? 0.5 * std::sqrt ( f1 / s2 ) : 0.0 ;
}
// ============================================================================
/* cosine on the angle between 2nd and 3rd particles in the  (3,1) rest frame
 *  \f$ \cos \theta_{23}^{R(3,1)}
 */
// ============================================================================
double Ostap::Kinematics::Dalitz::cos_23_R31 ( const double    s1    , const double s2 ) const 
{
  const double s3_ = s3 ( s1 , s2 ) ;
  //
  const double f1 = Ostap::Kinematics::triangle ( s () , s3_     , m2sq () ) ;
  const double f2 = Ostap::Kinematics::triangle ( s3_  , m3sq () , m1sq () ) ;
  //
  if ( 0 >= f1 || 0 >= f2 ) { return  -1 ; }
  //
  const double f  = ( s() - s3_ - m2sq () ) * ( s3_ + m3sq () - m1sq () ) 
    + 2 * s3_ * ( m2sq () + m3sq() - s2 ) ;
  //
  return f / std::sqrt ( f1 * f2 ) ;
}
// ============================================================================
/*  sine squared  on the angle between 2nd and 3rd particles in the  (3,1) rest frame
 *  \f$ \sin^2 \theta_{23}^{R(3,1)}
 */
// ============================================================================
double Ostap::Kinematics::Dalitz::sin2_23_R31 ( const double    s1    , const double s2 ) const 
{
  const double s3_ = s3 ( s1 , s2 ) ;
  //
  const double f1 = Ostap::Kinematics::triangle ( s () , s3_     , m2sq () ) ;
  const double f2 = Ostap::Kinematics::triangle ( s3_  , m3sq () , m1sq () ) ;  
  const double g  = Ostap::Kinematics::G ( s2 , s3_ , s () , m3sq () , m2sq () , m1sq () ) ;
  //
  return -4 * s3_ * g / ( f1 * f2 ) ;
}
// ============================================================================
/*  Dalitz plot density:
 *  \f$ R_3 = 
 *  \frac{1}{32s} \int \deriv s_1 \deriv s_2 \deriv Omega \deriv \phi_3 
 *  \Theta { -G \left( s_1, s_2 , s , m_2^2, m_1^2, m_3^2 \right) } \f$
 *  Here we take 
 *  \f$ \int \deriv \Omega = 4\pi\f$ and \f$ \int \deriv \phi_3 = 2\pi\f$. 
 */
// ============================================================================
double Ostap::Kinematics::Dalitz::density ( const double s1 , const double s2 ) const 
{
  static const double s_norm = 0.25 * M_PI * M_PI ;
  //
  if ( s1  < s1_min () || s1  > s1_max () ) { return 0 ; }
  if ( s2  < s2_min () || s2  > s2_max () ) { return 0 ; }
  const double s3_ =  s3 ( s1 , s2  ) ;
  if ( s3_ < s3_min () || s3_ > s3_max () ) { return 0 ; }
  //
  const double g = Ostap::Kinematics::G ( s1 , s2 , s () , m2sq () , m1sq () , m3sq () ) ;
  //
  return g < 0 ? s_norm / s ()  : 0.0 ;  
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
  const double f1 = Ostap::Kinematics::triangle ( s3 , s ()        , m2sq () ) ;
  if ( f1 < 0 ) {  return 0 ; }
  const double f2 = Ostap::Kinematics::triangle ( s3 , m_m3 * m_m3 , m1sq () ) ;  
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
  const double f1 = Ostap::Kinematics::triangle ( s1 , s ()        , m3sq () ) ;
  if ( f1 < 0 ) {  return 0 ; }
  const double f2 = Ostap::Kinematics::triangle ( s1 , m_m1 * m_m1 , m2sq () ) ;  
  if ( f1 < 0 ) {  return 0 ; }
  //
  static const double s_norm = 0.25 * M_PI * M_PI ;
  //
  return 0 < f1 && 0 < f2 ? s_norm * std::sqrt ( f1 * f2 ) / ( s () * s1 ) : 0.0 ; 
}
// ============================================================================
// Dalitz plot boundaries \f$ s_1^{min/max} ( s_2 ) \f$ 
// ============================================================================
std::pair<double,double>  
Ostap::Kinematics::Dalitz::s1_minmax_for_s2 
( const double s2 ) const 
{
  //
  const bool   m12_zero  = s_zero ( m_m1 ) && s_zero ( m_m2 ) ;
  const bool   m23_zero  = s_zero ( m_m2 ) && s_zero ( m_m3 ) ;
  const bool   m31_zero  = s_zero ( m_m3 ) && s_zero ( m_m1 ) ;
  //
  const bool   all_zero  = m12_zero        && m23_zero        ;
  //
  if ( all_zero )
  {
    //
    const double s_min = 0    ;
    const double s_max = s () ;
    //
    if      ( s_equal ( s2 , s_min ) ) { return std::make_pair ( s_min , s_max ) ; }
    else if ( s2 < s_min             ) { return std::make_pair ( 1.0   , -1.0  ) ; }
    else if ( s_equal ( s2 , s_max ) ) { return std::make_pair ( s_max ,  0.0  ) ; }
    else if ( s2 > s_max             ) { return std::make_pair ( 1.0   , -1.0  ) ; }
    //
    return std::make_pair ( 0.0 , s_max - s2 ) ; 
  }
  //
  if  ( m12_zero )
  {
    //
    const double s_min = m3sq () ;
    const double s_max = s () ;
    //
    if      ( s_equal ( s2 , s_min ) ) { return std::make_pair ( s_min ,  0.0 ) ; }
    else if ( s2 < s_min             ) { return std::make_pair ( 1.0   , -1.0 ) ; }
    else if ( s_equal ( s2 , s_max ) ) { return std::make_pair ( s_max ,  0.0 ) ; }
    else if ( s2 > s_max             ) { return std::make_pair ( 1.0   , -1.0 ) ; }
    //
    return  std::make_pair ( 0.0 , ( s() - s2 ) * ( s2 - m3sq () ) / s2 )  ;
  }
  //
  if ( m23_zero ) 
  {
    const double s_min = 0 ;
    const double s_max = ( sqs () - m_m1  ) * ( sqs () - m_m1  ) ;
    //
    if      ( s_equal ( s2 , s_min ) ) { return std::make_pair ( m1sq () , s()     ) ; }
    else if ( s2 < s_min             ) { return std::make_pair ( 1.0 , -1.0        ) ; }
    else if ( s_equal ( s2 , s_max ) ) { return std::make_pair ( m_m1  * sqs () , m_m1  * sqs () ) ; }
    else if ( s2 > s_max             ) { return std::make_pair ( 1.0 , -1.0        ) ; }
    //
  }
  //
  if ( m31_zero ) 
  {
    //
    const double s_min = m2sq () ;
    const double s_max = s () ;
    //
    if      ( s_equal ( s2 , s_min ) ) { return std::make_pair ( s_min , s_max ) ; }
    else if ( s2 < s_min             ) { return std::make_pair ( 1.0   , -1.0  ) ; }
    else if ( s_equal ( s2 , s_max ) ) { return std::make_pair ( s_max , s_min ) ; }
    else if ( s2 > s_max             ) { return std::make_pair ( 1.0   , -1.0  ) ; }
    //
    return  std::make_pair ( s_max * s_min / s2 , s_max +  s_min - s2 ) ; 
  }
  //
  const double s_a       = m1sq()  + m2sq()  ;
  //
  const double s_min     = s2_min () ;
  const double s_max     = s2_max () ;
  //
  if ( s2 < s_min || s2 > s_max ) { return std::make_pair ( 1.0 , -1.0 ) ; }
  //
  const double f1 = Ostap::Kinematics::triangle ( s2 , s ()    , m1sq () ) ;
  const double f2 = Ostap::Kinematics::triangle ( s2 , m2sq () , m3sq () ) ;
  //
  if ( f1 < 0 || f2 < 0 ) { return std::make_pair ( 1.0 , -1.0 ) ; }
  //
  const double b  = ( s2 - s() + m1sq() ) * ( s2 + m2sq() - m3sq () ) ;
  const double c  = std::sqrt ( f1 * f2 ) ;
  //
  const double q1 = s_a - ( b - c ) / ( 2 * s2 ) ;
  const double q2 = s_a - ( b + c ) / ( 2 * s2 ) ;
  //
  return std::make_pair ( q2 , q1 ) ;
}
// ============================================================================
// Dalitz plot boundaries \f$ s_2^{min/max} ( s_1 ) \f$ 
// ============================================================================
std::pair<double,double>  
Ostap::Kinematics::Dalitz::s2_minmax_for_s1 
( const double s1 ) const 
{
  //
  const bool   m12_zero  = s_zero ( m_m1 ) && s_zero ( m_m2 ) ;
  const bool   m23_zero  = s_zero ( m_m2 ) && s_zero ( m_m3 ) ;
  const bool   m31_zero  = s_zero ( m_m3 ) && s_zero ( m_m1 ) ;
  //
  const bool   all_zero  = m12_zero        && m23_zero        ;
  //
  if ( all_zero )
  {
    //
    const double s_min = 0    ;
    const double s_max = s () ;
    //
    if      ( s_equal ( s1 , s_min ) ) { return std::make_pair ( s_min , s_max ) ; }
    else if ( s1 < s_min             ) { return std::make_pair ( 1.0   , -1.0  ) ; }
    else if ( s_equal ( s1 , s_max ) ) { return std::make_pair ( s_max ,  0.0  ) ; }
    else if ( s1 > s_max             ) { return std::make_pair ( 1.0   , -1.0  ) ; }
    //
    return std::make_pair ( 0.0 , s_max - s1 ) ; 
  }
  //
  if  ( m23_zero )
  {
    //
    const double s_min = m1sq () ;
    const double s_max = s () ;
    //
    if      ( s_equal ( s1 , s_min ) ) { return std::make_pair ( s_min ,  0.0 ) ; }
    else if ( s1 < s_min             ) { return std::make_pair ( 1.0   , -1.0 ) ; }
    else if ( s_equal ( s1 , s_max ) ) { return std::make_pair ( s_max ,  0.0 ) ; }
    else if ( s1 > s_max             ) { return std::make_pair ( 1.0   , -1.0 ) ; }
    //
    return  std::make_pair ( 0.0 , ( s() - s1 ) * ( s1 - m1sq () ) / s1 )  ;
  }
  //
  if ( m12_zero ) 
  {
    const double s_min = 0 ;
    const double s_max = s1_max ()  ;
    //
    if      ( s_equal ( s1 , s_min ) ) { return std::make_pair ( m3sq() , s() ) ; }
    else if ( s1 < s_min             ) { return std::make_pair ( 1.0 , -1.0        ) ; }
    else if ( s_equal ( s1 , s_max ) ) { return std::make_pair ( m_m3  * sqs () , m_m3  * sqs () ) ; }
    else if ( s1 > s_max             ) { return std::make_pair ( 1.0 , -1.0        ) ; }
    //
  }
  if ( m31_zero ) 
  {
    //
    const double s_min = m2sq () ;
    const double s_max = s () ;
    //
    if      ( s_equal ( s1 , s_min ) ) { return std::make_pair ( s_min , s_max ) ; }
    else if ( s1 < s_min             ) { return std::make_pair ( 1.0   , -1.0  ) ; }
    else if ( s_equal ( s1 , s_max ) ) { return std::make_pair ( s_max , s_min ) ; }
    else if ( s1 > s_max             ) { return std::make_pair ( 1.0   , -1.0  ) ; }
    //
    return  std::make_pair ( s_max * s_min / s1 , s_max +  s_min - s1 ) ; 
  }
  //
  const double s_a       = m3sq () + m2sq()  ;
  //
  const double s_min     = s1_min () ;
  const double s_max     = s1_max () ;
  //
  if ( s1 < s_min || s1 > s_max ) { return std::make_pair ( 1.0 , -1.0 ) ; }
  //
  const double f1 = Ostap::Kinematics::triangle ( s1 , s ()    , m3sq () ) ;
  const double f2 = Ostap::Kinematics::triangle ( s1 , m2sq () , m1sq () ) ;
  //
  if ( f1 < 0 || f2 < 0 ) { return std::make_pair ( 1.0 , -1.0 ) ; }
  //
  const double b  = ( s1 - s() +  m3sq () ) * ( s1 + m2sq () - m1sq()  ) ;
  const double c  = std::sqrt ( f1 * f2 ) ;
  //
  const double q1 = s_a - ( b - c ) / ( 2 * s1 ) ;
  const double q2 = s_a - ( b + c ) / ( 2 * s1 ) ;
  //
  return std::make_pair ( q2 , q1 ) ; 
}
// ============================================================================
//                                                                      The END 
// ============================================================================
