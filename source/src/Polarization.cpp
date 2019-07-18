// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================
#include <cmath>
#include <climits>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Polarization.h"
#include "Ostap/Tensors.h"
#include "Ostap/Kinematics.h"
// ============================================================================
/** @file 
 *  Implementation file for functions from namespace Ostap::Polarization
 *  @date 2018-03-20 
 *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
 */
// ============================================================================
namespace 
{ 
  // ==========================================================================
  static_assert ( std::numeric_limits<float> ::is_specialized      ,
                  "std::numeric_limits<float>  is not specialized" ) ;
  /// large negative number 
  constexpr double s_INVALID = -0.9 * std::numeric_limits<float>::max () ;
  static_assert (  s_INVALID <  0   , "invalid negative number"    ) ;
  // ==========================================================================
}
// ============================================================================
namespace 
{
  typedef  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<long double> > LV ;
  // ==========================================================================
  Ostap::Math::Polarization::Frame  
  _frame_ ( const long double az      ,
            const long double bz      , 
            const LV&         AT      , 
            const LV&         BT      , 
            const LV&         P       , 
            const Ostap::Math::Polarization::UseMadisonConvention madison ) 
  {
    const LV Z ( az * AT + bz * BT ) ;
    //
    const long double atat = AT.M2  ( )    ;
    const long double atbt = AT.Dot ( BT ) ;
    const long double btbt = BT.M2  ( )    ;
    // 
    const long double ax = -az * atbt - bz * btbt ;
    const long double bx =  az * atat + bz * atbt ;
    //
    LV X ( ax * AT + bx * BT ) ;
    //
    const long double ixm = madison ? 1 / std::abs ( X.M() ) : -1 / std::abs ( X.M() ) ;  
    if   ( 0 < az * bx - ax * bz  ) { X *=  ixm ; }
    else                            { X *= -ixm ; }
    //
    // const short s =  0 < az * bx - ax * bz ? 1 : -1 ;  
    // const LV X ( ( s * ax ) * AT + ( s * bx ) * BT ) ;
    //
    const long double im = 1/P.M() ;
    const LV Y ( Ostap::Math::Tensors::Epsilon::epsilon( P , X , Z ) * -im  ) ;
    //
    typedef Ostap::LorentzVector VL ;
    return {{  VL( X ) , VL( Y ) , VL ( Z ) , 
          VL ( Ostap::Math::Tensors::Epsilon::epsilon ( X , Y , Z ) ) }} ;  
  }
  // ==========================================================================
}
// ============================================================================
/* get three  axes for polarization frame 
 *  @param f  the frame 
 *  @param P     4-momenta of the particle 
 *  @param beam1 4-momenta of the first colliding particle
 *  @param beam2 4-momenta of the second colliding particle
 *  @param madison use Madison convention?
 *  @return three polarization axes: x,y&z
 */
// ============================================================================
Ostap::Math::Polarization::Frame  
Ostap::Math::Polarization::frame
( const Ostap::Math::Polarization::Frames               f       , 
  const Ostap::LorentzVector&                           p       ,
  const Ostap::LorentzVector&                           beam1   ,
  const Ostap::LorentzVector&                           beam2   ,
  const Ostap::Math::Polarization::UseMadisonConvention madison ) 
{
  //
  const LV P1 ( beam1   ) ;
  const LV P2 ( beam2   ) ;
  
  const LV A  ( P1 + P2 ) ;
  const LV B  ( P1 - P2 ) ;
  
  const LV P  ( p )       ;
  
  ///  the particle mass  squared 
  const long double M2 = P.M2() ;  
  
  const long double AP = A.Dot ( P ) ;
  const long double BP = B.Dot ( P ) ;
  
  const LV AT ( A - ( AP / M2 ) * P ) ;
  const LV BT ( B - ( BP / M2 ) * P ) ;  
  //
  switch ( f ) 
  {
    // ========================================================================
  case Ostap::Math::Polarization::Frames::GottfriedJackson :
    {
      const long double az_GJ = -1/(AT+BT).M() ;             // ATTENTION HERE! 
      const long double bz_GJ = az_GJ          ;             // ATTENTION HERE!
      return _frame_ ( az_GJ , bz_GJ , AT , BT , P , madison ) ; 
    }
    // ========================================================================
  case Ostap::Math::Polarization::Frames::Target :
    {
      const long double az_TA =  1 / ( AT - BT ) . M() ;     // ATTENTION HERE!
      const long double bz_TA = -az_TA        ;              // ATTENTION HERE!
      return _frame_ ( az_TA , bz_TA , AT , BT , P , madison ) ; 
    }
    // ========================================================================
  case Ostap::Math::Polarization::Frames::CollinsSoper :
    {
      const long double c1 = 1/(AT-BT).M() ;
      const long double c2 = 1/(AT+BT).M() ;
      const long double a1 = c1 - c2       ;
      const long double a2 = c1 + c2       ;
      const long double d = 1 / ( a1 * AT - a2 * BT ).M() ;
      const long double az_CS =  -a1 * d ;                    // ATTENTION HERE  
      const long double bz_CS =   a2 * d ;                    // ATTENTION HERE 
      return _frame_ ( az_CS , bz_CS , AT , BT , P , madison ) ; 
    }
  case Ostap::Math::Polarization::Frames::Recoil : ;
  default : ;
  }
  // ==========================================================================
  /// helicity frame:
  const long double az_HX = 1 / AT.M() ;                     // ATTENTION HERE! 
  const long double bz_HX = 0          ;                     // ATTENTION HERE!
  //
  return _frame_ ( az_HX , bz_HX , AT , BT , P , madison ) ; 
}
// ============================================================================
/*  get the direction cosines of the particle direction 
 *  in the specified reference frame 
 *  @param p the particle
 *  @param f the frame 
 *  @return  direction  cosines 
 */
// ============================================================================
Ostap::Math::Polarization::Cosines 
Ostap::Math::Polarization::cosines 
( const Ostap::LorentzVector&             p , 
  const Ostap::Math::Polarization::Frame& f ) 
{
  const double irm = -1 / Ostap::Kinematics::restMomentum ( p , f[3] ) ;
  return {{  f[0].Dot(p)*irm , f[1].Dot(p)*irm ,f[2].Dot(p)*irm }} ;
}
// ============================================================================
/*  get the direction cosines of the particle direction 
 *  in the rest frame of particle m,  and the beam-momenta p1& p2  
 *  @param p the particle
 *  @param f the frame 
 *  @param m the particle that defined the frame 
 *  @param p1 4-momenta of the first colliding particle
 *  @param p2 4-momenta of the second colliding particle
 *  @param madison use Madison convention?
 *  @return  \$ (\cos \theta,\phi)\$ structure
 */
// ============================================================================
Ostap::Math::Polarization::Cosines 
Ostap::Math::Polarization::cosines 
( const Ostap::LorentzVector&                           p       , 
  const Ostap::Math::Polarization::Frames               f       , 
  const Ostap::LorentzVector&                           m       ,  
  const Ostap::LorentzVector&                           beam1   , 
  const Ostap::LorentzVector&                           beam2   ,
  const Ostap::Math::Polarization::UseMadisonConvention madison ) 
{ return cosines ( p , frame ( f , m , beam1 ,  beam2 , madison ) ) ; }
// ============================================================================
/*  get the angles \$ (\cos \theta,\phi)\$ for the particle 
 *  in the defined frame
 *  @param p the particle
 *  @param f the frame 
 *  @return  \$ (\cos \theta,\phi)\$ structure
 */
// ============================================================================
Ostap::Math::Polarization::Angles 
Ostap::Math::Polarization::angles
( const Ostap::LorentzVector&             p , 
  const Ostap::Math::Polarization::Frame& f ) 
{
  const Cosines cs = cosines ( p , f ) ;
  return { cs[2] , std::atan2 ( cs[1] , cs[0] ) } ;
}
// ============================================================================
/*  get the angles \$f(\cos \theta,\phi)\f$ for the particle 
 *  in the rest frame of particle m,  and the beam-momenta p1& p2  
 *  @param p the particle
 *  @param f the frame 
 *  @param m the particle that defined the frame 
 *  @param p1 4-momenta of the first colliding particle
 *  @param p2 4-momenta of the second colliding particle
 *  @param madison use Madison convention?
 *  @return  \f$ (\cos \theta,\phi)\f$ structure
 */
// ============================================================================
Ostap::Math::Polarization::Angles 
Ostap::Math::Polarization::angles
( const Ostap::LorentzVector&                           p       , 
  const Ostap::Math::Polarization::Frames               f       ,   
  const Ostap::LorentzVector&                           m       , 
  const Ostap::LorentzVector&                           beam1   , 
  const Ostap::LorentzVector&                           beam2   ,
  const Ostap::Math::Polarization::UseMadisonConvention madison ) 
{ return angles ( p , frame ( f , m , beam1 , beam2 , madison ) ) ; }
// ============================================================================
/*  get the polarization vectors from the frame 
 *  @param f the frame 
 *  @return polarization vectors (-1,0,+1)
 */
// ============================================================================
Ostap::Math::Polarization::PolVectors
Ostap::Math::Polarization::vectors 
( const Ostap::Math::Polarization::Frame& f ) 
{
  typedef Ostap::LorentzVector        VL  ;
  typedef Ostap::ComplexLorentzVector CLV ;
  const VL& ax = f[0] ;
  const VL& ay = f[1] ;
  const VL& az = f[2] ;
  //
  static const double               s_isq2 = 1/std::sqrt ( 2.0 ) ;
  static const std::complex<double> j { 0 , 1 } ;
  return {{ 
      CLV ( (  ax.X () - j * ay.X () ) * s_isq2 , 
            (  ax.Y () - j * ay.Y () ) * s_isq2 , 
            (  ax.Z () - j * ay.Z () ) * s_isq2 , 
            (  ax.T () - j * ay.T () ) * s_isq2 ) ,
        CLV ( az.X() ,  az.Y() ,  az.Z() ,  az.T () ) ,
        CLV ( ( -ax.X () - j * ay.X () ) * s_isq2 , 
              ( -ax.Y () - j * ay.Y () ) * s_isq2 , 
              ( -ax.Z () - j * ay.Z () ) * s_isq2 , 
              ( -ax.T () - j * ay.T () ) * s_isq2 ) }} ;
}
// ============================================================================


// ============================================================================
//                                                                      The END 
// ============================================================================
