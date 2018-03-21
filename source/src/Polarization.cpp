// ============================================================================
// Include files 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Polarization.h"
#include "Ostap/Tensors.h"
// ============================================================================
// ROOT
// ============================================================================
#include "Math/Boost.h"
// ============================================================================
/** @file 
 *  Implementation file for functions from namespace Ostap::Polarization
 *  @date 2018-03-20 
 *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
 */
// ============================================================================
/*  Boost LorentzVector into rest-frame of another Lorentz vector 
 *  @param what   the vextro to be bosted 
 *  @param frame  the 4-vector of the frame 
 *  @return boosted vector 
 */
// ============================================================================
Ostap::LorentzVector Ostap::Math::boost 
( const Ostap::LorentzVector& what  ,
  const Ostap::LorentzVector& frame )
{
  const ROOT::Math::Boost b { frame.BoostToCM() } ;
  return b ( what ) ;
}
// ============================================================================
namespace 
{
  typedef  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<long double> > LV ;
  // ==========================================================================
  Ostap::Math::Polarization::Frame  
  _frame_ ( const double az ,
            const double bz , 
            const LV&    AT , 
            const LV&    BT , 
            const LV&    P  ) 
  {
    const LV Z ( az * AT + bz * BT ) ;
    //
    const long double atbt = AT.Dot ( BT ) ;
    const long double atat = AT.M2  () ;
    const long double btbt = AT.M2  () ;
    // 
    const long double ax = -az * atbt + bz * btbt ;
    const long double bx =  az * atat + bz * atbt ;
    //
    LV X ( ax * AT + bx * BT ) ;
    if   ( 0 < az * bx - ax * bz  ) { X /=      std::abs ( X.M() ) ; }
    else                            { X /= -1 * std::abs ( X.M() ) ; }
    //
    const LV Y ( Ostap::Math::Tensors::Epsilon::epsilon( P , X , Z ) / P.M () ) ;
    //
    typedef Ostap::LorentzVector VL ;
    return {{ VL( X ) , VL( Y)  , VL ( Z ) }} ;  
  }
  // ==========================================================================
}
// ============================================================================
/* get three  axes for polarization frame 
 *  @param f  the frame 
 *  @param P     4-momenta of the particle 
 *  @param beam1 4-momenta of the first colliding particle
 *  @param beam2 4-momenta of the second colliding particle
 *  @return three polarization axes: x,y&z
 */
// ============================================================================
Ostap::Math::Polarization::Frame  
Ostap::Math::Polarization::frame
( Ostap::Math::Polarization::Frames f     , 
  const Ostap::LorentzVector&       p     ,
  const Ostap::LorentzVector&       beam1 ,
  const Ostap::LorentzVector&       beam2 )
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
      const long double az_GJ = -1/(AT+BT).M() ;         // ATTENTION HERE! 
      const long double bz_GJ = az_GJ          ;         // ATTENTION HERE!
      return _frame_ ( az_GJ , bz_GJ , AT , BT , P ) ; 
    }
    // ========================================================================
  case Ostap::Math::Polarization::Frames::Target :
    {
      const long double az_TA = 1/(AT-BT).M() ;         // ATTENTION HERE!
      const long double bz_TA = -az_TA        ;         // ATTENTION HERE!
      return _frame_ ( az_TA , bz_TA , AT , BT , P ) ; 
    }
    // ========================================================================
  case Ostap::Math::Polarization::Frames::CollinsSoper :
    {
      const long double d = 1/ ( -BP*AT+AP*BT).M() ;
      const long double az_CS =  BP*d ;                 // ATTENTION HERE  
      const long double bz_CS = -AP*d ;                 // ATTENTION HERE 
      return _frame_ ( az_CS , bz_CS , AT , BT , P ) ; 
    }
  case Ostap::Math::Polarization::Frames::Recoil : ;
  default : ;
  }
  // ==========================================================================
  /// helicity frame:
  const long double az_HX = 1/AT.M() ;                 // ATTENTION HERE! 
  const long double bz_HX = 0        ;                 // ATTENTION HERE!
  //
  return _frame_ ( az_HX , bz_HX , AT , BT , P ) ; 
}
// ============================================================================
/** get the angles \$ (\cos \theta,\phi)\$ for the particle 
 *  in the rest frame of particle m,  and the beam-momenta p1& p2  
 *  @param p the particle
 *  @param f the frame 
 *  @param m the particle that defined the frame 
 *  @param p1 4-momenta of the first colliding particle
 *  @param p2 4-momenta of the second colliding particle
 *  @return  \$ (\cos \theta,\phi)\$ structure
 */
// ============================================================================
Ostap::Math::Polarization::Angles 
Ostap::Math::Polarization::angles
( const Ostap::LorentzVector&       p     , 
  Ostap::Math::Polarization::Frames f     ,   
  const Ostap::LorentzVector&       m     , 
  const Ostap::LorentzVector&       beam1 , 
  const Ostap::LorentzVector&       beam2 )
{ return angles ( p , frame ( f , m , beam1 , beam2 ) ) ; }
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
  const long double cosx = -f[0].Dot ( p ) ;
  const long double cosy = -f[1].Dot ( p ) ;
  const long double cosz = -f[2].Dot ( p ) ;
  const long double r    = std::sqrt ( cosx * cosx + cosy * cosy + cosz * cosz ) ;
  //
  return {{ double ( cosz/r ) , double( std::atan2 ( cosy , cosx ) ) }} ;
}
// ============================================================================


// ============================================================================
//                                                                      The END 
// ============================================================================
