// Include files 
// ============================================================================
// STD&STL
// ============================================================================
#include <functional>
#include <tuple>
// ============================================================================
// local
// ============================================================================
#include "Ostap/Dalitz.h"
#include "Ostap/DalitzIntegrator.h"
// ============================================================================
// ROOT 
// ============================================================================
#include "Math/Integrator.h"
#include "Math/GSLIntegrator.h"
// ============================================================================
/** @file
 *  Implementation file for the class Ostap::Math::DalitzIntegrator
 *  @see class Ostap::Math::DalitzIntegral  
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-11-02
 */
// ============================================================================
// Constructor
// ============================================================================
Ostap::Math::DalitzIntegrator::DalitzIntegrator
( const double m1 , 
  const double m2 , 
  const double m3 ) 
  : m_m1 ( std::abs ( m1 ) ) 
  , m_m2 ( std::abs ( m2 ) ) 
  , m_m3 ( std::abs ( m3 ) ) 
{}
// ===========================================================================
namespace 
{
  // =========================================================================
  ROOT::Math::GSLIntegrator s_integrator_1D_1 
  ( ROOT::Math::Integration::kADAPTIVE ,
    ROOT::Math::Integration::kGAUSS51  , 1.e-6 , 1.e-6 , 20000 ) ; 
  // =========================================================================
  ROOT::Math::GSLIntegrator s_integrator_1D_2 
  ( ROOT::Math::Integration::kADAPTIVE ,
    ROOT::Math::Integration::kGAUSS51  , 1.e-6 , 1.e-6 , 20000 ) ; 
  // =========================================================================
  typedef std::tuple<Ostap::Math::DalitzIntegrator::function2*,
                     const Ostap::Kinematics::Dalitz*,double> DATA2 ;
  // =========================================================================
  double inner_function  ( const double s2 , void* param ) 
  {
    const DATA2* data     = static_cast<DATA2*> ( param ) ;
    auto*        function = std::get<0>(*data) ;
    const auto*  dalitz   = std::get<1>(*data) ;
    const double s1       = std::get<2>(*data) ;
    //
    if ( !dalitz->inside ( s1 , s2 ) ) { return 0 ; }
    //
    return (*function) ( s1 , s2 ) ;
  }
  // 
  double inner_integral  ( const double s1 , void* param ) 
  {
    DATA2*       data     = static_cast<DATA2*> ( param ) ;
    auto*        function = std::get<0>(*data) ;
    const auto*  dalitz   = std::get<1>(*data) ;
    //
    double s2min, s2max ;
    std::tie ( s2min , s2max ) = dalitz -> s2_minmax_for_s1 ( s1 ) ;
    if ( s2max <= s2min ) { return 0 ; }
    //
    std::get<2> ( *data ) = s1 ;
    // =========================================================================
    return s_integrator_1D_1.Integral ( &inner_function , 
                                        param           , 
                                        s2min           ,
                                        s2max           ) ;
  }
  // ===========================================================================
}
// ===========================================================================
/*  evaluate the integral over \f$s_1\f$ , \f$s_2\f$ variables 
 *  \f[ \int\int ds_1 ds_2 f(s, s_1,s_2) \f] 
 *  @param M the overall mass of the system \f$\sqrt{s} \f$
 *  @param fs the function \f$ f(s, s_1,s_2)\f$ 
 *  @return the integral over dalitz plot
 */
// ==========================================================================
double Ostap::Math::DalitzIntegrator::integrate_s1s2
( const double                             M , 
  Ostap::Math::DalitzIntegrator::function3 f ) const 
{
  auto      f1 = std::cref ( f      )  ;
  function2 f2 = std::bind ( f1 , M * M , std::placeholders::_1 , std::placeholders::_2 ) ;
  return integrate_s1s2 ( M , f2 ) ;  
}
// ===========================================================================
/*  evaluate the integral over \f$s_1\f$ , \f$s_2\f$ variables 
 *  \f[ \int\int ds_1 ds_2 f(s_1,s_2) \f] 
 *  @param M the overal mass of the system \f$ \sqrt{s} \f$
 *  @param fs the function \f$ f(s_1,s_2)\f$ 
 *  @return the integral over dalitz plot
 */
// ===========================================================================
double Ostap::Math::DalitzIntegrator::integrate_s1s2 
( const double M                           , 
  Ostap::Math::DalitzIntegrator::function2 f ) const 
{
  if ( M <= m_m1 + m_m2 + m_m3 ) { return 0.0 ; }
  // create the calculator 
  const Ostap::Kinematics::Dalitz dalitz { M , m_m1 , m_m2 , m_m3 } ;  
  // 
  DATA2 data = std::make_tuple( &f , &dalitz , 0.0 ) ;
  void* param = &data ;
  //
  return s_integrator_1D_2.Integral ( &inner_integral  , 
                                      param            ,
                                      dalitz.s1_min () , 
                                      dalitz.s1_max () ) ;
}  
// ============================================================================
/*  evaluate the integral over \f$e_2\f$ , \f$e_3\f$ variables 
 *  \f[ \int\int de_2 de_3 f(s,e_2,e_3) = 
 *  \int_{e_2^{min}}^{e_2^{max}} de_2 
 *  \int_{e_3^{min}(e_2)}^{e_3^{max}(e_2)} de_3 f(s,e_2,e_3)  \f] 
 *  @param  M the overall mass of the system \f$ \sqrt{s} \f$
 *  @param  f the function \f$ f(e_2,e_3)\f$ 
 *  @return the integral over dalitz plot
 */
// ==============================================================================
double Ostap::Math::DalitzIntegrator::integrate_e2e3 
( const double                             M , 
  Ostap::Math::DalitzIntegrator::function3 f ) const 
{
  auto      f1 = std::cref ( f      )  ;
  function2 f2 = std::bind ( f1 , M * M , std::placeholders::_1 , std::placeholders::_2 ) ;
  return integrate_e2e3 ( M , f2 ) ;  
} 
// ============================================================================
/*  evaluate the integral over \f$e_2\f$ , \f$e_3\f$ variables 
 *  \f[ \int\int de_2 de_3 f(e_2,e_3) = 
 *  \int_{e_2^{min}}^{e_2^{max}} de_2 
 *  \int_{e_3^{min}(e_2)}^{e_3^{max}(e_2)} de_3 f(e_2,e_3)  \f] 
 *  @param  M the overall mass of the system \f$ \sqrt{s} \f$
 *  @param  f the function \f$ f(e_2,e_3)\f$ 
 *  @return the integral over dalitz plot
 */
// ==============================================================================
double Ostap::Math::DalitzIntegrator::integrate_e2e3 
( const double                           M , 
  Ostap::Math::DalitzIntegrator::function2 f ) const 
{
  if ( M <= m_m1 + m_m2 + m_m3 ) { return 0.0 ; }
  // create the calculator 
  const Ostap::Kinematics::Dalitz dalitz { M , m_m1 , m_m2 , m_m3 } ;  
  // 
  auto         ff = std::cref ( f )  ;
  const double J  = 0.25 / ( M * M ) ; // jacobian 
  function2 fun = [&dalitz,&ff,J] ( const double s1 , const double s2 ) -> double 
    {
      const double e2 = dalitz.E2 ( s1 , s2 ) ;
      const double e3 = dalitz.E3 ( s1 , s2 ) ;
      //
      if ( e2 <= 0 || e3 <= 0 || !dalitz.inside ( s1 , s2 ) ) { return 0 ; }
      //
      return ff ( e2 , e3 ) * J ;
    } ;
  //
  return integrate_s1s2 ( M , fun ) ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================

