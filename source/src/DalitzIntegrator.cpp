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
// Local
// ============================================================================
#include "Integrator2D.h"
#include "local_math.h"
#include "local_gsl.h"
// ============================================================================
/** @file
 *  Implementation file for the class Ostap::Math::DalitzIntegrator
 *  @see class Ostap::Math::DalitzIntegral  
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-11-02
 */
// ============================================================================
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
  // ==========================================================================
}
// ============================================================================
namespace  
{
  // ==========================================================================
  /** @var s_CUBATURE 
   *  The actual cubatire-integrator 
   */
  const Ostap::Math::GSL::Integrator2D<Ostap::Math::DalitzIntegrator::function2> s_DI ;
  // ==========================================================================  
  /** @var s_MESSAGE  
   *  Error message fromm cubature 
   */
  const char s_MESSAGE[]  = "Integrate(Dalitz)" ;
  // ==========================================================================  
}
// ============================================================================
/*  evaluate the integral over \f$s_1\f$ , \f$s_2\f$ variables 
 *  \f[ \int\int ds_1 ds_2 f(s, s_1,s_2) \f] 
 *  @param M the overall mass of the system \f$\sqrt{s} \f$
 *  @param fs the function \f$ f(s, s_1,s_2)\f$ 
 *  @return the integral over dalitz plot
 */
// ==========================================================================
double Ostap::Math::DalitzIntegrator::integrate_s1s2
( Ostap::Math::DalitzIntegrator::function3 f3 ,
  const Ostap::Kinematics::Dalitz&         d  )
{
  auto      f1 = std::cref ( f3   ) ;
  function2 f2 = std::bind ( f1 , d .s() , std::placeholders::_1 , std::placeholders::_2 ) ;
  return integrate_s1s2 ( std::cref ( f2 ) , d.s() , d  ) ;  
}
// ============================================================================
/*  evaluate the integral over \f$s_1\f$ , \f$s_2\f$ variables 
 *  \f[ \int\int ds_1 ds_2 f(s, s_1,s_2) \f] 
 *  @param M the overall mass of the system \f$\sqrt{s} \f$
 *  @param fs the function \f$ f(s, s_1,s_2)\f$ 
 *  @return the integral over dalitz plot
 */
// ==========================================================================
double Ostap::Math::DalitzIntegrator::integrate_s1s2
( Ostap::Math::DalitzIntegrator::function3 f3 ,
  const double                             s  ,
  const Ostap::Kinematics::Dalitz0&        d  )
{
  if ( s  <= d.sqsumm() ) { return 0 ; }
  //
  auto      f1 = std::cref ( f3   ) ;
  function2 f2 = std::bind ( f1 , s , std::placeholders::_1 , std::placeholders::_2 ) ;
  return integrate_s1s2 ( std::cref ( f2 ) , s , d  ) ;  
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
( Ostap::Math::DalitzIntegrator::function2 f2 ,
  const double                             s  ,
  const Ostap::Kinematics::Dalitz0&        d  )
{
  if ( s <= d.sqsumm () ) { return 0 ; }
  //
  const double M = std::sqrt ( s  ) ;
  //
  // function2  fun = [s,&d,&f2] ( const double s1  , const double s2 ) -> double
  //                  {
  //                    return ( s < s1 + s2 || !d.inside ( s , s1 , s2 ) ) ? 0.0 : f2 ( s1 , s2 ) ;
  //                  };  
  // //
  // /// get integrator
  // const auto F = s_DI.make_function
  //   ( &fun ,
  //     d.s1_min () , d.s1_max ( M ) ,
  //     d.s2_min () , d.s2_max ( M ) ) ;
  //
  function2  fun = [s,M,&d,&f2] ( const double x1  , const double x2 ) -> double
                   {
                     double s1, s2 ;                     
                     std::tie ( s1 , s2 ) = d.x2s ( s  , x1  , x2 ) ;
                     const double J = d.J ( s  , s1  , s2 ) ;
                     return J <= 0 ? 0.0 : f2 ( s1 , s2 ) * J ;
                   };  
  //
  /// get integrator
  const auto F = s_DI.make_function
    ( &fun , -1  , 1 , d.s2_min () , d.s2_max ( M ) ) ;
  //
  int    ierror  =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  std::tie ( ierror , result , error ) = s_DI.cubature
    ( &F , 20000 , s_PRECISION , s_PRECISION , s_MESSAGE , __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
/* evaluate the integral over \f$s\f$ , \f$s_1\f$ variables 
 *  \f[ \int\int f( s, s_1,s_2) ds ds_1 \f]  
 *  @param f3 the function \f$  f(s,s_1,s_2) \f$
 *  @param s2 fixed value of s2 
 *  @param smax upper-edge fo inntegrato ove \f$ s \f$
 *  @param d  helper Dalitz-object 
 *  @return integral over \f$ s, s_1\f$
 */
// ============================================================================
double Ostap::Math::DalitzIntegrator::integrate_ss1
( function3                         f3   ,
  const double                      s2   ,
  const double                      smax ,
  const Ostap::Kinematics::Dalitz0& d    )
{
  if ( s2   <= d.s1_min () ||
       smax <= d.sqsumm () || 
       s2   >= d.s2_max ( smax ) ) { return 0 ; }
  //
  const double smin = d.sqsumm() + s2  - d.s2_min() ;
  if ( smax  <= smin ) { return 0 ; }
  //
  function2  fun = [&d,s2,&f3] ( const double y1  , const double y2 ) -> double
                   {
                     double s, s1 ;
                     std::tie ( s , s1 ) = d.y2s ( s2 , y1 , y2 ) ;
                     const double J = d.J ( s  , s1  , s2 ) ;
                     return J <= 0 ? 0.0 : f3 ( s , s1 , s2 ) * J ;
                   };
  
  /// get the  integrator
  const auto F = s_DI.make_function ( &fun , smin , smax  , -1  , 1 ) ;
  //
  int    ierror  =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  std::tie ( ierror , result , error ) = s_DI.cubature
    ( &F , 20000 , s_PRECISION , s_PRECISION , s_MESSAGE , __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
/*  evaluate the integral over \f$s_1\f$ , \f$s_2\f$ variables 
 *  \f[ \int\int ds_1 ds_2 f(s, s_1,s_2) \f] 
 *  @param M the overall mass of the system \f$\sqrt{s} \f$
 *  @param fs the function \f$ f(s, s_1,s_2)\f$ 
 *  @return the integral over dalitz plot
 */
// ==========================================================================
double Ostap::Math::DalitzIntegrator::integrate_s1s2
( const std::size_t                        tag , 
  Ostap::Math::DalitzIntegrator::function3 f3 ,
  const Ostap::Kinematics::Dalitz&         d  )
{
  auto      f1 = std::cref ( f3   ) ;
  function2 f2 = std::bind ( f1 , d .s() , std::placeholders::_1 , std::placeholders::_2 ) ;
  const std::size_t key = std::hash_combine ( tag , d.tag() ) ;
  return integrate_s1s2 ( key  , std::cref ( f2 ) , d.s() , d  ) ;  
}
// ============================================================================
/*  evaluate the integral over \f$s_1\f$ , \f$s_2\f$ variables 
 *  \f[ \int\int ds_1 ds_2 f(s, s_1,s_2) \f] 
 *  @param M the overall mass of the system \f$\sqrt{s} \f$
 *  @param fs the function \f$ f(s, s_1,s_2)\f$ 
 *  @return the integral over dalitz plot
 */
// ==========================================================================
double Ostap::Math::DalitzIntegrator::integrate_s1s2
( const std::size_t                        tag , 
  Ostap::Math::DalitzIntegrator::function3 f3 ,
  const double                             s  ,
  const Ostap::Kinematics::Dalitz0&        d  )
{
  if ( s  <= d.sqsumm() ) { return 0 ; }
  //
  auto      f1 = std::cref ( f3   ) ;
  function2 f2 = std::bind ( f1 , s , std::placeholders::_1 , std::placeholders::_2 ) ;
  const std::size_t key = std::hash_combine ( tag , s , d.tag() ) ;  
  return integrate_s1s2 ( key , std::cref ( f2 ) , s , d  ) ;  
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
( const std::size_t                        tag , 
  Ostap::Math::DalitzIntegrator::function2 f2  ,
  const double                             s  ,
  const Ostap::Kinematics::Dalitz0&        d   )
{
  if ( s <= d.sqsumm () ) { return 0 ; }
  //
  const double M = std::sqrt ( s  ) ;
  //
  function2  fun = [s,M,&d,&f2] ( const double x1  , const double x2 ) -> double
                   {
                     double s1, s2 ;                     
                     std::tie ( s1 , s2 ) = d.x2s ( s  , x1  , x2 ) ;
                     const double J = d.J ( s  , s1  , s2 ) ;
                     return J <= 0 ? 0.0 : f2 ( s1 , s2 ) * J ;
                   };  
  //
  /// get integrator
  const auto F = s_DI.make_function
    ( &fun , -1  , 1 , d.s2_min () , d.s2_max ( M ) ) ;
  //
  int    ierror  =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  const std::size_t key = std::hash_combine ( tag ,  s , d.tag() ) ;
  std::tie ( ierror , result , error ) = s_DI.cubature_with_cache 
    ( key , &F , 20000 , s_PRECISION , s_PRECISION , s_MESSAGE , __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
/*  evaluate the integral over \f$s\f$ , \f$s_1\f$ variables 
 *  \f[ \int\int f( s, s_1,s_2) ds ds_1 \f]  
 *  @param tag tag that indicate theuniquness of function
 *  @param f3 the function \f$  f(s,s_1,s_2) \f$
 *  @param s2 fixed value of s2 
 *  @param smax upper-edge fo inntegrato ove \f$ s \f$
 *  @param d  helper Dalitz-object 
 *  @return integral over \f$ s, s_1\f$
 */
double Ostap::Math::DalitzIntegrator::integrate_ss1
( const std::size_t                tag ,
  function3                         f3   ,
  const double                      s2   ,
  const double                      smax ,
  const Ostap::Kinematics::Dalitz0& d    )
{
  if ( s2   <= d.s1_min () ||
       smax <= d.sqsumm () || 
       s2   >= d.s2_max ( smax ) ) { return 0 ; }
  //
  const double smin = d.sqsumm() + s2  - d.s2_min() ;
  if ( smax  <= smin ) { return 0 ; }
  //
  function2  fun = [&d,s2,&f3] ( const double y1  , const double y2 ) -> double
                   {
                     double s, s1 ;
                     std::tie ( s , s1 ) = d.y2s ( s2 , y1 , y2 ) ;
                     const double J = d.J ( s  , s1  , s2 ) ;
                     return J <= 0 ? 0.0 : f3 ( s , s1 , s2 ) * J ;
                   };
  
  /// get the  integrator
  const auto F = s_DI.make_function ( &fun , smin , smax  , -1  , 1 ) ;
  //
  int    ierror  =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  //
  const std::size_t key = std::hash_combine ( tag , s2 , smax , d.tag() ) ;
  std::tie ( ierror , result , error ) = s_DI.cubature_with_cache 
    ( key , &F , 20000 , s_PRECISION , s_PRECISION , s_MESSAGE , __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
/*  evaluate the integral over \f$e_2\f$ , \f$e_3\f$ variables 
 *  \f[ \int\int de_2 de_3 f ( M ,e_2,e_3) = 
 *  \int_{e_2^{min}}^{e_2^{max}} de_2 
 *  \int_{e_3^{min}(e_2)}^{e_3^{max}(e_2)} de_3 f(M,e_2,e_3)  \f] 
 *  @param  f the function \f$ f(M,e_2,e_3)\f$ 
 *  @param d  helper Dalitz-object 
 *  @return the integral over Dalitz plot
 */
// ============================================================================
double Ostap::Math::DalitzIntegrator::integrate_e2e3
( Ostap::Math::DalitzIntegrator::function3 f3 ,
  const Ostap::Kinematics::Dalitz&         d  )
{
  const double M = d.M()  ;
  auto      f1 = std::cref ( f3 )  ;
  function2 f2 = std::bind ( f1 , d.M() , std::placeholders::_1 , std::placeholders::_2 ) ;
  return integrate_e2e3 ( std::cref ( f2 ) , d ) ;
}
// ============================================================================
/* evaluate the integral over \f$e_2\f$ , \f$e_3\f$ variables 
 *  \f[ \int\int de_2 de_3 f(e_2,e_3) = 
 *  \int_{e_2^{min}}^{e_2^{max}} de_2 
 *  \int_{e_3^{min}(e_2)}^{e_3^{max}(e_2)} de_3 f(e_2,e_3)  \f] 
 *  @param  f the function \f$ f(e_2,e_3)\f$ 
 *  @param d  helper Dalitz-object 
 *  @return the integral over dalitz plot
 */
// ============================================================================
double Ostap::Math::DalitzIntegrator::integrate_e2e3
( Ostap::Math::DalitzIntegrator::function2 f2 ,
  const Ostap::Kinematics::Dalitz&         d  ) 
{
  //
  const double M = d.M() ;
  auto         ff  = std::cref ( f2 )  ;
  const double J   = 0.25 / ( M * M ) ; // jacobian 
  function2    fun = [&d,&ff,J] ( const double s1 , const double s2 ) -> double
                     {
                       const double e2 = d.E2 ( s1 , s2 ) ;
                       const double e3 = d.E3 ( s1 , s2 ) ;
                       if ( e2 <= 0 || e3 <= 0 || !d.inside ( s1 , s2 ) ) { return 0 ; }
                       return ff ( e2 , e3 ) * J ;
                     } ;
  //
  return integrate_s1s2 ( std::cref ( fun ) , M * M , d ) ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================

