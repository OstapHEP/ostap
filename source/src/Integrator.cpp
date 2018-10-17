// =============================================================================
// Include files 
// =============================================================================
// STD&STL
// =============================================================================
#include <functional>
#include <iostream>
// =============================================================================
// Ostap
// =============================================================================
#include "Ostap/Integrator.h"
#include "Ostap/Workspace.h"
// =============================================================================
// Local
// =============================================================================
#include "Integrator1D.h"
#include "Integrator2D.h"
#include "local_hash.h"
// =============================================================================
/** @file 
 *  Implementation file for class Ostap::Math::Integrator
 *  @date 2018-10-11 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// =============================================================================
/*  calculate the integral 
 *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
 *  @param f1 the function 
 *  @param xmin lower integration edge 
 *  @param xmax uppr  integration edge
 *  @return the value of the integral 
 */
// =============================================================================
double Ostap::Math::Integrator::integrate
( Ostap::Math::Integrator::function1 f1   , 
  const double                       xmin , 
  const double                       xmax , 
  const Ostap::Math::WorkSpace&      ws   ) const 
{
  static const Ostap::Math::GSL::Integrator1D<function1> integrator {} ;
  auto F = integrator.make_function( &f1 ) ;
  //
  int    ierror ;
  double result ;
  double error  ;
  static const char s_message[] = "Ostap::Math::Integrator/integrate(1D)" ;
  std::tie ( ierror, result , error ) = 
    integrator.gaq_integrate 
    ( &F                , 
      xmin              ,   // lower integration edge  
      xmax              ,   // upper integration edge
      workspace ( ws )  ,   // workspace 
      s_PRECISION       ,   // absolute precision 
      s_PRECISION       ,   // relative precision 
      -1                ,   // limit 
      s_message         ,   // reason of failure 
      __FILE__          ,   // the file 
      __LINE__          ) ; // the line 
  //
  //
  std::cout << " I " << ";   result=" << result << std::endl ;
  //
  return result ;
}
// ============================================================================
/*  calculate the integral 
 *  \f[ r = \int_{x_{min}}^{x_{max}}f_2(x,y) dx \f]
 *  @param f2 the function 
 *  @param y parameter y
 *  @param xmin lower integration edge in x 
 *  @param xmax upper integration edge in x 
 *  @return the value of the integral 
 */
// ============================================================================
double Ostap::Math::Integrator::integrateX
( Ostap::Math::Integrator::function2 f2   , 
  const double                       y    , 
  const double                       xmin ,
  const double                       xmax ,
  const Ostap::Math::WorkSpace&      ws   ) const 
{
  auto f2_ = std::cref ( f2 ) ;
  auto f1  = std::bind ( f2_ , std::placeholders::_1 , y ) ;
  return integrate ( f1 , xmin , xmax , ws ) ;
} 
// ============================================================================
/** calculate the integral 
 *  \f[ r = \int_{y_{min}}^{y_{max}}f_2(x,y) dy \f]
 *  @param f2 the function 
 *  @param x parameter x
 *  @param ymin lower integration edge in y 
 *  @param ymax upper integration edge in y 
 *  @return the value of the integral 
 */
// ============================================================================
double Ostap::Math::Integrator::integrateY
( Ostap::Math::Integrator::function2 f2   , 
  const double                       x    , 
  const double                       ymin ,
  const double                       ymax ,
  const Ostap::Math::WorkSpace&      ws   ) const 
{
  auto f2_ = std::cref ( f2 ) ;
  auto f1  = std::bind ( f2_ , x , std::placeholders::_1 ) ;
  return integrate ( f1 , ymin , ymax , ws ) ;
}
// ==========================================================================
/** calculate the integral 
 *  \f[ r = \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}}f_2(x,y) dx dy \f]
 *  @param f2 the function 
 *  @param xmin lower integration edge in x 
 *  @param xmax upper integration edge in x 
 *  @param ymin lower integration edge in y 
 *  @param ymax upper integration edge in y 
 *  @param integration workspace (not used)
 *  @return the value of the integral 
 */
// ==========================================================================
double Ostap::Math::Integrator::integrate
( Ostap::Math::Integrator::function2 f2   , 
  const double                       xmin , 
  const double                       xmax ,
  const double                       ymin , 
  const double                       ymax , 
  const Ostap::Math::WorkSpace&     /* ws */ ) const 
{ return integrate ( std::cref ( f2 ) , xmin , xmax , ymin , ymax ) ; }
// ==========================================================================
/** calculate the integral 
 *  \f[ r = \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}}f_2(x,y) dx dy \f]
 *  @param f2 the function 
 *  @param xmin lower integration edge in x 
 *  @param xmax upper integration edge in x 
 *  @param ymin lower integration edge in y 
 *  @param ymax upper integration edge in y 
 *  @return the value of the integral 
 */
// ==========================================================================
double Ostap::Math::Integrator::integrate
( Ostap::Math::Integrator::function2 f2   , 
  const double                       xmin , 
  const double                       xmax ,
  const double                       ymin , 
  const double                       ymax ) const 
{ 
  //
  static const Ostap::Math::GSL::Integrator2D<function2> s_cubature{} ;
  static const char s_message[] = "Ostap::Math::Integrator/integrate(2D)" ;
  const auto F = s_cubature.make_function ( &f2 , xmin , xmax , ymin , ymax ) ;
  int     ierror =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  std::tie ( ierror , result , error ) = s_cubature.cubature 
    ( &F          ,   // the function  
      100000      ,   // limits  
      s_PRECISION ,   // absolute precision 
      s_PRECISION ,   // relative precision 
      //10000      ,   // limits  
      //1.e-5      , 
      //1.e-5      , 
      //
      s_message   ,   // message 
      __FILE__    ,   // the file name 
      __LINE__    ) ; // the line number 
  //
  return result ;
}
// =============================================================================


// =============================================================================
// Integration with cache 
// =============================================================================


// =============================================================================
/*  calculate the integral 
 *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
 *  @param f1 the function 
 *  @param xmin lower integration edge 
 *  @param xmax uppr  integration edge
 *  @return the value of the integral 
 */
// =============================================================================
double Ostap::Math::Integrator::integrate_with_cache
( const std::size_t                  tag  , 
  Ostap::Math::Integrator::function1 f1   , 
  const double                       xmin , 
  const double                       xmax , 
  const Ostap::Math::WorkSpace&      ws   ) const 
{
  static const Ostap::Math::GSL::Integrator1D<function1> integrator {} ;
  auto F = integrator.make_function( &f1 ) ;
  //
  int    ierror ;
  double result ;
  double error  ;
  static const char s_message[] = "Ostap::Math::Integrator/integrate(1Dc)" ;
  std::tie ( ierror, result , error ) = 
    integrator.gaq_integrate_with_cache   
    ( tag               , 
      &F                , 
      xmin              ,   // lower integration edge  
      xmax              ,   // upper integration edge
      workspace ( ws )  ,   // workspace 
      s_PRECISION       ,   // absolute precision 
      s_PRECISION       ,   // relative precision 
      -1                ,   // limit 
      s_message         ,   // reason of failure 
      __FILE__          ,   // the file 
      __LINE__          ) ; // the line 
  //
  std::cout << " I/tag= " << tag << ";   result=" << result << std::endl ;
  //
  return result ;
}
// ============================================================================
/*  calculate the integral 
 *  \f[ r = \int_{x_{min}}^{x_{max}}f_2(x,y) dx \f]
 *  @param f2 the function 
 *  @param y parameter y
 *  @param xmin lower integration edge in x 
 *  @param xmax upper integration edge in x 
 *  @return the value of the integral 
 */
// ============================================================================
double Ostap::Math::Integrator::integrateX_with_cache
( const std::size_t                  tag  , 
  Ostap::Math::Integrator::function2 f2   , 
  const double                       y    , 
  const double                       xmin ,
  const double                       xmax ,
  const Ostap::Math::WorkSpace&      ws   ) const 
{
  auto f2_ = std::cref ( f2 ) ;
  auto f1  = std::bind ( f2_ , std::placeholders::_1 , y ) ;
  return integrate_with_cache ( std::hash_combine ( tag , 'X' , y )  , f1 , xmin , xmax , ws ) ;
} 
// ============================================================================
/** calculate the integral 
 *  \f[ r = \int_{y_{min}}^{y_{max}}f_2(x,y) dy \f]
 *  @param f2 the function 
 *  @param x parameter x
 *  @param ymin lower integration edge in y 
 *  @param ymax upper integration edge in y 
 *  @return the value of the integral 
 */
// ============================================================================
double Ostap::Math::Integrator::integrateY_with_cache
( const std::size_t                  tag  , 
  Ostap::Math::Integrator::function2 f2   , 
  const double                       x    , 
  const double                       ymin ,
  const double                       ymax ,
  const Ostap::Math::WorkSpace&      ws   ) const 
{
  auto f2_ = std::cref ( f2 ) ;
  auto f1  = std::bind ( f2_ , x , std::placeholders::_1 ) ;
  return integrate_with_cache ( std::hash_combine ( tag , 'Y' , x ) , f1 , ymin , ymax , ws ) ;
}
// ==========================================================================
/** calculate the integral 
 *  \f[ r = \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}}f_2(x,y) dx dy \f]
 *  @param f2 the function 
 *  @param xmin lower integration edge in x 
 *  @param xmax upper integration edge in x 
 *  @param ymin lower integration edge in y 
 *  @param ymax upper integration edge in y 
 *  @param integration workspace (not used)
 *  @return the value of the integral 
 */
// ==========================================================================
double Ostap::Math::Integrator::integrate_with_cache
( const std::size_t                  tag  , 
  Ostap::Math::Integrator::function2 f2   , 
  const double                       xmin , 
  const double                       xmax ,
  const double                       ymin , 
  const double                       ymax , 
  const Ostap::Math::WorkSpace&     /* ws */ ) const 
{ return integrate_with_cache ( tag , std::cref ( f2 ) , xmin , xmax , ymin , ymax ) ; }
// ==========================================================================
/** calculate the integral 
 *  \f[ r = \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}}f_2(x,y) dx dy \f]
 *  @param f2 the function 
 *  @param xmin lower integration edge in x 
 *  @param xmax upper integration edge in x 
 *  @param ymin lower integration edge in y 
 *  @param ymax upper integration edge in y 
 *  @return the value of the integral 
 */
// ==========================================================================
double Ostap::Math::Integrator::integrate_with_cache
( const std::size_t                  tag  , 
  Ostap::Math::Integrator::function2 f2   , 
  const double                       xmin , 
  const double                       xmax ,
  const double                       ymin , 
  const double                       ymax ) const 
{ 
  //
  static const Ostap::Math::GSL::Integrator2D<function2> s_cubature{} ;
  static const char s_message[] = "Ostap::Math::Integrator/integrate(2Dc)" ;
  const auto F = s_cubature.make_function ( &f2 , xmin , xmax , ymin , ymax ) ;
  int     ierror =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  std::tie ( ierror , result , error ) = s_cubature.cubature_with_cache 
    ( tag         , 
      &F          ,   // the function  
      100000      ,   // limits  
      s_PRECISION ,   // absolute precision 
      s_PRECISION ,   // relative precision 
      // 20000      ,   // limits  
      // 1.e-5      , 
      // 1.e-5      , 
      //
      s_message   ,   // message 
      __FILE__    ,   // the file name 
      __LINE__    ) ; // the line number 
  //
  return result ;
}
// =============================================================================
//                                                                       The END 
// =============================================================================
