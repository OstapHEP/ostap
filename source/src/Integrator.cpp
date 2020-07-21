// =============================================================================
// Include files 
// =============================================================================
// STD&STL
// =============================================================================
#include <functional>
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
#include "local_math.h"
// =============================================================================
/** @file 
 *  Implementation file for class Ostap::Math::Integrator
 *  @date 2018-10-11 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// =============================================================================
// constructor with integration workspace size 
// =============================================================================
Ostap::Math::Integrator::Integrator 
( const std::size_t size ) 
  : m_workspace ( size )
{}
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
  const Ostap::Math::WorkSpace&      ws   ) 
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
  return result ;
}
// ============================================================================
/*  calculate the integral 
 *  \f[ r = \int_{-\infty}^{+\infty} f_1(x) dx \f]
 *  @param f1 the function 
 *  @param integration workspace 
 *  @return the value of the integral 
 */
// ============================================================================
double Ostap::Math::Integrator::integrate_infinity
( Ostap::Math::Integrator::function1 f1   , 
  const Ostap::Math::WorkSpace&      ws   ) 
{
  static const Ostap::Math::GSL::Integrator1D<function1> integrator {} ;
  auto F = integrator.make_function( &f1 ) ;
  //
  int    ierror ;
  double result ;
  double error  ;
  static const char s_message[] = "Ostap::Math::Integrator/integrate_infinity" ;
  std::tie ( ierror, result , error ) = 
    integrator.gaqi_integrate 
    ( &F                , 
      workspace ( ws )  ,   // workspace 
      s_PRECISION       ,   // absolute precision 
      s_PRECISION       ,   // relative precision 
      -1                ,   // limit 
      s_message         ,   // reason of failure 
      __FILE__          ,   // the file 
      __LINE__          ) ; // the line 
  //
  return result ;
}
// ============================================================================
/* calculate the integral 
 *  \f[ r = \int_{x_{min}}^{+\infty} f_1(x) dx \f]
 *  @param f1 the function 
 *  @param xmin lower integration edge 
 *  @param integration workspace 
 *  @return the value of the integral 
 */
// ============================================================================
double Ostap::Math::Integrator::integrate_to_infinity
( Ostap::Math::Integrator::function1 f1   , 
  const double                       xmin , 
  const Ostap::Math::WorkSpace&      ws   ) 
{
  static const Ostap::Math::GSL::Integrator1D<function1> integrator {} ;
  auto F = integrator.make_function( &f1 ) ;
  //
  int    ierror ;
  double result ;
  double error  ;
  static const char s_message[] = "Ostap::Math::Integrator/integrate_to_infinity" ;
  std::tie ( ierror, result , error ) = 
    integrator.gaqiu_integrate 
    ( &F                , 
      xmin              ,   // lower integration edge  
      workspace ( ws )  ,   // workspace 
      s_PRECISION       ,   // absolute precision 
      s_PRECISION       ,   // relative precision 
      -1                ,   // limit 
      s_message         ,   // reason of failure 
      __FILE__          ,   // the file 
      __LINE__          ) ; // the line 
  //
  return result ;
}
// ============================================================================
/* calculate the integral 
 *  \f[ r = \int_{-\infty}^{x_{max}} f_1(x) dx \f]
 *  @param f1 the function 
 *  @param xmax upper  integration edge
 *  @param integration workspace 
 *  @return the value of the integral 
 */
// ============================================================================
double Ostap::Math::Integrator::integrate_from_infinity
( Ostap::Math::Integrator::function1 f1   , 
  const double                       xmax , 
  const Ostap::Math::WorkSpace&      ws   ) 
{
  static const Ostap::Math::GSL::Integrator1D<function1> integrator {} ;
  auto F = integrator.make_function( &f1 ) ;
  //
  int    ierror ;
  double result ;
  double error  ;
  static const char s_message[] = "Ostap::Math::Integrator/integrate_to_infinity" ;
  std::tie ( ierror, result , error ) = 
    integrator.gaqil_integrate 
    ( &F                , 
      xmax              ,   // lower integration edge  
      workspace ( ws )  ,   // workspace 
      s_PRECISION       ,   // absolute precision 
      s_PRECISION       ,   // relative precision 
      -1                ,   // limit 
      s_message         ,   // reason of failure 
      __FILE__          ,   // the file 
      __LINE__          ) ; // the line 
  //
  return result ;
}
// ============================================================================
/** get Cauchy principal value integral 
 *  \f[ g(c) =  \int_{x_{min}}^{x_{max}}\frac{f(x)}{x-c}dx \f]
 *  @param f1 the function 
 *  @param c  the parameter 
 *  @param xmin lower integration edge 
 *  @param xmax upper  integration edge 
 *  @param integration workspace 
 */
// ============================================================================
double Ostap::Math::Integrator::cauchy_pv 
( Ostap::Math::Integrator::function1 f1   , 
  const double                       c    , 
  const double                       xmin , 
  const double                       xmax , 
  const Ostap::Math::WorkSpace&      ws   ) 
{
  if      ( s_equal ( xmax , xmin ) ) {return 0 ; }  
  else if (           xmax < xmin   ) 
  { return -1 * cauchy_pv ( std::cref ( f1 ) , c , xmax , xmin , ws ) ; }
  //
  // regular integration  
  if ( c < xmin || xmax < c ) 
  {
    auto f2 = std::cref ( f1 ) ;
    auto ff = [f2,c]  ( const double  x ) -> double 
      { return f2 ( x ) / ( x - c ) ; } ;
    //
    return integrate ( std::cref ( ff ) , xmin , xmax , ws ) ;
  }
  else if ( s_equal ( c , xmin ) ) 
  {
    auto f2 = std::cref ( f1 ) ;
    auto ff = [f2,xmin]  ( const double  x ) -> double 
      { return f2 ( x ) / ( x - xmin ) ; } ;
    //
    const double xc = xmin + 0.05 * ( xmax - xmin ) ;
    static const Ostap::Math::GSL::Integrator1D<function1> integrator {} ;
    auto F = integrator.make_function( &f1 ) ;
    //
    int    ierror ;
    double result ;
    double error  ;
    static const char s_message[] = "Ostap::Math::Integrator/integrate_to_infinity" ;
    std::tie ( ierror, result , error ) = 
      integrator.gaqp_integrate 
      ( &F                    , 
        xmin                  ,   // lower integration edge  
        xc                    ,   // high integration edge  
        std::vector<double>() ,   // other singular points 
        workspace ( ws )      ,   // workspace 
        s_PRECISION           ,   // absolute precision 
        s_PRECISION           ,   // relative precision 
        -1                    ,   // limit 
        s_message             ,   // reason of failure 
        __FILE__              ,   // the file 
        __LINE__              ) ; // the line 
    //
    return result + integrate ( std::cref ( ff ) , xc , xmax , ws ) ;
  }
  else if ( s_equal ( c , xmax ) ) 
  {
    auto f2 = std::cref ( f1 ) ;
    auto ff = [f2,xmax]  ( const double  x ) -> double 
      { return f2 ( x ) / ( x - xmax ) ; } ;
    //
    const double xc = xmax - 0.05 * ( xmax - xmin ) ;
    static const Ostap::Math::GSL::Integrator1D<function1> integrator {} ;
    auto F = integrator.make_function( &f1 ) ;
    //
    int    ierror ;
    double result ;
    double error  ;
    static const char s_message[] = "Ostap::Math::Integrator/integrate_to_infinity" ;
    std::tie ( ierror, result , error ) = 
      integrator.gaqp_integrate 
      ( &F                    , 
        xc                    ,   // lower integration edge  
        xmax                  ,   // high integration edge  
        std::vector<double>() ,   // other singular points 
        workspace ( ws )      ,   // workspace 
        s_PRECISION           ,   // absolute precision 
        s_PRECISION           ,   // relative precision 
        -1                    ,   // limit 
        s_message             ,   // reason of failure 
        __FILE__              ,   // the file 
        __LINE__              ) ; // the line 
    //
    return result + integrate ( std::cref ( ff ) , xmin , xc , ws ) ;
  }
  //
  // regular Cauchy integral
  //
  const double dx = std::min ( std::abs ( c - xmin ) , std::abs ( c - xmax ) ) / 2 ;
  //
  auto f2 = std::cref ( f1 ) ;
  auto ff = [f2,c]  ( const double  x ) -> double 
    { return f2 ( x ) / ( x - c ) ; } ;
  //
  static const Ostap::Math::GSL::Integrator1D<function1> integrator {} ;
  auto F = integrator.make_function( &f1 ) ;
  //
  int    ierror ;
  double result ;
  double error  ;
  static const char s_message[] = "Ostap::Math::Integrator/integrate_cauchy_pv" ;
  std::tie ( ierror, result , error ) = 
    integrator.gawc_integrate 
    ( &F                , 
      c - dx , c + dx   ,   // low and high integrtaion edges 
      c                 ,   // Cauchy's point 
      workspace ( ws )  ,   // workspace 
      s_PRECISION       ,   // absolute precision 
      s_PRECISION       ,   // relative precision 
      -1                ,   // limit 
      s_message         ,   // reason of failure 
      __FILE__          ,   // the file 
      __LINE__          ) ; // the line 
  //
  return result + 
    integrate ( std::cref ( ff ) , xmin   , c - dx , ws ) + 
    integrate ( std::cref ( ff ) , c + dx , xmax   , ws ) ;
  //
}
// ============================================================================
/*  get Cauchy principal value integral 
 *  \f[ g(c) =  \int_{x_{min}}^{+\infty}\frac{f(x)}{x-c}dx \f]
 *  @param f1 the function 
 *  @param c  the parameter 
 *  @param xmin lower integration edge 
 *  @param integration workspace 
 */
// ============================================================================
double Ostap::Math::Integrator::cauchy_pv_to_infinity 
( Ostap::Math::Integrator::function1 f1    , 
  const double                       c     , 
  const double                       xmin  , 
  const Ostap::Math::WorkSpace&      ws    ) 
{
  //
  if ( c < xmin ) 
  {
    auto f2 = std::cref ( f1 ) ;
    auto ff = [f2,c]  ( const double  x ) -> double 
      { return f2 ( x ) / ( x - c ) ; } ;
    //
    return integrate_to_infinity ( std::cref ( ff ) , xmin , ws ) ;
  }
  //
  double xx = c + ( c - xmin ) / 2 ;
  if ( s_equal ( xx , xmin ) ) { xx += 1 ; }
  //
  auto f2 = std::cref ( f1 ) ;
  auto ff = [f2,c]  ( const double  x ) -> double 
    { return f2 ( x ) / ( x - c ) ; } ;
  //
  return 
    cauchy_pv              ( std::cref ( f1 ) , c , xmin , xx , ws ) + 
    integrate_to_infinity  ( std::cref ( ff ) ,            xx , ws ) ;
}
// ============================================================================
/*  get Cauchy principal value integral 
 *  \f[ g(c) =  \int_{-\infty}^{+x_{max}}\frac{f(x)}{x-c}dx \f]
 *  @param f1 the function 
 *  @param c  the parameter 
 *  @param xmax upper integration edge 
 *  @param integration workspace 
 */
// ============================================================================
double Ostap::Math::Integrator::cauchy_pv_from_infinity 
( Ostap::Math::Integrator::function1 f1    , 
  const double                       c     , 
  const double                       xmax  , 
  const Ostap::Math::WorkSpace&      ws    ) 
{
  //
  if ( c > xmax ) 
  {
    auto f2 = std::cref ( f1 ) ;
    auto ff = [f2,c]  ( const double  x ) -> double 
      { return f2 ( x ) / ( x - c ) ; } ;
    //
    return integrate_from_infinity ( std::cref ( ff ) , xmax , ws ) ;
  }
  //
  double xx =  c - ( xmax - c ) / 2 ;
  if ( s_equal ( xx , xmax ) ) { xx -= 1 ; }
  //
  auto f2 = std::cref ( f1 ) ;
  auto ff = [f2,c]  ( const double  x ) -> double 
    { return f2 ( x ) / ( x - c ) ; } ;
  //
  return 
    integrate_from_infinity  ( std::cref ( ff ) ,            xx , ws ) +
    cauchy_pv                ( std::cref ( f1 ) , c , xx , xmax , ws ) ;
}
// ============================================================================
/*  Kramers-Kronig dispersion relation with n-subtractions 
 *  \f[ g(s) = \frac{s^n}{\pi} 
 *     \mathcal{P} \int_{x_{min}}^{+\infty} \frac{g(x)}{x^n(x-s)}dx \f] 
 *  @param f1    the function 
 *  @param s     s-parameter 
 *  @param xmin  lower integration range 
 *  @param n     number of subtracion
 *  @return value of the dispersion integral 
 */
// ============================================================================
double Ostap::Math::Integrator::kramers_kronig
( Ostap::Math::Integrator::function1 f1   , 
  const double                       s    , 
  const double                       xmin , 
  const unsigned short               n    , 
  const Ostap::Math::WorkSpace&      ws   ) 
{
  // 
  if ( 0 < n ) 
  {
    auto f2 = std::cref ( f1 ) ;
    auto ff = [f2,n] ( const double x ) -> double 
      { return f2 ( x ) / std::pow ( x , n ) ; } ;
    //
    return 
      std::pow ( s , n ) * 
      kramers_kronig ( std::cref ( ff ) , s , xmin , 0 , ws ) ;
  }
  // no subtractions 
  return cauchy_pv_to_infinity ( std::cref ( f1 ) , s , xmin , ws ) / M_PI ;
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
  const Ostap::Math::WorkSpace&      ws   ) 
{
  auto f2_ = std::cref ( f2 ) ;
  auto f1  = std::bind ( f2_ , std::placeholders::_1 , y ) ;
  return integrate ( std::cref ( f1 ) , xmin , xmax , ws ) ;
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
  const Ostap::Math::WorkSpace&      ws   )
{
  auto f2_ = std::cref ( f2 ) ;
  auto f1  = std::bind ( f2_ , x , std::placeholders::_1 ) ;
  return integrate ( std::cref ( f1 ) , ymin , ymax , ws ) ;
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
  const Ostap::Math::WorkSpace&     /* ws */ )  
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
  const double                       ymax ) 
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
double Ostap::Math::Integrator::integrate
( const std::size_t                  tag  , 
  Ostap::Math::Integrator::function1 f1   , 
  const double                       xmin , 
  const double                       xmax , 
  const Ostap::Math::WorkSpace&      ws   ) 
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
( const std::size_t                  tag  , 
  Ostap::Math::Integrator::function2 f2   , 
  const double                       y    , 
  const double                       xmin ,
  const double                       xmax ,
  const Ostap::Math::WorkSpace&      ws   ) 
{
  auto f2_ = std::cref ( f2 ) ;
  auto f1  = std::bind ( f2_ , std::placeholders::_1 , y ) ;
  return integrate ( std::hash_combine ( tag , 'X' , y )  , std::cref ( f1 ) , xmin , xmax , ws ) ;
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
( const std::size_t                  tag  , 
  Ostap::Math::Integrator::function2 f2   , 
  const double                       x    , 
  const double                       ymin ,
  const double                       ymax ,
  const Ostap::Math::WorkSpace&      ws   ) 
{
  auto f2_ = std::cref ( f2 ) ;
  auto f1  = std::bind ( f2_ , x , std::placeholders::_1 ) ;
  return integrate ( std::hash_combine ( tag , 'Y' , x ) , std::cref ( f1 ), ymin , ymax , ws ) ;
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
( const std::size_t                  tag  , 
  Ostap::Math::Integrator::function2 f2   , 
  const double                       xmin , 
  const double                       xmax ,
  const double                       ymin , 
  const double                       ymax , 
  const Ostap::Math::WorkSpace&     /* ws */ ) 
{ return integrate ( tag , std::cref ( f2 ) , xmin , xmax , ymin , ymax ) ; }
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
( const std::size_t                  tag  , 
  Ostap::Math::Integrator::function2 f2   , 
  const double                       xmin , 
  const double                       xmax ,
  const double                       ymin , 
  const double                       ymax ) 
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
