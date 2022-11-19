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
namespace 
{
  // ===========================================================================
  template <class FUNCTION>
  inline double fun_scale ( const FUNCTION&      fun     , 
                            const double         xmin    , 
                            const double         xmax    ,
                            const unsigned short rescale )
  { 
    if ( 0 == rescale ) { return 1.0 ; }
    const double dx = ( xmax - xmin ) / ( rescale + 1 ) ;
    double scale = 0 ;
    for ( unsigned short i = 1 ; i <= rescale ; ++i ) { scale += fun ( xmin + i * dx ) ; }
    return scale * ( xmin - xmax ) / rescale ;
  } 
  // ===========================================================================
}
// =============================================================================
// constructor with integration workspace size 
// =============================================================================
Ostap::Math::Integrator::Integrator 
( const std::size_t size ) 
  : m_workspace ( size )
  , m_abs_precision_gaq   ( s_APRECISION_GAQ   ) 
  , m_rel_precision_gaq   ( s_RPRECISION_GAQI  ) 
  , m_abs_precision_gaqi  ( s_APRECISION_GAQI  ) 
  , m_rel_precision_gaqi  ( s_RPRECISION_GAQI  ) 
  , m_abs_precision_gaqiu ( s_APRECISION_GAQIU ) 
  , m_rel_precision_gaqiu ( s_RPRECISION_GAQIU ) 
  , m_abs_precision_gaqil ( s_APRECISION_GAQIL ) 
  , m_rel_precision_gaqil ( s_RPRECISION_GAQIL ) 
  , m_abs_precision_gaqp  ( s_APRECISION_GAQP  ) 
  , m_rel_precision_gaqp  ( s_RPRECISION_GAQP  )
  , m_abs_precision_qawc  ( s_APRECISION_QAWC  ) 
  , m_rel_precision_qawc  ( s_RPRECISION_QAWC  )
  , m_abs_precision_cpv   ( s_APRECISION_QAWC  ) 
  , m_rel_precision_cpv   ( s_RPRECISION_QAWC  ) 
  , m_abs_precision_cpvi  ( s_APRECISION_QAWC  ) 
  , m_rel_precision_cpvi  ( s_RPRECISION_QAWC  ) 
  , m_abs_precision_kk    ( s_APRECISION_QAWC  ) 
  , m_rel_precision_kk    ( s_RPRECISION_QAWC  )
  , m_abs_precision_cube  ( s_APRECISION_CUBE  ) 
  , m_rel_precision_cube  ( s_RPRECISION_CUBE  ) 
{}
// =============================================================================
// set absolute/relative precision for GAG
// =============================================================================
void Ostap::Math::Integrator::set_precision_gaq   
( const double aprec ,
  const double rprec ) 
{
  m_abs_precision_gaq = 0 < aprec ? aprec : s_APRECISION_GAQ ;
  m_rel_precision_gaq = 0 < rprec ? rprec : s_RPRECISION_GAQ ;
}
// =============================================================================
// set absolute/relative precision for GAGI
// =============================================================================
void Ostap::Math::Integrator::set_precision_gaqi   
( const double aprec ,
  const double rprec ) 
{
  m_abs_precision_gaqi = 0 < aprec ? aprec : s_APRECISION_GAQI ;
  m_rel_precision_gaqi = 0 < rprec ? rprec : s_RPRECISION_GAQI ;
}
// =============================================================================
// set absolute/relative precision for GAGIL
// =============================================================================
void Ostap::Math::Integrator::set_precision_gaqil   
( const double aprec ,
  const double rprec ) 
{
  m_abs_precision_gaqil = 0 < aprec ? aprec : s_APRECISION_GAQIL ;
  m_rel_precision_gaqil = 0 < rprec ? rprec : s_RPRECISION_GAQIL ;
}
// =============================================================================
// set absolute/relative precision for GAGIU
// =============================================================================
void Ostap::Math::Integrator::set_precision_gaqiu
( const double aprec ,
  const double rprec ) 
{
  m_abs_precision_gaqiu = 0 < aprec ? aprec : s_APRECISION_GAQIU ;
  m_rel_precision_gaqiu = 0 < rprec ? rprec : s_RPRECISION_GAQIU ;
}
// =============================================================================
// set absolute/relative precision for GAQP
// =============================================================================
void Ostap::Math::Integrator::set_precision_gaqp
( const double aprec ,
  const double rprec ) 
{
  m_abs_precision_gaqp  = 0 < aprec ? aprec : s_APRECISION_GAQP ;
  m_rel_precision_gaqp  = 0 < rprec ? rprec : s_RPRECISION_GAQP ;
}
// =============================================================================
// set absolute/relative precision for QAWC
// =============================================================================
void Ostap::Math::Integrator::set_precision_qawc
( const double aprec ,
  const double rprec ) 
{
  m_abs_precision_qawc  = 0 < aprec ? aprec : s_APRECISION_QAWC ;
  m_rel_precision_qawc  = 0 < rprec ? rprec : s_RPRECISION_QAWC ;
}
// =============================================================================
// set absolute/relative precision for Cauchy PV integration 
// =============================================================================
void Ostap::Math::Integrator::set_precision_cpv
( const double aprec ,
  const double rprec ) 
{
  m_abs_precision_cpv  = 0 < aprec ? aprec : s_APRECISION_QAWC ;
  m_rel_precision_cpv  = 0 < rprec ? rprec : s_RPRECISION_QAWC ;
}
// =============================================================================
// set absolute/relative precision for Cauchy PV/inf integration 
// =============================================================================
void Ostap::Math::Integrator::set_precision_cpvi
( const double aprec ,
  const double rprec ) 
{
  m_abs_precision_cpvi = 0 < aprec ? aprec : s_APRECISION_QAWC ;
  m_rel_precision_cpvi = 0 < rprec ? rprec : s_RPRECISION_QAWC ;
}
// =============================================================================
// set absolute/relative precision for Kramers-Kronig integration
// =============================================================================
void Ostap::Math::Integrator::set_precision_kk
( const double aprec ,
  const double rprec ) 
{
  m_abs_precision_kk  = 0 < aprec ? aprec : s_APRECISION_QAWC ;
  m_rel_precision_kk  = 0 < rprec ? rprec : s_RPRECISION_QAWC ;
}
// =============================================================================
// set absolute/relative precision for cubature 
// =============================================================================
void Ostap::Math::Integrator::set_precision_cube
( const double aprec ,
  const double rprec ) 
{
  m_abs_precision_cube  = 0 < aprec ? aprec : s_APRECISION ;
  m_rel_precision_cube  = 0 < rprec ? rprec : s_RPRECISION ;
}
// =============================================================================
/*  calculate the integral 
 *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
 *  @param f1 the function 
 *  @param xmin lower integration edge 
 *  @param xmax upper integration edge
 *  @return value of the integral and the estimate of the uncertainty
 */
// =============================================================================
Ostap::Math::Integrator::result
Ostap::Math::Integrator::integrate_
( Ostap::Math::Integrator::function1 f1         , 
  const double                       xmin       , 
  const double                       xmax       , 
  const Ostap::Math::WorkSpace&      ws         , 
  const std::size_t                  tag        ,  
  const unsigned short               rescale    , 
  const double                       aprecision ,
  const double                       rprecision ) 
{
  //
  if ( s_equal ( xmin , xmax ) ) { return result ( 0 , 0 )  ; }
  //
  if ( 0 < rescale ) 
  {
    const double scale = fun_scale ( f1 , xmin , xmax , rescale ) ;
    if ( !s_zero ( scale ) && std::abs ( scale ) < 0.1 || 10 < std::abs ( scale ) ) 
    {
      const double iscale = 1.0 / scale ;
      auto f2 = std::cref ( f1 ) ;
      auto ff = [f2,iscale]  ( const double  x ) -> double { return f2 ( x ) * iscale ; } ;
      //
      const std::size_t ntag = 
        0 == tag ? tag : std::hash_combine ( tag , rescale , scale , iscale ) ;
      result r = integrate_ ( std::cref ( ff ) , xmin , xmax , ws , ntag , 0 ) ;
      return result ( scale * r.first , scale * r.second ) ;  
    }
  }
  //
  static const Ostap::Math::GSL::Integrator1D<function1> integrator {} ;
  auto F = integrator.make_function( &f1 ) ;
  //
  int    ierror ;
  double value  ;
  double error  ;
  static const char s_message[] = "Ostap::Math::Integrator/integrate(1D)" ;
  //
  std::tie ( ierror, value , error ) = 
    integrator.gaq_integrate    
    ( &F                , 
      xmin              ,   // lower integration edge  
      xmax              ,   // upper integration edge
      workspace ( ws )  ,   // workspace 
      0 < aprecision ? aprecision : s_APRECISION_GAQ  ,   // absolute precision 
      0 < rprecision ? rprecision : s_RPRECISION_GAQ  ,   // relative precision 
      -1                ,   // limit 
      s_message         ,   // reason of failure 
      __FILE__          ,   // the file 
      __LINE__          ,   // the line
      GSL_INTEG_GAUSS51 ,   // the rule 
      tag               ) ; // label/tag
  //
  return result ( value , error ) ;  
}
// ============================================================================
/*  calculate the integral 
 *  \f[ r = \int_{-\infty}^{+\infty} f_1(x) dx \f]
 *  @param f1 the function 
 *  @param integration workspace 
 *  @return value of the integral and the estimate of the uncertainty
 */
// ============================================================================
Ostap::Math::Integrator::result
Ostap::Math::Integrator::integrate_infinity_
( Ostap::Math::Integrator::function1 f1         , 
  const Ostap::Math::WorkSpace&      ws         ,
  const std::size_t                  tag        ,
  const double                       aprecision ,
  const double                       rprecision )
{
  static const Ostap::Math::GSL::Integrator1D<function1> integrator {} ;
  auto F = integrator.make_function( &f1 ) ;
  //
  int    ierror ;
  double value  ;
  double error  ;
  static const char s_message[] = "Ostap::Math::Integrator/integrate_infinity" ;
  std::tie ( ierror, value , error ) = 
    integrator.gaqi_integrate 
    ( &F                , 
      workspace ( ws )  ,   // workspace 
      0 < aprecision ? aprecision : s_APRECISION_GAQI ,   // absolute precision 
      0 < rprecision ? rprecision : s_RPRECISION_GAQI ,   // relative precision 
      -1                ,   // limit 
      s_message         ,   // reason of failure 
      __FILE__          ,   // the file 
      __LINE__          ,   // the line 
      tag               ) ; // tag/label 
  //
  return result ( value , error ) ;
}
// ============================================================================
/*  calculate the integral 
 *  \f[ r = \int_{x_{min}}^{+\infty} f_1(x) dx \f]
 *  @param f1 the function 
 *  @param xmin lower integration edge 
 *  @param integration workspace 
 *  @return value of the integral and the estimate of the uncertainty
 */
// ============================================================================
Ostap::Math::Integrator::result
Ostap::Math::Integrator::integrate_to_infinity_
( Ostap::Math::Integrator::function1 f1         , 
  const double                       xmin       , 
  const Ostap::Math::WorkSpace&      ws         ,
  const std::size_t                  tag        ,
  const double                       aprecision ,
  const double                       rprecision )
{
  static const Ostap::Math::GSL::Integrator1D<function1> integrator {} ;
  auto F = integrator.make_function( &f1 ) ;
  //
  int    ierror ;
  double value  ;
  double error  ;
  static const char s_message[] = "Ostap::Math::Integrator/integrate_to_infinity" ;
  std::tie ( ierror, value , error ) = 
    integrator.gaqiu_integrate 
    ( &F                 , 
      xmin               ,   // lower integration edge  
      workspace ( ws )   ,   // workspace 
      0 < aprecision ? aprecision : s_APRECISION_GAQIU ,   // absolute precision 
      0 < rprecision ? rprecision : s_RPRECISION_GAQIU ,   // relative precision 
      -1                 ,   // limit 
      s_message          ,   // reason of failure 
      __FILE__           ,   // the file 
      __LINE__           ,   // the line 
      tag                ) ; // tag/label 
  //
  return result ( value , error ) ;
}
// ============================================================================
/* calculate the integral 
 *  \f[ r = \int_{-\infty}^{x_{max}} f_1(x) dx \f]
 *  @param f1 the function 
 *  @param xmax upper  integration edge
 *  @param integration workspace 
 *  @return value of the integral and the estimate of the uncertainty
 */
// ============================================================================
Ostap::Math::Integrator::result
Ostap::Math::Integrator::integrate_from_infinity_
( Ostap::Math::Integrator::function1 f1         , 
  const double                       xmax       , 
  const Ostap::Math::WorkSpace&      ws         ,
  const std::size_t                  tag        ,
  const double                       aprecision ,
  const double                       rprecision )
{
  static const Ostap::Math::GSL::Integrator1D<function1> integrator {} ;
  auto F = integrator.make_function( &f1 ) ;
  //
  int    ierror ;
  double value  ;
  double error  ;
  static const char s_message[] = "Ostap::Math::Integrator/integrate_to_infinity" ;
  std::tie ( ierror , value , error ) = 
    integrator.gaqil_integrate 
    ( &F                 , 
      xmax               ,   // lower integration edge  
      workspace ( ws )   ,   // workspace 
      0 < aprecision ? aprecision : s_APRECISION_GAQIL ,   // absolute precision 
      0 < rprecision ? rprecision : s_RPRECISION_GAQIL ,   // relative precision 
      -1                 ,   // limit 
      s_message          ,   // reason of failure 
      __FILE__           ,   // the file 
      __LINE__           ,   // the line 
      tag                ) ; // tag/label 
  //
  return result  ( value , error ) ;
}
// ============================================================================
/** get Cauchy principal value integral 
 *  \f[ g(c) =  \int_{x_{min}}^{x_{max}}\frac{f(x)}{x-c}dx \f]
 *  @param f1 the function 
 *  @param c  the parameter 
 *  @param xmin lower integration edge 
 *  @param xmax upper  integration edge 
 *  @param integration workspace 
 *  @return value of the integral and the estimate of the uncertainty
 */
// ============================================================================
Ostap::Math::Integrator::result 
Ostap::Math::Integrator::cauchy_pv_
( Ostap::Math::Integrator::function1 f1         , 
  const double                       c          , 
  const double                       xmin       , 
  const double                       xmax       , 
  const Ostap::Math::WorkSpace&      ws         ,
  const std::size_t                  tag        ,
  const unsigned short               rescale    , 
  const double                       aprecision ,
  const double                       rprecision )
{
  if      ( s_equal ( xmax , xmin ) ) { return result ( 0 , 0 ) ; }
  //
  if (           xmax < xmin   ) 
  { 
    result r = cauchy_pv_ ( std::cref ( f1 ) , c , xmax , xmin , ws , tag , rescale , aprecision , rprecision ) ; 
    return result ( -r.first , r.second ) ;
  }
  //
  const double       scale = fun_scale ( f1 , xmin , xmax , rescale ) ;
  const double      iscale = s_zero    ( scale ) ? 1.0 : 1.0 / scale ;
  const std::size_t ntag   = 0 == tag ? tag : std::hash_combine ( tag , rescale , scale , iscale ) ;
  //
  auto  f2 = std::cref ( f1 ) ;
  //
  if ( s_equal ( c , xmin ) ) 
  {
    auto ff = [f2,c,xmin,iscale]  ( const double  x ) -> double 
      { return x <= xmin ? 0.0 : f2 ( x ) * iscale / ( x - c ) ; } ;
    //
    const double xlow = xmin - 0.1 * ( xmax - xmin ) ;
    result r = integrate_singular_ ( std::cref ( ff ) , xlow , xmax, { c } , ws , ntag , aprecision , rprecision ) ;
    return result ( r.first / iscale , r.second / iscale ) ;
  }
  //
  if ( s_equal ( c , xmax ) ) 
  {
    auto ff = [f2,c,xmax,iscale]  ( const double  x ) -> double 
      { return xmax <= x ? 0.0 : f2 ( x ) * iscale / ( x - c ) ; } ;
    //
    const double xhigh = xmax + 0.1 * ( xmax - xmin ) ;
    const std::size_t ntag = 
      0 == tag ? tag : std::hash_combine ( tag , rescale , scale , iscale ) ;
    result r = integrate_singular_ ( std::cref ( ff ) , xmin , xhigh , { c } , ws , ntag , aprecision , rprecision ) ;
    return result ( r.first / iscale , r.second / iscale ) ;
  }
  //
  auto  ff = [f2,c,iscale]  ( const double  x ) -> double 
    { return f2 ( x ) * iscale / ( x - c ) ; } ;
  //
  // regular integration  ?
  if ( c < xmin || c > xmax ) 
  { 
    result r = integrate_singular_ ( std::cref ( ff ) , xmin , xmax, { c } , ws , ntag , aprecision , rprecision ) ;
    return result ( r.first / iscale , r.second / iscale ) ;    
  }
  //
  // regular Cauchy integral
  //
  const double dx = std::min ( std::abs ( c - xmin ) , std::abs ( c - xmax ) ) / 2 ;
  //
  const double xlow  = c - dx ;
  const double xhigh = c + dx ;
  //
  const double  scale2 = fun_scale ( f1 , xlow , xhigh , rescale ) ;
  const double iscale2 = s_zero    ( scale2 ) ? 1.0 : 1.0 / scale2 ;
  function1 fs = [f2,iscale2] ( const double  x ) -> double { return f2 ( x ) * iscale2 ; } ;
  //
  const std::size_t n2tag   = 
    0 == tag ? tag : std::hash_combine ( tag , rescale , scale2 , iscale2 ) ;
  //
  static const Ostap::Math::GSL::Integrator1D<function1> integrator {} ;
  //
  auto F = s_equal ( iscale , 1.0 ) ? 
    integrator.make_function ( &f1 ) : 
    integrator.make_function ( &fs ) ;
  //
  int    ierror ;
  double value  ;
  double error  ;
  static const char s_message[] = "Ostap::Math::Integrator/integrate_cauchy_pv" ;
  std::tie ( ierror, value , error ) = 
    integrator.gawc_integrate 
    ( &F                , 
      xlow , xhigh      ,   // low and high integration edges 
      c                 ,   // Cauchy's point 
      workspace ( ws )  ,   // workspace 
      0 < aprecision ? aprecision : s_APRECISION_QAWC ,   // absolute precision 
      0 < rprecision ? rprecision : s_RPRECISION_QAWC ,   // relative precision 
      -1                ,   // limit 
      s_message         ,   // reason of failure 
      __FILE__          ,   // the file 
      __LINE__          ,   // the line 
      n2tag             ) ; // tag/label 
  //
  result r1 = integrate_singular_ ( std::cref ( ff ) , xmin   , c - dx , { c } , ws , ntag , aprecision , rprecision ) ;
  result r2 = integrate_singular_ ( std::cref ( ff ) , c + dx , xmax   , { c } , ws , ntag , aprecision , rprecision ) ;
  //
  return result ( value / iscale + r1.first  / iscale + r2.first  / iscale , 
                  error / iscale + r1.second / iscale + r2.second / iscale ) ;
  //
}
// ============================================================================
/*  get Cauchy principal value integral 
 *  \f[ g(c) =  \int_{x_{min}}^{+\infty}\frac{f(x)}{x-c}dx \f]
 *  @param f1 the function 
 *  @param c  the parameter 
 *  @param xmin lower integration edge 
 *  @param integration workspace 
 *  @return value of the integral and the estimate of the uncertainty
 */
// ============================================================================
Ostap::Math::Integrator::result 
Ostap::Math::Integrator::cauchy_pv_to_infinity_
( Ostap::Math::Integrator::function1 f1         , 
  const double                       c          , 
  const double                       xmin       , 
  const Ostap::Math::WorkSpace&      ws         , 
  const std::size_t                  tag        , 
  const unsigned short               rescale    ,
  const double                       aprecision ,
  const double                       rprecision )
{
  double xmax = 0 ;
  if ( s_equal ( c , xmin ) ) 
  {
    xmax = s_zero ( xmin ) ? 0.1 : 
      xmin < -1 ? 0.9 * xmin :
      xmin <  1 ? xmin + 0.1 : 1.1 * xmin ;
  }
  else if ( c < xmin ) { xmax = 2.2 * xmin - 1.2 * c    ; }
  else if ( c > xmin ) { xmax = 2.2 * c    - 1.2 * xmin ; }
  //
  auto f2 = std::cref ( f1 ) ;
  auto ff = [f2,c]  ( const double  x ) -> double { return f2 ( x ) / ( x - c ) ; } ;
  //
  result r1 = cauchy_pv_             ( std::cref ( f1 ) , c , xmin , xmax , ws , tag , rescale , aprecision , rprecision ) ;
  result r2 = integrate_to_infinity_ ( std::cref ( ff ) ,            xmax , ws , tag           , aprecision , rprecision ) ;
  //
  return result (  r1.first + r2.first , r1.second + r2.second ) ;
}
// ============================================================================
/*  get Cauchy principal value integral 
 *  \f[ g(c) =  \int_{-\infty}^{+x_{max}}\frac{f(x)}{x-c}dx \f]
 *  @param f1 the function 
 *  @param c  the parameter 
 *  @param xmax upper integration edge 
 *  @param integration workspace 
 *  @return value of the integral and the estimate of the uncertainty
 */
// ============================================================================
Ostap::Math::Integrator::result  
Ostap::Math::Integrator::cauchy_pv_from_infinity_
( Ostap::Math::Integrator::function1 f1         , 
  const double                       c          , 
  const double                       xmax       , 
  const Ostap::Math::WorkSpace&      ws         ,
  const std::size_t                  tag        ,
  const unsigned short               rescale    ,
  const double                       aprecision ,
  const double                       rprecision )
{
  //
  double xmin = 0 ;
  if ( s_equal ( c , xmax ) ) 
  {
    xmin = s_zero ( xmax ) ? -0.1 : 
      xmax >  1 ? 0.9 * xmax :
      xmin > -1 ? xmax - 0.1 : 1.1 * xmax ;
  }
  else if ( c < xmax ) { xmin = 2.2 * c    - 1.2 * xmax ; }
  else if ( c > xmax ) { xmin = 2.2 * xmax - 1.2 * c    ; }
  //
  auto f2 = std::cref ( f1 ) ;
  auto ff = [f2,c]  ( const double  x ) -> double { return f2 ( x ) / ( x - c ) ; } ;
  //
  result r1 = integrate_from_infinity_ ( std::cref ( ff ) ,            xmin , ws , tag           , aprecision , rprecision ) ;
  result r2 = cauchy_pv_               ( std::cref ( f1 ) , c , xmin , xmax , ws , tag , rescale , aprecision , rprecision ) ;
  //
  return result (  r1.first + r2.first , r1.second + r2.second ) ;
}
// ============================================================================
/*  Kramers-Kronig dispersion relation with n-subtractions 
 *  \f[ g(s) = \frac{s^n}{\pi} 
 *     \mathcal{P} \int_{x_{min}}^{+\infty} \frac{g(x)}{x^n(x-s)}dx \f] 
 *  @param f1    the function 
 *  @param s     s-parameter 
 *  @param xmin  lower integration range 
 *  @param n     number of subtracion
 *  @return value of the integral and the estimate of the uncertainty
 */
// ============================================================================
Ostap::Math::Integrator::result
Ostap::Math::Integrator::kramers_kronig_
( Ostap::Math::Integrator::function1 f1         , 
  const double                       s          , 
  const double                       xmin       , 
  const unsigned short               n          , 
  const Ostap::Math::WorkSpace&      ws         ,
  const std::size_t                  tag        , 
  const unsigned short               rescale    , 
  const double                       aprecision ,
  const double                       rprecision )
{
  // 
  if ( 0 < n ) 
  {
    auto f2 = std::cref ( f1 ) ;
    auto ff = [f2,n] ( const double x ) -> double { return f2 ( x ) / std::pow ( x , n ) ; } ;
    //
    const double scale = std::pow ( s , n ) ;
    result r = kramers_kronig_ ( std::cref ( ff ) , s , xmin , 0 , ws , tag , rescale , aprecision , rprecision ) ;
    //
    return result ( r.first * scale , r.second * scale ) ;
  }
  // no subtractions
  result r = cauchy_pv_to_infinity_ ( std::cref ( f1 ) , s , xmin , ws , tag , rescale , aprecision , rprecision ) ;
  return result ( r.first / M_PI , r.second / M_PI ) ;
}
// ============================================================================
/*  integration with known singular points 
 *  \f[ r = \int_{x_{min}}^{x_{max}}f_1(x) dx \f]
 *  @param  f1 the   function 
 *  @param xmin lower integration edge 
 *  @param xmax upper integration edge 
 *  @param points known singular points 
 *  @param ws integration workspace 
 *  @param tag unique tag/label 
 *  
 *  - Only singular poins between \f$ x_{min} \f$  and \f$ x_{max} \f$ are considered 
 *  - \f$ x_{min} \f$  and \f$ x_{max \f$ are also  considered as singular points 
 */
// ============================================================================
Ostap::Math::Integrator::result
Ostap::Math::Integrator::integrate_singular_
( function1                     f1         , 
  const double                  xmin       , 
  const double                  xmax       ,
  const std::vector<double>&    points     ,
  const Ostap::Math::WorkSpace& ws         , 
  const std::size_t             tag        , 
  const double                  aprecision ,
  const double                  rprecision )
{ 
  static const Ostap::Math::GSL::Integrator1D<function1> integrator {} ;
  auto F = integrator.make_function( &f1 ) ;
  //
  int    ierror ;
  double value  ;
  double error  ;
  static const char s_message[] = "Ostap::Math::Integrator/integrate_singular" ;
  std::tie ( ierror, value , error ) = 
    integrator.gaqp_integrate 
    ( &F                    , 
      xmin                  ,   // lower integration edge
      xmax                  ,   // high integration edge
      points                ,   // known singulatories 
      workspace ( ws )      ,   // workspace 
      0 < aprecision ? aprecision : s_APRECISION_GAQP     ,   // absolute precision 
      0 < rprecision ? rprecision : s_RPRECISION_GAQP     ,   // relative precision 
      -1                    ,   // limit 
      s_message             ,   // reason of failure 
      __FILE__              ,   // the file 
      __LINE__              ,   // the line 
      tag                   ) ; // tag/label 
  //
  return result ( value , error ) ;
}
// ============================================================================
/*  calculate the integral 
 *  \f[ r = \int_{x_{min}}^{x_{max}}f_2(x,y) dx \f]
 *  @param f2 the function 
 *  @param y parameter y
 *  @param xmin lower integration edge in x 
 *  @param xmax upper integration edge in x 
 *  @param tag unique label/tag 
 *  @param rescale rescale function for better numerical precision 
 *  @return the value of the integral 
 */
// ============================================================================
Ostap::Math::Integrator::result
Ostap::Math::Integrator::integrateX_
( Ostap::Math::Integrator::function2 f2         , 
  const double                       y          , 
  const double                       xmin       ,
  const double                       xmax       ,
  const Ostap::Math::WorkSpace&      ws         , 
  const std::size_t                  tag        ,
  const unsigned short               rescale    ,
  const double                       aprecision , 
  const double                       rprecision ) 
{
  auto f2_ = std::cref ( f2 ) ;
  auto f1  = std::bind ( f2_ , std::placeholders::_1 , y ) ;
  return integrate_ ( std::cref ( f1 ) , xmin , xmax , ws , tag , rescale , aprecision , rprecision ) ;
} 
// ============================================================================
/** calculate the integral 
 *  \f[ r = \int_{y_{min}}^{y_{max}}f_2(x,y) dy \f]
 *  @param f2 the function 
 *  @param x parameter x
 *  @param ymin lower integration edge in y 
 *  @param ymax upper integration edge in y 
 *  @param tag unique label/tag 
 *  @param rescale rescale function for better numerical precision 
 *  @return the value of the integral 
 */
// ============================================================================
Ostap::Math::Integrator::result
Ostap::Math::Integrator::integrateY_
( Ostap::Math::Integrator::function2 f2        , 
  const double                       x          , 
  const double                       ymin       ,
  const double                       ymax       ,
  const Ostap::Math::WorkSpace&      ws         ,
  const std::size_t                  tag        ,
  const unsigned short               rescale    ,
  const double                       aprecision , 
  const double                       rprecision ) 
{
  auto f2_ = std::cref ( f2 ) ;
  auto f1  = std::bind ( f2_ , x , std::placeholders::_1 ) ;
  return integrate_ ( std::cref ( f1 ) , ymin , ymax , ws , tag , rescale , aprecision , rprecision ) ;
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
 *  @return the value of the integral and the estimate for the error 
 */
// ==========================================================================
Ostap::Math::Integrator::result 
Ostap::Math::Integrator::integrate_
( Ostap::Math::Integrator::function2 f2         , 
  const double                       xmin       , 
  const double                       xmax       ,
  const double                       ymin       , 
  const double                       ymax       , 
  const Ostap::Math::WorkSpace&     /* ws */    , 
  const std::size_t                  tag        ,
  const double                       aprecision , 
  const double                       rprecision ) 
{ return integrate_ ( std::cref ( f2 ) , xmin , xmax , ymin , ymax , tag , aprecision , rprecision ) ; }
// ==========================================================================
/** calculate the integral 
 *  \f[ r = \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}}f_2(x,y) dx dy \f]
 *  @param f2 the function 
 *  @param xmin lower integration edge in x 
 *  @param xmax upper integration edge in x 
 *  @param ymin lower integration edge in y 
 *  @param ymax upper integration edge in y 
 *  @return the value of the integral and the estimate for the error 
 */
// ==========================================================================
Ostap::Math::Integrator::result
Ostap::Math::Integrator::integrate_
( Ostap::Math::Integrator::function2 f2         , 
  const double                       xmin       ,  
  const double                       xmax       ,
  const double                       ymin       , 
  const double                       ymax       ,
  const std::size_t                  tag        , 
  const double                       aprecision , 
  const double                       rprecision ) 
{ 
  //
  static const Ostap::Math::GSL::Integrator2D<function2> s_cubature{} ;
  static const char s_message[] = "Ostap::Math::Integrator/integrate(2D)" ;
  const auto F = s_cubature.make_function ( &f2 , xmin , xmax , ymin , ymax ) ;
  int     ierror =  0 ;
  double  value  =  1 ;
  double  error  = -1 ;
  std::tie ( ierror , value , error ) = s_cubature.cubature 
    ( &F          ,   // the function  
      100000      ,   // limits  
      0 < aprecision ? aprecision : s_APRECISION ,   // absolute precision 
      0 < rprecision ? rprecision : s_RPRECISION ,   // relative precision 
      //
      s_message   ,   // message 
      __FILE__    ,   // the file name 
      __LINE__    ,   // the line number 
      tag         ) ; // tag/label 
  //
  return result  ( value , error ) ;
}
// =============================================================================




// =============================================================================
//                                                                       The END 
// =============================================================================
