// =============================================================================
// Include files 
// =============================================================================
// STD&STL
// =============================================================================
#include <functional>
// =============================================================================
// Ostap
// =============================================================================
#include "Ostap/Hash.h"
#include "Ostap/Integrator.h"
#include "Ostap/Workspace.h"
// =============================================================================
// Local
// =============================================================================
#include "Integrator1D.h"
#include "Integrator2D.h"
#include "Integrator3D.h"
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
  inline double fun_scale 
  ( const FUNCTION&      fun     , 
    const double         xmin    , 
    const double         xmax    ,
    const unsigned short rescale )
  { 
    if ( 0 == rescale ) { return 1.0 ; }
    const double dx = ( xmax - xmin ) / ( rescale + 1 ) ;
    double scale = 0 ;
    // for ( unsigned short i = 1 ; i <= rescale ; ++i ) { scale += fun ( xmin + i * dx ) ; }
    for ( unsigned short i = 1 ; i <= rescale ; ++i ) 
    { scale += std::abs ( fun ( xmin + i * dx ) ) ; }
    return scale * ( xmin - xmax ) / rescale ;
  } 
  // ===========================================================================
}
// =============================================================================
// constructor with integration workspace the name 
// =============================================================================
Ostap::Math::Integrator::Integrator 
( const Ostap::Math::WorkSpace& ws ) 
  : m_qag_rule              ( GSL_INTEG_GAUSS61    )
  , m_abs_precision_qag     ( s_APRECISION_QAG     ) 
  , m_rel_precision_qag     ( s_RPRECISION_QAG     ) 
  , m_abs_precision_qagi    ( s_APRECISION_QAGI    ) 
  , m_rel_precision_qagi    ( s_RPRECISION_QAGI    ) 
  , m_abs_precision_qagiu   ( s_APRECISION_QAGIU   ) 
  , m_rel_precision_qagiu   ( s_RPRECISION_QAGIU   ) 
  , m_abs_precision_qagil   ( s_APRECISION_QAGIL   ) 
  , m_rel_precision_qagil   ( s_RPRECISION_QAGIL   ) 
  , m_abs_precision_qagp    ( s_APRECISION_QAGP    ) 
  , m_rel_precision_qagp    ( s_RPRECISION_QAGP    )
  , m_abs_precision_qawc    ( s_APRECISION_QAWC    ) 
  , m_rel_precision_qawc    ( s_RPRECISION_QAWC    )
  , m_abs_precision_cpv     ( s_APRECISION_QAWC    ) 
  , m_rel_precision_cpv     ( s_RPRECISION_QAWC    ) 
  , m_abs_precision_cpvi    ( s_APRECISION_QAWC    ) 
  , m_rel_precision_cpvi    ( s_RPRECISION_QAWC    ) 
  , m_abs_precision_kk      ( s_APRECISION_QAWC    ) 
  , m_rel_precision_kk      ( s_RPRECISION_QAWC    )
  , m_abs_precision_cquad   ( s_APRECISION_CQUAD   ) 
  , m_rel_precision_cquad   ( s_RPRECISION_CQUAD   ) 
  , m_abs_precision_romberg ( s_APRECISION_ROMBERG ) 
  , m_rel_precision_romberg ( s_RPRECISION_ROMBERG ) 
  , m_abs_precision_cube2   ( s_APRECISION_CUBE2D  ) 
  , m_rel_precision_cube2   ( s_RPRECISION_CUBE2D  ) 
  , m_abs_precision_cube3   ( s_APRECISION_CUBE3D  ) 
  , m_rel_precision_cube3   ( s_RPRECISION_CUBE3D  ) 
  , m_workspace ( ws )
{}
// =============================================================================
// constructor with integration workspace size & name 
// =============================================================================
Ostap::Math::Integrator::Integrator 
( const std::size_t    size         , 
  const unsigned short size_cquad   ,
  const unsigned short size_romberg ) 
  : Integrator ( Ostap::Math::WorkSpace ( size ,size_cquad , size_romberg ) )
{}
// ============================================================================
// constructor with the fictive name 
// ============================================================================
Ostap::Math::Integrator::Integrator 
( const std::string& /* name */ ) 
  : Integrator () 
{}
// ============================================================================
// set absolute/relative precision for GAG
// =============================================================================
void Ostap::Math::Integrator::set_precision_qag   
( const double aprec ,
  const double rprec ) 
{
  m_abs_precision_qag = 0 < aprec ? aprec : s_APRECISION_QAG ;
  m_rel_precision_qag = 0 < rprec ? rprec : s_RPRECISION_QAG ;
}
// =============================================================================
// set absolute/relative precision for GAGI
// =============================================================================
void Ostap::Math::Integrator::set_precision_qagi   
( const double aprec ,
  const double rprec ) 
{
  m_abs_precision_qagi = 0 < aprec ? aprec : s_APRECISION_QAGI ;
  m_rel_precision_qagi = 0 < rprec ? rprec : s_RPRECISION_QAGI ;
}
// =============================================================================
// set absolute/relative precision for GAGIL
// =============================================================================
void Ostap::Math::Integrator::set_precision_qagil   
( const double aprec ,
  const double rprec ) 
{
  m_abs_precision_qagil = 0 < aprec ? aprec : s_APRECISION_QAGIL ;
  m_rel_precision_qagil = 0 < rprec ? rprec : s_RPRECISION_QAGIL ;
}
// =============================================================================
// set absolute/relative precision for GAGIU
// =============================================================================
void Ostap::Math::Integrator::set_precision_qagiu
( const double aprec ,
  const double rprec ) 
{
  m_abs_precision_qagiu = 0 < aprec ? aprec : s_APRECISION_QAGIU ;
  m_rel_precision_qagiu = 0 < rprec ? rprec : s_RPRECISION_QAGIU ;
}
// =============================================================================
// set absolute/relative precision for QAGP
// =============================================================================
void Ostap::Math::Integrator::set_precision_qagp
( const double aprec ,
  const double rprec ) 
{
  m_abs_precision_qagp  = 0 < aprec ? aprec : s_APRECISION_QAGP ;
  m_rel_precision_qagp  = 0 < rprec ? rprec : s_RPRECISION_QAGP ;
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
// set absolute/relative precision for CQUAD
// =============================================================================
void Ostap::Math::Integrator::set_precision_cquad
( const double aprec ,
  const double rprec ) 
{
  m_abs_precision_cquad  = 0 < aprec ? aprec : s_APRECISION_CQUAD ;
  m_rel_precision_cquad  = 0 < rprec ? rprec : s_RPRECISION_CQUAD ;
}
// =============================================================================
// set absolute/relative precision for cubature 
// =============================================================================
void Ostap::Math::Integrator::set_precision_romberg
( const double aprec ,
  const double rprec ) 
{
  m_abs_precision_romberg  = 0 < aprec ? aprec : s_APRECISION_ROMBERG ;
  m_rel_precision_romberg  = 0 < rprec ? rprec : s_RPRECISION_ROMBERG ;
}
// =============================================================================
// set absolute/relative precision for 2d-cubature 
// =============================================================================
void Ostap::Math::Integrator::set_precision_cube2
( const double aprec ,
  const double rprec ) 
{
  m_abs_precision_cube2  = 0 < aprec ? aprec : s_APRECISION_CUBE2D ;
  m_rel_precision_cube2  = 0 < rprec ? rprec : s_RPRECISION_CUBE2D ;
}
// =============================================================================
// set absolute/relative precision for 3d-cubature 
// =============================================================================
void Ostap::Math::Integrator::set_precision_cube3
( const double aprec ,
  const double rprec ) 
{
  m_abs_precision_cube3  = 0 < aprec ? aprec : s_APRECISION_CUBE3D ;
  m_rel_precision_cube3  = 0 < rprec ? rprec : s_RPRECISION_CUBE3D ;
}
// =============================================================================
// set the QAG nitegrtaion rule 
// =============================================================================
void Ostap::Math::Integrator::set_qag_rule ( const int rule ) 
{
  m_qag_rule = 
    GSL_INTEG_GAUSS15 <= rule && rule <= GSL_INTEG_GAUSS61 <= rule ? rule :
    GSL_INTEG_GAUSS61 ;
} 
// =============================================================================
/*  calculate the integral 
 *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
 *  @param f1 the function 
 *  @param xmin lower integration edge 
 *  @param xmax upper  integration edge
 *  @param ws   integration workspace 
 *  @param tag   unique tag/label 
 *  @param rescale rescale function for better numerical precision  
 *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAG is used) 
 *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAG is used) 
 *  @param rule       the actual Gauss-Kronrod integration rule 
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
  const double                       rprecision , 
  const int                          rule       )
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
        0 == tag ? tag : Ostap::Utils::hash_combiner ( tag , rescale , scale , iscale ) ;
      result r = integrate_ ( std::cref ( ff ) , xmin , xmax , ws , ntag , 0 , aprecision , rprecision , rule ) ;
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
      0 < aprecision ? aprecision : s_APRECISION_QAG  ,   // absolute precision 
      0 < rprecision ? rprecision : s_RPRECISION_QAG  ,   // relative precision 
      -1                ,   // limit 
      s_message         ,   // reason of failure 
      __FILE__          ,   // the file 
      __LINE__          ,   // the line
      GSL_INTEG_GAUSS15 <= rule && 
      rule <= GSL_INTEG_GAUSS61 ? rule : GSL_INTEG_GAUSS61 , // the integrtaion rule
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
      0 < aprecision ? aprecision : s_APRECISION_QAGI ,   // absolute precision 
      0 < rprecision ? rprecision : s_RPRECISION_QAGI ,   // relative precision 
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
      0 < aprecision ? aprecision : s_APRECISION_QAGIU ,   // absolute precision 
      0 < rprecision ? rprecision : s_RPRECISION_QAGIU ,   // relative precision 
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
      0 < aprecision ? aprecision : s_APRECISION_QAGIL ,   // absolute precision 
      0 < rprecision ? rprecision : s_RPRECISION_QAGIL ,   // relative precision 
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
  const std::size_t ntag   = 0 == tag ? tag : Ostap::Utils::hash_combiner ( tag , rescale , scale , iscale ) ;
  //
  auto  f2 = std::cref ( f1 ) ;
  //
  if ( s_equal ( c , xmin ) ) 
  {
    auto ff = [f2,c,xmin,iscale]  ( const double  x ) -> double 
      { return x <= xmin ? 0.0 : f2 ( x ) * iscale / ( x - c ) ; } ;
    //
    const double xlow = xmin - 0.1 * ( xmax - xmin ) - 2 * std::abs ( c - xmin ) ;
    result r = integrate_singular_ ( std::cref ( ff ) , xlow , xmax, { c } , ws , ntag , aprecision , rprecision ) ;
    return result ( r.first / iscale , r.second / iscale ) ;
  }
  //
  if ( s_equal ( c , xmax ) ) 
  {
    auto ff = [f2,c,xmax,iscale]  ( const double  x ) -> double 
      { return xmax <= x ? 0.0 : f2 ( x ) * iscale / ( x - c ) ; } ;
    //
    const double xhigh = xmax + 0.1 * ( xmax - xmin ) + 2 * std::abs ( c - xmax ) ;
    const std::size_t ntag = 
      0 == tag ? tag : Ostap::Utils::hash_combiner ( tag , rescale , scale , iscale ) ;
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
    0 == tag ? tag : Ostap::Utils::hash_combiner ( tag , rescale , scale2 , iscale2 ) ;
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
  const double                       rprecision ,
  const double                       width      )
{
  const double xmax  = std::max ( xmin , cauchy_pv_b ( c , width ) ) ;
  //
  auto f2 = std::cref ( f1 ) ;
  auto fR = [f2,c,xmin]  ( const double  x ) -> double { return x <= xmin ? 0.0 : f2 ( x ) / ( x - c ) ; } ;
  //
  const double aprec = 0 < aprecision ? aprecision : s_APRECISION_QAWC ;
  const double rprec = 0 < aprecision ? aprecision : s_APRECISION_QAWC ;
  //
  result r1 = cauchy_pv_             ( std::cref ( f1 ) , c , xmin , xmax , ws , tag , rescale , aprec / 2 , rprec ) ;
  result r2 = integrate_to_infinity_ ( std::cref ( fR ) ,            xmax , ws , tag           , aprec / 2 , rprec ) ;
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
  const double                       rprecision ,
  const double                       width      ) 
{
  //
  const double xmin  = std::min ( xmax , cauchy_pv_a ( c , width ) ) ;
  //
  auto f2 = std::cref ( f1 ) ;
  auto fL = [f2,c,xmax]  ( const double  x ) -> double { return xmax <= x ? 0.0 : f2 ( x ) / ( x - c ) ; } ;
  //
  const double aprec = 0 < aprecision ? aprecision : s_APRECISION_QAWC ;
  const double rprec = 0 < aprecision ? aprecision : s_APRECISION_QAWC ;
  //
  result r1 = integrate_from_infinity_ ( std::cref ( fL ) ,            xmin , ws , tag           , aprec / 2 , rprec ) ;
  result r2 = cauchy_pv_               ( std::cref ( f1 ) , c , xmin , xmax , ws , tag , rescale , aprec / 2 , rprec ) ;
  //
  return result (  r1.first + r2.first , r1.second + r2.second ) ;
}
// ======================================================================================
/*  Cauchy principal value for the infinite range 
 *  \f[ g(c) = \mathcal{P} \int_{-\infty}^{+\infty}\frac{f(x)}{x-c}dx \f]
 *
 *  Integral is calculated as sum of three components
 *   - \f$ \int_{-\infty}^{a}\frac{f(x)}{x-c}dx \f$
 *   - \f$ \mathcal{P} \int_{a}^{b}\frac{f(x)}{x-c}dx \f$
 *   - \f$ int_{b}^{+\infty}\frac{f(x)}{x-c}dx \f$
 *
 *  @param f1      the function 
 *  @param c       the parameter 
 *  @param ws      integration workspace 
 *  @param tag     unuque label/tag 
 *  @param rescale rescale function for better numerical precision  
 *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAWC is used) 
 *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAWC is used) 
 *  @param width   width parameter
 *  @return value of the integral and the estimate of the uncertainty
 *
 *  @see Ostap::Math::Integrator::cauchy_pv
 *  @see Ostap::Math::Integrator::integrate_to_infinity 
 *  @see Ostap::Math::Integrator::integrate_from_infinity 
 *  @see Ostap::Math::Integrator::cauchy_pv_a
 *  @see Ostap::Math::Integrator::cauchy_pv_b
 */
// ======================================================================================
Ostap::Math::Integrator::result
Ostap::Math::Integrator::cauchy_pv_infinity_
( Ostap::Math::Integrator::function1 f1         , 
  const double                       c          , 
  const Ostap::Math::WorkSpace&      ws         ,
  const std::size_t                  tag        ,
  const unsigned short               rescale    ,
  const double                       aprecision , 
  const double                       rprecision ,
  const double                       width      ) 
{  
  //
  const double xmin  = cauchy_pv_a ( c , width ) ;
  const double xmax  = cauchy_pv_b ( c , width ) ;
  //
  const double aprec = 0 < aprecision ? aprecision : s_APRECISION_QAWC ;
  const double rprec = 0 < aprecision ? aprecision : s_APRECISION_QAWC ;
  //
  auto f2 = std::cref ( f1 ) ;
  auto fL = [f2,c,xmin]  ( const double  x ) -> double { return xmin <= x ? 0.0 : f2 ( x ) / ( x - c ) ; } ;
  auto fR = [f2,c,xmax]  ( const double  x ) -> double { return x <= xmax ? 0.0 : f2 ( x ) / ( x - c ) ; } ;
  //
  result rL = integrate_from_infinity_ ( std::cref ( fL ) ,            xmin , ws , tag           , aprec / 3 , rprec ) ;
  result rR = integrate_to_infinity_   ( std::cref ( fR ) ,     xmax        , ws , tag           , aprec / 3 , rprec ) ;
  result rC = cauchy_pv_               ( std::cref ( f1 ) , c , xmin , xmax , ws , tag , rescale , aprec / 2 , rprec ) ;
  //
  return result ( rL.first  + rR.first  + rC.first , rL.second + rR.second + rC.second ) ;
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
  const double                       rprecision ,
  const double                       width      ) 
{
  // 
  if ( 0 < n ) 
  {
    auto f2 = std::cref ( f1 ) ;
    auto ff = [f2,n] ( const double x ) -> double { return f2 ( x ) / std::pow ( x , n ) ; } ;
    //
    const double scale = std::pow ( s , n ) ;
    result r = kramers_kronig_ ( std::cref ( ff ) , s , xmin , 0 , ws , tag , rescale , aprecision , rprecision , width ) ;
    //
    return result ( r.first * scale , r.second * scale ) ;
  }
  // no subtractions
  result r = cauchy_pv_to_infinity_ ( std::cref ( f1 ) , s , xmin , ws , tag , rescale , aprecision , rprecision , width ) ;
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
      0 < aprecision ? aprecision : s_APRECISION_QAGP     ,   // absolute precision 
      0 < rprecision ? rprecision : s_RPRECISION_QAGP     ,   // relative precision 
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
Ostap::Math::Integrator::integrate2X_
( Ostap::Math::Integrator::function2 f2         , 
  const double                       y          , 
  const double                       xmin       ,
  const double                       xmax       ,
  const Ostap::Math::WorkSpace&      ws         , 
  const std::size_t                  tag        ,
  const unsigned short               rescale    ,
  const double                       aprecision , 
  const double                       rprecision ,
  const int                          rule       ) 
{
  auto f2_ = std::cref ( f2 ) ;
  auto f1  = std::bind ( f2_ , std::placeholders::_1 , y ) ;
  // new tag 
  static const std::string s_X { "integrateX"} ;
  const  std::size_t ntag = ( 0 == tag  ) ? tag : Ostap::Utils::hash_combiner ( tag , s_X , y ) ;
  return integrate_ ( std::cref ( f1 ) , xmin , xmax , ws , ntag , rescale , aprecision , rprecision , rule ) ;
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
Ostap::Math::Integrator::integrate2Y_
( Ostap::Math::Integrator::function2 f2        , 
  const double                       x          , 
  const double                       ymin       ,
  const double                       ymax       ,
  const Ostap::Math::WorkSpace&      ws         ,
  const std::size_t                  tag        ,
  const unsigned short               rescale    ,
  const double                       aprecision , 
  const double                       rprecision ,
  const int                          rule       ) 
{
  auto f2_ = std::cref ( f2 ) ;
  auto f1  = std::bind ( f2_ , x , std::placeholders::_1 ) ;
  // new tag 
  static const std::string s_Y { "integrateY"} ;
  const  std::size_t ntag = ( 0 == tag  ) ? tag : Ostap::Utils::hash_combiner ( tag , s_Y , x ) ;
  return integrate_ ( std::cref ( f1 ) , ymin , ymax , ws , ntag , rescale , aprecision , rprecision , rule ) ;
}
// =============================================================================
/* calculate the integral using double adaptive  CQUAD integrator 
 *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
 *  @param f1 the function 
 *  @param xmin lower integration edge 
 *  @param xmax upper integration edge
 *  @param ws   integration workspace 
 *  @param tag  unique tag/label 
 *  @param rescale rescale function for better numerical precision  
 *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAG is used) 
 *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAG is used) 
 *  @return value of the integral and the estimate of the uncertainty
 */
// =============================================================================
Ostap::Math::Integrator::result
Ostap::Math::Integrator::integrate_cquad_
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
      const std::size_t ntag = ( 0 == tag ) ? tag :
        Ostap::Utils::hash_combiner ( tag , rescale , scale , iscale ) ;
      result r = integrate_cquad_ ( std::cref ( ff ) , xmin , xmax , ws , ntag , 0 , aprecision , rprecision ) ;
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
  static const char s_message[] = "Ostap::Math::Integrator/integrate_cquad(1D)" ;
  //
  std::tie ( ierror, value , error ) = 
    integrator.cquad_integrate    
    ( &F                       , 
      xmin                     ,   // lower integration edge  
      xmax                     ,   // upper integration edge
      workspace_cquad ( ws )   ,   // workspace 
      0 < aprecision ? aprecision : s_APRECISION_CQUAD , // absolute precision 
      0 < rprecision ? rprecision : s_RPRECISION_CQUAD , // relative precision 
      s_message                ,   // reason of failure 
      __FILE__                 ,   // the file 
      __LINE__                 ,   // the line
      tag                      ) ; // label/tag
  //
  return result ( value , error ) ;  
}
// =============================================================================
/* calculate the integral using Romberg integrator 
 *  \f[ r = \int_{x_{min}}^{x_{max}} f_1(x) dx \f]
 *  @param f1 the function 
 *  @param xmin lower integration edge 
 *  @param xmax upper integration edge
 *  @param ws   integration workspace 
 *  @param tag  unique tag/label 
 *  @param rescale rescale function for better numerical precision  
 *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAG is used) 
 *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAG is used) 
 *  @return value of the integral and the estimate of the uncertainty
 */
// =============================================================================
Ostap::Math::Integrator::result
Ostap::Math::Integrator::integrate_romberg_
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
      const std::size_t ntag = ( 0 == tag ) ? tag :
        Ostap::Utils::hash_combiner ( tag , rescale , scale , iscale ) ;
      result r = integrate_romberg_ ( std::cref ( ff ) , xmin , xmax , ws , ntag , 0 , aprecision , rprecision ) ;
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
  static const char s_message[] = "Ostap::Math::Integrator/integrate_romberg(1D)" ;
  //
  std::tie ( ierror, value , error ) = 
    integrator.romberg_integrate    
    ( &F                       , 
      xmin                     ,   // lower integration edge  
      xmax                     ,   // upper integration edge
      workspace_romberg ( ws ) ,   // workspace 
      0 < aprecision ? aprecision : s_APRECISION_ROMBERG , // absolute precision 
      0 < rprecision ? rprecision : s_RPRECISION_ROMBERG , // relative precision 
      s_message                ,   // reason of failure 
      __FILE__                 ,   // the file 
      __LINE__                 ,   // the line
      tag                      ) ; // label/tag
  //
  return result ( value , error ) ;  
}
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
Ostap::Math::Integrator::integrate2_
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
  const auto F   = s_cubature.make_function ( &f2 , xmin , xmax , ymin , ymax ) ;
  int     ierror =  0 ;
  double  value  =  1 ;
  double  error  = -1 ;
  std::tie ( ierror , value , error ) = s_cubature.cubature 
    ( &F          ,   // the function  
      100000      ,   // limits  
      0 < aprecision ? aprecision : s_APRECISION_CUBE2D , // absolute precision 
      0 < rprecision ? rprecision : s_RPRECISION_CUBE2D , // relative precision 
      //
      s_message   ,   // message 
      __FILE__    ,   // the file name 
      __LINE__    ,   // the line number 
      tag         ) ; // tag/label 
  //
  return result  ( value , error ) ;
}
// =============================================================================
/* calculate the integral 
 *  \f[ r = \int_{x_{min}}^{x_{max}}
 *          \int_{y_{min}}^{y_{max}}
 *          \int_{z_{min}}^{z_{max}}f_3(x,y,z) dx dy dz\f]
 *  @param f3 the function 
 *  @param xmin lower integration edge in x 
 *  @param xmax upper integration edge in x 
 *  @param ymin lower integration edge in y 
 *  @param ymax upper integration edge in y 
 *  @param zmin lower integration edge in z 
 *  @param zmax upper integration edge in z 
 *  @param tag unique label/tag 
 *  @param aprecision absolute precision  (if non-positive s_APRECISION_CUBE is used) 
 *  @param aprecision relative precision  (if non-positive s_RPRECISION_CUBE is used) 
 *  @return the value of the integral and the estimate of the error 
 */
// =============================================================================
Ostap::Math::Integrator::result
Ostap::Math::Integrator::integrate3_
( Ostap::Math::Integrator::function3 f3         , 
  const double      xmin                 ,   
  const double      xmax                 ,
  const double      ymin                 , 
  const double      ymax                 ,
  const double      zmin                 , 
  const double      zmax                 ,
  const std::size_t tag                  ,
  const double      aprecision           , 
  const double      rprecision           ) 
{ 
  //
  static const Ostap::Math::GSL::Integrator3D<function3> s_cubature{} ;
  static const char s_message[] = "Ostap::Math::Integrator/integrate(3D)" ;
  const auto F   = s_cubature.make_function ( &f3 , xmin , xmax , ymin , ymax , zmin , zmax ) ;
  int     ierror =  0 ;
  double  value  =  1 ;
  double  error  = -1 ;
  std::tie ( ierror , value , error ) = s_cubature.cubature 
    ( &F          ,   // the function  
      500000      ,   // limits  
      0 < aprecision ? aprecision : s_APRECISION_CUBE3D ,   // absolute precision 
      0 < rprecision ? rprecision : s_RPRECISION_CUBE3D ,   // relative precision 
      //
      s_message   ,   // message 
      __FILE__    ,   // the file name 
      __LINE__    ,   // the line number 
      tag         ) ; // tag/label 
  //
  return result  ( value , error ) ;
}
// =============================================================================
// partial integration for 3D functions 
// ============================================================================
/* calculate the integral 
 *  \f[ r(z) = \int_{x_{min}}^{x_{max}}
 *             \int_{y_{min}}^{y_{max}} f_3(x,y,z) dx \f]
 *  @param f2 the function 
 *  @param z parameter z
 *  @param xmin lower integration edge in x 
 *  @param xmax upper integration edge in x 
 *  @param ymin lower integration edge in y 
 *  @param ymax upper integration edge in y 
 *  @param tag unique label/tag 
 *  @param aprecision absolute precision  (if non-positive s_APRECISION_CUBE2D is used) 
 *  @param aprecision relative precision  (if non-positive s_RPRECISION_CUBE2D is used) 
 *  @return the value of the integral and the estimate of the error 
 */
// ============================================================================
Ostap::Math::Integrator::result
Ostap::Math::Integrator::integrate3XY_
( Ostap::Math::Integrator::function3 f3 , 
  const double         z                , 
  const double         xmin             ,
  const double         xmax             ,
  const double         ymin             ,
  const double         ymax             ,
  const std::size_t    tag              , 
  const double         aprecision       , 
  const double         rprecision       ) 
{ 
  auto f3_ = std::cref ( f3 ) ;
  auto f2  = std::bind ( f3_ , std::placeholders::_1 , std::placeholders::_2 , z ) ;
  // new tag 
  static const std::string s_XY { "integrateXY"} ;
  const  std::size_t ntag = ( 0 == tag  ) ? tag : Ostap::Utils::hash_combiner ( tag , s_XY , z ) ;
  return integrate2_ ( std::cref ( f2 ) , 
                       xmin , xmax , 
                       ymin , ymax , 
                       ntag , aprecision , rprecision ) ;
}
// ============================================================================
/** calculate the integral 
 *  \f[ r(y) = \int_{x_{min}}^{x_{max}}
 *             \int_{z_{min}}^{z_{max}} f_3(x,y,z) dx \f]
 *  @param f3 the function 
 *  @param y parameter z
 *  @param xmin lower integration edge in x 
 *  @param xmax upper integration edge in x 
 *  @param zmin lower integration edge in z 
 *  @param zmax upper integration edge in z 
 *  @param tag unique label/tag 
 *  @param aprecision absolute precision  (if non-positive s_APRECISION_CUBE2D is used) 
 *  @param aprecision relative precision  (if non-positive s_RPRECISION_CUBE2D is used) 
 *  @return the value of the integral and the estimate of the error 
 */
// ============================================================================
Ostap::Math::Integrator::result
Ostap::Math::Integrator::integrate3XZ_
( Ostap::Math::Integrator::function3 f3 , 
  const double         y                , 
  const double         xmin             ,
  const double         xmax             ,
  const double         zmin             ,
  const double         zmax             ,
  const std::size_t    tag              , 
  const double         aprecision       , 
  const double         rprecision       )
{ 
  auto f3_ = std::cref ( f3 ) ;
  auto f2  = std::bind ( f3_ , std::placeholders::_1 , y , std::placeholders::_2 ) ;
  // new tag 
  static const std::string s_XZ { "integrateXZ"} ;
  const  std::size_t ntag = ( 0 == tag  ) ? tag : Ostap::Utils::hash_combiner ( tag , s_XZ , y ) ;
  return integrate2_ ( std::cref ( f2 ) , 
                       xmin , xmax , 
                       zmin , zmax , 
                       ntag , aprecision , rprecision ) ;
 }
// ============================================================================
/** calculate the integral 
 *  \f[ r(z) = \int_{y_{min}}^{y_{max}}
 *             \int_{z_{min}}^{z_{max}} f_3(x,y,z) dydz \f]
 *  @param f3 the function 
 *  @param x parameter x
 *  @param ymin lower integration edge in y 
 *  @param ymax upper integration edge in y 
 *  @param zmin lower integration edge in z 
 *  @param zmax upper integration edge in z 
 *  @param tag unique label/tag 
 *  @param aprecision absolute precision  (if non-positive s_APRECISION_CUBE2D is used) 
 *  @param aprecision relative precision  (if non-positive s_RPRECISION_CUBE2D is used) 
 *  @return the value of the integral and the estimate of the error 
 */
// ============================================================================
Ostap::Math::Integrator::result
Ostap::Math::Integrator::integrate3YZ_
( Ostap::Math::Integrator::function3 f3 , 
  const double         x                , 
  const double         ymin             ,
  const double         ymax             ,
  const double         zmin             ,
  const double         zmax             ,
  const std::size_t    tag              , 
  const double         aprecision       , 
  const double         rprecision       ) 
{
  auto f3_ = std::cref ( f3 ) ;
  auto f2  = std::bind ( f3_ , x , std::placeholders::_1 , std::placeholders::_2 ) ;
  // new tag 
  static const std::string s_YZ { "integrateYZ"} ;
  const  std::size_t ntag = ( 0 == tag  ) ? tag : Ostap::Utils::hash_combiner ( tag , s_YZ , x ) ;
  return integrate2_ ( std::cref ( f2 ) , 
                       ymin , ymax , 
                       zmin , zmax , 
                       ntag , aprecision , rprecision ) ;
}
// ============================================================================
/*  calculate the integral 
 *  \f[ r(y,z) = \int_{x_{min}}^{x_{max}}f_3(x,y,z) dx  \f]
 *  @param f2 the function 
 *  @param y parameter y
 *  @param z parameter z
 *  @param xmin lower integration edge in x 
 *  @param xmax upper integration edge in x 
 *  @param integration workspace 
 *  @param tag unique label/tag 
 *  @param rescale rescale function for better numerical precision 
 *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAG is used) 
 *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAG is used) 
 *  @return the value of the integral and the estimate of the error 
 */
// ============================================================================
Ostap::Math::Integrator::result
Ostap::Math::Integrator::integrate3X_
( Ostap::Math::Integrator::function3 f3         , 
  const double                       y          , 
  const double                       z          , 
  const double                       xmin       ,
  const double                       xmax       ,
  const Ostap::Math::WorkSpace&      ws         , 
  const std::size_t                  tag        , 
  const unsigned short               rescale    ,
  const double                       aprecision , 
  const double                       rprecision , 
  const int                          rule       ) 
{
  auto f3_ = std::cref ( f3 ) ;
  auto f1  = std::bind ( f3_ , std::placeholders::_1 , y , z ) ;
  // new tag 
  static const std::string s_X { "integrateX"} ;
  const  std::size_t ntag = ( 0 == tag  ) ? tag : Ostap::Utils::hash_combiner ( tag , s_X , y , z ) ;
  return integrate_ ( std::cref ( f1 ) , xmin , xmax , ws , ntag , rescale , aprecision , rprecision , rule ) ;
}
// ============================================================================
/*  calculate the integral 
 *  \f[ r(x,z) = \int_{y_{min}}^{y_{max}}f_3(x,y,z) dx  \f]
 *  @param f3 the function 
 *  @param x parameter x
 *  @param z parameter z
 *  @param ymin lower integration edge in y 
 *  @param ymax upper integration edge in y 
 *  @param integration workspace 
 *  @param tag unique label/tag 
 *  @param rescale rescale function for better numerical precision 
 *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAG is used) 
 *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAG is used) 
 *  @return the value of the integral and the estimate of the error 
 */
// ============================================================================
Ostap::Math::Integrator::result
Ostap::Math::Integrator::integrate3Y_
( Ostap::Math::Integrator::function3 f3         , 
  const double                       x          , 
  const double                       z          , 
  const double                       ymin       ,
  const double                       ymax       ,
  const Ostap::Math::WorkSpace&      ws         , 
  const std::size_t                  tag        , 
  const unsigned short               rescale    ,
  const double                       aprecision , 
  const double                       rprecision , 
  const int                          rule       ) 
{
  auto f3_ = std::cref ( f3 ) ;
  auto f1  = std::bind ( f3_ , x , std::placeholders::_1 , z ) ;
  // new tag 
  static const std::string s_Y { "integrateY"} ;
  const  std::size_t ntag = ( 0 == tag  ) ? tag : Ostap::Utils::hash_combiner ( tag , s_Y , x , z ) ;
  return integrate_ ( std::cref ( f1 ) , ymin , ymax , ws , ntag , rescale , aprecision , rprecision , rule ) ;
}
// ============================================================================
/*  calculate the integral 
 *  \f[ r(x,y) = \int_{z_{min}}^{z_{max}}f_3(x,y,z) dz  \f]
 *  @param f3 the function 
 *  @param x parameter x
 *  @param y parameter y
 *  @param zmin lower integration edge in z 
 *  @param zmax upper integration edge in z 
 *  @param integration workspace 
 *  @param tag unique label/tag 
 *  @param rescale rescale function for better numerical precision 
 *  @param aprecision absolute precision  (if non-positive s_APRECISION_QAG is used) 
 *  @param aprecision relative precision  (if non-positive s_RPRECISION_QAG is used) 
 *  @return the value of the integral and the estimate of the error 
 */
// ============================================================================
Ostap::Math::Integrator::result
Ostap::Math::Integrator::integrate3Z_
( Ostap::Math::Integrator::function3 f3         , 
  const double                       x          , 
  const double                       y          , 
  const double                       zmin       ,
  const double                       zmax       ,
  const Ostap::Math::WorkSpace&      ws         , 
  const std::size_t                  tag        , 
  const unsigned short               rescale    ,
  const double                       aprecision , 
  const double                       rprecision , 
  const int                          rule       ) 
{
  auto f3_ = std::cref ( f3 ) ;
  auto f1  = std::bind ( f3_ , x , y , std::placeholders::_1 ) ;
  // new tag 
  static const std::string s_Z { "integrateZ"} ;
  const  std::size_t ntag = ( 0 == tag  ) ? tag : Ostap::Utils::hash_combiner ( tag , s_Z , x , y ) ;
  return integrate_ ( std::cref ( f1 ) , zmin , zmax , ws , ntag , rescale , aprecision , rprecision , rule ) ;
}
// ============================================================================




// ============================================================================
// get low limit/cutoff for Cauchy-PV integration cutoff  
// ============================================================================
double Ostap::Math::Integrator::cauchy_pv_a  
( const double c     ,
  const double width ) 
{ 
  const double ac = std::abs ( c ) ;
  return ( 0 < width ) ? ( c - width ) : ac < 1 ? -2.0 : c - std::max ( 1.0 , 0.001 * ac ) ;
}
// ============================================================================
// get high limit/cutoff  for Cauchy-PV integrtaion cutoff  
// ============================================================================
double Ostap::Math::Integrator::cauchy_pv_b 
( const double c     ,
  const double width )
{
  const double ac = std::abs ( c ) ;
  return ( 0 < width ) ? ( c + width ) : ac < 1 ? +2.0 : c + std::max ( 1.0 , 0.001 * ac ) ;
}
// ============================================================================


// =============================================================================
//                                                                       The END 
// =============================================================================
