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
#include "Ostap/Workspace.h"
// ============================================================================
// Local
// ============================================================================
#include "Integrator1D.h"
#include "Integrator2D.h"
#include "local_math.h"
#include "local_hash.h"
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
  // ==========================================================================
  /** @var s_DI21
   *  The actual 2D-integrator 
   */
  const Ostap::Math::GSL::Integrator2D<Ostap::Math::DalitzIntegrator::function2> s_DI21 ;
  // ==========================================================================
  /** @var s_DI22
   *  The actual 2D-integrator 
   */
  const Ostap::Math::GSL::Integrator2D<Ostap::Math::DalitzIntegrator::function2> s_DI22 ;
  // ==========================================================================
  /** @var s_DI1
   *  The actual 1D-integrator 
   */
  const Ostap::Math::GSL::Integrator1D<Ostap::Math::DalitzIntegrator::function1> s_DI1 ;
  // ==========================================================================  
  /** @var s_MESSAGE1  
   *  Error message from integration
   */
  const char s_MESSAGE1[]  = "Integrate1(Dalitz)" ;
  // ==========================================================================  
  /** @var s_MESSAGE2
   *  Error message from cubature 
   */
  const char s_MESSAGE2[]  = "Integrate2(Dalitz)" ;
  // ==========================================================================  
  /** @var s_MAXCALLS1 
   *  Maximal number of (uncached)function calls 
   */
  const unsigned int s_MAXCALLS1 = 150000 ;
  // ==========================================================================
  /** @var s_MAXCALLS2 
   *  Maximal number of (cached)function calls 
   */
  const unsigned int s_MAXCALLS2 = 250000 ;
  // ==========================================================================  
  /** @var local integration workspace  
   */
  const Ostap::Math::WorkSpace  s_WS {} ;
  // ==========================================================================  
}
// ============================================================================

// ============================================================================
// 1D-integration (no workspace)
// ============================================================================

// ============================================================================
/*  evaluate integral over \f$s\f$ for \f$ f(s,s_1,s_2) \f$
 *  \f[ F(s_1,s_2)  = \int_{s_{mim}}^{s_{max}} ds f(s,s_1,s_2) \f]
 *  @param f3  the fun§ãtion \f$  f(s,s_1,s_2)\f$
 *  @param s1    value of \f$ s_1\f$
 *  @param s2    value of \f$ s_2\f$
 *  @param smax  upper inntegration limit for  \f$s\f$
 *  @param d     helper Dalitz-object 
 */
// ============================================================================
double Ostap::Math::DalitzIntegrator::integrate_s
( Ostap::Math::DalitzIntegrator::function3 f3   ,
  const double                             s1   ,
  const double                             s2   ,
  const double                             smax , 
  const Ostap::Kinematics::Dalitz0&        d    ) 
{ return integrate_s ( std::cref ( f3 ) , s1 , s2 , smax , d , s_WS ) ; }
// ============================================================================
/*  evaluate integral over \f$s_1\f$ for \f$ f(a,s_1,s_2) \f$
 *  \f[ F(s,s_2)  = \int  ds_1 f(s,s_1,s_2) \f]
 *  @param f3  the funсtion \f$  f(s,s_1,s_2)\f$
 *  @param s     value of \f$ s\f$
 *  @param s2    value of \f$ s_2\f$
 *  @param d     helper Dalitz-object 
 */
// ============================================================================
double Ostap::Math::DalitzIntegrator::integrate_s1
( Ostap::Math::DalitzIntegrator::function3 f3 ,
  const double                             s  ,
  const double                             s2 ,
  const Ostap::Kinematics::Dalitz0&        d  ) 
{ return integrate_s1 ( std::cref ( f3 ) , s , s2 , d , s_WS ) ; }
// ============================================================================
/*  evaluate integral over \f$s_1\f$ for \f$ f(s_1,s_2) \f$
 *  \f[ F(s_2)  = \int  ds_1 f(s_1,s_2) \f]
 *  @param f2  the funсtion \f$  f(s_1,s_2)\f$
 *  @param s2    value of \f$ s_2\f$
 *  @param d     helper Dalitz-object 
 */
// ============================================================================
double Ostap::Math::DalitzIntegrator::integrate_s1
( Ostap::Math::DalitzIntegrator::function2 f2 ,
  const double                             s2 ,
  const Ostap::Kinematics::Dalitz&         d  )
{ return integrate_s1  ( std::cref  ( f2 ) , s2 , d , s_WS ) ; }
// ============================================================================
/*  evaluate integral over \f$s_1\f$ for \f$ f(s_1,s_2) \f$
 *  \f[ F(s_2)  = \int  ds_1 f(s_1,s_2) \f]
 *  @param f2    the funсtion \f$  f(s_1,s_2)\f$
 *  @param s     value of \f$ s \f$
 *  @param s2    value of \f$ s_2\f$
 *  @param d     helper Dalitz-object 
 */
// ============================================================================
double Ostap::Math::DalitzIntegrator::integrate_s1
( Ostap::Math::DalitzIntegrator::function2 f2 ,
  const double                             s  ,
  const double                             s2 ,
  const Ostap::Kinematics::Dalitz0&        d  )
{ return integrate_s1 ( std::cref ( f2 ) , s , s2 , d , s_WS ) ; }
// ============================================================================


// ============================================================================
// 1D-integration with workspace 
// ============================================================================

// ============================================================================
/*  evaluate integral over \f$s\f$ for \f$ f(s,s_1,s_2) \f$
 *  \f[ F(s_1,s_2)  = \int_{s_{mim}}^{s_{max}} ds f(s,s_1,s_2) \f]
 *  @param f3  the funсtion \f$  f(s,s_1,s_2)\f$
 *  @param s1    value of \f$ s_1\f$
 *  @param s2    value of \f$ s_2\f$
 *  @param smax  upper inntegration limit for  \f$s\f$
 *  @param d     helper Dalitz-object 
 *  @param ws    integration workspace  
 */
// ============================================================================
double Ostap::Math::DalitzIntegrator::integrate_s
( Ostap::Math::DalitzIntegrator::function3 f3 ,
  const double                             s1   ,
  const double                             s2   ,
  const double                             smax , 
  const Ostap::Kinematics::Dalitz0&        d    ,
  const Ostap::Math::WorkSpace&            ws   )
{
  if  ( s1 <= d.s1_min (      ) ||
        s2 <= d.s1_min (      ) ||
        s1 >= d.s1_max ( smax ) ||
        s2 >= d.s2_max ( smax ) ) { return 0 ; }
  
  const double smin = s1 + s2 + d.s3_min() - d.summ2 () ;
  if ( smax  <= smin ) { return 0 ; }
  //
  auto ff = std::cref ( f3 ) ;
  function1 fun = [&ff,s1,s2,&d] ( const double s ) -> double
    { return d.inside ( s , s1 , s2 ) ? ff ( s , s1 , s2 ) : 0.0 ; } ;
  
  /// get integrator
  const auto F = s_DI1.make_function ( &fun ) ;
  //
  int    ierror  =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  std::tie ( ierror , result , error ) =
    s_DI1.gaq_integrate 
    ( &F                , 
      smin              ,   // lower integration edge  
      smax              ,   // upper integration edge
      workspace ( ws )  ,   // workspace 
      s_PRECISION       ,   // absolute precision 
      s_PRECISION       ,   // relative precision 
      -1                ,   // limit 
      s_MESSAGE1        ,   // reason of failure 
      __FILE__          ,   // the file 
      __LINE__          ) ; // the line 
  //
  return result ;  
}


// ============================================================================
/*  evaluate integral over \f$s_1\f$ for \f$ f(a,s_1,s_2) \f$
 *  \f[ F(s,s_2)  = \int  ds_1 f(s,s_1,s_2) \f]
 *  @param f3  the funсtion \f$  f(s,s_1,s_2)\f$
 *  @param s     value of \f$ s\f$
 *  @param s2    value of \f$ s_2\f$
 *  @param d     helper Dalitz-object 
 *  @param ws    integration workspace  
 */
// ============================================================================
double Ostap::Math::DalitzIntegrator::integrate_s1
( Ostap::Math::DalitzIntegrator::function3 f3 ,
  const double                             s  ,
  const double                             s2 ,
  const Ostap::Kinematics::Dalitz0&        d  ,
  const Ostap::Math::WorkSpace&            ws )
{
  if ( s < d.sqsumm () || s2 <= d.s2_min() || s2 >= d.s2_max ( s ) ) { return 0 ; }
  //
  auto ff = std::cref ( f3 ) ;
  function2 fun = [&ff,s,&d] ( const double s_1  , const double s_2 ) -> double
    { return d.inside ( s , s_1 , s_2 ) ? ff ( s , s_1 , s_2 ) : 0.0 ; } ;
  //
  return integrate_s1 ( std::cref ( fun ) , s , s2 , d , ws ) ;
}
// ============================================================================
/*  evaluate integral over \f$s_1\f$ for \f$ f(s_1,s_2) \f$
 *  \f[ F(s_2)  = \int  ds_1 f(s_1,s_2) \f]
 *  @param f2  the funсtion \f$  f(s_1,s_2)\f$
 *  @param s2    value of \f$ s_2\f$
 *  @param d     helper Dalitz-object 
 *  @param ws    integration workspace  
 */
// ============================================================================
double Ostap::Math::DalitzIntegrator::integrate_s1
( Ostap::Math::DalitzIntegrator::function2 f2 ,
  const double                             s2 ,
  const Ostap::Kinematics::Dalitz&         d  ,
  const Ostap::Math::WorkSpace&            ws )
{ return integrate_s1  ( std::cref  ( f2 ) , d.s() , s2  , d  , ws  ) ; }
// ============================================================================
/*  evaluate integral over \f$s_1\f$ for \f$ f(s_1,s_2) \f$
 *  \f[ F(s_2)  = \int  ds_1 f(s_1,s_2) \f]
 *  @param f2    the funсtion \f$  f(s_1,s_2)\f$
 *  @param s     value of \f$ s \f$
 *  @param s2    value of \f$ s_2\f$
 *  @param d     helper Dalitz-object 
 *  @param ws    integration workspace  
 */
// ============================================================================
double Ostap::Math::DalitzIntegrator::integrate_s1
( Ostap::Math::DalitzIntegrator::function2 f2 ,
  const double                             s  ,
  const double                             s2 ,
  const Ostap::Kinematics::Dalitz0&        d  ,
  const Ostap::Math::WorkSpace&            ws )
{
  if ( s < d.sqsumm () || s2 <= d.s2_min() || s2 >= d.s2_max ( s ) ) { return 0 ; }
  //
  auto ff = std::cref ( f2 ) ;
  function1 fun = [&ff,s,s2,&d] ( const double s1 ) -> double
    { return d.inside ( s , s1 , s2 ) ? ff ( s1 , s2 ) : 0.0 ; } ;
  
  double s1mn , s1mx ;
  std::tie ( s1mn , s1mx ) = d.s1_minmax_for_s_s2 ( s , s2 ) ;
  if ( s1mx  <= s1mn ) { return 0 ; }
  //
  
  /// get integrator
  const auto F = s_DI1.make_function ( &fun ) ;
  //
  int    ierror  =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  std::tie ( ierror , result , error ) =
    s_DI1.gaq_integrate 
    ( &F                , 
      s1mn              ,   // lower integration edge  
      s1mx              ,   // upper integration edge
      workspace ( ws )  ,   // workspace 
      s_PRECISION       ,   // absolute precision 
      s_PRECISION       ,   // relative precision 
      -1                ,   // limit 
      s_MESSAGE1        ,   // reason of failure 
      __FILE__          ,   // the file 
      __LINE__          ) ; // the line 
  //
  return result ;
}
// ============================================================================



// ============================================================================
// 1D-integration with cache 
// ============================================================================


// ============================================================================
/*  evaluate integral over \f$s\f$ for \f$ f(s,s_1,s_2) \f$
 *  \f[ F(s_1,s_2)  = \int_{s_{mim}}^{s_{max}} ds f(s,s_1,s_2) \f]
 *  @param tag tag that indicate the uniquness of function
 *  @param f3  the funсtion \f$  f(s,s_1,s_2)\f$
 *  @param s1    value of \f$ s_1\f$
 *  @param s2    value of \f$ s_2\f$
 *  @param smax  upper inntegration limit for  \f$s\f$
 *  @param d     helper Dalitz-object 
 *  @param ws    integration workspace  
 */
// ============================================================================
double Ostap::Math::DalitzIntegrator::integrate_s
( const std::size_t                        tag  ,
  Ostap::Math::DalitzIntegrator::function3 f3   ,
  const double                             s1   ,
  const double                             s2   ,
  const double                             smax , 
  const Ostap::Kinematics::Dalitz0&        d    ,
  const Ostap::Math::WorkSpace&            ws   )
{
  if  ( s1 <= d.s1_min (      ) ||
        s2 <= d.s1_min (      ) ||
        s1 >= d.s1_max ( smax ) ||
        s2 >= d.s2_max ( smax ) ) { return 0 ; }
  
  const double smin = s1 + s2 + d.s3_min() - d.summ2 () ;
  if ( smax  <= smin ) { return 0 ; }
  //
  auto ff = std::cref ( f3 ) ;
  function1 fun = [&ff,s1,s2,&d] ( const double s ) -> double
                  {
                    return d.inside ( s , s1 , s2 ) ? ff ( s , s1 , s2 ) : 0.0 ;
                  } ;
  
  /// get integrator
  const auto F = s_DI1.make_function ( &fun ) ;
  //
  int    ierror  =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  //
  std::tie ( ierror , result , error ) =
    s_DI1.gaq_integrate_with_cache  
    ( tag               ,
      &F                , 
      smin              ,   // lower integration edge  
      smax              ,   // upper integration edge
      workspace ( ws )  ,   // workspace 
      s_PRECISION       ,   // absolute precision 
      s_PRECISION       ,   // relative precision 
      -1                ,   // limit 
      s_MESSAGE1        ,   // reason of failure 
      __FILE__          ,   // the file 
      __LINE__          ) ; // the line 
  //
  return result ;  
}
// ============================================================================


// ============================================================================
/*  evaluate integral over \f$s_ы\f$ for \f$ f(a,s_1,s_2) \f$
 *  \f[ F(s,s_2)  = \int  ds_1 f(s,s_1,s_2) \f]
 *  @param tag tag that indicate the uniquness of function
 *  @param f3  the funсtion \f$  f(s,s_1,s_2)\f$
 *  @param s     value of \f$ s\f$
 *  @param s2    value of \f$ s_2\f$
 *  @param d     helper Dalitz-object 
 *  @param ws    integration workspace  
 */
// ============================================================================
double Ostap::Math::DalitzIntegrator::integrate_s1
( const std::size_t                        tag ,
  Ostap::Math::DalitzIntegrator::function3 f3  ,
  const double                             s   ,
  const double                             s2  ,
  const Ostap::Kinematics::Dalitz0&        d   ,
  const Ostap::Math::WorkSpace&            ws  )
{
  if ( s < d.sqsumm () || s2 <= d.s2_min() || s2 >= d.s2_max ( s ) ) { return 0 ; }
  //
  auto ff = std::cref ( f3 ) ;
  function2 fun = [&ff,s,&d] ( const double s_1  , const double s_2 ) -> double
                  {
                    return d.inside ( s , s_1 , s_2 ) ? ff ( s , s_1 , s_2 ) : 0.0 ;
                  } ;
  //
  return integrate_s1 ( tag , std::cref ( fun ) , s , s2 , d , ws ) ;
}
// ============================================================================
/*  evaluate integral over \f$s_ы\f$ for \f$ f(s_1,s_2) \f$
 *  \f[ F(s_2)  = \int  ds_1 f(s_1,s_2) \f]
 *  @param tag tag that indicate the uniquness of function
 *  @param f2  the funсtion \f$  f(s_1,s_2)\f$
 *  @param s2    value of \f$ s_2\f$
 *  @param d     helper Dalitz-object 
 *  @param ws    integration workspace  
 */
// ============================================================================
double Ostap::Math::DalitzIntegrator::integrate_s1
( const std::size_t                        tag ,
  Ostap::Math::DalitzIntegrator::function2 f2  ,
  const double                             s2  ,
  const Ostap::Kinematics::Dalitz&         d   ,
  const Ostap::Math::WorkSpace&            ws  )
{ return integrate_s1  ( tag , std::cref  ( f2 ) , d.s() , s2  , d  , ws  ) ; }
// ============================================================================
/*  evaluate integral over \f$s_ы\f$ for \f$ f(s_1,s_2) \f$
 *  \f[ F(s_2)  = \int  ds_1 f(s_1,s_2) \f]
 *  @param tag tag that indicate the uniquness of function
 *  @param f2    the funсtion \f$  f(s_1,s_2)\f$
 *  @param s     value of \f$ s \f$
 *  @param s2    value of \f$ s_2\f$
 *  @param d     helper Dalitz-object 
 *  @param ws    integration workspace  
 */
// ============================================================================
double Ostap::Math::DalitzIntegrator::integrate_s1
( const std::size_t                        tag ,
  Ostap::Math::DalitzIntegrator::function2 f2  ,
  const double                             s   ,
  const double                             s2  ,
  const Ostap::Kinematics::Dalitz0&        d   ,
  const Ostap::Math::WorkSpace&            ws  )
{
  if ( s < d.sqsumm () || s2 <= d.s2_min() || s2 >= d.s2_max ( s ) ) { return 0 ; }
  //
  auto ff = std::cref ( f2 ) ;
  function1 fun = [&ff,s,s2,&d] ( const double s1 ) -> double
                  {
                    return d.inside ( s , s1 , s2 ) ? ff ( s1 , s2 ) : 0.0 ;
                   } ;
  
  double s1mn , s1mx ;
  std::tie ( s1mn , s1mx ) = d.s1_minmax_for_s_s2 ( s , s2 ) ;
  if ( s1mx  <= s1mn ) { return 0 ; }
  //
  
  /// get integrator
  const auto F = s_DI1.make_function ( &fun ) ;
  //
  int    ierror  =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  std::tie ( ierror , result , error ) =
    s_DI1.gaq_integrate_with_cache  
    ( tag               ,
      &F                , 
      s1mn              ,   // lower integration edge  
      s1mx              ,   // upper integration edge
      workspace ( ws )  ,   // workspace 
      s_PRECISION       ,   // absolute precision 
      s_PRECISION       ,   // relative precision 
      -1                ,   // limit 
      s_MESSAGE1        ,   // reason of failure 
      __FILE__          ,   // the file 
      __LINE__          ) ; // the line 
  //
  return result ;
}
// ============================================================================




// ============================================================================
// 2D-integration
// ============================================================================

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
  const double x2_min = d.s2_min (   ) ;
  const double x2_max = d.s2_max ( M ) ;
  //
  double f_avg = 0 ;
  { 
    const unsigned short nx1 { 4 } ;
    const unsigned short nx2 { 5 } ;
    const double         dx1 = 2.0                 / nx1 ;
    const double         dx2 = ( x2_max - x2_min ) / nx2 ;
    //
    double f_value = 0.0 ;
    for ( unsigned short ix2 = 0 ; ix2 < nx2 ; ++ix2 ) 
    {
      const double x2 = x2_min + ( 0.5 + ix2 ) * dx2 ;
      for ( unsigned short ix1 = 0 ; ix1 < nx1 ; ++ix1 ) 
      {
        const double x1 = -1.0  + ( 0.5 + ix1 ) * dx1 ;
        double s1 , s2 ;
        std::tie ( s1 , s2 ) = d.x2s ( s , x1 , x2 ) ;
        f_avg += f2 ( s1 , s2 ) ;
      }
    }
    f_avg /= ( nx1 * nx2 ) ;
  }
  //
  const double f_norm = s_zero ( f_avg ) ? 1.0 : 1.0 / f_avg ;
  //
  function2 fun = [s,M,&d,&f2,f_norm] ( const double x1  , const double x2 ) -> double
    {
      double s1, s2 ;                     
      std::tie ( s1 , s2 ) = d.x2s ( s , x1 , x2 ) ;
      const double J = d.J ( s  , s1  , s2 ) ;
      return J <= 0 ? 0.0 : f2 ( s1 , s2 ) * J * f_norm  ;
    };  
  //
  /// get integrator
  const auto F = s_DI21.make_function ( &fun , -1  , 1 , x2_min , x2_max ) ;
  //
  int    ierror  =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  std::tie ( ierror , result , error ) = s_DI21.cubature
    ( &F , s_MAXCALLS1 , s_PRECISION , s_PRECISION , s_MESSAGE2 , __FILE__ , __LINE__ ) ;
  //
  return result / f_norm ;
}
// ============================================================================
/* evaluate the integral over \f$s\f$ , \f$s_1\f$ variables 
 *  \f[ \int\int f( s, s_1,s_2) ds ds_1 \f]  
 *  @param f3 the function \f$  f(s,s_1,s_2) \f$
 *  @param s2 fixed value of s2 
 *  @param smax upper-edge for integration over \f$ s \f$
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
  const double smin = d.sqsumm() + s2 - d.s2_min() ;
  if ( smax  <= smin ) { return 0 ; }
  //
  const double y1_min = smin ;
  const double y1_max = smax ;
  //
  double f_avg = 0 ;
  { 
    const unsigned short ny1 { 5 } ;
    const unsigned short ny2 { 4 } ;
    const double         dy1 = ( y1_max - y1_min ) / ny1 ;
    const double         dy2 = 2.0                 / ny2 ;
    //
    double f_value = 0.0 ;
    for ( unsigned short iy2 = 0 ; iy2 < ny2 ; ++iy2 ) 
    {
      const double y2 = -1.0  + ( 0.5 + iy2 ) * dy2 ;
      for ( unsigned short iy1 = 0 ; iy1 < ny1 ; ++iy1 ) 
      {
        const double y1 = y1_min  + ( 0.5 + iy1 ) * dy1 ;
        double s , s1 ;
        std::tie ( s , s1 ) = d.y2s ( s2 , y1 , y2 ) ;
        f_avg += f3 ( s , s1 , s2 ) ;
      }
    }
    f_avg /= ( ny1 * ny2 ) ;
  }
  //
  const double f_norm = s_zero ( f_avg ) ? 1.0 : 1.0 / f_avg ;
  //
  //
  function2  fun = [&d,s2,&f3,f_norm] ( const double y1  , const double y2 ) -> double
    {
      double s, s1 ;
      std::tie ( s , s1 ) = d.y2s ( s2 , y1 , y2 ) ;
      const double J = d.J ( s  , s1  , s2 ) ;
      return J <= 0 ? 0.0 : f3 ( s , s1 , s2 ) * J  * f_norm;
    };
  
  /// get the  integrator
  const auto F = s_DI22.make_function ( &fun , y1_min , y1_max  , -1  , 1 ) ;
  //
  int    ierror  =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  std::tie ( ierror , result , error ) = s_DI22.cubature
    ( &F , s_MAXCALLS1 , s_PRECISION , s_PRECISION , s_MESSAGE2 , __FILE__ , __LINE__ ) ;
  //
  return result / f_norm ;
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
  const double x2_min = d.s2_min (   ) ;
  const double x2_max = d.s2_max ( M ) ;
  //
  double f_avg = 0 ;
  { 
    const unsigned short nx1 { 4 } ;
    const unsigned short nx2 { 5 } ;
    const double         dx1 = 2.0                 / nx1 ;
    const double         dx2 = ( x2_max - x2_min ) / nx2 ;
    //
    double f_value = 0.0 ;
    for ( unsigned short ix2 = 0 ; ix2 < nx2 ; ++ix2 ) 
    {
      const double x2 = x2_min + ( 0.5 + ix2 ) * dx2 ;
      for ( unsigned short ix1 = 0 ; ix1 < nx1 ; ++ix1 ) 
      {
        const double x1 = -1.0  + ( 0.5 + ix1 ) * dx1 ;
        double s1 , s2 ;
        std::tie ( s1 , s2 ) = d.x2s ( s , x1 , x2 ) ;
        f_avg += f2 ( s1 , s2 ) ;
      }
    }
    f_avg /= ( nx1 * nx2 ) ;
  }
  //
  const double f_norm = s_zero ( f_avg ) ? 1.0 : 1.0 / f_avg ;
  //
  function2  fun = [s,M,&d,&f2,f_norm] ( const double x1  , const double x2 ) -> double
    {
      double s1, s2 ;                     
      std::tie ( s1 , s2 ) = d.x2s ( s  , x1  , x2 ) ;
      const double J = d.J ( s  , s1  , s2 ) ;
      return J <= 0 ? 0.0 : f2 ( s1 , s2 ) * J * f_norm ;
    };  
  //
  /// get integrator
  const auto F = s_DI21.make_function ( &fun , -1  , 1 , x2_min , x2_max ) ;
  //
  int    ierror  =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  const std::size_t key = std::hash_combine ( tag ,  s , d.tag() ) ;
  std::tie ( ierror , result , error ) = s_DI21.cubature_with_cache 
    ( key , &F , s_MAXCALLS2 , s_PRECISION , s_PRECISION , s_MESSAGE2 , __FILE__ , __LINE__ ) ;
  //
  return result / f_norm ;
}
// ============================================================================
/*  evaluate the integral over \f$s\f$ , \f$s_1\f$ variables 
 *  \f[ \int\int f( s, s_1,s_2) ds ds_1 \f]  
 *  @param tag tag that indicate theuniquness of function
 *  @param f3 the function \f$  f(s,s_1,s_2) \f$
 *  @param s2 fixed value of s2 
 *  @param smax upper-edge for inntegration over \f$ s \f$
 *  @param d  helper Dalitz-object 
 *  @return integral over \f$ s, s_1\f$
 */
// ============================================================================
double Ostap::Math::DalitzIntegrator::integrate_ss1
( const std::size_t                 tag ,
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
  const double y1_min = smin ;
  const double y1_max = smax ;
  //
  double f_avg = 0 ;
  { 
    const unsigned short ny1 { 5 } ;
    const unsigned short ny2 { 4 } ;
    const double         dy1 = ( y1_max - y1_min ) / ny1 ;
    const double         dy2 = 2.0                 / ny2 ;
    //
    double f_value = 0.0 ;
    for ( unsigned short iy2 = 0 ; iy2 < ny2 ; ++iy2 ) 
    {
      const double y2 = -1.0  + ( 0.5 + iy2 ) * dy2 ;
      for ( unsigned short iy1 = 0 ; iy1 < ny1 ; ++iy1 ) 
      {
        const double y1 = y1_min  + ( 0.5 + iy1 ) * dy1 ;
        double s , s1 ;
        std::tie ( s , s1 ) = d.y2s ( s2 , y1 , y2 ) ;
        f_avg += f3 ( s , s1 , s2 ) ;
      }
    }
    f_avg /= ( ny1 * ny2 ) ;
  }
  //
  const double f_norm = s_zero ( f_avg ) ? 1.0 : 1.0 / f_avg ;
  //
  function2  fun = [&d,s2,&f3,f_norm] ( const double y1  , const double y2 ) -> double
    {
      double s, s1 ;
      std::tie ( s , s1 ) = d.y2s ( s2 , y1 , y2 ) ;
      const double J = d.J ( s  , s1  , s2 ) ;
      return J <= 0 ? 0.0 : f3 ( s , s1 , s2 ) * J * f_norm  ;
    };
  
  /// get the  integrator
  const auto F = s_DI22.make_function ( &fun , smin , smax  , -1  , 1 ) ;
  //
  int    ierror  =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  //
  const std::size_t key = std::hash_combine ( tag , s2 , smax , d.tag() ) ;
  std::tie ( ierror , result , error ) = s_DI22.cubature_with_cache 
    ( key , &F , s_MAXCALLS2 , s_PRECISION , s_PRECISION , s_MESSAGE2 , __FILE__ , __LINE__ ) ;
  //
  return result / f_norm ;
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

