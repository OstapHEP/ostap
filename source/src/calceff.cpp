

//$Id:%
// ============================================================================
// STD& STL
// ============================================================================
#include <cmath>
#include <limits>
#include <tuple>
// ============================================================================
// GSL 
// ============================================================================
#include "gsl/gsl_roots.h"
#include "gsl/gsl_multiroots.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Binomial.h"
#include "Ostap/MoreMath.h"
// ============================================================================
// Local 
// ============================================================================
#include "local_math.h"
// ============================================================================
/** @file
 *  implementation of functions from file Ostap/Binomial.h
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2015-09-17
 */
// ============================================================================
namespace
{
  // // ==========================================================================
  // typedef std::tuple<unsigned long,unsigned long,double> PARAMS ;
  // // ==========================================================================
  // double  B1Df ( double x , void* p )
  // {
  //   const PARAMS* params = static_cast<PARAMS*> ( p ) ;
  //   const unsigned long A = std::get<0> ( params ) ;
  //   const unsigned long B = std::get<1> ( params ) ;
  //   const double        P = std::get<2> ( params ) ;
  //   return Ostap::Math::beta_inc ( 1.0 * A , 1.0 * B , x ) - P ;
  // }
  // double B1Ddf ( double x , void* p )
  // {
  //   const PARAMS* params = static_cast<PARAMS*> ( p ) ;
  //   const unsigned long A = std::get<0> ( params ) ;
  //   const unsigned long B = std::get<1> ( params ) ;
  //   return Ostap::Math::dbeta_inc ( 1.0 * A , 1.0 * B , x ) ;
  // }
  // void   B1Dfdf
  // ( double  x  ,
  //   void*   p  , 
  //   double* f  ,
  //   double* df )
  // {
  //   const PARAMS* params = static_cast<PARAMS*> ( p ) ;
  //   const unsigned long A = std::get<0> ( params ) ;
  //   const unsigned long B = std::get<1> ( params ) ;
  //   const double        P = std::get<2> ( params ) ;
  //   *f  = Ostap::Math::beta_inc  ( 1.0 * A , 1.0 * B , x ) - P ;
  //   *df = Ostap::Math::dbeta_inc ( 1.0 * A , 1.0 * B , x )     ;
  // }
  // // ================================================================================
  
  // int bayes    ( gsl_vector* x , void* p , gsl_vector* f )
  // {
  //   const double alpha  = gsl_vector_get ( x , 0 ) ; // parameter alpha 
  //   const double beta   = gsl_vector_get ( x , 1 ) ; // parameter   beta 
  //   const double lambda = gsl_vector_get ( x , 2 ) ; // Lagrange  muliplier 
  //   //
  //   std::array<double,3> params = static_cast<std::array<double,3> > ( p ) ;
  //   //
  //   const double a = params [ 0 ] ;
  //   const double b = params [ 1 ] ;
  //   const double P = params [ 2 ] ;
  //   //
  //   gsl_vector_set ( f , 0 , 1 + lambda * Ostap::Math::dbeta_inc ( a , b , alpha ) ) ;
  //   gsl_vector_set ( f , 1 , 1 + lambda * Ostap::Math::dbeta_inc ( a , b , beta  ) ) ;
  //   gsl_vector_set ( f , 2 ,
  //                    Ostap::Math::beta_inc ( a , b , beta  ) -
  //                    Ostap::Math::beta_inc ( a , b , alpha ) - P ) ;
  //   //
  //   return GSL_SUCCESS ;
  // }
  // // ==========================================================================
  // int bayes_fd ( gsl_vector* x , void* p , gsl_matrix* J )
  // {
    
  //   const double alpha  = gsl_vector_get ( x , 0 ) ; // parameter alpha 
  //   const double beta   = gsl_vector_get ( x , 1 ) ; // parameter   beta 
  //   const double lambda = gsl_vector_get ( x , 2 ) ; // Lagrange  muliplier 
  //   //
  //   std::array<double,3> params = static_cast<std::array<double,3> > ( p ) ;
  //   //
  //   const double a = params [ 0 ] ;
  //   const double b = params [ 1 ] ;
  //   const double P = params [ 2 ] ;
  //   //

  // }
    
                 
  // ==========================================================================  
}
// ============================================================================
/*  Bayes' theorem-based interval
 *  @see M.Paterno, "Calculating efficiencies and their uncertainties", 
 *                   FERMILAB-TM-2286-CD
 *  @see https://inspirehep.net/literature/669498
 *  @see DOI: 10.2172/15017262
 */
// ============================================================================
std::pair<double,double>
Ostap::Math::bayes_interval
( const unsigned long accepted  ,
  const unsigned long rejected  ,
  const double        conflevel )
{
  if ( 0 == accepted && 0 == rejected ) { return std::make_pair ( 0.0 , 1.0 ) ; }
  //
  if ( 1 <= conflevel ) { return std::make_pair ( 0.0 , 1.0 ) ; }
  const double a = accepted ;
  const double p = a / ( a + rejected ) ;
  if ( 0 >= conflevel  ) { return std::make_pair ( p , p )    ; }
  //
  
  if ( !accepted ) { return std::make_pair ( 0.0 , 1 - std::pow ( 1 - conflevel , 1.0 / ( rejected + 1 ) )       )   ; }
  
  if ( !rejected ) { return std::make_pair (           std::pow ( 1 - conflevel , 1.0 / ( accepted + 1 ) ) , 1.0 ) ; }  
  
  // if ( !accepted )
  //   {        
  //     //
  //     PARAMS params { accepted , rejected , !accepted ? ( 1 - conflevel ) : conflevel }  ;
  //     // brent 
  //     gsl_function  F ;
  //     F  .f  = &B1D     ;
      
  //     // Brent 
  //     const gsl_root_fsolver_type* TB ;
  //     gsl_root_fsolver*            sb ;      
  //     TB = gsl_root_fsolver_brent (    ) ;
  //     sb = gsl_root_fsolver_alloc ( TB ) ;
  //     gsl_root_fsolver_set ( sb , &F , 0 , 1 ) ; 
  //     //

  //     // ======================================================================
  //     // several iterations of Brent
  //     // ======================================================================
  //     int biter   = 0 ;
  //     int bstatus = 0 ;
  //     double xlow = 0 ;
  //     double xigh = 1 ;
      
  //     for ( int biter = 0 ; biter < 5 ; ++biter )
  //       {
  //         int bstatus = gsl_root_fsolver_iterate ( sb ) ;
      
  //       }
      
  //     // switch to Newton 
      
  //     // Newton 
  //     const gsl_root_fdfsolver_type* T2 ;
  //     gsl_root_fdfsolver*            s2 ;      
  //     T2 = gsl_root_fdfsolver_newton (    ) ;
  //     s2 = gsl_root_fdfsolver_alloc  ( T2 ) ;
  //     //
  //     gsl_root_fdfsolver_set ( s , &F , 0 , 1 ) ; 


      
  //     // newtom 
  //     gsl_function_fdf FDF ;
  //     FDF.f   = &B1Df   ;
  //     FDF.df  = &B1Ddf  ;
  //     FDF.fdf = &B1Dfdf ;
  //     //
  //     F  .params = &params ;
  //     FDF.params = &params ;
  //     //
  //     //
      
  //   }
  
  return std::make_pair ( p , p  ) ;  
}

// ============================================================================
//                                                                      The END 
// ============================================================================


