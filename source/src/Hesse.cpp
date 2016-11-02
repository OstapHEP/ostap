// $Id$
// ============================================================================
// Include files
// ============================================================================
// STD& STL
// ============================================================================
#include <cmath>
// ============================================================================
// GSL 
// ============================================================================
#include "gsl/gsl_math.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_linalg.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatusCode.h"
#include "Ostap/Hesse.h"
// ============================================================================
// Local
// ============================================================================
#include "GSL_sentry.h"
// ============================================================================
namespace 
{
  // ==========================================================================
  typedef Ostap::Math::GSL::GSL_Error_Handler Sentry ;
  // ==========================================================================
  typedef Ostap::Math::GSL::Hesse::function function ;
  // ==========================================================================
  /// 1/sqrt(2) 
  const double s_SQRT2i = 1.0 / std::sqrt ( 2.0 ) ;                // 1/sqrt(2) 
  // ==========================================================================
  // finite difference coefficients 
  // ==========================================================================
  const double s_8 [] = {
    -   1. / 560 ,
    +   8. / 315 , 
    -   1. /   5 , 
    +   8. /   5 , 
    - 205. /  72 , 
    +   8. /   5 , 
    -   1. /   5 , 
    +   8. / 315 , 
    -   1. / 560 
  } ;
  // ==========================================================================
  // finite difference coefficients 
  // ==========================================================================
  const double s_6 [] = {
    -   0.       , 
    +   1. /  90 ,
    -   3. /  20 , 
    +   3. /   2 , 
    -  49. /  18 , 
    +   3. /   2 , 
    -   3. /  20 ,  
    +   1. /  90 ,
    -   0.       
  } ;
  // ==========================================================================
  inline double dot     
  ( const double v1 [9] , 
    const double v2 [9] )
  {
    double result = 0 ;
    for ( unsigned int i = 0 ; i < 9 ; ++i ) 
    { result +=  v1[i] * v2[i] ; }
    return result ;
  }
  // ==========================================================================
  inline double absdot     
  ( const double  v1 [9] , 
    const double  v2 [9] )
  {
    double result = 0 ;
    for ( unsigned int i = 0 ; i < 8 ; ++i ) 
    { result += std::abs ( v1[i] * v2[i] ) ; }
    return result ;
  }
  // ==========================================================================
  inline double deriv2_8 
  ( const double   values [9]   , 
    const double   x            , 
    const double   h0           ,
    double*        abserr_round , 
    double*        abserr_trunc )
  {
    //
    const double h02 =   h0 * h0 ;
    //
    const double o8  =    dot ( values , s_8 ) ;
    const double o6  =    dot ( values , s_6 ) ;
    //
    const double o8a = absdot ( values , s_8 ) ;
    const double o6a = absdot ( values , s_6 ) ;
    //
    const double e6  = o6a * GSL_DBL_EPSILON ;
    const double e8  = o8a * GSL_DBL_EPSILON ;
    //
    const double r8  = o8 / h02 ;
    const double r6  = o6 / h02 ;
    //
    const double dy = 6 * 
      GSL_MAX ( std::abs ( r8 )  , std::abs ( r6 ) ) * 
      std::abs ( x / h0 ) * GSL_DBL_EPSILON ;
    //
    *abserr_trunc = std::abs  ( o6 - o8 ) / h02        ;
    *abserr_round = std::fabs ( e8 + e6 ) / h02  + dy  ;
    // 
    return r8 ;
  }
  // ===========================================================================
  /** calculate 
   *  a = x+h*d 
   *  @param x (INPUT)  vector x 
   *  @param d (INPUT)  vector d 
   *  @param h (INPUT)  scalar h 
   *  @param a (UPDATE) output vector a 
   */
  inline 
  void 
  _update_  
  ( const gsl_vector* x , 
    const gsl_vector* d ,
    const double      h , 
    gsl_vector*       a ) 
  {
    // fill the step vector 
    for ( std::size_t i = 0 ; i < x->size ; ++i ) 
    {
      const double xi = gsl_vector_get ( x , i ) ;
      const double di = gsl_vector_get ( d , i ) ; 
      gsl_vector_set   ( a , i , xi + h * di  ) ;  
    }
  }
  // ==========================================================================
  // get the second derivative along the direction d
  // ==========================================================================
  double deriv2_8
  ( const gsl_multimin_function* f  , 
    const gsl_vector*            x  ,
    const gsl_vector*            d  ,
    double                       h  ,
    double*                  round  , 
    double*                  trunc  , 
    gsl_vector*                  a ) 
  {
    //
    double values[9] ;
    const double h0 = 0.25 * std::abs ( h ) ;
    //
    
    for ( std::size_t j = 0 ; j < 9 ; ++j ) 
    {
      const double s = ( j - 4.0 )  * h0 ;      
      _update_ ( x , d , s , a ) ;
      values [ j ] = GSL_MULTIMIN_FN_EVAL_F( f , a  ) ;
    }
    //
    double xx = 0 ;
    for ( unsigned int j = 0 ; j < x->size ; ++j ) 
    {
      const double dj = gsl_vector_get ( d , j ) ; 
      if ( 0 == dj ) { continue ; }
      const double xj = gsl_vector_get ( x , j ) ;
      xx = std::max ( xx , std::abs ( xj ) + std::abs ( dj ) ) ;
    }
    //
    return deriv2_8 ( values , 
                      xx     , 
                      h0     , 
                      round  , 
                      trunc  ) ;
  }
  // ==========================================================================
  // get the second derivative along the direction d
  // ==========================================================================
  double deriv2
  ( const gsl_multimin_function* f , 
    const gsl_vector*            x ,
    const gsl_vector*            d ,
    double                       h ,
    double*                      e , 
    gsl_vector*                  a ) 
  {
    //
    double round ;
    double trunc ;
    //
    double result = deriv2_8 ( f  , x  , d  , h  , &round , &trunc , a  ) ;
    //
    *e = round + trunc ;
    //
    double         h_new = h ;
    double         res   = result ;
    unsigned short iter  = 0 ;
    while ( 0 < trunc && 
            0 < round &&
            ++iter <  16 )
    {
      //
      const double rt    = round / trunc / 2 ;
      const double corr  = 
        1 < rt  && rt < 1000 ? 
        pow ( rt , 1.0    ) :
        pow ( rt , 1. / 7 ) ;
      //
      if ( std::abs ( corr - 1 )  < 0.25 ) { break ; }   // BREAK   
      //
      if ( 0 != corr ) { h_new *= corr ; }
      else             { h_new /= 2    ; }
      //
      res   = deriv2_8 ( f , x , d , h_new , &round , &trunc , a ) ;
      //
      if ( iter > 2 && ( *e ) < ( trunc + round) ) { break    ; } // BREAK 
      //
      if (             ( *e ) < ( trunc + round) ) { continue ; } // CONTINUE 
      //
      result = res ;
      *e     = trunc + round ;
      //
    }
    //
    return result ;
  }
  // ==========================================================================
  // get the second derivative along pseudo-axis 
  // ==========================================================================
  double deriv2 
  ( const gsl_multimin_function* f , 
    const gsl_vector*            x ,
    const unsigned short         i , 
    const unsigned short         j ,
    double                       h ,
    double*                  error , 
    gsl_vector*                  a ,  // helper vector (workspace)
    gsl_vector*                  b )  // helper vector (workspace)
  {
    //
    gsl_vector_set_zero ( b ) ;    
    //
    // define the direction 
    if ( i == j  ) 
    { 
      gsl_vector_set ( b , i , 1 ) ; 
    }
    else 
    { 
      gsl_vector_set ( b , i ,  s_SQRT2i ) ;
      gsl_vector_set ( b , j ,  s_SQRT2i ) ;
    }
    //
    return deriv2  ( f , x , b , h , error , a ) ;
    // ========================================================================
  }
  // ==========================================================================
} // end of namespace 
// ============================================================================
// HESSE 
// ============================================================================
/*  constructor with all parameters 
 *  @param f      the function to be used 
 *  @param x      the point for hessian to be evaluated 
 *  @param params the parameters for the function 
 *  @param h the step-size (guess)
 */
// ============================================================================
Ostap::Math::GSL::Hesse::Hesse
( function          f      ,
  const gsl_vector* x      ,
  void*             params , 
  const double      h      ) 
//
  : m_func   ( f      ) 
  , m_x      ( x      ) 
  , m_params ( params ) 
//
  , m_h      ( std::abs ( h ) )
//
  , m_hesse  ( 0      )
  , m_aux    ( 0      )
  , m_cov2   ( 0      )
//
  , m_a      ( 0 ) 
  , m_b      ( 0 )
{
  m_a     = gsl_vector_calloc ( x -> size ) ;
  m_b     = gsl_vector_calloc ( x -> size ) ;
  //
}
// ============================================================================
/// destrictor 
// ============================================================================
Ostap::Math::GSL::Hesse::~Hesse ()
{
  // free matrices 
  if ( 0 != m_cov2  ) { gsl_matrix_free ( m_cov2  ) ; m_cov2  = 0 ; }
  if ( 0 != m_aux   ) { gsl_matrix_free ( m_aux   ) ; m_aux   = 0 ; }
  if ( 0 != m_hesse ) { gsl_matrix_free ( m_hesse ) ; m_hesse = 0 ; }
  // free vectors 
  if ( 0 != m_a     ) { gsl_vector_free ( m_a     ) ; m_a     = 0 ; }
  if ( 0 != m_b     ) { gsl_vector_free ( m_b     ) ; m_b     = 0 ; }
}
// ============================================================================
Ostap::StatusCode Ostap::Math::GSL::Hesse::calcHesse ()
{
  //
  // if ( 0 == m_x    ) { return InvalidPoint    ; }
  // if ( 0 == m_func ) { return InvalidFunction ; }
  //
  if ( 0 != m_hesse ) { gsl_matrix_free ( m_hesse ) ; m_hesse = 0 ; }
  if ( 0 != m_aux   ) { gsl_matrix_free ( m_aux   ) ; m_cov2  = 0 ; }
  //
  // allocate new matrix 
  //
  m_hesse = gsl_matrix_calloc ( size () , size () )  ;
  m_aux   = gsl_matrix_calloc ( size () , size () )  ;
  //
  double error = 0 ;
  //
  gsl_multimin_function F ;
  F.f      = m_func      ;
  F.params = m_params    ;
  F.n      = m_x -> size ;
  //
  // fill auxillary Hesse matrix 
  for ( std::size_t i = 0 ; i < m_aux->size1 ; ++i ) 
  {
    //
    for ( std::size_t j = 0 ; j < m_aux->size2 ; ++j )
    {
      //
      const double hij = deriv2  ( &F       , 
                                   m_x      , 
                                   i        , 
                                   j        , 
                                   m_h      , 
                                   &error   , 
                                   m_a      ,  
                                   m_b      ) ;
      //
      gsl_matrix_set ( m_aux , i , j , hij ) ;
      //
    }
    //
  }
  //
  // adjust hesse matrix 
  for ( std::size_t i = 0 ; i < m_hesse->size1 ; ++i ) 
  {
    const double hii = gsl_matrix_get ( m_aux , i , i ) ;
    //
    gsl_matrix_set ( m_hesse , i , i , hii ) ;
    //
    for ( std::size_t j = i + 1 ; j < m_hesse->size2 ; ++j )
    {
      //
      const double hjj = gsl_matrix_get ( m_aux , j , j ) ;
      //
      const double hij = gsl_matrix_get ( m_aux , i , j ) ;
      const double hji = gsl_matrix_get ( m_aux , j , i ) ;
      //
      const double h = 0.5 * ( hij + hji ) - 0.5 * ( hii + hjj ) ;
      //
      gsl_matrix_set ( m_hesse , i , j , h ) ;      
      gsl_matrix_set ( m_hesse , j , i , h ) ;
      //
    }
  }
  //   aux <--- hesse 
  gsl_matrix_memcpy ( m_aux , m_hesse ) ;
  //
  return StatusCode::SUCCESS ;
}
// ============================================================================
Ostap::StatusCode Ostap::Math::GSL::Hesse::calcCov2 ()
{
  //
  Ostap::Math::GSL::GSL_Error_Handler sentry ;
  //
  if ( 0 != m_cov2  ) { gsl_matrix_free ( m_cov2  ) ; m_cov2  = 0 ; }
  //
  if ( 0 == m_hesse || 0 == m_aux ) 
  {
    StatusCode sc = calcHesse () ;
    if  ( sc.isFailure() ) { return sc ; }
  }
  //
  m_cov2 = gsl_matrix_calloc ( size() , size() )  ;
  //
  // copy hesse into aux 
  gsl_matrix_memcpy ( m_aux , m_hesse ) ;;
  //
  // use LU decomposition 
  StatusCode sc = invert_LU_1 ( m_aux , m_cov2  ) ;
  if ( sc.isFailure() ) 
  {
    gsl_matrix_free ( m_cov2 ) ; m_cov2 =  0 ;
    return sc ;
  }
  //
  gsl_matrix_scale ( m_cov2 , 2 ) ;
  //
  return StatusCode::SUCCESS ;
}
// ============================================================================
/*  invert the matrix using LU decomposition
 *  @param matrix (UPDATE) the matrix to be inverted 
 *  @param result (UPDATE) the result 
 *  @return status code 
 *  @attention the input matrix will be screwed up!
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2012-05-28
 */
// ============================================================================
Ostap::StatusCode 
Ostap::Math::GSL::invert_LU_1 
( gsl_matrix* matrix , 
  gsl_matrix* result ) 
{
  //
  if ( 0 == matrix )                        { return 10 ; }
  if ( 0 == result )                        { return 11 ; }
  if ( matrix -> size1 != matrix -> size2 ) { return 12 ; }
  if ( result -> size1 != result -> size2 ) { return 13 ; }
  if ( result -> size1 != matrix -> size1 ) { return 14 ; }
  if ( result -> size2 != matrix -> size2 ) { return 15 ; }
  //
  // Make LU decomposition of input matrix 
  //
  // permutations :
	gsl_permutation * perm = gsl_permutation_alloc ( matrix -> size1  );
  //
	int sigdet ;
  int ierror ;
  //
	ierror = gsl_linalg_LU_decomp ( matrix , perm, &sigdet );
  if  ( ierror  ) 
  {
    gsl_permutation_free ( perm ) ;
    gsl_error ( "Error from LU_decomp " , __FILE__ , __LINE__ , ierror ) ; 
    return Ostap::StatusCode ( 100 + ierror ) ;
  }
  //
	// Invert the matrix m_aux 
	ierror = gsl_linalg_LU_invert ( matrix , perm, result  ) ;
  if  ( ierror  ) 
  {
    gsl_permutation_free ( perm ) ;
    gsl_error ( "Error from LU_invert " , __FILE__ , __LINE__ , ierror ) ; 
    return Ostap::StatusCode ( 200 + ierror ) ;
  }
  //
  gsl_permutation_free ( perm ) ;
  //
  return Ostap::StatusCode::SUCCESS ;
}
// ======================================================================
/*  invert the matrix using LU decomposition
 *  @param matrix (INPUT) the matrix to be inverted 
 *  @param result (UPDATE) the result 
 *  @return status code 
 *  @attention the input matrix will be preserved 
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2012-05-28
 */
// ======================================================================
Ostap::StatusCode 
Ostap::Math::GSL::invert_LU_2
( const gsl_matrix* matrix ,
  gsl_matrix*       result ) 
{
  //
  if ( 0 == matrix )                        { return 10 ; }
  if ( 0 == result )                        { return 11 ; }
  if ( matrix -> size1 != matrix -> size2 ) { return 12 ; }
  if ( result -> size1 != result -> size2 ) { return 13 ; }
  if ( result -> size1 != matrix -> size1 ) { return 14 ; }
  if ( result -> size2 != matrix -> size2 ) { return 15 ; }
  //
  // make the intermediaet matrix 
  //
  gsl_matrix* aux = gsl_matrix_alloc ( matrix->size1 , matrix->size2 ) ;
  //
  gsl_matrix_memcpy ( aux , matrix ) ;
  //
  Ostap::StatusCode sc = invert_LU_1 ( aux , result ) ;
  //
  // delete the auxillary matrix ;
  gsl_matrix_free ( aux ) ;
  //
  return sc ;
}  
// ============================================================================
// The END 
// ============================================================================
