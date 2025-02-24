// $Id$
// =============================================================================
// Include files 
// =============================================================================
// STD & STL
// =============================================================================
#include <cmath>
#include <limits>
#include <sstream>
// =============================================================================
// GSL
// =============================================================================
#include "gsl/gsl_errno.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_multimin.h"
// =============================================================================
// Ostap
// =============================================================================
#include "Ostap/StatusCode.h"
#include "Ostap/ValueWithError.h"
#include "Ostap/Hesse.h"
#include "Ostap/LinAlg.h"
#include "Ostap/Chi2Fit.h"
// =============================================================================
// Local
// =============================================================================
#include "GSL_sentry.h"
// =============================================================================
/** @file  
 *  Implementation file for class Gaudi::Math::Chi2Fit
 *
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2012-05-26
 */
// =============================================================================
namespace 
{
  // ==========================================================================
  const double s_Inf = 0.5 * std::numeric_limits<double>::max() ;
  // ==========================================================================
  class Chi2 
  {
    // ========================================================================
    friend  class Ostap::Math::Chi2Fit ;
    // ========================================================================
  public:
    // ========================================================================
    typedef Ostap::Math::Chi2Fit::VE    VE   ;
    typedef Ostap::Math::Chi2Fit::DATA  DATA ;
    typedef Ostap::Math::Chi2Fit::CMPS  CMPS ;
    // ========================================================================
  protected : // protected interface 
    // ========================================================================
    Chi2 ( const Ostap::Math::Chi2Fit& fit ) ;
    // ========================================================================
    ~Chi2 () ; // destructor
    // ========================================================================
  protected:
    // ========================================================================
    /// perform the fit 
    Ostap::StatusCode fit () ;
    // ========================================================================
  private:    // disabled methods 
    // ========================================================================
    /// the default constructor is disabled 
    Chi2 () ; // the default constructor is disabled 
    /// the copy constructor is disabled 
    Chi2 ( const Chi2& ) ; // the copy constructor is disabled 
    /// the assigmenet operator is disabled 
    Chi2& operator=( const Chi2& ) ; // the assigmenet operator is disabled 
    // ========================================================================
  public : // GSL-specific stuff
    // ======================================================================
    /// calculate the function value 
    double f      ( const gsl_vector* x ) const ; // function
    /// calcualate the gradient 
    void   df     ( const gsl_vector* x , 
                    gsl_vector*       g ) const ; // gradient 
    /// calculate the function and gradient 
    double fdf    ( const gsl_vector* x , 
                    gsl_vector*       g ) const ; // f & df 
    /// calculate the function, gradient and hessian 
    double fdfddf ( const gsl_vector* x , 
                    gsl_vector*       g , 
                    gsl_matrix*       h ) const ; // f & df & ddf 
    // ========================================================================
  public:
    // ========================================================================
    const gsl_vector* solution   () const { return m_solution   ; }
    // const gsl_matrix* hessian    () const { return m_hessian    ; }    
    const gsl_matrix* covariance () const { return m_covariance ; }    
    // ========================================================================
  public:
    // ========================================================================
    double      chi2         () const { return m_chi2         ; }    
    std::size_t calls_f      () const { return m_calls_f      ; }
    std::size_t calls_df     () const { return m_calls_df     ; }
    std::size_t calls_fdf    () const { return m_calls_fdf    ; }
    std::size_t calls_fdfddf () const { return m_calls_fdfddf ; }
    std::size_t points       () const { return m_points       ; }
    std::size_t iters        () const { return m_iter         ; }
    // ========================================================================
  private: 
    // ========================================================================
    const Ostap::Math::Chi2Fit* m_fit         ;
    //
    mutable std::size_t m_calls_f      ;
    mutable std::size_t m_calls_df     ;
    mutable std::size_t m_calls_fdf    ;
    mutable std::size_t m_calls_fdfddf ;
    mutable std::size_t m_points       ;
    // ========================================================================
    /// the solution 
    gsl_vector* m_solution   ;  // the current solution 
    /// the hessian 
    gsl_matrix* m_hessian    ;  // the hessian 
    /// the covariance 
    gsl_matrix* m_covariance ;   // the covariance 
    // ========================================================================
    double      m_chi2 ;
    std::size_t m_iter ;
    // ========================================================================
  } ;  
  // ==========================================================================
  // GSL specific stuff 
  // ==========================================================================
  /// get the function value 
  double chi2_f ( const gsl_vector *x, void * params )
  {
    // 
    const Chi2* fit = (const Chi2*) params ;
    // 
    return fit->f ( x )  ;
  }
  // ==========================================================================
  /// get the fradient 
  void chi2_df ( const gsl_vector *x, void *params, gsl_vector *df )
  {
    // 
    const Chi2* fit = (const Chi2*) params ;
    // 
    fit->df ( x , df )  ;
  }
  // ==========================================================================
  /// get the function and gradient 
  void chi2_fdf ( const gsl_vector *x, void *params, double *f, gsl_vector *df) 
  {
    //
    //
    const Chi2* fit = (const Chi2*) params ;
    // 
    *f = fit->fdf ( x , df )  ;
    //
  }
  // ===========================================================================
  // implementation of class Chi2 
  // ===========================================================================
  /// contribution to chi2  
  inline Chi2::VE term_chi2 
  ( Chi2::VE             res   ,    
    const Chi2::CMPS&    cmps  , 
    const unsigned short index , 
    const gsl_vector*    x     ) 
  {
    //
    for ( unsigned int j = 0 ; j < cmps.size() ; ++j ) 
    {
      // get the parameter from GSL 
      const double    aj = gsl_vector_get ( x , j ) ;
      // get the component value 
      const Chi2::VE& cj = cmps [ j ] [ index ] ;
      //
      res  -= aj * cj ;
    }
    //
    return res ;
  }
  // ========================================================================
  /// get the contribution to the gradient 
  inline 
  void term_grad 
  ( const Chi2::VE&      res   , 
    const Chi2::CMPS&    cmps  , 
    const unsigned short index , 
    const gsl_vector*    x     ,
    gsl_vector*          g     )
  {
    //
    const double rv     = res.value ()  ;
    const double rv2    = rv * rv       ;
    const double cov2   = res.cov2  ()  ;
    //
    // 3b. the second loop over components 
    for ( unsigned int j = 0 ; j < cmps.size() ; ++j ) 
    {
      // get the parameter from GSL 
      const double    aj = gsl_vector_get ( x , j ) ;
      // get the component value 
      const Chi2::VE& cj  = cmps[j][index] ;
      //
      // d(res)   / d(l_j)
      const double drdj =         - cj.value () ;
      // d(cov2)  / d(l_j)
      const double dsdj =  2 * aj * cj.cov2  () ;
      //
      // d(chi2)  / d(l_j) 
      //
      double dcdj   = 2 * rv  / cov2        * drdj ;
      dcdj         += -   rv2 / cov2 / cov2 * dsdj ;
      //
      // update the gradient
      const double gold = gsl_vector_get ( g , j ) ;
      gsl_vector_set ( g , j , gold + dcdj  ) ; // update gradient
    }
    //
  }
  // ==========================================================================
  /// get the contributions to the hessian 
  inline 
  void term_hesse
  ( const Chi2::VE&      res   , 
    const Chi2::CMPS&    cmps  , 
    const unsigned short index , 
    const gsl_vector*    x     ,
    gsl_matrix*          h     )
  {
    //
    const double rv     = res.value ()  ;
    const double rv2    = rv * rv       ;
    const double cov2   = res.cov2  ()  ;
    //
    // 3b. the second loop over components 
    for ( unsigned int j = 0 ; j < cmps.size() ; ++j ) 
    {
      // get the parameter from GSL 
      const double    aj = gsl_vector_get ( x , j ) ;
      // get the component value 
      const Chi2::VE& cj = cmps [ j ][ index ] ;
      //
      // d(res)   / d(l_j)
      const double drdj =         - cj.value () ;
      // d(cov2)  / d(l_j)
      const double dsdj =  2 * aj * cj.cov2  () ;
      //
      // next step: hessian 
      for ( unsigned  int k = j ; k < cmps.size() ; ++k )
      {
        // get the parameter from GSL 
        const double    ak = gsl_vector_get ( x , k ) ;
        // get the component value 
        const Chi2::VE& ck = cmps [ k ][ index ] ;
        //
        // d(res)   / d(l_k)
        const double drdk  =         - ck.value () ;
        // d(cov2)  / d(l_k)
        const double dsdk  =  2 * ak * ck.cov2  () ;
        //
        const double d2sjk = ( j == k ) ? 2 * cj.cov2() : 0.0 ;
        //
        //
        double dh = 0                                      ;
        dh += -     rv2 / cov2 / cov2        * d2sjk       ; 
        dh +=   2              / cov2        * drdj * drdk ; 
        dh += - 2 * rv  / cov2 / cov2        * drdj * dsdk ;
        dh += - 2 * rv  / cov2 / cov2        * drdk * dsdj ;
        dh += - 2 * rv2 / cov2 / cov2 / cov2 * dsdj * dsdk ;
        //
        const double hjk = gsl_matrix_get ( h , j , k ) ;
        //
        if ( j == k ) 
        { 
          gsl_matrix_set ( h , j , k , hjk + dh ) ; 
        }
        else 
        { 
          // enforce the symmetry:
          gsl_matrix_set ( h , j , k , hjk + dh ) ; 
          gsl_matrix_set ( h , k , j , hjk + dh ) ;
        }
        //
      }
      //
    }
    //
  }
  // ===========================================================================
  // calculate the function value 
  // ===========================================================================
  double Chi2::f ( const gsl_vector* x ) const // function
  {
    //
    ++m_calls_f   ;
    //
    const DATA& data = m_fit -> data () ;
    const CMPS& cmps = m_fit -> cmps () ;
    //
    m_points    = 0 ;
    double chi2 = 0 ;
    // the loop over the data-points 
    for ( std::size_t i = 0 ; i < data.size() ; ++ i ) 
    {
      // calculate the chi2-term 
      VE res = term_chi2 ( data[i] , cmps , i , x ) ;
      //
      if ( 0 >= res.cov2() ) { continue ; } // CONTINUE 
      //
      const double rv2 = res.value() * res.value() ;
      //
      chi2 += rv2 / res.cov2() ; // accumulate 
      //
      ++m_points ;
      //
    }
    //
    return chi2 ;
  } 
  // ==========================================================================
  // calculate the gradient
  // ==========================================================================
  void Chi2::df ( const gsl_vector* x , gsl_vector* g ) const 
  { 
    ++m_calls_df  ;
    fdf ( x , g )  ; 
    --m_calls_fdf ;
  }
  // ==========================================================================
  // calculate the function and gradient 
  // ==========================================================================
  double Chi2::fdf ( const gsl_vector* x , gsl_vector *g ) const
  {
    ++m_calls_fdf  ;
    //
    m_points = 0 ;
    // 1. reset the gradient 
    gsl_vector_set_zero ( g ) ;
    //    
    const DATA& data = m_fit -> data () ;
    const CMPS& cmps = m_fit -> cmps () ;
    //
    double chi2 = 0 ;
    // 2. the loop over the data points 
    for ( std::size_t i = 0 ; i < data.size() ; ++ i ) 
    {
      //
      VE res = term_chi2 ( data[i] , cmps , i , x ) ;
      //
      if ( 0 >= res.cov2() ) { continue ; } // CONTINUE 
      //
      const double rv2 = res.value() * res.value() ;
      //
      chi2 += rv2 / res.cov2() ; // accumulate 
      //
      ++m_points;
      //
      // update the gradient 
      term_grad  ( res , cmps , i , x , g ) ;
      //
    }
    //
    return chi2 ;
  }
  // ==========================================================================
  // calculate the function, gradient and hesse-matrix
  // ==========================================================================
  double Chi2::fdfddf ( const gsl_vector* x , 
                        gsl_vector*       g , 
                        gsl_matrix*       h ) const
  {
    ++m_calls_fdfddf ;
    //
    m_points = 0 ;
    // 1. reset the gradient & hessian  
    if ( 0 != g ) { gsl_vector_set_zero ( g ) ; }
    gsl_matrix_set_zero ( h ) ;
    //    
    const DATA& data = m_fit -> data () ;
    const CMPS& cmps = m_fit -> cmps () ;
    //
    double chi2 = 0 ;
    // 2. the loop over the data points 
    for ( std::size_t i = 0 ; i < data.size() ; ++ i ) 
    {
      //
      VE res = term_chi2 ( data[i] , cmps , i , x ) ;
      //
      if ( 0 >= res.cov2() ) { continue ; } // CONTINUE 
      //
      const double rv2 = res.value() * res.value() ;
      //
      chi2 += rv2 / res.cov2() ; // accumulate 
      //
      ++m_points;
      //
      // update the gradient
      if ( 0 != g ) { term_grad  ( res , cmps , i , x , g ) ; } 
      //
      // update the hessian
      term_hesse ( res , cmps , i , x , h ) ;
      //
    }
    //
    return chi2 ;
  }
  // ==========================================================================
  // constructor 
  // ==========================================================================
  Chi2::Chi2 ( const Ostap::Math::Chi2Fit& fit ) 
    : m_fit          ( &fit ) 
    , m_calls_f      ( 0 )
    , m_calls_df     ( 0 ) 
    , m_calls_fdf    ( 0 ) 
    , m_calls_fdfddf ( 0 ) 
    , m_points       ( 0 ) 
    , m_solution     ( 0 ) 
    , m_hessian      ( 0 )
    , m_covariance   ( 0 )
    , m_chi2         ( s_Inf )
    , m_iter         ( 0 ) 
  {} 
  // ==========================================================================
  /// destructor
  // ==========================================================================
  Chi2::~Chi2 ()  // destructor
  {
    /// free allocated vectors & matrices 
    if ( 0 != m_solution   ) 
    { gsl_vector_free ( m_solution   ) ; m_solution   = 0 ; }
    if ( 0 != m_hessian    ) 
    { gsl_matrix_free ( m_hessian    ) ; m_hessian    = 0 ; }
    if ( 0 != m_covariance ) 
    { gsl_matrix_free ( m_covariance ) ; m_covariance = 0 ; }
    // ========================================================================
  }
  // ==========================================================================
  // perform the fit 
  // ==========================================================================
  Ostap::StatusCode Chi2::fit () 
  {
    //
    Ostap::Math::GSL::GSL_Error_Handler sentry ;
    //
    const std::size_t size = m_fit->cmps().size() ;
    /// the function  
    gsl_multimin_function_fdf F ;
    F.n      = size        ;
    F.f      = &chi2_f     ;
    F.df     = &chi2_df    ;
    F.fdf    = &chi2_fdf   ;
    F.params = (void*)this ;
    // ========================================================================
    /// the minimizer 
    const gsl_multimin_fdfminimizer_type *T = 
      // gsl_multimin_fdfminimizer_conjugate_fr;
      gsl_multimin_fdfminimizer_vector_bfgs2;
    //
    gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc ( T , size );
    // 
    gsl_vector* x = gsl_vector_calloc ( size ) ;
    //
    /// 
    double step = 0  ;
    const DATA& init = m_fit->init() ;
    for ( unsigned int i =0 ; i < x->size ; ++i ) 
    {
      const VE& f = init[i] ;
      gsl_vector_set  ( x , i , 0.5 * f.value()  ) ;
      step = std::max ( step , f.error () ) ;
    }
    //
    gsl_multimin_fdfminimizer_set  ( s , &F , x , 0.1 * step , 0.05 );
    //
    m_iter = 0 ;
    int          status    ;
    do 
    {
      ++m_iter ;
      //
      status = gsl_multimin_fdfminimizer_iterate ( s );
      if ( status ) 
      { 
        gsl_error ( "Chi2Fit error from minimization" , 
                    __FILE__ , __LINE__ , status ) ;
        /// break ;                                                  // BREAK
        continue  ;                                                  // BREAK
      }  
      //
      status = gsl_multimin_test_gradient ( s->gradient , 1.e-4  );
      if ( GSL_SUCCESS == status )
      { 
        // std::cout << " Converge!!! " << std::endl ; 
      }
      //
      // std::cout
      // << " ITER  : " << m_iter 
      //  << " POINT : " << s->f << " " << (*s->x) << std::endl ;  
      // std::cout << " GRAD : "  << (*s->gradient) << std::endl ;
    }
    while ( status == GSL_CONTINUE && m_iter < 5000 ) ; 
    //
    //
    if ( GSL_SUCCESS != status && 27 != status ) { return 300 + status ; }
    //
    
    // get chi2 
    m_chi2 = s->f ;
    //
    // get the solution 
    m_solution = gsl_vector_alloc ( s->x->size ) ;
    gsl_vector_memcpy ( m_solution , s->x ) ;
    // std::cout << " solution! " << (*m_solution) << std::endl ;
    
    // std::cout 
    //<< " #f      " << calls_f      () 
    //<< " #df     " << calls_df     () 
    //<< " #fdf    " << calls_fdf    ()   
    //<< " #fdfddf " << calls_fdfddf ()   
    //<< " #points " << points       ()   
    //<< std::endl ;

    
    //
    /// calcualte hessian:
    //
    m_hessian = gsl_matrix_alloc ( s->x->size , s->x->size );
    fdfddf ( s->x , 0 , m_hessian ) ;
    
    // std::cout << " ANALYTICAL " << std::endl ;
    // std::cout << (*m_hessian)   << std::endl ;
    
    // calculate the covariance  matrix 
    m_covariance = gsl_matrix_alloc ( s->x->size , s->x->size );
    
    Ostap::StatusCode sc = Ostap::Math::GSL::invert_LU_2 ( m_hessian , m_covariance  ) ;
    gsl_matrix_scale ( m_covariance , 2 ) ;
    
    //
    gsl_multimin_fdfminimizer_free ( s ) ;
    gsl_vector_free                ( x ) ;
    //
    return sc ;
  }
  // ==========================================================================
} //                                                  end of anonymous namespace 
// ============================================================================
// constructor 
// ============================================================================
namespace 
{
  Ostap::Math::Chi2Fit::VE
  _adjust_ ( const Ostap::Math::Chi2Fit::DATA& data ,
             Ostap::Math::Chi2Fit::DATA&       comp ,
             const std::size_t                 ncmp )  
  {
    //
    typedef Ostap::Math::Chi2Fit::VE VE ;
    //
    const std::size_t ds = data.size() ;
    const std::size_t cs = comp.size() ;
    //
    if      ( ds > cs ) 
    { comp.insert ( comp.end   () ,    ds - cs   , VE ()      ) ; }
    else if ( cs > ds )   
    { comp.erase  ( comp.begin () +  ( cs - ds ) , comp.end() ) ; }
    //
    VE sumd = Ostap::Math::sum     ( data ) ;
    VE sumc = Ostap::Math::sum     ( comp ) ;
    //
    if ( 0 != sumc.value() ) 
    {
      //
      VE f = sumd/sumc ;
      //
      //
      f.setValue ( f . value () / ncmp        ) ;
      f.setCov2  ( f . cov2  () * ncmp * ncmp ) ;
      //
      return f ; 
    }
    //
    VE asumd = Ostap::Math::abssum ( data ) ;
    VE asumc = Ostap::Math::abssum ( data ) ;
    //
    if ( 0 != asumc ) 
    {
      //
      VE f = 0 < sumd.value() ? asumd / asumc : -asumd / asumc ;  
      //
      f.setValue  ( f . value () / ncmp        ) ;
      f.setCov2   ( f . cov2  () * ncmp * ncmp ) ;
      //
      return f ; 
    }
    // 
    return VE ( 1 , 1 ) ;
  }
}
// ============================================================================
Ostap::Math::Chi2Fit::Chi2Fit 
( const Ostap::Math::Chi2Fit::DATA& data , 
  const Ostap::Math::Chi2Fit::DATA& cmps )
//
  : m_data (     data )
  , m_cmps ( 1 , cmps ) 
  , m_init () 
//
  , m_code ( Ostap::StatusCode::SUCCESS ) 
  , m_solu ( 0 ) 
  , m_cov2 ( 0 ) 
//
  , m_chi2   ( s_Inf )
  , m_calls  ( 0 )
  , m_iters  ( 0 )
  , m_points ( 0 )
{
  //
  // adjust the components and get the init values  
  for ( CMPS::iterator icmp = m_cmps.begin() ; m_cmps.end () != icmp ; ++icmp ) 
  { m_init.push_back ( _adjust_ ( m_data , *icmp , m_cmps.size() ) ) ; }
  //  
  Chi2 c2 ( *this ) ;
  //
  m_code = c2.fit () ;
  //
  if ( m_code.isSuccess() ) 
  {
    m_chi2   = c2.chi2    () ;
    m_calls  =
      c2.calls_f      () + 
      c2.calls_df     () + 
      c2.calls_fdf    () +
      c2.calls_fdfddf () ;
    m_iters  = c2.iters   () ;
    m_points = c2.points  () ;
    //
    if ( 0 != c2.covariance () )
    {
      gsl_matrix* cov2 = gsl_matrix_alloc ( size() , size() ) ;
      gsl_matrix_memcpy ( cov2 , c2.covariance () ) ;
      m_cov2 = cov2 ;
    }
    //
    if ( 0 != c2.solution() ) 
    {
      gsl_vector* solu  = gsl_vector_alloc ( size() ) ;
      gsl_vector_memcpy ( solu , c2.solution () ) ;      
      m_solu = solu ;
    }
    //
  }  
}
// ============================================================================
Ostap::Math::Chi2Fit::Chi2Fit 
( const Ostap::Math::Chi2Fit::DATA& data , 
  const Ostap::Math::Chi2Fit::CMPS& cmps )
  : m_data ( data )
  , m_cmps ( cmps ) 
  , m_init () 
//
  , m_code ( Ostap::StatusCode::SUCCESS ) 
  , m_solu ( 0 ) 
  , m_cov2 ( 0 ) 
//
  , m_chi2   ( s_Inf )
  , m_calls  ( 0 )
  , m_iters  ( 0 )
  , m_points ( 0 )
{
  //
  // adjust the components and get the init values  
  //
  for ( CMPS::iterator icmp = m_cmps.begin() ; m_cmps.end () != icmp ; ++icmp ) 
  { m_init.push_back ( _adjust_ ( m_data , *icmp , m_cmps.size () ) ) ; }
  //
  Chi2 c2 ( *this ) ;
  //
  m_code  = c2.fit   () ;
  //
  if ( m_code.isSuccess() ) 
  {
    m_chi2   = c2.chi2    () ;
    m_calls  = c2.calls_f () + c2.calls_df () + c2.calls_fdf () ;
    m_iters  = c2.iters   () ;
    m_points = c2.points  () ;
    //
    if ( 0 != c2.covariance () )
    {
      gsl_matrix* cov2 = gsl_matrix_alloc ( size() , size() ) ;
      gsl_matrix_memcpy ( cov2 , c2.covariance () ) ;
      m_cov2 = cov2 ;
    }
    //
    if ( 0 != c2.solution() ) 
    {
      gsl_vector* solu  = gsl_vector_alloc ( size() ) ;
      gsl_vector_memcpy ( solu , c2.solution () ) ;      
      m_solu = solu ;
    }
    //
  }  
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Math::Chi2Fit::~Chi2Fit () // destructor 
{
  //
  if ( 0 != m_solu )
  {
    gsl_vector* v = (gsl_vector*) m_solu     ;
    gsl_vector_free ( v ) ;       m_solu = 0 ; 
  }
  //
  if ( 0 != m_cov2 )
  {
    gsl_matrix* m = (gsl_matrix*) m_cov2     ;
    gsl_matrix_free ( m ) ;       m_cov2 = 0 ; 
  }
  //
}
// ============================================================================
// get the parameter 
// ============================================================================
Ostap::Math::Chi2Fit::VE          
Ostap::Math::Chi2Fit::param  ( const unsigned int index ) const 
{
  // 
  if ( m_code.isFailure() ) { return VE ( -s_Inf , -s_Inf ) ; }
  if ( 0 == m_cov2        ) { return VE ( -s_Inf , -s_Inf ) ; }
  if ( 0 == m_solu        ) { return VE ( -s_Inf , -s_Inf ) ; }
  if ( size() <= index    ) { return VE ( -s_Inf , -s_Inf ) ; }
  //
  gsl_vector* v = (gsl_vector*) m_solu ;
  gsl_matrix* m = (gsl_matrix*) m_cov2 ;
  //
  if ( size() != v->size  ) { return VE ( -s_Inf , -s_Inf ) ; }
  if ( size() != m->size1 ) { return VE ( -s_Inf , -s_Inf ) ; }
  if ( size() != m->size2 ) { return VE ( -s_Inf , -s_Inf ) ; }
  //
  return VE ( gsl_vector_get ( v , index         ) , 
              gsl_matrix_get ( m , index , index ) ) ;
  //
}
// ============================================================================
// get all parameters at once 
// ============================================================================
Ostap::Math::Chi2Fit::DATA
Ostap::Math::Chi2Fit::params () const            // get all parameters at once
{
  Ostap::Math::Chi2Fit::DATA pars ;
  //
  if ( m_code.isFailure() ) { return pars ; }
  if ( 0 == m_cov2        ) { return pars ; }
  if ( 0 == m_solu        ) { return pars ; }
  //
  gsl_vector* v = (gsl_vector*) m_solu ;
  gsl_matrix* m = (gsl_matrix*) m_cov2 ;
  //
  if ( size() != v->size  ) { return pars ; }
  if ( size() != m->size1 ) { return pars ; }
  if ( size() != m->size2 ) { return pars ; }
  //
  for ( unsigned int i = 0 ; i < size() ; ++i ) 
  { pars.push_back ( param ( i ) ) ; }
  //
  return pars ;
}
// ============================================================================
// get the covariance matrix elements
// ============================================================================
double Ostap::Math::Chi2Fit::cov2   ( const unsigned int i1 , 
                                      const unsigned int i2 ) const 
{
  // 
  if ( m_code.isFailure() ) { return -s_Inf ; }
  if ( 0 == m_cov2        ) { return -s_Inf ; }
  if ( 0 == m_solu        ) { return -s_Inf ; }
  if ( size() <= i1       ) { return -s_Inf ; }
  if ( size() <= i2       ) { return -s_Inf ; }
  //
  gsl_vector* v = (gsl_vector*) m_solu ;
  gsl_matrix* m = (gsl_matrix*) m_cov2 ;
  //
  if ( size() != v->size  ) { return -s_Inf ; }
  if ( size() != m->size1 ) { return -s_Inf ; }
  if ( size() != m->size2 ) { return -s_Inf ; }
  //
  return gsl_matrix_get ( m , i1 , i2 )  ;
  //
}
// ============================================================================
// output 
// ============================================================================
std::ostream& 
Ostap::Math::Chi2Fit::fillStream ( std::ostream& s ) const 
{
  s << "  Status    : " << m_code             << std::endl
    << "  Chi2      : " << chi2   ()          << std::endl
    << "  #Calls    : " << ncalls ()          << std::endl
    << "  #iters    : " << niters ()          << std::endl
    << "  #points   : " << points ()          << std::endl
    << "  #size     : " << size   ()          << std::endl
    << "  #dof      : " << points () - size() << std::endl ;
  
  if ( m_code.isSuccess () && 
       0 !=  m_solu        && 
       0 !=  m_cov2          )
  {
    gsl_vector* v = (gsl_vector*) m_solu ;
    gsl_matrix* m = (gsl_matrix*) m_cov2 ;
    //
    if ( v->size == m->size1 && m->size1 == m->size2 ) 
    {
      for  (unsigned int i = 0 ; i < size() ; ++i ) 
      {
        s << "  Parameter : " << param(i) << std::endl ;
      } 
      s << "  Solution  : "   << (*v ) << std::endl 
        << "  Covariance: \n" << (*m ) << std::endl;
      //
    }
  }
  return s ;
}
// ============================================================================
// output 
// ============================================================================
std::string Ostap::Math::Chi2Fit::toString () const 
{
  std::ostringstream s;
  fillStream ( s ) ;
  return s.str();
}
// ============================================================================
// The END 
// ============================================================================
