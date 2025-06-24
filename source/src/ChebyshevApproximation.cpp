// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <algorithm>
// ============================================================================
// GSL
// ============================================================================
#include "gsl/gsl_math.h"
#include "gsl/gsl_chebyshev.h"
// ============================================================================
// ROOT 
// ============================================================================
#include "Math/GSLFunctionAdapter.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/ChebyshevApproximation.h"
#include "Ostap/PyCallable.h"
#include "Ostap/Polynomials.h"
#include "Ostap/StatusCode.h"
// ============================================================================
// Local
// ============================================================================
#include "local_math.h"
#include "status_codes.h"
// ============================================================================
/** @file
 *  Implementation file for class Ostap::Math::ChebyshevApproximation
 *  @date 2019-09-25 
 *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  const char s_ERROR   [] = "Invalid gsl_cheb_series object!"                 ;
  const char s_METHOD1 [] = "Ostap::Math::ChebyshevApproximation::evaluate"   ;
  const char s_METHOD2 [] = "Ostap::Math::ChebyshevApproximation::derivative" ;
  const char s_METHOD3 [] = "Ostap::Math::ChebyshevApproximation::integral"   ;
  const char s_METHOD4 [] = "Ostap::Math::ChebyshevApproximation::operator+"  ;
  const char s_METHOD5 [] = "Ostap::Math::ChebyshevApproximation::operator*"  ;
  const char s_METHOD6 [] = "Ostap::Math::ChebyshevApproximation::polynomial" ;
  const Ostap::StatusCode s_SC =  Ostap::StatusCode::FAILURE                  ;
  // ==========================================================================
}
// ============================================================================
/*  constructor from the function, low/high-limits and the approximation order 
 *  @param func the function 
 *  @param a the low-limit 
 *  @param b the high-limit 
 *  @param N the approximation order 
 */
// ============================================================================
Ostap::Math::ChebyshevApproximation::ChebyshevApproximation
( std::function<double(double)> func , 
  const double                  a    , 
  const double                  b    , 
  const unsigned short          N    ) 
  : m_a ( std::min ( a , b ) )
  , m_b ( std::max ( a , b ) )
  , m_N ( N )
  , m_chebyshev ( nullptr )
{
  //
  ROOT::Math::GSLFunctionAdapter< std::function<double(double)> > adapter ;
  const void* p = &func ;
  //
  gsl_function F ;
  F.function = &adapter.F ;
  F.params   = const_cast<void*> (p) ;
  //
  gsl_cheb_series* ns = gsl_cheb_alloc ( m_N ) ;  
  gsl_cheb_init ( ns , &F , m_a , m_b ) ;
  //
  m_chebyshev = (char*) ns ;
}
// ============================================================================
/* constructor from the function, low/high-limits and the approximation order 
 *  @param func the function 
 *  @param a the low-limit 
 *  @param b the high-limit 
 *  @param N the approximation order 
 */
// ============================================================================
Ostap::Math::ChebyshevApproximation::ChebyshevApproximation
( const Ostap::Functions::PyCallable& func , 
  const double             a    , 
  const double             b    , 
  const unsigned short     N    ) 
  : ChebyshevApproximation ( std::function<double(double)> ( std::cref ( func ) ) , a ,  b , N )
{}
// ============================================================================

// ============================================================================
// default (protected) constructor 
// ============================================================================
Ostap::Math::ChebyshevApproximation::ChebyshevApproximation()
  : m_a ( 0 ) 
  , m_b ( 1 ) 
  , m_N ( 0 ) 
  , m_chebyshev ( nullptr ) 
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Math::ChebyshevApproximation::ChebyshevApproximation 
( const Ostap::Math::ChebyshevApproximation& right ) 
  : m_a         ( right.m_a )
  , m_b         ( right.m_b )
  , m_N         ( right.m_N )
  , m_chebyshev ( nullptr   )
{
  gsl_cheb_series* ns = gsl_cheb_alloc ( m_N ) ;
  gsl_cheb_series* os = (gsl_cheb_series*) right.m_chebyshev ;
  //
  ns -> a        = os -> a        ;
  ns -> b        = os -> b        ;
  ns -> order    = os -> order    ;
  ns -> order_sp = os -> order_sp ;
  //
  // copy the coefficients and the function values 
  //
  std::copy ( os -> c , os-> c + ( m_N + 1 ) , ns -> c ) ;
  std::copy ( os -> f , os-> f + ( m_N + 1 ) , ns -> f ) ;
  //
  m_chebyshev = (char*) ns ;
}
// ============================================================================
// move constructor 
// ============================================================================
Ostap::Math::ChebyshevApproximation::ChebyshevApproximation 
( Ostap::Math::ChebyshevApproximation&& right ) 
  : m_a         ( right.m_a )
  , m_b         ( right.m_b )
  , m_N         ( right.m_N )
  , m_chebyshev ( right.m_chebyshev )
{
  right.m_chebyshev = nullptr ;
}
// ============================================================================
//  destructor
// ============================================================================
Ostap::Math::ChebyshevApproximation::~ChebyshevApproximation() 
{
  if ( m_chebyshev ) 
  {
    gsl_cheb_series* cs = (gsl_cheb_series*) m_chebyshev ;
    gsl_cheb_free ( cs )  ;
    m_chebyshev = nullptr ;
  } 
}
// ============================================================================
// the main method: evaluate the approximation sum 
// ============================================================================
double Ostap::Math::ChebyshevApproximation::evaluate 
( const double x ) const 
{
  //
  Ostap::Assert ( m_chebyshev ,
                  s_ERROR     ,
                  s_METHOD6   ,
                  INVALID_CHEBYSHEV , __FILE__ , __LINE__  ) ;
  //
  Ostap::Assert ( m_chebyshev , s_ERROR , s_METHOD1 , s_SC ) ;
  //
  gsl_cheb_series* cs = (gsl_cheb_series*) m_chebyshev ;
  return x < m_a ? 0.0 : x > m_b ? 0.0 : gsl_cheb_eval ( cs , x ) ;
}
// ============================================================================
/*  the main method: evaluate the approximation sum, 
 *  using at most <code>n</code> terms 
 */
// ============================================================================
double Ostap::Math::ChebyshevApproximation::evaluate 
( const double         x , 
  const unsigned short n ) const 
{
  //
  Ostap::Assert ( m_chebyshev ,
                  s_ERROR     ,
                  s_METHOD6   ,
                  INVALID_CHEBYSHEV , __FILE__ , __LINE__  ) ;
  //
  gsl_cheb_series* cs = (gsl_cheb_series*) m_chebyshev ;
  return x < m_a ? 0.0 : x > m_b ? 0.0 : gsl_cheb_eval_n ( cs , n , x ) ;
}
// ============================================================================
/*  the main method: evaluate the approximation sum 
 *  @return the approximation with the error estimate 
 */
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ChebyshevApproximation::eval_err 
( const double         x ) const 
{
  //
  Ostap::Assert ( m_chebyshev ,
                  s_ERROR     ,
                  s_METHOD6   ,
                  INVALID_CHEBYSHEV , __FILE__ , __LINE__  ) ;
  //
  if ( x < m_a || x > m_b ) { return 0 ; }
  //
  gsl_cheb_series* cs = (gsl_cheb_series*) m_chebyshev ;
  //
  double result =  0 ;
  double error  =  0 ;
  //
  gsl_cheb_eval_err ( cs , x , &result , &error ) ;
  //
  return Ostap::Math::ValueWithError ( result , error * error ) ;
}
// ============================================================================
/*  the main method: evaluate the approximation sum
 *  using at most <code>n</code> terms 
 *  @retutn the approximation with the error estimate 
 */
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ChebyshevApproximation::eval_err 
( const double         x ,
  const unsigned short n ) const 
{
  Ostap::Assert ( m_chebyshev ,
                  s_ERROR     ,
                  s_METHOD6   ,
                  INVALID_CHEBYSHEV , __FILE__ , __LINE__  ) ;
  //
  if ( x < m_a || x > m_b ) { return 0 ; }
  //
  gsl_cheb_series* cs = (gsl_cheb_series*) m_chebyshev ;
  //
  double result =  0 ;
  double error  =  0 ;
  //
  gsl_cheb_eval_n_err ( cs , n , x , &result , &error ) ;
  //
  return Ostap::Math::ValueWithError ( result , error * error ) ;
}
// ============================================================================
// get a derivative  
// ============================================================================
Ostap::Math::ChebyshevApproximation
Ostap::Math::ChebyshevApproximation::derivative () const 
{
  //
  Ostap::Assert ( m_chebyshev ,
                  s_ERROR     ,
                  s_METHOD6   ,
                  INVALID_CHEBYSHEV , __FILE__ , __LINE__  ) ;
  //
  gsl_cheb_series*       ds = gsl_cheb_alloc ( m_N ) ;
  const gsl_cheb_series* cs = (const gsl_cheb_series*) m_chebyshev ;
  //
  gsl_cheb_calc_deriv ( ds , cs ) ;
  //
  ChebyshevApproximation deriv ;
  deriv.m_a         = m_a ;
  deriv.m_b         = m_b ;
  deriv.m_N         = m_N ;
  deriv.m_chebyshev = (char*) ds ;
  //
  return deriv ;
}
// ============================================================================
// get an integral: \f$ F(x) \equiv \int_a^{z} f(t) \deriv t  + C \f$ 
// ============================================================================
Ostap::Math::ChebyshevApproximation
Ostap::Math::ChebyshevApproximation::integral 
( const double C ) const 
{
  //
  Ostap::Assert ( m_chebyshev ,
                  s_ERROR     ,
                  s_METHOD6   ,
                  INVALID_CHEBYSHEV , __FILE__ , __LINE__  ) ;
  //
  gsl_cheb_series*       ds = gsl_cheb_alloc ( m_N ) ;
  const gsl_cheb_series* cs = (const gsl_cheb_series*) m_chebyshev ;
  //
  gsl_cheb_calc_integ ( ds , cs ) ;
  //
  ds->c[0] += 2*C ;
  //
  ChebyshevApproximation integ ;
  integ.m_a         = m_a ;
  integ.m_b         = m_b ;
  integ.m_N         = m_N ;
  integ.m_chebyshev = (char*) ds ;
  //
  return integ ;
}
// ============================================================================
// shift by a constant 
// ============================================================================
Ostap::Math::ChebyshevApproximation&
Ostap::Math::ChebyshevApproximation::operator+= ( const double a )
{
  Ostap::Assert ( m_chebyshev ,
                  s_ERROR     ,
                  s_METHOD6   ,
                  INVALID_CHEBYSHEV , __FILE__ , __LINE__  ) ;
  //
  gsl_cheb_series* cs = (gsl_cheb_series*) m_chebyshev ;
  //
  cs -> c[0] += 2 * a ;
  //
  return *this ;
}
// ============================================================================
// scale by a constant 
// ============================================================================
Ostap::Math::ChebyshevApproximation&
Ostap::Math::ChebyshevApproximation::operator*= ( const double a )
{
  Ostap::Assert ( m_chebyshev ,
                  s_ERROR     ,
                  s_METHOD6   ,
                  INVALID_CHEBYSHEV , __FILE__ , __LINE__  ) ;
  //
  gsl_cheb_series* cs = (gsl_cheb_series*) m_chebyshev ;
  std::transform ( cs->c                 ,
                   cs->c + cs->order + 1 ,
                   cs->c                 ,
                   [a]( double v ) -> double { return v * a ; } );
  //
  return *this ;
}
// ============================================================================
/**  convert it to pure chebyshev sum 
 *  suppressing the coefficients that are small enough 
 *   - it is numerically zero: \f$ c_k \approx 0 \f% 
 *   - or if epsilon  > 0: \f$ \left| c_k \right| \le \epsilon \f$ 
 *   - or if scale   != 0: \f$ \left| s \right| + \left| c_k \right| \approx \left| s \right| \f$        
 */
 // ============================================================================
Ostap::Math::ChebyshevSum
Ostap::Math::ChebyshevApproximation::polynomial
( const double epsilon ,
  const double scale   ) const
{
  //
  Ostap::Assert ( m_chebyshev ,
                  s_ERROR     ,
                  s_METHOD6   ,
                  INVALID_CHEBYSHEV , __FILE__ , __LINE__  ) ;
  //
  gsl_cheb_series* cs = (gsl_cheb_series*) m_chebyshev ;
  Ostap::Math::ChebyshevSum cp { cs->c , cs->c + cs->order + 1 , cs->a , cs->b } ;
  cp.setPar ( 0 , 0.5 * cp.par ( 0 ) ) ;
  //
  cp.remove_noise ( epsilon , scale ) ; 
  return cp ;
}
// ============================================================================
/** Get sum of all small/neglected terms 
 *  Coefficients are small if 
 *   - it is numerically zero: \f$ c_k \approx 0 \f% 
 *   - or if epsilon  > 0: \f$ \left| c_k \right| \le \epsilon \f$ 
 *   - or if scale   != 0: \f$ \left| s \right| + \left| c_k \right| \approx \left| s \right| \f$        
 */
// ============================================================================
Ostap::Math::ChebyshevSum
Ostap::Math::ChebyshevApproximation::noise 
( const double epsilon ,
  const double scale   ) const
{
  //
  Ostap::Assert ( m_chebyshev ,
                  s_ERROR     ,
                  s_METHOD6   ,
                  INVALID_CHEBYSHEV , __FILE__ , __LINE__  ) ;
  //
  gsl_cheb_series* cs = (gsl_cheb_series*) m_chebyshev ;
  Ostap::Math::ChebyshevSum cp { cs->c , cs->c + cs->order + 1 , cs->a , cs->b } ;
  cp.setPar ( 0 , 0.5 * cp.par ( 0 ) ) ;
  //
  const bool   eps    = 0 < epsilon        ;
  const bool   sca    = scale              ;
  const double ascale = std::abs ( scale ) ;
  //
  const std::size_t NN { cp.npars() } ;
  for ( std::size_t k = 0 ; k < NN ; ++k )
    {
      const double absp = std::abs ( cp.par ( k ) ) ;
      //
      if      ( s_zero ( absp )                           ) {} 
      else if ( eps && absp <= epsilon                    ) {} 
      else if ( sca && s_equal ( ascale + absp , ascale ) ) {}
      else    { cp.setPar ( k , 0.0 , true ) ;  }
    }
  //
  return cp ;
}

// ============================================================================
//                                                                      The END 
// ============================================================================
