// ============================================================================
// Incldue files 
// ============================================================================
// STD&STL
// ============================================================================
#include <memory>
#include <map> 
#include <cstring>
#include <numeric>
#include <algorithm>
#include <tuple>
#include <array>
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/Hash.h"
#include "Ostap/ECDF.h"
#include "Ostap/Power.h"
#include "Ostap/MoreMath.h"
#include "Ostap/StatEntity.h"
#include "Ostap/WStatEntity.h"
#include "Ostap/Moments.h"
// ============================================================================
// Local
// ============================================================================
#include "Exception.h"
#include "status_codes.h"
#include "local_math.h"
#include "local_hash.h"
#include "syncedcache.h"
// ============================================================================
/** @file 
 *  Implementation file for class Ostap::Math::ECDF
 *  @see Ostap::Math::ECDF
 *  @date 2024-09-16 
 *  @author Vanya BELYAEV 
 */
// ============================================================================
namespace
{
  // ==========================================================================
  /// actual kernel function type 
  typedef double (*FUNPTR)(double) ; 
  typedef std::tuple<
    Ostap::Math::DensityEstimator::Kernel,
    FUNPTR,
    double, 
    double, 
    double, 
    double> RECORD ;
  //
  typedef std::array<RECORD,Ostap::Math::DensityEstimator::Last+1> TABLE ;  
  const TABLE s_TABLE { 
    RECORD ( Ostap::Math::DensityEstimator::Uniform,     &Ostap::Math::k_uniform     ,1,1./3   ,1./2                   ,0.929) ,
    RECORD ( Ostap::Math::DensityEstimator::Triangular,  &Ostap::Math::k_triangular  ,1,1./6   ,2./3                   ,0.986) ,
    RECORD ( Ostap::Math::DensityEstimator::Epanechnikov,&Ostap::Math::k_epanechnikov,1,1./5   ,3./5                   ,1.000) ,
    RECORD ( Ostap::Math::DensityEstimator::Quartic,     &Ostap::Math::k_quartic     ,1,1./7   ,5./7                   ,0.994) ,
    RECORD ( Ostap::Math::DensityEstimator::Triweight,   &Ostap::Math::k_triweight   ,1,1./9   ,350./429               ,0.987) ,
    RECORD ( Ostap::Math::DensityEstimator::Tricube,     &Ostap::Math::k_tricube     ,1,35./243,175./447               ,0.998) ,
    RECORD ( Ostap::Math::DensityEstimator::Gaussian,    &Ostap::Math::k_gaussian    ,0,1      ,0.5/std::sqrt(M_PI)    ,0.951) ,
    RECORD ( Ostap::Math::DensityEstimator::Cosine,      &Ostap::Math::k_cosine      ,1,1-8/(M_PI*M_PI) , M_PI*M_PI/16 ,0.999) ,
    RECORD ( Ostap::Math::DensityEstimator::Logistic,    &Ostap::Math::k_logistic    ,0,M_PI*M_PI/3 , 1.6              ,0.887) ,
    RECORD ( Ostap::Math::DensityEstimator::Sigmoid,     &Ostap::Math::k_sigmoid     ,0,M_PI*M_PI/4 , 2/(M_PI*M_PI)    ,0.843) ,
  } ;
  // ==========================================================================
  /// cached version of normalized incomplete beta-function
  double beta_inc  
  ( const double alpha      ,  
    const double beta       , 
    const double z          ,
    std::size_t  N  = 20000 ) 
  {
   /// cache for normalized incomplete beta-function;
   typedef std::map<std::size_t,double>  MAP   ;
   typedef SyncedCache<MAP>              CACHE ;
   /// the actual integration cache
   static CACHE   s_cache {}      ; // integration cache
   ///
   const std::size_t key { Ostap::Utils::hash_combiner ( z , alpha , beta ) }; 
   ///
   { 
    CACHE::Lock lock { s_cache.mutex () } ;
    auto it = s_cache->find ( key ) ;
    if ( s_cache->end () != it ) { return it->second ; } 
   }
   // calcluate the result 
   const double result = Ostap::Math::beta_inc ( alpha , beta , z ) ;
   {
    CACHE::Lock lock { s_cache.mutex () } ;
    if ( N < s_cache->size() ) { s_cache->clear() ; }
    s_cache->insert ( std::make_pair ( key , result ) ) ; 
   }
   //
   return result ;
  }
  /// cached version of normalized incomplete beta-function
  double beta_log  
  ( const double alpha      ,  
    const double beta       ,
    const double t          , 
    const double dt         ,
    std::size_t  N  = 20000 ) 
  {
   /// cache for normalized incomplete beta-function;
   typedef std::map<std::size_t,double>  MAP   ;
   typedef SyncedCache<MAP>              CACHE ;
   /// the actual integration cache
   static CACHE   s_cache {}      ; // integration cache
   ///
   const std::size_t key { Ostap::Utils::hash_combiner ( t , alpha , beta , t , dt ) }; 
   ///
   { 
    CACHE::Lock lock { s_cache.mutex () } ;
    auto it = s_cache->find ( key ) ;
    if ( s_cache->end () != it ) { return it->second ; } 
   }
  // calcluate the result
  double result  = std::log  ( dt ) ; 
  result        += ( alpha - 1 ) * std::log ( t       ) ;  
  result        += ( beta  - 1 ) * std::log ( 1 - t   ) ; 
  result        -= Ostap::Math::lnbeta ( alpha , beta ) ;
  result         = std::exp ( result ) ;   
  {
    CACHE::Lock lock { s_cache.mutex () } ;
    if ( N < s_cache->size() ) { s_cache->clear() ; }
    s_cache->insert ( std::make_pair ( key , result ) ) ; 
   }
   //
   return result ;
  }
}
// ============================================================================
double Ostap::Math::DensityEstimator::kernel 
( const double                                u , 
  const Ostap::Math::DensityEstimator::Kernel k )
{
  //
  Ostap::Assert ( Uniform <= k && k <=  Last , 
		  "Invalid Kernel!" , 
		  "Ostap::Math::DensityEstimator::kernel"  ,
		  INVALID_KERNEL , __FILE__, __LINE__ ) ;
  //
  return 
    Gaussian     == k ? Ostap::Math::k_gaussian     ( u ) :
    Epanechnikov == k ? Ostap::Math::k_epanechnikov ( u ) :
    Triangular   == k ? Ostap::Math::k_triangular   ( u ) :
    Quartic      == k ? Ostap::Math::k_quartic      ( u ) : 
    Triweight    == k ? Ostap::Math::k_triweight    ( u ) :
    Tricube      == k ? Ostap::Math::k_tricube      ( u ) : 
    Cosine       == k ? Ostap::Math::k_cosine       ( u ) : 
    Logistic     == k ? Ostap::Math::k_logistic     ( u ) : 
    Sigmoid      == k ? Ostap::Math::k_sigmoid      ( u ) : Ostap::Math::k_uniform ( u ) ;  
}
// ============================================================================

// ============================================================================
// Cosine kernel 
// ============================================================================
double Ostap::Math::k_cosine       ( const double u ) 
{ return std::abs ( u ) <= 1 ? 0.25 * M_PI * std::cos ( 0.5 * M_PI * u ) : 0.0 ;  } 
// ============================================================================
// Logistic Kernel
// ============================================================================
double Ostap::Math::k_logistic     ( const double u )
{ return std::abs ( u ) < s_EXP_OVERFLOW ? 1 / ( 2.0 * std::cosh ( u ) + 2 ) :  0.0 ;  }
// ============================================================================
// sigmoid kernel 
// ============================================================================
double Ostap::Math::k_sigmoid      ( const double u ) 
{ return std::abs ( u ) < s_EXP_OVERFLOW ? s_PIi / std::cosh ( u ) : 0.0 ; }
// ============================================================================
// Gaussian kernel 
// ============================================================================
double Ostap::Math::k_gaussian    ( const double u ) { return gauss_pdf ( u ) ; }
// ============================================================================

// ============================================================================
// get the "optimal" value for the smoothing parameter 
// ============================================================================
double Ostap::Math::DensityEstimator::hopt
( const  Ostap::Math::ECDF& data )
{
  if ( 2 > data.N() ) { return -1 ; }
  //
  Ostap::StatEntity cnt {} ;
  for ( const auto v : data ) { cnt += v ; } 
  //
  double value = cnt.rms() ;
  //
  if ( 4 <= data.N () )
    {
      // the first * third  quartiles
      const std::size_t i1 = static_cast<std::size_t> ( 0.25 * data.N () ) ;
      const std::size_t i3 = static_cast<std::size_t> ( 0.75 * data.N () ) ;
      const double x1 = data [ i1 ] ;
      const double x3 = data [ i3 ] ;
      if ( x1 < x3 ) { value = std::min ( value , ( x3 - x1 ) / 1.34 ) ;} 
    }
  //
  return 0.9 * value * std::pow ( data.N() , -0.2 );
}
// ============================================================================
// get the "optimal" value for the smoothing parameter 
// ============================================================================
double Ostap::Math::DensityEstimator::hopt
( const  Ostap::Math::WECDF& data )
{
  if ( 2 > data.N() || 2 > data.nEff() ) { return -1 ; }
  //
  Ostap::WStatEntity cnt {} ;
  for ( const auto v : data ) { cnt.add ( v.first , v.second ) ; } 
  //
  double value = cnt.rms() ;
  //
  if ( 4 <= data.N () && 4 <= data.nEff() )
    {      
      // the first * third  quartiles, well... not really... @todo FIX ME!
      const std::size_t i1 = static_cast<std::size_t> ( 0.25 * data.N () ) ;
      const std::size_t i3 = static_cast<std::size_t> ( 0.75 * data.N () ) ;
      const double x1 = data [ i1 ].first ;
      const double x3 = data [ i3 ].first ;
      if ( x1 < x3 ) { value = std::min ( value , ( x3 - x1 ) / 1.34 ) ;} 
    }
  //
  return 0.9 * value * std::pow ( data.nEff () , -0.2 );
}
// ============================================================================
// Standard constructor from  data
// ============================================================================
Ostap::Math::ECDF::ECDF
( const Ostap::Math::ECDF::Data&  data          ,
  const bool                      complementary )
  : Ostap::Math::ECDF::ECDF ( data.begin() , data.end() , complementary )
{}
// ============================================================================
// constructor to create complementary/oridnary ECDF
// ============================================================================
Ostap::Math::ECDF::ECDF
( const Ostap::Math::ECDF & right         , 
  const bool                complementary ) 
  : ECDF ( right ) 
{
  m_complementary = complementary ; 
}
// ============================================================================

// ============================================================================
// copy assignement 
// ============================================================================
Ostap::Math::ECDF&
Ostap::Math::ECDF::operator=
( const Ostap::Math::ECDF&  right )
{
  if ( this == &right ) { return *this ; }
  m_data          = right.m_data          ;
  m_complementary = right.m_complementary ;
  return *this ;
}
// ============================================================================
// move assignement 
// ============================================================================
Ostap::Math::ECDF&
Ostap::Math::ECDF::operator=
( Ostap::Math::ECDF&&  right )
{
  if ( this == &right ) { return *this ; }
  this->swap ( right ) ;
  return *this ; 
}
// ============================================================================
// check that ECDF is OK: there are some entries 
// ============================================================================
Ostap::Math::ECDF&
Ostap::Math::ECDF::check_me ()
{
  Data::iterator remove = std::remove_if
    ( m_data.begin () ,
      m_data.end   () ,
      [] ( const double x ) -> bool
      { return !std::isfinite ( x ) ;} ) ;
  //
  m_data.erase ( remove , m_data.end() ) ; 
  //
  Ostap::Assert ( !m_data.empty()                     ,
                  "No data for Empirical CDF"         ,
                  "Ostap::Math::CDF"                  ,
                  INVALID_DATA   , __FILE__ ,__LINE__ ) ;
  //
  return *this ;
}
// ============================================================================
// swap two objects 
// ============================================================================
void Ostap::Math::ECDF::swap ( Ostap::Math::ECDF& right )
{
  std::swap ( m_data          , right.m_data          ) ;
  std::swap ( m_complementary , right.m_complementary ) ;
}
// ============================================================================
// the main method 
// ============================================================================
double Ostap::Math::ECDF::evaluate   ( const double x ) const
{
  if      ( x < m_data.front () ) { return m_complementary ? 1.0 : 0.0 ; } 
  else if ( x > m_data.back  () ) { return m_complementary ? 0.0 : 1.0 ; } 
  //
  const double result = double ( rank ( x ) ) / m_data.size () ;
  return m_complementary ? ( 1 - result ) : result ; 
}
// ============================================================================
// the main method 
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ECDF::estimate ( const double x ) const
{
  const std::size_t NN = m_data.size() ;
  //
  if      ( x < m_data.front () ) { return Ostap::Math::binomEff ( m_complementary ? NN : 0u , NN ) ; }
  else if ( x > m_data.back  () ) { return Ostap::Math::binomEff ( m_complementary ? 0u : NN , NN ) ; }
  //
  const std::size_t success  =
    std::upper_bound ( m_data.begin () , m_data.end () , x ) - m_data.begin() ;
  //
  return Ostap::Math::binomEff ( m_complementary ? NN - success : success , NN ) ;
}
// ============================================================================
// add a value to data container  
// ============================================================================
Ostap::Math::ECDF& 
Ostap::Math::ECDF::add
( const double value  )
{
  if ( !std::isfinite ( value ) ) { return *this ; } 
  auto where = std::upper_bound ( m_data.begin () , m_data.end   () , value ) ;
  m_data.insert ( where , value ) ;
  return *this ;
}
// ============================================================================
// add values to data container  
// ============================================================================
Ostap::Math::ECDF&
Ostap::Math::ECDF::add
( const Ostap::Math::ECDF& values )
{
  /// prepare the output 
  Data tmp  ( values.size() + m_data.size () ) ;
  /// merge two sorted containers 
  std::merge ( values.m_data.begin () ,                     
               values.m_data.end   () ,
               m_data.begin        () ,
               m_data.end          () ,                     
               tmp.begin           () ) ;
  /// swap the merged result  with own data 
  std::swap ( m_data , tmp ) ;
  return *this ;
}
// ============================================================================
// add a value to data container  
// ============================================================================
Ostap::Math::ECDF& 
Ostap::Math::ECDF::add
( const Ostap::Math::ECDF::Data& values )
{ return add ( values.begin() , values.end() ) ; }
// ============================================================================
Ostap::Math::ECDF
Ostap::Math::ECDF::__add__  ( const double                    x ) const 
{ ECDF c { *this } ; c.add ( x ) ; return c ; }
// ============================================================================
Ostap::Math::ECDF
Ostap::Math::ECDF::__add__  ( const Ostap::Math::ECDF&        x ) const 
{ ECDF c { *this } ; c.add ( x ) ; return c ; }
// ============================================================================
Ostap::Math::ECDF
Ostap::Math::ECDF::__add__  ( const Ostap::Math::ECDF::Data&  x ) const 
{ ECDF c { *this } ; c.add ( x ) ; return c ; }
// ============================================================================
// ============================================================================
// get ranks of the elements in the (pooled) sample
// ============================================================================
Ostap::Math::ECDF::Indices
Ostap::Math::ECDF::ranks
( const Ostap::Math::ECDF& sample ) const
{
  const Data::size_type N  = size() ;
  // fill output array with N
  Indices result ( sample.size () ,  N ) ;
  Data::size_type NS= sample.size() ;
  for ( Data::size_type i = 0 ; i < NS ; ++i )
    {
      Data::size_type r = rank ( sample.m_data [ i ] ) ;
      result [ i ] = r ;
      // try to be a bit more efficient, the rest of array is filled with N 
      if ( N <= r ) { break ; }
    }
  return result ;
}
// ============================================================================
/*  assuming that x comes from the same distribution 
 *  return transformed value  \f$ g = f(x) \f$, such that 
 *  \f$ g \f$  has Gaussian distribution
 */
// ============================================================================
double Ostap::Math::ECDF::gauss   ( const double x ) const
{ return Ostap::Math::probit ( uniform ( x ) ) ; }
// ============================================================================
/*  assuming that x comes from the same distribution 
 *  return transformed value  \f$ u = f(x) \f$, such that 
 *  \f$ u \f$  has uniform distribution for \f$ 0 \le  u \le 1 \f$ 
 */
// ============================================================================
double Ostap::Math::ECDF::uniform ( const double x ) const
{
  return 
    ( ( x < xmin () ) ? (       1.0 / size () ) :
      ( x > xmax () ) ? ( 1.0 - 1.0 / size () ) : 
      ( std::lower_bound ( m_data.begin () , m_data.end () , x ) - m_data.begin() ) * 1.0 / size () ) ;
}
// =============================================================================
/* get p-quantile of distributtion: \f$ 1 \le p \le1  \f$
 *  @param p      (INPUT) quantile
 *  @param alphap (INPUT) parameter alphap \f$ 0 \le\alpha_p \le 1 \f$ 
 *  @param abetap (INPUT) parameter betap \f$ 0 \le\beta_p \le 1 \f$ 
*/    
// =============================================================================
double Ostap::Math::ECDF::quantile
( const double p      ,  
  const double alphap , 
  const double betap  ) const 
{
  if      ( !p     || s_zero  ( p )     ) { return xmin () ; }
  else if ( 1 == p || s_equal ( p , 1 ) ) { return xmax () ; } 
  //
  Ostap::Assert ( 0 <= p && p <= 1 , 
		  "Invalid quantile!" , 
		  "Ostap::Math::ECDF::quantile" , 
		  INVALID_QUANTILE , __FILE__  , __LINE__ ) ;
  Ostap::Assert ( 0 <= alphap && alphap <= 1 , 
		  "Invalid alphap!" , 
		  "Ostap::Math::ECDF::quantile" , 
		  INVALID_QUANTILE , __FILE__  , __LINE__ ) ;
  Ostap::Assert ( 0 <= betap && betap <= 1 , 
		  "Invalid betap!" , 
		  "Ostap::Math::ECDF::quantile" , 
		  INVALID_QUANTILE , __FILE__  , __LINE__ ) ;
  //
  const double          m = alphap + p * ( 1 - alphap - betap ) ;
  const Data::size_type n = N () ; 
  const double          a = p * n  + m ; 
  const int             j = static_cast<int> ( std::floor ( a ) ) ; 
  //
  if       ( j <  0     ) { return m_data.front() ; } 
  else if  ( n <= j + 1 ) { return m_data.back () ; } 
  //
  const double          g = a - j ;
  return ( 1 - g ) * m_data [ j ] + g * m_data [ j + 1 ] ; 
}
// ============================================================================
/* get p-quantile
 *  @see https://arxiv.org/abs/2304.07265
 *  @see Andrey Akinshin, "Weighted quantile estimators", arXiv:2304.07265     
 */
// ============================================================================
double Ostap::Math::ECDF::quantile_HF
( const double                   p ,
  const Ostap::Math::ECDF::QType t )
{
  if      ( !p     || s_zero  ( p )     ) { return xmin () ; }
  else if ( 1 == p || s_equal ( p , 1 ) ) { return xmax () ; } 
  //
  Ostap::Assert ( 0 <= p && p <= 1 , 
		  "Invalid p-quantile!" , 
		  "Ostap::Math::ECDF::quantile_HF" , 
		  INVALID_QUANTILE , __FILE__  , __LINE__ ) ;
  Ostap::Assert ( One <= t && t <= Nine , 
		  "Invalid quantile estimator type!" , 
		  "Ostap::Math::ECDF::quantile_HF" , 
		  INVALID_QUANTILE , __FILE__  , __LINE__ ) ;
  //
  const Data::size_type n = N ()  ;
  //
  auto index = [ n ] ( const double h ) -> std::size_t
  {
    const long j = Ostap::Math::round ( h ) ;
    return  ( j < 0 ) ? 0 :  ( n <= j ) ? ( n - 1 ) : std::size_t ( j ) ;  
  } ;
  //
  auto QF = [this,index] ( const double h ) ->double
  {
    const std::size_t jf = index ( std::floor ( h ) ) ;
    const std::size_t jc = index ( std::ceil  ( h ) ) ;
    const double xf = this->data ( jf ) ;
    const double xc = this->data ( jc ) ;
    //
    return xf + ( h - jf ) * ( xc - xf ) ;
  } ;
  //
  switch ( t ) 
    {
    case One :
      {
	const double h = n * p ;
	const std::size_t j = index ( std::ceil ( h ) ) ; 
	return data ( j ) ;
      }
    case Two :
      {
	const double h = n * p + 0.5 ;
	const std::size_t j1 = index ( std::ceil ( h + 0.5 ) ) ; 
	const std::size_t j2 = index ( std::ceil ( h - 0.5 ) ) ; 
	const double x1 = data ( j1 ) ;
	const double x2 = data ( j2 ) ;
	return  0.5 * ( x1 + x2 ) ;
      }
    case Three :
      {
	const double h = n * p ;
	const std::size_t j = index ( h  ) ; 
	return data ( j ) ;
      }
    case Four :
      {
	const double h = n * p ;
	return QF ( h ) ;
      }
    case Five :
      {
	const double h = n * p + 0.5 ;
	return QF ( h ) ;
      }
    case Six :
      {
	const double h = ( n + 1 ) * p ;
	return QF ( h ) ;
      }
    case Seven :
      {
	const double h = ( n - 1 ) * p + 1 ;
	return QF ( h ) ;
      }
    case Eight :
      {
	const double h = ( n + 1./3  ) * p + 1./3  ;
	return QF ( h ) ;
      }
    case Nine :
      {
	const double h = ( n + 0.25 ) * p + 0.375 ;
	return QF ( h ) ;
      }
    default :
      return 0 ;
    } ;
  //
  return 0 ;
}
// ============================================================================
/* Get Harrel-Davis estimator for quantile function
 *  @param p  (INPUT) quantile 
 *  @retiurn Harrel-Davis quantile estimator
 *  @see https://doi.org/10.1093/biomet/69.3.635
 *  @see F.E. Harrel and C.E.Davis,  "A new distribution-free quantile estimator",
 *       Biometrika 63.9 (Dec 1982), pp. 635-640
 */
// ============================================================================
double Ostap::Math::ECDF::quantile_HD
( const double p ) const 
{
  if      ( !p     || s_zero  ( p     ) ) { return m_data.front () ; }
  else if ( 1 == p || s_equal ( p , 1 ) ) { return m_data.back  () ; } 
  //
  Ostap::Assert ( 0 <= p && p <= 1     , 
		  "Invalid quantile!"              , 
		  "Ostap::Math::ECDF::quantile_HD" , 
		  INVALID_QUANTILE , __FILE__  , __LINE__ ) ;
  //
  double result = 0 ;
  const std::size_t n = N () ;
  const double alpha = ( n + 1 ) * p ;
  const double beta  = ( n + 1 ) * ( 1 - p ) ;
  for ( std::size_t i = 0 ; i < n ; ++i )
    {
      const double ti    = ( i + 1.0 ) / n ;
      const double value = m_data [ i ] ; 
      if ( !value ) { continue ; } 
      // 
      if (  n < 100 )
	{
	  result += value * ::beta_inc ( alpha , beta , ti          , 20 * n ) ; 
	  result -= value * ::beta_inc ( alpha , beta , i * 1.0 / n , 20 * n ) ; 
	}
      else 
	{
	  const double t  =  ( i + 0.5 ) / n ;
	  result += value  * ::beta_log ( alpha, beta , t , 1.0 / n , 20 * n ) ; 
	}
    }
  //
  return result ;  
}
// ============================================================================
// statistics (as Counter)
// ============================================================================
Ostap::StatEntity Ostap::Math::ECDF::counter() const 
{
  Ostap::StatEntity cnt {} ;
  for ( auto v : m_data ) { cnt.add ( v ) ; } 
  return cnt ;
}  
// ============================================================================
/* statistics (as statistics)
 * @param stat (UPDATE) input statistic object
 * @return updated statistics object
 */
// =============================================================================
Ostap::Math::Statistic& 
Ostap::Math::ECDF::statistics
( Ostap::Math::Statistic& stat ) 
{
  for ( auto v : m_data ) { stat.update ( v ) ; } 
  return stat ;
}
// ============================================================================
// For weighted data 
// ============================================================================
/*  check that WECDF is OK: 
 *  - there are some entries
 *  - sum of weigths is positive 
 *  - sum of squaed weigths is positive 
 *  - remove elments <code>!std::isfinite</code>
 */
// ============================================================================
Ostap::Math::WECDF&
Ostap::Math::WECDF::check_me ()
{
  Data::iterator remove = std::remove_if
    ( m_data.begin () ,
      m_data.end   () ,
      [] ( const Data::value_type& item ) -> bool
      { return
          !std::isfinite ( item.first  ) ||
          !std::isfinite ( item.second ) || !item.second ; } ) ;
  // 
  m_data.erase ( remove , m_data.end() ) ; 
  //
  Ostap::Assert ( !m_data.empty()                     ,
                  "No data for Empirical CDF"         ,
                  "Ostap::Math::WCDF"                 ,
                  INVALID_DATA   , __FILE__ ,__LINE__ ) ;
  //
  m_sumw  = calc_sumw  () ;
  m_sumw2 = calc_sumw2 () ;
  Ostap::Assert ( 0 < m_sumw && 0 < m_sumw2      ,
                  "Non-positive sum of weights!" ,
                  "Ostap::Math::WECDF"           , 
                  INVALID_DATA   , __FILE__ ,__LINE__ ) ;
  //
  return *this ;
}
// ============================================================================
/* Constructor from  data
 *  data must be non-empty!
 */ 
// ============================================================================
Ostap::Math::WECDF::WECDF
( const Ostap::Math::WECDF::Data& data          ,
  const bool                      complementary )
  : m_data          ( data )
  , m_sumw          ( 0    ) 
  , m_sumw2         ( 0    ) 
  , m_complementary ( complementary ) 
{
  //
  if ( !std::is_sorted ( m_data.begin() , m_data.end () , COMPARE() ) )
    { std::sort ( m_data.begin() , m_data.end() , COMPARE () ) ; }
  //
  this->check_me() ;
}
// ============================================================================
/*  Constructor from data
 *  data must be non-empty!
 */ 
// ============================================================================
Ostap::Math::WECDF::WECDF
( const Ostap::Math::ECDF::Data&  data          ,
  const Ostap::Math::ECDF::Data&  weights       ,
  const bool                      complementary )
  : m_data          (   ) 
  , m_sumw          ( 0 ) 
  , m_sumw2         ( 0 ) 
  , m_complementary ( complementary )
{
  //
  const std::size_t nd = data   .size() ;
  const std::size_t nw = weights.size() ;
  //
  m_data.reserve ( nd ) ;
  for ( std::size_t i = 0 ; i < nd ; ++i )
    {
      const double value = data [ i ] ;
      const double weight = ( i < nw ) ? weights [ i ] : 1.0 ;
      if ( !std::isfinite ( value ) || !std::isfinite ( weight ) || !weight ) { continue ; } 
      m_data.emplace_back ( value , weight ) ;
    }
  //
  if ( !std::is_sorted ( m_data.begin() , m_data.end() , COMPARE () ) ) 
    { std::sort ( m_data.begin() , m_data.end() , COMPARE () ) ; } 
  //
  this->check_me() ;
}
// ============================================================================
// Standard constructor from  data
// ============================================================================
Ostap::Math::WECDF::WECDF
( const Ostap::Math::ECDF::Data&  data          ,
  const bool                      complementary  )
  :  WECDF ( data , Ostap::Math::ECDF::Data() , complementary )
{}
// ============================================================================
Ostap::Math::WECDF::WECDF
( const Ostap::Math::WECDF& right         ,
  const bool                complementary )
  : WECDF ( right )
{
  m_complementary = complementary ;
}
// ============================================================================
Ostap::Math::WECDF::WECDF
( const Ostap::Math::ECDF& right         ,
  const bool               complementary )
  : WECDF ( right )
{
  m_complementary = complementary ;
}
// ============================================================================
Ostap::Math::WECDF::WECDF
( const Ostap::Math::ECDF& right            )
  : m_data          () 
  , m_sumw          ( right.size ()         ) 
  , m_sumw2         ( right.size ()         ) 
  , m_complementary ( right.complementary() ) 
{
  m_data.reserve ( right.size() ) ;
  for ( auto d : right.data() ) { m_data.emplace_back ( d , 1.0 ) ; }
}
// ============================================================================
// swap two objects 
// ============================================================================
void Ostap::Math::WECDF::swap ( Ostap::Math::WECDF& right )
{
  std::swap ( m_data          , right.m_data          ) ;
  std::swap ( m_weights       , right.m_weights       ) ;
  std::swap ( m_sumw          , right.m_sumw          ) ;
  std::swap ( m_sumw2         , right.m_sumw2         ) ;
  std::swap ( m_complementary , right.m_complementary ) ;
}
// ============================================================================
// add a value to data container  
// ============================================================================
Ostap::Math::WECDF&
Ostap::Math::WECDF::add
( const Ostap::Math::WECDF::Entry& entry ) 
{
  const double value  = entry.first ;
  const double weight = entry.second ;
  //
  if ( !std::isfinite ( value ) || !std::isfinite ( weight ) || !weight) { return *this ; }  
  //
  auto where = std::upper_bound ( m_data.begin () , m_data.end () , entry , COMPARE () ) ;
  m_data.insert ( where , entry ) ;
  //
  m_sumw  += weight ;
  m_sumw2 += weight * weight ;
  //
  Ostap::Assert ( 0 < m_sumw && 0 <= m_sumw2         ,
                  "Non-positive sum of weights!"     ,
                  "Ostap::Math::WECDF::add"          ,
                  INVALID_DATA , __FILE__ , __LINE__ ) ;
  //
  return *this ;
}
// ============================================================================
// add values to data container  
// ============================================================================
Ostap::Math::WECDF&
Ostap::Math::WECDF::add
( const Ostap::Math::WECDF& values ) 
{
  Data tmp   ( m_data.size() + values.size() ) ;
  /// merge two sorted containers 
  std::merge ( m_data.begin        () ,
               m_data.end          () ,
               values.m_data.begin () ,
               values.m_data.end   () ,
               tmp.begin           () ,
	       COMPARE             () ) ;
  ///
  std::swap ( m_data , tmp  ) ;
  ///
  m_sumw  += values.m_sumw  ;
  m_sumw2 += values.m_sumw2 ;
  ///
  return *this ;
}
// ============================================================================
// add values to data container  
// ============================================================================
Ostap::Math::WECDF&
Ostap::Math::WECDF::add
( const Ostap::Math::ECDF& values ) 
{
  Data aux {} ; aux.reserve ( values.size() ) ; 
  for ( auto d : values.data() ) { aux.emplace_back ( d , 1.0 ) ; }
  //
  Data tmp   ( m_data.size() + aux.size() ) ;
  /// merge two sorted containers 
  std::merge ( m_data.begin () ,
               m_data.end   () ,
               aux.begin    () ,
               aux.end      () ,
               tmp.begin    () ,
	       COMPARE      () ) ;
  ///
  std::swap ( m_data , tmp  ) ;
  ///
  m_sumw  += values.size() ;
  m_sumw2 += values.size() ;
  ///
  return *this ;
}// ============================================================================
// add values to data container  
// ============================================================================
Ostap::Math::WECDF&
Ostap::Math::WECDF::add
( const Ostap::Math::ECDF::Data& values ) 
{
  Data aux {} ; aux.reserve ( values.size() ) ; 
  for ( auto d : values ) { aux.emplace_back ( d , 1.0 ) ; }
  //
  return this->add ( aux ) ;
}  
// ============================================================================
// add values to data container  
// ============================================================================
Ostap::Math::WECDF&
Ostap::Math::WECDF::add
( const Ostap::Math::WECDF::Data& values ) 
{
  const Data* input = &values ;
  //
  Data values2 ;  
  if ( !std::is_sorted ( values.begin() , values.end() , COMPARE () ) ) 
    {
      values2 = values ;
      std::sort ( values2.begin() , values2.end  () , COMPARE () ) ;
      input = &values2 ;
    }
  /// temporary  dataset 
  Data tmp   ( m_data.size() + input -> size() ) ;  
  /// merge two sorted containers 
  std::merge ( m_data.begin  () ,
               m_data.end    () ,
               input->begin  () ,
               input->end    () ,
               tmp.begin     () ,
	       COMPARE       () ) ;
  //
  std::swap ( m_data , tmp ) ;
  //
  return this->check_me() ;
}
// ============================================================================

// ============================================================================
// the main method 
// ============================================================================
double Ostap::Math::WECDF::evaluate   ( const double x ) const
{
  if      ( x < xmin () ) { return m_complementary ? 1.0 : 0.0 ; } 
  else if ( x > xmax () ) { return m_complementary ? 0.0 : 1.0 ; } 
  //
  const Entry entry { x , 1.0 } ;
  // NB: note the comparison criteria! 
  auto found = std::upper_bound ( m_data.begin () , m_data.end () , entry , COMPARE () ) ;
  //
  const double wsum   = calc_sumw ( found - m_data.begin() ) ;
  const double result = wsum / m_sumw ;
  //
  return m_complementary ? ( 1 - result ) : result ; 
}
// ============================================================================
// the main method 
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::WECDF::estimate ( const double x ) const
{
  const std::size_t NN = m_data.size() ;
  //
  typedef Ostap::Math::ValueWithError  VE ;
  //
  const VE all  { m_sumw , m_sumw2 } ;
  //
  if      ( x < xmin () )
    {
      const VE none { 0 , std::pow ( m_data.front().second , 2 ) } ;
      return m_complementary ?
        Ostap::Math::binomEff2 ( all  , none ) :
        Ostap::Math::binomEff2 ( none , all  ) ; }
  else if ( x > xmax () )
    {
      const VE none { 0 , std::pow ( m_data.back().second , 2 ) } ;
      return m_complementary ?
        Ostap::Math::binomEff2 ( none , all  ) :
        Ostap::Math::binomEff2 ( all  , none ) ;
    }
  //
  const Entry entry { x , 1.0 } ;
  // NB: note the comparison criteria! 
  auto found = std::upper_bound ( m_data.begin () , m_data.end () , entry , COMPARE () ) ;
  //
  const double wsum   = calc_sumw  ( found - m_data.begin() ) ;
  const double w2sum  = calc_sumw2 ( found - m_data.begin() ) ;
  //  
  const VE  acc {          wsum ,           w2sum } ;
  const VE  rej { m_sumw - wsum , m_sumw2 - w2sum } ;
  //
  return m_complementary ?
    Ostap::Math::binomEff2 ( acc , rej ) : 
    Ostap::Math::binomEff2 ( rej , acc ) ;
}
// ============================================================================
// get ranks of the elements in the (pooled) sample
// ============================================================================
Ostap::Math::WECDF::Indices
Ostap::Math::WECDF::ranks
( const Ostap::Math::ECDF& sample ) const
{
  const Data::size_type N  = size() ;
  // fill outptut array with N
  Indices result ( sample.size () ,  N ) ;
  Data::size_type NS = sample.size() ;
  for ( Data::size_type i = 0 ; i < NS ; ++i )
    {
      Data::size_type r = rank ( sample.data ( i ) ) ;
      result [ i ] = r ;
      // try to be a bit more efficient, the rest of array is filled with N 
      if ( N <= r ) { break ; }
    }
  return result ;
}
// ============================================================================
// get ranks of the elements in the (pooled) sample
// ============================================================================
Ostap::Math::WECDF::Indices
Ostap::Math::WECDF::ranks
( const Ostap::Math::WECDF& sample ) const
{
  const Data::size_type N  = size() ;
  // fill output array with N
  Indices result ( sample.size () ,  N ) ;
  Data::size_type NS = sample.size() ;
  for ( Data::size_type i = 0 ; i < NS ; ++i )
    {
      Data::size_type r = rank ( sample.data ( i ) ) ;
      result [ i ] = r ;
      // try to be a bit more efficient, the rest of array is filled with N 
      if ( N <= r ) { break ; }
    }
  return result ;
}
// ============================================================================  
// statistics (as Counter)
// ============================================================================
Ostap::WStatEntity 
Ostap::Math::WECDF::counter() const 
{
  Ostap::WStatEntity cnt{} ;
  for ( const auto& v : m_data ) { cnt.add( v.first , v.second ) ; } 
  return cnt ;
}  
// ============================================================================
// copy assignement 
// ============================================================================
Ostap::Math::WECDF&
Ostap::Math::WECDF::operator=
( const Ostap::Math::WECDF&  right )
{
  if ( this == &right ) { return *this ; }
  m_data          = right.m_data          ;
  m_weights       = right.m_weights       ;
  m_sumw          = right.m_sumw          ;
  m_sumw2         = right.m_sumw2         ;
  m_complementary = right.m_complementary ;
  return *this ;
}
// ============================================================================
// move assignement 
// ============================================================================
Ostap::Math::WECDF&
Ostap::Math::WECDF::operator=
( Ostap::Math::WECDF&&  right )
{
  if ( this == &right ) { return *this ; }
  this->swap ( right ) ;
  return *this ; 
}


// ============================================================================
/* statistics (as statistics)
 * @param stat (UPDATE) input statistic object
 * @return updated statistics object
 */
// =============================================================================
Ostap::Math::WStatistic& 
Ostap::Math::WECDF::statistics
( Ostap::Math::WStatistic& stat ) 
{
  for ( const auto& v : m_data ) { stat.update ( v.first , v.second  ) ; } 
  return stat ;
}
// ============================================================================
/* Get Harrel-Davis estimator for quantile function
 *  @param p  (INPUT) quantile 
 *  @retiurn Harrel-Davis quantile estimator
 *
 *  @see https://arxiv.org/abs/2304.07265
 *  @see Andrey Akinshin, "Weighted quantile estimators", arXiv:2304.07265
 *   
 *  @see https://doi.org/10.1093/biomet/69.3.635
 *  @see F.E. Harrel and C.E.Davis,  "A new distribution-free quantile estimator",
 *       Biometrika 63.9 (Dec 1982), pp. 635-640
 */
// ============================================================================
double Ostap::Math::WECDF::quantile_HD ( const double p ) const 
{
  if      ( !p     || s_zero  ( p     ) ) { return xmin () ; }
  else if ( 1 == p || s_equal ( p , 1 ) ) { return xmax () ; } 
  //
  Ostap::Assert ( 0 <= p && p <= 1     , 
      "Invalid quantile!"              , 
      "Ostap::Math::WECDF::quantile_HD" , 
      INVALID_QUANTILE , __FILE__  , __LINE__ ) ;
//
const double sw_inv { 1./sumw () } ;
const double nstar  { nEff () } ; 
const double alpha = ( nstar + 1 ) * p ;
const double beta  = ( nstar + 1 ) * ( 1 - p ) ;
//
double result = 0 ;
const std::size_t n = N () ;
double wsum         = 0    ;  
for ( std::size_t i = 0 ; i < n ; ++i )
  {
    const double value  = m_data [ i ].first  ;
    const double weight = m_data [ i ].second ;
    if ( !value || !weight ) { continue ; }
    //
    const double tp = sw_inv * ( wsum          ) ;
    const double ti = sw_inv * ( wsum + weight ) ;
    wsum += weight ;
    // 
    result += value * ::beta_inc ( alpha , beta , ti , 20 * N () ) ; 
    result -= value * ::beta_inc ( alpha , beta , tp , 20 * N () ) ; 
  }
  //
  return result ;  
}
// ============================================================================
/* create the empirical PDF from empirical CDF 
 *  @attention data are not copied!
 */
// ============================================================================
Ostap::Math::EPDF::EPDF
( const Ostap::Math::ECDF&                    cdf ,
  const Ostap::Math::DensityEstimator::Kernel k   ,
  const double                                h   ) 
  : m_cdf   ( cdf )
  , m_k     ( k   )
  , m_h     ( h   )
{
  Ostap::Assert ( Ostap::Math::DensityEstimator::Uniform <= k &&
		  Ostap::Math::DensityEstimator::Last    >= k && k < s_TABLE.size () , 
		  "Invalid Kernel!"   , 
		  "Ostap::Math::EPDF" ,
		  INVALID_KERNEL , __FILE__, __LINE__ ) ;
  // check the smoothing parameter 
  if ( m_h <= 0 ) { m_h = Ostap::Math::DensityEstimator::hopt ( m_cdf ) ; }
  Ostap::Assert ( 0 < m_h             ,
		  "Invalid smoothing parameter" ,
		  "Ostap::Math::EPDF" ,
		  INVALID_SMOOTH , __FILE__ , __LINE__ ) ;  
}
// =============================================================================
// get the PDF
// =============================================================================
double Ostap::Math::EPDF::evaluate ( const double x ) const
{
  FUNPTR k = std::get<1> ( s_TABLE[m_k] ) ;
  double s = std::get<2> ( s_TABLE[m_k] ) ;
  if  ( s <= 0 ) { s = 5 ; } ;
  //
  const double xmn = x - s * m_h ;  
  const double xmx = x + s * m_h ;
  //
  if ( m_cdf.xmax () < xmn || m_cdf.xmin () > xmx ) { return 0 ; }
  Ostap::Math::ECDF::iterator imin = std::lower_bound ( m_cdf.begin () , m_cdf.end () , xmn ) ;
  if ( imin == m_cdf.end() ) { return 0 ; }
  Ostap::Math::ECDF::iterator imax = std::upper_bound ( imin + 1 , m_cdf.end () , xmx ) ;
  double       value = 0 ;
  const double ih    = 1 / m_h ;
  for ( Ostap::Math::ECDF::iterator i = imin ; i < imax ; ++i )
    {
      const double u = ih * ( x - *i ) ;
      value += (*k) ( u ) ;
    }
  return ih * value / m_cdf.N() ;
}
// ============================================================================
/* create the emppirical PDF from empirical CDF 
 *  @attention data are not copied!
 */
// ============================================================================
Ostap::Math::WEPDF::WEPDF
( const Ostap::Math::WECDF&                   cdf ,
  const Ostap::Math::DensityEstimator::Kernel k   ,
  const double                                h   ) 
  : m_cdf ( cdf )
  , m_k   ( k   )
  , m_h   ( h   )
{
  Ostap::Assert ( Ostap::Math::DensityEstimator::Uniform <= k &&
		  Ostap::Math::DensityEstimator::Last    >= k && k < s_TABLE.size () , 
		  "Invalid Kernel!"    , 
		  "Ostap::Math::WEPDF" ,
		  INVALID_KERNEL , __FILE__, __LINE__ ) ;
  // check the smoothing parameter 
  if ( m_h <= 0 ) { m_h = Ostap::Math::DensityEstimator::hopt ( m_cdf ) ; }
  Ostap::Assert ( 0 < m_h             ,
		  "Invalid smoothing parameter" ,
		  "Ostap::Math::WEPDF" ,
		  INVALID_SMOOTH , __FILE__ , __LINE__ ) ;  
}
// =============================================================================
// get the PDF
// =============================================================================
double Ostap::Math::WEPDF::evaluate ( const double x ) const
{
  FUNPTR k = std::get<1> ( s_TABLE[m_k] ) ;
  double s = std::get<2> ( s_TABLE[m_k] ) ;
  if  ( s <= 0 ) { s = 5 ; } ;
  //
  const double xmn = x - s * m_h ;  
  const double xmx = x + s * m_h ;
  //
  if ( m_cdf.xmax () < xmn || m_cdf.xmin () > xmx ) { return 0 ; }
  //
  const Ostap::Math::WECDF::COMPARE cmp {} ;
  //
  Ostap::Math::WECDF::iterator imin = std::lower_bound ( m_cdf.begin () , m_cdf.end () , xmn , cmp ) ;
  if ( imin == m_cdf.end() ) { return 0 ; }
  Ostap::Math::WECDF::iterator imax = std::upper_bound ( imin + 1 , m_cdf.end () , xmx , cmp ) ;
  double       value = 0 ;
  const double ih    = 1 / m_h ;
  for ( Ostap::Math::WECDF::iterator i = imin ; i < imax ; ++i )
    {
      const double w = i->second ;
      if ( !w ) { continue ; } 
      const double u = ih * ( x - i->first ) ;
      value += w * (*k) ( u ) ;
    }
  return ih * value / m_cdf.sumw ();
}
// ============================================================================
// update smoothing parameters
// ============================================================================
bool Ostap::Math::EPDF::setH      ( const double h )
{
  if ( ( 0 < h ) && s_equal ( m_h , h ) ) { return false ; }
  m_h = h ;
  // check the smoothing parameter 
  if ( m_h <= 0 ) { m_h = Ostap::Math::DensityEstimator::hopt ( m_cdf ) ; }
  return true ;
}
// ============================================================================
// update kernel kernel
// ============================================================================
bool Ostap::Math::EPDF::setKernel
( const Ostap::Math::DensityEstimator::Kernel k )
{
  Ostap::Assert ( Ostap::Math::DensityEstimator::Uniform <= k &&
		  Ostap::Math::DensityEstimator::Last    >= k && k < s_TABLE.size () , 
		  "Invalid Kernel!"    , 
		  "Ostap::Math::EPDF::setKernel" ,
		  INVALID_KERNEL , __FILE__, __LINE__ ) ;
  if ( k == m_k ) { return false ; }
  m_k = k ;
  return true ;
}
// ============================================================================
// update smoothing parameters
// ============================================================================
bool Ostap::Math::WEPDF::setH      ( const double h )
{
  if ( ( 0 < h ) && s_equal ( m_h , h ) ) { return false ; }
  m_h = h ;
  // check the smoothing parameter 
  if ( m_h <= 0 ) { m_h = Ostap::Math::DensityEstimator::hopt ( m_cdf ) ; }
  return true ;
}
// ============================================================================
// update kernel kernel
// ============================================================================
bool Ostap::Math::WEPDF::setKernel
( const Ostap::Math::DensityEstimator::Kernel k )
{
  Ostap::Assert ( Ostap::Math::DensityEstimator::Uniform <= k &&
		  Ostap::Math::DensityEstimator::Last    >= k && k < s_TABLE.size () , 
		  "Invalid Kernel!"    , 
		  "Ostap::Math::WEPDF::setKernel" ,
		  INVALID_KERNEL , __FILE__, __LINE__ ) ;
  if ( k == m_k ) { return false ; }
  m_k = k ;
  return true ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================

