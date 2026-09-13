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
#include "Ostap/Quantiles.h"
#include "Ostap/StatusCode.h"
// ============================================================================
// ROOT 
// ============================================================================
#include "TH1.h"
// ============================================================================
// Local
// ============================================================================
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
  /// cached version of normalized incomplete beta-function
  double beta_inc  
  ( const double alpha      ,  
    const double beta       , 
    const double z          ,
    std::size_t  N  = 50000 ) 
  {
   /// cache for normalized incomplete beta-function;
   typedef std::map<std::size_t,double>  MAP   ;
   typedef SyncedCache<MAP>              CACHE ;
   /// the actual cache
   static CACHE   s_cache {}      ; // integration cache
   ///
   const std::size_t key { Ostap::Utils::hash_combiner ( z , alpha , beta , z ) }; 
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
    std::size_t  N  = 50000 ) 
  {
   /// cache for normalized incomplete beta-function;
   typedef std::map<std::size_t,double>  MAP   ;
   typedef SyncedCache<MAP>              CACHE ;
   /// the actual cache
   static CACHE   s_cache {}      ; // integration cache
   ///
   const std::size_t key { Ostap::Utils::hash_combiner ( alpha , beta , t , dt ) }; 
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
  // ==========================================================================
}
// ============================================================================

// ============================================================================
// Standard constructor from  data
// ============================================================================
Ostap::Math::ECDF::ECDF
( const Ostap::Math::ECDF::Data&  data          ,
  const bool                      complementary )
  : Ostap::Math::ECDF::ECDF ( data.begin() , data.end() , complementary )
{}
// ============================================================================
// Standard constructor from  data
// ============================================================================
Ostap::Math::ECDF::ECDF
( const double  data          ,
  const bool    complementary )
  : Ostap::Math::ECDF::ECDF ( std::vector<double>( 1 , data ) , complementary )
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
// check that ECDF is OK: there are some entries 
// ============================================================================
Ostap::Math::ECDF&
Ostap::Math::ECDF::cleanup ()
{
  // find invalid elements 
  Data::iterator remove = std::remove_if
    ( m_data.begin () ,
      m_data.end   () ,
      [] ( const double x ) -> bool
      { return !std::isfinite ( x ) ; } ) ;
  // remove them! 
  m_data.erase ( remove , m_data.end() ) ; 
  //
  return *this ;
}
// ============================================================================
// merge two objects 
// ============================================================================
Ostap::Math::ECDF
Ostap::Math::ECDF::merge
( const Ostap::Math::ECDF& right ) const 
{ ECDF result { *this } ; result.add ( right ) ; return result  ; }
// ============================================================================
// swap two objects 
// ============================================================================
void Ostap::Math::ECDF::swap
( Ostap::Math::ECDF& right )
{
  std::swap ( m_data          , right.m_data          ) ;
  std::swap ( m_complementary , right.m_complementary ) ;
  m_counter.swap ( right.m_counter ) ;
}
// ============================================================================
// the main method 
// ============================================================================
double Ostap::Math::ECDF::evaluate   ( const double x ) const
{
  //
  if      (     m_data.empty () ) { return m_complementary ? 1.0 : 0.0 ; }              
  else if ( x < m_data.front () ) { return m_complementary ? 1.0 : 0.0 ; } 
  else if ( x > m_data.back  () ) { return m_complementary ? 0.0 : 1.0 ; } 
  //
  const double result = double ( rank ( x ) ) / m_data.size () ;
  return m_complementary ? ( 1 - result ) : result ; 
}
// ============================================================================
// the main method 
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ECDF::estimate
( const double x ) const
{
  static const Ostap::Math::ValueWithError s_VE_zero { 0.0 , 0.0 } ;
  static const Ostap::Math::ValueWithError s_VE_one  { 1.0 , 0.0 } ;
  //
  const std::size_t NN = m_data.size() ;
  //
  if      (     m_data.empty () ) { return m_complementary ? s_VE_one : s_VE_zero ; } 
  else if ( x < m_data.front () ) { return Ostap::Math::binomEff ( m_complementary ? NN : 0u , NN ) ; }
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
  m_data   .insert ( where , value ) ;
  m_counter.add    ( value ) ;
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
  Data tmp   ( values.size() + m_data.size () ) ;
  /// merge two sorted containers 
  std::merge ( values.m_data.begin () ,                     
               values.m_data.end   () ,
               m_data.begin        () ,
               m_data.end          () ,                     
               tmp.begin           () ) ;
  /// swap the merged result  with own data 
  std::swap ( m_data , tmp ) ;
  //
  m_counter += values.m_counter ;
  // 
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
/*  get p-quantile uaing Hyndman-Fan estimator 
 *  @see https://arxiv.org/abs/2304.07265
 *  @see Andrey Akinshin, "Weighted quantile estimators", arXiv:2304.07265     
 */
// =============================================================================
double Ostap::Math::ECDF::quantile
( const double                               p ,
  const Ostap::QuantileTypes::HyndmanFanType t ) const 
{
  //
  if      ( !p     || s_zero  ( p     ) || p < 0 ) { return low_edge  () ; }
  else if ( 1 == p || s_equal ( p , 1 ) || p > 1 ) { return high_edge () ; }
  //
  Ostap::Assert ( 0 < p && p < 1     , 
                  "Invalid probability"                      , 
                  "Ostap::Math::ECDF::quantile"              , 
                  INVALID_PROBABILITY , __FILE__  , __LINE__ ) ;
  //
  if ( 1 == size() ) { return m_data.front() ; }
  //
  const Ostap::Math::HyndmanFan estimator { t , false } ;
  return estimator.quantile ( m_data.begin () , m_data.end   () , p ) ;
}
// =============================================================================
/*  Get Harrell-Davis estimator for quantile function
 *  @param p  (INPUT) quantile 
 *  @retiurn Harrel-Davis quantile estimator
 *  @see https://doi.org/10.1093/biomet/69.3.635
 *  @see F.E. Harrel and C.E.Davis,  "A new distribution-free quantile estimator",
 *       Biometrika 63.9 (Dec 1982), pp. 635-640
 */
// =============================================================================
double Ostap::Math::ECDF::quantile
( const double                                     p ,
  const Ostap::QuantileTypes::HarrellDavisType& /* t */ ) const
{
  //
  if      ( !p     || s_zero  ( p     ) || p < 0 ) { return low_edge  () ; }
  else if ( 1 == p || s_equal ( p , 1 ) || p > 1 ) { return high_edge () ; }
  //
  Ostap::Assert ( 0 < p && p < 1     , 
                  "Invalid probability"                      , 
                  "Ostap::Math::ECDF::quantile"              , 
                  INVALID_PROBABILITY , __FILE__  , __LINE__ ) ;
  //
  if ( 1 == size() ) { return m_data.front() ; }
  //
  
  const Ostap::Math::HarrellDavis estimator { false } ;
  return estimator.quantile ( m_data.begin () , m_data.end   () , p ) ;
}
// =============================================================================
/*  get p-quantile of distribution: \f$ 1 \le p \le1  \f$
 *  @see scipy.stats.mstats
 *
 *  Typical values for alphap, betap are:
 * 
 * - (0,1) : p(k) = k/n : linear interpolation of cdf (R type 4)
 * - (.5,.5) : p(k) = (k - 1/2.)/n : piecewise linear function (R type 5)
 * - (0,0) : p(k) = k/(n+1) : (R type 6)
 * - (1,1) : p(k) = (k-1)/(n-1): p(k) = mode[F(x[k])]. (R type 7, R default)
 * - (1/3,1/3): p(k) = (k-1/3)/(n+1/3): Then p(k) ~ median[F(x[k])].
 *   The resulting quantile estimates are approximately median-unbiased 
 *   regardless of the distribution of x. (R type 8)
 * - (3/8,3/8): p(k) = (k-3/8)/(n+1/4): Blom.
 *   The resulting quantile estimates are approximately 
 *   unbiased if x is normally distributed (R type 9)
 * - 0(.4,.4) : approximately quantile unbiased (Cunnane)
 * - (.35,.35): APL, used with PWM
 *
 *  @param p      (INPUT) quantile
 *  @param alphap (INPUT) parameter alphap \f$ 0 \le\alpha_p \le 1 \f$ 
 *  @param abetap (INPUT) parameter betap \f$ 0 \le\beta_p \le 1 \f$ 
 */    
// =============================================================================
double Ostap::Math::ECDF::quantile
( const double                                p  ,
  const Ostap::QuantileTypes::ABQuantileType& ab ) const
{
  //
  if      ( !p     || s_zero  ( p     ) || p < 0 ) { return low_edge  () ; }
  else if ( 1 == p || s_equal ( p , 1 ) || p > 1 ) { return high_edge () ; }
  //
  Ostap::Assert ( 0 < p && p < 1     , 
                  "Invalid probability"                      , 
                  "Ostap::Math::ECDF::quantile"              , 
                  INVALID_PROBABILITY , __FILE__  , __LINE__ ) ;
  //
  if ( 1 == size() ) { return m_data.front() ; }
  //
  const Ostap::Math::ABQuantile estimator { ab , false } ;
  return estimator.quantile ( m_data.begin () , m_data.end   () , p ) ;  
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
// project data content into 1D histogram
// ============================================================================
void Ostap::Math::ECDF::project 
( TH1& histo ) const 
{
  Ostap::Assert ( 1 == histo.GetDimension()      , 
                  "Invalild histogram dimension" , 
                  "Ostap::Math::ECDF::project"   , 
                  INVALID_TH1                    ) ; 
  histo.Reset() ;
  if ( !histo.GetSumw2() ) { histo.Sumw2() ; }
  histo.FillN ( m_data.size () , m_data.data () , nullptr ) ;
}
// ============================================================================
// low_edge <= xmin 
// ============================================================================
double Ostap::Math::ECDF::low_edge  () const
{ return m_data.empty() ? -std::numeric_limits<double>::max () : xmin () ; } 
// ============================================================================
///  xmax < high_edge 
// ============================================================================
double Ostap::Math::ECDF::high_edge () const
{
  return m_data.empty() ?
    +std::numeric_limits<double>::max () :
    Ostap::Math::high_edge ( xmax () ) ;
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
Ostap::Math::WECDF::cleanup ()
{
  Data::iterator remove = std::remove_if
  ( m_data.begin () ,
    m_data.end   () ,
    [] ( const Data::value_type& item ) -> bool
    { return
        !std::isfinite ( item.first  ) ||
        !std::isfinite ( item.second ) || !item.second ; } ) ;
  m_data.erase ( remove , m_data.end() ) ; 
  return *this ;
}
// ============================================================================
// constructor 
// ============================================================================
Ostap::Math::WECDF::WECDF
( const bool                      complementary )
  : m_data          ( )
  , m_complementary ( complementary ) 
  , m_counter       ( ) 
{} 
// ============================================================================
/* Constructor from  data
 *  data must be non-empty!
 */ 
// ============================================================================
Ostap::Math::WECDF::WECDF
( const Ostap::Math::WECDF::Data& data          ,
  const bool                      complementary )
  : m_data          ( data )
  , m_complementary ( complementary ) 
  , m_counter       ( ) 
{
  // (1) remove bad elements 
  cleanup () ;
  // (2) sort it
  std::sort ( m_data.begin() , m_data.end() , COMPARE () ) ;  
  // (3) update counter
  for ( const auto& v : m_data ) { m_counter.add ( v.first , v.second ) ; } 
}
// ============================================================================
/*  Constructor from data
 *  data must be non-empty!
 */ 
// ============================================================================
Ostap::Math::WECDF::WECDF
( const Ostap::Math::ECDF::Data&  data          ,
  const Ostap::Math::ECDF::Data&  weights       ,
  const double                    weight        , 
  const bool                      complementary )
  : WECDF ( complementary ) 
{
  //
  add ( data , weights , weight ) ;
}
// ============================================================================
// Standard constructor from  data
// ============================================================================
Ostap::Math::WECDF::WECDF
( const Ostap::Math::ECDF::Data&  data          ,
  const double                    weight        , 
  const bool                      complementary )
  : WECDF ( complementary ) 
{
  add ( data , weight ) ;
}
// ============================================================================
Ostap::Math::WECDF::WECDF
( const Ostap::Math::ECDF& right         ,
  const double             weight        , 
  const bool               complementary )
  : WECDF ( complementary ) 
{
  add ( right , weight ) ;
}
// ============================================================================
Ostap::Math::WECDF::WECDF
( const Ostap::Math::WECDF& right         ,
  const double              weight        , 
  const bool                complementary )
  : WECDF ( complementary ) 
{
  ///
  if ( std::isfinite ( weight ) && weight )
  {
    /// copy data 
    m_data = right.m_data ;
    /// scale data 
    if ( 1 == weight || s_equal ( 1 , weight ) ) { m_counter = right.m_counter ; }
    else
    {
      for ( auto& entry : m_data )
      {
        entry.second *= weight ;
        m_counter.add ( entry.first , entry.second ) ;
      }
    }
  }
}
// ============================================================================
// merge two objects 
// ============================================================================
Ostap::Math::WECDF
Ostap::Math::WECDF::merge
( const Ostap::Math::WECDF& right ) const 
{ WECDF result { *this } ; result.add ( right ) ; return result  ; }
// ============================================================================
// merge two objects 
// ============================================================================
Ostap::Math::WECDF
Ostap::Math::WECDF::merge
( const Ostap::Math::ECDF&  right ) const 
{ WECDF result { *this } ; result.add ( right ) ; return result  ; }
// ============================================================================
// swap two objects 
// ============================================================================
void Ostap::Math::WECDF::swap
( Ostap::Math::WECDF& right )
{
  std::swap ( m_data          , right.m_data          ) ;
  std::swap ( m_complementary , right.m_complementary ) ;
  m_counter.swap ( right.m_counter ) ;
}
// ============================================================================
/** scale all weights by a factor s
 *  - no action is s if not finite or one 
 *  - reset if s is zero
 */
// ============================================================================
Ostap::Math::WECDF&
Ostap::Math::WECDF::scale
( const double s )
{
  /// no action 
  if ( !std::isfinite ( s ) || 1 == s || s_equal ( s , 1 ) ) { return *this ; }
  //
  /// clear/reset/erase 
  if ( !s ) { m_data.clear () ; m_counter.reset() ; return *this ; }
  //
  /// scale
  m_counter.reset () ;
  for ( auto& entry : m_data )
  {
    entry.second *= s ;
    m_counter.add ( entry.first , entry.second ) ;
  }
  //
  return *this ;
}
// ============================================================================

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
  if ( !std::isfinite ( value ) || !std::isfinite ( weight ) || !weight ) { return *this ; }  
  //
  auto where = std::upper_bound ( m_data.begin () , m_data.end () , entry , COMPARE () ) ;
  m_data.insert ( where , entry  ) ;
  // update the counter 
  m_counter.add ( value , weight ) ;
  return *this ;
}
// ============================================================================
// add values to data container  
// ============================================================================
Ostap::Math::WECDF&
Ostap::Math::WECDF::add
( const Ostap::Math::WECDF& values , 
  const double              weight )    
{
  if ( !std::isfinite ( weight ) || !weight ) { return *this ; }
  
  // allocate temporary storage 
  Data tmp ( m_data.size() + values.size() ) ;
  //
  /// apply additional weight?
  if ( 1 != weight && !s_equal ( 1 , weight ) )
  { 
    Data scaled_data { values.m_data } ;
    for ( auto& v : scaled_data ) { v.second *= weight ; }
    //
    std::merge ( m_data.begin      () ,
                 m_data.end        () ,
                 scaled_data.begin () ,
                 scaled_data.end   () ,
                 tmp.begin         () ,
                 COMPARE           () ) ;
    // update counter 
    for ( const auto& v : scaled_data) { m_counter.add ( v.first, v.second ) ; }
    //
    std::swap ( m_data , tmp  ) ;
    return *this ; 
  }
  //
  /// simple merge 
  std::merge ( m_data.begin        () ,
               m_data.end          () ,
               values.m_data.begin () ,
               values.m_data.end   () ,
               tmp.begin           () ,
               COMPARE             () ) ;
  // update counter 
  m_counter += values.m_counter ; 
  ///
  std::swap ( m_data , tmp  ) ;
  return *this ;
}

// ============================================================================
// add values to data container  
// ============================================================================
Ostap::Math::WECDF&
Ostap::Math::WECDF::add
( const Ostap::Math::ECDF& values , 
  const double             weight )  
{
  if ( !std::isfinite ( weight ) || !weight ) { return *this ; } 
  //  
  Data aux {} ; aux.reserve ( values.size() ) ; 
  for ( auto d : values.data() ) { aux.emplace_back ( d , weight ) ; }
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
  /// update counter 
  for ( const auto& v : aux ) { m_counter.add ( v.first , v.second ) ; }
  return *this ;
}
// ============================================================================
// add values to data container  
// ============================================================================
Ostap::Math::WECDF&
Ostap::Math::WECDF::add
( const Ostap::Math::ECDF::Data& values ,
  const double                   weight ) 
{
  if ( !std::isfinite ( weight ) || !weight ) { return *this ; }
  //
  /// (1) create temporary array 
  Data values2 {} ; values2.reserve ( values.size() ) ;
  /// (2) fill it with on-fly removal of bad values  
  for ( auto d : values ) { if ( std::isfinite ( d ) ) { values2.emplace_back ( d , weight ) ; } }
  /// (3) sort it 
  std::sort  ( values2.begin() , values2.end  () , COMPARE () ) ;        
  /// (4) one more temporary  dataset 
  Data tmp   ( m_data.size() + values2.size() ) ;  
  /// (5) merge two sorted containers 
  std::merge ( m_data.begin  () ,
               m_data.end    () ,
               values2.begin () ,
               values2.end   () ,
               tmp.begin     () ,
               COMPARE       () ) ;
  /// (6) swap containers 
  std::swap ( m_data , tmp ) ;
  /// (7) update counters 
  for ( const auto& v : values2 ) { m_counter.add ( v.first , v.second ) ; }
  return *this ;
}  
// ============================================================================
// add values to data container  
// ============================================================================
Ostap::Math::WECDF&
Ostap::Math::WECDF::add
( const Ostap::Math::WECDF::Data& values ,
  const double                    weight ) 
{
  // ATTENTION!
  if ( !std::isfinite ( weight ) || !weight ) { return *this ; } 
  //
  /// (1) make local copy 
  Data values2 { values } ;
  /// (2) scale it if needed
  if ( 1 != weight && !s_equal ( 1 , weight ) )
  { for ( auto& entry : values2 ) { entry.second *= weight; } }
  ///
  /// (3) remove bad elements 
  Data::iterator remove = std::remove_if
  ( values2.begin () ,
    values2.end   () ,
    [] ( const Data::value_type& item ) -> bool
    { return
        !std::isfinite ( item.first  ) ||
        !std::isfinite ( item.second ) || !item.second ; } ) ;
  if ( values2.end() != remove ) { values2.erase ( remove , values2.end() ) ; } 
  /// (4) sort it 
  std::sort ( values2.begin() , values2.end  () , COMPARE () ) ;      
  ///
  /// (5) temporary  dataset 
  Data tmp   ( m_data.size() + values2.size() ) ;  
  /// (6) merge two sorted containers 
  std::merge ( m_data.begin  () ,
               m_data.end    () ,
               values2.begin () ,
               values2.end   () ,
               tmp.begin     () ,
               COMPARE       () ) ;
  /// (7) swap containers 
  std::swap ( m_data , tmp ) ;
  /// (8) update counters 
  for ( const auto& v : values2 ) { m_counter.add ( v.first , v.second ) ; }
  return *this ;
}
// ============================================================================
// add values to data container  
// ============================================================================
Ostap::Math::WECDF&
Ostap::Math::WECDF::add
( const Ostap::Math::ECDF::Data& values  ,
  const Ostap::Math::ECDF::Data& weights , 
  const double                   weight  ) 
{
  //
  // ATTENTION!
  if ( !std::isfinite ( weight ) || !weight ) { return *this ; } 
  //
  Ostap::Assert ( values.size() == weights.size() , 
                  "Mismatch values/weights size"  , 
                  "Ostap::Math::WECDF::add"       ,
                  INVALID_VALUES_WEIGHTS_MATCH    , __FILE__ , __LINE__ ) ;
  //
  /// (1)  merge/zip two containers 
  Data values2 ( values.size () ) ;
  std::transform ( values .begin () ,
                   values .end   () ,
                   weights.begin () ,
                   values2.begin () ,
                   [weight] ( const double v , const double w )
                   { return std::make_pair ( v , w * weight ) ; } ) ;

  /// (2) remove bad elements 
  Data::iterator remove = std::remove_if
  ( values2.begin () ,
    values2.end   () ,
    [] ( const Data::value_type& item ) -> bool
    { return
        !std::isfinite ( item.first  ) ||
        !std::isfinite ( item.second ) || !item.second ; } ) ;
  if ( values2.end() != remove ) { values2.erase ( remove , values2.end() ) ; } 
  ///
  /// (4) sort it 
  std::sort ( values2.begin() , values2.end  () , COMPARE () ) ;      
  /// (5) temporary  dataset 
  Data tmp   ( m_data.size() + values2.size() ) ;  
  /// (6) merge two sorted containers 
  std::merge ( m_data.begin  () ,
               m_data.end    () ,
               values2.begin () ,
               values2.end   () ,
               tmp.begin     () ,
               COMPARE       () ) ;
  /// (7) swap containers 
  std::swap ( m_data , tmp ) ;
  /// (8) update counters 
  for ( const auto& v : values2 ) { m_counter.add ( v.first , v.second ) ; }
  return *this ;
}
// ============================================================================
// the main method 
// ============================================================================
double Ostap::Math::WECDF::evaluate   ( const double x ) const
{
  //
  if      ( m_data.empty () ) { return m_complementary ? 1.0 : 0.0 ; }              
  else if ( x < xmin ()     ) { return m_complementary ? 1.0 : 0.0 ; } 
  else if ( x > xmax ()     ) { return m_complementary ? 0.0 : 1.0 ; } 
  //
  const Entry entry { x , 1.0 } ;
  // NB: note the comparison criteria! 
  auto found = std::upper_bound ( m_data.begin () , m_data.end () , entry , COMPARE () ) ;
  //
  const double wsum   = ( found == m_data.end() ) ? sumw() : sumw ( found - m_data.begin() ) ;
  const double result = wsum / sumw ()  ;
  //
  return m_complementary ? ( 1 - result ) : result ; 
}
// ============================================================================
// the main method 
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::WECDF::estimate ( const double x ) const
{
  //
  static const Ostap::Math::ValueWithError s_VE_zero { 0.0 , 0.0 } ;
  static const Ostap::Math::ValueWithError s_VE_one  { 1.0 , 0.0 } ;
  //
  if  ( m_data.empty () ) { return m_complementary ? s_VE_one : s_VE_zero ; } 
  //
  const std::size_t NN = m_data.size() ;
  //
  typedef Ostap::Math::ValueWithError  VE ;
  //
  const VE all  { sumw () , sumw2 () } ;
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
  const double wsum  = ( found == m_data.end() ) ? sumw  () : sumw  ( found - m_data.begin() ) ;
  const double w2sum = ( found == m_data.end() ) ? sumw2 () : sumw2 ( found - m_data.begin() ) ;
  //  
  const VE acc {           wsum ,           w2sum } ;
  const VE rej { sumw () - wsum , sumw2() - w2sum } ;
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
double Ostap::Math::WECDF::quantile
( const double                                     p    ,
  const Ostap::QuantileTypes::HarrellDavisType& /* t */ ) const 
{
  //
  if      ( !p     || s_zero  ( p     ) || p < 0 ) { return low_edge  () ; }
  else if ( 1 == p || s_equal ( p , 1 ) || p > 1 ) { return high_edge () ; }
  //
  Ostap::Assert ( 0 < p && p < 1     , 
                  "Invalid probability"                      , 
                  "Ostap::Math::WECDF::quantile"             , 
                  INVALID_PROBABILITY , __FILE__  , __LINE__ ) ;
  //
  const double w  = sumw  () ;
  const double w2 = sumw2 () ;
  //
  Ostap::Assert (  0 < w   , 
                   "Invalid sum of weights!"            , 
                   "Ostap::Math::WECDF::quantile"       , 
                   INVALID_SUMW  , __FILE__  , __LINE__ ) ;
  //
  Ostap::Assert (  0 < w2  , 
                   "Invalid sum of squared weights!"    , 
                   "Ostap::Math::WECDF::quantile"       , 
                   INVALID_SUMW2 , __FILE__  , __LINE__ ) ;
  //
  if ( 1 == size() ) { return m_data.front().first ; }
  //
  const Ostap::Math::WHarrellDavis estimator { false } ;
  return estimator.quantile ( m_data.begin () , m_data.end () , p , w , w2 ) ;
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
double Ostap::Math::WECDF::quantile
( const double                                    p     ) const
{
  const Ostap::QuantileTypes::HarrellDavisType t {} ;
  return quantile ( p , t ) ;
}
// ============================================================================
// project data content into 1D histogram[
// ============================================================================
void Ostap::Math::WECDF::project 
( TH1& histo ) const 
{
  Ostap::Assert ( 1 == histo.GetDimension()      , 
                  "Invalild histogram dimension" , 
                  "Ostap::Math::ECDF::project"   , 
                  INVALID_TH1                    ) ; 
  histo.Reset() ;
  if ( !histo.GetSumw2() ) { histo.Sumw2() ; }
  for ( const auto& item : m_data )
  { histo.Fill ( item.first , item.second ) ; }
}
// ============================================================================
// low_edge <= xmin 
// ============================================================================
double Ostap::Math::WECDF::low_edge  () const
{ return m_data.empty() ? -std::numeric_limits<double>::max () : xmin () ; } 
// ============================================================================
//  xmax < high_edge 
// ============================================================================
double Ostap::Math::WECDF::high_edge () const
{
  return m_data.empty() ?
    +std::numeric_limits<double>::max () :
    Ostap::Math::high_edge ( xmax () ) ;
}
  

// ============================================================================
//                                                                      The END 
// ============================================================================

