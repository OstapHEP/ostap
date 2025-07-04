// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// =============================================================================
#include <algorithm>
#include <cmath> 
#include <map>
// =============================================================================
// Ostap
// =============================================================================
#include "Ostap/Hash.h"
#include "Ostap/Math.h"
#include "Ostap/MoreMath.h"
#include "Ostap/QuantileTypes.h"
#include "Ostap/Quantile.h"
#include "Ostap/Quantiles.h"
#include "Ostap/StatusCode.h"
// =============================================================================
// Local
// =============================================================================
#include "local_math.h"
#include "local_hash.h"
#include "syncedcache.h"
#include "status_codes.h"
// =============================================================================
namespace 
{
  // ===========================================================================
  const std::size_t s_MAX_CACHE = 100000 ;
  // ===========================================================================
  /** calculate \f$ \ I_{t_1}(\alpha,\beta) - I_{t_2} ( \alpha, \beta) f$, where 
   *  \f$ I_z(x,y) \f$ is normalized inncomplete beta function   
   *  @see Ostap::Math::beta_inc  
   *  - protection is added 
   *  - caching is applied 
   */
  double _WHD_  
  ( const double alpha , 
    const double beta  , 
    const double t1    ,
    const double t2    )  
  {
    /// ;
    typedef std::map<std::size_t,double>  MAP   ;
    typedef SyncedCache<MAP>              CACHE ;
    /// the actual integration cache
    static CACHE   s_cache {}      ; // integration cache
    ///
    if ( s_equal ( t1 , t2 ) ) { return 0 ; }
    //
    const std::size_t key { Ostap::Utils::hash_combiner ( alpha , beta , t1 , t2 ) } ; 
    /// lookup in cache 
    { 
      CACHE::Lock lock { s_cache.mutex () } ;
      auto it = s_cache->find ( key ) ;
      if ( s_cache->end () != it ) { return it->second ; }
    }
    /// calculate it 
    const double pmax = std::max ( alpha , beta ) ;
    const double t    = 0.5 *    ( t1 + t2 ) ;
    const double dt   = std::abs ( t2 - t1 ) ;  
    ///
    const auto fun1 = [alpha,beta,t1,t2]() -> double
    { return Ostap::Math::beta_inc ( alpha , beta , t1 ) - 
	Ostap::Math::beta_inc ( alpha , beta , t2 ) ; } ;
    // 
    const auto fun2 = [alpha,beta,t1,t2]() -> double
    {  
      const double t       = 0.5 *    ( t1 + t2 ) ;
      const double dt      = std::abs ( t2 - t1 ) ;  
      double       result  = std::log  ( dt ) ;
      result              += ( alpha - 1 ) * std::log ( t       ) ;
      result              += ( beta  - 1 ) * std::log ( 1 - t   ) ;
      result              -= Ostap::Math::lnbeta ( alpha , beta ) ;
      result               = std::exp ( result ) ;
      return t2 < t1 ? result : -result ;
    } ; 
    ///
    const double result = pmax < 100 ? fun1 () : fun2() ; 
    //
    /// add result to the cache 
    {
      CACHE::Lock lock { s_cache.mutex () } ;
      if ( s_MAX_CACHE < s_cache->size() ) { s_cache->clear() ; }
      s_cache->insert ( std::make_pair ( key , result ) ) ;
    }
    //
    return result ;
  } ; 
  // ============================================================================
} //                                               The end of anynymous namesapce  
// ==============================================================================
// constructor
// =============================================================================
Ostap::QuantileTypes::ABQuantileType::ABQuantileType 
( const double alpha ,
  const double beta  )
  : m_alpha ( alpha )
  , m_beta  ( beta  )
{
  if ( s_zero  ( m_alpha     ) ) { m_alpha = 0 ; }
  if ( s_zero  ( m_beta      ) ) { m_beta  = 0 ; }
  if ( s_equal ( m_alpha , 1 ) ) { m_alpha = 1 ; }
  if ( s_equal ( m_beta  , 1 ) ) { m_beta  = 1 ; }
  //
  Ostap::Assert ( 0 <= m_alpha && m_alpha <= 1            ,
                  "Invalid alpha!"                        ,
                  "Ostap::Math::ABQ"                      ,
                  INVALID_QUANTILE , __FILE__  , __LINE__ ) ;
  Ostap::Assert ( 0 <= m_beta && m_beta <= 1              ,
                  "Invalid beta!"                         ,
                  "Ostap::Math::ABQ"                      ,
                  INVALID_QUANTILE , __FILE__  , __LINE__ ) ;
}
// ==============================================================================
// constructor 
// ==============================================================================
Ostap::QuantileTypes::HarrellDavisType::HarrellDavisType() {} 
// ==============================================================================
// constructor 
// ==============================================================================
Ostap::QuantileTypes::P2QuantileType::P2QuantileType() {} 
// ==============================================================================

// ==============================================================================
Ostap::Math::QCheck::QCheck
( const  bool check )
  :  m_check ( check)
{}
// ==============================================================================
void Ostap::Math::QCheck::throw_exception 
( const char* message , 
  const char* f       , 
  const long  l       ) const 
{
  Ostap::throwException 
    ( message , 
      "Ostap::Math::QCheck" , 
      INVALID_DATA , f ? f : __FILE__  , 0 <= l ? l : __LINE__ ) ;
}
// =============================================================================
// constructor
// =============================================================================
Ostap::Math::HyndmanFan::HyndmanFan
( const Ostap::QuantileTypes::HyndmanFanType  t     , 
  const bool                                  check ) 
  : m_t      ( t    ) 
  , m_check ( check ) 
{
  Ostap::Assert ( Ostap::QuantileTypes::HyndmanFanType::One <= t  &&
                  t <= Ostap::QuantileTypes::HyndmanFanType::Nine  , 
                  "Invalid QuantileType!" , 
                  "Ostap::Math::HyndmanFan" , 
                  INVALID_QUANTILE  , __FILE__ , __LINE__ ) ;
} 
// ==============================================================================
// constructor
// =============================================================================
Ostap::Math::ABQuantile::ABQuantile
( const double alpha , 
  const double beta  ,  
  const bool   check ) 
  : m_abq   ( alpha , beta )
  , m_check ( check ) 
{}
// ==============================================================================
// constructor
// =============================================================================
Ostap::Math::ABQuantile::ABQuantile
( const Ostap::QuantileTypes::ABQuantileType& abq   , 
  const bool                                  check ) 
  : m_abq   ( abq   ) 
  , m_check ( check ) 
{}
// =============================================================================
// constructor
// =============================================================================
Ostap::Math::HarrellDavis::HarrellDavis 
( const bool check ) 
  : m_check  ( check )
{} ;
// =============================================================================
/*  calculate \f$ \ I_{t_1}(\alpha,\beta) - I_{t_2} ( \alpha, \beta) f$, where 
 *  \f$ I_z(x,y) \f$ is normalized incomplete beta function   
 *  @see Ostap::Math::beta_inc  
 *  - protection is added 
 *  - cachinng is applie 
 */
// =============================================================================
double Ostap::Math::HarrellDavis::WHD 
( const std::size_t N  , 
  const double      p  , 
  const double      t1 ,
  const double      t2 ) 
{
  const double alpha = ( N + 1 ) * p ;
  const double beta  = ( N + 1 ) * ( 1 - p ) ;
  return _WHD_  ( alpha , beta , t1  , t2  ) ;
}
// ============================================================================= 


// =============================================================================
// constructor
// =============================================================================
Ostap::Math::WHarrellDavis::WHarrellDavis 
( const bool check ) 
  : m_check  ( check )
{} ;
// =============================================================================
/*  calculate \f$ \ I_{t_1}(\alpha,\beta) - I_{t_2} ( \alpha, \beta) f$, where 
 *  \f$ I_z(x,y) \f$ is normalized incomplete beta function   
 *  @see Ostap::Math::beta_inc  
 *  - protection is added 
 *  - cachinng is applie 
 */
// =============================================================================
double Ostap::Math::WHarrellDavis::WHD
( const double alpha , 
  const double beta  , 
  const double t1    ,
  const double t2    ) 
{
  return _WHD_  ( alpha , beta , t1  , t2  ) ;
}
// ============================================================================= 

// =============================================================================
//                                                                       The END
// ============================================================================= 
