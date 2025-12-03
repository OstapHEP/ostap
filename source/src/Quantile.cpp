// ============================================================================
// Include files 
// ============================================================================
//  STD&STL
// ============================================================================
#include <algorithm>
#include <limits>
#include <cmath> 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Quantile.h"
#include "Ostap/Quantiles.h"
#include "Ostap/StatusCode.h"
// ============================================================================
// local
// ============================================================================
#include "status_codes.h"
#include "local_math.h" 
// ============================================================================
/** @file  
 *  Implementation file for class Ostap::Math::P2Quantile
 *  @date 2020-05-21 
 *  @author Vanya BELYAEV Ivan/Belyaev@itep.ru
 */
// ============================================================================
namespace
{
  // ==========================================================================
  /// Parabolic interpolator for P2 algorithm 
  template <class VN, class VQ>
  inline double Parabolic
  ( const VN&         n ,
    const VQ&         q ,
    const std::size_t i , 
    const double      d )
  {
    //
    const double n2  = n [ i + 1 ] * 1.0 - n [ i - 1 ] ;
    const double nm1 = n [ i     ] * 1.0 - n [ i - 1 ] ; 
    const double np1 = n [ i + 1 ] * 1.0 - n [ i     ] ; 
    const double qp1 = q [ i + 1 ] - q [ i     ] ; 
    const double qm1 = q [ i     ] - q [ i - 1 ] ;
    //
    return q [ i ] + d / n2 * ( ( nm1 + d ) * qp1 / np1 + 
                                ( np1 - d ) * qm1 / nm1 ) ; 
  }
  // ==========================================================================
  /// Linear interpolator for P2 algorithm 
  template <class VN, class VQ>
  inline double Linear
  ( const VN&         n ,
    const VQ&         q ,
    const std::size_t i , 
    const int         d )
  {
    const double dq = q [ i + d ]       - q [ i ] ;
    const double dn = n [ i + d ] * 1.0 - n [ i ] ;
    //
    return q [ i ] + d * dq / dn ;
  }
  // =========================================================================
  /// Adjust N&Q
  template <class VN, class VQ, class VS>
  inline void Adjust 
  ( VN&               n  ,
    VQ&               q  ,
    const VS&         ns ,
    const std::size_t i  ) 
  {  
    const double dd = ns [ i ] - n [ i ] ;
    //
    const long long dnp = n [ i + 1 ] * 1L - n [ i ] ;
    const long long dnm = n [ i - 1 ] * 1L - n [ i ] ;
    //
    if ( (  1 <= dd &&  1 < dnp ) || ( -1 >= dd && -1 > dnm )  )  
      {
        //
        const int di =  ( 0 <= dd ) ? 1 : -1 ;
        //
        const double qp = Parabolic ( n , q , i , di ) ;
        if ( std::isfinite ( qp ) && q [ i - 1 ] < qp && qp < q [ i + 1 ] ) { q [ i ] = qp                        ; }
        else                                                                { q [ i ] = Linear ( n , q , i , di ) ; }
        //
        n [ i ] += di ;
      } 
  } ;
  // ========================================================================
  template <class VP, class VS>
  inline void UpdateNs
  ( const VP& p  ,
    VS&       ns ,
    const Ostap::Math::Quantile::size_type maxIndex )
  {
    // number of probabilities 
    const std::size_t M = p .size() ; 
    const std::size_t N = ns.size() ; 
    //
    // Principal markers
    ns.front () = 0;
    for ( std::size_t i = 0 ; i < M ; ++i )
      { ns [ i * 2 + 2 ] = maxIndex * p [ i ] ; }
    ns.back() = maxIndex ;
    //
    // Middle markers    
    ns [ 1 ] = maxIndex * p.front () * 0.5 ;
    for ( std::size_t i = 1; i < M ; ++i )
      {
	const double pr = 0.5 * ( p [ i - 1 ] + p [ i ] ) ;
	ns [2 * i + 1 ] = maxIndex * pr ;
      }
    ns [ N - 2 ] = maxIndex * ( 1 + p.back() ) * 0.5 ;
  } ;
  // =========================================================================
} //                                            The end of anonymois namespace
// ===========================================================================

// ===========================================================================
// P2 algorithm for (approximate) quantile estimamtion  
// ============================================================================

// ============================================================================
// contrcuctor 
// =============================================================================
Ostap::Math::Quantile::Quantile
( const double                                p , 
  const Ostap::Math::Quantile::Initialization s ) 
  : m_init    ( s ) 
  , m_p       ( p ) 
  , m_N       ( 0 )
  , m_q       {}
  , m_ns      {}
  , m_n       { 0 , 1 , 2 , 3 , 4 }
  , m_counter {}  
{
  Ostap::Assert ( Classic  == s || Adaptive == s          ,
                  "Invalid Initialization strategy"       ,  
                  "Ostap::Math::Quantile"                 ,
                  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //		  
  Ostap::Assert ( 0 < m_p && m_p < 1                         ,  
                  "Invalid quantile probability"             ,  
                  "Ostap::Math::Quantile"                    , 
                  INVALID_PROBABILITY , __FILE__ , __LINE__  ) ;
  //
  // Classic initialization strategy 
  m_ns = std::array<double,5> { 0 , 2 * m_p , 4 * m_p , 2 * m_p + 2 , 4 } ;
  //
}
// =============================================================================
// reset quantile 
// =============================================================================
void Ostap::Math::Quantile::reset  ()
{
  m_N = 0 ;
  m_counter.reset() ;
} 
// =============================================================================
// update the counter 
// =============================================================================
Ostap::Math::Quantile&
Ostap::Math::Quantile::add  
( const double value ) 
{
  /// (1) skip non-finite values 
  if ( !std::isfinite ( value ) ) { return *this ; }
  /// (2) update the helper counter 
  m_counter.add ( value ) ;
  //
  if ( m_N  < 5 )
    {
      m_q [ m_N ] = value ;
      ++m_N ;
      /// explicitely sort them 
      std::sort ( m_q.begin () , m_q.begin() + m_N ) ;
      if ( 5 == m_N && Adaptive == m_init )
        {
          //
          std::copy ( m_q.begin() , m_q.end() , m_ns.begin() ) ;
          //
          m_n [ 1 ] = Ostap::Math::round ( 2 * m_p     ) ;
          m_n [ 2 ] = Ostap::Math::round ( 4 * m_p     ) ;
          m_n [ 3 ] = Ostap::Math::round ( 2 * m_p + 2 ) ;
          //
          m_q [ 1 ] = m_ns [ m_n [ 1 ] ] ;
          m_q [ 2 ] = m_ns [ m_n [ 2 ] ] ;
          m_q [ 3 ] = m_ns [ m_n [ 3 ] ] ;
          // 
          m_ns = std::array<double,5> { 0 , 2 * m_p , 4 * m_p , 2 * m_p + 2 , 4  } ;
        }      
      //
      return *this ; 
    }
  //
  const std::size_t k =     
    value < m_q [ 0 ] ? 0 :
    value < m_q [ 1 ] ? 0 :
    value < m_q [ 2 ] ? 1 :
    value < m_q [ 3 ] ? 2 : 3 ;
  //
  /// adjust minimal & maximal 
  if      ( 0 == k ) { m_q .front () = std::min ( m_q.front () , value ) ; } 
  else if ( 3 <= k ) { m_q .back  () = std::max ( m_q.back  () , value ) ; } 
  //
  /// increment 
  for ( std::size_t i = k + 1 ; i < 5 ; ++i ) { ++m_n [ i ] ; }
  //
  /// adjust ns 
  m_ns [ 1 ] = m_N *       m_p   / 2 ;
  m_ns [ 2 ] = m_N *       m_p       ;
  m_ns [ 3 ] = m_N * ( 1 + m_p ) / 2 ;
  m_ns [ 4 ] = m_N                   ;
  //
  // ==============================================================
  // Adjust the counter s
  // =================================================================  
  // Adjust 
  if ( 0.5 <= m_p )
    {
      Adjust ( m_n , m_q , m_ns , 1 ) ;
      Adjust ( m_n , m_q , m_ns , 2 ) ;
      Adjust ( m_n , m_q , m_ns , 3 ) ;
    }
  else
    {
      Adjust ( m_n , m_q , m_ns , 3 ) ;
      Adjust ( m_n , m_q , m_ns , 2 ) ;
      Adjust ( m_n , m_q , m_ns , 1 ) ;
    }
  //
  ++m_N ;
  //
  return *this ;
}
// ============================================================================
// swap two objects
// ============================================================================
void Ostap::Math::Quantile::swap
( Ostap::Math::Quantile& right )
{
  std::swap ( m_init , right.m_init ) ;
  std::swap ( m_p    , right.m_p    ) ;
  std::swap ( m_N    , right.m_N    ) ;
  std::swap ( m_q    , right.m_q    ) ;
  std::swap ( m_ns   , right.m_ns   ) ;
  std::swap ( m_n    , right.m_n    ) ;
  //
  m_counter.swap  ( right.m_counter ) ;
}
// ============================================================================
// get the quantile value
// ============================================================================
double Ostap::Math::Quantile::quantile () const
{
  // regular case 
  if ( 5 < m_N ) { return m_q [ 2 ] ; }  
  //
  /// if there are not enougth entry, use the explicit quantile  estimator
  const Ostap::Math::HarrellDavis hd {} ;
  return hd ( m_q.begin() , m_q.begin() + m_N , m_p ) ; 
}
// ============================================================================
// get the triplet: (min. quantile, max)
// ============================================================================
std::array<double,3> Ostap::Math::Quantile::quantiles () const
{ return {{ this->min () , this->quantile() , this->max() }} ; }
// ============================================================================

// ============================================================================
// Extended P2 algorithm for (approximate) quantile estimamtion
// ============================================================================

// ============================================================================
// Constructor from vector of probabiliies 
// ============================================================================
Ostap::Math::Quantiles::Quantiles
( const std::vector<double>& p )
  : m_p  (   )
  , m_N  ( 0 )
  , m_q  (   )
  , m_ns (   )
  , m_n  (   )    
{
  //
  // use only good/valid probabilities
  //
  m_p.reserve  ( p.size  () ) ;
  std::copy_if ( p.begin () ,
                 p.end   () ,
                 std::back_inserter ( m_p ) , 
                 []( const double v ) -> bool
                 { return 0 < v && v < 1 ; } ) ;
  //
  std::sort ( m_p.begin () , m_p.end () );
  m_p.erase ( std::unique ( m_p.begin () , m_p.end () ) , m_p.end () ) ;
  //
  Ostap::Assert ( !m_p.empty ()                              ,
                  "Empty array of probabilities"             ,
                  "Ostap::Math::Quantiles"                   ,
                  INVALID_PROBABILITIES, __FILE__ , __LINE__ ) ;
  //
  const std::size_t MC = 2 * m_p.size() + 3 ;
  //
  m_q  . resize ( MC , 0.0 ) ;
  m_ns . resize ( MC , 0.0 ) ;
  m_n  . resize ( MC , 0   ) ;
  //
}
// ============================================================================
// from quantile index 
// ============================================================================
Ostap::Math::Quantiles::Quantiles
( const std::size_t N  )
  : m_p  (   )
  , m_N  ( 0 )
  , m_q  (   )
  , m_ns (   )
  , m_n  (   )    
  , m_counter ()
{
  // =========================================================================
  Ostap::Assert ( 1 <= N ,
                  "Invalid quantile index "                     ,
                  "Ostap::Math::Quantiles"                      ,
                  INVALID_QUANTILE_INDEX  , __FILE__ , __LINE__ ) ;
  //
  m_p.resize ( N ) ;
  for ( std::size_t i = 0 ; i < N ; ++ i )
    { m_p [ i ] = ( i + 1.0 ) / ( N + 1 ) ; }
  //
  const std::size_t MC = 2 * m_p.size() + 3 ;
  //
  m_q  . resize ( MC , 0.0 ) ;
  m_ns . resize ( MC , 0.0 ) ;
  m_n  . resize ( MC , 0   ) ;
  // ===========================================================================
}
// =============================================================================
// reset quantile 
// =============================================================================
void Ostap::Math::Quantiles::reset  ()
{
  m_N = 0 ;
  m_counter.reset() ;
}
// =============================================================================
// add the value 
// =============================================================================
Ostap::Math::Quantiles&
Ostap::Math::Quantiles::add
( const double value )
{
  /// (1) skip non-finite values; 
  if ( !std::isfinite ( value ) ) { return *this ; }  
  /// (2) update internal helper counter (mainly for debugging) 
  m_counter.add ( value ) ;
  //
  const std::size_t MC  = m_ns.size() ;
  //
  if ( m_N < MC )
    {
      m_q [ m_N ] = value ;
      ++m_N ;
      // keep it sorted 
      std::sort ( m_q.begin() , m_q.begin() + m_N  ) ;
      //
      if ( m_N == MC )
        {
          //
          UpdateNs ( m_p , m_ns , MC - 1 ) ;
          //
          for ( std::size_t i = 0 ; i < MC  ; ++i )
            { m_n [ i ] = Ostap::Math::round ( m_ns [ i ] ) ; }
          //
          std::copy ( m_q.begin() , m_q.end() , m_ns.begin() ) ;
          //
          for ( std::size_t i = 0 ; i < MC ; ++i)
            { m_q [ i ] = m_ns [ m_n [ i ] ] ; }
          //
          UpdateNs ( m_p , m_ns , MC - 1 ) ;	  
        }
      //
      return *this ;
      // ==============================================================
    }
  // ==================================================================
  const auto ik = std::upper_bound ( m_q.begin() , m_q.end() , value ) ;
  
  const std::size_t k =
    m_q.begin () == ik ? 0      :
    m_q.end   () == ik ? MC - 2 :
    ( ik - m_q.begin() ) - 1 ;
  
  //
  if      ( 0     == k  ) { m_q.front () = std::min ( m_q.front () , value ) ; } 
  else if ( k + 2 >= MC ) { m_q.back  () = std::max ( m_q.back  () , value ) ; } 
  //
  for ( std::size_t i = k + 1 ; i < MC ; ++i ) { ++m_n [ i ] ; }
  //
  UpdateNs ( m_p , m_ns , m_N ) ;  
  //
  std::size_t left  = 1      ;
  std::size_t right = MC - 2 ; 
  while ( left <= right )
    {
      const double l = std::abs ( m_ns [ left  ] / m_N - 0.5 ) ;
      const double r = std::abs ( m_ns [ right ] / m_N - 0.5 ) ;
      //
      const std::size_t index = ( l <= r ) ? left++ : right-- ;
      Adjust ( m_n , m_q , m_ns , index ) ;
    }
  //
  ++m_N ;
  //
  return *this ;
}
// ======================================================================
// swap two objects
// ============================================================================
void Ostap::Math::Quantiles::swap
( Ostap::Math::Quantiles& right )
{
  std::swap ( m_p    , right.m_p    ) ;
  std::swap ( m_N    , right.m_N    ) ;
  std::swap ( m_q    , right.m_q    ) ;
  std::swap ( m_ns   , right.m_ns   ) ;
  std::swap ( m_n    , right.m_n    ) ;
  //
  m_counter.swap ( right.m_counter ) ;
}
// ============================================================================
// get all quantiles 
// ============================================================================
std::vector<double>
Ostap::Math::Quantiles::quantiles () const
{
  std::vector<double> result {} ;
  result.reserve   ( m_p.size() + 2 ) ; // ATTENTION!!! 
  result.push_back ( this->min()    ) ; // ADD MIN
  //
  if (  m_N <= m_q.size() )
    {
      Ostap::Math::HarrellDavis hd { false } ;
      for ( const auto p : m_p )
        { result.push_back ( hd ( m_q.begin() , m_q.begin() + m_N , p ) ) ; }
      //
      result.push_back ( this->max() ) ; // ADD MAX 
      return result ; 
    }
  //
  const std::size_t M = NP () ;
  for ( std::size_t i = 0 ; i < M ; ++i  )
    {
      const std::size_t index = i * 2 + 2 ; 
      result.push_back ( m_q [ index ] ) ;
    }
  //
  result.push_back ( this->max() ) ; // ADD MAX
  return result ; 
}
// ============================================================================
// get the quantile value
// ============================================================================
double Ostap::Math::Quantiles::quantile
( const std::size_t index ) const
{
  if ( m_p.size () <= index ) { return max () ; }
  //  
  /// regular case  
  if ( m_q.size () < m_N    ) { return m_q [ index * 2 + 2 ] ; }
  /// if there are not enought quantiles, use the explicit algorithms 
  const Ostap::Math::HarrellDavis hd { false } ;
  return hd ( m_q.begin() , m_q.begin() + m_N , m_p [ index ] ) ;  
}
// ============================================================================
//                                                                      The END  
// ============================================================================
