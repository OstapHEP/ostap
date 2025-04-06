// ============================================================================
// Incldue files
// ============================================================================
// STD & STL
// ===========================================================================
#include <array>
#include <limits>
#include <map>
#include <vector>
// ===========================================================================
// Ostap
// ===========================================================================
#include "Ostap/Choose.h"
#include "Ostap/Clenshaw.h"
#include "Ostap/StatusCode.h"
#include "Ostap/Bernulli.h"
#include "Ostap/Polynomials.h"
// ===========================================================================
// Local
// ===========================================================================
#include "status_codes.h"
#include "syncedcache.h"
#include "local_math.h"
// ===========================================================================
/** @file 
 *  Implementation of functions and classes from namespace Ostap::Math::Bernulli
 *  @authr Vanya BELYAEV Ivam.Belyave@cern,.ch
 *  @date 2025-04-05
 */
// ===========================================================================
namespace Ostap
{
  // =========================================================================
  namespace Math
  {
    // =======================================================================
    /** Riemann's Zeta function \f$ n\ne 1\f$:
     *  \f$ \zeta ( n ) = \sum_k k^{-n}\f$ 
     */
    double zeta ( const int    n ) ;
    // =======================================================================
  }
  // ========================================================================
}
// ===========================================================================
namespace
{
  // =========================================================================
  static_assert ( std::numeric_limits<long double>::is_specialized             ,
                  "std::numeric_limits<long double> is not specialized"        ) ;
  static_assert ( std::numeric_limits<long long>::is_specialized               ,
                  "std::numeric_limits<long long> is not specialized"          ) ;
  static_assert ( std::numeric_limits<unsigned long long>::is_specialized      ,
                  "std::numeric_limits<unsigned long long> is not specialized" ) ;
  // =========================================================================
  const std::array<long double ,36> s_B { 
    +1 *1.0L / 1, 
    -1 *1.0L / 2, 
    +1 *1.0L / 6, 
    +0 *1.0L / 1, 
    -1 *1.0L / 30, 
    +0 *1.0L / 1, 
    +1 *1.0L / 42, 
    +0 *1.0L / 1, 
    -1 *1.0L / 30, 
    +0 *1.0L / 1, 
    +5 *1.0L / 66, 
    +0 *1.0L / 1, 
    -691 *1.0L / 2730, 
    +0 *1.0L / 1, 
    +7 *1.0L / 6, 
    +0 *1.0L / 1, 
    -3617 *1.0L / 510, 
    +0 *1.0L / 1, 
    +43867 *1.0L / 798, 
    +0 *1.0L / 1, 
    -174611 *1.0L / 330, 
    +0 *1.0L / 1, 
    +854513 *1.0L / 138, 
    +0 *1.0L / 1, 
    -236364091 *1.0L / 2730, 
    +0 *1.0L / 1, 
    +8553103 *1.0L / 6, 
    +0 *1.0L / 1, 
    -23749461029 *1.0L / 870, 
    +0 *1.0L / 1, 
    +8615841276005 *1.0L / 14322, 
    +0 *1.0L / 1, 
    -7709321041217 *1.0L / 510, 
    +0 *1.0L / 1, 
    +2577687858367 *1.0L / 6, 
    +0 *1.0L / 1 } ;
  // ==========================================================================
  // Coefficients of Bernulli polynoilas of order N 
  const std::vector<long double>&
  bernulli_poly
  ( const unsigned short N )
  {
    // ========================================================================
    typedef std::vector<long double>        VALUES ; 
    typedef std::map<unsigned short,VALUES> MAP    ;
    typedef SyncedCache<MAP>                STORE  ;
    // ========================================================================
    /// add six first polynomials explicitly 
    static STORE s_store { MAP { 
        { 0 , { 1.0L ,       } } ,
        { 1 , { 1.0L , -0.5L } } , 
        { 2 , { 1.0L , -1.0L , +1.0L/6 } } ,      
        { 3 , { 1.0L , -1.5L , +0.5L     , +0.0L } } ,      
        { 4 , { 1.0L , -2.0L , +1.0L     , +0.0L , -1.0L / 30 } } ,
        { 5 , { 1.0L , -2.5L , +5.0L / 3 , +0.0L , -1.0L / 6  , +0.0L } } ,
        { 6 , { 1.0L , -3.0L , +2.5L     , +0.0L , -0.5L      , +0.0L , +1.0L / 42 } } } } ; 
    // ========================================================================
    { // look into the store
      // ======================================================================
      STORE::Lock lock ( s_store.mutex() ) ;
      auto it = s_store->find ( N ) ;
      if  ( s_store->end() != it ) { return it->second ; }
      // ======================================================================
    } // ======================================================================
    /// calculate the coefficients:
    VALUES values {} ; values.reserve ( N + 1 ) ;
    for ( std::size_t k = 0 ; k <= N ; ++k )
      {
        long double value = Ostap::Math::bernulli( k ) ;
        if ( k < 30 ) { value *= Ostap::Math::choose        ( N , k ) ; }
        else          { value *= Ostap::Math::choose_double ( N , k ) ; }
        values.push_back ( value ) ;
      }
    // =======================================================================
    { // update the store
      // ======================================================================
      STORE::Lock lock ( s_store.mutex() ) ;
      s_store->insert ( std::make_pair ( N , values ) ) ;
      // ======================================================================
    } // ======================================================================
    // ========================================================================
    { // ======================================================================
      // look into the store
      // ======================================================================
      STORE::Lock lock ( s_store.mutex() ) ;
      auto it = s_store->find ( N ) ;
      Ostap::Assert ( s_store->end() != it    ,
                      "Invalid SyncedCache!"  ,
                      "Ostap::Math::Bernulli" ,
                      INVALID_CACHE , __FILE__ , __LINE__ ) ;
      return it->second ;
    }
    // ========================================================================
  } // The end of helper function 
  // ==========================================================================
} // The end of anonymous namespace
// ============================================================================
/** Get Bernulli number 
 *  - \f$ B_0 = 1 \f$ 
 *  - \f$ B_1 = -\frac{1}{2} \f$ 
 *  - \f$ B_{2k+1) = 0 \f$ 
 */
// ============================================================================
double Ostap::Math::bernulli ( const unsigned short n )
{
  if ( n < s_B.size() ) { return s_B [ n ] ; }
  if ( 1 == n % 2     ) { return 0.0       ; }
  // =
  typedef std::map<unsigned short,double> MAP   ;
  typedef SyncedCache<MAP>                STORE ;
  /// the actual store 
  static  STORE s_store {} ; // the actual store 
  // ==========================================================================
  { // look into the store
    // ========================================================================
    STORE::Lock lock ( s_store.mutex() ) ;
    auto it = s_store->find ( n ) ;
    if  ( s_store->end() != it ) { return it->second ; }
    // ========================================================================
  } // ========================================================================
  // ==========================================================================
  const double result = -n * Ostap::Math::zeta ( 1 - n ) ;
  // ==========================================================================
  { // update the store
    // ========================================================================
    STORE::Lock lock ( s_store.mutex() ) ;
    s_store->insert ( std::make_pair ( n , result ) ) ;
    // ========================================================================
  } // ========================================================================
  // ==========================================================================
  return result ;
  // ==========================================================================
}
// ============================================================================
// constructor from the order 
// ============================================================================
Ostap::Math::Bernulli::Bernulli
( const unsigned short N )
  : m_N ( N )
{}
// ============================================================================
// evaluate the polynomial
// ============================================================================
double Ostap::Math::Bernulli::evaluate
( const double x ) const
{
  const std::vector<long double>& C = bernulli_poly ( m_N ) ;
  return Ostap::Math::Clenshaw::monomial_sum ( C.begin() , C.end() , x ).first ;
}
// ============================================================================
// derivative 
// ============================================================================
double Ostap::Math::Bernulli::derivative
( const double x ) const
{
  if ( 0 == m_N ) { return 0 ; }
  const Bernulli D ( m_N - 1 ) ;
  return m_N * D ( x ) ;
}
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::Bernulli::integral
( const double xmin ,
  const double xmax ) const
{
  if ( s_equal ( xmin , xmax ) ) { return 0 ; }
  const Bernulli I ( m_N + 1 ) ;
  return ( I ( xmax ) - I ( xmin ) ) / ( m_N + 1 ) ;
}
// ============================================================================



// ============================================================================
// convert Bernulli polynomila into regular polynomial
// ============================================================================
Ostap::Math::Polynomial::Polynomial
( const Ostap::Math::Bernulli& bp )
  : Polynomial ( bp.degree() , -1.0 , 1.0 )
{
  const std::vector<long double>& C = bernulli_poly ( bp.degree () ) ;
  Ostap::Assert ( m_pars.size() == C.size()           ,
                  "Invalid static sstructure"         ,
                  "Ostap::Math::Bernulli"             , 
                  INVALID_CACHE , __FILE__ , __LINE__ ) ;
  //
  std::copy ( C.rbegin () , C.rend () , m_pars.begin () ) ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================
