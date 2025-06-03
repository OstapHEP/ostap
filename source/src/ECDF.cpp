// ============================================================================
// Incldue files 
// ============================================================================
// STD&STL
// ============================================================================
#include <memory>
#include <cstring>
#include <numeric>
#include <algorithm>
#include <tuple>
#include <array>
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/ECDF.h"
#include "Ostap/Power.h"
#include "Ostap/MoreMath.h"
#include "Ostap/WStatEntity.h"
// ============================================================================
// Local
// ============================================================================
#include "Exception.h"
#include "status_codes.h"
#include "local_math.h"
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
// check that ECDF is OK: there are some entreiss 
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

