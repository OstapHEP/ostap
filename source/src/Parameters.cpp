// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================ \
#include <algorithm>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/Parameters.h"
// ============================================================================
/** @file 
 *  Implementation file for class Ostap::Math::Parameters
 *  @see Ostap::Math::Parameters
 *  @author Vanya BELYAEV
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  static_assert ( std::numeric_limits<double>::is_specialized           , 
                  "mumeric_limits are not specialized for doubles"      ) ;
  // ==========================================================================
  /// equality criteria for doubles
  const Ostap::Math::Equal_To<double> s_equal {} ;       // equality criteria for doubles
  /// zero for doubles  
  const Ostap::Math::Zero<double>     s_zero  {} ;       // zero for doubles
  /// zero fo vectors 
  const Ostap::Math::Zero< std::vector<double> > s_vzero {} ; // zero for vectors
  // ==========================================================================
}
// ============================================================================
// PARAMETERS
// ============================================================================
// constructor from number of parameters 
// ============================================================================
Ostap::Math::Parameters::Parameters 
( const std::size_t np ) 
  : m_pars ( np , 0.0 ) 
{}
// ============================================================================
// constructor from  the list of parameters 
// ============================================================================
Ostap::Math::Parameters::Parameters 
( const std::vector<double>&  pars   ) 
  : m_pars ( pars ) 
{}
// ============================================================================
// constructor from the list of parameters 
// ============================================================================
Ostap::Math::Parameters::Parameters 
(       std::vector<double>&& pars   ) 
  : m_pars ( std::forward<std::vector<double> > ( pars ) ) 
{}
// ============================================================================
// all zero ?
// ============================================================================
bool Ostap::Math::Parameters::zero  () const { return s_vzero ( m_pars ) ; }
// ============================================================================
// set k-parameter
// ============================================================================
bool Ostap::Math::Parameters::_setPar 
( const std::size_t k     , 
  const double      value ,
  const bool        force ) 
{
  if ( m_pars.size() <= k                         ) { return false ; }
  if ( s_equal ( m_pars [ k ] , value ) && !force ) { return false ; }
  m_pars [ k ] = value ;
  return true ;
}
// ============================================================================
// rest all parameters to zero 
// ============================================================================
void Ostap::Math::Parameters::reset () 
{ std::fill ( m_pars.begin() , m_pars.end() ,  0.0 ) ; }
// ============================================================================
// swap two objects 
// ============================================================================
void Ostap::Math::Parameters::swap (  Ostap::Math::Parameters& right ) 
{ std::swap ( m_pars ,  right.m_pars ); }
// ============================================================================
/*  filter out very small terms/parameters 
 *  the term is considered to be very small if
 *   - it is numerically zero: \f$ c_k \approx 0 \f% 
 *   - or if epsilon  > 0: \f$ \left| c_k \right| \le \epsilon \f$ 
 *   - or if scale   != 0: \f$ \left| s \right| + \left| c_k \right| \approx \left| s \right| \f$ 
 *  @param  epsilon  parameter to define "smalness" of terms
 *  @param  scale    parameter to define "smalness" of terms
 *  @return number of nullified terms
 */
// ============================================================================
std::size_t
Ostap::Math::Parameters::remove_noise
( const double epsilon , 
  const double scale   )
{
  std::size_t  num    = 0                  ;
  const bool   eps    = 0 <  epsilon       ;
  const bool   sca    = scale              ;
  const double ascale = std::abs ( scale ) ;
  //
  const std::size_t NN { m_pars.size() } ;
  for ( std::size_t k = 0 ; k < NN ; ++k )
    {
      const double absp = std::abs ( m_pars [ k ] ) ;
      if      ( s_zero ( absp )                           ) { m_pars [ k ] = 0 ; ++num ; }
      else if ( eps && absp <= epsilon                    ) { m_pars [ k ] = 0 ; ++num ; }
      else if ( sca && s_equal ( ascale + absp , ascale ) ) { m_pars [ k ] = 0 ; ++num ; }
    }
  //
  return num ;
}    
// ============================================================================


// ============================================================================
// join two vectors togather 
// ============================================================================
std::vector<double>
Ostap::Math::Parameters::join
( const std::vector<double>& a ,
  const std::vector<double>& b )
{
  //
  if      ( a.empty() ) { return b ; } 
  else if ( b.empty() ) { return a ; } 
  //
  std::vector<double> ab ( a.size () + b.size() ) ;
  //
  std::copy ( a.begin () , a.end () , ab.begin()            ) ;
  std::copy ( b.begin () , b.end () , ab.begin() + a.size() ) ;
  //
  return ab ;
}
// ============================================================================
// join scalar and vector together
// ============================================================================
std::vector<double>
Ostap::Math::Parameters::join
( const double               a ,
  const std::vector<double>& b )
{
  //
  if ( b.empty() ) { return { a } ; }
  //
  std::vector<double> ab ( 1 + b.size() ) ;
  //
  ab[0] = a ;
  std::copy ( b.begin () , b.end () , ab.begin() + 1 ) ;
  //
  return ab ;
}
// ============================================================================
// join 2 scalars and vector together
// ============================================================================
std::vector<double>
Ostap::Math::Parameters::join
( const double               a1 ,
  const double               a2 ,
  const std::vector<double>& b  )
{
  //
  if ( b.empty() ) { return { a1 , a2 } ; }
  //
  std::vector<double> ab ( 2 + b.size() ) ;
  //
  ab [ 0 ] = a1 ;
  ab [ 1 ] = a2  ;
  //
  std::copy ( b.begin () , b.end () , ab.begin() + 2 ) ;
  //
  return ab ;
}
// ============================================================================
// join 3 scalars and vector together
// ============================================================================
std::vector<double>
Ostap::Math::Parameters::join
( const double               a1 ,
  const double               a2 ,
  const double               a3 ,
  const std::vector<double>& b  )
{
  //
  if ( b.empty() ) { return { a1 , a2 , a3 } ; }
  //
  std::vector<double> ab ( 3 + b.size() ) ;
  //
  ab [ 0 ] = a1 ;
  ab [ 1 ] = a2  ;
  ab [ 2 ] = a3  ;
  //
  std::copy ( b.begin () , b.end () , ab.begin() + 3 ) ;
  //
  return ab ;
}
// ============================================================================
// join 4 scalars and vector together
// ============================================================================
std::vector<double>
Ostap::Math::Parameters::join
( const double               a1 ,
  const double               a2 ,
  const double               a3 ,
  const double               a4 ,
  const std::vector<double>& b  )
{
  //
  if ( b.empty() ) { return { a1 , a2 , a3 , a4 } ; }
  //
  std::vector<double> ab ( 4 + b.size() ) ;
  //
  ab [ 0 ] = a1 ;
  ab [ 1 ] = a2  ;
  ab [ 2 ] = a3  ;
  ab [ 3 ] = a4  ;
  //
  std::copy ( b.begin () , b.end () , ab.begin() + 4 ) ;
  //
  return ab ;
}
// ============================================================================
//  join vector & scalar together
// ============================================================================
std::vector<double>
Ostap::Math::Parameters::join
( const std::vector<double>& a , 
  const double               b )
{
  //
  if ( a.empty() ) { return { b } ; }
  //
  std::vector<double> ab ( a.size() + 1 ) ;
  //
  std::copy ( a.begin () , a.end () , ab.begin() ) ;
  ab [ a.size() ] = b ;
  //
  return ab ;
}
// ============================================================================
//  join vector & 2 scalars together
// ============================================================================
std::vector<double>
Ostap::Math::Parameters::join
( const std::vector<double>& a  , 
  const double               b1 , 
  const double               b2 ) 
{
  //
  if ( a.empty() ) { return { b1 , b2 } ; }
  //
  std::vector<double> ab ( a.size() + 2 ) ;
  //
  std::copy ( a.begin () , a.end () , ab.begin() ) ;
  ab [ a.size()     ] = b1 ;
  ab [ a.size() + 1 ] = b2 ;
  //
  return ab ;
}
// ============================================================================
//  join vector & 3 scalars together
// ============================================================================
std::vector<double>
Ostap::Math::Parameters::join
( const std::vector<double>& a  , 
  const double               b1 , 
  const double               b2 , 
  const double               b3 ) 
{
  //
  if ( a.empty() ) { return { b1 , b2 , b3 } ; }
  //
  std::vector<double> ab ( a.size() + 3 ) ;
  //
  std::copy ( a.begin () , a.end () , ab.begin() ) ;
  ab [ a.size()     ] = b1 ;
  ab [ a.size() + 1 ] = b2 ;
  ab [ a.size() + 2 ] = b3 ;
  //
  return ab ;
}
// ============================================================================
//  join vector & 4 scalars together
// ============================================================================
std::vector<double>
Ostap::Math::Parameters::join
( const std::vector<double>& a  , 
  const double               b1 , 
  const double               b2 , 
  const double               b3 , 
  const double               b4 ) 
{
  //
  if ( a.empty() ) { return { b1 , b2 , b3 , b4 } ; }
  // 
  std::vector<double> ab ( a.size() + 4 ) ;
  //
  std::copy ( a.begin () , a.end () , ab.begin() ) ;
  ab [ a.size()     ] = b1 ;
  ab [ a.size() + 1 ] = b2 ;
  ab [ a.size() + 2 ] = b3 ;
  ab [ a.size() + 4 ] = b4 ;
  //
  return ab ;
}




// ============================================================================
//                                                                      The END 
// ============================================================================
