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
//                                                                      The END 
// ============================================================================
