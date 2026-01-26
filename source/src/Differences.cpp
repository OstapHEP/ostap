// $Id$
// ============================================================================
// Include files 
// ============================================================================
// STD& STL 
// ============================================================================
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <memory>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Differences.h"
// ============================================================================
/** @file
 *  Implementation file for functiond form namesapce Ostap::Math::Differences
 *  @see Ostap::Math::Differences
 *  @date 2011-09-27 
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 */
// ============================================================================
/*  Evaluate N-th forward difference of function <code>fun</code>
 *  @param fun the function 
 *  @param N   the order of the difference 
 *  @param x   the point 
 *  @param h   the step 
 *  @return N-th forward dirrerence 
 */
// ============================================================================
double Ostap::Math::Differences::forward 
( std::function<double(double)> fun , 
  const unsigned short          N   , 
  const double                  x   , 
  const double                  h   )
{ return Ostap::Math::Differences::forward_ ( fun , N , x , h ) ; }
// ============================================================================
/*  Evaluate N-th backward difference of function <code>fun</code>
 *  @param fun the function 
 *  @param N   the order of the difference 
 *  @param x   the point 
 *  @param h   the step 
 *  @return N-th backward dirrerence 
 */
// ============================================================================
double Ostap::Math::Differences::backward 
( std::function<double(double)> fun , 
  const unsigned short          N   , 
  const double                  x   , 
  const double                  h   )
{ return Ostap::Math::Differences::backward_ ( fun , N , x , h ) ; }
// ============================================================================
/*  Evaluate N-th central difference of function <code>fun</code>
 *  @param fun the function 
 *  @param N   the order of the difference 
 *  @param x   the point 
 *  @param h   the step 
 *  @return N-th central dirrerence 
 */
// ============================================================================
double Ostap::Math::Differences::central
( std::function<double(double)> fun , 
  const unsigned short          N   , 
  const double                  x   , 
  const double                  h   )
{ return Ostap::Math::Differences::central_ ( fun , N , x , h ) ; }
// ============================================================================

// ============================================================================
/* construcructor for the order 
 *  @param N  the order of the finite differenece 
 */
// ============================================================================
Ostap::Math::Differences::FiniteDifference::FiniteDifference
( const unsigned short N )
  : m_N ( N )
{}
// ============================================================================
// get forward difference 
// ============================================================================
double Ostap::Math::Differences::FiniteDifference::forward
( std::function<double(double)> f , 
  const double                  x ,
  const double                  h ) const
{
  return Ostap::Math::Differences::forward_ ( f , m_N , x , h ) ;
}
// ============================================================================
// get backward difference 
// ============================================================================
double Ostap::Math::Differences::FiniteDifference::backward 
( std::function<double(double)> f , 
  const double                  x ,
  const double                  h ) const 
{
  return Ostap::Math::Differences::backward_ ( f , m_N , x , h ) ;
}
// ============================================================================
// get central difference 
// ============================================================================
double Ostap::Math::Differences::FiniteDifference::central 
( std::function<double(double)> f , 
  const double                  x ,
  const double                  h ) const 
{
  return Ostap::Math::Differences::central_ ( f , m_N , x , h ) ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================
