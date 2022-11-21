// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================
#include <functional>
#include <algorithm>
#include <vector>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Piecewise.h"
// ============================================================================
// local
// ============================================================================
#include "Exception.h"
// ============================================================================
/** @file
 *  Implementation file for class Ostap::Math::Piecewise
 *  @date 2020-06-29 
 *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
 */
// ============================================================================
namespace
{
  // ==========================================================================
  const std::string s_ERROR   { "Cannot insert new value!"    } ;
  const std::string s_METHOD  { "Ostap::Math::Piecewise::add" } ;
  const std::string s_METHOD3 { "Ostap::Math::Piecewise(3)"   } ;
  const std::string s_METHOD4 { "Ostap::Math::Piecewise(4)"   } ;
  // ==========================================================================
}
// ============================================================================
// default constructor : create constant function
// ============================================================================
Ostap::Math::Piecewise::Piecewise 
( const double value ) 
  : m_edges() 
  , m_funcs ( 1 , FPAIR ( []( const double /* x */ ) -> double { return 1 ; } , value ) )
{}
// ============================================================================
/*  add new function, defined for  \f$ x\ge x_i\$ 
 *  @attention xi must be larger than any previosly aded ranges!
 *  @param xi value of \f$ x_i \f$
 *  @param fi the function to be used for \f$x\ge x_i\f$ 
 */
// ============================================================================
void Ostap::Math::Piecewise::add_ 
( const double                  xi , 
  std::function<double(double)> fi ,
  const double                  si ) 
{
  // remove the entries that are not actual anymore 
  while ( 1 <= m_edges.size() && xi < m_edges.back() ) 
  {
    m_edges.pop_back() ;
    m_funcs.pop_back() ;
  }
  //
  Ostap::Assert ( m_edges.empty () || xi > m_edges.back() , s_ERROR , s_METHOD ) ;
  //
  m_edges.push_back ( xi ) ; 
  m_funcs.push_back ( FPAIR ( fi , si ) ) ;
}
// ============================================================================
// scale the function 
// ============================================================================
Ostap::Math::Piecewise&
Ostap::Math::Piecewise:: operator*= ( const double value ) 
{
  for ( auto& p : m_funcs ) { p.second *= value ; }
  return *this ;
}
// ============================================================================
// scale the function 
// ============================================================================
Ostap::Math::Piecewise&
Ostap::Math::Piecewise:: operator/= ( const double value ) 
{ (*this) *= ( 1.0 / value ) ; return *this ; }


// ============================================================================
//                                                                      The END 
// ============================================================================

