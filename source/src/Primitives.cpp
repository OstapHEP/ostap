// ============================================================================
// Include files 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Primitives.h"
// ============================================================================
// Local
// ============================================================================
#include "Exception.h"
#include "local_math.h"
// ============================================================================
/** @file  
 *  Implementaiton file for functions from the file   Ostap/Primitives.h 
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2020-07-25
 */
// ============================================================================
Ostap::Math::Moebius::Moebius
( const double a ,
  const double b ,
  const double c ,
  const double d )
  : Moebius ( a , b, c , d , Id() )
{}
// ==========================================================================
bool Ostap::Math::Moebius::check() const
{
  Ostap::Assert ( !s_zero ( m_a * m_d - m_b * m_c )    ,
                  "invalid parameters!"  , 
                  "Ostap::Math::Moebius" ) ;
  return true ;
}
// ===========================================================================
Ostap::Math::Step::Step
( const double a ,
  const double b )
  : Step ( a , b, Ostap::Math::Id() )
{}


// ============================================================================
//                                                                      The END 
// ============================================================================
