// Include files
// ============================================================================
//  STD&STL
// ============================================================================
#include <ostream> 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatusCode.h"
// ============================================================================
/** Implementation file for class Ostap::StatusCode
 *  @see Ostap::StatusCode 
 */
// ============================================================================
// printout 
std::ostream&
Ostap::operator<<
( std::ostream&            os ,
  const Ostap::StatusCode& sc )
{
  if      ( sc.isSuccess     ()                     ) { return os << "SUCCESS"     ; }
  else if ( sc.isRecoverable ()                     ) { return os << "RECOVERABLE" ; }
  else if ( Ostap::StatusCode::FAILURE == sc.code() ) { return os << "FAILRUE"     ; }
  //
  return os << "FAILURE(" << sc.code() << ")" ;
}
// ============================================================================
//                                                                      The END  
// ============================================================================
