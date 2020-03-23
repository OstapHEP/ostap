// Include files
// ============================================================================
//  STD&STL
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
std::ostream& Ostap::operator<< ( std::ostream& s , const Ostap::StatusCode& sc )
{
  if      ( sc.isSuccess     ()                     ) { return s << "SUCCESS"     ; }
  else if ( sc.isRecoverable ()                     ) { return s << "RECOVERABLE" ; }
  else if ( Ostap::StatusCode::FAILURE == sc.code() ) { return s << "FAILRUE"     ; }
  //
  return s << "FAILURE(" << sc.getCode() << ")" ;
}
// ============================================================================
// The END  
// ============================================================================
