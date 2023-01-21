#ifndef OSTAP_EXCEPTION2_H
#define OSTAP_EXCEPTION2_H
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <string>
#include <exception>
#include <memory>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatusCode.h"
#include "Ostap/Exception.h"
// ============================================================================
namespace Ostap 
{
  // ==========================================================================
  inline bool Assert
  ( const bool         assertion                        ,
    const std::string& message                          ,
    const std::string& tag       = "Ostap"              ,
    const StatusCode&  sc        =  StatusCode::FAILURE )
  { return assertion ? true : throwException ( message , tag , sc ).isSuccess() ; }
  // ===========================================================================
  template < unsigned int N1, unsigned int N2>
  inline bool Assert
  ( const bool        assertion                            ,
    const char        (&message)[N1]                       ,
    const char        (&tag)    [N2]                       ,
    const StatusCode& sc             = StatusCode::FAILURE )
  { return assertion ? true :
      throwException ( std::string ( message , message + N1 ) ,
                       std::string ( tag     , tag     + N2 ) , sc ).isSuccess() ; }
  // ===========================================================================
} //                                                      end of namespace Ostap
// =============================================================================
//                                                                       The END 
// =============================================================================
#endif  // OSTAP_EXCEPTION2_H
// =============================================================================
