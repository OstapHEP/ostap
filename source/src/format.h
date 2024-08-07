// ============================================================================
#ifndef OSTAP_FORMAT_H 
#define OSTAP_FORMAT_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <string>
// ============================================================================
namespace Ostap
{
  // ============================================================================
  /** very ugly way to avoid boost::format 
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   */
  // ============================================================================
  /// format single number ..
  std::string format ( const std::string& fmt    ,
                       double             value  ) ;
  /// format single number ..
  std::string format ( const std::string& fmt    ,
                       long               value  ) ;
  /// format two numbers
  std::string format ( const std::string& fmt    , 
                       double             value1 ,
                       double             value2 ) ;
  /// format three numbers
  std::string format ( const std::string& fmt    , 
                       double             value1 , 
                       double             value2 ,
                       double             value3 ) ;
  /// format four numbers
  std::string format ( const std::string& fmt    ,
                       double             value1 , 
                       double             value2 ,
                       double             value3 ,
                       double             value4 ) ; 
  /// format five numbers
  std::string format ( const std::string& fmt    , 
                       double             value1 , 
                       double             value2 ,
                       double             value3 ,
                       double             value4 , 
                       double             value5 ) ; 
// ============================================================================
}
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_FORMAT_H
// ============================================================================
