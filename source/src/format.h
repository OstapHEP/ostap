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
  std::string format
  ( const std::string&  fmt    ,
    const double        value  ) ;
  /// format single number ..
  std::string format
  ( const std::string&  fmt    ,
    const long          value  ) ;
  /// format single number ..
  std::string format
  ( const std::string&  fmt    ,
    const unsigned long value  ) ;
  /// format two numbers
  std::string format
  ( const std::string&  fmt    , 
    const double        value1 ,
    const double        value2 ) ;
  /// format three numbers
  std::string format
  ( const std::string&  fmt    , 
    const double        value1 , 
    const double        value2 ,
    const double        value3 ) ;
  /// format four numbers
  std::string format
  ( const std::string&  fmt    ,
    const double        value1 , 
    const double        value2 ,
    const double        value3 ,
    const double        value4 ) ; 
  /// format five numbers
  std::string format
  ( const std::string&  fmt    , 
    const double        value1 , 
    const double        value2 ,
    const double        value3 ,
    const double        value4 , 
    const double        value5 ) ; 
  // ============================================================================
}
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_FORMAT_H
// ============================================================================
