// Include files 
// ============================================================================
//  STD& STL
// ============================================================================
#include <cstdio>
// ============================================================================
// local
// ============================================================================
#include "format.h"
// ============================================================================
/** @file
 *  very ugly way to avoid boost::format 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// ============================================================================
// format single number ..
// ============================================================================
std::string Ostap::format 
( const std::string& fmt    ,
  const double       value1 ) 
{
  const  std::size_t s_LEN = 1024 ;
  static char        s_buffer[ s_LEN ] ;
  const int result = snprintf ( s_buffer    , 
                                s_LEN       , 
                                fmt.c_str() , 
                                value1      ) ;
  return 
    0 <= result && (unsigned int) result < s_LEN ?
    std::string ( s_buffer , s_buffer + result ) : 
    fmt + std::to_string ( value1 ) ;
}
// ============================================================================
// format single number ..
// ============================================================================
std::string Ostap::format 
( const std::string& fmt    ,
  const long         value1 ) 
{
  const  std::size_t  s_LEN  = 1024 ;
  static char         s_buffer[ s_LEN ] ;
  const  int          result = snprintf ( s_buffer    , 
                                          s_LEN       ,
                                          fmt.c_str() , 
                                          value1      ) ;
  return 
    0 <= result && (unsigned int) result < s_LEN ?
    std::string ( s_buffer , s_buffer + result ) : 
    fmt + std::to_string ( value1 ) ;
}
// ============================================================================
// format single number ..
// ============================================================================
std::string Ostap::format 
( const std::string&  fmt    ,
  const unsigned long value1 ) 
{
  const  std::size_t  s_LEN  = 1024 ;
  static char         s_buffer[ s_LEN ] ;
  const  int          result = snprintf ( s_buffer    , 
                                          s_LEN       ,
                                          fmt.c_str() , 
                                          value1      ) ;
  return 
    0 <= result && (unsigned int) result < s_LEN ?
    std::string ( s_buffer , s_buffer + result ) : 
    fmt + std::to_string ( value1 ) ;
}
// ============================================================================
// format two numbers ..
// ============================================================================
std::string Ostap::format 
( const std::string& fmt    ,
  const double       value1 , 
  const double       value2 ) 
{
  const  std::size_t s_LEN  = 1024 ;
  static char        s_buffer [ s_LEN ] ;
  const  int         result = snprintf ( s_buffer    , 
                                         s_LEN       , 
                                         fmt.c_str() , 
                                         value1      , 
                                         value2      ) ;
  return 
    0 < result && (unsigned int) result < s_LEN ?
    std::string ( s_buffer , s_buffer + result ) : 
    fmt + std::to_string ( value1 ) + " " + std::to_string ( value2 ) ;
}
// ============================================================================
// format three numbers ..
// ============================================================================
std::string Ostap::format
( const std::string& fmt    ,
  const double       value1 , 
  const double       value2 ,
  const double       value3 ) 
{
  const  std::size_t s_LEN  = 1024 ;
  static char        s_buffer [ s_LEN ] ;
  const  int         result = snprintf ( s_buffer    , 
                                         s_LEN       , 
                                         fmt.c_str() , 
                                         value1      , 
                                         value2      , 
                                         value3      ) ;
  return 
    0 < result && (unsigned int) result < s_LEN ?
    std::string ( s_buffer , s_buffer + result ) : 
    ( fmt + std::to_string ( value1 ) + 
      " " + std::to_string ( value2 ) +
      " " + std::to_string ( value3 ) ) ;
}
// ============================================================================
// format four numbers ..
// ============================================================================
std::string Ostap::format
( const std::string& fmt    ,
  const double       value1 , 
  const double       value2 ,
  const double       value3 ,
  const double       value4 ) 
{
  const  std::size_t s_LEN = 1024 ;
  static char        s_buffer[ s_LEN ] ;
  const  int         result = snprintf ( s_buffer    , 
                                         s_LEN       , 
                                         fmt.c_str() , 
                                         value1      , 
                                         value2      , 
                                         value3      , 
                                         value4      ) ;
  //
  return 
    0 < result && (unsigned int) result < s_LEN ?
    std::string ( s_buffer , s_buffer + result ) : 
    ( fmt + std::to_string ( value1 ) + 
      " " + std::to_string ( value2 ) +
      " " + std::to_string ( value3 ) + 
      " " + std::to_string ( value4 ) ) ;
}
// ============================================================================
// format five numbers ..
// ============================================================================
std::string Ostap::format 
( const std::string& fmt    ,
  const double       value1 , 
  const double       value2 ,
  const double       value3 ,
  const double       value4 ,
  const double       value5 ) 
{
  const  std::size_t s_LEN  = 1024 ;
  static char        s_buffer[ s_LEN ] ;
  const  int         result = snprintf ( s_buffer    , 
                                         s_LEN       , 
                                         fmt.c_str() , 
                                         value1      , 
                                         value2      , 
                                         value3      , 
                                         value4      , 
                                         value5      ) ;
  //
  return 
    0 < result && ( unsigned int) result < s_LEN ?
    std::string ( s_buffer , s_buffer + result ) :
    ( fmt + std::to_string ( value1 ) + 
      " " + std::to_string ( value2 ) +
      " " + std::to_string ( value3 ) + 
      " " + std::to_string ( value4 ) + 
      " " + std::to_string ( value5 ) ) ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================
