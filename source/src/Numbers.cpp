// ============================================================================
// Include files
// ============================================================================
// STD&STL 
// ============================================================================
#include <cstdint>
#include <limits>
#include <array>
// ============================================================================
/// local namespace to hide all techinical symbols
namespace 
{  
  // ==========================================================================
  // Limits? 
  // ==========================================================================
  static_assert ( std::numeric_limits<float> ::is_specialized                ,
                  "std::numeric_limits<float>       is not specialized"      ) ;
  static_assert ( std::numeric_limits<double>::is_specialized                ,
                  "std::numeric_limits<double>      is not specialized"      ) ;
  static_assert ( std::numeric_limits<long double>::is_specialized           ,
                  "std::numeric_limits<long double> is not specialized"      ) ;
  // ==========================================================================  
  static_assert ( std::numeric_limits<float>::has_infinity                   ,
                  "std::numeric_limits<float> does not have infinity"        ) ;
  static_assert ( std::numeric_limits<double>::has_infinity                  ,
                  "std::numeric_limits<double> does not have infinity"       ) ;
  static_assert ( std::numeric_limits<long double>::has_infinity             ,
                  "std::numeric_limits<long double> does not have infinity"  ) ;
  // ==========================================================================    
  static_assert ( std::numeric_limits<float>::has_quiet_NaN                  ,
                  "std::numeric_limits<float> does not have quiet NaN"       ) ;
  static_assert ( std::numeric_limits<double>::has_quiet_NaN                 ,
                  "std::numeric_limits<double> does not have quiet NaN"      ) ;
  static_assert ( std::numeric_limits<long double>::has_quiet_NaN            ,
                  "std::numeric_limits<long double> does not have quiet NaN" ) ;
  // ==========================================================================
  
  // ==========================================================================
  static_assert ( std::numeric_limits<unsigned short> ::is_specialized         ,
                  "std::numeric_limits<unsigned short>  is not specialized"    ) ;
  static_assert ( std::numeric_limits<unsigned int>   ::is_specialized         ,
                  "std::numeric_limits<unsigned int>  is not specialized"      ) ;
  static_assert ( std::numeric_limits<unsigned long>  ::is_specialized         ,
                  "std::numeric_limits<unsigned long> is not specialized"      ) ;
  static_assert ( std::numeric_limits<unsigned long long>  ::is_specialized    ,
                  "std::numeric_limits<unsigned long long> is not specialized" ) ;
  // ==========================================================================
  
  // ==========================================================================
  static_assert ( std::numeric_limits<signed short> ::is_specialized         ,
                  "std::numeric_limits<signed short>  is not specialized"    ) ;
  static_assert ( std::numeric_limits<signed int>   ::is_specialized         ,
                  "std::numeric_limits<signed int>  is not specialized"      ) ;
  static_assert ( std::numeric_limits<signed long>  ::is_specialized         ,
                  "std::numeric_limits<signed long> is not specialized"      ) ;
  static_assert ( std::numeric_limits<signed long long>  ::is_specialized    ,
                  "std::numeric_limits<signed long long> is not specialized" ) ;
  // ==========================================================================
  
  // ==========================================================================
  static_assert ( std::numeric_limits<unsigned char> ::is_specialized        ,
                  "std::numeric_limits<unsigned char> is not specialized"    ) ;
  static_assert ( std::numeric_limits<signed char>    ::is_specialized         ,
                  "std::numeric_limits<signed char>   is not specialized"    ) ;
  static_assert ( std::numeric_limits<char>           ::is_specialized         ,
                  "std::numeric_limits<char>          is not specialized"    ) ;
  // ==========================================================================
  
  // ==========================================================================
  static_assert ( std::numeric_limits<std::int8_t>  ::is_specialized  && 
                  std::numeric_limits<std::int8_t>  ::is_integer      &&   
                  std::numeric_limits<std::int8_t>  ::is_signed       ,
                  "std::numeric_limits<std::int8_t> is not OK"        ) ;
  // ==========================================================================
  static_assert ( std::numeric_limits<std::int16_t>  ::is_specialized  && 
                  std::numeric_limits<std::int16_t>  ::is_integer      &&   
                  std::numeric_limits<std::int16_t>  ::is_signed       ,
                  "std::numeric_limits<std::int16_t> is not OK"        ) ;
  // ==========================================================================
  static_assert ( std::numeric_limits<std::int32_t>  ::is_specialized  && 
                  std::numeric_limits<std::int32_t>  ::is_integer      &&   
                  std::numeric_limits<std::int32_t>  ::is_signed       ,
                  "std::numeric_limits<std::int32_t> is not OK"        ) ;
  // ==========================================================================
  static_assert ( std::numeric_limits<std::int64_t>  ::is_specialized  && 
                  std::numeric_limits<std::int64_t>  ::is_integer      &&   
                  std::numeric_limits<std::int64_t>  ::is_signed       ,
                  "std::numeric_limits<std::int64_t> is not OK"        ) ;
  // ==========================================================================
  static_assert ( std::numeric_limits<std::intmax_t>  ::is_specialized  && 
                  std::numeric_limits<std::intmax_t>  ::is_integer      &&   
                  std::numeric_limits<std::intmax_t>  ::is_signed       ,
                  "std::numeric_limits<std::intmax_t> is not OK"        ) ;
  // ==========================================================================
  
  // ==========================================================================
  static_assert ( std::numeric_limits<std::uint8_t>  ::is_specialized  && 
                  std::numeric_limits<std::uint8_t>  ::is_integer      &&   
                  !std::numeric_limits<std::uint8_t> ::is_signed       ,
                  "std::numeric_limits<std::uint8_t> is not OK"        ) ;
  // ==========================================================================
  static_assert ( std::numeric_limits<std::uint16_t>  ::is_specialized  && 
                  std::numeric_limits<std::uint16_t>  ::is_integer      &&   
                  !std::numeric_limits<std::uint16_t> ::is_signed       ,
                  "std::numeric_limits<std::uint16_t> is not OK"        ) ;
  // ==========================================================================
  static_assert ( std::numeric_limits<std::uint32_t>  ::is_specialized  && 
                  std::numeric_limits<std::uint32_t>  ::is_integer      &&   
                  !std::numeric_limits<std::uint32_t> ::is_signed       ,
                  "std::numeric_limits<std::uint32_t> is not OK"        ) ;
  // ==========================================================================
  static_assert ( std::numeric_limits<std::uint64_t>  ::is_specialized  && 
                  std::numeric_limits<std::uint64_t>  ::is_integer      &&   
                  !std::numeric_limits<std::uint64_t> ::is_signed       ,
                  "std::numeric_limits<std::uint64_t> is not OK"        ) ;
  // ==========================================================================
  static_assert ( std::numeric_limits<std::uintmax_t>  ::is_specialized  && 
                  std::numeric_limits<std::uintmax_t>  ::is_integer      &&   
                  !std::numeric_limits<std::uintmax_t> ::is_signed       ,
                  "std::numeric_limits<std::uintmax_t> is not OK"        ) ;
  // ==========================================================================
  
  // ==========================================================================
  static_assert ( std::numeric_limits<char>  ::is_specialized  && 
                  std::numeric_limits<char>  ::is_integer      &&   
                  "std::numeric_limits<char> is not OK"        ) ;
  // ==========================================================================
  static_assert ( std::numeric_limits<signed char>  ::is_specialized  && 
                  std::numeric_limits<signed char>  ::is_integer      &&   
                  std::numeric_limits<signed char>  ::is_signed       ,
                  "std::numeric_limits<signed char> is not OK"        ) ;
  // ==========================================================================
  static_assert ( std::numeric_limits<unsigned char>  ::is_specialized  && 
                  std::numeric_limits<unsigned char>  ::is_integer      &&   
                  !std::numeric_limits<unsigned char> ::is_signed       ,
                  "std::numeric_limits<unsigned char> is not OK"        ) ;
  // ==========================================================================

  // ==========================================================================
  static_assert ( std::numeric_limits<std::size_t>  ::is_specialized  && 
                  std::numeric_limits<std::size_t>  ::is_integer      &&   
                  !std::numeric_limits<std::size_t> ::is_signed       ,
                  "std::numeric_limits<std::size_t> is not OK"        ) ;
  // ==========================================================================
  
  // ==========================================================================
  static_assert ( std::numeric_limits<std::intptr_t>  ::is_specialized  && 
                  std::numeric_limits<std::intptr_t>  ::is_integer      &&   
                  std::numeric_limits<std::intptr_t> ::is_signed        ,
                  "std::numeric_limits<std::intptr_t> is not OK"        ) ;
  // ==========================================================================

  // ==========================================================================
  static_assert ( std::numeric_limits<std::uintptr_t>  ::is_specialized  && 
                  std::numeric_limits<std::uintptr_t>  ::is_integer      &&   
                  !std::numeric_limits<std::uintptr_t> ::is_signed       ,
                  "std::numeric_limits<std::uintptr_t> is not OK"        ) ;
  // ==========================================================================
  
  
  
}


// ============================================================================
//                                                                      The END 
// ============================================================================
