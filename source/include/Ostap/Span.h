// ============================================================================
#ifndef OSTAP_SPAN_H 
#define OSTAP_SPAN_H 1
// ============================================================================
//
#if defined(__has_include)
#if __has_include(<version>)
#include <version>
#endif
#endif
// 
#if defined(__cpp_lib_span) && __cpp_lib_span >= 202002L
#include <span>
#define OSTAP_HAS_STD_SPAN 1
#else
#define OSTAP_HAS_STD_SPAN 0
#endif
//
// ============================================================================
#endif // OSTAP_SPAN_H 
// ============================================================================
//                                                                      The END
// ============================================================================
