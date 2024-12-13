// ============================================================================
#ifndef OSTAP_LOCAL_HASH_H 
#define OSTAP_LOCAL_HASH_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <functional>
#include <complex>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Utils.h"
// ============================================================================
namespace std
{
  // ==========================================================================
  template <>
  struct hash<const char*> 
  {
    inline size_t operator() ( const char *str ) const 
    {
      using Ostap::Utils::hash_combine  ;
      std::size_t seed =  0 ;
      char c;
      while ( ( c = *str++ ) ) { hash_combine ( seed , c ) ; }
      return seed ;
    }
  } ;
  // ==========================================================================
  template <typename T, int N>
  struct hash<T(&)[N]> 
  {
    inline size_t operator() ( const T(&s)[N] ) const 
    { 
      using Ostap::Utils::hash_range ;
      return hash_range ( &s , &s+N ); 
    }
  } ;
  // ==========================================================================
  template <typename T, int N>
  struct hash<const T(&)[N]> 
  {
    inline size_t operator() ( const T(&s)[N] ) const 
    { 
      using Ostap::Utils::hash_range ;
      return hash_range ( &s , &s+N ); 
    }
  } ;
  // ==========================================================================
  template <class T, class A>
  struct hash<std::vector<T,A> >
  {
    size_t operator() ( const std::vector<T,A>& v ) const 
    { 
      using Ostap::Utils::hash_range ;
      return hash_range ( v.begin () , v.end () ) ;
    }
  } ;
  // ===========================================================================
  template<class T>
  struct std::hash< std::complex<T> >
  {
    std::size_t operator()(const std::complex<T>& s ) const noexcept
    {
      const std::size_t hr = std::hash<T>{} ( std::real ( s ) ) ;
      const std::size_t hi = std::hash<T>{} ( std::imag ( s ) ) ;
      using Ostap::Utils::hash_combiner ;
      return hash_combiner ( hr , hi ) ;
    }
    // =======================================================================
  }; 
  // ==========================================================================
} //
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_LOCAL_HASH_H
// ============================================================================
