// ============================================================================
#ifndef OSTAP_LOCAL_HASH_H 
#define OSTAP_LOCAL_HASH_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <functional>
// ============================================================================
namespace std
{
  // ==========================================================================
  ///  @see http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2018/p0814r2.pdf
  // ==========================================================================
  template<typename T>
  void _hash_combine ( size_t& seed, const T& val)
  {
    seed ^= std::hash<T>()(val) + 0x9e3779b9 + (seed<<6) + (seed>>2);
  }
  // ============================================================================
  // template<typename RT = size_t, typename... T>
  // RT hash_combine ( const T&... args )
  // {
  //   size_t seed = 0;
  //   (_hash_combine(seed,args) , ... ); // create hash value with seed over all args
  //   return seed;
  // }
  // ===========================================================================
  inline void _hash_combine ( size_t& seed ) {}
  template <typename T, typename... Types>
  void _hash_combine ( std::size_t& seed , const T& val, const Types&... args)
  {
    _hash_combine(seed,val);
    _hash_combine(seed,args...);
  }
  template<typename RT = size_t, typename... T>
  RT hash_combine ( const T&... args )
  {
    size_t seed = 0;
    _hash_combine ( seed, args...) ;
    return seed ;
  }
  // ==========================================================================
  template <typename RT=size_t , typename T1, typename ...T>
  RT hash_combine  ( const T&...   args , const T1& t1 ) 
  {
    size_t seed = 0 ;
    _hash_combine ( seed , t1 ) ;
    _hash_combine ( seed , _hash_combine ( args... ) ) ;
    return seed ; 
  }
  // ==========================================================================
  template <>
  struct hash<const char*> 
  {
    size_t operator() ( const char *str ) const 
    {
      std::size_t seed =  0 ;
      char c;
      while ( c = *str++) { _hash_combine ( seed , c ) ; }
      return seed ;
    }
  } ;
  // ==========================================================================
} //
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_LOCAL_HASH_H
// ============================================================================
