// ============================================================================
#ifndef OSTAP_HASH_H
#define OSTAP_HASH_H 1
// ============================================================================
// STD&STL
// ============================================================================
#include <functional>
// ============================================================================
/** @file Ostap/Hash.h
 *  few minor utilities for hashing 
 *  - hash_combine 
 *  - hash_combiner
 *  - hash_range 
 *  @author Vanya Belyaev
 *  @date   2018-03-23
 */
namespace Ostap
{
  // ==========================================================================
  namespace Utils 
  {
    // ========================================================================
    ///  @see http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2018/p0814r2.pdf
    // ==========================================================================
    inline void hash_combine ( size_t& /* seed */ ) {}
    // ========================================================================
    template<typename T>
    void hash_combine ( std::size_t& seed, const T& val)
    { seed ^= std::hash<T>()(val) + 0x9e3779b9 + (seed<<6) + (seed>>2); }
    // ========================================================================
    /// hash for sequence of objects 
    template <typename T, typename... Types>
    void hash_combine ( std::size_t& seed , const T& val, const Types&... args)
    {
      hash_combine ( seed , val     ) ;
      hash_combine ( seed , args... ) ;
    }
    // ========================================================================
    /// combine hash for sequence of objects 
    template<typename RT = std::size_t, typename... T>
    RT hash_combiner ( const T&... args )
    {
      std::size_t seed = 0;
      hash_combine ( seed, args...) ;
      return seed ;
    }
    // ========================================================================
    /// hash the range 
    template <class   ITERATOR>
    std::size_t hash_range ( ITERATOR i1 , ITERATOR i2 )
    {
      std::size_t seed = 0 ;
      for ( ; i1 != i2 ; ++i1 ) { hash_combine ( seed , *i1 ) ; }
      return seed ;
    }
    // ========================================================================
    /// hash the range 
    template <class CONTAINER >
    std::size_t hash_range ( const CONTAINER& cnt )
    { return hash_range ( cnt.begin () , cnt.end() ) ; }
    // ========================================================================
  } //                                        The end of namespace Ostap::Utils 
  // ==========================================================================
} //                                                The  end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_HASH_H
// ============================================================================
