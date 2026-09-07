// ============================================================================
#ifndef OSTAP_EPSILON_H 
#define OSTAP_EPSILON_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL 
// ============================================================================
#include <type_traits> // for std::decay_t
#include <limits>      // for std::numeric_limits
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    /** Use epsilon from numeric_limits as compile-time constant
     *  @code
     *  template <typename T>
     *  void my func ( T a ,
     *                 T e = Epsilon_v<T> ) ;
     *  @endcode
     */
    template <class T>
    constexpr std::decay_t<T> epsilon_v = []()
    { using CleanT = std::decay_t<T>;
      static_assert(std::numeric_limits<CleanT>::is_specialized, "Type has no numeric_limits!");
      return std::numeric_limits<CleanT>::epsilon(); }();
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_EPSILON_H
// ============================================================================
