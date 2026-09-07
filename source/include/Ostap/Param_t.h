// ============================================================================
#ifndef OSTAP_PARAM_T_H 
#define OSTAP_PARAM_T_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL 
// ============================================================================
#include <type_traits> 
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    namespace detail
    {
      // ======================================================================
      /**
       * @tparam T The input type to analyze and optimize for parameter passing.
       * @brief Metafunction to determine the optimal parameter passing type.
       * @details Automatically evaluates whether a type is small enough to be passed 
       *          by value or requires a const reference, while also enforcing 
       *          pointer-to-const safety for raw pointers.
       */
      template <class T>
      struct param_traits
      {
        /// Decay the type to remove top-level references, const-qualifiers, 
        /// and to convert arrays/functions to pointers for uniform analysis.
        using clean_type = std::decay_t<T>;
        
        /// Determines if the type is "small" and cheap to pass by value.
        /// Evaluates to true if the type is trivially copyable and its size 
        /// does not exceed twice the machine word size (typically up to 16 bytes on 64-bit).
        static constexpr bool is_small = 
          std::is_trivially_copyable_v<clean_type> && 
          (sizeof(clean_type) <= 2 * sizeof(void*));
        
        /// Normalizes raw pointers into pointers-to-const (e.g., T* becomes const T*) 
        /// to prevent accidental modification of the pointee inside the function.
        /// Preserves other types as-is.
        using base_type = std::conditional_t<
          std::is_pointer_v<clean_type>,
          std::add_pointer_t<std::add_const_t<std::remove_pointer_t<clean_type>>>,
          clean_type
          >;
        
        /// The final optimized parameter type:
        /// - Passed by value (base_type) if it is small.
        /// - Passed by reference-to-const (const base_type&) if it is large.
        using type = std::conditional_t<is_small, base_type, const base_type&>;
      };
      // ======================================================================
    }
    // ========================================================================
    /**
     * @tparam T The input type.
     * @brief Convenient alias template for cleaner syntax in function signatures.
     */
    template <class T>
    using param_t = typename detail::param_traits<T>::type;
    // ======================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_PARAM_T_H
// ============================================================================
