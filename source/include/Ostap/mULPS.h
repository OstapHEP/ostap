// ============================================================================
#ifndef OSTAP_MULPS_H 
#define OSTAP_MULPS_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL 
// ============================================================================
#include <type_traits> // for std::decay_t
#include <complex>     
#include <array>     
#include <vector>     
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    /** @var mULPS_float
     *  "tolerance" parameter for "Lomont"-compare of floating point numbers.
     *  It corresponds to relative ("Knuth/GLS") tolerance of about ~6*10^-6
     *  for values in excess of 10^-37.
     *
     *  @see Ostap::Math::Lomont::compare_float 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2010-01-02
     */
    constexpr unsigned short mULPS_float = 100 ;
    // ========================================================================
    /** @var mULPS_float_low
     *  "Low-tolerance" parameter for "Lomont"-compare of floating point numbers.
     *  It corresponds to relative ("Knuth/GLS") tolerance of about ~6*10^-5
     *  for values in excess of 10^-37.
     *
     *  @see Ostap::Math::Lomont::compare_float 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2010-01-02
     */
    constexpr unsigned short mULPS_float_low = 1000 ;
    // =========================================================================
    /** @var mULPS_double
     *  "tolerance" parameter for "Lomont"-compare of floating point numbers.
     *  It corresponds to relative ("Knuth/GLS") tolerance of about ~6*10^-13
     *  for values in excess of 10^-304.
     *  @see Ostap::Math::Lomont::compare_double
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2010-01-02
     */
    constexpr unsigned int mULPS_double = 1000 ;
    // ========================================================================
    namespace detail
    {
      // ======================================================================
      /// Default fallback: return mULPS_double for any generic type
      template <class T>
      struct ulps_deducer
      {
        static constexpr auto value = mULPS_double ;
      };      
      // ======================================================================
      /// Specialization for float types
      template <>
      struct ulps_deducer<float>
      {
        static constexpr auto value = mULPS_float;
      };      
      // ======================================================================
      /// Specialization for complex types: delegate to the underlying scalar type T
      template <class T>
      struct ulps_deducer<std::complex<T>>
      {
        static constexpr auto value = ulps_deducer<std::decay_t<T>>::value;
      } ;
      // ======================================================================
      /// Specialization for std::vector: delegate to its element type
      template <class TYPE, class ALLOCATOR>
      struct ulps_deducer<std::vector<TYPE, ALLOCATOR>>
      {
        static constexpr auto value = ulps_deducer<std::decay_t<TYPE>>::value;
      };
      // ======================================================================      
      /// Specialization for std::array: delegate to its element type
      template <class TYPE, std::size_t N>
      struct ulps_deducer<std::array<TYPE, N>>
      {
        static constexpr auto value = ulps_deducer<std::decay_t<TYPE>>::value;
      };
      // =====================================================================
    } 
    // =======================================================================
    /// Compile-time selection of mULPs supporting raw scalars and complex types
    template <class T>
    inline constexpr auto mULPS = detail::ulps_deducer<std::decay_t<T>>::value;
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_MULPS_H
// ============================================================================
