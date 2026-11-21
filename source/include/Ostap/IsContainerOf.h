// ============================================================================
#ifndef OSTAP_ISCONTAINEROF_H
#define OSTAP_ISCONTAINEROF_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <iterator>     // std::iterator_traits
#include <type_traits>  // std::enable_if, std::enable_if_t, std::is_convertible, std::is_convertible_v
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  /** @struct is_container_of
   *  Type trait checking whether a range/container @p CONTAINER supports
   *  iteration via ADL/std::begin and std::end and has an element value type
   *  implicitly convertible to @p T.
   *
   *  @tparam CONTAINER Candidate container or range type.
   *  @tparam T         Target element type to convert each element to.
   */
  template < class CONTAINER , class T , class = void >
  struct is_container_of : std::false_type {} ;
  
  // ==========================================================================
  /** Specialization for types exposing valid begin/end and iterator_traits.
   *
   *  Uses std::void_t with SFINAE so non-iterable types or types without a
   *  valid iterator value_type cleanly fall back to std::false_type rather
   *  than triggering hard compile errors during overload resolution.
   */
  template < class CONTAINER , class T >
  struct is_container_of
    < CONTAINER ,
      T ,
      std::void_t
      < decltype ( std::begin ( std::declval < const CONTAINER& > () ) ) ,
        decltype ( std::end   ( std::declval < const CONTAINER& > () ) ) ,
        typename std::iterator_traits
        < decltype ( std::begin ( std::declval < const CONTAINER& > () ) ) >::value_type
      > >
    : std::bool_constant
      < std::is_convertible_v
        < typename std::iterator_traits
          < decltype ( std::begin ( std::declval < const CONTAINER& > () ) ) >::value_type ,
          T
        > > {} ;

  // ==========================================================================
  /** Convenience variable template for Ostap::is_container_of. */
  template < class CONTAINER , class T >
    inline constexpr bool is_container_of_v =
    is_container_of < CONTAINER , T > :: value ;
  // ==========================================================================
} //                                                     End of namespace Ostap
// ============================================================================
#endif // OSTAP_ISCONTAINEROF_H
// ============================================================================
//                                                                      The END
// ============================================================================
