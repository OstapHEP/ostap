// ============================================================================
#ifndef OSTAP_TYPES_H
#define OSTAP_TYPES_H 1
// ============================================================================
// Include files  
// ============================================================================
// STD&STL
// ============================================================================
#if defined(__has_include)
#if __has_include(<version>)
#include <version>
#endif
#endif
// ============================================================================
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <limits>
#include <string>
#include <vector>
#include <map> 
// ============================================================================
#if defined(__cpp_lib_byte) && __cpp_lib_byte >= 201603L
#include <cstddef>
#define OSTAP_HAS_STD_BYTE 1
#else
#define OSTAP_HAS_STD_BYTE 0
#endif
// ============================================================================
/** @file Ostap/Types.h
 *  Helper file with definition of few useful types  
 *  @author  Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2025-06-15
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  /** the type for Event  Index (for TTree/RooAbsData looping
   * It should be in agreement with TTree::kMaxEntries 
   * @see TTree::kMaxEntries
   * @see TVirtualTreePlayer::kMaxEntries 
   */   
  using EventIndex = unsigned long ;
  // ==========================================================================
  /** @var FirstEvent
   *     Index for the first event in the loop
   */
  constexpr EventIndex FirstEvent { 0 } ;
  // ==========================================================================
  /** @var LastEvent
   *  Index for the last (exclusive) event in the loop
   */   
  constexpr EventIndex LastEvent { std::numeric_limits<EventIndex>::max() } ; 
  // ==========================================================================
  /// the data type for ranges 
  typedef double DataType ;
  // ========================================================================== 
  /** @var MinValue 
   *  minimal value for various ranges 
   */
  constexpr DataType MinValue { -std::numeric_limits<DataType>::max () } ; 
  // ==========================================================================
  /** @var MaxValue 
   *  maximal value for various ranges 
   */
  constexpr DataType MaxValue {  std::numeric_limits<DataType>::max () } ; 
  // ==========================================================================
  /// Types for keys 
  using Key        = std::string ; 
  /// Type fot names 
  using Name       = Key ;
  /// Dictionary type with string keys 
  template <typename Value>
  using Dict       = std::map<Key,Value>  ;
  /// the dictorinnaty 
  using Dictionary = Dict<std::string> ; 
  /// vector of strings/keys  
  using Strings    = std::vector<Key>     ;
  /// vector of strings/keys  
  using Keys       =  std::vector<Key>    ; 
  /// Ditto 
  using Names      = std::vector<Name>    ;
  /// vector of doubles 
  using Doubles    = std::vector<double>  ;  
  // =========================================================================
  
  // =========================================================================
  /// Numeric Type ? 
  template <typename T>
  struct is_numeric
  {
    // =======================================================================
  private:
    // =======================================================================
    using CleanT = typename std::remove_cv<T>::type;
    // =======================================================================
  public:
    // =======================================================================
    static constexpr bool value = std::is_arithmetic<CleanT>::value ;
    // =======================================================================
  };
  // =========================================================================  
  /// Numeric Type or Bytes (for buffers) 
  template <typename T>
  struct is_numeric_or_byte
  {
  private:
    // =======================================================================
    using CleanT = typename std::remove_cv<T>::type;
    // =======================================================================
  public:
    // =======================================================================
    static constexpr bool value = std::is_arithmetic<CleanT>::value
#if defined(OSTAP_HAS_STD_BYTE) && OSTAP_HAS_STD_BYTE      
      || std::is_same<CleanT, std::byte>::value
#endif 
      ;
    // =======================================================================
  }; 
  // =========================================================================
  /// convertibe to double ? iterators for numerical sequences 
  template <typename T>
  struct is_convertible_to_double 
  {
    // =======================================================================
  private:
    // =======================================================================
    using CleanT = typename std::remove_cv<T>::type;
    // =======================================================================
  public:
    // =======================================================================
    static constexpr bool value = std::is_convertible<CleanT,double>::value ;
    // =======================================================================
  };
  // =========================================================================

  // =========================================================================
  template <typename T>
  constexpr bool is_numeric_v               = is_numeric<T>::value;
  // =========================================================================
  template <typename T>
  constexpr bool is_numeric_or_byte_v       = is_numeric_or_byte<T>::value;
  // =========================================================================
  template <typename T>
  constexpr bool is_convertible_to_double_v = is_convertible_to_double<T>::value;
  // =========================================================================

  // =========================================================================
} //                                                The end of namespace Ostap
// ===========================================================================
#endif // OSTAP_TYPES_H  
// ===========================================================================
//                                                                     The END
// ===========================================================================


