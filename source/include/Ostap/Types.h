// ============================================================================
#ifndef OSTAP_TYPES_H
#define OSTAP_TYPES_H 1
// ============================================================================
// Inclide files  
// ============================================================================
// STD&STL
// ============================================================================
#include <cstdint>
#include <limits>
#include <string>
#include <vector>
#include <map> 
// ============================================================================
/** @file OStap/Typees.h
 *  Helper file with efinitini of few useful types  
 *  @author  Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2025-06-15
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  /** the tupe for Event  Index (for TTtee/RooAbsData looping
   * IT should be in agreement with TTree::kMaxEntries 
   * @see TTree::kMaxEntries
   * @see TVirtualTreePlayer::kMaxEntries 
   */   
  using EventIndex = unsigned long ;
  // ==========================================================================
  /** @var FirstEvent
   *     Index for the first event in the loop
   */
  const EventIndex FirstEvent { 0 } ;
  // ==========================================================================
  /** @var LastEvent
   *  Index for the last (exclusive) event in the loop
   */   
  const EventIndex LastEvent { std::numeric_limits<EventIndex>::max() } ; 
  // ==========================================================================
  /// the data type for ranges 
  typedef double DataType ;
  // ==========================================================================
  /** @var MinValue 
   *  minimal value for varioys ranges 
   */
  const DataType MinValue { -std::numeric_limits<DataType>::max () } ; 
  // ==========================================================================
  /** @var MaxValue 
   *  minimal value for varioys ranges 
   */
  const DataType MaxValue {  std::numeric_limits<DataType>::max () } ; 
  // ==========================================================================
  /// Types for keys 
  using Key  = std::string ; 
  /// Type fot names 
  using Name = Key ;
  /// Disctionanry type with string keys 
  template <typename Value>
  using Dict    = std::map<Key,Value> ;
  /// vector of strings/keys  
  using Strings = std::vector<Key>    ;
  /// vector of strings/keys  
  using Keys    =  std::vector<Key>   ; 
  /// Ditto 
  using Names   = std::vector<Name>   ;
  /// vector of doubles 
  using Doubles = std::vector<double> ;  
  // =========================================================================
}  //                                               The end of namespace Ostap
// ===========================================================================
#endif // OSTAP_TYPES_H  
// ===========================================================================
//                                                                     The END
// ===========================================================================


