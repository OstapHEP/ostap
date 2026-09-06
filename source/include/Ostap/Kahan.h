// ============================================================================
#ifndef OSTAP_KAHAN_H 
#define OSTAP_KAHAN_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL 
// ============================================================================
#include <cstdef>
#include <cstdint>
#include <iterator>
#include <type_traits>
#include <vector>
#include <array>
// ============================================================================
/** @file Ostap/Kahan.h
 *  Kahan's summation
 */
// ===========================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    /** Kahan summation 
     *  @see https://en.wikipedia.org/wiki/Kahan_summation_algorithm
     *  \f$ r = \sum_i x_i \f$ 
     *  @code
     *  // pseudocode 
     *  function KahanSum(input)
     *    var sum = 0.0
     *    var c = 0.0                 
     *    // A running compensation for lost low-order bits.
     *    for i = 1 to input.length do
     *       var y = input[i] - c     
     *    // So far, so good: c is zero.
     *       var t = sum + y          
     *    // Alas, sum is big, y small, so low-order digits of y are lost.
     *       c = (t - sum) - y        
     *    // (t - sum) cancels the high-order part of y; subtracting y recovers negative (low part of y)
     *       sum = t                  
     *    // Algebraically, c should always be zero. Beware overly-aggressive optimizing compilers!
     *    next i                      
     *    // Next time around, the lost low part will be added to y in a fresh attempt.
     *  return sum
     *  @endcode 
     *  @param begin (INPUT) begin-iterator for the input data 
     *  @param end   (INPUT) end-iterator for the input data 
     */
    template <class ITERATOR>
    inline auto sum_kahan
    ( ITERATOR begin , 
      ITERATOR end   )
    {
      // Automatically choose return type: preserve long double precision if requested, default to double otherwise
      using VALUE      = typename std::iterator_traits<ITERATOR>::value_type  ;
      using ReturnType = std::conditional_t<std::is_same_v<std::decay_t<VALUE>, long double>, long double, double>;
      
      long double sum = 0 ;
      long double c   = 0 ;
      for ( ; begin != end ; ++begin ) 
      {
        volatile const long double y = (*begin) - c ;
        volatile const long double t = sum      + y ;
        c        = ( t - sum ) - y ;
        sum      =   t             ;
      }
      return static_cast<RetrunType> ( sum ) ;
    }    
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
#endif // OSTAP_KAHAN_H
// ============================================================================
//                                                                      The END 
// ============================================================================
