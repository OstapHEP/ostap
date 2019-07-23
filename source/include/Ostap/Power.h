// ============================================================================
#ifndef OSTAP_POWER_H 
#define OSTAP_POWER_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL 
// ============================================================================
#include <cmath>
// ============================================================================
/** @file Ostap/Power.h
 *  This file is taken from the LoKi project - 
 *    "C++ ToolKit  for Smart and Friendly Physics Analysis"
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2001-01-23 
 */
// ============================================================================
namespace Ostap
{ 
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    // using namespace std ; ?
    // ========================================================================
#ifdef __INTEL_COMPILER         // Disable ICC remark
#pragma warning(disable:2259) //  non-pointer conversion may lose significant bits
#pragma warning(push)
#endif
    // ========================================================================
    /** Simple utility for efficient "pow".
     *  It works only for positive integer powers.
     *
     *  @code 
     *   
     *   const double result = Ostap::Math::pow ( value , 10 ) ;
     *
     *  @endcode 
     *
     *  The actual code is copied from 
     *     std::__cmath_power bits/cmath.tcc 
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@lapp.in2p3.fr
     *  @date 2005-04-09 
     */
    template<typename TYPE>
    inline TYPE pow ( TYPE __x , unsigned long __n )
    {
      //
      TYPE __y = __n % 2 ? __x : 1;
      //
      while ( __n >>= 1 )
      {
        __x = __x * __x;
        if ( __n % 2) { __y = __y * __x; }
      } 
      //
      return __y ;
    }
    // ========================================================================
#ifdef __INTEL_COMPILER         // Disable ICC remark
#pragma warning(push)
#endif
    // ========================================================================
  } //                                                    end of namespace Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
// The END 
// ============================================================================
#endif // OSTAP_POWER_H
// ============================================================================
