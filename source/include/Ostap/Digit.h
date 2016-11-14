// ============================================================================
#ifndef OSTAP_DIGIT_H 
#define OSTAP_DIGIT_H 1
// ============================================================================
// Include files
// ============================================================================
//  STD& STL 
// ============================================================================
#include <limits>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/TypeWrapper.h"
#include "Ostap/IPower.hpp"
#include "Ostap/Power.h"
// ============================================================================
/** @file
 *
 *  The collection of useful utilities for evaluation of 
 *  decimal digits for the unsigned integral values 
 *
 *  The utilities for evaluation of 
 *  Nth decimal digit for the unsigned integral values:
 *
 *   -  Ostap::Math::IDigit 
 *      The most efficient compile-time evaluator.  
 *      It is applicable if both N and value are compile-time constanst 
 *   -  Ostap::Math::Digit 
 *      Rather efficient fuctor, N is compile-time constanst 
 *   -  Ostap::Math::digit
 *      The regular function, the least efficient evaluator 
 *  
 *  The utilities for evaluation of 
 *  decimal N1->N2 decimal digits for the unsigned integral values:
 *
 *   -  Ostap::Math::IDigits 
 *      The most efficient compile-time evaluator.  
 *      It is applicable if both N1/N2 and value are compile-time constanst 
 *   -  Ostap::Math::Digits 
 *      Rather efficient fuctor, N1&N2 are compile-time constanst 
 *   -  Ostap::Math::digits
 *      The regular function, the least efficient evaluator 
 *
 *  @attention the least significant decimal digit is numbered as #0
 *
 *  @see Ostap::Math::IDigit 
 *  @see Ostap::Math::Digit 
 *  @see Ostap::Math::digit
 *
 *  @see Ostap::Math::IDigits 
 *  @see Ostap::Math::Digits 
 *  @see Ostap::Math::digits
 *
 *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
 *  @date 2008-08-01
 */
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
      /** @struct Check10 
       *  Simple Helper structure to check if the type is able to 
       *  contain the certain number of decimal digits 
       *
       *   - "safe" 
       *   - "value"
       *
       *  @attention the least significat decimal digir is numbered as #0
       *
       *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
       *  @date 2008-08-01
       */
      template <class TYPE,unsigned int N>
      struct Check10 
      {
        // ====================================================================
        static_assert (  std::numeric_limits<TYPE>::is_specialized && 
                         std::numeric_limits<TYPE>::is_integer     && 
                         !std::numeric_limits<TYPE>::is_signed      , 
                         "Check10: inappropriate type" ) ;
        // ====================================================================
        enum {
          /// Nth digit is "save"
          safe  = N <  (unsigned int) std::numeric_limits<TYPE>::digits10 ,
          /// Nth digit is still OK
          value = N <= (unsigned int) std::numeric_limits<TYPE>::digits10 
          } ;
        // ====================================================================
      };
      // ======================================================================
      /** @struct _IDigit
       *  Simple structure form compile-time evaluation 
       *  of the decimal 
       *  digit for the given number 
       *
       *  @attention the least significat decimal digir is numbered as #0
       *
       *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
       *  @date 2008-08-01
       */
      template <class TYPE,
                typename Ostap::Math::TypeWrapper<TYPE>::value_type I , 
                unsigned int N>
      struct _IDigit
      {
        // ====================================================================
        static_assert (  std::numeric_limits<TYPE>::is_specialized && 
                         std::numeric_limits<TYPE>::is_integer     && 
                         !std::numeric_limits<TYPE>::is_signed      , 
                         "_IDigit: inappropriate type" ) ;
        // ====================================================================
        static_assert ( N <= std::numeric_limits <TYPE>::digits10 , 
                        "_IDigit: invalid index" )  ;
        // ====================================================================
        enum 
          { 
            value = 
            Check10<TYPE,N>::safe ? 
            (I/Ostap::Math::IPower<TYPE              ,10,N>::value )%10 :
            (I/Ostap::Math::IPower<unsigned long long,10,N>::value )%10 
          } ;
        // ====================================================================
      } ;
    // ======================================================================
      /** @struct _IDigits 
       *  Helper structure for compile-time evaluation of 
       *  range N1->N2 of decimal digits from the integral type 
       *
       *  @attention the least significat decimal digit is numbered as #0
       *
       *  @author Vanya BELYAEV Ivan.Belyaev@nilkhef.nl
       *  @date 2008-07-31
       */
      template <class TYPE,
                typename Ostap::Math::TypeWrapper<TYPE>::value_type I , 
                unsigned int N1 ,
                unsigned int N2 >  
      struct _IDigits
      {
        // ====================================================================
        static_assert (  std::numeric_limits<TYPE>::is_specialized && 
                         std::numeric_limits<TYPE>::is_integer     && 
                         !std::numeric_limits<TYPE>::is_signed      , 
                        "_IDigits: inappropriate type"              ) ;
        // ====================================================================
        static_assert ( N1 < N2  
                        && N1 <= std::numeric_limits<TYPE>::digits10  
                        && N2 <= std::numeric_limits<TYPE>::digits10 +1 , 
                        "_IDigits: invalid indices" ) ;
        // ====================================================================
        enum { 
          value = 
          Check10<TYPE,N1>::safe && Check10<TYPE,N2-N1>::safe     ?
          (I/Ostap::Math::IPower<TYPE              ,10,N1>::value )%
          Ostap::Math::IPower<TYPE,10,N2-N1>::value :
          (I/Ostap::Math::IPower<unsigned long long,10,N1>::value )%
          Ostap::Math::IPower<unsigned long long,10,N2-N1>::value 
        };
        // ====================================================================`
      } ;
      // ======================================================================
      /** @struct _Digit 
       *  Helper structure to get Nth decimal digits from the integral type 
       *
       *  @attention the least significant decimal digit is numbered as #0
       *
       *  @author Vanya BELYAEV Ivan.Belyaev@nilkhef.nl
       *  @date 2008-07-31
       */
      template <class TYPE, unsigned int N>
      struct _Digit : public std::unary_function<TYPE,int>
      {
        // ====================================================================
        static_assert (  std::numeric_limits<TYPE>::is_specialized && 
                         std::numeric_limits<TYPE>::is_integer     && 
                         !std::numeric_limits<TYPE>::is_signed      , 
                         "_Digit: inappropriate type"                ) ;                         
        // ====================================================================
        static_assert ( N <= std::numeric_limits<TYPE>::digits10   , 
                        "_Digit: invaild index"                      ) ;
        // ====================================================================
        enum { value = Ostap::Math::IPower<unsigned long long,10,N>::value } ;
        // ====================================================================
        inline int operator () ( const TYPE v ) const 
        { return ( v / value ) % 10 ; }
        // ====================================================================
      } ;
      // ======================================================================
      template <class TYPE, unsigned int N, bool OK>
      struct __Dig10 
      { enum { value = Ostap::Math::IPower<TYPE,10,N>::value } ; } ;
      // ======================================================================
      template <class TYPE, unsigned int N>
      struct __Dig10<TYPE,N,false> 
      { enum { value = Ostap::Math::IPower<unsigned long long,10,N>::value } ; } ;
      // ======================================================================
      /** @struct _Digits 
       *  Helper structure for evaluation of 
       *  range N1->N2 of decimal digits from the integral type 
       *
       *  @attention the least significat decimal digit is numbered as #0
       *
       *  @author Vanya BELYAEV Ivan.Belyaev@nilkhef.nl
       *  @date 2008-07-31
       */
      template <class TYPE, 
                unsigned int N1,
                unsigned int N2>
      struct _Digits : public std::unary_function<TYPE,TYPE>
      {
      private:
        // ====================================================================
        static_assert (  std::numeric_limits<TYPE>::is_specialized && 
                         std::numeric_limits<TYPE>::is_integer     && 
                         !std::numeric_limits<TYPE>::is_signed      , 
                         "_Digits: inappropriate type"                ) ;
        // ====================================================================
        static_assert ( N1 < N2  
                        && N1 <= std::numeric_limits<TYPE>::digits10     
                        && N2 <= std::numeric_limits<TYPE>::digits10 + 1 , 
                        "_Digits: invalid index"  ) ;
        // ====================================================================
        enum {
          val1 = __Dig10<TYPE,N1,
          Check10<TYPE,N1>::safe && Check10<TYPE,N2-N1>::safe>:: value ,
          val2 = __Dig10<TYPE,N2-N1,
          Check10<TYPE,N1>::safe && Check10<TYPE,N2-N1>::safe>:: value 
        } ;
        // ====================================================================        
      public:
        // ====================================================================
        /// the only on eessential method 
        inline TYPE operator() ( const TYPE v ) const 
        { return static_cast<TYPE> (  ( v / val1 ) % val2 ) ; }
        // ====================================================================
      } ;  
    }
    // ========================================================================
    /** @struct IDigit
     *  Simple structure form compile-time evaluation 
     *  of the Nth decimal digit for the given number 
     
     *  @code
     * 
     *  const int dig5 = IDigit<unsigned int,4362736,5>::value ;
     * 
     *  @endcode 
     *
     *  @attention the least significant decimal digit is numbered as #0
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@nilkhef.nl
     *  @date 2008-07-31
     */ 
    // ========================================================================
    template <class TYPE, 
              typename Ostap::Math::TypeWrapper<TYPE>::value_type I , 
              unsigned int N>
    struct IDigit : public detail::_IDigit<TYPE,I,N> {};
    // ========================================================================
    /** @struct IDigits 
     *  Simple structr efor compile-time evaluation of 
     *  the range of decimal digits N1->N2 (N2 is excluded) 
     *  for the given integral type 
     *  
     *  @code 
     *
     *  const int dig02 = IDigits<unsigned int,4362736,0,2>::value ;
     *
     *  @endcode 
     *
     *  @attention the least significant decimal digit is numbered as #0
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@nilkhef.nl
     *  @date 2008-07-31
     */  
    template <class TYPE, 
              typename Ostap::Math::TypeWrapper<TYPE>::value_type I , 
              unsigned int N1 , 
              unsigned int N2 >
    struct IDigits : public detail::_IDigits<TYPE,I,N1,N2> {};
    // ========================================================================
    /** @struct Digit 
     *  simple structure for evaluation of Nth digit for the integral type 
     * 
     *  @code
     *   
     *   const unsigned int value  = ... ;
     * 
     *   Digit<unsigned int,5> digit5 ;
     *  
     *   const int dig5 = digit5 ( value ) ;
     *
     *  @endcode 
     *
     *
     *  @attention the least significat decimal digit is numbered as #0
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
     *  @date 2008-07-31 
     */
    template <class TYPE, unsigned int N>
    struct Digit : public detail::_Digit<TYPE,N> {};
    // ========================================================================
    /** @struct Digits 
     *  simple structure for evaluation of 
     *  range of decomal digits N1->N2 (N2 is excluded) for the integral type 
     * 
     *  @code
     *   
     *   const unsigned int value  = ... ;
     * 
     *   Digits<unsigned int,0,5> digits05 ;
     *  
     *   const int dig05 = digit05 ( value ) ;
     *
     *  @endcode 
     *
     *  @attention the least significat decimal digit is numbered as #0
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
     *  @date 2008-07-31 
     */
    template <class        TYPE , 
              unsigned int N1   , 
              unsigned int N2   >
    struct Digits : public detail::_Digits<TYPE,N1,N2> {};
    // ========================================================================
    /** simple function which evaluate N-th decimal digit for the 
     *  integral value 
     *
     *  @code
     *
     *   unsigned int value = ... ;
     *   unsigned int N     = ... ;
     * 
     *   const int digN = Ostap::Math::digit ( value , N ) ;
     *
     *  @endcode 
     * 
     *  @attention the least significat decimal digit is numbered as #0
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
     *  @date 2008-07-09
     */
    template <class TYPE>
    inline TYPE digit ( const TYPE value , const unsigned int N  ) 
    {
      // ======================================================================
      static_assert (  std::numeric_limits<TYPE>::is_specialized && 
                       std::numeric_limits<TYPE>::is_integer     && 
                       !std::numeric_limits<TYPE>::is_signed      , 
                       "digit: inappropriate type"                 ) ;
      // ======================================================================
      if      (  N > (unsigned int) std::numeric_limits<TYPE>::digits10 ) { return 0 ; } // RETURN 
      else if (  N < (unsigned int) std::numeric_limits<TYPE>::digits10 ) 
        {
        // ====================================================================
        const TYPE ten = 10 ;
        const TYPE aux = Ostap::Math::pow ( ten , N ) ;
        return static_cast<TYPE> ( ( value/aux ) % 10 ) ;               // RETURN 
        // ====================================================================
      }
      // ======================================================================
      const unsigned long long val = value ;
      const unsigned long long ten = 10    ;
      const unsigned long long aux = Ostap::Math::pow ( ten , N ) ;
      return static_cast<TYPE> (  (val/aux) % 10 ) ;                  // RETURN
      // ======================================================================
    }  
    // ========================================================================
    /** simple function which evaluates the range of decimal digits 
     *  N1->N2 (N2 is excluded) for the integral values 
     *
     *  @code
     *
     *   unsigned int value = ... ;
     *   unsigned int N1 = ... ;
     *   unsigned int N2 = ... ;
     * 
     *   const bool digN1N2 = Ostap::Math::digits( value , N1 , N2  ) ;
     *
     *  @endcode 
     * 
     *  @attention the least significat decimal digit is numbered as #0
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
     *  @date 2008-07-09
     */
    template <class TYPE>
    inline TYPE digits ( const TYPE         value , 
                         const unsigned int N1    ,
                         const unsigned int N2    ) 
    {
      // ======================================================================
      static_assert (  std::numeric_limits<TYPE>::is_specialized && 
                       std::numeric_limits<TYPE>::is_integer     && 
                       !std::numeric_limits<TYPE>::is_signed      , 
                       "digits: inappropriate type"                ) ;
      // ======================================================================
      if      (  N2 >  1 + std::numeric_limits<TYPE>::digits10 )
      { return digits ( value , N1 , 1 + std::numeric_limits<TYPE>::digits10 ) ; }
      // ======================================================================
      if      (  N1 >= N2 || 
                 N1 > (unsigned int) std::numeric_limits<TYPE>::digits10 )
      { return 0 ; }                                                  // RETURN 
      //
      if ( N1      < (unsigned int) std::numeric_limits<TYPE>::digits10 && 
           N2 - N1 < (unsigned int) std::numeric_limits<TYPE>::digits10 ) 
      {
        // ====================================================================
        const TYPE ten = 10 ;
        const TYPE aux1 = Ostap::Math::pow ( ten ,      N1 ) ;
        const TYPE aux2 = Ostap::Math::pow ( ten , N2 - N1 ) ;
        return static_cast<TYPE> ( (value/aux1)%aux2 ) ;              // RETURN 
        // ====================================================================
      }
      // ======================================================================
      const unsigned long long val  = value ;
      const unsigned long long ten  = 10 ;
      const unsigned long long aux1 = Ostap::Math::pow ( ten , N1      ) ;
      const unsigned long long aux2 = Ostap::Math::pow ( ten , N2 - N1 ) ;
      return static_cast<TYPE>( (val/aux1)%aux2 ) ;                   // RETURN
      // ======================================================================
    }
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_DIGIT_H
// ============================================================================
