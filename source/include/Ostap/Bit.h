// $Id$
// ============================================================================
#ifndef OSTAP_BIT_H 
#define OSTAP_BIT_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL 
// ============================================================================
#include <functional>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/TypeWrapper.h"
// ============================================================================
// Boost
// ============================================================================
#include "boost/integer_traits.hpp"
#include "boost/static_assert.hpp"
// ============================================================================
/** @file
 *
 *  The collection of utilities for evaluation of bits for 
 *  the unsigned integral values 
 *
 *  The utilities for evaluation of Nth bit for 
 *  the unsigned integral values: 
 * 
 *   -  Ostap::Math::IBit 
 *      The most efficient compile-time evaluator.
 *      It is applicable if both N and value are compile-time constanst 
 *   -  Ostap::Math::Bit
 *      Rather efficient fuctor, N is compile-time constanst 
 *   -  Ostap::Math::bit
 *      The regular function, the least efficient evaluator 
 *  
 *  The utilities for evaluation of range of bits [N1,N2) for 
 *  the unsigned integral values. N2 is *NOT* included:
 *
 *   -  Ostap::Math::IBits 
 *      The most efficient compile-time evaluator.
 *      It is applicable if both N1&N2 and value are compile-time constanst 
 *   -  Ostap::Math::Bits
 *      Rather efficient fuctor, N1&N2 are compile-time constanst 
 *   -  Ostap::Math::bits
 *      The regular function, the least efficient evaluator 
 * 
 *  @attention the least significant bit is numbered as #0
 *
 *  @see Ostap::Math::IBit
 *  @see Ostap::Math::Bit 
 *  @see Ostap::Math::bit
 *
 *  @see Ostap::Math::IBits
 *  @see Ostap::Math::Bits 
 *  @see Ostap::Math::bits
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
      /** @struct Check 
       *  Simple structure to check the if the type has sufficient 
       *  length to address Nth bit 
       *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
       *  @date 2008-08-01
       */
      template <class TYPE,unsigned int N>
      struct Check 
      {
        // ====================================================================
        static_assert ( boost::integer_traits<TYPE>::is_specialized
                        && boost::integer_traits<TYPE>::is_integral 
                        &&!boost::integer_traits<TYPE>::is_signed   , 
                        "Check: Inappropriate type") ;      
        // ====================================================================
        enum { value =  N < (unsigned int) boost::integer_traits <TYPE>::digits } ;
        // ====================================================================
      };
      // ======================================================================
      /** @struct _IBit 
       *  Simple helper structure to extract Nth bit of the integral type 
       *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
       *  @date 2008-08-01
       */
      template <class TYPE, 
                typename Ostap::Math::TypeWrapper<TYPE>::value_type I, 
                unsigned int N>
      struct _IBit
      {
      private:
        // ====================================================================
        static_assert ( boost::integer_traits<TYPE>::is_specialized 
                        &&!boost::integer_traits<TYPE>::is_signed   , 
                        "_IBit: Inappropriate type"                 ) ;      
        // ====================================================================
      public:
        // ====================================================================
        enum { value = 0 < ( I & ( static_cast<TYPE>(1) << N ) ) } ;
        // ====================================================================
      };
      // ======================================================================
      /** @struct _IBits 
       *  Simple helper structure to extract N1->N2th 
       *  bits of the integral type 
       *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
       *  @date 2008-08-01
       */
      template <class TYPE,
                typename Ostap::Math::TypeWrapper<TYPE>::value_type I, 
                unsigned int N1,
                unsigned int N2>
      struct _IBits
      {
      private:
        // ====================================================================
        static_assert ( boost::integer_traits<TYPE>::is_specialized 
                        &&!boost::integer_traits<TYPE>::is_signed   , 
                        "_IBits: Inappropriate type"                 ) ;      
        // ====================================================================
      public:
      // ====================================================================
        enum 
          { 
            value = 
            ( I & ( static_cast<TYPE> ( -1 )
                    << ( boost::integer_traits<TYPE>::digits + N1 - N2 ) 
                    >> ( boost::integer_traits<TYPE>::digits      - N2 ) ) ) 
            >> N1
          } ;
        // ====================================================================
      };
      // ======================================================================
      /** @struct _Bit 
       *  Simple helper structure to extract Nth bit of the integral type 
       *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
       *  @date 2008-08-01
       */
      template <class TYPE, unsigned int N>
      struct _Bit : public std::unary_function<TYPE,bool>
      {
      private:
        // ====================================================================
        static_assert ( boost::integer_traits<TYPE>::is_specialized
                        && boost::integer_traits<TYPE>::is_integral 
                        &&!boost::integer_traits<TYPE>::is_signed   , 
                        "_Bit: Inappropriate type"                  ) ;      
        // ====================================================================
      public:
        // ====================================================================
        inline bool operator() ( const TYPE value ) const
        {
          static const TYPE s_mask = static_cast<TYPE> ( 1 ) << N  ;
          return ( (value & s_mask) == 0 ? false:true );
        } 
        // ====================================================================
      };
      // ======================================================================
      /** @struct _Bits 
       *  Simple helper structure to extract N1->N2 bits of the integral type 
       *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
       *  @date 2008-08-01
       */
      template <class TYPE,unsigned int N1,unsigned int N2>
      struct _Bits : public std::unary_function<TYPE,TYPE> 
      {
      private:
        // ====================================================================
        static_assert ( boost::integer_traits<TYPE>::is_specialized
                        && boost::integer_traits<TYPE>::is_integral 
                        &&!boost::integer_traits<TYPE>::is_signed   , 
                        "_Bits: Inappropriate type"                  ) ;      
        // ====================================================================
        static_assert ( N1 < N2  &&
                        N1 <  (unsigned int) boost::integer_traits<TYPE>::digits && 
                        N2 <= (unsigned int) boost::integer_traits<TYPE>::digits , 
                        "_Bits: invalid indices"                     ) ;      
        // ====================================================================
      public:
        // ====================================================================
        inline TYPE operator() ( const TYPE value ) const
        {
          // ==================================================================
          static const TYPE s_mask = 
            ( static_cast<TYPE> ( -1 ) 
              << ( boost::integer_traits<TYPE>::digits + N1 - N2 ) )  
            >> ( boost::integer_traits<TYPE>::digits      - N2 ) ;
          // =================================================================
          return ( value & s_mask ) >> N1 ;
          // ==================================================================
        }
        // ====================================================================
      } ;
      // ======================================================================
    } // end of namespace detail
    // ========================================================================
    /** @struct IBit 
     *  (Compile-time) evaluation of the Nth bit for unsigned integral I:
     *
     *  @code 
     * 
     *  const bool bit0 = Ostap::Math::IBit<111,0>::value ;
     *  const bool bit1 = Ostap::Math::IBit<111,1>::value ;
     *  const bool bit2 = Ostap::Math::IBit<111,2>::value ;
     *  const bool bit3 = Ostap::Math::IBit<111,3>::value ;
     *
     *  @endcode 
     *
     *  @attention the least significant bit is numbered as bit#0 
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
     *  @date 2008-07-09
     */
    template <class TYPE, 
              typename TypeWrapper<TYPE>::value_type I, 
              unsigned int  N>
    struct IBit : public detail::_IBit<TYPE,I,N> { };
    // ========================================================================
    /** @struct IBits 
     *
     *  (Compile-time) Get the N1th to N2th bits (N2 is not included) 
     *  from the unsigned 
     *  integral
     *
     *  @code
     *
     *   // gets bits [0,5) from the number 100 
     *   const unsigned int ibits = IBits<unsigned int,100,0,5>::value ;
     *
     *  @endcode 
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
     *  @date 2008-08-01
     */    
    template <class TYPE,
              typename TypeWrapper<TYPE>::value_type I, 
              unsigned int N1,
              unsigned int N2>
    struct IBits : public detail::_IBits<TYPE,I,N1,N2> {};
    // ========================================================================
    /** @struct Bit
     *  Simple strcuture which evaluates N-th bit for the 
     *  integral value, where N is compile-time constant  
     *     
     *  
     *  @code
     * 
     *  TYPE value = ... ;
     *
     *  Ostap::Math::Bit<TYPE,10> bit ; 
     *  const bool bit10 = bit ( value ) ; // check bit #10 
     *
     *  @endcode 
     *
     *  @attention the least significant bit is numbered as bit#0 
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
     *  @date 2008-07-09
     */
    template <class TYPE, unsigned int N>
    struct Bit : public detail::_Bit<TYPE,N> {};
    // ========================================================================
    /** @struct Bits 
     *
     *  Get the N1th to N2th bits (N2 is not included) from the unsigned 
     *  integral
     *
     *  @code
     * 
     *   const unsigned int value = ... ;
     *  
     *   Bits<unsigned int,0,5> bits05 ;
     * 
     *   const bits = bits05 ( value ) ; 
     *
     *  @endcode 
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
     *  @date 2008-08-01
     */
    template <class TYPE, unsigned int N1, unsigned int N2>
    struct Bits : public detail::_Bits<TYPE,N1,N2> {};
    // ========================================================================
    /** simple function which evaluates N-th bit for the 
     *  integral value 
     *
     *  @code
     * 
     *  TYPE value = ... ;
     * 
     *  const bool bit10 = Ostap::Math::bit ( value , 10 ) ; // check bit #10 
     *
     *  @endcode 
     *
     *  @attention the least significant bit is numbered as bit#0 
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
     *  @date 2008-07-09
     */
    template <class TYPE>
    inline bool bit ( const TYPE value , const unsigned int N  ) 
    {
      // ======================================================================
      static_assert ( boost::integer_traits<TYPE>::is_specialized
                      && boost::integer_traits<TYPE>::is_integral 
                      &&!boost::integer_traits<TYPE>::is_signed   , 
                      "bit: inappropriate type"                   ) ;
      // ======================================================================
      static const TYPE one = static_cast<TYPE>(1) ;
      return N < (unsigned int) boost::integer_traits <TYPE>::digits ? 
        ( value & ( one << N ) )!=0 : 0 ;
      // ======================================================================
    }  
    // ========================================================================
    template <class TYPE>
    inline TYPE bits ( const TYPE         value , 
                       const unsigned int N1    , 
                       const unsigned int N2    )
    {
      // ======================================================================
      static_assert ( boost::integer_traits<TYPE>::is_specialized
                      && boost::integer_traits<TYPE>::is_integral 
                      &&!boost::integer_traits<TYPE>::is_signed   , 
                      "bits: inappropriate type"                   ) ;
      // ======================================================================
      if ( N2 >  (unsigned int) boost::integer_traits<TYPE>::digits ) 
      { return bits ( value , N1 , boost::integer_traits<TYPE>::digits ) ; } 
      // 
      if ( N1 >= N2 ||
           N1 >= (unsigned int) boost::integer_traits<TYPE>::digits ) { return 0 ; }
      // ======================================================================
      const TYPE mask = 
        ( static_cast<TYPE> ( -1 ) 
          << ( boost::integer_traits<TYPE>::digits + N1 - N2 ) )
        >> ( boost::integer_traits<TYPE>::digits      - N2 ) ;
      // ======================================================================
      return ( value & mask ) >> N1 ;
      // ======================================================================
    }
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
// The END 
// ============================================================================
#endif // OSTAP_BIT_H
// ============================================================================
