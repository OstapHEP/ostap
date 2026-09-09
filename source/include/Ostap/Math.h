// ============================================================================
#ifndef OSTAP_MATH_H 
#define OSTAP_MATH_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL 
// ============================================================================
#include <cstdint>
#include <complex>
#include <cmath>
#include <limits>
#include <utility>
#include <iterator>
#include <type_traits>
#include <algorithm>
#include <functional>
#include <vector>
#include <array>
#include <numeric>
#if defined ( __cplusplus ) && defined ( __cpp_lib_math_constants ) && ( 201907L <= __cpp_lib_math_constants )
#include <numbers>
#endif
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Lomont.h"
#include "Ostap/Power.h"
#include "Ostap/Epsilon.h"
#include "Ostap/Param_t.h"
#include "Ostap/mULPS.h"
#include "Ostap/Comparisons.h"
// ============================================================================
/** @file Ostap/Math.h
 *  collection of generic math functions and classes 
 */
// ===========================================================================
namespace Ostap
{
  // ==========================================================================
  /** @namespace Ostap::Math Math.h
   *  collection of generic math functions and classes 
   */
  namespace Math 
  {
    // ========================================================================
    // Parameters for numerical calculations (M.Needham)
    // ========================================================================
    /// High tolerance
    constexpr double hiTolerance    = 1e-40 ;
    /// Low  tolerance
    constexpr double lowTolerance   = 1e-20 ;
    /// Loose tolerance
    constexpr double looseTolerance = 1e-6  ;
    // ========================================================================
    
    // ======================================================================
    /// get \f$ s(n) = (-1)^{n} \f$  
    constexpr inline std::int8_t sign ( const int n ) { return ( 1 == n % 2 ) ? -1 : +1 ; }
    // ======================================================================
    
    // ========================================================================
    /** return "min_by_abs"
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2007-08-17
     */        
    template <class TYPE> 
    inline TYPE absMin
    ( param_t<TYPE> v1 ,
      param_t<TYPE> v2 ) 
    { return std::min ( std::abs ( v1 ) , std::abs ( v2 ) ) ; }
    // ========================================================================
    /** return  "max_by_abs"
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2007-08-17
     */
    template <class TYPE> 
    inline TYPE absMax
    ( param_t<TYPE> v1 ,
      param_t<TYPE> v2 ) 
    { return std::max ( std::abs ( v1 ) , std::abs ( v2 ) ) ; }
    // ========================================================================

    // ========================================================================
    template <class TYPE> struct    Zero ;
    template <class TYPE> struct NonZero ;
    // ========================================================================    
    
    // ========================================================================
    /// partial specialisation for the complex values 
    template <class TYPE>
    struct Tiny< std::complex<TYPE> >
    {
    public:
      // ======================================================================
      // constructor 
      Tiny ( TYPE  b ) : m_tiny ( b ) {} ;
      // ======================================================================
      /// comparison
      inline bool operator() ( const std::complex<TYPE>& v ) const 
      { return m_tiny ( v.real() ) && m_tiny ( v.imag () ) ; }
      // ======================================================================
    private :
      // ======================================================================
      /// default consyruictor is disabled 
      Tiny () ; // default consyruictor is disabled 
      // ======================================================================      
    private:
      // ======================================================================
      // the comparison criteria 
      Tiny <TYPE> m_tiny ;
      // ======================================================================
    } ;

    // ========================================================================
    /** round to nearest integer, rounding half a  way from  zero 
     *  @see std::llround 
     *  @see https://en.cppreference.com/w/cpp/numeric/math/round
     *  @author Vanya BELYAEV Ivan.Belyaev
     */
    std::intmax_t round ( const double x ) ;
    // ========================================================================
    /** round to nearest integer, rounding half a  way from  zero 
     *  @see std::llround 
     *  @see https://en.cppreference.com/w/cpp/numeric/math/round
     *  @author Vanya BELYAEV Ivan.Belyaev
     */
    std::intmax_t round ( const long double x ) ; 
    // ========================================================================
    /** round to nearest integer, rounding half a way from zero 
     *  @author Vanya BELYAEV Ivan.Belyaev
     */
    inline std::intmax_t round ( const float  x ) { return round ( 1.0 * x ) ; }
    // ========================================================================
    /** "Round" the complex value: 
     *  round both the real and imaginary components
     *  @code
     *  std::complex<double> x = ...
     *  std::complex<double> r = round ( x ) ;      
     *  @endcode
     *  @see Ostap::Math::round
     */
    template <class TYPE>
    inline std::complex<TYPE> round
    ( const std::complex<TYPE>& z )
    {
      /// round real & imaginary parts 
      const TYPE re_z { round ( z.real () ) } ;
      const TYPE im_z { round ( z.imag () ) } ;
      //
      return std::complex<TYPE> ( re_z , im_z ) ;
    }
    // ========================================================================
    /** Round  down
     *  @see std::floor 
     *  @see  https://en.wikipedia.org/wiki/Rounding 
     */
    inline std::intmax_t round_down ( const double v )
    { return round ( std::floor ( v ) ) ; } 
    // ========================================================================
    /** Round  up 
     *  @see std::ceil 
     *  @see  https://en.wikipedia.org/wiki/Rounding 
     */
    inline std::intmax_t round_up  ( const double v )
    { return round ( std::ceil  ( v ) ) ; } 
    // =========================================================================
    /** Round  toward zero
     *  @see std::trunc  
     *  @see  https://en.wikipedia.org/wiki/Rounding 
     */
    inline std::intmax_t round_toward_zero ( const double v )
    { return round ( std::trunc ( v ) ) ; } 
    // ==========================================================================
    /** Round  away zero
     *  @see std::ceil    
     *  @see  https://en.wikipedia.org/wiki/Rounding 
     */
    inline std::intmax_t round_away_zero ( const double v )
    { 
      const double av = std::ceil ( std::abs ( v ) ) ;
      return round ( 0 <= v ? av : -av ) ; 
    } 
    // ==========================================================================
    /** Round half-up 
     *  @see std::floor    
     *  @see  https://en.wikipedia.org/wiki/Rounding 
     */
    inline std::intmax_t round_half_up ( const double v )
    { return round ( std::floor ( v + 0.5 ) ) ;  } 
    // ==========================================================================
    /** Round half-down 
     *  @see std::ceil    
     *  @see  https://en.wikipedia.org/wiki/Rounding 
     */
    inline std::intmax_t round_half_down ( const double v )
    { return round ( std::ceil ( v - 0.5 ) ) ;  }
    // ========================================================================
    /** Round half-toward sero 
     *  @see std::floor    
     *  @see std;:ceil 
     *  @see  https://en.wikipedia.org/wiki/Rounding 
     */
    inline std::intmax_t round_half_toward_zero ( const double v )
    { return round ( 0 <= v ? std::ceil ( v - 0.5 ) : std::floor ( v + 0.5 ) ) ;  } 
    // ========================================================================
    /** Round half-away sero 
     *  @see std::floor    
     *  @see std;:ceil 
     *  @see  https://en.wikipedia.org/wiki/Rounding 
     */
    inline std::intmax_t round_half_away_zero ( const double v )
    { return round ( 0 <= v ? std::floor ( v + 0.5 ) : std::ceil ( v - 0.5 ) ) ;  } 
    // ========================================================================
    /** Round half to even  
     *  Aka: 
     *  -  convergent rounding
     *  - statistician's rounding, 
     *  - Dutch rounding, 
     *  - Gaussian rounding, 
     *  - odd–even rounding,
     *  - bankers' rounding.
     *  @see std::lrint
     */    
    std::intmax_t round_half_even ( const double v ) ; 
    // =======================================================================
    /** Round half to even  
     *  Aka: 
     *  - convergent rounding
     *  - statistician's rounding, 
     *  - Dutch rounding, 
     *  - Gaussian rounding, 
     *  - odd–even rounding,
     *  - bankers' rounding.
     *  @see std::lrint
     */    
    inline std::intmax_t banker ( const double v ) { return round_half_even ( v ) ; }  
    // ========================================================================
    /** get mantissa and (decimal) exponent 
     *  similar to std::frexp, but radix=10)
     *  @param x  INPUT  value 
     *  @param e  UPDATE exponent 
     *  @return  mantissa     (0.1<=m<1)
     *  @author Vanya BELYAEV Ivan.Belyaev       
     *  @date 2015-07-21
     */
    double frexp10 ( const double x , int& e ) ;
    // ========================================================================
    /** get mantissa and (decimal) exponent 
     *  similar to std::frexp, but radix=10)
     *  @param x  INPUT  value 
     *  @param e  UPDATE exponent 
     *  @return  mantissa    (0.1<=m<1) 
     *  @author Vanya BELYAEV Ivan.Belyaev       
     *  @date 2015-07-21
     */
    float frexp10 ( const float x , int& e ) ;
    // ========================================================================
    /** get mantissa and (decimal) exponent 
     *  similar to std::frexp, but radix=10)
     *  @param x  INPUT  value 
     *  @return   pair of mantissa (0.1<=m<1) and (decimal) exponent 
     *  @author Vanya BELYAEV Ivan.Belyaev       
     *  @date 2015-07-21
     */
    std::pair<double,int>
    frexp10 ( const double x ) ;
    // ========================================================================
    /** get mantissa and binary exponent 
     *  similar to std::frexp
     *  @param x  INPUT  value 
     *  @return   pair of mantissa (0.5<=m<1) and (binary) exponent 
     *  @author Vanya BELYAEV Ivan.Belyaev       
     *  @date 2015-07-21
     */
    std::pair<double,int>
    frexp2  ( const double x ) ;
    // ========================================================================
    /** Multiplies a floating point value num by the number 10(radix)
     *  raised to the exp power. 
     *  similar to std::lxep, but radix is 10 
     *  @param num  input value 
     *  @param exp  the 10-base exponent 
     *  @return     propertly scaled input value 
     *  @attention  It is not very efficient, but OK for our purpose
     */
    double ldexp10
    ( const double value ,
      const short  expo  ) ;
    // ========================================================================
    /** Round to N-significant digits 
     *  @param x  INPUT  input value 
     *  @param n  INPUT  number of significnat digits 
     *  @return rounded value 
     *  @author Vanya BELYAEV Ivan.Belyaev       
     *  @date 2015-07-21
     */
    double round_N ( const double x , const unsigned short n ) ;
    // ========================================================================
    /** Round to N-significant digits 
     *  @param x  INPUT  input value 
     *  @param n  INPUT  number of significnat digits 
     *  @return rounded value 
     *  @author Vanya BELYAEV Ivan.Belyaev       
     *  @date 2015-07-21
     */
    float round_N ( const float x , const unsigned short n ) ;
    // ========================================================================
    /** Is the value actually long ?
     *  @author Vanya BELYAEV Ivan.Belyaev       
     *  @date 2011-07-18
     */
    bool islong ( const long double x ) ;
    // ========================================================================
    /** Is the value actually long ?
     *  @author Vanya BELYAEV Ivan.Belyaev       
     *  @date 2011-07-18
     */
    bool islong ( const double x ) ;
    // ========================================================================
    /** Is the value actually long ?
     *  @author Vanya BELYAEV Ivan.Belyaev       
     *  @date 2011-07-18
     */
    bool islong ( const float  x ) ;
    // ========================================================================
    /** Is the value actually int ?
     *  @author Vanya BELYAEV Ivan.Belyaev       
     *  @date 2011-07-18
     */
    bool isint  ( const long double x ) ;
    // ========================================================================
    /** Is the value actually int ?
     *  @author Vanya BELYAEV Ivan.Belyaev       
     *  @date 2011-07-18
     */
    bool isint  ( const double x ) ;
    // ========================================================================    
    /** Is the value actually int ?
     *  @author Vanya BELYAEV Ivan.Belyaev       
     *  @date 2011-07-18
     */
    bool isint  ( const float  x ) ;
    // ========================================================================
    /** Is the value actually unsigned int ?
     *  @author Vanya BELYAEV Ivan.Belyaev       
     *  @date 2011-07-18
     */
    bool isuint     ( const double x ) ;
    // =========================================================================
    /// Is it a short value ?
    bool isshort    ( const double x ) ;
    // =========================================================================
    /// Is it an unsigned  short value ?
    bool isushort   ( const double x ) ;
    // =========================================================================
    /** Is the value actually unsigned long ?
     *  @author Vanya BELYAEV Ivan.Belyaev       
     *  @date 2011-07-18
     */
    bool isulong    ( const double x ) ;
    // =========================================================================
    /// Is the value actually long long?
    bool islonglong ( const double x ) ;
    // ========================================================================
    /// Is the value actually unsigned  long long?
    bool isulonglong ( const double x ) ;
    // ========================================================================
    /// Is floating value actually std::int8_t ?
    bool isint8 ( const double x ) ;
    // ========================================================================
    /// Is floating value actually std::int16_t ?
    bool isint16 ( const double x ) ;
    // ========================================================================
    /// Is floating value actually std::int32_t ?
    bool isint32 ( const double x ) ;
    // ========================================================================
    /// Is floating value actually std::int64_t ?
    bool isint64 ( const double x ) ;
    // ========================================================================
    /// Is floating value actually std::uint8_t ?
    bool isuint8 ( const double x ) ;
    // ========================================================================
    /// Is floating value actually std::uint16_t ?
    bool isuint16 ( const double x ) ;
    // ========================================================================
    /// Is floating value actually std::uint32_t ?
    bool isuint32 ( const double x ) ;
    // ========================================================================
    /// Is floating value actually std::uint64_t ?
    bool isuint64 ( const double x ) ;
    // ========================================================================
    /// Is floating value actually char  ?
    bool ischar ( const double x ) ;
    // ========================================================================
    /// Is floating value actually signed char  ?
    bool isschar ( const double x ) ;
    // ========================================================================
    /// Is floating value actually unsigned char  ?
    bool isuchar ( const double x ) ;
    // ========================================================================
    /// Is floating value actually std::intmax_t ? 
    bool isintmax  ( const double x ) ; 
    /// Is floating value actually std::uintmax_t ? 
    bool isuintmax ( const double x ) ; 
    /// Is floating value actually std::size_t  ? 
    bool issize    ( const double x ) ; 
    // ========================================================================
    /// Is floating value actually std::intptr_t 
    bool isintptr  ( const double x ) ; 
    // ========================================================================
    /// Is floating value actually std::uintptr_t 
    bool isuintptr ( const double x ) ; 
    // ========================================================================
    
    // ========================================================================    
    /** check if the double value is actually equal to the integer value  
     *  @param val value to be compared with the integer 
     *  @param ref the reference integer number 
     *  @param mULPS the precision 
     *  @see Ostap::Math::Lomont::compare_double 
     *  @see Ostap::Math::mULPS_double
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2008-09-17
     */
    bool equal_to_int 
    ( const double       val                  , 
      const int          ref                  , 
      const unsigned int mULPS = mULPS_double ) ;
    // ========================================================================
    /** check if the double value is actually equal to the integer value  
     *  @param ref the reference integer  number 
     *  @param val value to be compared with the integer 
     *  @param mULPS the precision 
     *  @see Ostap::Math::Lomont::compare_double 
     *  @see Ostap::Math::mULPS_double
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2008-09-17
     */        
    inline bool equal_to_int 
    ( const int          ref                  , 
      const double       val                  , 
      const unsigned int mULPS = mULPS_double ) 
    { return equal_to_int ( val , ref , mULPS ) ; }
    // ========================================================================
    /** check if the double value is actually equal to the unsigned integer value  
     *  @param val value to be compared with the unsigned integer 
     *  @param ref the reference unsigned integer number 
     *  @param mULPS the precision 
     *  @see Ostap::Math::Lomont::compare_double 
     *  @see Ostap::Math::mULPS_double
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2008-09-17
     */        
    bool equal_to_uint 
    ( const double       val                  , 
      const unsigned int ref                  , 
      const unsigned int mULPS = mULPS_double ) ;
    // ========================================================================
    /** check if the double value is actually equal to the integer value  
     *  @param val value to be compared with the unsigned integer 
     *  @param ref the reference unsigned integer number 
     *  @param mULPS the precision 
     *  @see Ostap::Math::Lomont::compare_double 
     *  @see Ostap::Math::mULPS_double
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2008-09-17
     */        
    inline bool equal_to_uint 
    ( const unsigned int ref                  , 
      const double       val                  ,
      const unsigned int mULPS = mULPS_double ) 
    { return equal_to_uint ( val , ref , mULPS ) ; }
    // ========================================================================
    /// signed sqrt 
    inline double signed_sqrt ( const double value ) 
    { return 0 <= value ? std::sqrt ( value ) : -std::sqrt( std::abs ( value ) ) ; }  
    // ========================================================================
    /** make dot-multiplication of two sequences using std::fma 
     *  \f$ r = \sum_i  x_i y_i \f$
     *  @param begin  (INPUT) start of the first sequence 
     *  @param end    (INPUT) end   of the first sequence 
     *  @param begin2 (INPUT) start iterator of the second sequence
     *  @return   "dot" product of two sequences 
     *
     *  @see http://en.cppreference.com/w/cpp/numeric/math/fma
     *  "...the function std::fma evaluates faster 
     *   (in addition to being more precise) than the expression x*y+z for 
     *   float, double, and long double arguments, respectively. "
     */
    template <class ITERATOR1, class ITERATOR2>
    inline double dot_fma
    ( ITERATOR1 begin  , 
      ITERATOR1 end    ,
      ITERATOR2 begin2 ) 
    {
      long double dot = 0 ;
      for ( ; begin != end ; ++begin, ++begin2 ) 
      {
        const long double x1 = *begin  ;
        const long double x2 = *begin2 ;        
        dot = std::fma ( x1 , x2 , dot ) ; 
      }
      return dot  ;  
    }
    // ========================================================================
    /** make dot-multiplication of two sequences using std::fma 
     *  \f$ r = \sum_i  x_i y_i \f$
     *  @param x     (INPUT) the first sequence 
     *  @param begin (INPUT) start iterator of the second sequence
     *  @return   "dot" product of two sequences 
     */
    template <unsigned int N, class TYPE, class ITERATOR>
    inline double dot_fma
    ( const std::array<TYPE,N>& x     , 
      ITERATOR                  begin )  
    { return dot_fma ( x.begin() , x.end() , begin ) ; }
    // ========================================================================
    /** make dot-multiplication of two sequences using std::fma 
     *  \f$ r = \sum_i  x_i y_i \f$
     *  @param x (INPUT) the first sequence 
     *  @param y (INPUT) the second sequence
     *  @return   "dot" product of two sequences 
     */
    template <unsigned int N, class TYPE1, class TYPE2>
    inline double dot_fma
    ( const std::array<TYPE1,N>& x , 
      const std::array<TYPE2,N>& y )  
    { return dot_fma ( x.begin() , x.end() , y.begin() ) ; }
    // ========================================================================
    /** make dot-multiplication of two sequences using std::fma 
     *  \f$ r = \sum_i  x_i y_i \f$
     *  @param x (INPUT) the first sequence 
     *  @param y (INPUT) begin-iterator for the second sequence 
     *  @return   "dot" product of two sequences 
     */
    template <class TYPE,unsigned int N, class ITERATOR> 
    inline double dot_fma 
    ( TYPE(&x)[N] , 
      ITERATOR y  ) { return dot_fma ( x , x + N , y ) ; }
    // ========================================================================
    /** make dot-multiplication of two sequences using std::fma 
     *  \f$ r = \sum_i  x_i y_i \f$
     *  @param x (INPUT) the first sequence 
     *  @param y (INPUT) the second sequence 
     *  @return   "dot" product of two sequences 
     */
    template <class TYPE1, class TYPE2, unsigned int N> 
    inline double dot_fma 
    ( TYPE1(&x)[N] , 
      TYPE2(&y)[N] ) { return dot_fma ( x , x + N , y ) ; }
    // ========================================================================
    /** make dot-multiplication of two sequences using std::fma 
     *  \f$ r = \sum_i  x_i y_i \f$
     *  @param N (INPUT) length of the sequences 
     *  @param x (INPUT) the first sequence 
     *  @param y (INPUT) the second sequence 
     *  @param skipx skipx  skip first <code> skipx</code> elements from x 
     *  @param skipy skipy  skop first <code> skipy</code> elements from y 
     *  @return   "dot" product of two sequences 
     */
    inline double dot_fma_ 
    ( const unsigned int N         , 
      const double*      x         , 
      const double*      y         , 
      const unsigned int skipx     ,  
      const unsigned int skipy     ) 
    { return dot_fma ( x + skipx     , 
                       x + skipx + N ,
                       y + skipy     ) ; }    
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
    inline double sum_kahan
    ( ITERATOR begin , 
      ITERATOR end   )
    {
      long double sum = 0 ;
      long double c   = 0 ;
      for ( ; begin != end ; ++begin ) 
      {
        volatile const long double y = (*begin) - c ;
        volatile const long double t = sum      + y ;
        c        = ( t - sum ) - y ;
        sum      =   t             ;
      }
      return sum ;
    }
    // ========================================================================
    /** make dot-multiplication of two sequences based on Kahan summation 
     *  \f$ r = \sum_i  x_i y_i \f$
     *  @param xbegin (INPUT) begin-iterator for the first sequence  
     *  @param xend   (INPUT) end-iterator for the first sequence 
     *  @param ybegin (INPUT) begin-iterator for the second sequence  
     */
    template <class ITERATOR1, class ITERATOR2>
    inline double dot_kahan
    ( ITERATOR1 xbegin , 
      ITERATOR1 xend   , 
      ITERATOR2 ybegin )
    {
      long double sum = 0 ;
      long double c   = 0 ;
      for ( ; xbegin != xend ; ++xbegin , ++ybegin ) 
      {
        const          long double v = (*xbegin ) * (*ybegin) ;
        volatile const long double y = v   - c ;
        volatile const long double t = sum + y ;
        c   = ( t - sum ) - y ;
        sum =   t             ;
      }
      return sum ;
    }
    // ========================================================================
    /** make dot-multiplication of two sequences using Kahan summation
     *  \f$ r = \sum_i  x_i y_i \f$
     *  @param x     (INPUT) the first sequence 
     *  @param begin (INPUT) start iterator of the second sequence
     *  @return   "dot" product of two sequences 
     */
    template <unsigned int N, class TYPE, class ITERATOR>
    inline double dot_kahan
    ( const std::array<TYPE,N>& x     , 
      ITERATOR                  begin )  
    { return dot_kahan ( x.begin() , x.end() , begin ) ; }
    // ========================================================================
    /** make dot-multiplication of two sequences using Kahan summation 
     *  \f$ r = \sum_i  x_i y_i \f$
     *  @param x (INPUT) the first sequence 
     *  @param y (INPUT) the second sequence
     *  @return   "dot" product of two sequences 
     */
    template <unsigned int N, 
              class TYPE1 , 
              class TYPE2 , 
              typename std::enable_if<std::is_convertible<TYPE1,long double>::value,bool>::type = true , 
              typename std::enable_if<std::is_convertible<TYPE2,long double>::value,bool>::type = true >
    inline double dot_kahan
    ( const std::array<TYPE1,N>& x , 
      const std::array<TYPE2,N>& y )  
    { return dor_kahan ( x.begin() , x.end() , y.begin() ) ; }
    // ========================================================================    
    /** make dot-multiplication of two sequences using Kahan summation
     *  \f$ r = \sum_i  x_i y_i \f$
     *  @param x (INPUT) the first sequence 
     *  @param y (INPUT) begin-iterator for the second sequence 
     *  @return   "dot" product of two sequences 
     */
    template <class TYPE     ,
              std::size_t  N , 
              class ITERATOR , 
              typename value_type = typename std::iterator_traits<ITERATOR>::value_type     ,
              typename std::enable_if<std::is_convertible<TYPE,long double>::value,bool>::type       = true  , 
              typename std::enable_if<std::is_convertible<value_type,long double>::value,bool>::type = true>
    inline double dot_kahan 
    ( TYPE(&x)[N] , 
      ITERATOR y  )
    { return dot_kahan ( x , x + N , y ) ; }
   // ========================================================================
    /** make dot-multiplication of two sequences using Kahan summation
     *  \f$ r = \sum_i  x_i y_i \f$
     *  @param x (INPUT) the first sequence 
     *  @param y (INPUT) the second sequence 
     *  @return   "dot" product of two sequences 
     */
    template <class TYPE1, 
              class TYPE2, 
              std::size_t  N ,
              typename std::enable_if<std::is_convertible<TYPE1,long double>::value,bool>::type = true  , 
              typename std::enable_if<std::is_convertible<TYPE2,long double>::value,bool>::type = true  > 
    inline double dot_kahan
    ( TYPE1(&x)[N] , 
      TYPE2(&y)[N] )
    { return dot_kahan ( x , x + N , y ) ; }
    // ========================================================================
    /** make dot-multiplication of two sequences using Kahan summation
     *  \f$ r = \sum_i  x_i y_i \f$
     *  @param N (INPUT) length of the sequences 
     *  @param x (INPUT) the first sequence 
     *  @param y (INPUT) the second sequence 
     *  @param skipx skipx  skip first <code> skipx</code> elements from x 
     *  @param skipy skipy  skop first <code> skipy</code> elements from y 
     *  @return   "dot" product of two sequences 
     */
    inline double dot_kahan_
    ( const unsigned int N     , 
      const double*      x     , 
      const double*      y     , 
      const std::size_t  skipx ,
      const std::size_t  skipy )
    { return dot_kahan ( x + skipx     ,
                         x + skipx + N , 
                         y + skipy     ) ; }
    // ========================================================================
    /// simple scaling of elements of non-constant sequence        
    template <class ITERATOR  , 
              typename SCALAR ,
              typename value_type = typename std::iterator_traits<ITERATOR>::value_type     ,
              typename std::enable_if<std::is_convertible<value_type,long double>::value,bool>::type = true , 
              typename std::enable_if<std::is_convertible<SCALAR,long double>::value,bool>::type     = true >
    void scale
    ( ITERATOR first  ,
      ITERATOR last   , 
      SCALAR   factor )
    { for ( ; first != last ; ++first ) { (*first) *= factor ; } }
    // ========================================================================
    /// shift all elements of non-constant sequence        
    template <class    ITERATOR, 
              typename SCALAR  ,
              typename value_type = typename std::iterator_traits<ITERATOR>::value_type     ,
              typename std::enable_if<std::is_convertible<value_type,long double>::value,bool>::type = true , 
              typename std::enable_if<std::is_convertible<SCALAR,long double>::value,bool>::type     = true >
    void shift 
    ( ITERATOR first  ,
      ITERATOR last   , 
      SCALAR   factor )
    { for ( ; first != last ; ++first ) { (*first) += factor ; } }
    // ========================================================================
    /// simple scale and shift of elements of non-constant sequence        
    template <class    ITERATOR, 
              typename SCALAR  ,
              typename value_type = typename std::iterator_traits<ITERATOR>::value_type     ,
              typename std::enable_if<std::is_convertible<value_type,long double>::value,bool>::type = true , 
              typename std::enable_if<std::is_convertible<SCALAR,long double>::value,bool>::type     = true >
    void scale_and_shift 
    ( ITERATOR first ,
      ITERATOR last  , 
      SCALAR   scale ,
      SCALAR   shift )
    { for ( ; first != last ; ++first ) 
      { (*first) = std::fma ( *first , scale , shift ) ; } }
    // ========================================================================
    /// scale & shift all elements of vector 
    template <class TYPE     , 
              class ALLOCATOR, 
              class SCALAR   , 
              typename std::enable_if<std::is_convertible<TYPE  ,long double>::value,bool>::type = true , 
              typename std::enable_if<std::is_convertible<SCALAR,long double>::value,bool>::type = true >    
    inline void scale_and_shift
    ( std::vector<TYPE,ALLOCATOR>& vct , 
      SCALAR scale , 
      SCALAR shift ) 
    { scale_and_shift ( vct.begin() , vct.end() , scale , shift ) ; }
    // ========================================================================
    /// simple scaling of exponents for all elements of non-constant sequence        
    template <class ITERATOR> 
    void scale_exp2
    ( ITERATOR    first ,
      ITERATOR    last  , 
      const short iexp  )
    { if ( 0 != iexp ) { for ( ; first != last ; ++first ) { (*first) = std::ldexp ( *first , iexp ) ; } } }
    // ========================================================================
    /// scale all elements of vector 
    template <class    TYPE      ,
              class    ALLOCATOR , 
              typename SCALAR    >
    void scale
    ( std::vector<TYPE,ALLOCATOR>& vct , SCALAR factor ) 
    { scale    ( vct.begin() , vct.end () , factor ) ; }
    // ========================================================================
    /// scale all elements of vector 
    template <class TYPE, class ALLOCATOR>    
    inline void scale_exp2
    ( std::vector<TYPE,ALLOCATOR>& vct , const int iexp  ) 
    { scale_exp2 ( vct.begin() , vct.end () , iexp ) ; }
    // ========================================================================
    /// scale all elements of vector  by 2**s
    template <class TYPE, class ALLOCATOR>    
    inline
    std::vector<TYPE,ALLOCATOR> 
    ldexp ( std::vector<TYPE,ALLOCATOR> vct , const short iexp )
    {
      if ( 0 != iexp ) { scale_exp2 ( vct , iexp ) ; }
      return vct ;
    }
    // ========================================================================
    /// shift all elements of vector 
    template <class TYPE , class ALLOCATOR, typename SCALAR>
    void shift ( std::vector<TYPE, ALLOCATOR>& vct , SCALAR factor ) 
    { shift    ( vct.begin() , vct.end () , factor ) ; }
    // ========================================================================
    template <class ITERATOR> 
    void negate ( ITERATOR first , ITERATOR last ) 
    { for ( ; first != last ; ++first ) { (*first) = -(*first) ; } }
    // ========================================================================
    template <class TYPE, class ALLOCATOR> 
    void negate ( std::vector<TYPE,ALLOCATOR>& vct ) 
    { negate ( vct.begin() , vct.end() ) ; }
    // ========================================================================
    
    // ========================================================================
    /// check that all elements osf sequence are finite        
    template <class ITERATOR, 
              typename = std::enable_if_t<std::is_arithmetic_v<typename std::iterator_traits<ITERATOR>::value_type>>>
    inline bool isfinite
    ( ITERATOR begin ,
      ITERATOR end   )
    { return std::all_of ( begin ,
                           end   ,
                           [](const auto& val) -> bool
                           { return std::isfinite ( val ) ; } ) ; }
    
    // ========================================================================
    /// sum of all vector elements \f$ \Sum v_i \f$
    template <class ITERATOR>
    inline long double sum
    ( ITERATOR      begin , 
      ITERATOR      end   )
    { return std::reduce ( begin , end , 0.0L ) ; }
    // ========================================================================
    /// sum of all absolute values of vector elements \f$ \Sum \left| v_i  \right| \f$
    template <class ITERATOR>
    inline long double sum1
    ( ITERATOR      begin , 
      ITERATOR      end   )      
    { return std::transform_reduce
        ( begin ,
          end   ,
          0.0L  ,
          std::plus<>() ,
          []( long double v ) { return std::abs ( v ) ; } ) ; }
    // ========================================================================
    /// sum of all squared values of vector elements \f$ \Sum v_i^2 \f$
    template <class ITERATOR>
    inline long double sum2
    ( ITERATOR      begin , 
      ITERATOR      end   )      
    { return std::transform_reduce
        ( begin       ,
          end         ,
          0.0L        ,
          std::plus<>() ,
          []( const long double v ) { return v * v ; } ) ; }
    // ========================================================================
    /// sum of all powered values of vector elements \f$ \Sum \left| v_i\right|^p \f$
    template <class ITERATOR>
    inline long double sum_pow
    ( ITERATOR      begin , 
      ITERATOR      end   ,
      const double  p     ) 
    { return std::transform_reduce
        ( begin       ,
          end         ,
          0.0L        ,
          std::plus<>() ,
          [ p ]( const long double v ) { return std::pow ( std::abs ( v ) , p ) ; } ) ; } 
           
    // ========================================================================
    /** Calculate p-norm for the vector 
     *  \f$ |v|_{p} \equiv = \left( \sum_i  \left| v_i\right|^{p} \right)^{1/p}\f$ 
     *  Few special cases:
     *  - p==1        : sum of absolute values 
     *  - p==infinity : the maximal absolute value 
     *  @param begin begin-itetator for the sequnce of coefficients 
     *  @param end   end-iterator for the sequnce of coefficients 
     *  @param pinv  (1/p)
     *  @return p-norm of the vector
     */
    template <class ITERATOR>
    long double p_norm 
    ( ITERATOR      begin , 
      ITERATOR      end   , 
      const double  pinv  )  //  1/p
    {
      /// check i/p
      const double ip = pinv < 0 ? 0 : pinv > 1 ? 1  : pinv ;
      ///
      long  double r  = 0 ;
      /// few "easy" cases:  treat explicitely
      /// 1) (p==1)        : sum of absolute values 
      if      ( 1 == ip ) { return sum1 ( begin , end ) ; } 
      /// 2) (p==infinity) : maximal  absolute value 
      else if ( 0 == ip )    // p = infinity
      {
        for ( ; begin != end ; ++begin ) 
        {
          const long double c = *begin ;
          r = std::max ( r , std::abs ( c ) ) ; 
        }
        return r ;                                                      // RETURN 
      }
      /// 3) (p==2) : frequent case 
      else if ( 0.5 == ip )  // p = 2 : frequent case 
      { return std::sqrt ( sum2 ( begin , end ) ) ; }
      /// 4) not very large integer 
      else if (  ( 0.05 < ip ) && Ostap::Math::isushort ( 1 / ip ) ) 
      {
        const unsigned short p = Ostap::Math::round ( 1 / ip ) ;
        for ( ; begin != end ; ++begin ) 
        {
          const long double c = *begin ;
          r += Ostap::Math::POW ( std::abs ( c ) , p ) ; 
        }
        return std::pow ( r , ip ) ;                                    // RETURN 
      }
      /// 5) generic case
      const long double p = 1.0L/ip ;
      //
      return std::pow ( sum_pow ( begin , end , p ) , ip ) ;
    }
    // ========================================================================
    /** Calculate p-norm for the vector 
     *  \f$ |v|_{p} \equiv = \left( \sum_i  \left| v_i\right|^{p} \right)^{1/p}\f$ 
     *  @param vct  the vector 
     *  @param pinv  (1/p)
     *  @return p-norm of the vector
     */
    template <class TYPE, class ALLOCATOR>
    long double p_norm 
    ( const std::vector<TYPE, ALLOCATOR>& vct  , 
      const double                        pinv )  //  1/p
    { return p_norm ( vct.begin() , vct.end() , pinv ) ; }
    // ========================================================================
    /** sign of the number 
     *  @see https://stackoverflow.com/a/4609795
     */
    template <typename T> 
    inline constexpr std::int8_t signum ( T x , std::false_type /* is_signed */ ) 
    { return   T ( 0 ) < x ; }    
    template <typename T> 
    inline constexpr std::int8_t signum ( T x , std::true_type  /* is_signed */ ) 
    { return ( T ( 0 ) < x ) - ( x < T ( 0 ) ) ; }
    template <typename T>
    inline constexpr std::int8_t signum ( T x ) 
    { return signum ( x , std::is_signed<T>() ) ; }
    // ========================================================================
    /// positive number ?
    template <typename T>
    inline constexpr bool is_positive ( T x )
    { return 0 < signum ( x ) ; }
    /// negative number ?
    template <typename T>
    inline constexpr bool is_negative ( T x )
    { return 0 > signum ( x ) ; }
    /// numbers of same sign ?
    template <typename T>
    inline constexpr bool same_sign   ( T  x , T y )
    { return 0 < signum ( x ) * signum ( y ) ; }
    // ========================================================================
    /** number of (strickt) sign-variations in the sequence
     *  @param first begin-iterator for the input sequence 
     *  @param last  end-iterator for the  input sequence 
     *  @return number if strickt sign variations 
     */
    template <class ITERATOR, class ZERO>
    std::size_t sign_changes
    ( ITERATOR first , 
      ITERATOR last  , 
      ZERO     zero  )
    {
      while ( first != last && zero ( *first ) ) { ++first ; }
      //
      if ( first == last ) { return 0 ; }         //   RETURN
      //
      signed   int si = signum ( *first ) ;
      std::size_t  nc = 0 ;
      for ( ITERATOR j = first + 1 ; j != last ; ++j )
      {
        if ( zero ( *j ) )  { continue ;  }         // CONTINUE
        const signed int sj  = signum ( *j ) ;
        if ( 0 <= si * sj ) { continue ;  }         // CONTINUE 
        nc +=1  ;
        si = sj ;
      }
      return nc ;
    }
    // ========================================================================
    /** reduce the argument to the desired range ( low , high ) 
     *  @code
     *  double xx = 123.45 ;
     *  double x , n ;
     *  std::tie ( x , n ) = reduce ( xx , -5. , 5. ) ;
     *  @endcode  
     *  - Suitable for reduction of 
     *    periodic functions with period \f$ L = high-low \f4 
     */
    template <class TYPE=double>
    std::pair<TYPE,TYPE>
    reduce
    ( const TYPE x    ,
      const TYPE low  ,
      const TYPE high )
    {
      const TYPE xmin = std::min ( low , high ) ;
      const TYPE xmax = std::max ( low , high ) ;
      const TYPE L    = xmax - xmin ;
      const TYPE N    = std::floor ( ( x - xmin ) / L ) ;
      return std::make_pair ( x - N * L , N )   ;
    }
    // ========================================================================
    
    // ========================================================================
    /// valid permutation of indices? 
    bool        is_valid_permutation ( const std::vector<std::size_t>& indices ) ;      
    /// valid permutation of indices? 
    inline bool valid_permutation    ( const std::vector<std::size_t>& indices )
    { return is_valid_permutation ( indices ) ; } 
    // ========================================================================
    
    // ========================================================================
    template <class T>
    inline const T& min ( const T& a ) { return a ; }
    template <class T>
    inline const T& min ( const T& a , const T& b ) { return std::min ( a, b ) ; }
    template <class T, typename... Ts>
    inline const T& min ( const T& a , const T& b , const Ts&... ts ) 
    { return min ( min ( a , b ) , std::forward<Ts> ( ts )...  ) ; }
    // 
    template <class T>
    inline const T& max ( const T& a ) { return a ; }
    template <class T>
    inline const T& max ( const T& a , const T& b ) { return std::max ( a , b ) ; }
    template <class T, typename... Ts>
    inline const T& max ( const T& a , const T& b , const Ts&... ts ) 
    { return max ( max ( a , b ) , std::forward<Ts> ( ts )...  ) ; }
    //
    template <class T>
    inline std::pair<const T&,const T&> minmax ( const T& a ) 
    { return std::pair<const T&,const T&> ( a , a ) ; }
    template <class T>
    inline std::pair<const T&,const T&> minmax ( const T& a , const T&b) 
    { return a <= b ? std::pair<const T&,const T&> ( a , b ) : std::pair<const T&,const T&> ( b , a ) ; }
    template <class T, typename... Ts>
    inline std::pair<const T&,const T&> minmax ( const T& a , const T&b, const Ts&... ts) 
    { 
      std::pair<const T&,const T&> m { minmax ( b , std::forward<Ts>(ts)... ) } ;
      return std::pair<const T&,const T&> ( min ( a , m.first ) , max ( a , m.second ) ) ;
    }
    //
    template <class T>
    inline T absmin ( const T& a ) { return std::abs ( a ) ; }
    template <class T>
    inline T absmin ( const T& a , const T& b ) { return std::min ( std::abs ( a ) , std::abs ( b ) ) ; }
    template <class T, typename... Ts>
    inline T absmin ( const T& a , const T& b , Ts&&... ts ) 
    { return absmin ( absmin ( a , b ) , std::forward<Ts> ( ts )...  ) ; }
    //
    template <class T>
    inline T absmax ( const T& a ) { return std::abs ( a ) ; }
    template <class T>
    inline T absmax ( const T& a , const T& b ) { return std::max ( std::abs ( a ) , std::abs ( b ) ) ; }
    template <class T, typename... Ts>
    inline T absmax ( const T& a , const T& b , Ts&&... ts ) 
    { return absmax ( absmax ( a , b ) , std::forward<Ts> ( ts )...  ) ; }
    //
    
    // ========================================================================
    // some vector sums 
    // ========================================================================    
    
    // ========================================================================
    /// Nolume of N-ball 
    template <unsigned short>
    class NBallVolume_ ;
    // ========================================================================
    /// Nolume of N-ball 
    template <>
    class NBallVolume_<0>
    {
    public:
      // ======================================================================
      /// volume of 0-ball is 1 independently on the radius        
      static double volume ( const double /* radius */ ) { return unit_volume ; }
      /// volume of 0-ball of unit radius 
      static double volume () { return unit_volume   ; }
      /// volume of unit 0-ball is 1 
      static constexpr double unit_volume = 1 ;
      // ======================================================================
    } ;
    // ========================================================================
    /// Nolume of 1-ball 
    template <>
    class NBallVolume_<1>
    {
    public:
      // ======================================================================
      /// volume of 1-ball is 2*R 
      static double volume ( const double radius )
      { return unit_volume * radius ; }
      /// volume of 1-ball of unit radius 
      static double volume () { return unit_volume   ; }
      /// volume of unit 2-ball is 2  
      static constexpr double unit_volume = 2 ;
      // ======================================================================
    } ;
    // ========================================================================
    /// Volume of n-ball 
    template <unsigned short N>
    class NBallVolume_
    {
    public: 
      // ======================================================================
      /// volume of n-ball
      static double volume ( const double radius  )
      { return unit_volume * std::pow ( radius , N ) ; }
      /// volume of n-ball of unit radius 
      static double volume () { return unit_volume   ; }
      /// volume of unit n-ball is : 2*pi/n*V(n-2) 
      static constexpr double unit_volume =
#if defined ( __cplusplus ) && defined ( __cpp_lib_math_constants ) && ( 201907L <= __cpp_lib_math_constants )
        ( 2 * std::numbers::pi_v<long double>         * NBallVolume_<N-2>::unit_volume ) / N  ;
#else 
        ( 2 * 3.141592653589793238462643383279502884L * NBallVolume_<N-2>::unit_volume ) / N  ;
#endif 
      // ======================================================================
    } ;
    // ========================================================================
    /** Volume of n-ball of given radius r  
     *  @param n dimension
     *  @param r radius 
     */
    double nball_volume
    ( const unsigned short n       ,
      const double         r = 1.0 ) ;
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_MATH_H
// ============================================================================
