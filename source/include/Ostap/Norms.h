// ============================================================================
#ifndef OSTAP_NORMS_H
#define OSTAP_NORMS_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL 
// ============================================================================
#include <cstddef>     
#include <cmath>       
#include <iterator>    
#include <type_traits> 
#include <limits>      
#include <algorithm>
#include <numeric>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Math.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    
    // ========================================================================
    // Basic vector norms
    // ========================================================================

    /** @brief Compute L0 "norm" (count of non-zero elements) in a sequence 
     *  ||v||_0 = count(|v_i| > eps)
     *  @param begin (INPUT) begin-iterator 
     *  @param end   (INPUT) end-iterator 
     *  @param eps   (INPUT) Absolute tolerance threshold below which elements are considered zero
     *  @return Number of non-zero elements
     */
    template <typename ITERATOR , 
              typename VALUE    = typename std::iterator_traits<ITERATOR>::value_type            ,
              typename          = std::enable_if_t<std::is_convertible_v<VALUE, double>> >
    inline std::size_t norm_L0
    ( ITERATOR     begin ,
      ITERATOR     end   , 
      const double eps   = std::numeric_limits<double>::epsilon () )
      { return std::count_if
        ( begin ,
          end   , 
          [ eps ] ( const double value ) -> bool 
          { return std::isfinite ( value ) && value && ( eps < std::abs ( value ) ) ; } ) ; }

   // =======================================================================
    /** @brief Compute L1 norm (Manhattan norm / sum of absolute values)
     *  ||v||_1 = sum(|v_i|)
     *  @param begin (INPUT) begin-iterator 
     *  @param end   (INPUT) end-iterator 
     *  @return L1 norm value
     */
    template <typename ITERATOR , 
              typename VALUE    = typename std::iterator_traits<ITERATOR>::value_type            ,
              typename          = std::enable_if_t<std::is_convertible_v<VALUE, double>> >
    inline auto norm_L1
    ( ITERATOR     begin ,
      ITERATOR     end   ) 
    {
      // Automatically choose return type: preserve long double precision if requested, default to double otherwise
      using ReturnType = std::conditional_t<std::is_same_v<std::decay_t<VALUE>, long double>, long double, double>;
      return std::transform_reduce
        ( begin         ,
          end           ,
          ReturnType{0} ,
          std::plus<>() ,
          []( const auto x ) -> ReturnType 
          { return static_cast<ReturnType> ( std::abs ( x ) ) ; } ) ;
    }

    
    // ========================================================================
    /** @brief Compute fast L2 norm 
     *  ||v||_2 = sqrt(v . v)
     *  @param[in] v Input vector
     *  @param begin (INPUT) begin-iterator 
     *  @param end   (INPUT) end-iterator 
     *  @return L_2 norm value
     */
    template <typename ITERATOR , 
              typename VALUE    = typename std::iterator_traits<ITERATOR>::value_type    ,
              typename          = std::enable_if_t<std::is_convertible_v<VALUE, double>> >
    inline auto norm_L2
    ( ITERATOR     begin ,
      ITERATOR     end   ) 
    {
      // Automatically choose return type: preserve long double precision if requested, default to double otherwise
      using ReturnType = std::conditional_t<std::is_same_v<std::decay_t<VALUE>, long double>, long double, double>;
      //
      return std::sqrt
      ( std::transform_reduce
        ( begin         ,
          end           ,
          ReturnType{0} ,
          std::plus<>() ,
          []( const auto x )
          { const auto val = static_cast<ReturnType> ( x ) ; return val * val ; } ) );
    }

    // ==========================================================================
    /** @brief Compute stable L2 norm (Euclidean norm) avoiding overflow/underflow
     *  @param begin (INPUT) begin-iterator 
     *  @param end   (INPUT) end-iterator 
     *  @return L_2 norm value
     */
    template <typename ITERATOR , 
              typename VALUE    = typename std::iterator_traits<ITERATOR>::value_type    ,
              typename          = std::enable_if_t<std::is_convertible_v<VALUE, double>> >
    inline auto norm_L2_safe 
    ( ITERATOR     begin ,
      ITERATOR     end   ) 
    {
      //
      using ReturnType = std::conditional_t<std::is_same_v<std::decay_t<VALUE>, long double>, long double, double>;
      if ( begin == end) { return ReturnType(0); }
      //
      // Step 1: Find the maximum absolute value in the sequence to use as a scaling factor
      ReturnType scale = 0;
      for ( auto it = begin; it != end; ++it)
      { scale = std::max ( scale , static_cast<ReturnType> ( std::abs ( *it ) ) ) ; }

      // If all elements are zero (or vector is effectively zero)
      if ( 0 == scale ) { return ReturnType{0} ; } 
    
      // Step 2: Sum the squares of scaled elements: sum((v_i / scale)^2)
      ReturnType sum_sq = std::transform_reduce
      ( begin         ,
        end           ,
        ReturnType{0} ,
       std::plus<>()  ,
       [ scale ]( const auto x )
        { const auto val = static_cast<ReturnType>(x) / scale; return val * val; } ) ;

      // Step 3: Scale back the result
      return scale * std::sqrt ( sum_sq ) ;
    }

    // ========================================================================
    /** @brief Compute L_infinity norm (Chebyshev norm / maximum absolute element)
     *  ||v||_inf = max(|v_i|)
     *  @param begin (INPUT) begin-iterator 
     *  @param end   (INPUT) end-iterator 
     *  @return L_inf norm value: maxmimm absolute value
     */
    template <typename ITERATOR , 
              typename VALUE    = typename std::iterator_traits<ITERATOR>::value_type            ,
              typename          = std::enable_if_t<std::is_convertible_v<VALUE, double>> >
    inline auto norm_Linf
    ( ITERATOR     begin ,
      ITERATOR     end   ) 
    {
      // Automatically choose return type: preserve long double precision if requested, default to double otherwise
      using ReturnType = std::conditional_t<std::is_same_v<std::decay_t<VALUE>, long double>, long double, double>;
      return std::transform_reduce
      ( begin         ,
        end           ,
        ReturnType{0} ,
        []( auto a , auto b ) { return std::max ( a , b ) ; } ,
        []( const auto x    ) { return static_cast<ReturnType> ( std::abs ( x ) ) ; } ) ;
    }

    // ========================================================================
    /** @brief Compute L_infinity norm (Chebyshev norm / maximum absolute element)
     *  ||v||_inf = max(|v_i|)
     *  @param begin (INPUT) begin-iterator 
     *  @param end   (INPUT) end-iterator 
     *  @return L_inf norm value: maxmimm absolute value
     */
    template <typename ITERATOR , 
              typename VALUE    = typename std::iterator_traits<ITERATOR>::value_type            ,
              typename          = std::enable_if_t<std::is_convertible_v<VALUE, double>> >
    inline auto norm_max
    ( ITERATOR     begin ,
      ITERATOR     end   ) 
     { return norm_Linf ( begin , end ) ; } 

    
    // ========================================================================
    /** @brief Compute generalized Lp-norm
     *  ||v||_p = (sum(|v_i|^p))^(1/p)
     *  @param begin (INPUT) begin-iterator 
     *  @param end   (INPUT) end-iterator 
     *  @param p     (INPUT) Power parameter \f$ p \ge 0 \f$
     *  @return Lp norm value
     */
    template <typename ITERATOR , 
              typename VALUE    = typename std::iterator_traits<ITERATOR>::value_type            ,
              typename          = std::enable_if_t<std::is_convertible_v<VALUE, double>> >
    inline auto norm_Lp
    ( ITERATOR     begin ,
      ITERATOR     end   ,
      const double p     ) 
    {
      using ReturnType = std::conditional_t<std::is_same_v<std::decay_t<VALUE>, long double>, long double, double>;

      if ( begin == end ) { return ReturnType { 0 } ; }

      static const Ostap::Math::Equal_To<double> s_equal {} ;
      static const Ostap::Math::Zero    <double> s_zero  {} ;
      
      if      ( std::isinf ( p )            ) { return norm_Linf    ( begin , end ) ; }
      else if ( 1 == p || s_equal ( 1 , p ) ) { return norm_L1      ( begin , end ) ; }
      else if ( 2 == p || s_equal ( 2 , p ) ) { return norm_L2_safe ( begin , end ) ; }
      else if ( 0 >= p || s_zero  ( p     ) ) { return static_cast<ReturnType> ( norm_L0 ( begin , end , std::numeric_limits<double>::epsilon () ) ) ; }
        
      // Step 1: Find the maximum absolute value to scale elements and prevent overflow/underflow
      ReturnType scale = 0;
      for ( auto it = begin; it != end; ++it)
      { scale = std::max ( scale , static_cast<ReturnType> ( std::abs ( *it ) ) ) ; }

      if ( 0 == scale ) { return ReturnType { 0 } ; }

      // Step 2: Sum the p-th powers of scaled elements: sum((|v_i| / scale)^p)
      ReturnType sum_p = std::transform_reduce
      ( begin          ,
        end            ,
        ReturnType {0} ,
        std::plus<>()  ,
        [p, scale]( const auto x )
        { const auto val = static_cast<ReturnType>(std::abs(x)) / scale ; return std::pow ( val , p ) ; } ) ;

      // Step 3: Scale back and apply 1/p power
      return scale * std::pow ( sum_p , ReturnType{1} / p ) ;
    }
        
    // ========================================================================
  } //                                                    end of namespace Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
#endif // OSTAP_NORMS_H
// ============================================================================
//                                                                      The END 
// ============================================================================
