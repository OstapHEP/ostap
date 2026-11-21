// ============================================================================
#ifndef OSTAP_QUANTILES_H 
#define OSTAP_QUANTILES_H 1
// ============================================================================
// STD & STL
// ============================================================================
#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <type_traits>
#include <utility>
#include <vector>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/MakeArray.h"
#include "Ostap/Math.h"
#include "Ostap/QuantileTypes.h"
#include "Ostap/WeightedEntry.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /** @class QuantileMixin Ostap/Quantiles.h
     *  @CRTP mixin class providing common quantile evaluation interfaces.
     *
     *  Exposes uniform interfaces over sorted input sequences [first, last)
     *  for:
     *  - a single probability @c p,
     *  - a vector of probabilities @c probs (`std::vector<double>`),
     *  - equidistant @c N-quantiles (`quantiles_<N>` and named wrappers such as
     *    @c minmax, @c median, @c quartiles, @c deciles, @c percentiles).
     *
     *  Derived estimators implement @c quantile(first, last, p) and may
     *  override @c quantiles(first, last, probs) or @c quantiles_<N>(first, last)
     *  when single-pass batch evaluation is advantageous.
     *
     *  @tparam QUANTILE Concrete derived quantile estimator class (CRTP).
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
     *  @date   2025-06-20
     */
    template <class QUANTILE>
    class QuantileMixin
    {
    public:
      // ======================================================================
      /** Evaluate quantile for a single probability @c p on sorted data.
       *
       *  @tparam ITERATOR   Iterator type over sorted values convertible to @c double.
       *  @param  first      Begin iterator of the sorted input sequence.
       *  @param  last       End iterator of the sorted input sequence.
       *  @param  p          Probability level in @f$[0, 1]@f$.
       *  @return Estimated quantile value corresponding to @c p.
       *  @attention Input sequence @c [first, last) is assumed to be sorted.
       */
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true>
      inline double operator()
      ( ITERATOR     first , 
        ITERATOR     last  , 
        const double p     ) const
      {
        const QUANTILE& q = static_cast<const QUANTILE&>( *this ) ;
        return q.quantile ( first , last , p ) ;
      }
      // ======================================================================
      /** Evaluate quantiles for a vector of probabilities on sorted data.
       *
       *  @tparam ITERATOR   Iterator type over sorted values convertible to @c double.
       *  @param  first      Begin iterator of the sorted input sequence.
       *  @param  last       End iterator of the sorted input sequence.
       *  @param  probs      Vector of probability levels.
       *  @return Vector of estimated quantiles matching @c probs order.
       *  @attention Input sequence @c [first, last) is assumed to be sorted.
       */
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true>
      inline std::vector<double> quantiles
      ( ITERATOR                   first ,
        ITERATOR                   last  ,
        const std::vector<double>& probs ) const
      {
        const QUANTILE& q = static_cast<const QUANTILE&>( *this ) ;
        std::vector<double> result ;
        result.reserve ( probs.size () ) ;
        for ( const double p : probs ) { result.push_back ( q.quantile ( first , last , p ) ) ; }
        return result ;
      }
      // ======================================================================
      /** @brief Evaluate @c N+1 equidistant grid quantiles (@c p = k / N, @c k = 0..N).
       *
       *  Includes boundary points corresponding to @c p = 0 and @c p = 1.
       *
       *  @tparam N          Number of equal probability intervals (@c N >= 1).
       *  @tparam ITERATOR   Iterator type over sorted values convertible to @c double.
       *  @param  first      Begin iterator of the sorted input sequence.
       *  @param  last       End iterator of the sorted input sequence.
       *  @return Array of @c N+1 quantile estimates.
       *  @attention Input sequence @c [first, last) is assumed to be sorted.
       */
      template <unsigned int N, class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value&&1<=N,bool>::type = true >
      inline std::array<double,N+1>
      quantiles_
      ( ITERATOR first , 
        ITERATOR last  ) const 
      {
        const auto q = [this,first,last] ( const std::size_t k ) -> double
        { return (*this) ( first , last , k * 1.0 / N ) ; } ;
        return Ostap::Math::make_array ( q , std::make_index_sequence<N+1>() ) ;
      } ;
      // ======================================================================
      /// 1-quantiles (minimum and maximum): returns array of size 2.
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true >
      inline std::array<double,2>
      minmax ( ITERATOR first , ITERATOR last ) const
      { return this->template quantiles_<1> ( first , last ) ; }
      // ======================================================================
      /// 2-quantiles (min, median, max): returns array of size 3.
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true >
      inline std::array<double,3>
      median ( ITERATOR first , ITERATOR last ) const
      { return this->template quantiles_<2> ( first , last ) ; }
      // ======================================================================
      /// 3-quantiles (terciles + boundaries): returns array of size 4.
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true >
      inline std::array<double,4>
      terciles ( ITERATOR first , ITERATOR last ) const
      { return this->template quantiles_<3> ( first , last ) ; }
      // ======================================================================
      /// 4-quantiles (quartiles + boundaries): returns array of size 5.
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true >
      inline std::array<double,5>
      quartiles ( ITERATOR first , ITERATOR last ) const
      { return this->template quantiles_<4> ( first , last ) ; }
      // ======================================================================
      /// 5-quantiles (quintiles + boundaries): returns array of size 6.
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true >
      inline std::array<double,6>
      quintiles ( ITERATOR first , ITERATOR last ) const
      { return this->template quantiles_<5> ( first , last ) ; }
      // ======================================================================
      /// 6-quantiles (sextiles + boundaries): returns array of size 7.
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true >
      inline std::array<double,7>
      sextiles ( ITERATOR first , ITERATOR last ) const
      { return this->template quantiles_<6> ( first , last ) ; }
      // ======================================================================
      /// 7-quantiles (septiles + boundaries): returns array of size 8.
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true >
      inline std::array<double,8>
      septiles ( ITERATOR first , ITERATOR last ) const
      { return this->template quantiles_<7> ( first , last ) ; }
      // ======================================================================
      /// 8-quantiles (octiles + boundaries): returns array of size 9.
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true >
      inline std::array<double,9>
      octiles ( ITERATOR first , ITERATOR last ) const
      { return this->template quantiles_<8> ( first , last ) ; }
      // ======================================================================
      /// 10-quantiles (deciles + boundaries): returns array of size 11.
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true >
      inline std::array<double,11>
      deciles ( ITERATOR first , ITERATOR last ) const
      { return this->template quantiles_<10> ( first , last ) ; }
      // ======================================================================
      /// 20-quantiles (ventiles + boundaries): returns array of size 21.
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true >
      inline std::array<double,21>
      ventiles ( ITERATOR first , ITERATOR last ) const
      { return this->template quantiles_<20> ( first , last ) ; }
      // ======================================================================
      /// 100-quantiles (percentiles + boundaries): returns array of size 101.
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true >
      inline std::array<double,101>
      percentiles ( ITERATOR first , ITERATOR last ) const
      { return this->template quantiles_<100> ( first , last ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class QCheck Ostap/Quantiles.h
     *  Helper checker validating non-empty and (optionally) sorted ranges.
     */
    class QCheck
    {
    public:
      // ====================================================================
      /** @brief Construct checker with optional sorting validation.
       *  @param check If @c true, verifies that input range is sorted.
       */
      QCheck ( const bool check = false ) ;
      // ====================================================================
      /** @brief Validate input range @c [first, last).
       *  @throws Throws if range is empty or if sorting check is enabled and fails.
       */
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true  >
      inline void check
      ( ITERATOR first ,
        ITERATOR last  ) const 
      {
        if ( first == last )
          { this->throw_exception ( "Input data cannot be empty!" , __FILE__ , __LINE__ ) ; }
        if ( this->m_check && !std::is_sorted ( first , last ) )
          { this->throw_exception ( "Input data must be sorted!"  , __FILE__ , __LINE__ ) ; }
      }
      // ====================================================================
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,WeightedEntry>::value
                                        && !std::is_convertible<value_type,double>::value,bool>::type = true>
      inline void check
      ( ITERATOR first ,
        ITERATOR last  ) const
      {
        if ( first == last )
          { this->throw_exception ( "Input data cannot be empty!" , __FILE__ , __LINE__ ) ; }
        if ( this->m_check &&
             !std::is_sorted ( first , last ,
                               [] ( const WeightedEntry& a , const WeightedEntry& b ) -> bool
                               { return a.value () < b.value () ; } ) )
          { this->throw_exception ( "Input data must be sorted!"  , __FILE__ , __LINE__ ) ; }
      }
      // ====================================================================
    private:
      // ====================================================================
      void throw_exception 
      ( const char* message            ,  
        const char* file     = nullptr , 
        const long  line     = -1      ) const ;
      // ====================================================================
    private:
      // ====================================================================
      /// @brief Flag controlling whether @c std::is_sorted is checked.
      bool m_check { false } ;
      // ====================================================================
    } ;
    // ========================================================================
    /** @class HyndmanFan Ostap/Quantiles.h
     *  @brief Hyndman-Fan taxonomy of sample quantile estimators (Types 1..9).
     *  @see https://en.wikipedia.org/wiki/Quantile
     *  @see https://doi.org/10.2307%2F2684934
     */
    class HyndmanFan : public QuantileMixin<HyndmanFan>
    {
    public:
      // ======================================================================
      /** @brief Construct Hyndman-Fan estimator of given type.
       *  @param t     Hyndman-Fan estimator variant (default: Type 8).
       *  @param check Enable range sorting validation when @c true.
       */
      HyndmanFan 
      ( const Ostap::QuantileTypes::HyndmanFanType t = Ostap::QuantileTypes::HyndmanFanType::Eight , 
        const bool check                             = false ) ; 
      // ======================================================================
      /** @brief Compute @c p-quantile for sorted range @c [first, last).
       *  @param first Begin iterator of sorted sequence.
       *  @param last  End iterator of sorted sequence.
       *  @param p     Probability in @f$[0, 1]@f$.
       *  @return Estimated quantile value.
       */
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true  >
      inline double quantile
      ( ITERATOR      first ,
        ITERATOR      last  , 
        const double  p     ) const 
      {
        m_check.check ( first , last ) ;
        if ( 0 >= p ) { return *first ; } 
        //
        const std::size_t N = static_cast<std::size_t>( std::distance ( first , last ) ) ;
        if ( 1 == N ) { return *first ; }
        if ( 1 <= p ) { std::advance ( first , N - 1 ) ; return *first ; }
        //
        const double h = 
          ( Ostap::QuantileTypes::HyndmanFanType::One   == m_t ) ?   N * p                 :
          ( Ostap::QuantileTypes::HyndmanFanType::Two   == m_t ) ?   N * p + 0.5           :
          ( Ostap::QuantileTypes::HyndmanFanType::Three == m_t ) ?   N * p - 0.5           : 
          ( Ostap::QuantileTypes::HyndmanFanType::Four  == m_t ) ?   N * p                 :
          ( Ostap::QuantileTypes::HyndmanFanType::Five  == m_t ) ?   N * p + 0.5           : 
          ( Ostap::QuantileTypes::HyndmanFanType::Six   == m_t ) ?   N * p - 0.5           :
          ( Ostap::QuantileTypes::HyndmanFanType::Seven == m_t ) ? ( N - 1    ) * p + 1    :
          ( Ostap::QuantileTypes::HyndmanFanType::Eight == m_t ) ? ( N + 1./3 ) * p + 1./3 :
          ( N + 0.25 ) * p + 0.375 ;
        //
        // Clamp to valid 0-based index range [0, N - 1].
        const double hh = std::clamp ( h - 1 , 0.0 , N - 1.0 ) ;
        //
        if ( Ostap::Math::islong ( hh ) )
        {
          const std::size_t nn = static_cast<std::size_t> ( hh ) ;
          if ( nn ) { std::advance ( first , nn ) ; }             
          return *first ;
        }
        //
        // Handle non-interpolating / step-wise Hyndman-Fan types (1..3) explicitly.
        if ( Ostap::QuantileTypes::HyndmanFanType::One == m_t )
        {
          const std::size_t nn = Ostap::Math::round_up ( hh ) ;
          if ( nn ) { std::advance ( first , nn ) ; } 
          return *first ;
        }
        else if ( Ostap::QuantileTypes::HyndmanFanType::Two == m_t )
        {
          const std::size_t n1 = Ostap::Math::round_half_down ( hh ) ;
          const std::size_t n2 = Ostap::Math::round_half_up   ( hh ) ;
          if ( n1      ) { std::advance ( first , n1      ) ; } 
          const double v1 = *first ; 
          if ( n1 != n2 ) { std::advance ( first , n2 - n1 ) ; } 
          const double v2 = *first ;
          return 0.5 * ( v1 + v2 ) ;
        } 
        else if ( Ostap::QuantileTypes::HyndmanFanType::Three == m_t )
        {
          const std::size_t nn = Ostap::Math::banker ( hh ) ;
          if ( nn ) { std::advance ( first , nn ) ; }          
          return *first ;
        }
        //
        // Continuous linear interpolation between floor and ceiling order statistics.
        const std::size_t nf = Ostap::Math::round_down ( hh ) ;
        const std::size_t nc = Ostap::Math::round_up   ( hh ) ; 
        //
        if ( nf       ) { std::advance ( first , nf      ) ; } 
        const double vf = *first ;
        if ( nc != nf ) { std::advance ( first , nc - nf ) ; } 
        const double vc = *first ;
        //
        return vf + ( hh - std::floor ( hh ) ) * ( vc - vf ) ;
      }
      // ======================================================================
    private:
      // ======================================================================
      /// Selected Hyndman-Fan quantile estimator type.
      Ostap::QuantileTypes::HyndmanFanType m_t     { Ostap::QuantileTypes::HyndmanFanType::Eight } ;
      /// Range validator.
      QCheck                               m_check { false } ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class ABQuantile Ostap/Quantiles.h
     *  Two-parameter (@c alpha, @c beta) linear interpolation quantile estimator.
     */
    class ABQuantile : public QuantileMixin<ABQuantile>
    {
    public:
      // ====================================================================
      /// Construct from explicit @c alpha and @c beta parameters.
      ABQuantile
      ( const double alpha = 0.4   , 
        const double beta  = 0.4   , 
        const bool   check = false ) ;
      /// Construct from an @c ABQuantileType descriptor.
      ABQuantile
      ( const Ostap::QuantileTypes::ABQuantileType& abq          , 
        const bool                                  check = false ) ;
      // ====================================================================
      /// Compute @c p-quantile for sorted range @c [first, last).
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true  >
      inline double quantile 
      ( ITERATOR      first ,
        ITERATOR      last  , 
        const double  p     ) const 
      {
        m_check.check ( first , last ) ; 
        if ( 0 >= p ) { return *first ; } 
        // 
        const std::size_t N = static_cast<std::size_t>( std::distance ( first , last ) ) ;
        if ( 1 == N ) { return *first ; }
        if ( 1 <= p ) { std::advance ( first , N - 1 ) ; return *first ; } 
        //
        const double m = m_abq.m ( p ) ;
        const double a = p * N + m ;
        const int    j = static_cast<int> ( std::floor ( a ) ) ;
        //
        if      ( 0 > j      ) { return *first ; }
        else if ( N <= j + 1 ) { std::advance ( first , N - 1 ) ; return *first ; }
        //
        const double g = a - j ;
        std::advance ( first , j ) ;
        const double v1 = *first ;
        std::advance ( first , 1 ) ;
        const double v2 = *first ;
        //
        return ( 1 - g ) * v1 + g * v2 ;
      }
      // ====================================================================
      /// Return parameter @c alpha.
      inline double alpha () const { return m_abq.alpha () ; }
      /// Return parameter @c beta.
      inline double beta  () const { return m_abq.beta  () ; }
      // ====================================================================
    private:
      // ====================================================================
      /// Parameter pair (@c alpha, @c beta).
      Ostap::QuantileTypes::ABQuantileType m_abq   {} ;
      /// Range validator.
      QCheck                               m_check { false } ; 
      // ====================================================================
    } ; 
    // ========================================================================
    /** @class HarrellDavis Ostap/Quantiles.h
     *  Harrell-Davis distribution-free quantile estimator for unweighted data.
     *  @see https://doi.org/10.1093/biomet/69.3.635
     *  @attention Computationally intensive for large sample sizes @c N.
     */
    class HarrellDavis : public QuantileMixin<HarrellDavis>
    {
    public:
      // ======================================================================
      /** Construct Harrell-Davis estimator.
       *  @param check Enable range sorting validation when @c true.
       */
      HarrellDavis ( const bool check = false ) ;  
      // ======================================================================
      /** Compute a single @c p-quantile over sorted range @c [first, last).
       *  @param first Begin iterator of sorted sequence.
       *  @param last  End iterator of sorted sequence.
       *  @param p     Probability level.
       *  @return Estimated quantile value.
       */
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true  >
      inline double quantile 
      ( ITERATOR      first ,
        ITERATOR      last  , 
        const double  p     ) const 
      {
        m_check.check ( first , last ) ; 
        if ( 0 >= p ) { return *first ; } 
        // 
        const std::size_t N = static_cast<std::size_t>( std::distance ( first , last ) ) ;
        if ( 1 == N ) { return *first ; }
        if ( 1 <= p ) { std::advance ( first , N - 1 ) ; return *first ; } 
        //
        double result = 0 ;
        for ( std::size_t i = 0 ; first != last ; ++first , ++i )
          {
            const double ti    = ( i + 1.0 ) / N ;
            const double tp    =   i * 1.0   / N ;
            const double value = *first ;
            if ( !value ) { continue ; }
            //
            result += value * WHD ( N , p , ti , tp ) ; 
          }
        return result ;
      }
      // ======================================================================
      /** Compute quantiles for a vector of probabilities in a single pass.
       *
       *  Traverses @c [first, last) once and accumulates contributions for all
       *  requested probabilities @c probs simultaneously.
       *
       *  @param first Begin iterator of sorted sequence.
       *  @param last  End iterator of sorted sequence.
       *  @param probs Vector of probability levels.
       *  @return Vector of estimated quantiles corresponding to @c probs.
       */
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true  >
      inline std::vector<double> quantiles
      ( ITERATOR                   first ,
        ITERATOR                   last  ,
        const std::vector<double>& probs ) const
      {
        std::vector<double> result ( probs.size () , 0.0 ) ;
        this->eval_probs ( first , last , probs.begin () , probs.end () , result.begin () ) ;
        return result ;
      }
      // ======================================================================
      /** Compute @c N+1 equidistant grid quantiles in a single pass.
       *
       *  Traverses @c [first, last) once for @c p = k / N (@c k = 0..N).
       *
       *  @tparam N     Number of equal probability intervals (@c N >= 1).
       *  @param  first Begin iterator of sorted sequence.
       *  @param  last  End iterator of sorted sequence.
       *  @return Array of @c N+1 quantile estimates.
       */
      template <unsigned int N, class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value&&1<=N,bool>::type = true >
      inline std::array<double,N+1>
      quantiles_
      ( ITERATOR first ,
        ITERATOR last  ) const
      {
        std::array<double,N+1> probs  {} ;
        std::array<double,N+1> result {} ;
        for ( std::size_t k = 0 ; k <= N ; ++k )
          { probs [ k ] = k * 1.0 / N ; }
        this->eval_probs ( first , last , probs.begin () , probs.end () , result.begin () ) ;
        return result ;
      }
      // ======================================================================
    private:
      // ======================================================================
      /** Single-pass evaluation helper over @c [first, last) for arbitrary
       *  probability sequence @c [p_first, p_last).
       */
      template <class ITERATOR, class PITER, class OITER>
      inline void eval_probs
      ( ITERATOR first   ,
        ITERATOR last    ,
        PITER    p_first ,
        PITER    p_last  ,
        OITER    out     ) const
      {
        m_check.check ( first , last ) ;
        const std::size_t M = static_cast<std::size_t>( std::distance ( p_first , p_last ) ) ;
        if ( first == last || 0 == M ) { return ; }
        //
        const std::size_t count = static_cast<std::size_t>( std::distance ( first , last ) ) ;
        if ( 1 == count )
        {
          std::fill_n ( out , M , static_cast<double>( *first ) ) ;
          return ;
        }
        //
        std::fill_n ( out , M , 0.0 ) ;
        ITERATOR it = first ;
        for ( std::size_t i = 0 ; it != last ; ++it , ++i )
        {
          const double value = *it ;
          const double ti    = ( i + 1.0 ) / count ;
          const double tp    =   i * 1.0   / count ;
          //
          PITER pit = p_first ;
          OITER oit = out     ;
          for ( ; pit != p_last ; ++pit , ++oit )
          {
            const double p = *pit ;
            if      ( 0.0 >= p ) { if ( 0 == i )         { *oit = value ; } }
            else if ( 1.0 <= p ) { if ( count - 1 == i ) { *oit = value ; } }
            else if ( value )
              { *oit += value * WHD ( count , p , ti , tp ) ; }
          }
        }
      }
      // ======================================================================
    private:
      // ======================================================================
      /// Range validator.
      QCheck m_check { false } ;
      // ======================================================================
      /** Compute @f$ I_{t_1}(\alpha,\beta) - I_{t_2}(\alpha,\beta) @f$
       *  using normalized incomplete beta function with caching.
       */
      static double WHD 
      ( const std::size_t N  , 
        const double      p  , 
        const double      t1 ,
        const double      t2 ) ; 
      // ======================================================================
    } ;
    // ========================================================================
    /** @class WHarrellDavis Ostap/Quantiles.h
     *  Weighted Harrell-Davis quantile estimator over @c WeightedEntry ranges.
     *  @see https://doi.org/10.1093/biomet/69.3.635
     *  @attention Computationally intensive for large datasets.
     */
    class WHarrellDavis : public QuantileMixin<WHarrellDavis>
    {
    public:
      // ======================================================================
      /// The actual Entry type storing value and weight (`WeightedEntry`).
      typedef WeightedEntry Entry ;
      // ======================================================================
      /** Construct weighted Harrell-Davis estimator.
       *  @param check Enable range sorting validation when @c true.
       */
      WHarrellDavis ( const bool check = false ) ;  
      // ======================================================================
      /** Compute single @c p-quantile using precalculated @c sumw and @c sumw2.
       *  @param first Begin iterator of sorted weighted sequence.
       *  @param last  End iterator of sorted weighted sequence.
       *  @param p     Probability level in @f$(0, 1)@f$.
       *  @param sumw  Precalculated sum of weights (must be positive).
       *  @param sumw2 Precalculated sum of squared weights (must be positive).
       *  @return Weighted Harrell-Davis quantile estimate.
       */
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,Entry>::value,bool>::type = true>
      inline double quantile 
      ( ITERATOR     first ,
        ITERATOR     last  , 
        const double p     ,
        const double sumw  ,
        const double sumw2 ) const
      {
        // Inverse total weight and effective sample size n* = (sumw)^2 / sumw2.
        const double sw_inv = 1.0 / sumw ;
        const double nstar  = sumw * sumw / sumw2 ;
        const double alpha  = ( nstar + 1.0 ) *       p       ;
        const double beta   = ( nstar + 1.0 ) * ( 1.0 - p ) ;        
        //
        double wsum   = 0.0 ;
        double result = 0.0 ; 
        for ( ; first != last ; ++first )
        {
          const Entry  entry  { *first } ;
          const double value  = entry.value  () ;
          const double weight = entry.weight () ;
          if ( !weight ) { continue ; } 
          //
          const double tp = std::min ( std::max ( sw_inv * wsum , 0.0 ) , 1.0 ) ;            
          wsum += weight ;
          const double ti = std::min ( std::max ( sw_inv * wsum , 0.0 ) , 1.0 ) ;            
          //
          result += value * WHD ( alpha , beta , ti , tp ) ; 
        }
        return result ; 
      }
      // ======================================================================
      /** Compute single @c p-quantile, computing @c sumw and @c sumw2 first.
       *  @param first Begin iterator of sorted weighted sequence.
       *  @param last  End iterator of sorted weighted sequence.
       *  @param p     Probability level in @f$(0, 1)@f$.
       *  @return Weighted Harrell-Davis quantile estimate.
       */
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,Entry>::value,bool>::type = true>
      inline double quantile 
      ( ITERATOR     first ,
        ITERATOR     last  , 
        const double p     ) const 
      {
        double sumw  = 0.0 ;
        double sumw2 = 0.0 ;        
        this->sum_weights ( first , last , sumw , sumw2 ) ;
        return this->quantile ( first , last , p , sumw , sumw2 ) ;
      }
      // ======================================================================
      /** Compute quantiles for @c probs using precalculated @c sumw and @c sumw2
       *         in a single pass over @c [first, last).
       *  @param first Begin iterator of sorted weighted sequence.
       *  @param last  End iterator of sorted weighted sequence.
       *  @param probs Vector of probability levels.
       *  @param sumw  Precalculated sum of weights (must be positive).
       *  @param sumw2 Precalculated sum of squared weights (must be positive).
       *  @return Vector of quantile estimates corresponding to @c probs.
       */
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,Entry>::value,bool>::type = true>
      inline std::vector<double> quantiles
      ( ITERATOR                   first ,
        ITERATOR                   last  ,
        const std::vector<double>& probs ,
        const double               sumw  ,
        const double               sumw2 ) const
      {
        std::vector<double> result ( probs.size () , 0.0 ) ;
        this->eval_probs ( first , last , probs.begin () , probs.end () , result.begin () , sumw , sumw2 ) ;
        return result ;
      }
      // ======================================================================
      /** Compute quantiles for @c probs, calculating @c sumw and @c sumw2 first.
       *  @param first Begin iterator of sorted weighted sequence.
       *  @param last  End iterator of sorted weighted sequence.
       *  @param probs Vector of probability levels.
       *  @return Vector of quantile estimates corresponding to @c probs.
       */
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,Entry>::value,bool>::type = true>
      inline std::vector<double> quantiles
      ( ITERATOR                   first ,
        ITERATOR                   last  ,
        const std::vector<double>& probs ) const
      {
        double sumw  = 0.0 ;
        double sumw2 = 0.0 ;
        this->sum_weights ( first , last , sumw , sumw2 ) ;
        return this->quantiles ( first , last , probs , sumw , sumw2 ) ;
      }
      // ======================================================================
      /** Compute @c N+1 equidistant grid quantiles using precalculated
       *         @c sumw and @c sumw2 in a single pass over @c [first, last).
       *  @tparam N     Number of equal probability intervals (@c N >= 1).
       *  @param  first Begin iterator of sorted weighted sequence.
       *  @param  last  End iterator of sorted weighted sequence.
       *  @param  sumw  Precalculated sum of weights (must be positive).
       *  @param  sumw2 Precalculated sum of squared weights (must be positive).
       *  @return Array of @c N+1 quantile estimates.
       */
      template <unsigned int N, class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,Entry>::value&&1<=N,bool>::type = true >
      inline std::array<double,N+1>
      quantiles_
      ( ITERATOR     first ,
        ITERATOR     last  ,
        const double sumw  ,
        const double sumw2 ) const
      {
        std::array<double,N+1> probs  {} ;
        std::array<double,N+1> result {} ;
        for ( std::size_t k = 0 ; k <= N ; ++k ) { probs [ k ] = k * 1.0 / N ; }
        this->eval_probs ( first , last , probs.begin () , probs.end () , result.begin () , sumw , sumw2 ) ;
        return result ;
      }
      // ======================================================================
      /** Compute @c N+1 equidistant grid quantiles, calculating @c sumw
       *         and @c sumw2 first.
       *  @tparam N     Number of equal probability intervals (@c N >= 1).
       *  @param  first Begin iterator of sorted weighted sequence.
       *  @param  last  End iterator of sorted weighted sequence.
       *  @return Array of @c N+1 quantile estimates.
       */
      template <unsigned int N, class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,Entry>::value&&1<=N,bool>::type = true >
      inline std::array<double,N+1>
      quantiles_
      ( ITERATOR first ,
        ITERATOR last  ) const
      {
        double sumw  = 0.0 ;
        double sumw2 = 0.0 ;
        this->sum_weights ( first , last , sumw , sumw2 ) ;
        return this->template quantiles_<N> ( first , last , sumw , sumw2 ) ;
      }
      // ======================================================================
    private:
      // ======================================================================
      /** Accumulate sum of weights (@c sumw) and sum of squared weights (@c sumw2)
       *  over @c [first, last).
       */
      template <class ITERATOR>
      inline void sum_weights
      ( ITERATOR first ,
        ITERATOR last  ,
        double&  sumw  ,
        double&  sumw2 ) const
      {
        sumw  = 0.0 ;
        sumw2 = 0.0 ;
        std::for_each ( first ,
                        last  ,
                        [&sumw,&sumw2] ( const Entry& entry ) -> void
                        {
                          const double weight = entry.weight () ;
                          if ( !weight ) { return ; } 
                          sumw  +=          weight ;
                          sumw2 += weight * weight ;
                        } ) ;        
      }
      // ======================================================================
      /** Single-pass weighted Harrell-Davis evaluation over @c [first, last)
       *  for arbitrary probability sequence @c [p_first, p_last).
       */
      template <class ITERATOR, class PITER, class OITER>
      inline void eval_probs
      ( ITERATOR     first   ,
        ITERATOR     last    ,
        PITER        p_first ,
        PITER        p_last  ,
        OITER        out     ,
        const double sumw    ,
        const double sumw2   ) const
      {
        m_check.check ( first , last ) ;
        const std::size_t M = static_cast<std::size_t>( std::distance ( p_first , p_last ) ) ;
        if ( first == last || 0 == M ) { return ; }
        //
        const double sw_inv = 1.0 / sumw ;
        const double nstar  = sumw * sumw / sumw2 ;
        //
        // Precompute Beta parameters (alpha, beta) for each requested probability.
        std::vector<double> alphas ( M , 0.0 ) ;
        std::vector<double> betas  ( M , 0.0 ) ;
        std::size_t idx = 0 ;
        for ( PITER pit = p_first ; pit != p_last ; ++pit , ++idx )
        {
          const double p = *pit ;
            alphas [ idx ] = ( nstar + 1.0 ) *       p       ;
            betas  [ idx ] = ( nstar + 1.0 ) * ( 1.0 - p ) ;
        }
        //
        std::fill_n ( out , M , 0.0 ) ;
        double wsum = 0.0 ;
        for ( ; first != last ; ++first )
        {
          const Entry  entry  { *first } ;
          const double value  = entry.value  () ;
          const double weight = entry.weight () ;
          if ( !weight ) { continue ; }
          //
          const double tp = std::min ( std::max ( sw_inv * wsum , 0.0 ) , 1.0 ) ;
          wsum += weight ;
          const double ti = std::min ( std::max ( sw_inv * wsum , 0.0 ) , 1.0 ) ;
          //
          PITER pit = p_first ;
          OITER oit = out     ;
          for ( std::size_t m = 0 ; pit != p_last ; ++pit , ++oit , ++m )
          {
            const double p = *pit ;
            if      ( 0.0 >= p ) { if ( 0.0 == wsum ) { *oit = value ; } }
            else if ( 1.0 <= p ) { *oit = value ; }
            else                 { *oit += value * WHD ( alphas [ m ] , betas [ m ] , ti , tp ) ; }
          }
        }
      }
      // ======================================================================
    private:
      // ======================================================================
      /// Range validator.
      QCheck m_check { false } ;
      // ======================================================================
      /** @brief Compute @f$ I_{t_1}(\alpha,\beta) - I_{t_2}(\alpha,\beta) @f$
       *         using normalized incomplete beta function with caching.
       */
      static double WHD 
      ( const double alpha , 
        const double beta  , 
        const double t1    ,
        const double t2    ) ; 
      // ======================================================================
    } ;
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                The  end of namespace Ostap 
// ============================================================================
#endif // OSTAP_QUANTILES_H
// ============================================================================
//                                                                      The END 
// ============================================================================
