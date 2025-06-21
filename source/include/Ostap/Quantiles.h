// ============================================================================
#ifndef OSTAP_QUANTILES_H 
#define OSTAP_QUANTILES_H 1
// ============================================================================
// STD&STL
// ============================================================================
#include <algorithm>
#include <cmath>
#include <array>
// ============================================================================
// Ostap
// =============================================================================
#include "Ostap/Math.h"
#include "Ostap/MakeArray.h"
// =============================================================================
namespace Ostap
{
  // ===========================================================================
  namespace  Math
  {
    // =========================================================================
    /** get array of N-quantiles  (+min&max values)
     *  @retur array of lenth N+1 with N-quanties
     *  - first element is minaml value 
     *  - laselement if maximal values 
     */
    template <unsigned short N, class QUANTILE, class ITERATOR,
              typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
              typename = std::enable_if<std::is_convertible<value_type,double>::value&&1<=N> >
    inline std::array<double,N+1>
    quantiles_   
    ( const QUANTILE& quantile ,  
      ITERATOR        first     ,
      ITERATOR        last      )
    {
        auto q = [&quantile,first,last] ( unsigned short k ) -> double
        { return quantile ( first , last , k * 1.0 / ( N + 1 ) ) ; } ;
        return Ostap::Math::make_array ( q , std::make_index_sequence<N+1>() )  ;
    }
    // ========================================================================
    /** get array of N-quantiles  (+min&max values)
     *  @retur array of lenth N+1 with N-quanties
     *  - first element is minaml value 
     *  - laselement if maximal values 
     */
    template <unsigned short N, class QUANTILE, class CONTAINER,
              typename value_type = typename CONTAINER::value_type,
              typename = std::enable_if<std::is_convertible<value_type,double>::value&&1<=N> >
    inline std::array<double,N+1>
    quantiles_   
    ( const QUANTILE&   quantile  ,  
      const CONTAINER&  container )
    {
      return quantiles_<N> ( quantile , container.begin () , container.end() ) ;
    } 
    // ========================================================================
    /** @class QBase
     *  Helper base class for quantile estimators 
     */
    class QBase
    {
      public:
        // ====================================================================
        QBase ( const bool check = false ) ;
        // ====================================================================
      protected :
        // ====================================================================
        template <class ITERATOR,
              typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
              typename = std::enable_if<std::is_convertible<value_type,double>::value> >
      inline void check
      ( ITERATOR first ,
        ITERATOR last  ) const 
        {
          if  ( first == last ) 
          { this->throw_exception ( "Input data cannot be empty!" , __FILE__ , __LINE__ ) ; }
          if ( this->m_check && !std::is_sorted ( first , last ) )
          { this->throw_exception ( "Input data must be sorted!"  , __FILE__ , __LINE__ ) ; }
        }
        // ====================================================================
      protected :
        // ====================================================================
        void throw_exception 
        ( const char* message            ,  
          const char* file     = nullptr , 
          const long  line     = -1      ) const ;
        // ====================================================================
      private:
        // ====================================================================
        /// chech data ?
        bool m_check { false } ; //  check data 
        // ====================================================================
    };
    // ========================================================================
    /** @class HyndmanFan
      * Taxomomy of quantile estimators
      * @see https://en.wikipedia.org/wiki/Quantile
      * @see https://doi.org/10.2307%2F2684934
      */
     class HyndmanFan : protected QBase 
     {
      //========================================================================
      public:
      // =======================================================================
      enum QuantileType {
        One = 1 ,
        Two     ,
        Three   ,
        Four    ,
        Five    ,
        Six     ,
        Seven   ,
        Eight   ,
        Nine    , 
      } ;
      // =======================================================================
      public:
      // =======================================================================
      /// constructor
      HyndmanFan 
        ( const QuantileType t     = Eight , 
          const bool         check = false ) ; 
      // =======================================================================
      public:
      // =======================================================================
      /** get quanitile 
       *  @param firts (INPUT) begin-iteratir for the input sequence
       *  @param last  (INPUT) end-iterator for the innut sequence 
       *  @param p     (INPUT) probability   0<p<1 
       *  @return p-quantile
       *  @attention  the sequence is assumed to be sorted 
       *  @attention  exception if throw for the empty sequence
       */
      template <class ITERATOR,
              typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
              typename = std::enable_if<std::is_convertible<value_type,double>::value> >
      inline double quantile  
      ( ITERATOR      first       ,
        ITERATOR      last        , 
        const double  p     = 0.5 ) const 
        {
          // check input 
          this -> check ( first , last ) ;
          if  ( 0 >= p ) { return *first      ; } 
          // 
          const std::size_t  N  = std::distance ( first , last ) ;
          if  ( 1 == N ) { return *first ; }
          if  ( 1 <= p ) { std::advance ( first , N - 1 ) ; return *first ;} 
          //   
         const double h = 
            ( One   == m_t ) ?   N * p                  :
            ( Two   == m_t ) ?   N * p + 0.5            :
            ( Three == m_t ) ?   N * p - 0.5            : 
            ( Four  == m_t ) ?   N * p                  :
            ( Five  == m_t ) ?   N * p + 0.5            : 
            ( Six   == m_t ) ?   N * p - 0.5            :
            ( Seven == m_t ) ? ( N - 1    ) * p + 1     :
            ( Eight == m_t ) ? ( N + 1./3 ) * p + 1./3  :
            ( N + 0.25 ) * p + 0.375 ;

          /// clumped and adjusted for zero-indexed  
          const double hh = std::clamp ( h - 1 , 0.0 , N - 1.0 ) ;

          if ( Ostap::Math::islong  ( hh ) )
            {
              const std::size_t nn = static_cast<std::size_t> ( hh ) ;
              if ( nn ) { std::advance ( first , nn ) ; } 
              return *first ;
            }
          //
          // first three cases explicitley 
          //
          if ( One == m_t )
          {
            const std::size_t nn = Ostap::Math::round_up ( hh ) ;
            std::advance ( first , nn ) ;
            return *first ; 
          }
          else if ( Two == m_t )
          {
            const std::size_t n1 = Ostap::Math::round_half_down ( hh ) ;
            const std::size_t n2 = Ostap::Math::round_half_up   ( hh ) ;
            if ( n1       ) { std::advance ( first , n1      ) ; } 
            const double v1 = *first ; 
            if ( n1 != n2 ) { std::advance ( first , n2 - n1 ) ; } 
            const double v2 = *first ; 
            return 0.5 * ( v1 + v2 ) ;
          } 
          else if( Three == m_t )
          {
            const std::size_t nn = Ostap::Math::banker ( hh ) ;
            if ( nn ) { std::advance ( first ,  nn ) ; } 
            return *first ;
          }
          //
          const std::size_t nf = Ostap::Math::round_down ( hh ) ;
          const std::size_t nc = Ostap::Math::round_up   ( hh ) ; 
          //
          if ( nf        ) { std::advance ( first , nf      ) ; } 
          const double vf = *first ;
          if ( nc != nf  ) { std::advance ( first , nc - nf ) ; } 
          const double vc = *first ;
          //
          return vf + ( hh - std::floor ( hh ) ) * ( vc - vf ) ; 
        }
      // ======================================================================
      private: 
      // ======================================================================
      /// Algorithm to use 
      QuantileType m_t     { Eight } ; // algorithm to use
      /// check input sequence 
      bool         m_check { false } ; // check input sequence  
      // ======================================================================
     } ;
    // ========================================================================
   /** @class Quanitle2
    *  Variant of Hyndman-Fan with lnear innterplaton using two parameteers
    *  - alpha  \f$ 0 \le \alpha \le 1 \f$ 
    *  - beta   \f$ 0 \le \beta  \le 1 \f$ 
    *
    * Typical values for alphap, betap are:
    *
    * - (0,1) : p(k) = k/n : linear interpolation of cdf (R type 4)
    * - (.5,.5) : p(k) = (k - 1/2.)/n : piecewise linear function (R type 5)
    * - (0,0) : p(k) = k/(n+1) : (R type 6)
    * - (1,1) : p(k) = (k-1)/(n-1): p(k) = mode[F(x[k])]. (R type 7, R default)
    * - (1/3,1/3): p(k) = (k-1/3)/(n+1/3): Then p(k) ~ median[F(x[k])].
    *   The resulting quantile estimates are approximately median-unbiased
    *   regardless of the distribution of x. (R type 8)
    * - (3/8,3/8): p(k) = (k-3/8)/(n+1/4): Blom.
    *   The resulting quantile estimates are approximately
    *   unbiased if x is normally distributed (R type 9)
    * - (0.4,0.4) : approximately quantile unbiased (Cunnane)
    * - (.35,.35): APL, used with PWM
    */  
   class ABQuantile  : protected QBase
    {
      public:
        // ====================================================================
        ABQuantile
          ( const double alpha = 0.4   , 
            const double beta  = 0.4   , 
            const bool   check = false ) ;
        // ====================================================================
      public:
        // ====================================================================
        template <class ITERATOR,
              typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
              typename = std::enable_if<std::is_convertible<value_type,double>::value> >
        inline double quantile  
        ( ITERATOR      first       ,
          ITERATOR      last        , 
          const double  p     = 0.5 ) const 
        {
          this->check ( first , last ) ; 
          if  ( 0 >= p ) { return *first      ; } 
          // 
          const std::size_t  N  = std::distance ( first , last ) ;
          if  ( 1 == N ) { return *first ; }
          if  ( 1 <= p ) { std::advance ( first , N - 1 ) ; return *first ;} 
          //
          const double m = m_alpha + p * ( 1 - m_alpha - m_beta ) ;
          const double a = p * N  + m ;
          const int    j = static_cast<int> ( std::floor ( a ) ) ;
          //
          if      ( 0 > j      ) { return *first ; }
          else if ( N <= j + 1 ) { std::advance ( first , N - 1 ) ; return *first ; }
  
          const double          g = a - j ;
          std::advance ( first , j ) ;
          const double v1 = *first ;
          std::advance ( first , 1 ) ;
          const double v2 = *first ;
          //
          return ( 1 - g ) * v1  + g * v2  ;
        }
        // ====================================================================
        private:
        // ====================================================================
        /// parameter alpha
        double m_alpha { 0.4 } ;
        /// parameter beta 
        double m_beta  { 0.4 } ;
        // ====================================================================
    } ; 
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                The  end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_QUANTILES_H
// ============================================================================
