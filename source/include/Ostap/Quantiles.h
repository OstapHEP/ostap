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
#include "Ostap/QuantileTypes.h"
// =============================================================================
namespace Ostap
{
  // ===========================================================================
  namespace  Math
  {
    // =========================================================================
    /** @class QuanitleMixin
     *  Mixin class to evaluate varuous quantiles 
     *  @author Vanya BELYAEV Iban.Belyaev@cern.ch
     *  @date 2025-06-20
     */
    template <class QUANTILE>
    class QuantileMixin
    {
    public:
      // ==========================================================================
      /** the main method: calculate quantile for SORTED data
       *  @param p     (INPUT) probability    \f$ 0 \le p \le 1 \f$ 
       *  @param first (INPUT begin-iterator for SORTED data 
       *  @param last  (INPUT end-iterator for SORTED data 
       *  @reutrn quantile value correspondg to probability p 
       *  @attention  data s assumed to be SORTED!
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
      // ==========================================================================
      /** get N-quantiles 
       *  - p=0 and p=1 2quantiels are included!
       *  @param first (INPUT begin-iterator for SORTED data 
       *  @param last  (INPUT end-iterator for SORTED data 
       *  @return N-quantiles
       *  @attention  data s assumed to be SORTED!
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
        { return (*this) ( first , last , k * 1.0 / N  ) ; } ;
        return Ostap::Math::make_array ( q , std::make_index_sequence<N+1>() )  ;
      } ;
      // 1-quantiles:   min, max
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true >
      inline std::array<double,2>
      minmax 
      ( ITERATOR first , 
        ITERATOR last  ) const
      { return this->template quantiles_<1> ( first  , last ) ;	 }
      // 2-quantiles: min,median & max
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true >
      inline std::array<double,3>
      median
      ( ITERATOR first , 
        ITERATOR last  ) const
      { return this->template quantiles_<2> ( first  , last ) ;}
      /// 3-quantiles: min,q1,q2,max 
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true  >
      inline std::array<double,4>
      terciles
      ( ITERATOR first , 
        ITERATOR last  ) const
      { return this->template quantiles_<3> ( first  , last ) ;}
      /// 4-quantiles: min,q1,q2,q3,max
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true  >
      inline std::array<double,5>
      quartiles
      ( ITERATOR first , 
        ITERATOR last  ) const
      { return this->template quantiles_<4> ( first  , last ) ;}
      /// 5-quantiles:
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true  >
      inline std::array<double,6>
      quintiles 
      ( ITERATOR first , 
        ITERATOR last  ) const
      { return this->template quantiles_<5> ( first  , last ) ;}
      /// 6-quantiles:
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true  >
      inline std::array<double,7>
      sextiles 
      ( ITERATOR first , 
        ITERATOR last  ) const
      { return this->template quantiles_<6> ( first  , last ) ;}
      /// 7-quantiles:
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true  >
      inline std::array<double,8>
      septiles 
      ( ITERATOR first , 
        ITERATOR last  ) const
      { return this->template quantiles_<7> ( first  , last ) ;}
      /// 8-quantiles:
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true  >
      inline std::array<double,9>
      octiles 
      ( ITERATOR first , 
        ITERATOR last  ) const
      { return this->template quantiles_<8> ( first  , last ) ;}
      /// 10-quantiles:
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true  >
      inline std::array<double,11>
      deciles 
      ( ITERATOR first , 
        ITERATOR last  ) const
      { return this->template quantiles_<10> ( first  , last ) ;}
      /// 20-quantiles:
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true  >
      inline std::array<double,21>
      ventiles 
      ( ITERATOR first , 
        ITERATOR last  ) const
      { return this->template quantiles_<20> ( first  , last ) ;}
      /// 100-quantiles:
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true  >
      inline std::array<double,101>
      percentiles  
      ( ITERATOR first , 
        ITERATOR last  ) const
      { return this->template quantiles_<100> ( first  , last ) ;}
      // ======================================================================
      template <unsigned short N,typename = std::enable_if<(1<=N)> >
      inline std::array<double,N+1> 
      quantiles_ ( const std::vector<double>& data ) const
      { return this -> template quantiles_<N> ( data.begin(), data.end () ) ; }
      // ======================================================================
      inline std::array<double,2>
      minmax  ( const std::vector<double>& data ) const
      { return this -> minmax ( data.begin() , data.end () ) ; }
      // ======================================================================
      inline std::array<double,3>
      median ( const std::vector<double>& data ) const
      { return this -> median ( data.begin() , data.end () ) ; }
      // ======================================================================
      inline std::array<double,4>
      terciles ( const std::vector<double>& data ) const
      { return this -> terciles  ( data.begin() , data.end () ) ; }
      // ======================================================================
      inline std::array<double,5>
      quartiles ( const std::vector<double>& data ) const
      { return this -> quartiles ( data.begin() , data.end () ) ; }
      // ======================================================================
      inline std::array<double,6>
      quintiles ( const std::vector<double>& data ) const
      { return this -> quintiles ( data.begin() , data.end () ) ; }
      // ======================================================================
      inline std::array<double,7>
      sextiles ( const std::vector<double>& data ) const
      { return this -> sextiles ( data.begin() , data.end () ) ; }
      // ======================================================================
      inline std::array<double,8>
      septiles ( const std::vector<double>& data ) const
      { return this -> septiles ( data.begin() , data.end () ) ; }
      // ======================================================================
      inline std::array<double,9>
      octiles ( const std::vector<double>& data ) const
      { return this -> octiles ( data.begin() , data.end () ) ; }
      // ======================================================================
      inline std::array<double,11>
      deciles ( const std::vector<double>& data ) const
      { return this -> deciles ( data.begin() , data.end () ) ; }
      // ======================================================================
      inline std::array<double,21>
      ventiles ( const std::vector<double>& data ) const
      { return this -> ventiles ( data.begin() , data.end () ) ; }
      // ======================================================================
      inline std::array<double,101>
      percentiles ( const std::vector<double>& data ) const
      { return this -> percentiles ( data.begin() , data.end () ) ; }
      // ======================================================================
      inline double operator()
      ( const std::vector<double>& data , 
        const double               p    ) const
      { return (*this) ( data.begin () , data.end() , p ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class QCheck
     *  Helper checker for quantile estimators 
     */
    class QCheck
    {
    public:
      // ====================================================================
      QCheck ( const bool check = false ) ;
      // ====================================================================
    public : 
      // ====================================================================
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true  >
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
    private : 
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
    class HyndmanFan : public QuantileMixin<HyndmanFan>
    {
      // =======================================================================
    public:      
      // =======================================================================
      /// constructor
      HyndmanFan 
      ( const Ostap::QuantileTypes::HyndmanFanType t =
        Ostap::QuantileTypes::HyndmanFanType::Eight , 
        const bool check = false ) ; 
      // =======================================================================
    public:
      // =======================================================================
      /** get quantiile 
       *  @param firts (INPUT) begin-iteratir for the input sequence
       *  @param last  (INPUT) end-iterator for the innut sequence 
       *  @param p     (INPUT) probability   0<p<1 
       *  @return p-quantile
       *  @attention  the sequence is assumed to be sorted 
       *  @attention  exception if throw for the empty sequence
       */
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true  >
      inline double quantile
      ( ITERATOR      first ,
        ITERATOR      last  , 
        const double  p     ) const 
      {
        // check input 
        m_check.check ( first , last ) ;
        if  ( 0 >= p ) { return *first ; } 
        // 
        const std::size_t  N  = std::distance ( first , last ) ;
        // 
        if  ( 1 == N ) { return *first ; }
        //
        if  ( 1 <= p ) { std::advance ( first , N - 1 ) ; return *first ; }
        //
        const double h = 
          ( Ostap::QuantileTypes::HyndmanFanType::One   == m_t ) ?   N * p                  :
          ( Ostap::QuantileTypes::HyndmanFanType::Two   == m_t ) ?   N * p + 0.5            :
          ( Ostap::QuantileTypes::HyndmanFanType::Three == m_t ) ?   N * p - 0.5            : 
          ( Ostap::QuantileTypes::HyndmanFanType::Four  == m_t ) ?   N * p                  :
          ( Ostap::QuantileTypes::HyndmanFanType::Five  == m_t ) ?   N * p + 0.5            : 
          ( Ostap::QuantileTypes::HyndmanFanType::Six   == m_t ) ?   N * p - 0.5            :
          ( Ostap::QuantileTypes::HyndmanFanType::Seven == m_t ) ? ( N - 1    ) * p + 1     :
          ( Ostap::QuantileTypes::HyndmanFanType::Eight == m_t ) ? ( N + 1./3 ) * p + 1./3  :
          ( N + 0.25 ) * p + 0.375 ;
        //
        /// clamped and adjusted for zero-indexed  
        const double hh = std::clamp ( h - 1 , 0.0 , N - 1.0 ) ;
        // 
        if ( Ostap::Math::islong  ( hh ) )
          {
            const std::size_t nn = static_cast<std::size_t> ( hh ) ;
            if ( nn ) { std::advance ( first , nn ) ; }	    
            return *first ;
          }
        //
        // first three cases explicitley 
        //
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
            if ( n1       ) { std::advance ( first , n1      ) ; } 
            const double v1 = *first ; 
            if ( n1 != n2 ) { std::advance ( first , n2 - n1 ) ; } 
            const double v2 = *first ;
            //
            return 0.5 * ( v1 + v2 ) ;
          } 
        else if ( Ostap::QuantileTypes::HyndmanFanType::Three == m_t )
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
    private :
      // ======================================================================
      /// quantile type
      Ostap::QuantileTypes::HyndmanFanType  m_t { Ostap::QuantileTypes::HyndmanFanType:: Eight } ;
      /// check ?
      QCheck       m_check { false } ; // check 
      // ======================================================================
    };
    // ========================================================================
    /** @class ABQuantile
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
    class ABQuantile: public QuantileMixin<ABQuantile>
    {
    public:
      // ====================================================================
      ABQuantile
      ( const double alpha = 0.4   , 
        const double beta  = 0.4   , 
        const bool   check = false ) ;
      // ====================================================================
      ABQuantile
      ( const Ostap::QuantileTypes::ABQuantileType& abq          , 
        const bool                                 check = false ) ;
      // ====================================================================
    public:
      // ====================================================================
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true  >
      inline double quantile 
      ( ITERATOR      first ,
        ITERATOR      last  , 
        const double  p     ) const 
      {
        m_check.check ( first , last ) ; 
        if  ( 0 >= p ) { return *first      ; } 
        // 
        const std::size_t  N  = std::distance ( first , last ) ;
        if  ( 1 == N ) { return *first ; }
        if  ( 1 <= p ) { std::advance ( first , N - 1 ) ; return *first ;} 
        //
        // const double m = alpha ()  + p * ( 1 - alpha () - beta () ) ;
        const double m = m_abq.m ( p ) ;
        //
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
    public:
      // ====================================================================
      /// get alpha 
      inline double alpha () const { return m_abq.alpha () ; }
      /// get beta 
      inline double beta  () const { return m_abq.beta  () ; }
      // ====================================================================
    private:
      // ====================================================================
      /// alpha/beta keeper 
      Ostap::QuantileTypes::ABQuantileType m_abq   {} ;
      /// check 
      QCheck m_check { false } ; 
      // ====================================================================
    } ; 
    // ========================================================================
    /** @class HarrellDavis 
     *  Harrell-Davis quatile etimator 
     *  @attention it cna e CPU espensive for large data sets
     */
    class HarrellDavis : public QuantileMixin<HarrellDavis>
    {
      // ======================================================================
    public:
      // ======================================================================
      HarrellDavis ( const bool check = false ) ;  
      // ======================================================================
    public:
      // ======================================================================
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,double>::value,bool>::type = true  >
      inline double quantile 
      ( ITERATOR      first ,
        ITERATOR      last  , 
        const double  p     ) const 
      {
        m_check.check ( first , last ) ; 
        if  ( 0 >= p ) { return *first      ; } 
        // 
        const std::size_t  N  = std::distance ( first , last ) ;
        if  ( 1 == N ) { return *first ; }
        if  ( 1 <= p ) { std::advance ( first , N - 1 ) ; return *first ;} 
        //
        double result = 0 ;
        for ( std::size_t i = 0 ; first != last ; ++first , ++i )
          {
            const double ti    = ( i + 1.0 ) / N ;
            const double tp    =   i * 1.0   / N ;
            const double value = *first ;
            if ( !value ) { continue ; }
            //
            result += value * WHD ( N  , p , ti , tp ) ; 
          }
        return result ;
      } ;  
      // ======================================================================
    private :
      // ======================================================================
      /// check ? 
      QCheck  m_check { false } ; // check ? 
      // ======================================================================
    private:
      // ======================================================================
      /** calculate 
       *  \f$ \ I_{t_1}(\alpha,\beta) - I_{t_2} ( \alpha, \beta) f$, where 
       *  \f$ I_z(x,y) \f$ is normalized inncomplete beta function   
       *  @see Ostap::Math::beta_inc  
       *  - protection is added 
       *  - caching is applied 
       */
      static double WHD 
      ( const std::size_t N  , 
        const double      p  , 
        const double      t1 ,
        const double      t2 ) ; 
      // ======================================================================
    } ; //                        The end of the class  Ostap::Math::Harrelavis
    // ========================================================================     
    /** @class WHarrellDavis 
     *  Harrell-Davis quatile  etimator for weighted datas 
     *  @attention it can be  CPU espensive for large data sets
     */
    class WHarrellDavis : public QuantileMixin<WHarrellDavis>
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor 
      WHarrellDavis ( const bool check = false ) ;  
      // ======================================================================
    public:
      // ======================================================================
      struct Entry : public std::pair <double,double> 
      {
        typedef std::pair<double,double> Pair ;
        /// constructor from the value and weight 
        Entry
        ( const double value  = 0 ,
          const double weight = 1 )
          : Pair ( value , weight ) {} 
        /// constructor from the pair         
        Entry ( Pair entry  ) : Pair( entry ){}
      } ;
      // ======================================================================
    public : 
      // ======================================================================
      /** calculate the quantile using the preecalculated sum of weights and sum of squared weigths          
       *  @param first (INPUT) start of (ordered) data-sequence 
       *  @param last  (INPUT) end of   (ordered) data-sequence 
       *  @param p     (input) probability 0 < p < 1 
       *  @param sumw  (input) sum of weights (must be positive) 
       *  @param sumw2 (input) sum of squared weights (must be positive) 
       *  @return Harrell-Davis quantile 
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
        /// inverse total weight 
        const double sw_inv = 1/sumw ;
        /// the number of effectibe entries 
        const double nstar  = sumw * sumw / sumw2 ;
        /// effective alphas & betas
        const double alpha  = ( nstar + 1 ) *       p   ;
        const double beta   = ( nstar + 1 ) * ( 1 - p ) ;        
        //
        /// loop over all (ordered) entries 
        double wsum   =  0 ;
        double result =  0 ; 
        for ( ; first != last ; ++first )
          {
            const Entry  entry { *first } ;
            const double value  = entry.first   ;
            const double weight = entry.second  ;
            //
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
      /** calculate the quantile 
       *  @param first (INPUT) start of (ordered) data-sequence 
       *  @param last  (INPUT) end of   (ordered) data-sequence 
       *  @param p     (input) probability 0 < p < 1 
       *  @return Harrell-Davis quantile 
       */
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,Entry>::value,bool>::type = true>
      inline double quantile 
      ( ITERATOR     first ,
        ITERATOR     last  , 
        const double p     ) const 
      {
        /// (1) calculate  sum of weights and sum of squred weights
        double sumw  = 0 ;
        double sumw2 = 0 ;        
        std::for_each ( first ,
                        last  ,
                        [&sumw,&sumw2] ( Entry entry ) -> void
                        {
                          const double weight = entry.second ;
                          if ( !weight ) { return ; } 
                          sumw  +=          weight ;
                          sumw2 += weight * weight ;
                        } ) ;        
        /// (2) calculate the quantile
        return this->quantile ( first , last , p , sumw , sumw2 ) ;
      }
      // ======================================================================
    private :
      // ======================================================================
      /// check ? 
      QCheck  m_check { false } ; // check ? 
      // ======================================================================
    private:
      // ======================================================================
      /** calculate 
       *  \f$ \ I_{t_1}(\alpha,\beta) - I_{t_2} ( \alpha, \beta) f$, where 
       *  \f$ I_z(x,y) \f$ is normalized inncomplete beta function   
       *  @see Ostap::Math::beta_inc  
       *  - protection is added 
       *  - caching is applied 
       */
      static double WHD 
      ( const double alpha , 
        const double beta  , 
        const double t1    ,
        const double t2    ) ; 
      // ======================================================================
    } ; //                      The end of the class  Ostap::Math::WHarrelavis
    // ========================================================================             
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                The  end of namespace Ostap 
// ============================================================================
#endif // OSTAP_QUANTILES_H
// ============================================================================
//                                                                      The END 
// ============================================================================
