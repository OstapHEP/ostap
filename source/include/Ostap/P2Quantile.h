// ============================================================================
#ifndef OSTAP_P2QUANTILE_H 
#define OSTAP_P2QUANTILE_H 1
// ============================================================================
// STD&STL
// ============================================================================
#include <utility>
#include <algorithm>
#include <array>
// ============================================================================
//  GSL
// ============================================================================
#include "gsl/gsl_rstat.h"
// ============================================================================
// Ostap
// =============================================================================
#include "Ostap/Statistic.h"
// =============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace  Math
  {
    // ========================================================================
    /** @class QantileP2 
     *  P2 algorithm for (approximate quantile estimamtion)
     *  @see https://aakinshin.net/posts/p2-quantile-estimator-intro/
     */
    class  QuantileP2 : public Ostap::Math::Statistic
    {
      public :
        // =======================================================================
        /** constructor from  p
         *  @parameter p (INPUT) p: \f$ 0 \le p \le 1 \f$  
         */  
        QuantileP2 ( const double p  = 0.5) ;
        // ======================================================================
      public :
        // ======================================================================
        /// add one more measurable, update quantiles    
        QuantileP2& add ( const double v ) ;
        /// add group of measurements 
        // =====================================================================
        // Generic counter interface 
        void update ( const double value ) override { add  ( value ) ; }
        // ======================================================================
      public:
        // ======================================================================
        /// sample size 
        inline std::size_t N    () const { return m_N  ; }
        /// sample size 
        inline std::size_t size () const { return N () ; }
        // ======================================================================
      private:
        // ======================================================================
        /// quantile 
        double                    m_p   { 0.5 } ; // quantile 
        /// sample size 
        std::size_t               m_N   { 0   } ;  // sample size 
        std::array<double,5>      m_q   {     } ;
        std::array<std::size_t,5> m_n   {     } ;
        std::array<double,5>      m_np  {     } ;
        std::array<double,5>      m_dnp {     } ;
      // =======================================================================
    };
    // ========================================================================
    namespace GSL
    {
      // ======================================================================
      /** class P2Quantile Ostap/P2Quantile.h
       *  Get running (approximate) quantile using P2-algorithm 
       *  @see https://dl.acm.org/doi/10.1145/4372.4378
       */
      class P2Quantile
      {                                         \
      public:
        // ====================================================================
        /** Standard constructor for quantile
         *  @param p quatile \f$ 0 < p > 1 \f$
         */
        P2Quantile  ( const double p ) ;
        // ===================================================================
        // copy contructor 
        P2Quantile  ( const P2Quantile&  right ) ;
        // move constructor
        P2Quantile  (       P2Quantile&& right ) ;
        // no default constructor 
        P2Quantile  () = delete ;
        // destructor
        ~P2Quantile () ;
        // ====================================================================
        // copy assignement
        P2Quantile& operator=( const P2Quantile&  right ) ;
        // move assignement
        P2Quantile& operator=(       P2Quantile&& right ) ;
        // ====================================================================
        // add single measurement 
        void add ( const double x )
        {
          if ( m_ws == nullptr ) { m_ws = gsl_rstat_quantile_alloc  ( m_p ) ; }
          gsl_rstat_quantile_add ( x , m_ws ) ;
        }
        // add everal measurements
        template <class ITERATOR>
        void add ( ITERATOR begin ,
                   ITERATOR end   )
        {
          if ( m_ws == nullptr ) { m_ws = gsl_rstat_quantile_alloc  ( m_p ) ; }
          for ( ; begin != end ; ++begin )
          { gsl_rstat_quantile_add ( *begin , m_ws ) ; }
        }
        // ====================================================================
        /// get number of measurements
        std::size_t  n () const { return m_ws == nullptr ? 0 : m_ws->n ; }
        // ====================================================================
        // get the quantile value
        double value() const ;
        // ====================================================================
        // get the quantile value via conversion operator 
        operator double () const { return value () ; }
        // ====================================================================
        /// get tht quantile definition
        inline double p ()  const { return m_p ; }
        // ====================================================================
        // swapping two quantile counters 
        void swap ( P2Quantile& right )
        {
          std::swap  ( m_ws  , right.m_ws ) ;
          std::swap  ( m_p   , right.m_p  ) ;        
        }
        // ====================================================================
      private:
        // ====================================================================
        //  GSL workspace  
        gsl_rstat_quantile_workspace*  m_ws { nullptr } ;
        //  nominal quantile
        double m_p ;
        // ====================================================================
      } ;    
      // ======================================================================
      /// swap two quantile counters 
      inline void awap ( P2Quantile& a , P2Quantile& b ) { a.swap ( b ) ; }
      // ======================================================================
    } //                                 The end of namespace  Ostap::Math::GSL
    // ========================================================================
  } //                                         The end od namesapce Ostap::Math
  // ==========================================================================
} //                                                The  end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_P2QUANTILE_H
// ============================================================================
