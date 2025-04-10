// ============================================================================
#ifndef OSTAP_GOF_H 
#define OSTAP_GOF_H 1
// ============================================================================
// Incldue files
// ============================================================================
// STD STL 
// ============================================================================
#include <iterator>
#include <algorithm>
#include <cmath>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Buffer.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  /** @namespace GoF
   *  Collection of utilities for Goodness-of-Fit estimators 
   *  @see https://doi.org/10.1111/1467-9868.00337
   *  @author Vanya BELYAEV Ivan.BEelyaev@cern.ch 
   */
  namespace GoF
  {
    // ========================================================================
    /** Kolmogorov-Smirnov estimator for input sorted array \f$ F( x_i) \f$ 
     * - (sorted) input array of CDF ( x_i ) 
     *  @see https://doi.org/10.1111/1467-9868.00337
     *  @param begin begin-iterator for input sorted arrat of CDF( x(i) ) 
     *  @param end   end-iterator for input sorted arrat of CDF( x(i) ) 
     *  @return value of Kolmogorov-Smirnov estimator 
     */
    template <class ITERATOR>
    inline double
    kolmogorov_smirnov
    ( ITERATOR begin , 
      ITERATOR end   )
      {
        const std::size_t N = std::distance ( begin , end ) ;
        double result = - std::numeric_limits<double>::max() ;
        for ( std::size_t i = 0 ; begin !=end ; ++i , ++begin )
          {
            const double f = *begin ;
            result = std::max ( result , std::max ( ( i + 1.0 ) / N  - f , f - i * 1.0 / N ) ) ;
          }
        return result ; // * result ;
      } ;
    // ========================================================================
    /** Kuiper estimator for input sorted array \f$ F( x_i) \f$ 
     * - (sorted) input array of CDF ( x_i ) 
     *  @see https://doi.org/10.1111/1467-9868.00337
     *  @param begin begin-iterator for input sorted arrat of CDF( x(i) ) 
     *  @param end   end-iterator for input sorted arrat of CDF( x(i) ) 
     *  @return value of Kuiper estimator 
     */
    template <class ITERATOR>
    inline double
    kuiper
    ( ITERATOR begin , 
      ITERATOR end   )
      {
        const std::size_t N = std::distance ( begin , end ) ;
        double d1 = - std::numeric_limits<double>::max() ;
        double d2 = - std::numeric_limits<double>::max() ;
        for ( std::size_t i = 0 ; begin !=end ; ++i , ++begin )
          {
            const double f = *begin ;
            d1 = std::max ( d1 , ( i + 1.0 ) / N - f ) ;
            d2 = std::max ( d2 , f - i * 1.0 / N     ) ;
          }
        return d1 + d2 ;
      } ;
    // ========================================================================
    /** Anderson-Sarling estimator for input sorted array \f$ F( x_i) \f$ 
     * - (sorted) input array of CDF ( x_i ) 
     *  @see https://doi.org/10.1111/1467-9868.00337
     *  @param begin begin-iterator for input sorted arrat of CDF( x(i) ) 
     *  @param end   end-iterator for input sorted arrat of CDF( x(i) ) 
     *  @return value of Andersen-Darling estimator 
     */
    template <class ITERATOR>
    inline double
    anderson_darling
    ( ITERATOR begin , 
      ITERATOR end   )
      {
        const std::size_t N = std::distance ( begin , end ) ;
        double result = 0 ;
        for ( std::size_t i = 0 ; begin !=end ; ++i, ++begin )
          {
            const double f = *begin ;
            result += ( i + 0.5 ) * std::log ( f ) + ( N - i + 0.5 ) * std::log ( 1.0 - f ) ;
          }
        return  -2 * result / N - N ;
      } ;
    // ========================================================================
    /** Cramer-von Mises estimator for input sorted array \f$ F( x_i) \f$ 
     * - (sorted) input array of CDF ( x_i ) 
     *  @see https://doi.org/10.1111/1467-9868.00337
     *  @param begin begin-iterator for input sorted arrat of CDF( x(i) ) 
     *  @param end   end-iterator for input sorted arrat of CDF( x(i) ) 
     *  @return value of Cramer-von Mises estimator 
     */
    template <class ITERATOR>
    inline double
    cramer_von_mises
    ( ITERATOR begin , 
      ITERATOR end   )
      {
        const std::size_t N = std::distance ( begin , end ) ;
        double result = 0 ;
        for ( std::size_t i = 0 ; begin !=end ; ++i, ++begin )
          {
            const double f = *begin ;
            result += std::pow ( f - ( i + 0.5 ) / N , 2 ) ;
          }
        return result + 1 / ( 12.0 * N ) ;
      } ;
    // ========================================================================
    /** Zhang's Z_K estimator for input sorted array \f$ F( x_i) \f$ 
     * - (sorted) input array of CDF ( x_i ) 
     *  @see https://doi.org/10.1111/1467-9868.00337
     *  @param begin begin-iterator for input sorted arrat of CDF( x(i) ) 
     *  @param end   end-iterator for input sorted arrat of CDF( x(i) ) 
     *  @return value of Zhang's Z_K estimator 
     */
    template <class ITERATOR>
    inline double
    ZK 
    ( ITERATOR begin , 
      ITERATOR end   )
      {
        const std::size_t N = std::distance ( begin , end ) ;
        double result = - std::numeric_limits<double>::max() ;
        for ( std::size_t i = 0 ; begin !=end ; ++i, ++begin )
          {
            const double f  = *begin ;
            const double i1 = i     + 0.5 ;
            const double i2 = N - i - 0.5 ;           
            const double r1 = i2 * std::log ( i1 / ( N *        f   ) ) ;
            const double r2 = i2 * std::log ( i2 / ( N *  ( 1 - f ) ) ) ;            
            result = std::max ( result , r1 + r2 ) ;
          }
        return result ;
      } ;
    // ========================================================================
    /** Zhang's Z_A estimator for input sorted array \f$ F( x_i) \f$ 
     * - (sorted) input array of CDF ( x_i ) 
     *  @see https://doi.org/10.1111/1467-9868.00337
     *  @param begin begin-iterator for input sorted array of CDF( x(i) ) 
     *  @param end   end-iterator for input sorted array of CDF( x(i) ) 
     *  @return value of Zhang's Z_A estimator 
     */
    template <class ITERATOR>
    inline double
    ZA 
    ( ITERATOR begin , 
      ITERATOR end   )
      {
        const std::size_t N = std::distance ( begin , end ) ;
        double result = 0 ;
        for ( std::size_t i = 0 ; begin !=end ; ++i, ++begin )
          {
            const double f  = *begin ;
            const double i1 = i     + 0.5 ;
            const double i2 = N - i - 0.5 ;
            result -= std::log ( f ) / i2 + std::log ( 1 - f ) / i1 ;
          }
        return result ;
      } ;
    // ========================================================================
    /** Zhang's Z_C estimator for input sorted array \f$ F( x_i) \f$ 
     * - (sorted) input array of CDF ( x_i ) 
     *  @see https://doi.org/10.1111/1467-9868.00337
     *  @param begin begin-iterator for input sorted array of CDF( x(i) ) 
     *  @param end   end-iterator for input sorted array of CDF( x(i) ) 
     *  @return value of Zhang's Z_C estimator
     */
    template <class ITERATOR>
    inline double
    ZC 
    ( ITERATOR begin , 
      ITERATOR end   )
      {
        const std::size_t N = std::distance ( begin , end ) ;
        double result = 0 ;
        for ( std::size_t i = 0 ; begin !=end ; ++i, ++begin )
          {
            const double f  = *begin ;
            const double ni = ( N - 0.5 ) / ( i + 0.25 ) - 1.0 ;
            result += std::pow ( std::log ( ( 1.0/f - 1.0 ) / ni ) , 2 ) ;
          }
        return result ;
      } ;
    // ========================================================================

    // ========================================================================
    // overloads with buffers 
    // ========================================================================
    
    // ========================================================================
    /** Kolmogorov-Smirnov estimator for input sorted array \f$ F( x_i) \f$ 
     * - (sorted) input array of CDF ( x_i ) 
     *  @see https://doi.org/10.1111/1467-9868.00337
     *  @param buffer (sorted) input array of CDF ( x_i ) 
     *  @return value of Kolmogorov-Smirnov estimator
     */
    template <class DATA>
    inline double
    kolmogorov_smirnov  
    ( const Ostap::Utils::Buffer<DATA>& buffer )
    { return kolmogorov_smirnov ( buffer.begin() , buffer.end () ) ; }
    // ========================================================================
    /** Kuiper estimator for input sorted array \f$ F( x_i) \f$ 
     * - (sorted) input array of CDF ( x_i ) 
     *  @see https://doi.org/10.1111/1467-9868.00337
     *  @param begin begin-iterator for input sorted arrat of CDF( x(i) ) 
     *  @param end   end-iterator for input sorted arrat of CDF( x(i) ) 
     *  @return value of Kuiper estimator 
     */
    template <class DATA>
    inline double
    kuiper
    ( const Ostap::Utils::Buffer<DATA>& buffer )
    { return kuiper ( buffer.begin() , buffer.end () ) ; }
    // ========================================================================    
    /** Anderson-Sarling estimator for input sorted array \f$ F( x_i) \f$ 
     * - (sorted) input array of CDF ( x_i ) 
     *  @see https://doi.org/10.1111/1467-9868.00337
     *  @param buffer (sorted) input array of CDF ( x_i ) 
     *  @return value of Andersen-Darling estimator
     */
    template <class DATA>
    inline double
    anderson_darling
    ( const Ostap::Utils::Buffer<DATA>& buffer )
    { return anderson_darling ( buffer.begin() , buffer.end () ) ; }
    // =========================================================================
    /** Cramer-von Mises estimator for input sorted array \f$ F( x_i) \f$ 
     * - (sorted) input array of CDF ( x_i ) 
     *  @see https://doi.org/10.1111/1467-9868.00337
     *  @param buffer (sorted) input array of CDF ( x_i ) 
     *  @return value of Cramer-von Mises estimator 
     */
    template <class DATA>
    inline double
    cramer_von_mises
    ( const Ostap::Utils::Buffer<DATA>& buffer )
    { return cramer_von_mises ( buffer.begin() , buffer.end () ) ; }
    // =========================================================================
    /** Zhang's Z_A estimator for input sorted array \f$ F( x_i) \f$ 
     * - (sorted) input array of CDF ( x_i ) 
     *  @see https://doi.org/10.1111/1467-9868.00337
     *  @param buffer (sorted) input array of CDF ( x_i ) 
     *  @return value of Zhang's Z_A estimator 
     */
    template <class DATA>
    inline double
    ZA 
    ( const Ostap::Utils::Buffer<DATA>& buffer )
    { return ZA ( buffer.begin() , buffer.end () ) ; }
    // =========================================================================
    /** Zhang's Z_K estimator for input sorted array \f$ F( x_i) \f$ 
     * - (sorted) input array of CDF ( x_i ) 
     *  @see https://doi.org/10.1111/1467-9868.00337
     *  @param buffer (sorted) input array of CDF ( x_i ) 
     *  @param begin begin-iterator for input sorted arrat of CDF( x(i) ) 
     *  @param end   end-iterator for input sorted arrat of CDF( x(i) ) 
     *  @return value of Zhang's Z_K estimator 
     */
    template <class DATA>
    inline double
    ZK 
    ( const Ostap::Utils::Buffer<DATA>& buffer )
    { return ZK ( buffer.begin() , buffer.end () ) ; }
    // =========================================================================
    /** Zhang's Z_C estimator for input sorted array \f$ F( x_i) \f$ 
     * - (sorted) input array of CDF ( x_i ) 
     *  @see https://doi.org/10.1111/1467-9868.00337
     *  @param buffer (sorted) input array of CDF ( x_i ) 
     *  @return valeu of Zhang's Z_C estimator
     */
    template <class DATA>
    inline double
    ZC 
    ( const Ostap::Utils::Buffer<DATA>& buffer )
    { return ZC ( buffer.begin() , buffer.end () ) ; }
    // ========================================================================
  } //                                          The END of namespace Ostap::GoF 
  // ==========================================================================
} //                                                  The END of namespaceOstap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_GOF_H
// ============================================================================
