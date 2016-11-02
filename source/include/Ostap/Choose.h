// $Id$
// ============================================================================ 
#ifndef OSTAP_CHOOSE_H 
#define OSTAP_CHOOSE_H 1
// ============================================================================
// Include files
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {    
    // ========================================================================
    /** calculate the binomial coefficient C(n,k) = n!/((n-k)!*k!)
     *  the result is exact for all n,k<=67
     *  @warning In case of overflow std::numeric_limits<unsigned long long>::max is returned 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-03-08
     */
    unsigned long long  choose ( const unsigned short n ,
                                 const unsigned short k ) ;
    // ========================================================================
    /** calculate the logarithm of binomial coefficient
     *  \f$ \log C^n_k \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-03-08
     */
    double log_choose ( const unsigned short n ,
                        const unsigned short k ) ;
    // ========================================================================
    /** calculate the binomial coefficient C(k,n) = n!/((n-k)!*k!)
     *  @author Vanya BELYAEV Ivan.Belyaev@irep.ru
     *  @date 2015-03-08
     */
    double choose_double       ( const unsigned short n , 
                                 const unsigned short k ) ;
    // ========================================================================
    /** calculate the generalized binomial coefficient C(a,k) 
     *  \f$C(\alpha,k) = \frac{\alpha}{k}\frac{\alpha-1}{k-1}...\f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-03-08
     */
    double gen_choose ( const double         a ,
                        const unsigned short k ) ;
    // ========================================================================
    /** calculate the generalized binomial coefficient C(n/2,k) 
     *  \f$C(n,k) = \frac{n/2}{k}\frac{n/2-1}{k-1}...\f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-03-08
     */
    double choose_half ( const int            n ,
                         const unsigned short k ) ;
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Gaudi
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_CHOOSE_H
// ============================================================================
