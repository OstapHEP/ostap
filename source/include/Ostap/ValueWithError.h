// ============================================================================
#ifndef OSTAP_VALUEWITHERRORS_H 
#define OSTAP_VALUEWITHERRORS_H 1
// ============================================================================
// Include files
// ============================================================================
// STD& STL 
// ============================================================================
#include <iosfwd>
#include <utility>
#include <vector>
#include <string>
#include <cmath>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatusCode.h"
// ============================================================================
/** @file Ostap/ValueWithError.h
 *  Collection of useful objects with associated "covariances".
 *  The concept has been stollen from Wouter Hulsbergen's lines 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2009-06-03
 */
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    /** @class ValueWithError 
     *  The most simple representation of "value with associated error"
     *  The concept has been stollen from Wouter Hulsbergen's lines 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 20090603
     */
    class ValueWithError 
    {
    public:
      // ======================================================================
      typedef double                                               Value      ;
      typedef double                                               Covariance ;
      // ======================================================================
    public:
      // ======================================================================
      typedef std::vector<ValueWithError>                          Vector     ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the value and covariance 
      ValueWithError ( const double value      = 0 , 
                       const double covariance = 0 ) ;
      // ======================================================================
      /** constructor from the (value,error)-pair 
       *   - first  element is "value" 
       *   - second element is "error" 
       *  @param p  (value,error)-pair 
       */
      explicit ValueWithError ( const std::pair<double,double>& p ) ;
      // ======================================================================
    public: // trivial accessors 
      // ======================================================================
      /// get the value 
      double value      () const { return m_value ; }
      /// get the covariance 
      double cov2       () const { return m_cov2  ; }
      /// get the covariance 
      double covariance () const { return m_cov2  ; }
      /** get the error 
       *  @attention negative erorr is returned for invalid covariance 
       *  @return the error estimate 
       */
      double error      () const ;
      // ======================================================================
    public: // setters 
      // ======================================================================
      void setValue      ( const double v ) {  m_value = v   ; }
      void setCov2       ( const double c ) {  m_cov2  = c   ; }
      void setCovariance ( const double c ) { setCov2  ( c ) ; }
      /** set the error 
       *  @attention invalid covariance is set for negative error 
       */
      void setError      ( const double e ) ; 
      // ======================================================================
    public: // finally it is just a value 
      // ======================================================================
      /// cast 
      operator double () const { return value () ; }                    // cast 
      // ======================================================================
    public:  // operators
      // ======================================================================      
      /// += 
      ValueWithError& operator+= ( const ValueWithError& right ) ;        // += 
      /// -= 
      ValueWithError& operator-= ( const ValueWithError& right ) ;        // -= 
      /// *= 
      ValueWithError& operator*= ( const ValueWithError& right ) ;        // *= 
      /// /= 
      ValueWithError& operator/= ( const ValueWithError& right ) ;        // /= 
      // ======================================================================
    public: 
      // ======================================================================
      /// += 
      ValueWithError& operator+= ( const double v ) ;                     // += 
      /// -= 
      ValueWithError& operator-= ( const double v ) ;                     // -= 
      /// *= 
      ValueWithError& operator*= ( const double v ) ;                     // *= 
      /// /= 
      ValueWithError& operator/= ( const double v ) ;                     // /= 
      // ======================================================================
    public: // unary - 
      // ======================================================================
      /// unary - 
      ValueWithError operator-() const ;                             // unary - 
      // ======================================================================
    public: // as pair:
      // ======================================================================
      /// get as std::pair 
      std::pair<double,double> asPair () const 
      { return std::pair<double,double> ( value () , error () ) ; }
      /// conversion to std::pair
      operator std::pair<double,double>() const 
      { return std::pair<double,double> ( value () , error () ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the mean value 
      ValueWithError mean ( const ValueWithError& right ) const ;
      /// get chi2 distance 
      double chi2 ( const ValueWithError& right ) const ;
      /// get chi2 distance 
      double chi2 ( const double          right ) const ;
      /** get ``residual''
       *  defined as  signed \f$\sqrt \chi^2 \f$
       */
      double residual ( const ValueWithError& right ) const ;
      /** get ``residual''
       *  defined as  signed \f$\sqrt \chi^2 \f$
       */
      double residual ( const double          right ) const ;
      /** get Kullback-Liebler divergency 
       *  @return KL-divergency for valid arguments, -1 otherwise
       */
      double kullback   ( const ValueWithError& right ) const ;
      /** get (squared) Hellinger distance
       *  @see https://en.wikipedia.org/wiki/Hellinger_distance
       *  @return heilinger distance for valid arguments, -1 otherwise
       */
      double hellinger2 ( const ValueWithError& right ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /** evaluate the "fraction" \f$  \frac{a}{a+b} \f$ 
       *  @param  b the parameter "b" for the fraction 
       *  @return a/(a+b) 
       */
      ValueWithError frac ( const ValueWithError& b ) const ;
      /** evaluate the "fraction" \f$  \frac{a}{a+b} \f$ 
       *  @param  b the parameter "b" for the fraction 
       *  @return a/(a+b) 
       */
      ValueWithError frac ( const double          b ) const ;
      /** evaluate the "asymmetry" \f$  \frac{a-b}{a+b} \f$ 
       *  @param  b the parameter "b" for the fraction 
       *  @return (a-b)/(a+b) 
       */
      ValueWithError asym ( const ValueWithError& b ) const ;
      /** evaluate the "asymmetry" \f$  \frac{a-b}{a+b} \f$ 
       *  @param  b the parameter "b" for the fraction 
       *  @return (a-b)/(a+b) 
       */
      ValueWithError asym ( const double          b ) const ;
      // ======================================================================
    public: // NaN and Finite 
      // ======================================================================
      ///  finite ? 
      bool isfinite () const ;      
      ///  normal ? 
      bool isnormal () const ;      
      ///  check for NaN
      bool isnan    () const ;
      ///  check for inf
      bool isinf    () const ;
      // ======================================================================
    public: // good ?
      // ======================================================================
      /// check for goodness: finite values and non-negative covariance 
      bool isgood   () const ;
      /// check for goodness: finite values and non-negative covariance 
      bool good     () const { return isgood   () ; }
      // ======================================================================
    public: // helper functions for Python:
      // ======================================================================
      ///    a + right 
      ValueWithError __add__       ( const ValueWithError& right ) const ;
      ///    a + right 
      ValueWithError __add__       ( const double          right ) const ;
      ///        right + a  
      ValueWithError __radd__      ( const double          right ) const { return __add__ ( right ) ; }
      ///    a - right  
      ValueWithError __sub__       ( const ValueWithError& right ) const ;
      ///    a - right  
      ValueWithError __sub__       ( const double          right ) const ;
      ///        right - a   
      ValueWithError __rsub__      ( const double          right ) const ;
      ///    a * right  
      ValueWithError __mul__       ( const ValueWithError& right ) const ;
      ///    a * right  
      ValueWithError __mul__       ( const double          right ) const ;
      ///        right * a   
      ValueWithError __rmul__      ( const double          right ) const { return __mul__ ( right ) ; }
      ///    a / right  
      ValueWithError __truediv__   ( const ValueWithError& right ) const ;
      ///    a / right  
      ValueWithError __truediv__   ( const double          right ) const ;
      ///        right / a   
      ValueWithError __rtruediv__  ( const double          right ) const ;
      ///    a / right  
      ValueWithError __div__       ( const ValueWithError& right ) const { return __truediv__  ( right ) ; }
      ///    a / right  
      ValueWithError __div__       ( const double          right ) const { return __truediv__  ( right ) ; }
      ///        right / a   
      ValueWithError __rdiv__      ( const double          right ) const { return __rtruediv__ ( right ) ; }
      ///
      ValueWithError& __iadd__     ( const ValueWithError& right ) { (*this) += right ; return *this ; }
      ValueWithError& __imul__     ( const ValueWithError& right ) { (*this) *= right ; return *this ; }
      ValueWithError& __isub__     ( const ValueWithError& right ) { (*this) -= right ; return *this ; }
      ValueWithError& __itruediv__ ( const ValueWithError& right ) { (*this) /= right ; return *this ; }
      ValueWithError& __iadd__     ( const double          right ) { (*this) += right ; return *this ; }
      ValueWithError& __imul__     ( const double          right ) { (*this) *= right ; return *this ; }
      ValueWithError& __isub__     ( const double          right ) { (*this) -= right ; return *this ; }
      ValueWithError& __itruediv__ ( const double          right ) { (*this) /= right ; return *this ; }
      //
      ValueWithError& __idiv__     ( const ValueWithError& right ) { return __itruediv__ ( right ) ; }
      ValueWithError& __idiv__     ( const double          right ) { return __itruediv__ ( right ) ; }
      ///
      ///  abs ( a )   
      ValueWithError __abs__    () const ;
      /// -me
      ValueWithError __neg__    () const ;
      /// +me  (no effect) 
      ValueWithError __pos__    () const ; 
      /// me**e 
      ValueWithError __pow__    ( const int             e ) const ;
      /// me**e 
      ValueWithError __pow__    ( const double          e ) const ;
      /// me**e 
      ValueWithError __pow__    ( const ValueWithError& e ) const ;
      /// e**me 
      ValueWithError __rpow__   ( const int             e ) const ;
      /// e**me 
      ValueWithError __rpow__   ( const double          e ) const ;
      /// exp(me) 
      ValueWithError __exp__    () const ;
      /// exp2(me) 
      ValueWithError __exp2__   () const ;
      /// expm1(me) 
      ValueWithError __expm1__  () const ;
      /// log(me) 
      ValueWithError __log__    () const ;
      /// log2(me) 
      ValueWithError __log2__   () const ;
      /// log10(me) 
      ValueWithError __log10__  () const ;
      /// log1p(me) 
      ValueWithError __log1p__  () const ;
      /// sqrt(me) 
      ValueWithError __sqrt__   () const ;
      /// sqrt(me) 
      ValueWithError __cbrt__   () const ;
      /// sin(me) 
      ValueWithError __sin__    () const ;
      /// cos(me) 
      ValueWithError __cos__    () const ;
      /// tan(me) 
      ValueWithError __tan__    () const ;
      /// sinh(me) 
      ValueWithError __sinh__   () const ;
      /// cosh(me) 
      ValueWithError __cosh__   () const ;
      /// tanh(me) 
      ValueWithError __tanh__   () const ;
      /// erf 
      ValueWithError __erf__    () const ;
      /// erfc
      ValueWithError __erfc__   () const ;
      /// asin 
      ValueWithError __asin__   () const ;
      /// acos
      ValueWithError __acos__   () const ;
      /// atan
      ValueWithError __atan__   () const ;
      /// asinh 
      ValueWithError __asinh__  () const ;
      /// acosh
      ValueWithError __acosh__  () const ;
      /// atanh
      ValueWithError __atanh__  () const ;
      /// sinc
      ValueWithError __sinc__   () const ;
      /// tgamma 
      ValueWithError __tgamma__ () const ;
      /// lgamma 
      ValueWithError __lgamma__ () const ;
      /// igamma 
      ValueWithError __igamma__ () const ;
      /// pdigamma
      ValueWithError __psi__    () const ;
      /// pdigamma
      ValueWithError __psi__    ( const unsigned short  n ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// printout 
      std::ostream& fillStream ( std::ostream& s  ) const ;          // printout 
      /// printout using format 
      std::ostream& fillStream 
      ( std::ostream&      s    , 
        const std::string& fmt  ) const ;          // printout 
      /// conversion to the string 
      std::string   toString   ( ) const ;
      /// conversion to the string using format 
      std::string   toString   ( const std::string& fmt ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// swap the content with another objects  
      inline void swap ( ValueWithError& right ) ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual value 
      double m_value ;                             //          the actual value 
      /// the associated covariance
      double m_cov2 ;                              // the associated covariance
      // ======================================================================
    } ;
    // ========================================================================
    /// addition 
    inline ValueWithError 
    operator+ ( const ValueWithError& v1 , const ValueWithError& v2 ) 
    { return v1.__add__ ( v2 ) ; }
    /// addition 
    inline ValueWithError 
    operator+ ( const ValueWithError& v1 , const double          v2 ) 
    { return v1.__add__ ( v2 ) ; }
    /// addition 
    inline ValueWithError 
    operator+ ( const double          v1 , const ValueWithError& v2 ) 
    { return v2 + v1 ; }
    // ========================================================================
    /// subtraction 
    inline ValueWithError 
    operator- ( const ValueWithError& v1 , const ValueWithError& v2 ) 
    { return v1.__sub__ ( v2 ) ; }
    /// subtraction 
    inline ValueWithError 
    operator- ( const ValueWithError& v1 , const double          v2 ) 
    { return v1.__sub__ ( v2 ) ; }
    /// subtraction 
    inline ValueWithError 
    operator- ( const double          v1 , const ValueWithError& v2 ) 
    { return v2.__rsub__ ( v1 ) ; }
    // ========================================================================
    /// multiplication 
    inline ValueWithError 
    operator* ( const ValueWithError& v1 , const ValueWithError& v2 ) 
    { return v1.__mul__ ( v2 ) ; }
    /// multiplication 
    inline ValueWithError 
    operator* ( const ValueWithError& v1 , const double          v2 ) 
    { return v1.__mul__ ( v2 ) ; }
    /// multiplication 
    inline ValueWithError 
    operator* ( const double          v1 , const ValueWithError& v2 ) 
    { return v2 * v1 ; }
    // ========================================================================
    /// division
    inline ValueWithError 
    operator/ ( const ValueWithError& v1 , const ValueWithError& v2 ) 
    { return v1.__div__ ( v2 ) ; }
    /// subtraction 
    inline ValueWithError 
    operator/ ( const ValueWithError& v1 , const double          v2 ) 
    { return v1.__div__ ( v2 ) ; }
    /// subtraction 
    inline ValueWithError 
    operator/ ( const double          v1 , const ValueWithError& v2 ) 
    { return v2.__rdiv__ ( v1 ) ; }
    // ========================================================================
    inline std::ostream& operator<< 
      ( std::ostream& s , const ValueWithError& v ) 
    { return v.fillStream ( s ) ; }
    // ========================================================================
    /// evaluate chi2 
    inline double chi2
    ( const ValueWithError& a , 
      const ValueWithError& b ) { return a.chi2 ( b ) ; }
    /// evaluate chi2 
    inline double chi2
    ( const ValueWithError& a , 
      const double          b ) { return a.chi2 ( b ) ; }
    /// evaluate chi2 
    inline double chi2
    ( const double          b ,
      const ValueWithError& a ) { return a.chi2 ( b ) ; }
    // ========================================================================
    /// evaluate the mean of a and b 
    inline ValueWithError mean
    ( const ValueWithError& a , 
      const ValueWithError& b ) { return a.mean ( b ) ; }
    // =======================================================================
    /** evaluate the mean of a and b 
     *  taking into account correlation coefficient <code>-1<=rho<=1</code>
     *  @param a (INPUT) the first argument 
     *  @param b (INPUT) the second argument 
     *  @param rho (INPUT) correlation coefficient \f$-1\le\rho\le1\f$
     */
    ValueWithError mean
    ( const ValueWithError& a   , 
      const ValueWithError& b   , 
      const double          rho ) ; 
    // ========================================================================
    /** get Kullback-Liebler divergency 
     *  @see Ostap::Math::ValueWithError::kullback
     *  @return the divergency for valid arguments, -1 otherwise
     */
    inline double kullback 
    ( const ValueWithError& a , 
      const ValueWithError& b ) { return a.kullback ( b ) ; }
    // ========================================================================
    /** get (squared) Hellinger distance for two (gaussian) variables 
     *  @see Ostap::Math::ValueWithError::hellinger
     *  @see https://en.wikipedia.org/wiki/Hellinger_distance
     *  @return squared Hellinger distance for two (gaussian) variables 
     */
    inline double hellinger2
    ( const ValueWithError& a ,
      const ValueWithError& b ) { return a.hellinger2 ( b ) ; }
    // ========================================================================
    /// evaluate the "fraction":  \f$ \frac{a}{a+b} \f$ 
    inline ValueWithError frac 
    ( const ValueWithError& a , 
      const ValueWithError& b ) { return a.frac ( b ) ; }
    /** \overload evaluate the "fraction": \f$ \frac{a}{a+b} \f$ 
     */
    inline ValueWithError frac 
    ( const ValueWithError& a , 
      const double          b ) { return a.frac ( b ) ; }
    /** \overload evaluate the "fraction": \f$ \frac{a}{a+b} \f$ 
     */
    inline ValueWithError frac 
    ( const double          a , 
      const ValueWithError& b ) { return frac ( ValueWithError ( a ) , b ) ; }
    // ========================================================================
    /// evaluate the "asymmetry":   \f$ \frac{a-b}{a+b} \f$ 
    inline ValueWithError asym 
    ( const ValueWithError& a , 
      const ValueWithError& b ) { return a.asym ( b ) ; }
    /** \overload evaluate the "asymmetry": \f$ \frac{a-b}{a+b} \f$  
     */
    inline ValueWithError asym 
    ( const ValueWithError& a , 
      const double          b ) { return a.asym ( b ) ; }
    /** \overload  evaluate the "asymmetry": \f$ \frac{a-b}{a+b} \f$ 
     */
    inline ValueWithError asym 
    ( const double          a , 
      const ValueWithError& b ) { return asym ( ValueWithError ( a ) , b ) ; }
    // ========================================================================
    /** Evaluate <code>abs(a)</code> 
     *  @param a (INPUT) the value 
     *  @return the absolute value 
     */
    ValueWithError abs 
    ( const ValueWithError& a ) ;
    // ========================================================================
    /** evaluate the binomial efficiency for Bernulli scheme 
     *  @param n_success (INPUT) number of 'success' 
     *  @param N_total   (INPUT) total number 
     *  @return the binomial efficiency 
     */
    ValueWithError binomEff   
    ( const size_t n_success , 
      const size_t N_total   ) ;
    // ========================================================================
    /** evaluate the binomial efficiency interval using Wilson's prescription
     *  @param n_success (INPUT) number of 'success' 
     *  @param N_total   (INPUT) total number 
     *  @return the binomial efficiency 
     */
    ValueWithError wilsonEff   
    ( const size_t n_success , 
      const size_t N_total   ) ;
    // ========================================================================
    /** evaluate the binomial efficiency interval 
     *  using Agresti-Coull's prescription
     *  @param n_success (INPUT) number of 'success' 
     *  @param N_total   (INPUT) total number 
     *  @return the binomial efficiency 
     */
    ValueWithError agrestiCoullEff   
    ( const size_t n_success , 
      const size_t N_total   ) ;
    // ========================================================================
    /** simple evaluation of efficiency from statistically independend
     * "exclusive" samples "accepted" and "rejected"
     *  \f[ \varepsilon = \frac{1}{ 1 + \frac{N_{rejected}}{N_{accepted}}}\f] 
     *  @param accepted  (IN) accepted sample 
     *  @param rejected  (IN) rejected sample 
     *  @return the binomial efficiency 
     *  @see Ostap::Math::binomEff2
     */
    ValueWithError exclusiveEff
    ( const ValueWithError& accepted , 
      const ValueWithError& rejected ) ;
    // ========================================================================
    /** evaluate the binomial efficiency for Bernulli scheme with weights 
     *  \f[ R = \frac{N_{acc}}{N_{tot}} = \frac{N_{acc}}{N_{acc}+N_{rej}} = 
     *          \left( 1 + \frac{N_{rej}}{N_{acc}}\right)^{-1} \f]
     *  @param nAccepted (INPUT) number of accepted (weighted) events 
     *  @param nRejected (INPUT) number of rejected (weighted) events 
     *  @return the binomial efficiency 
     *  @see Ostap::Math::exclusiveEff 
     */
    ValueWithError binomEff2   
    ( const ValueWithError& nAccepted , 
      const ValueWithError& nRejected ) ;
    // ========================================================================
    /** simple evaluation of efficiency using Zech's prescription 
     *  @param n_success (IN) accepted sub-sample 
     *  @param N_total   (IN) total     sample 
     */
    ValueWithError zechEff
    ( const ValueWithError& n_success , 
      const ValueWithError& N_total   ) ;
    // ========================================================================
    /** calculate the ratio of weighted to unweighted statistic 
     *  \f[ R = \frac{N_w}{N}  = \frac{ \sum_1^{N} w_i }{N} \f] 
     *  using Jackknife' method:
     *  \f[ \sigma^2(R) = \left( \sum_1^N w_i^2 - NR^2 \right) / (N-1)^2 \f] 
     *  - thanks to Wouter Hulsbergen 
     *  @see http://en.wikipedia.org/wiki/Jackknife_%28statistics%29
     *  The result has proper behaviour : 
     *  uncertainty in R goes to zero if 
     *  dispersion on weights go to zero.
     *  @param   nWeighted (input) statistic of weighted sample 
     *  @param   n         (input) size      of origial sample 
     *  @return  ratio R with the proper uncertaities 
     *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
     *  @date 2014-01-25
     */
    ValueWithError effJackknife  
    ( const ValueWithError& nWeighted , 
      const unsigned long   n         ) ;  
    // ========================================================================
    /** \overload evaluate <code>pow(a,b)</code>:   \f$ a^b \f$ 
     *  @param a (INPUT) the base 
     *  @param b (INPUT) the exponent 
     *  @return the <c>a</c> raised to power <c>b</c> 
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError pow 
    ( const ValueWithError& a , 
      const int             b ) ;
    // ========================================================================    
    /** \overload evaluate <code> pow(a,b)<code>:  \f$ a^b \f$ 
     *  @param a (INPUT) the base 
     *  @param b (INPUT) the exponent 
     *  @return the <c>a</c> raised to power <c>b</c> 
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError pow 
    ( const ValueWithError& a , 
      const double          b ) ;
    // ========================================================================
    /** \overload evaluate <code>pow(a,b)</code>:  \f$ a^b \f$ 
     *  @param a (INPUT) the base 
     *  @param b (INPUT) the exponent 
     *  @return the <c>a</c> raised to power <c>b</c> 
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError pow 
    ( const int             a , 
      const ValueWithError& b ) ;
    // ========================================================================    
    /** \overload evaluate <code>pow(a,b)</code>: \f$ a^b \f$ 
     *  @param a (INPUT) the base 
     *  @param b (INPUT) the exponent 
     *  @return the <c>a</c> raised to power <c>b</c> 
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError pow 
    ( const double          a , 
      const ValueWithError& b ) ;
    // ========================================================================    
    /** \overload evaluate <code>pow(a,b)</code>: \f$ a^b \f$ 
     *  @param a (INPUT) the base 
     *  @param b (INPUT) the exponent 
     *  @return the <c>a</c> raised to power <c>b</c> 
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError pow 
    ( const ValueWithError& a , 
      const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate <code>exp(b)</code>: \f$ \mathrm{e}^b \f$ 
     *  @param b (INPUT) the exponent 
     *  @return the <c>e</c> raised to power <c>b</c> 
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError exp
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate <code>exp2(b)<code>: \f$ 2^b \f$ 
     *  @param b (INPUT) the exponent 
     *  @return the <c>2</c> raised to power <c>b</c> 
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError exp2
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate <code>expm1(b)</code>: \f$ \mathrm{e}^b-1 \f$ 
     *  @param b (INPUT) the exponent 
     *  @return  expm1(b) 
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError expm1
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate <code>log(b)</code>: \f$ \log b \f$
     *  @param b (INPUT) the parameter 
     *  @return logarithm
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError log
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate <code>log2(b)</code>: \f$ \log_2 b \f$  
     *  @param b (INPUT) the parameter 
     *  @return logarithm on base 2 
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError log2
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate <code>log10(b)</code>: \f$ \log_{10} b \f$  
     *  @param b (INPUT) the parameter 
     *  @return logarithm on base 10 
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError log10
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate <code>log1p(b)</code>:  \f$ \log (1+b) \f$
     *  @param b (INPUT) the parameter 
     *  @return  log1p(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError log1p
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate <code>sqrt(b)</code>: \f$ \sqrt{b}\f$
     *  @param b (INPUT) the parameter 
     *  @return  sqrt(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError sqrt
    ( const ValueWithError& b ) ;
    // ========================================================================
    /** \overload evaluate "signed sqrt": \f$ \sign(b)\sqrt{\left|b\right|}\f$ 
     *  @param a (INPUT) the value 
     *  @return the signed-sqrt  
     */
    ValueWithError signed_sqrt 
    ( const ValueWithError& a ) ;
    // ========================================================================    
    /** evaluate <code>cbrt(b)</code>: \f$ \sqrt[3]{b}\f$
     *  @param b (INPUT) the parameter 
     *  @return  cbrt(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError cbrt
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate <code>sin(b)</code>:  \f$ \sin{b}\f$
     *  @param b (INPUT) the parameter 
     *  @return  sin(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError sin 
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate <code>cos(b)</code>: \f$ \cos{b}\f$
     *  @param b (INPUT) the parameter 
     *  @return  cos(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError cos 
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate <code>tan(b)</code>: \f$ \tan{b}\f$
     *  @param b (INPUT) the parameter 
     *  @return  tan(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError tan 
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate <code>sinh(b)</code>:  \f$ \sinh{b}\f$
     *  @param b (INPUT) the parameter 
     *  @return  sinh(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError sinh 
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate <code>cosh(b)</code>:  \f$ \cosh{b}\f$
     *  @param b (INPUT) the parameter 
     *  @return  cosh(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError cosh 
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate <code>tanh(b)</code>:   \f$ \tanh{b}\f$
     *  @param b (INPUT) the parameter 
     *  @return  tanh(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError tanh 
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** \overload evaluate <code>sech(b)</code>: \f$ \frac{1}{\cosh b}\f$
     *  @param b (INPUT) the parameter 
     *  @return  sech(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError sech 
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** \overload evaluate <code>erf(b)</code>, the error function 
     *  @param b (INPUT) the parameter 
     *  @return  erf(b)
     *  @warning invalid and small covariances are ignored 
     *  @see https://en.wikipedia.org/wiki/Error_function
     */
    ValueWithError erf
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** \overload evaluate <code>erfc(b)</code>, complementary error function 
     *  @param b (INPUT) the parameter 
     *  @return  erfc(b)
     *  @warning invalid and small covariances are ignored 
     *  @see https://en.wikipedia.org/wiki/Error_function
     */
    ValueWithError erfc
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** \overload evaluate <code>erfi(b)</code>,  imaginary error function 
     *  @param b (INPUT) the parameter 
     *  @return  erfc(b)
     *  @warning invalid and small covariances are ignored 
     *  @see https://en.wikipedia.org/wiki/Error_function
     */
    ValueWithError erfi
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** \overload evaluate <code>erfcx(b)</code>, complementary scaled error function 
     *  @param b (INPUT) the parameter 
     *  @return  erfc(b)
     *  @warning invalid and small covariances are ignored 
     *  @see https://en.wikipedia.org/wiki/Error_function
     */
    ValueWithError erfcx
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** \overload evaluate <code>probit(b)</code> 
     *  @param b (INPUT) the parameter 
     *  @return  erfc(b)
     *  @see https://en.wikipedia.org/wiki/Probit
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError probit
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate <code>asin(b)</code>: \f$ \asin  b \f$ 
     *  @param b (INPUT) the parameter 
     *  @return  asin(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError asin
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate <code>acos(b)</code>:  \f$ \acos b \f$ 
     *  @param b (INPUT) the parameter 
     *  @return  acos(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError acos
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate <code>atan(b)</code>:  \f$  \atan b \f$
     *  @param b (INPUT) the parameter 
     *  @return  atan(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError atan
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate <code>atan2(y,x)</code> 
     *  @param y    (INPUT) the parameter 
     *  @param x    (INPUT) the parameter 
     *  @param corr (INPUT) the correlation coefficient   -1<=corr<=1 
     *  @return  atan2(y,x)
     *  @warning invalid and small covariances are ignored
     *  @warning no error estimate for x=y=0 point 
     */
    ValueWithError atan2
    ( const ValueWithError& y        ,
      const ValueWithError& x        , 
      const double          corr = 0 ) ;
    // ========================================================================    
    /** evaluate <code>asinh(b)</code>  
     *  @param b (INPUT) the parameter 
     *  @return  asinh(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError asinh
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate <code>acosh(b)</code>
     *  @param b (INPUT) the parameter 
     *  @return  acosh(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError acosh
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate <code>atanh(b)</code>
     *  @param b (INPUT) the parameter 
     *  @return  atanh(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError atanh
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate  <code>sinc(x)</code>: \f$ \frac{\sin x}{x}\f$
     *  @param x (INPUT) parameter
     *  @return sinc(x)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError sinc 
    ( const ValueWithError& b ) ;     
    // ========================================================================    
    /** evaluate <code>tgamma(b)</code>:  \f$ \Gamma(b) \f$ 
     *  @param b (INPUT) the parameter 
     *  @return  Gamma(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError tgamma
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate <code>lgamma(b)</code>: \f$ \log\Gamma (b)\f$
     *  @param b (INPUT) the parameter 
     *  @return  log(Gamma(b))
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError lgamma
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate <code>igamma(b)</code>:  \f$ \frac{1}{\Gamma(b)}\f$ 
     *  @param b (INPUT) the parameter 
     *  @return  1/Gamma(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError igamma
    ( const ValueWithError& b ) ;
    // ========================================================================
    /** \overload evaluate Pochhammer symbol 
     *  \f[ (x)^n = x ( x + 1) ( x + 1 ) ... ( x + n - 1 ) = \Pi^{k-1}_{k=0} (x + k) \f] 
     *  @see https://en.wikipedia.org/wiki/Falling_and_rising_factorials
     *  @param x (INPUT) the parameter 
     *  @param n (INPUT) the parameter 
     *  @return  pochhammer  symbol 
     *  @warning invalid and small covariances are ignored 
     *  @see Ostap::Math::rising_factorial
     *  @see Ostap::Math::falling_factorial
     *  @see Ostap::Math::pochhammer 
     */
    ValueWithError pochhammer 
    ( const ValueWithError& x ,  
      const unsigned short  n ) ;
    // ========================================================================
    /** Complete elliptic integral \f$ K(k) \f$  
     *  \[ K(k) \equiv F ( \frac{\pi}{2}, k ) \f] 
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     */
    ValueWithError elliptic_K ( const ValueWithError& k ) ;
    // ========================================================================
    /**  Complete elliptic integral \f$ E(k) \f$  
     *  \[ E(k) \equiv E ( \frac{\pi}{2}, k ) \f] 
     *  @see https://en.wikipedia.org/wiki/Elliptic_integral
     */
    ValueWithError elliptic_E ( const ValueWithError& k ) ;
    // ========================================================================
    /** Triangle, aka Kallen function 
     *  \f[ \lambda ( a , b, c ) = a^2 + b^2 + c^2 - 2ab - 2bc - 2ca \f]
     *  @see see https://en.wikipedia.org/wiki/K%C3%A4ll%C3%A9n_function          
     *  @see Ostap::Kinematics::triangle 
     *  @see Ostap::Kinematics::kallen 
     *  @see Ostap::Math::kallen
     */
    ValueWithError triangle 
    ( const ValueWithError& x , 
      const double          y , 
      const double          z ) ;                           
    // ========================================================================    
    /** Triangle, aka Kallen function 
     *  \f[ \lambda ( a , b, c ) = a^2 + b^2 + c^2 - 2ab - 2bc - 2ca \f]
     *  @see see https://en.wikipedia.org/wiki/K%C3%A4ll%C3%A9n_function          
     *  @see Ostap::Kinematics::triangle 
     *  @see Ostap::Kinematics::kallen 
     *  @see Ostap::Math::triangle 
     */
    inline ValueWithError kallen 
    ( const ValueWithError& x , 
      const double          y , 
      const double          z ) { return triangle ( x , y , z ) ; }
    // ========================================================================
    /** momentum in the rest frame for the two-body decay  \f$ n\rightarrow m_1 m_2 \f$ 
     *  \f[ q ( m , m_1 , m_2 )  \equiv 
     *  \frac{\lambda^{1/2}\left( m^2, m_1^2, m_2^2\right)}{2m} \f]
     *  @see  Ostap::Kinematics::q 
     */
    ValueWithError q
    ( const ValueWithError& m  , 
      const double          m1 ,
      const double          m2 ) ;
    // ========================================================================
    /** \overload evaluate standard Gauss PDF 
     *  @param x the value 
     *  @return valeu of the standard gaussian PDF  
     */
    ValueWithError gauss_pdf 
    ( const ValueWithError& x         , 
      const double          mu    = 0 , 
      const double          sigma = 1 ) ;
    // ========================================================================
    /** \overload evaluate standard Gauss CDF 
     *  @param x the value 
     *  @return value of the standard gaussian CDF  
     */
    ValueWithError gauss_cdf  
    ( const ValueWithError& x         ,
      const double          mu    = 0 , 
      const double          sigma = 1 ) ;
    // ========================================================================    
    /** evaluate <code>hypot(x,y)</code>: \f$ \sqrt( x^2 + y^2 ) \f$
     *  @param x (INPUT) the first parameter
     *  @param y (INPUT) the second parameter
     *  @param c (INPUT) the correlation coefficient  (-1<=c<=1)
     *  @return the value of <code>hypot</code> function
     */
    ValueWithError  hypot
    ( const ValueWithError& x     , 
      const ValueWithError& y     , 
      const double          c = 0 ) ;
    // ========================================================================    
    /** evaluate <code>hypot(x,y)</code>: \f$ \sqrt( x^2 + y^2 ) \f$
     *  @param x (INPUT) the first parameter
     *  @param y (INPUT) the second parameter
     *  @param c (INPUT) the correlation coefficient  (-1<=c<=1)
     *  @return the value of <code>hypot</code> function
     */
    ValueWithError  hypot
    ( const ValueWithError& x     , 
      const double          y     ) ;
    // ========================================================================
    /** evaluate <code>hypot(x,y)</code>: \f$ \sqrt( x^2 + y^2 ) \f$
     *  @param x (INPUT) the first parameter
     *  @param y (INPUT) the second parameter
     *  @param c (INPUT) the correlation coefficient  (-1<=c<=1)
     *  @return the value of <code>hypot</code> function
     */
    inline ValueWithError  hypot
    ( const double          x     , 
      const ValueWithError& y     ) { return hypot ( y , x ) ; }
    // ========================================================================    
    /** evaluate beta-function \f$ \Beta(x,y) \f$ 
     *  @param x (INPUT) the first parameter
     *  @param y (INPUT) the second parameter
     *  @param c (INPUT) the correlation coefficient  (-1<=c<=1)
     *  @return the value of beta function 
     */
    ValueWithError beta 
    ( const ValueWithError& x     , 
      const ValueWithError& y     , 
      const double          c = 0 ) ;
    // =========================================================================
    /** evaluate beta-function \f$ \Beta(x,y) \f$ 
     *  @param x (INPUT) the first parameter
     *  @param y (INPUT) the second parameter
     *  @return the value of beta function 
     */
    ValueWithError beta 
    ( const ValueWithError& x     , 
      const double          y     ) ;
    // ========================================================================
    /** evaluate beta-function \f$ \Beta(x,y) \f$ 
     *  @param x (INPUT) the first parameter
     *  @param y (INPUT) the second parameter
     *  @return the value of beta function 
     */
    inline ValueWithError beta 
    ( const double          x  , 
      const ValueWithError& y  ) { return beta ( y , x ) ; }
    // =========================================================================
    /** evaluate log(beta-function) \f$ \log \Beta(x,y) \f$ 
     *  @param x (INPUT) the first parameter
     *  @param y (INPUT) the second parameter
     *  @param c (INPUT) the correlation coefficient  (-1<=c<=1)
     *  @return the value of log of the beta function 
     */
    ValueWithError lnbeta 
    ( const ValueWithError& x     , 
      const ValueWithError& y     , 
      const double          c = 0 ) ;
    // =========================================================================
    /** evaluate log(beta-function) \f$ \log \Beta(x,y) \f$ 
     *  @param x (INPUT) the first parameter
     *  @param y (INPUT) the second parameter
     *  @return the value of log of the beta function 
     */
    ValueWithError lnbeta 
    ( const ValueWithError& x     , 
      const double          y     ) ;
    // ========================================================================
    /** evaluate log(beta-function) \f$ \log \Beta(x,y) \f$ 
     *  @param x (INPUT) the first parameter
     *  @param y (INPUT) the second parameter
     *  @return the value of log of the beta function 
     */
    inline ValueWithError lnbeta 
    ( const double          x     ,
      const ValueWithError& y     ) { return lnbeta ( y , x ) ; }
    // =========================================================================
    /** calculate psi/digamma function 
     *  @see Ostap::Math::digamma 
     *  @see Ostap::Math::polygamma
     *  @see Ostap::Math::psi 
     *  @return the value of digamma function 
     */
    ValueWithError psi 
    ( const ValueWithError& x ) ;
    // ========================================================================
    /** calculate psi/polygamma function 
     *  @see Ostap::Math::digamma 
     *  @see Ostap::Math::polygamma
     *  @see Ostap::Math::psi 
     *  @return the value of polygamma function 
     */
    ValueWithError psi 
    ( const ValueWithError& x , 
      const unsigned short  n ) ;
    // ========================================================================
    /** calculate psi/digamma function 
     *  @see Ostap::Math::digamma 
     *  @see Ostap::Math::polygamma
     *  @see Ostap::Math::psi 
     *  @return the value of digamma function 
     */
    inline ValueWithError digamma 
    ( const ValueWithError& x ) { return psi ( x ) ; }
    // ========================================================================
    /** calculate trigamma function 
     *  @see Ostap::Math::digamma 
     *  @see Ostap::Math::polygamma
     *  @see Ostap::Math::psi 
     *  @return the value of digamma function 
     */
    inline ValueWithError trigamma 
    ( const ValueWithError& x ) { return psi ( x , 1 ) ; }
    // ========================================================================
    /** calculate psi/polygamma function 
     *  @see Ostap::Math::digamma 
     *  @see Ostap::Math::polygamma
     *  @see Ostap::Math::psi 
     *  @return the value of polygamma function 
     */
    inline ValueWithError polygamma 
    ( const ValueWithError& x , 
      const unsigned short  n ) { return psi ( x , n ) ; }
    // ========================================================================    
    /** evaluate <code>fma(x,y,z)</code>: \f$ xy+z \f$  
     *  @param y    (INPUT) the parameter 
     *  @param x    (INPUT) the parameter 
     *  @param z    (INPUT) the parameter 
     *  @param cxy  (INPUT) the correlation coefficient   -1<=c_xy<=1 
     *  @param cxz  (INPUT) the correlation coefficient   -1<=c_xz<=1 
     *  @param cyz  (INPUT) the correlation coefficient   -1<=c_yz<=1 
     *  @return  fma(x,y,z)
     *  @warning invalid and small covariances are ignored
     */
    ValueWithError fma 
    ( const ValueWithError& x        ,
      const ValueWithError& y        , 
      const ValueWithError& z        , 
      const double          cxy = 0  ,
      const double          cxz = 0  ,
      const double          cyz = 0  ) ;
    // ========================================================================
    /// Regular Bessel function \f$ J_n(s)\f$
    ValueWithError bessel_Jn 
    ( const int             n  , 
      const ValueWithError& x  ) ;
    // ========================================================================
    /// Regular Bessel function \f$ J_{\nu}(s)\f$
    ValueWithError bessel_Jnu 
    ( const double          nu , 
      const ValueWithError& x  ) ;
    // ========================================================================
    /// Irregular Bessel function \f$ Y_n(s)\f$
    ValueWithError bessel_Yn 
    ( const int             n  , 
      const ValueWithError& x  ) ;
    // ========================================================================
    /// Irregular Bessel function \f$ Y_{\nu}(s)\f$
    ValueWithError bessel_Ynu 
    ( const double          nu , 
      const ValueWithError& x  ) ;
    // ========================================================================
    /// Modified Bessel function \f$ I_n(s)\f$
    ValueWithError bessel_In
    ( const int             n  , 
      const ValueWithError& x  ) ;
    // ========================================================================
    /// Modified Bessel function \f$ I_{\nu}(s)\f$
    ValueWithError bessel_Inu 
    ( const double          nu , 
      const ValueWithError& x  ) ;
    // ========================================================================
    /// Modified Bessel function \f$ K_n(s)\f$
    ValueWithError bessel_Kn 
    ( const int             n  , 
      const ValueWithError& x  ) ;
    // ========================================================================
    /// Modified Bessel function \f$ K_{\nu}(s)\f$
    ValueWithError bessel_Knu 
    ( const double          nu , 
      const ValueWithError& x  ) ;
    // ========================================================================
    /// check for NaN
    inline bool isnan    ( const ValueWithError& v ) { return v.isnan    () ; }
    /// finite ?
    inline bool isfinite ( const ValueWithError& v ) { return v.isfinite () ; }
    /// infinte ? 
    inline bool isinf    ( const ValueWithError& v ) { return v.isinf    () ; }    
    /// normal ?
    inline bool isnormal ( const ValueWithError& v ) { return v.isnormal () ; }    
    /// check for goodness 
    inline bool isgood   ( const ValueWithError& v ) { return v.isgood   () ; }    
    /// check for goodness 
    inline bool good     ( const ValueWithError& v ) { return v.good     () ; }    
    // ========================================================================    
    /** Does this object represent natural number?
     *  - non-negative integer value 
     *  - cov2 == value  or cov2 == 0 
     */
    bool natural_number 
    ( const ValueWithError& value ) ;
    // ========================================================================
    /** Does this object represent natural entry in histogram
     *  - non-negative integer value 
     *  - cov2 == value or ( 0 == value && 1 == cov2 ) 
     */
    bool natural_entry 
    ( const ValueWithError& value ) ;
    // ========================================================================
    /** make a sum two elements taking into account the 
     *  correlation coefficient  
     *  @param a  (input) the first value 
     *  @param b  (input) the second value 
     *  @param c  (input) the correlation coefficient
     *  @return a+b 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2012-11-09
     */
    ValueWithError  sum 
    ( const ValueWithError& a     , 
      const ValueWithError& b     , 
      const double          c = 0 ) ;
    // ========================================================================
    /** make a sum two elements taking into account the 
     *  correlation coefficient  
     *  @param a  (input) the first value 
     *  @param b  (input) the second value 
     *  @param c  (input) the correlation coefficient
     *  @return a+b 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2012-11-09
     */
    ValueWithError  sum2 
    ( const ValueWithError& a     , 
      const ValueWithError& b     , 
      const double          c = 0 ) ;
    // ========================================================================
    /** make a subtraction of two elements taking into account the 
     *  correlation coefficient  
     *  @param a  (input) the first value 
     *  @param b  (input) the second value 
     *  @param c  (input) the correlation coefficient
     *  @return a-b 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2012-11-09
     */
    ValueWithError  subtract
    ( const ValueWithError& a     , 
      const ValueWithError& b     , 
      const double          c = 0 ) ;
    // ========================================================================
    /** make a  multiplication two elements taking into account the 
     *  correlation coefficient  
     *  @param a  (input) the first value 
     *  @param b  (input) the second value 
     *  @param c  (input) the correlation coefficient
     *  @return a*b 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2012-11-09
     */
    ValueWithError  multiply
    ( const ValueWithError& a     , 
      const ValueWithError& b     , 
      const double          c = 0 ) ;
    // ========================================================================
    /** make a division of two elements taking into account the 
     *  correlation coefficient  
     *  @param a  (input) the first value 
     *  @param b  (input) the second value 
     *  @param c  (input) the correlation coefficient
     *  @return a/b 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2012-11-09
     */
    ValueWithError  divide
    ( const ValueWithError& a     , 
      const ValueWithError& b     , 
      const double          c = 0 ) ;
    // ========================================================================
    /** calculate "fraction" of two elements \f$ \alpha = \frac{a}{a+b} \f$ 
     *  taking into account the correlation coefficient  
     *  @param a  (input) the first value 
     *  @param b  (input) the second value 
     *  @param c  (input) the correlation coefficient
     *  @return a/(a+b)
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2012-11-09
     */
    ValueWithError  fraction
    ( const ValueWithError& a     , 
      const ValueWithError& b     , 
      const double          c = 0 ) ;
    // ========================================================================
    /** calculate the "effective" background-to-signal ratio from the valeu 
     *  and its uncertainty using the identity
     *  \f$ \frac{\sigma(S)}{S} = \frac{\sqrt{S}}{S}\sqrt{1+\frac{B}{S}}\f$.
     *  From this identity one has
     *  \f$ \left.\frac{B}{S}\right|_{\mathrm{eff}} \equiv \frac{\sigma^2(S)}{S} -1 \f$
     *  @param v the value 
     *  @return the effective backround-to-signal ratio or -1 
     */
    ValueWithError b2s   
    ( const ValueWithError& v ) ;
    // ========================================================================
    /** calculate the "effective purity" ratio using the identity
     *  \f[ p_{\mathrm{eff}} = \frac{S}{S+B} = \frac{1}{1+\frac{B}{S}}\f]
     *  and the effective "background-to-signal" ratio is estimated as 
     *  \f[ \left.\frac{B}{S}\right|_{\mathrm{eff}} = \frac{\sigma^2(S)}{S} -1 \f], 
     *  finally one gets 
     *  \f[ p_{\mathrm{eff}} \equiv \frac{S}{\sigma^2(S)}\f]
     *  @see Ostap::Math::b2s 
     *  @param v the value 
     *  @return the effective purity or -1 
     */
    ValueWithError purity 
    ( const ValueWithError& v ) ;
    // ========================================================================
    /** calculate "asymmetry" of two elements: \f$ \kappa = \frac{a-b}{a+b}\f$ 
     *  taking into account the correlation coefficient  
     *  @param a  (input) the first value 
     *  @param b  (input) the second value 
     *  @param c  (input) the correlation coefficient
     *  @return (a-b)/(a+b)
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2012-11-09
     */
    ValueWithError  asymmetry
    ( const ValueWithError& a     , 
      const ValueWithError& b     , 
      const double          c = 0 ) ;
    // ========================================================================
    /** simple linear interpolation 
     *  @param  x  the point to evaluate the function 
     *  @param  x0 the abscissa for the first  point
     *  @param  y0 the function value for the first  point
     *  @param  x1 the abscissa for the second point
     *  @param  y1 the function value for the second point
     *  @return linear interpolation at point x
     */
    ValueWithError interpolate_1D 
    ( const double          x  , 
      const double          x0 ,
      const ValueWithError& y0 , 
      const double          x1 ,
      const ValueWithError& y1 ) ;    
    // ========================================================================
    /** simple (bi)linear interpolation 
     *  @param x  the x-coordiate to evaluate the function 
     *  @param y  the y-coordiate to evaluate the function 
     *  @param x0 the x-coordinate for the first  pair of points
     *  @param x1 the x-coordinate for the second pair of points
     *  @param y0 the y-coordinate for the first  pair of points
     *  @param y1 the y-coordinate for the second pair of points
     *  @param v00 the function value 
     *  @param v01 the function value 
     *  @param v10 the function value 
     *  @param v11 the function value 
     *  @return bilinear interpolation at point (x,y)
     */
    ValueWithError interpolate_2D 
    ( const double          x   , 
      const double          y   , 
      const double          x0  ,
      const double          x1  ,
      const double          y0  ,
      const double          y1  ,
      const ValueWithError& v00 , 
      const ValueWithError& v01 , 
      const ValueWithError& v10 , 
      const ValueWithError& v11 ) ;
    // ========================================================================
    /** simple interpolation 
     *  - if vector of y is larger  than vector of x, extra values are ignored 
     *  - if vector of y is shorter than vector of x, missing entries assumed to be zero 
     *  @param y_i        INPUT  vector of yi 
     *  @param x_i        INPUT  vector of xi
     *  @param x          INPUT  the point where the function to be evaluated 
     *  @param correlated INPUT  correlated uncertaties in yi?
     */
    ValueWithError interpolate 
    ( const std::vector<ValueWithError>& y_i                , 
      const std::vector<double>&         x_i                ,
      const double                       x                  , 
      const bool                         correlated = false ) ;
    // ========================================================================
    /** simple interpolation 
     *  - if vector of y is larger  than vector of x, extra values are ignored 
     *  - if vector of y is shorter than vector of x, missing entries assumed to be zero 
     *  @param y_i        INPUT  vector of yi 
     *  @param x_i        INPUT  vector of xi
     *  @param x          INPUT  the point where the function to be evaluated 
     */
    ValueWithError interpolate 
    ( const std::vector<double>&         y_i                , 
      const std::vector<double>&         x_i                ,
      const ValueWithError&              x                  ) ;
    // ========================================================================
    /** simple interpolation 
     *  - if vector of y is larger  than vector of x, extra values are ignored 
     *  - if vector of y is shorter than vector of x, missing entries assumed to be zero 
     *  @param y_i        INPUT  vector of yi 
     *  @param x_i        INPUT  vector of xi
     *  @param x          INPUT  the point where the function to be evaluated 
     *  @param correlated INPUT  correlated uncertaties in yi?
     */
    ValueWithError interpolate 
    ( const std::vector<ValueWithError>& y_i                , 
      const std::vector<double>&         x_i                ,
      const ValueWithError&              x                  , 
      const bool                         correlated = false ) ;
    // ========================================================================
    /** get the sum of the vector 
     *  @param vct the vector
     *  @param ini the intial value 
     *  @return the sum over the vector 
     */
    ValueWithError sum 
    ( const std::vector<ValueWithError>& vct                    , 
      ValueWithError                     ini = ValueWithError() ) ;
    // ========================================================================    
    /** get the sum of the vector 
     *  @param vct the vector
     *  @param ini the intial value 
     *  @return the sum over the vector 
     */
    ValueWithError accumulate 
    ( const std::vector<ValueWithError>& vct                    , 
      ValueWithError                     ini = ValueWithError() ) ;
    // ========================================================================    
    /** get the sum of the absolute values
     *  @param vct the vector
     *  @return the sum over the vector 
     */
    ValueWithError abssum
    ( const std::vector<ValueWithError>& vct ) ;
    // ========================================================================    
    /** get the sum of the absolute values
     *  @param vct the vector
     *  @return the sum over the vector 
     */
    inline 
    ValueWithError sumabs 
    ( const std::vector<ValueWithError>& vct )
    { return abssum ( vct )  ; }
    // ========================================================================    
    /** evaluate polynomial
     *  \f[ f(x) = Чаa_0 + a_1x + a_2x^2 + ... + a_{n-1}x^{n-1} + a_nx^n\f]
     *  such as \f$f(0) = a_0 \f$      
     *  using Horner rule
     *  @param poly  INPUT the coefficients
     *  @param x     INPUT argument 
     *  @return value of polynomial
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError horner_a0 
    ( const std::vector<double>& poly ,
      const ValueWithError&      x    ) ;
    // ======================================================================
    /** evaluate polynomial
     *  \f[ f(x) = a_0x^n + a_1x^{n-1}+ ... + a_{n-1}x + a_n\f], 
     *  such as \f$f(0) = a_n \f$      
     *  using Horner rule
     *  @param poly  INPUT the coefficients 
     *  @param x     INPUT argument 
     *  @return value of polynomial
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError horner_aN
    ( const std::vector<double>& poly ,
      const ValueWithError&      x    ) ;
    // ========================================================================    
    /// the output operator for the vector 
    std::ostream& operator<<
    ( std::ostream& s , 
      const std::vector<Ostap::Math::ValueWithError>& v ) ;
    // ========================================================================
    /// swap the content with another object
    inline void ValueWithError::swap 
    ( ValueWithError& right ) 
    {
      std::swap ( m_value , right.m_value ) ;
      std::swap ( m_cov2  , right.m_cov2  ) ;
    }
    // ========================================================================
    /// swap the content of two objects 
    inline void swap 
    ( ValueWithError& left , 
      ValueWithError& right ) { left.swap ( right ) ; }
    // ========================================================================
    // converison to string
    inline std::string to_string 
    ( const ValueWithError& v ) 
    { return v.toString() ; }
    // ========================================================================
  } //                                         The end of namespace Ostap::Math 
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_VALUEWITHERRORS_H
// ============================================================================
