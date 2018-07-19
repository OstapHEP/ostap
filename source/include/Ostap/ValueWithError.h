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
/** @file 
 *  Collection of useful objects with associated "covariances".
 *  The concept has been stollen from Wouter Hulsbergen's lines 
 *  @author Vanya BELYAEV Ivane.Belyaev@itep.ru
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
     *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
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
      ValueWithError __add__  ( const ValueWithError& right ) const ;
      ///    a + right 
      ValueWithError __add__  ( const double          right ) const ;
      ///        right + a  
      ValueWithError __radd__ ( const double          right ) const { return __add__ ( right ) ; }
      ///    a - right  
      ValueWithError __sub__  ( const ValueWithError& right ) const ;
      ///    a - right  
      ValueWithError __sub__  ( const double          right ) const ;
      ///        right - a   
      ValueWithError __rsub__ ( const double          right ) const ;
      ///    a * right  
      ValueWithError __mul__  ( const ValueWithError& right ) const ;
      ///    a * right  
      ValueWithError __mul__  ( const double          right ) const ;
      ///        right * a   
      ValueWithError __rmul__ ( const double          right ) const { return __mul__ ( right ) ; }
      ///    a / right  
      ValueWithError __div__  ( const ValueWithError& right ) const ;
      ///    a / right  
      ValueWithError __div__  ( const double          right ) const ;
      ///        right / a   
      ValueWithError __rdiv__ ( const double          right ) const ;
      ///  abs ( a )   
      ValueWithError __abs__  () const ;
      /// -me
      ValueWithError __neg__  () const ;
      /// +me  (no effect) 
      ValueWithError __pos__  () const ; 
      /// me**e 
      ValueWithError __pow__  ( const int             e ) const ;
      /// me**e 
      ValueWithError __pow__  ( const double          e ) const ;
      /// me**e 
      ValueWithError __pow__  ( const ValueWithError& e ) const ;
      /// e**me 
      ValueWithError __rpow__ ( const int             e ) const ;
      /// e**me 
      ValueWithError __rpow__ ( const double          e ) const ;
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
      /// tgamma 
      ValueWithError __tgamma__ () const ;
      /// lgamma 
      ValueWithError __lgamma__ () const ;
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
     *  return the divergency for valid arguments, -1 otherwise
     */
    inline double kullback 
    ( const ValueWithError& a , 
      const ValueWithError& b ) { return a.kullback ( b ) ; }
    // ========================================================================
    /** get (squared) Hellinger distance for two (gaussian) varibales 
     *  @see https://en.wikipedia.org/wiki/Hellinger_distance
     *  return Hellinger distance for two (gaussian) varibales 
     */
    inline double hellinger2
    ( const ValueWithError& a ,
      const ValueWithError& b ) { return a.hellinger2 ( b ) ; }
    // ========================================================================
    /// evaluate the "fraction"  a/(a+b)
    inline ValueWithError frac 
    ( const ValueWithError& a , 
      const ValueWithError& b ) { return a.frac ( b ) ; }
    /// evaluate the "fraction"  a/(a+b)
    inline ValueWithError frac 
    ( const ValueWithError& a , 
      const double          b ) { return a.frac ( b ) ; }
    /// evaluate the "fraction"  a/(a+b)
    inline ValueWithError frac 
    ( const double          a , 
      const ValueWithError& b ) { return frac ( ValueWithError ( a ) , b ) ; }
    // ========================================================================
    /// evaluate the "asymmetry"  (a-b)/(a+b)
    inline ValueWithError asym 
    ( const ValueWithError& a , 
      const ValueWithError& b ) { return a.asym ( b ) ; }
    /// evaluate the "asymmetry"  (a-b)/(a+b)
    inline ValueWithError asym 
    ( const ValueWithError& a , 
      const double          b ) { return a.asym ( b ) ; }
    /// evaluate the "asymmetry"  (a-b)/(a+b)
    inline ValueWithError asym 
    ( const double          a , 
      const ValueWithError& b ) { return asym ( ValueWithError ( a ) , b ) ; }
    // ========================================================================
    /** evaluate abs(a) 
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
     *  \f$ \varepsilon = \frac{1}{ 1 + \frac{N_{rejected}}{N_accepted}}\f$ 
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
     *  \f[ R = \frac{N_acc}{N_tot} = \frac{N_acc}{N_acc+N_rej} = 
     *          \left( 1 + \frac{N_rej}{N_acc}\right)^{-1} \f]
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
     *  using jackknife method:
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
    /** evaluate pow(a,b)
     *  @param a (INPUT) the base 
     *  @param b (INPUT) the exponent 
     *  @return the <c>a</c> raised to power <c>b</c> 
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError pow 
    ( const ValueWithError& a , 
      const int             b ) ;
    // ========================================================================    
    /** evaluate pow(a,b)
     *  @param a (INPUT) the base 
     *  @param b (INPUT) the exponent 
     *  @return the <c>a</c> raised to power <c>b</c> 
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError pow 
    ( const ValueWithError& a , 
      const double          b ) ;
    // ========================================================================
    /** evaluate pow(a,b)
     *  @param a (INPUT) the base 
     *  @param b (INPUT) the exponent 
     *  @return the <c>a</c> raised to power <c>b</c> 
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError pow 
    ( const int             a , 
      const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate pow(a,b)
     *  @param a (INPUT) the base 
     *  @param b (INPUT) the exponent 
     *  @return the <c>a</c> raised to power <c>b</c> 
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError pow 
    ( const double          a , 
      const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate pow(a,b)
     *  @param a (INPUT) the base 
     *  @param b (INPUT) the exponent 
     *  @return the <c>a</c> raised to power <c>b</c> 
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError pow 
    ( const ValueWithError& a , 
      const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate exp(b)
     *  @param b (INPUT) the exponent 
     *  @return the <c>e</c> raised to power <c>b</c> 
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError exp
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate exp2(b)
     *  @param b (INPUT) the exponent 
     *  @return the <c>2</c> raised to power <c>b</c> 
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError exp2
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate expm1(b)
     *  @param b (INPUT) the exponent 
     *  @return  expm1(b) 
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError expm1
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate log(b)
     *  @param b (INPUT) the parameter 
     *  @return logarithm
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError log
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate log2(b)
     *  @param b (INPUT) the parameter 
     *  @return logarithm on base 2 
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError log2
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate log10(b)
     *  @param b (INPUT) the parameter 
     *  @return logarithm on base 10 
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError log10
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate log1p(b)
     *  @param b (INPUT) the parameter 
     *  @return  log1p(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError log1p
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate sqrt(b)
     *  @param b (INPUT) the parameter 
     *  @return  sqrt(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError sqrt
    ( const ValueWithError& b ) ;
    // ========================================================================
    /** evaluate 'signed sqrt' 
     *  @param a (INPUT) the value 
     *  @return the signed-sqrt  
     */
    ValueWithError signed_sqrt 
    ( const ValueWithError& a ) ;
    // ========================================================================    
    /** evaluate cbrt(b)
     *  @param b (INPUT) the parameter 
     *  @return  cbrt(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError cbrt
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate sin(b)
     *  @param b (INPUT) the parameter 
     *  @return  sin(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError sin 
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate cos(b)
     *  @param b (INPUT) the parameter 
     *  @return  cos(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError cos 
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate tan(b)
     *  @param b (INPUT) the parameter 
     *  @return  tan(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError tan 
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate sinh(b)
     *  @param b (INPUT) the parameter 
     *  @return  sinh(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError sinh 
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate cosh(b)
     *  @param b (INPUT) the parameter 
     *  @return  cosh(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError cosh 
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate tanh(b)
     *  @param b (INPUT) the parameter 
     *  @return  tanh(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError tanh 
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate sech(b)
     *  @param b (INPUT) the parameter 
     *  @return  sech(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError sech 
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate erf(b) the error function 
     *  @param b (INPUT) the parameter 
     *  @return  erf(b)
     *  @warning invalid and small covariances are ignored 
     *  @see https://en.wikipedia.org/wiki/Error_function
     */
    ValueWithError erf
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate erfc(b) complementary error function 
     *  @param b (INPUT) the parameter 
     *  @return  erfc(b)
     *  @warning invalid and small covariances are ignored 
     *  @see https://en.wikipedia.org/wiki/Error_function
     */
    ValueWithError erfc
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate erfi(b)  imaginary error function 
     *  @param b (INPUT) the parameter 
     *  @return  erfc(b)
     *  @warning invalid and small covariances are ignored 
     *  @see https://en.wikipedia.org/wiki/Error_function
     */
    ValueWithError erfi
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate erfcx(b) complementary scaled error function 
     *  @param b (INPUT) the parameter 
     *  @return  erfc(b)
     *  @warning invalid and small covariances are ignored 
     *  @see https://en.wikipedia.org/wiki/Error_function
     */
    ValueWithError erfcx
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate probit(b)
     *  @param b (INPUT) the parameter 
     *  @return  erfc(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError probit
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate asin(b)
     *  @param b (INPUT) the parameter 
     *  @return  asin(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError asin
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate acos(b)
     *  @param b (INPUT) the parameter 
     *  @return  acos(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError acos
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate atan(b)
     *  @param b (INPUT) the parameter 
     *  @return  atan(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError atan
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate atan2(y,x)
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
    /** evaluate asinh(b)
     *  @param b (INPUT) the parameter 
     *  @return  asinh(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError asinh
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate acosh(b)
     *  @param b (INPUT) the parameter 
     *  @return  acosh(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError acosh
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate atanh(b)
     *  @param b (INPUT) the parameter 
     *  @return  atanh(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError atanh
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate tgamma(b)
     *  @param b (INPUT) the parameter 
     *  @return  Gamma(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError tgamma
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate lgamma(b)
     *  @param b (INPUT) the parameter 
     *  @return  log(Gamma(b))
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError lgamma
    ( const ValueWithError& b ) ;
    // ========================================================================    
    /** evaluate igamma(b)
     *  @param b (INPUT) the parameter 
     *  @return  1/Gamma(b)
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError igamma
    ( const ValueWithError& b ) ;
    // ========================================================================
    /** evaluate Pochhammer symbol 
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
      const unsigned short n ) ;
    // ========================================================================
    /** evaluate standard Gauss PDF 
     *  @param x the value 
     *  @return valeu of the standard gaussian PDF  
     */
    ValueWithError gauss_pdf ( const ValueWithError& x ) ;
    // ========================================================================
    /** evaluate standard Gauss CDF 
     *  @param x the value 
     *  @return value of the standard gaussian CDF  
     */
    ValueWithError gauss_cdf  ( const ValueWithError& x ) ;    
    // ========================================================================    
    /** evaluate <code>hypot(x,y)</code>
     *  \f$ \sqrt( x^2 + y^2 ) \f$
     *   @param x (INPUT) the first parameter
     *   @param y (INPUT) the second parameter
     *   @param c (INPUT) the correlation coefficient  (-1<=c<=1)
     *   @return the value of <code>hypot</code> function
     */
    ValueWithError  hypot
    ( const ValueWithError& x     , 
      const ValueWithError& y     , 
      const double          c = 0 ) ;
    // ========================================================================    
    /** evaluate fma(x,y,z) = x*y+x 
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
    bool natural_number ( const ValueWithError& value ) ;
    // ========================================================================
    /** Does this object represent natural entry in histogram
     *  - non-negative integer value 
     *  - cov2 == value or ( 0 == value && 1 == cov2 ) 
     */
    bool natural_entry  ( const ValueWithError& value ) ;
    // ========================================================================
    /** make a sum two elements taking into account the 
     *  correlation coefficient  
     *  @param a  (input) the first value 
     *  @param b  (input) the second value 
     *  @param c  (input) the correlation coefficient
     *  @return a+b 
     *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
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
     *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
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
     *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
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
     *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
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
     *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
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
     *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
     *  @date 2012-11-09
     */
    ValueWithError  fraction
    ( const ValueWithError& a     , 
      const ValueWithError& b     , 
      const double          c = 0 ) ;
    // ========================================================================
    /** calculate "asymmetry" of two elements  \f$ \kappa = \frac{a-b}{a+b}\f$ 
     *  taking into account the correlation coefficient  
     *  @param a  (input) the first value 
     *  @param b  (input) the second value 
     *  @param c  (input) the correlation coefficient
     *  @return (a-b)/(a+b)
     *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
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
    ValueWithError abssum ( const std::vector<ValueWithError>& vct ) ;
    // ========================================================================    
    /** get the sum of the absolute values
     *  @param vct the vector
     *  @return the sum over the vector 
     */
    inline 
    ValueWithError sumabs ( const std::vector<ValueWithError>& vct )
    { return abssum ( vct )  ; }
    // ========================================================================    
    /** evaluate polynomial
     *  \f$f(x) = a_0 + a_1x + a_2x^2 + ... + a_{n-1}x^{n-1} + a_nx^n\f$
     *  such as \f$f(0) = a_0 \f$      
     *  using Horner rule
     *  @param poly  INPUT the coefficients
     *  @param x     INPUT argument 
     *  @return value of polynomial
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError horner_a0 ( const std::vector<double>& poly ,
                               const ValueWithError&      x    ) ;
    // ======================================================================
    /** evaluate polynomial
     *  \f$f(x) = a_0x^n + a_1x^{n-1}+ ... + a_{n-1}x + a_n\f$, 
     *  such as \f$f(0) = a_n \f$      
     *  using Horner rule
     *  @param poly  INPUT the coefficients 
     *  @param x     INPUT argument 
     *  @return value of polynomial
     *  @warning invalid and small covariances are ignored 
     */
    ValueWithError horner_aN ( const std::vector<double>& poly ,
                               const ValueWithError&      x    ) ;
    // ========================================================================    
    /// the output operator for the vector 
    std::ostream& operator<<
    ( std::ostream& s , 
      const std::vector<Ostap::Math::ValueWithError>& v ) ;
    // ========================================================================
    /// swap the content with another object
    inline void ValueWithError::swap ( ValueWithError& right ) 
    {
      std::swap ( m_value , right.m_value ) ;
      std::swap ( m_cov2  , right.m_cov2  ) ;
    }
    // ========================================================================
    /// swap the content of two objects 
    inline void swap ( ValueWithError& left , 
                       ValueWithError& right ) { left.swap ( right ) ; }
    // ========================================================================
    // converison to string
    inline std::string to_string ( const ValueWithError& v ) 
    { return v.toString() ; }
    // ========================================================================
  } //                                             end of namespace Ostap::Math 
  // ==========================================================================
} //                                                    end of namespace Ostap
// ============================================================================
// The END 
// ============================================================================
#endif // OSTAP_VALUEWITHERRORS_H
// ============================================================================
