// $Id$
// ============================================================================
#ifndef OSTAP_BERNSTEIN_H
#define OSTAP_BERNSTEIN_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <functional>
#include <vector>
#include <complex>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/NSphere.h"
#include "Ostap/Polynomials.h"
// ============================================================================
/** @file Ostap/Bernstein.h
 *  Set of useful math-functions, related to Bernstein polynomials
 *
 *  @see http://en.wikipedia.org/wiki/Bernstein_polynomial
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2010-04-19
 *
 *                    $Revision$
 *  Last modification $Date$
 *                 by $author$
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /// forward declaration
    class LegendreSum  ; // forward declaration
    class ChebyshevSum ; // forward declaration
    class Polynomial   ; // forward declaration
    // ========================================================================
    /** @class Bernstein
     *  The sum of bernstein's polynomial of order N
     *  \f$f(x) = \sum_i a_i B^n_i(x)\f$, where
     *  \f$ B^n_k(x) = C^n_k x^k(1-x)^{n-k}\f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     */
    class Bernstein : public Ostap::Math::PolySum
    {
      // ======================================================================
    public:
      // ======================================================================
      /** helper structure to denote the basic Bernstein polynomials B(k,N)
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       */
      class Basic
      {
      public:
        // ====================================================================
        explicit
          Basic ( const unsigned short k = 0 ,
                  const unsigned short N = 0 )
          : m_k ( k )
          , m_N ( N )
        {}
        // ====================================================================
      public :
        // ====================================================================
        unsigned short k () const { return m_k ; }
        unsigned short N () const { return m_N ; }
        // ====================================================================
      private:
        // ====================================================================
        unsigned short m_k ;
        unsigned short m_N ;
        // ====================================================================
      } ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      Bernstein ( const unsigned short        N     = 0 ,
                  const double                xmin  = 0 ,
                  const double                xmax  = 1 ) ;
      // ======================================================================
      /// constructor from N+1 coefficients
      Bernstein ( const std::vector<double>&  pars      ,
                  const double                xmin  = 0 ,
                  const double                xmax  = 1 ) ;
      // ======================================================================
      /// construct the basic bernstein polinomial  B(k,N)
      Bernstein  ( const Basic&              basic     ,
                   const double              xmin  = 0 ,
                   const double              xmax  = 1 ) ;
      // ======================================================================
      /// template constructor from sequence of parameters
      template <class ITERATOR>
        Bernstein ( ITERATOR                 first ,
                    ITERATOR                 last  ,
                    const double             xmin  ,
                    const double             xmax  )
        : Ostap::Math::PolySum ( first , last )
        , m_xmin ( std::min ( xmin, xmax ) )
        , m_xmax ( std::max ( xmin, xmax ) )
      {}
      // ======================================================================
      /// constructor  from Bernstein polynomial from *different* domain
      Bernstein ( const Bernstein& poly ,
                  const double     xmin ,
                  const double     xmax ) ;
      // ======================================================================
      /** construct Bernstein interpolant
       *  @param x    vector of abscissas
       *  @param y    vector of function values
       *  @param xmin low  edge for Bernstein polynomial
       *  @param xmin high edge for Bernstein polynomial
       *  - if vector of y is longer  than vector x, extra values are ignored
       *  - if vector of y is shorter than vector x, missing entries are assumed to be zero
       *  It relies on Newton-Bernstein algorithm
       *  @see http://arxiv.org/abs/1510.09197
       *  @see Mark Ainsworth and Manuel A. Sanches,
       *       "Computing of Bezier control points of Largangian interpolant
       *       in arbitrary dimension", arXiv:1510.09197 [math.NA]
       *  @see http://adsabs.harvard.edu/abs/2015arXiv151009197A
       */
      Bernstein ( const std::vector<double>& x         ,
                  const std::vector<double>& y         ,
                  const double               xmin  = 0 ,
                  const double               xmax  = 1 ) ;
      // ======================================================================
      /// copy
      Bernstein ( const Bernstein&  ) = default ;
      /// move
      Bernstein (       Bernstein&& ) = default ;
      // ======================================================================
      /** constructor from Legendre polynomial
       *  @see http://www.sciencedirect.com/science/article/pii/S0377042700003769 eq.20
       */
      explicit Bernstein ( const LegendreSum&  poly ) ;
      // ======================================================================
      /// constructor from Chebyshev polynomial
      explicit Bernstein ( const ChebyshevSum& poly ) ;
      // ======================================================================
      /// constructor from simple monomial form
      explicit Bernstein ( const Polynomial&   poly ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value of polynomial
      double evaluate ( const double x ) const ;
      // ======================================================================
      /// get the value
      double operator () ( const double x ) const
      { return x < m_xmin ? 0 : x > m_xmax ? 0 : evaluate ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get lower edge
      double xmin () const { return m_xmin ; }
      /// get upper edge
      double xmax () const { return m_xmax ; }
      // ======================================================================
    public:
      // ======================================================================
      /// is it a decreasing function?
      bool   decreasing    () const ;
      /// is it a increasing function?
      bool   increasing    () const ;
      /// is it a monothonical function?
      bool   monothonic    () const { return increasing() || decreasing() ; }
      /// is it a constant function?
      bool   constant      () const ;
      // ======================================================================
    public: // convert from local to global variables
      // ======================================================================
      double x ( const double t ) const
      { return       m_xmin   + ( m_xmax - m_xmin ) * t ; }
      double t ( const double x ) const
      { return ( x - m_xmin ) / ( m_xmax - m_xmin )     ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between xmin and xmax
      double    integral   () const ;
      /// get the integral between low and high
      double    integral   ( const double low , const double high ) const ;
      /** get indefinite integral  as function object
       *  \f$ I(x) = \int^{x}_{x_{min}} B(t) dt + C \f$
       *  @param C the integration constant
       */
      Bernstein indefinite_integral ( const double C = 0 ) const ;
      /// get the derivative at point "x"
      double    derivative          ( const double x     ) const ;
      /// get derivative as function object
      Bernstein derivative          () const ;
      // ======================================================================
    public :
      // ======================================================================
      /** elevate it:
       *  represent as Bernstein polynomial of order N+r
       *  @param r  INPUT increase of degree
       *  @return new polynomial of order N+r
       */
      Bernstein elevate  ( const unsigned short r ) const ;
      // ======================================================================
      /** reduce it:
       *  represent as Bernstein polynomial of order N-r
       *  @param r  INPUT decrease of degree
       *  @return new polynomial of order N-r
       */
      Bernstein reduce  ( const unsigned short r ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /** calculate ``nearest'' polynomial (in the sense of q-norm) of lower degree,
       *  where q-norm is defined as:
       *  \f$ \left| f \right|_{q} = \left( \sum_i \left|c_i\right|^q\right)^{\frac{1}{q}} \f$
       *
       *  - q_inv = 0.0 ->  \f$ max_k    \left|c_k\right|  \f$
       *  - q_inv = 0.5 ->  \f$ \sqrt{ \sum_k  c_k^2 }     \f$
       *  - q_inv = 1.0 ->  \f$ \sum_k \left| c_k \right|  \f$
       *  @see  N.Rezvani and R.M. Corless,
       *       "The Nearest Polynomial With A Given Zero, Revisited"
       *        ACM SIGSAM Bulletin, Vol. 39, No. 3, September 2005
       *  @see http://dl.acm.org/citation.cfm?doid=1113439.1113442
       */
      Bernstein nearest ( const double q_inv = 0 ) const ;
      // ======================================================================
      /** calculate q-norm of the polynomial
       *  where q-norm is defined as:
       *  \f$ \left| f \right|_{q} = \left( \sum_i \left|c_i\right|^q\right)^{\frac{1}{q}} \f$
       *
       *  - q_inv = 0.0 ->  \f$ max_k    \left|c_k\right|  \f$
       *  - q_inv = 0.5 ->  \f$ \sqrt{ \sum_k  c_k^2 }     \f$
       *  - q_inv = 1.0 ->  \f$ \sum_k \left| c_k \right|  \f$
       *
       *  q-norm for Bernstein polynomials is "equivalent" to common
       *  \f$L^2\f$-norm for generic functions,
       *  in the sense that
       *  \f$   \gamma_n \left| f \right|_{\inf}
       *             \le \left| f \right|_{\int}
       *             \le \left| f \right|_{\inf} \f$,
       *  where \f$  \left| f \right|_{\int}^2 = \int  \left| f \right|^2 dx \f$ and
       *  \f$ \gamma_n = \min_{j} \int \left|  B^n_j(x) \right|^2 dx
       *            \sim  (\pi n)^{-\frac{3}{4}} \f$
       *  similarly \f$ \left| f \right|_I \le \left| f \right|_1 \f$,
       */
      double    norm   ( const double q_inv = 0 ) const ;
      // ======================================================================
      /** filter out very small terms
       *  the term is considered to be very small if
       *   - it is numerically zero
       *   - or if epsilon > 0,
       *          abs ( c(k) * C(n,k) * k^k(n-k)^(n-k)/n^n ) < epsilon
       *  Since the maximum value for each term of
       *  \f$ c_k C^n_k \frac{ k^k (n-k)^{n-k}{ n^n}}\f$
       *  @param epsilon  parameter to define "smalness" of terms
       *  @returm number of nullified terms
       */
      unsigned short remove_noise ( const double epsilon = 0 ) ;
      // ======================================================================
      /** how close are two polynomials in q-norm?
       *  where q-norm is defined as:
       *  \f$ \left| f \right|_{q} = \left( \sum_i \left|c_i\right|^q\right)^{\frac{1}{q}} \f$
       *
       *  - q_inv = 0.0 -> \f$ max_k    \left|c_k\right|  \f$
       *  - q_inv = 0.5 -> \f$ \sqrt{ \sum_k  c_k^2 }     \f$
       *  - q_inv = 1.0 -> \f$ \sum_k \left| c_k \right|  \f$
       */
      double distance ( const Bernstein& other , const double q_inv = 0 ) const ;
      // ======================================================================
    public:  // polynomial division
      // ======================================================================
      /** polynomial division
       *  \f$  f(x) = q(z)*g(x) + r(x) \f$
       *  @return the pair q(x),r(x)
       */
      std::pair<Bernstein,Bernstein> divmod   ( const Bernstein& g ) const ;
      // ======================================================================
      /** polynomial division
       *  \f$  f(x) = q(z)*g(x) + r(x) \f$
       *  @return the quotient q(x)
       */
      // ======================================================================
      Bernstein                      quotient ( const Bernstein& g ) const ;
      /** polynomial division
       *  \f$  f(x) = q(z)*g(x) + r(x) \f$
       *  @return the reminder r(x)
       */
      Bernstein                      reminder ( const Bernstein& g ) const ;
      // ======================================================================
      /** get Greatest Common Divisor
       */
      Bernstein  gcd ( const Bernstein& b ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /** get the leading power coefficient
       *  \f$ f(x) = hx^n + ....\f$
       *  @return the coefficient at x^n
       */
      double    head  () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// simple  manipulations with polynoms: shift it!
      Bernstein& operator += ( const double a ) ;
      /// simple  manipulations with polynoms: shift it!
      Bernstein& operator -= ( const double a ) ;
      /// simple  manipulations with polynoms: scale it!
      Bernstein& operator *= ( const double a ) ;
      /// simple  manipulations with polynoms: scale it!
      Bernstein& operator /= ( const double a ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// negate it!
      Bernstein  operator-() const ;
      // ======================================================================
    public: // a bit of operators for python
      // ======================================================================
      /// Sum of Bernstein polynomial and a constant
      Bernstein __add__   ( const double value ) const ;
      /// Sum of Bernstein polynomial and a constant
      Bernstein __radd__  ( const double value ) const ;
      /// Product of Bernstein polynomial and a constant
      Bernstein __mul__   ( const double value ) const ;
      /// Product of Bernstein polynomial and a constant
      Bernstein __rmul__  ( const double value ) const ;
      /// Subtract a constant from Benrstein polynomial
      Bernstein __sub__   ( const double value ) const ;
      /// Constant minus Bernstein polynomial
      Bernstein __rsub__  ( const double value ) const ;
      /// Divide Benrstein polynomial by a constant
      Bernstein __div__   ( const double value ) const ;
      /// Negate Bernstein polynomial
      Bernstein __neg__   () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// Sum of   Bernstein polynomials
      Bernstein __add__   ( const Bernstein& other ) const ;
      /// Subtract Bernstein polynomials
      Bernstein __sub__   ( const Bernstein& other ) const ;
      /// Multiply Bernstein polynomials
      Bernstein __mul__   ( const Bernstein& other ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the underlying Bernstein polynomial  (self here)
      const Ostap::Math::Bernstein& bernstein () const { return *this ; }
      // ======================================================================
    public:
      // ======================================================================
      /// the sum two Bernstein polynomials
      Bernstein  sum      ( const Bernstein&        other ) const ;
      /// subtract Bernstein polynomials
      Bernstein  subtract ( const Bernstein&        other ) const ;
      /// multiply Bernstein polynomials
      Bernstein  multiply ( const Bernstein&        other ) const ;
      /// multiply Bernstein polynomials with the basic bernstein polynomial
      Bernstein  multiply ( const Bernstein::Basic& other ) const ;
      /** multiply Bernstein polynomial with
       *  \f$ (x-x_{min})^i(x_{max}-x)^j \f$
       */
      Bernstein  multiply ( const unsigned short i ,
                            const unsigned short j ) const ;
      /// power function
      Bernstein  pow      ( const unsigned short i ) const ;
      // ======================================================================
    public:  // various assignements
      // ======================================================================
      /// copy assignement
      Bernstein& operator=( const Bernstein&  right ) ;
      /// move assignement
      Bernstein& operator=(       Bernstein&& right ) ;
      /// assignement from the constant
      Bernstein& operator=( const double      right ) ;
      // ======================================================================
    private:  // internal data
      // ======================================================================
      /// the left edge of interval
      double m_xmin  ;                             // the left edge of interval
      /// the right edge of interval
      double m_xmax  ;                             // the right edge of interval
      // ======================================================================
    };
    // ========================================================================
    ///  Bernstein plus      constant
    inline Bernstein operator+( const Bernstein& p , const double v )
    { return Bernstein ( p ) += v ; } //  Bernstein plus constant
    ///  Bernstein multiply  constant
    inline Bernstein operator*( const Bernstein& p , const double v )
    { return Bernstein ( p ) *= v ; } //  Bernstein plus constant
    ///  Bernstein minus constant
    inline Bernstein operator-( const Bernstein& p , const double v )
    { return Bernstein ( p ) -= v ; } //  Bernstein plus constant
    ///  Bernstein divide constant
    inline Bernstein operator/( const Bernstein& p , const double v )
    { return Bernstein ( p ) /= v ; } //  Bernstein plus constant
    ///  Constant plus  Bernstein
    inline Bernstein operator+( const double v , const Bernstein& p ) { return p +   v  ; }
    ///  Constant times Bernstein
    inline Bernstein operator*( const double v , const Bernstein& p ) { return p *   v  ; }
    ///  Constant minus Bernstein
    inline Bernstein operator-( const double v , const Bernstein& p ) { return v + (-p) ; }
    // ========================================================================
    ///  Bernstein plus     Bernstein
    inline Bernstein operator+( const Bernstein& a , const Bernstein& b  )
    { return a.sum      ( b ) ; } //  Bernstein plus     Bernstein
    ///  Bernstein minus    Bernstein
    inline Bernstein operator-( const Bernstein& a , const Bernstein& b  )
    { return a.subtract ( b ) ; } //  Bernstein subtract Bernstein
    ///  Bernstein multiply Bernstein
    inline Bernstein operator*( const Bernstein& a , const Bernstein& b  )
    { return a.multiply ( b ) ; } //  Bernstein multiply  Bernstein
    ///  polynomial division: quotient
    inline Bernstein operator/( const Bernstein& a , const Bernstein& b  )
    { return a.quotient ( b ) ; } // polynomial division: quotient
    ///  polynomial division: reminder
    inline Bernstein operator%( const Bernstein& a , const Bernstein& b  )
    { return a.reminder ( b ) ; } // polynomial division: reminder
    // ========================================================================
    /** polynomial division
     *  Return pair of polynomials q and r, such as
     *  \f$ a = qb + r \f$
     *  @param a
     *  @param b
     */
    inline
    std::pair<Bernstein,Bernstein> divmod
    ( const Bernstein& a ,
      const Bernstein& b ) { return a.divmod ( b ) ; }
    // ========================================================================
    /** get the integral between low and high for a product of Bernstein
     *  polynom and the exponential function with the exponent tau
     *  \f[  \int_{a}^{b} \mathcal{B} e^{\tau x } \mathrm{d}x \f]
     *  @param poly  bernstein polynomial
     *  @param tau   slope parameter for exponential
     *  @param a     low  integration range
     *  @param b     high integration range
     */
    double integrate
    ( const Ostap::Math::Bernstein& poly ,
      const double                  tau  ,
      const double                  a    ,
      const double                  b    ) ;
    // ========================================================================
    /** get the integral between 0 and 1 for a product of basic  Bernstein
     *  polynom and the exponential function with the exponent tau
     *  \f[  \int_{0}^{1} \mathcal{B} e^{\tau x } \mathrm{d}x \f]
     *  @param b     basic bernstein polynomial
     *  @param tau   slope parameter for exponential
     */
    double integrate
    ( const Ostap::Math::Bernstein::Basic& b    ,
      const double                         tau  ) ;
    // =======================================================================
    /** get the integral between \f$x_{min}\f$ and \f$x_{max}\f$ for
     *  a product of Bernstein polynom and the exponential function
     *   with the exponent tau
     *  \f[  \int_{x_{min}}^{x_{max}} \mathcal{B} e^{\tau x } \mathrm{d}x \f]
     *  @param poly  bernstein polynomial
     *  @param tau   slope parameter for exponential
     */
    double integrate
    ( const Ostap::Math::Bernstein& poly ,
      const double                  tau  ) ;
    // ========================================================================
    /** get the integral between 0 and 1 for a product of basic  Bernstein
     *  polynom and monomial or degree m
     *  \f[  \int_{0}^{1} \mathcal{B} \frac{x^m}{m!} \mathrm{d}x \f]
     *  @param b     basic bernstein polynomial
     *  @param m     degree of monomial
     */
    double integrate_poly
    ( const Ostap::Math::Bernstein::Basic& b ,
      const unsigned short                 m ) ;
    // =======================================================================
    /** get the integral between xmin and xmax Bernstein
     *  polynom and monomial or degree m
     *  \f[  \int_{x_min}^{x_max} \mathcal{B} \frac{(x-x_min)^m}{m!} \mathrm{d}x \f]
     *  @param b     basic bernstein polynomial
     *  @param m     degree of monomial
     */
    double integrate_poly
    ( const Ostap::Math::Bernstein& b ,
      const unsigned short          m ) ;
    // ========================================================================
    /** get the integral between xmin and xmax Bernstein
     *  polynom and monomial or degree m
     *  \f[  \int_{low}^{high} \mathcal{B} \frac{(x-x_min)^m}{m!} \mathrm{d}x \f]
     *  @param b     basic bernstein polynomial
     *  @param m     degree of monomial
     *  @param low   low  integration limit
     *  @param high  high integtation limit
     */
    double integrate_poly
    ( const Ostap::Math::Bernstein& b    ,
      const unsigned short          m    ,
      const double                  low  ,
      const double                  high ) ;
    // =======================================================================
    /** de Casteljau algorithm for summation of Bernstein polynomials
     *  \f$ f(x) = \sum_i p_i B_ik(x) \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-02-10
     */
    double casteljau
    ( const std::vector<double>& pars ,
      const double               x    ) ;
    // ========================================================================
    /** Dual basic Bernstein function
     *  The dual basic functions \f$ d^n_j(x)\f$ are dedeined as
     *   \f$  \int_{x_{min}}^{x_{max}}   b^n_k(x) d^n_j(x) = \delta_{kj}\f$,
     *   where \f$b^n_k(x)\f$ is basic Bernstein polynomial
     *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
     *  @date 2016-07-03
     */
    class BernsteinDualBasis : public std::unary_function<double,double>
    {
      // ======================================================================
    public :
      // ======================================================================
      unsigned short k () const { return m_k                  ; }
      unsigned short N () const { return m_bernstein.degree() ; }
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      BernsteinDualBasis ( const unsigned short N     = 0 ,
                           const unsigned short k     = 0 ) ;
      /// copy constructor
      BernsteinDualBasis  ( const BernsteinDualBasis&  right ) ;
      /// cconstructor
      BernsteinDualBasis  (       BernsteinDualBasis&& right ) ;
      /// destructor
      ~BernsteinDualBasis () ;
      // ======================================================================
    public:
      // ======================================================================
      /// calculate the value of dual bernstein function
      double operator() ( const double x ) const
      { return m_k <= N() ? m_bernstein ( x ) : 0.0 ; }
      // ======================================================================
    public:
      // ======================================================================
      ///  get the parameters
      double par ( const unsigned short i ) const { return m_bernstein.par( i ) ; }
      ///  get all parameters
      const std::vector<double>& pars()     const { return m_bernstein.pars() ; }
      // ======================================================================
    private:
      // ======================================================================
      /// the index
      unsigned short   m_k         ; // the index
      /// the actual bernstein polynomial
      Bernstein        m_bernstein ; // the actual bernstein polynomial
      // ======================================================================
    };
    // ========================================================================
    /** @class BernsteinEven
     *  A special case of BErnstein polynomial with symmetry:
     *  \f$ f( \frac{x_{max}+x_{min}}{2} - x ) \equiv  \frac{x_{max}+x_{min}}{2} + x ) \f$
     *  @see Ostap::Math::Bernstein
     *  @author Vanya Belyaev Ivan.Belyaev@iep.ru
     *  @date 2016-10-02
     */
    class BernsteinEven: public std::unary_function<double,double>
    {
    public:
      // ======================================================================
      /** constructor
       *  the actual degree of polynomial will be 2*N
       *  @param N  parameter that defiend the order of polynomial (2*N)
       *  @param xmin low edge
       *  @param xmax high edge
       */
      BernsteinEven ( const unsigned short N     = 0 ,
                      const double         xmin  = 0 ,
                      const double         xmax  = 1 ) ;
      /** constructor from list of coefficients
       *  @param xmin low edge
       *  @param xmax high edge
       */
      BernsteinEven ( const std::vector<double>& pars      ,
                      const double               xmin  = 0 ,
                      const double               xmax  = 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value of polynomial
      double evaluate ( const double x ) const { return m_bernstein.evaluate ( x ) ; }
      // ======================================================================
      /// get the value
      double operator () ( const double x ) const { return m_bernstein ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// the effective degree of polynomial
      unsigned short degree() const { return 2*m_N   ; }
      /// number of parameters
      unsigned short npars () const { return m_N + 1 ; }
      /// all zero ?
      bool                 zero  () const { return m_bernstein.zero() ; }
      /** set k-parameter
       *  @param k index
       *  @param value new value
       *  @return true if parameter is actually changed
       */
      bool setPar          ( const unsigned short k , const double value ) ;
      /** set k-parameter
       *  @param k index
       *  @param value new value
       *  @return true iof parameter is actually changed
       */
      bool setParameter    ( const unsigned short k , const double value )
      { return setPar      ( k , value ) ; }
      /// get the parameter value
      double  par          ( const unsigned short k ) const
      { return  k < m_N  ? m_bernstein.par ( k )  : 0.0 ; }
      /// get the parameter value
      double  parameter    ( const unsigned short k ) const { return par ( k ) ; }
      /// get all parameters (by value!!! COPY!!)
      std::vector<double> pars () const ;
      // ======================================================================
    public: // convert from local to global variables
      // ======================================================================
      double x ( const double t ) const { return m_bernstein.x ( t ) ; }
      double t ( const double x ) const { return m_bernstein.t ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get lower edge
      double xmin () const { return m_bernstein.xmin () ; }
      /// get upper edge
      double xmax () const { return m_bernstein.xmax () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// simple  manipulations with polynoms: shift it!
      BernsteinEven& operator += ( const double a ) { m_bernstein += a ; return *this ; }
      /// simple  manipulations with polynoms: shift it!
      BernsteinEven& operator -= ( const double a ) { m_bernstein -= a ; return *this ; }
      /// simple  manipulations with polynoms: scale it!
      BernsteinEven& operator *= ( const double a ) { m_bernstein *= a ; return *this ; }
      /// simple  manipulations with polynoms: scale it!
      BernsteinEven& operator /= ( const double a ) { m_bernstein /= a ; return *this ; }
      // ======================================================================
    public: // a bit of operators for python
      // ======================================================================
      /// Sum of Bernstein polynomial and a constant
      BernsteinEven __add__   ( const double value ) const ;
      /// Sum of Bernstein polynomial and a constant
      BernsteinEven __radd__  ( const double value ) const ;
      /// Product of Bernstein polynomial and a constant
      BernsteinEven __mul__   ( const double value ) const ;
      /// Product of Bernstein polynomial and a constant
      BernsteinEven __rmul__  ( const double value ) const ;
      /// Subtract a constant from Benrstein polynomial
      BernsteinEven __sub__   ( const double value ) const ;
      /// Constant minus Bernstein polynomial
      BernsteinEven __rsub__  ( const double value ) const ;
      /// Divide Benrstein polynomial by a constant
      BernsteinEven __div__   ( const double value ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between xmin and xmax
      double    integral   () const { return m_bernstein.integral() ; }
      /// get the integral between low and high
      double    integral   ( const double low , const double high ) const
      { return m_bernstein.integral( low , high ) ; }
      /** get indefinite integral  as function object
       *  \f$ I(x) = \int^{x}_{x_{min}} B(t) dt + C \f$
       *  @param C the integration constant
       */
      Bernstein indefinite_integral ( const double C = 0 ) const
      { return m_bernstein.indefinite_integral( C ) ; }
      /// get the derivative at point "x"
      double    derivative          ( const double x     ) const
      { return m_bernstein.derivative ( x ) ; }
      /// get derivative as function object
      Bernstein derivative          () const
      { return m_bernstein.derivative () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// convert to normal bernstein polynomial
      const Bernstein& bernstein() const { return m_bernstein ; }
      /// convert to normal bernstein polynomial
      operator const Bernstein&()  const { return m_bernstein ; }
      // ======================================================================
    private:
      // ======================================================================
      /// the half-order
      unsigned short m_N         ;
      /// the actual Bernstein polynomial
      Bernstein      m_bernstein ; // the actual Bernstein polynomial
      // ======================================================================
    } ;
    ///  Bernstein plus      constant
    inline BernsteinEven operator+( const BernsteinEven& p , const double v )
    { return BernsteinEven ( p ) += v ; } //  Bernstein plus constant
    ///  Bernstein multiply  constant
    inline BernsteinEven operator*( const BernsteinEven& p , const double v )
    { return BernsteinEven ( p ) *= v ; } //  Bernstein plus constant
    ///  Bernstein minus constant
    inline BernsteinEven operator-( const BernsteinEven& p , const double v )
    { return BernsteinEven ( p ) -= v ; } //  Bernstein plus constant
    ///  Bernstein divide constant
    inline BernsteinEven operator/( const BernsteinEven& p , const double v )
    { return BernsteinEven ( p ) /= v ; } //  Bernstein plus constant
    ///  Constant plus  Bernstein
    inline BernsteinEven operator+( const double v , const BernsteinEven& p )
    { return p +   v  ; }
    ///  Constant times Bernstein
    inline BernsteinEven operator*( const double v , const BernsteinEven& p )
    { return p *   v  ; }
    ///  Constant minus Bernstein
    inline BernsteinEven operator-( const double v , const BernsteinEven& p )
    { return v + -1*p; }
    // ========================================================================
    /** @class Positive
     *  The "positive" polynomial of order N
     *  Actually it is a sum of basic bernstein polynomials with
     *  non-negative coefficients
     *  \f$ f(x) = \sum_i  \alpha^2_i B^n_i(x)\f$, where
     *  \f$  \sum_i \alpha^2_i = 1\f$ and
     *  parameters \f$ \alpha_i(p_0,p_1,....p_{n-1})\f$ form
     *  n-dimension sphere
     */
    class Positive : public std::unary_function<double,double>
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      Positive ( const unsigned short        N     =  1 ,
                 const double                xmin  =  0 ,
                 const double                xmax  =  1 ) ;
      // ======================================================================
      /// constructor from N phases
      Positive ( const std::vector<double>&  phases     ,
                 const double                xmin  =  0 ,
                 const double                xmax  =  1 ) ;
      /// constructor from the sphere with coefficients
      Positive ( const Ostap::Math::NSphere& sphere    ,
                 const double                xmin = 0  ,
                 const double                xmax = 0  ) ;
      // ======================================================================
      /// copy
      Positive ( const Positive&  right ) ;
      /// move
      Positive (       Positive&& right ) ;
      // ======================================================================
      virtual ~Positive () ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double operator () ( const double x ) const { return m_bernstein ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get number of parameters
      std::size_t npars () const { return m_sphere.nPhi () ; }
      /// set k-parameter
      bool setPar       ( const unsigned short k , const double value ) ;
      /// set k-parameter
      bool setParameter ( const unsigned short k , const double value )
      { return setPar   ( k , value ) ; }
      /// get the parameter value
      double  par       ( const unsigned short k ) const
      { return m_sphere.par ( k ) ; }
      /// get all parameters (phases on sphere)
      const std::vector<double>& pars  () const { return m_sphere   .pars () ; }
      /// get bernstein coefficients
      const std::vector<double>& bpars () const { return m_bernstein.pars () ; }
      // ======================================================================
    public:  // some characteristics
      // ======================================================================
      /// degree
      unsigned short degree      () const { return m_bernstein.degree() ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the parameter value
      double  parameter ( const unsigned short k ) const { return par ( k ) ; }
      /// get lower edge
      double xmin () const { return m_bernstein.xmin () ; }
      /// get upper edge
      double xmax () const { return m_bernstein.xmax () ; }
      /// transform variables
      double x ( const double t ) const { return m_bernstein. x ( t )  ; }
      double t ( const double x ) const { return m_bernstein. t ( x )  ; }
      // ======================================================================
    public:
      //
      bool decreasing  () const { return m_bernstein.decreasing()     ; }
      bool increasing  () const { return m_bernstein.increasing()     ; }
      bool monothonic  () const { return increasing() || decreasing() ; }
      bool constant    () const { return m_bernstein.constant()       ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between xmin and xmax
      double integral   () const ;
      /// get the integral between low and high
      double integral   ( const double low , const double high ) const ;
      /// get the derivative
      double derivative ( const double x ) const
      { return m_bernstein.derivative ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the underlying Bernstein polynomial
      const Ostap::Math::Bernstein& bernstein () const { return m_bernstein ; }
      /// get the parameter sphere
      const Ostap::Math::NSphere&   sphere    () const { return m_sphere    ; }
      /// get the indefinite integral
      Bernstein indefinite_integral ( const double C = 0 ) const
      { return m_bernstein.indefinite_integral ( C ) ; }
      /// get the derivative
      Bernstein derivative          () const
      { return m_bernstein.derivative          () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// copy assignement
      Positive& operator=( const Positive&  right ) ;
      /// move assignement
      Positive& operator=(       Positive&& right ) ;
      // ======================================================================
    public:  /// basic operations  for python
      // ======================================================================
      /// Sum of Bernstein polynomial and a constant
      Bernstein __add__   ( const double value ) const { return m_bernstein + value ; }
      /// Sum of Bernstein polynomial and a constant
      Bernstein __radd__  ( const double value ) const { return m_bernstein + value ; }
      /// Product of Bernstein polynomial and a constant
      Bernstein __mul__   ( const double value ) const { return m_bernstein * value ; }
      /// Product of Bernstein polynomial and a constant
      Bernstein __rmul__  ( const double value ) const { return m_bernstein * value ; }
      /// Subtract a constant from Benrstein polynomial
      Bernstein __sub__   ( const double value ) const { return m_bernstein - value ; }
      /// Constant minus Bernstein polynomial
      Bernstein __rsub__  ( const double value ) const { return value - m_bernstein ; }
      /// Divide Bernstein polynomial by a constant
      Bernstein __div__   ( const double value ) const { return m_bernstein / value ; }
      /// Negate Bernstein polynomial
      Bernstein __neg__   () const { return -m_bernstein ; }
      // ======================================================================
    protected:
      // ======================================================================
      /// update bernstein coefficiencts
      virtual bool updateBernstein () ;
      // ======================================================================
    protected:
      // ======================================================================
      /// the actual bernstein polynomial
      Ostap::Math::Bernstein m_bernstein ; // the actual bernstein polynomial
      /// arameters sphere
      Ostap::Math::NSphere   m_sphere    ;
      // ======================================================================
    } ;
    // ========================================================================
    ///  Positive plus      constant
    inline Bernstein operator+( const Positive& p , const double v )
    { return p.bernstein() + v ; }
    ///  Positive multiply  constant
    inline Bernstein operator*( const Positive& p , const double v )
    { return p.bernstein() * v ; }
    ///  Positive minus     constant
    inline Bernstein operator-( const Positive& p , const double v )
    { return p.bernstein() - v ; }
    ///  Positive divide constant
    inline Bernstein operator/( const Positive& p , const double v )
    { return p.bernstein() / v ; }
    ///  Constant plus  Positive
    inline Bernstein operator+( const double v , const Positive& p ) { return p + v  ; }
    ///  Constant times Positive
    inline Bernstein operator*( const double v , const Positive& p ) { return p * v  ; }
    ///  Constant minus Positive
    inline Bernstein operator-( const double v , const Positive& p )
    { return v - p.bernstein() ; }
    // ========================================================================
    /** @class PositiveEven
     *  The "positive" polynomial of order N, symmetric as
     *  \f$ f( \frac{x_{max}+x_{min}}{2} - x ) \equiv  \frac{x_{max}+x_{min}}{2} + x ) \f$
     *  Actually it is a sum of basic bernstein polynomials with
     *  non-negative coefficients
     *  \f$ f(x) = \sum_i  \alpha^2_i B^n_i(x)\f$, where
     *  \f$  \sum_i \alpha^2_i = 1\f$ and
     *  parameters \f$ \alpha_i(p_0,p_1,....p_{n-1})\f$ form
     *  n-dimension sphere
     *  @see Ostap::Math::Bernstein
     *  @see Ostap::Math::BernsteinEven
     *  @see Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2016-10-02
     */
    class PositiveEven : public std::unary_function<double,double>
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      PositiveEven ( const unsigned short        N     =  1 ,
                     const double                xmin  =  0 ,
                     const double                xmax  =  1 ) ;
      // ======================================================================
      /// constructor from N phases
      PositiveEven ( const std::vector<double>&  phases     ,
                     const double                xmin  =  0 ,
                     const double                xmax  =  1 ) ;
      /// constructor from the sphere with coefficients
      PositiveEven ( const Ostap::Math::NSphere& sphere    ,
                     const double                xmin = 0  ,
                     const double                xmax = 0  ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double operator () ( const double x ) const { return m_even           ( x ) ; }
      /// get the value
      double evaluate    ( const double x ) const { return m_even.evaluate  ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get number of parameters
      std::size_t npars () const { return m_sphere.nPhi () ; }
      /// set k-parameter
      bool setPar       ( const unsigned short k , const double value ) ;
      /// set k-parameter
      bool setParameter ( const unsigned short k , const double value )
      { return setPar   ( k , value ) ; }
      /// get the parameter value
      double  par       ( const unsigned short k ) const
      { return m_sphere.par ( k ) ; }
      /// get all parameters (phases on sphere)
      const std::vector<double>& pars  () const { return m_sphere   .pars () ; }
      /// get bernstein coefficients (by value, copy)
      std::vector<double>        bpars () const { return m_even     .pars () ; }
      // ======================================================================
    public:  // some characteristics
      // ======================================================================
      /// degree
      unsigned short degree      () const { return m_even.degree() ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the parameter value
      double  parameter ( const unsigned short k ) const { return par ( k ) ; }
      /// get lower edge
      double xmin () const { return m_even.xmin () ; }
      /// get upper edge
      double xmax () const { return m_even.xmax () ; }
      /// transform variables
      double x ( const double t ) const { return m_even. x ( t )  ; }
      double t ( const double x ) const { return m_even. t ( x )  ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between xmin and xmax
      double integral   () const ;
      /// get the integral between low and high
      double integral   ( const double low , const double high ) const ;
      /// get the derivative
      double derivative ( const double x ) const
      { return m_even.derivative ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the underlying Bernstein Even polynomial
      const Ostap::Math::BernsteinEven& bernsteinEven () const { return m_even ; }
      const Ostap::Math::BernsteinEven& even          () const { return m_even ; }
      /// get the underlying Bernstein polynomial
      const Ostap::Math::Bernstein& bernstein () const { return m_even.bernstein() ; }
      /// get the parameter sphere
      const Ostap::Math::NSphere&   sphere    () const { return m_sphere    ; }
      /// get the indefinite integral
      Bernstein indefinite_integral ( const double C = 0 ) const
      { return m_even.indefinite_integral ( C ) ; }
      /// get the derivative
      Bernstein derivative          () const
      { return m_even.derivative          () ; }
      // ======================================================================
    public:  /// basic operations  for python
      // ======================================================================
      /// Sum of Bernstein polynomial and a constant
      BernsteinEven __add__   ( const double value ) const { return m_even + value  ; }
      /// Sum of Bernstein polynomial and a constant
      BernsteinEven __radd__  ( const double value ) const { return m_even + value  ; }
      /// Product of Bernstein polynomial and a constant
      BernsteinEven __mul__   ( const double value ) const { return m_even * value  ; }
      /// Product of Bernstein polynomial and a constant
      BernsteinEven __rmul__  ( const double value ) const { return m_even * value  ; }
      /// Subtract a constant from Benrstein polynomial
      BernsteinEven __sub__   ( const double value ) const { return m_even - value  ; }
      /// Constant minus Bernstein polynomial
      BernsteinEven __rsub__  ( const double value ) const { return value  - m_even ; }
      /// Divide Bernstein polynomial by a constant
      BernsteinEven __div__   ( const double value ) const { return m_even / value  ; }
      // ======================================================================
    protected:
      // ======================================================================
      /// update bernstein coefficients
      bool updateBernstein () ;
      // ======================================================================
    protected:
      // ======================================================================
      /// the actual bernstein polynomial
      Ostap::Math::BernsteinEven m_even   ; // the actual bernstein polynomial
      /// arameters sphere
      Ostap::Math::NSphere       m_sphere ;
      // ======================================================================
    } ;
    ///  Positive plus      constant
    inline BernsteinEven operator+( const PositiveEven& p , const double v )
    { return p.even() + v ; }
    ///  Positive multiply  constant
    inline BernsteinEven operator*( const PositiveEven& p , const double v )
    { return p.even() * v ; }
    ///  Positive minus     constant
    inline BernsteinEven operator-( const PositiveEven& p , const double v )
    { return p.even() - v ; }
    ///  Positive divide constant
    inline BernsteinEven operator/( const PositiveEven& p , const double v )
    { return p.even() / v ; }
    ///  Constant plus  Positive
    inline BernsteinEven operator+( const double v , const PositiveEven& p )
    { return p + v  ; }
    ///  Constant times Positive
    inline BernsteinEven operator*( const double v , const PositiveEven& p )
    { return p * v  ; }
    ///  Constant minus Positive
    inline BernsteinEven operator-( const double v , const PositiveEven& p )
    { return v - p.even() ; }
    // ========================================================================
    /** @class Monothonic
     *  The "positive" monothonic polynomial of order N
     *  Actually it is a sum of basic bernstein polynomials with
     *  non-negative coefficients
     */
    class Monothonic : public Ostap::Math::Positive
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      Monothonic
        ( const unsigned short       N          =    1 ,
          const double               xmin       =    0 ,
          const double               xmax       =    1 ,
          const bool                 increasing = true ) ;
      // ======================================================================
      /// constructor from N phases
      Monothonic
        ( const std::vector<double>& pars              ,
          const double               xmin       =    0 ,
          const double               xmax       =    1 ,
          const bool                 increasing = true ) ;
      // ======================================================================
      /// constructor positive spline
      Monothonic
        ( const Positive&            poly              ,
          const bool                 increasing        ) ;
      // ======================================================================
      /// copy  constructor
      Monothonic ( const Monothonic&  right ) ;
      /// move
      Monothonic (       Monothonic&& right ) = default ;
      // ======================================================================
      virtual ~Monothonic() ;
      // ======================================================================
    public:
      // ======================================================================
      /// increasing ?
      bool increasing () const { return  m_increasing  ; }
      /// decreasing ?
      bool decreasing () const { return !increasing () ; }
      /// monothonic
      bool monothonic () const { return  true  ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the minimal value of function
      double fun_min () const ; // get the minimal value of function
      /// get the maximal value of function
      double fun_max () const ; // get the maximal value of function
      // ======================================================================
    protected:
      // ======================================================================
      /// update bernstein coefficients
      bool updateBernstein () override;
      // ======================================================================
    protected:
      // ======================================================================
      /// increasing ?
      bool                   m_increasing ; // increasing ?
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Convex
     *  The "positive" polynomial of order N with
     *  fixed sign of first and second derivatives
     *  Actually it is a sum of basic bernstein polynomials with
     *  non-negative coefficients
     */
    class Convex : public Ostap::Math::Monothonic
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      Convex
        ( const unsigned short       N          =    1 ,
          const double               xmin       =    0 ,
          const double               xmax       =    1 ,
          const bool                 increasing = true ,
          const bool                 convex     = true ) ;
      // ======================================================================
      /// constructor from N phases
      Convex
        ( const std::vector<double>& pars              ,
          const double               xmin       =    0 ,
          const double               xmax       =    1 ,
          const bool                 increasing = true ,
          const bool                 convex     = true ) ;
      // ======================================================================
      /// constructor from polynom
      Convex
        ( const Positive&            poly       ,
          const bool                 increasing ,
          const bool                 convex     ) ;
      // ======================================================================
      /// constructor from polynom
      Convex
        ( const Monothonic&          poly      ,
          const bool                 convex    ) ;
      // ======================================================================
      /// copy constructor
      Convex ( const Convex&         right     ) ;
      /// move
      Convex (       Convex&&        right     ) = default ;
      // ======================================================================
      virtual ~Convex() ;
      // ======================================================================
    public:
      // ======================================================================
      /// convex     ?
      bool   convex    () const { return  m_convex    ; }
      /// convex     ?
      bool   concave   () const { return   !convex () ; }
      // ======================================================================
    protected:
      // ======================================================================
      /// update bernstein coefficients
      bool updateBernstein () override;
      // ======================================================================
    protected:
      // ======================================================================
      /// convex ?
      bool                   m_convex     ; // iconvex ?
      // ======================================================================
    } ;
    // ========================================================================

    // ========================================================================
    /** @class ConvexOnly
     *  The "positive" polynomial of order N with
     *  fixed sign the second derivatives
     *  Actually it is a sum of basic bernstein polynomials with
     *  non-negative coefficients
     */
    class ConvexOnly : public Ostap::Math::Positive
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      ConvexOnly
        ( const unsigned short       N          =    1 ,
          const double               xmin       =    0 ,
          const double               xmax       =    1 ,
          const bool                 convex     = true ) ;
      // ======================================================================
      /// constructor from N phases
      ConvexOnly
        ( const std::vector<double>& pars              ,
          const double               xmin       =    0 ,
          const double               xmax       =    1 ,
          const bool                 convex     = true ) ;
      // ======================================================================
      /// constructor from polynom
      ConvexOnly
        ( const Positive&            poly       ,
          const bool                 convex     ) ;
      // ======================================================================
      /// copy constructor
      ConvexOnly ( const ConvexOnly&   right     ) ;
      /// move
      ConvexOnly (       ConvexOnly&&  right     ) = default ;
      // ======================================================================
      virtual ~ConvexOnly() ;
      // ======================================================================
    public:
      // ======================================================================
      /// convex     ?
      bool   convex    () const { return  m_convex    ; }
      /// convex     ?
      bool   concave   () const { return   !convex () ; }
      // ======================================================================
    protected:
      // ======================================================================
      /// update bernstein coefficients
      bool updateBernstein () override;
      // ======================================================================
    protected:
      // ======================================================================
      /// convex ?
      bool                   m_convex     ; // iconvex ?
      // ======================================================================
    } ;
    // ========================================================================

    // ========================================================================
    // 2D-models
    // ========================================================================
    /** @class Bernstein2D
     *  The Bernstein's polynomial of order Nx*Ny
     */
    class Bernstein2D : public std::binary_function<double,double,double>
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      Bernstein2D ( const unsigned short       nX    =  1 ,
                    const unsigned short       nY    =  1 ,
                    const double               xmin  =  0 ,
                    const double               xmax  =  1 ,
                    const double               ymin  =  0 ,
                    const double               ymax  =  1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double operator () ( const double x ,
                           const double y ) const ;
      // ======================================================================
    public: // setters
      // ======================================================================
      /// set k-parameter
      bool setPar       ( const unsigned int   k     ,
                          const double         value ) ;
      /// set k-parameter
      bool setParameter ( const unsigned int   k     ,
                          const double         value )
      { return ( k < m_pars.size() ) && setPar ( k , value ) ; }
      /// set (l,m)-parameter
      bool setPar       ( const unsigned short l     ,
                          const unsigned short m     ,
                          const double         value ) ;
      /// set (l,m)-parameter
      bool setParameter ( const unsigned short l     ,
                          const unsigned short m     ,
                          const double         value )
      { return setPar   ( l , m  , value ) ; }
      // ======================================================================
    public: // getters
      // ======================================================================
      /// get (l,m)-parameter
      double  par       ( const unsigned short l ,
                          const unsigned short m ) const ;
      /// get (l,m)-parameter
      double  parameter ( const unsigned short l ,
                          const unsigned short m ) const { return par (  l , m  ) ; }
      /// get k-parameter
      double  par       ( const unsigned int k ) const
      { return k < m_pars.size() ? m_pars[k] : 0.0 ; }
      /// get k-parameter
      double  parameter ( const unsigned int k ) const { return par ( k ) ; }
      /// get all parameters at once
      const std::vector<double>& pars() const { return m_pars ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the actual number of parameters
      std::size_t npars () const { return m_pars.size() ; }
      /// get lower edge
      double xmin () const { return m_xmin ; }
      /// get upper edge
      double xmax () const { return m_xmax ; }
      /// get lower edge
      double ymin () const { return m_ymin ; }
      /// get upper edge
      double ymax () const { return m_ymax ; }
      /// get the polynomial order (X)
      unsigned short nX () const { return m_nx ; }
      /// get the polynomial order (Y)
      unsigned short nY () const { return m_ny ; }
      // ======================================================================
    public:  // transformations
      // ======================================================================
      double x  ( const double tx ) const
      { return xmin ()  + ( xmax () - xmin () ) * tx ; }
      double y  ( const double ty ) const
      { return ymin ()  + ( ymax () - ymin () ) * ty ; }
      double tx ( const double x ) const
      { return  ( x - xmin () ) / ( xmax () - xmin () )      ; }
      double ty ( const double y ) const
      { return  ( y - ymin () ) / ( ymax () - ymin () )      ; }
      // ======================================================================
    public: // general integration
      // ======================================================================
      /** get the integral over 2D-region
       *  \f[ \int_{x_{low}}^{x_{high}}\int_{y_{low}}^{y_{high}}
       *      \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y \f]
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integral   ( const double xlow , const double xhigh ,
                          const double ylow , const double yhigh ) const ;
      /** integral over x-dimension
       *  \f[ \int_{x_{low}}^{x_{high}} \mathcal{B}(x,y) \mathrm{d}x\f]
       *  @param x     variable
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integrateX ( const double y    ,
                          const double xlow , const double xhigh ) const ;
      /** integral over y-dimension
       *  \f[ \int_{y_{low}}^{y_{high}} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param y     variable
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       */
      double integrateY ( const double x    ,
                          const double ylow , const double yhigh ) const ;
      // ======================================================================
    public: // special cases
      // ======================================================================
      /** get the integral over 2D-region
       *  \f[  x_min < x < x_max, y_min< y< y_max\f]
       */
      double integral   () const ;
      /** integral over x-dimension
       *  \f[ \int_{x_{min}}^{x_{max}} \mathcal{B}(x,y) \mathrm{d}x\f]
       *  @param x     variable
       */
      double integrateX ( const double y    ) const ;
      /** integral over y-dimension
       *  \f[ \int_{y_{min}}^{y_{max}} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param y     variable
       */
      double integrateY ( const double x    ) const ;
      // ======================================================================
    public: // few helper functions to expose internals
      // ======================================================================
      /// evaluate the basic polynomials
      double basicX ( const unsigned short i , const double         x ) const
      { return ( i > m_nx || x < m_xmin || x < m_xmax ) ? 0.0 : m_bx[i](x) ; }
      /// evaluate the basic polynomials
      double basicY ( const unsigned short i , const double         y ) const
      { return ( i > m_ny || y < m_ymin || y < m_ymax ) ? 0.0 : m_by[i](y) ; }
      /// expose some internals
      const Bernstein& basicX ( const unsigned short i ) const { return m_bx[i] ; }
      /// expose some internals
      const Bernstein& basicY ( const unsigned short i ) const { return m_by[i] ; }
      // ======================================================================
    private:
      // ======================================================================
      // polynom order in x-dimension
      unsigned short m_nx ; // polynom order in x-dimension
      // polynom order in y-dimension
      unsigned short m_ny ; // polynom order in y-dimension
      /// the list of parameters
      std::vector<double>  m_pars ;                // the list of parameters
      /// the left edge of interval
      double m_xmin  ;                             // the left edge of interval
      /// the right edge of interval
      double m_xmax  ;                             // the right edge of interval
      /// the left edge of interval
      double m_ymin  ;                             // the left edge of interval
      /// the right edge of interval
      double m_ymax  ;                             // the right edge of interval
      // ======================================================================
    private:
      // ======================================================================
      ///  vectors of basic  Bernstein polynomials
      typedef std::vector<Bernstein>  VB ;
      ///  vector  of basic  Bernstein polynomials
      VB m_bx ; //  vector  of basic  Bernetin polynomials
      ///  vector  of basic  Bernstein polynomials
      VB m_by ; //  vector  of basic  Bernetin polynomials
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Positive2D
     *  The "positive" 2D-polynomial of order Nx*Ny
     *  Actually it is a sum of basic bernstein 2D-polynomials with
     *  non-negative coefficients
     */
    class Positive2D : public std::binary_function<double,double,double>
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      Positive2D ( const unsigned short       Nx    =  1 ,
                   const unsigned short       Ny    =  1 ,
                   const double               xmin  =  0 ,
                   const double               xmax  =  1 ,
                   const double               ymin  =  0 ,
                   const double               ymax  =  1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double operator () ( const double x , const double y ) const
      { return m_bernstein ( x , y ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get number of parameters
      std::size_t npars () const { return m_sphere.nPhi () ; }
      /// set k-parameter
      bool setPar       ( const unsigned int k , const double value ) ;
      /// set k-parameter
      bool setParameter ( const unsigned int k , const double value )
      { return setPar   ( k , value ) ; }
      /// get the parameter value
      double  par       ( const unsigned int k ) const ;
      /// get the parameter value
      double  parameter ( const unsigned int k ) const { return par ( k ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get lower/upper edges
      double         xmin () const { return m_bernstein.xmin () ; }
      double         xmax () const { return m_bernstein.xmax () ; }
      double         ymin () const { return m_bernstein.ymin () ; }
      double         ymax () const { return m_bernstein.ymax () ; }
      // polynom order
      unsigned short nX   () const { return m_bernstein.nX   () ; }
      unsigned short nY   () const { return m_bernstein.nY   () ; }
      // ======================================================================
    public:
      // ======================================================================
      // transform variables
      double tx ( const double  x ) const { return m_bernstein.tx (  x ) ; }
      double ty ( const double  y ) const { return m_bernstein.ty (  y ) ; }
      double  x ( const double tx ) const { return m_bernstein. x ( tx ) ; }
      double  y ( const double ty ) const { return m_bernstein. y ( ty ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /** get the integral over 2D-region
       *  \f[ \int_{x_low}^{x_high}\int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f]
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integral   ( const double xlow , const double xhigh ,
                          const double ylow , const double yhigh ) const ;
      /** integral over x-dimension
       *  \f[ \int_{x_{low}}^{x_{high}} \mathcal{B}(x,y) \mathrm{d}x\f]
       *  @param x     variable
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       */
      double integrateX ( const double y    ,
                          const double xlow , const double xhigh ) const
      { return m_bernstein.integrateX ( y , xlow , xhigh ) ; }
      /** integral over y-dimension
       *  \f[ \int_{y_{low}}^{y_{high}} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param x     variable
       *  @param ylow  low  edge in x
       *  @param yhigh high edge in x
       */
      double integrateY ( const double x    ,
                          const double ylow , const double yhigh ) const
      { return m_bernstein.integrateY ( x , ylow , yhigh ) ; }
      // ======================================================================
    public: // specific
      // ======================================================================
      /** get the integral over 2D-region
       *  \f[ \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}}
       *        \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f]
       */
      double integral   () const ;
      /** integral over x-dimension
       *  \f[ \int_{x_{min}}^{x_{max}} \mathcal{B}(x,y) \mathrm{d}x\f]
       *  @param x     variable
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       */
      double integrateX ( const double y    ) const
      { return m_bernstein.integrateX ( y ) ; }
      /** integral over y-dimension
       *  \f[ \int_{y_{min}}^{y_{max}} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param x     variable
       *  @param ylow  low  edge in x
       *  @param yhigh high edge in x
       */
      double integrateY ( const double x    ) const
      { return m_bernstein.integrateY ( x ) ; }
      // ======================================================================
    public: // ingeredients
      // =====================================================================
      // get the bernstein polinomial in 2D
      const  Ostap::Math::Bernstein2D& bernstein () const
      { return m_bernstein ; }
      /// get the parameter sphere
      const  Ostap::Math::NSphere&     sphere    () const
      { return m_sphere ; }
      // ======================================================================
    private:
      // ======================================================================
      /// update bernstein coefficients
      bool updateBernstein () ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual bernstein polynomial
      Ostap::Math::Bernstein2D m_bernstein ; // the actual bernstein polynomial
      /// the external parameter sphere
      Ostap::Math::NSphere     m_sphere    ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Bernstein2DSym
     *  The symmetric Bernstein's polynomial of order N*N
     */
    class Bernstein2DSym : public std::binary_function<double,double,double>
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      Bernstein2DSym ( const unsigned short       n     =  1 ,
                       const double               xmin  =  0 ,
                       const double               xmax  =  1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double operator () ( const double x ,
                           const double y ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get number of parameters
      std::size_t npars () const { return m_pars.size() ; }
      /// set k-parameter
      bool setPar       ( const unsigned int   k     ,
                          const double         value ) ;
      /// set k-parameter
      bool setParameter ( const unsigned int   k     ,
                          const double         value )
      { return ( k < m_pars.size() ) && setPar ( k , value ) ; }
      /// set (l,m)-parameter
      bool setPar       ( const unsigned short l     ,
                          const unsigned short m     ,
                          const double         value ) ;
      /// set (l,m)-parameter
      bool setParameter ( const unsigned short l     ,
                          const unsigned short m     ,
                          const double         value )
      { return setPar   ( l , m  , value ) ; }
      /// get (l,m)-parameter
      double  par       ( const unsigned short l ,
                          const unsigned short m ) const ;
      /// get (l,m)-parameter value
      double  parameter ( const unsigned short l ,
                          const unsigned short m ) const { return par (  l , m  ) ; }
      /// get k-parameter
      double  par       ( const unsigned int   k ) const
      { return k < m_pars.size() ? m_pars [k] : 0.0 ; }
      /// get k-parameter
      double  parameter ( const unsigned int   k ) const { return par ( k ) ; }
      /// get all parameters at once
      const std::vector<double>& pars() const { return m_pars ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get lower edge
      double xmin () const { return m_xmin    ; }
      /// get upper edge
      double xmax () const { return m_xmax    ; }
      /// get lower edge
      double ymin () const { return   xmin () ; }
      /// get upper edge
      double ymax () const { return   xmax () ; }
      // ======================================================================
      unsigned short n  () const { return m_n  ; }
      unsigned short nX () const { return n () ; }
      unsigned short nY () const { return n () ; }
      // ======================================================================
    public:
      // ======================================================================
      double x  ( const double tx ) const
      { return xmin ()  + ( xmax () - xmin () ) * tx ; }
      double y  ( const double ty ) const
      { return ymin ()  + ( ymax () - ymin () ) * ty ; }
      double tx ( const double x ) const
      { return  ( x - xmin () ) / ( xmax () - xmin () ) ; }
      double ty ( const double y ) const
      { return  ( y - ymin () ) / ( ymax () - ymin () ) ; }
      // ======================================================================
    public: // generic integrals
      // ======================================================================
      /** get the integral over 2D-region
       *  \f[ \int_{x_{low}}^{x_{high}}\int_{y_{low}}^{y_{high}}
       *  \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f]
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integral   ( const double xlow , const double xhigh ,
                          const double ylow , const double yhigh ) const ;
      /** integral over x-dimension
       *  \f[ \int_{x_{low}}^{x_{high}} \mathcal{B}(x,y) \mathrm{d}x\f]
       *  @param x     variable
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       */
      double integrateX ( const double y    ,
                          const double xlow , const double xhigh ) const ;
      /** integral over y-dimension
       *  \f[ \int_{y_{low}}^{x_{high}} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param y     variable
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       */
      double integrateY ( const double x    ,
                          const double ylow , const double yhigh ) const ;
      // ======================================================================
    public: // specific integrals
      // ======================================================================
      /** get the integral over 2D-region
       *  \f[ \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}}
       *  \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f]
       */
      double integral   () const ;
      /** integral over x-dimension
       *  \f[ \int_{x_{min}}^{x_{max}} \mathcal{B}(x,y) \mathrm{d}x\f]
       *  @param x     variable
       */
      double integrateX ( const double y ) const ;
      /** integral over y-dimension
       *  \f[ \int_{y_{min}}^{x_{max}} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param y     variable
       */
      double integrateY ( const double x ) const ;
      // ======================================================================
    public: // few helper functions to expose internals
      // ======================================================================
      /// evaluate the basic polynomials
      double basic  ( const unsigned short i , const double         x ) const
      { return ( i > m_n || x < m_xmin || x < m_xmax ) ? 0.0 : m_b[i](x) ; }
      /// expose some internals
      const Bernstein& basic ( const unsigned short i ) const { return m_b[i] ; }
      // ======================================================================
    private:
      // ======================================================================
      // polynom order
      unsigned short m_n  ; // polynom order in x-dimension
      /// the list of parameters
      std::vector<double>  m_pars ;                // the list of parameters
      /// the left edge of interval
      double m_xmin  ;                             // the left edge of interval
      /// the right edge of interval
      double m_xmax  ;                             // the right edge of interval
      // ======================================================================
    private:
      // ======================================================================
      ///  vectors of basic  Bernstein polynomials
      typedef std::vector<Bernstein>  VB ;
      ///  vector  of basic  Bernetin polynomials
      VB m_b  ; //  vector  of basic  Bernstein polynomials
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Positive2DSym
     *  The "positive" symmetrical polynomial of order Nx*Ny
     *  Actually it is a sum of basic bernstein 2D-polynomials with
     *  non-negative coefficients
     */
    class Positive2DSym : public std::binary_function<double,double,double>
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      Positive2DSym ( const unsigned short       Nx    =  1 ,
                      const double               xmin  =  0 ,
                      const double               xmax  =  1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double operator () ( const double x , const double y ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get number of parameters
      std::size_t npars () const { return m_sphere.nPhi () ; }
      /// set k-parameter
      bool setPar       ( const unsigned int k , const double value ) ;
      /// set k-parameter
      bool setParameter ( const unsigned int k , const double value )
      { return setPar   ( k , value ) ; }
      /// get the parameter value
      double  par       ( const unsigned int k ) const ;
      /// get the parameter value
      double  parameter ( const unsigned int k ) const { return par ( k ) ; }
      /// get lower/upper edges
      double         xmin () const { return m_bernstein.xmin () ; }
      double         xmax () const { return m_bernstein.xmax () ; }
      double         ymin () const { return m_bernstein.ymin () ; }
      double         ymax () const { return m_bernstein.ymax () ; }
      // polynom order
      unsigned short n    () const { return m_bernstein.n    () ; }
      unsigned short nX   () const { return m_bernstein.nX   () ; }
      unsigned short nY   () const { return m_bernstein.nY   () ; }
      // ======================================================================
    public:
      // ======================================================================
      double tx ( const double  x ) const { return m_bernstein.tx (  x ) ; }
      double ty ( const double  y ) const { return m_bernstein.ty (  y ) ; }
      double  x ( const double tx ) const { return m_bernstein. x ( tx ) ; }
      double  y ( const double ty ) const { return m_bernstein. y ( ty ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /** get the integral over 2D-region
       *  \f[ \int_{x_low}^{x_high}\int_{y_low}^{y_high}
       *   \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y \f]
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integral   ( const double xlow , const double xhigh ,
                          const double ylow , const double yhigh ) const ;
      // ======================================================================
      /** integral over x-dimension
       *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param x     variable
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integrateX ( const double y    ,
                          const double xlow , const double xhigh ) const ;
      // ======================================================================
      /** integral over x-dimension
       *  \f[ \int_{x_low}^{x_high} \mathcal{B}(x,y) \mathrm{d}x\f]
       *  @param y     variable
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       */
      double integrateY ( const double x    ,
                          const double ylow , const double yhigh ) const ;
      // ======================================================================
    public: // specific
      // ======================================================================
      /** get the integral over 2D-region
       *  \f[ \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}}
       *   \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y \f]
       */
      double integral   () const ;
      /** integral over x-dimension
       *  \f[ \int_{x_{min}}^{x_{max}} \mathcal{B}(x,y) \mathrm{d}x\f]
       *  @param x     variable
       */
      double integrateX ( const double y ) const ;
      /** integral over y-dimension
       *  \f[ \int_{y_{min}}^{y_{max}} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param y     variable
       */
      double integrateY ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the bernstein 2D polynom
      const Ostap::Math::Bernstein2DSym& bernstein() const
      { return m_bernstein ; }
      /// get the parameter sphere
      const  Ostap::Math::NSphere&       sphere   () const
      { return m_sphere ; }
      // ======================================================================
    private:
      // ======================================================================
      /// update bernstein coefficients
      bool updateBernstein () ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual bernstein polynomial
      Ostap::Math::Bernstein2DSym m_bernstein ; // the actual bernstein polynomial
      /// Parameter sphere
      Ostap::Math::NSphere        m_sphere ;
      // ======================================================================
    } ;
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_BERNSTEIN_H
// ============================================================================
