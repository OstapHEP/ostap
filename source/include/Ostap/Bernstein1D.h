// ============================================================================
#ifndef OSTAP_BERNSTEIN1D_H
#define OSTAP_BERNSTEIN1D_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <cmath>
#include <array>
#include <iterator>
#include <vector>
#include <complex>
#include <algorithm>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/NSphere.h"
#include "Ostap/Bernstein.h"
// ============================================================================
/** @file Ostap/Bernstein1D.h
 *  Set of useful math-functions, related to Bernstein polynomials
 *  @see http://en.wikipedia.org/wiki/Bernstein_polynomial
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2010-04-19
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    // Forward declarations of important classes 
    // ========================================================================
    class Positive   ;
    class Monotonic  ;
    class Convex     ;
    class ConvexOnly ;
    // ========================================================================
    /** @class BernsteinEven
     *  A special case of Bernstein polynomial with symmetry relatiev to midpoints
     *  \f$ f( \frac{x_{max}+x_{min}}{2} - x ) \equiv  \frac{x_{max}+x_{min}}{2} + x ) \f$
     *  @see Ostap::Math::Bernstein
     *  @author Vanya Belyaev Ivan.Belyaev@iep.ru
     *  @date 2016-10-02
     */
    class BernsteinEven
    {
    public:
      // ======================================================================
      /** constructor
       *  @param N    degree of even Bernstein polynomial 
       *  @param xmin low edge
       *  @param xmax high edge
       */
      BernsteinEven
      ( const unsigned short N     = 0 ,
        const double         xmin  = 0 ,
        const double         xmax  = 1 ) ;
      // ======================================================================
      /** constructor from list of parameters 
       *  @param pars vector of parameters 
       *  @param xmin low edge
       *  @param xmax high edge
       */
      BernsteinEven 
      ( const std::vector<double>& pars      ,
        const double               xmin  = 0 ,
        const double               xmax  = 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value of the polynomial
      double evaluate ( const double x ) const { return m_bernstein.evaluate ( x ) ; }
      // ======================================================================
      /// get the value
      double operator () ( const double x ) const { return m_bernstein ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// the degree of polynomial
      unsigned short degree () const { return m_bernstein.degree ()     ; }
      /// number of parameters
      unsigned short npars  () const { return m_bernstein.npars  () / 2 ; }
      /// all zero ?
      bool           zero   () const { return m_bernstein.zero()        ; }
      // ======================================================================
      /** set k-parameter
       *  @param k index
       *  @param value new value
       *  @return true if parameter is actually changed
       */
      bool setPar  
      ( const unsigned short k     , 
        const double         value ) ;
      // ======================================================================
      /** set k-parameter
       *  @param k index
       *  @param value new value
       *  @return true if parameter is actually changed
       */
      bool setParameter  
      ( const unsigned short k     , 
        const double         value ) 
      { return setPar      ( k , value ) ; }
      // =====================================================================
      /** set several/all parameters at once 
       *  @param begin  start itertaor for the sequence of coefficients 
       *  @param end    end   itertaor for the sequence of coefficients 
       *  @return true if at least one parameter is actually changed 
       */
      template<typename ITERATOR,
               typename value_type = typename std::iterator_traits<ITERATOR>::value_type ,
               typename = std::enable_if<std::is_convertible<value_type,long double>::value> >
      inline bool setPars 
      ( ITERATOR begin , 
        ITERATOR end   ) 
      {
        bool updated = false ;
        const unsigned short npb = m_bernstein.npars  () ;
        for ( unsigned short k   = 0 ; k < npb && begin != end ; ++k, ++begin )
        {
          const bool updated1 = m_bernstein.setPar (       k    , *begin ) ;
          const bool updated2 = m_bernstein.setPar ( npb - k -1 , *begin ) ;
          updated = updated1 || updated2 ? true : updated ;
        }
        return updated ;
      }
      /// get the parameter value
      double  par          ( const unsigned short k ) const
      { return m_bernstein.par ( k ) ; }
      /// get the parameter value
      double  parameter    ( const unsigned short k ) const { return par ( k ) ; }
      /// get all parameters (by value!!! COPY!!)
      std::vector<double> pars () const ;
      // ======================================================================
    public: // convert from local to global variables
      // ======================================================================
      /// local to global
      double x ( const double t ) const { return m_bernstein.x ( t ) ; }
      /// gloal to local 
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
      BernsteinEven __add__     ( const double value ) const ;
      /// Sum of Bernstein polynomial and a constant
      BernsteinEven __radd__    ( const double value ) const ;
      /// Product of Bernstein polynomial and a constant
      BernsteinEven __mul__     ( const double value ) const ;
      /// Product of Bernstein polynomial and a constant
      BernsteinEven __rmul__    ( const double value ) const ;
      /// Subtract a constant from Benrstein polynomial
      BernsteinEven __sub__     ( const double value ) const ;
      /// Constant minus Bernstein polynomial
      BernsteinEven __rsub__    ( const double value ) const ;
      /// Divide Benrstein polynomial by a constant
      BernsteinEven __truediv__ ( const double value ) const ;
      BernsteinEven __div__     ( const double value ) const { return __truediv__ ( value ) ; }
      // ======================================================================
      Bernstein __add__  ( const Bernstein& a ) const { return m_bernstein + a ; }
      Bernstein __radd__ ( const Bernstein& a ) const { return a + m_bernstein ; }
      Bernstein __sub__  ( const Bernstein& a ) const { return m_bernstein - a ; }
      Bernstein __rsub__ ( const Bernstein& a ) const { return a - m_bernstein ; }
      Bernstein __mul__  ( const Bernstein& a ) const { return m_bernstein * a ; }
      Bernstein __rmul__ ( const Bernstein& a ) const { return a * m_bernstein ; }
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
      /// get the unique tag 
      std::size_t tag () const { return m_bernstein.tag () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// swap two objects 
      void swap ( BernsteinEven& right ) 
      { Ostap::Math::swap ( m_bernstein , right.m_bernstein ) ; } 
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
      /// the actual underlying "regular" Bernstein polynomial
      Bernstein      m_bernstein ; // the actual Bernstein polynomial
      // ======================================================================
    } ;
    // ========================================================================
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
    inline Bernstein operator+ ( const BernsteinEven& a , const Bernstein&     b ) 
    { return a.bernstein () + b ; }
    inline Bernstein operator+ ( const Bernstein&     a , const BernsteinEven& b ) 
    { return a + b.bernstein () ; }
    inline Bernstein operator- ( const BernsteinEven& a , const Bernstein&     b ) 
    { return a.bernstein () - b ; }
    inline Bernstein operator- ( const Bernstein&     a , const BernsteinEven& b ) 
    { return a - b.bernstein () ; }
    inline Bernstein operator* ( const BernsteinEven& a , const Bernstein&     b ) 
    { return a.bernstein () * b ; }
    inline Bernstein operator* ( const Bernstein&     a , const BernsteinEven& b  ) 
    { return a * b.bernstein () ; }
    // ========================================================================
    /// Swapping function for even bernstein polynomials 
    inline void swap 
    ( BernsteinEven& a , 
      BernsteinEven& b ) { a.swap ( b ) ; }
    // ========================================================================
    /** @class Positive
     *  The "positive" polynomial of order N.
     *
     *  Positive polynomials are described according to  Karlin and Shapley 
     *  @see S.Karlin and L.S. Shapley,"Geometry of Moment Space",
     *       Memoirs of the Amer.Math.Soc., 12, 1953 
     *  @see https://bookstore.ams.org/memo-1-12/
     *
     *  For \f$n=2m\f$ the polynomial, non-negative at the \f$ 0\le x \le 1\f$ 
     *  interval is written as 
     *
     *  \f[ P_{2m}(x) = 
     *  \alpha        A\sum_{j=1}^{m}   \left(x-x_{2j-1} \right)^2} + 
     *  \beta x (1-x) B\sum_{j=1}^{m-1} \left(x-x_{2j}   \right)^2} \f] 
     *
     *  and for \f$n=2m+1\f$ the polynomial is :
     *
     *  \f[ P_{2m+1}(x) = 
     *  \alpha (1-x) A\sum_{j=1}^{m} \left(x-x_{2j-1} \right)^2} + 
     *  \beta     x  B\sum_{j=1}^{m} \left(x-x_{2j}   \right)^2} \f] 
     *  where \f$0\le x_{1} \le ... \le x_{n-1} \le 1\f$ and 
     *  \f$ 0 < \alpha \f$, \f$ 0 < \beta \f$
     *
     *  - normalization constants \f$ A\f$ and \f$ B \f$ are chosen such 
     *    that corresponding polynomial terms have unit integrals 
     *  - constants \f$ \alpha \f$ and \f$ \\beta \f$ are parameterized
     *    using the phase \f$ \phi_0\f$, such that 
     *    \f$ \alpha = \cos^2 \phi_0 \f$ and \f$ \beta  = \sin^2 \phi_0 \f$
     *  - The ordered parameters $0 \le ... \le x_i \le x_{i+1} \le ... \le 1 $ 
     *    are parameterized in terms of multidimensional sphere, 
     *    using \f$ n-1\f$ phases \f$ \phi_i \f$, \f$ 1 \le i < n\f$ 
     * 
     *  The special cases:  
     *  - For \f$ n = 0 \f$ polynom is \f$ P_0(x) \equiv 1  \f$
     *  - For \f$ n = 1 \f$ polynom is \f$ P_1(x) \equiv \cos^2 \phi_0 ( 1 - x ) + \sin^2 \phi_0 x \f$
     *
     *  With such choice of parameters \f$ \phi_i \f$, \f$ 0 \le i < n \f$ one gets 
     *  - positivity      \f$ 0 \le P_{n}(x)         \f$ for \f$ 0 \le x \le 1 \f$ 
     *  - normalisation:  \f$ \int_0^1 P_n(x) dx = 1 \f$ 
     *
     *  @see Ostap::Math::Bernstein
     *  @author Vanya BELYAEV Ivan.Belayev@itep.ru
     *  @date 2010-04-19
     */
    class Positive 
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from the order
       *  @param N    degree of polynomial
       *  @param xmin low-edge 
       *  @param xmax high-edge 
       */
      Positive
      ( const unsigned short        N     =  1 ,
        const double                xmin  =  0 ,
        const double                xmax  =  1 ) ;
      // ======================================================================
      /** constructor from N parameters/phases 
       *  @param phases  list of parameters/phase 
       *  @param xmin low-edge 
       *  @param xmax high-edge 
       */
      Positive
      ( const std::vector<double>&  phases     ,
        const double                xmin  =  0 ,
        const double                xmax  =  1 ) ;
      // ======================================================================
      /** constructor from the sequence of   parameters 
       *  @param begin  start-iterator for sequence of coefficients 
       *  @param end     start-iterator for sequence of coefficients 
       *  @param xmin low-edge 
       *  @param xmax high-edge 
       */
      template<typename ITERATOR,
               typename value_type = typename std::iterator_traits<ITERATOR>::value_type ,
               typename = std::enable_if<std::is_convertible<value_type,long double>::value> >
      Positive
      ( ITERATOR     begin    , 
        ITERATOR     end      , 
        const double xmin = 0 ,
        const double xmax = 1 ) 
        : Positive ( std::distance ( begin , end ) , xmin , xmax ) 
      { setPars ( begin , end ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double operator () ( const double x ) const { return m_bernstein ( x ) ; }
      // ======================================================================
    public: // PAR-interface 
      // ======================================================================
      /// get number of parameters (==degree)
      unsigned short npars () const { return m_bernstein.degree () ; }
      // =======================================================================
      /// get the parameter value
      double  par       ( const unsigned short k ) const
      { 
        const unsigned short nA = m_sphereA.npars () ;
        const unsigned short nR = m_sphereR.npars () ;
        return 
          k < nA      ? m_sphereA.par ( k      ) :
          k < nA + nR ? m_sphereR.par ( k - nA ) : 0.0 ;
      }
      // ======================================================================
      /// get the paramete  values 
      double parameter ( const unsigned short k ) const { return par ( k ) ; }
      /// get all parameters (phases on sphere) BY VALUE!
      std::vector<double> pars  () const ;
      // ======================================================================
      /// set k-parameter
      bool setPar     
      ( const unsigned short k     , 
        const double         value ) 
      {
        const unsigned short nA = m_sphereA.npars () ;
        const unsigned short nR = m_sphereR.npars () ;
        const bool update = 
          k < nA      ? m_sphereA.setPhase ( k      , value ) :
          k < nA + nR ? m_sphereR.setPhase ( k - nA , value ) : false ;
        return update ? updateBernstein () : false ;
      }
      // ======================================================================
      /// set k-parameter
      bool setParameter ( const unsigned short k , const double value )
      { return setPar   ( k , value ) ; }
      // ======================================================================
      /** set several/all parameters at once 
       *  @param begin  start itertaor for the sequence of coefficients 
       *  @param end    end   itertaor for the sequence of coefficients 
       *  @return true if at least one parameter is actually changed 
       */
      template<typename ITERATOR,
               typename value_type = typename std::iterator_traits<ITERATOR>::value_type ,
               typename = std::enable_if<std::is_convertible<value_type,long double>::value> >
      inline bool setPars 
      ( ITERATOR begin , 
        ITERATOR end   ) 
      {
        const auto NN = std::distance ( begin  , end ) ;
        const unsigned short nA = m_sphereA.npars () ;
        const bool updatedA =           m_sphereA.setPars ( begin      , end ) ;
        const bool updatedR = nA < NN ? m_sphereR.setPars ( begin + nA , end ) : false ;
        return updatedA || updatedR ? updateBernstein ( ) : false ;
      }
      // ======================================================================
      /** set several/all parameters at once 
       *  @param pars (INPUT) vector of parameters 
       *  @return true if at least one parameter is actually changed 
       */
      inline bool setPars ( const std::vector<double>& pars ) 
      { return setPars ( pars.begin () , pars.end () ) ; }
      // ======================================================================
    public:
      // ======================================================================
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
      /// get lower edge
      double xmin () const { return m_bernstein.xmin () ; }
      /// get upper edge
      double xmax () const { return m_bernstein.xmax () ; }
      /// transform variables
      double x ( const double t ) const { return m_bernstein. x ( t )  ; }
      double t ( const double x ) const { return m_bernstein. t ( x )  ; }
      // ======================================================================
    public:
      // ======================================================================
      /// constant ?
      bool constant    () const { return m_bernstein.constant   () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between xmin and xmax
      double integral   () const { return 1 ; }
      /// get the integral between low and high
      double integral   ( const double low , const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the underlying Bernstein polynomial
      const Ostap::Math::Bernstein& bernstein () const { return m_bernstein ; }
      /// get the parameter sphere
      const Ostap::Math::NSphere&   asphere   () const { return m_sphereA   ; }
      /// get the parameter sphere
      const Ostap::Math::NSphere&   rsphere   () const { return m_sphereR   ; }
      /// get the indefinite integral
      Bernstein indefinite_integral ( const double C = 0 ) const
      { return m_bernstein.indefinite_integral ( C ) ; }
      /// get the derivative
      Bernstein derivative () const { return m_bernstein.derivative () ; }
      /// get the derivative at point x 
      double derivative             ( const double x ) const
      { return m_bernstein.derivative ( x  ) ; }
      // ======================================================================
    public:  /// basic operations  for python
      // ======================================================================
      /// Sum of Bernstein polynomial and a constant
      Bernstein __add__     ( const double value ) const { return m_bernstein + value   ; }
      /// Sum of Bernstein polynomial and a constant
      Bernstein __radd__    ( const double value ) const { return m_bernstein + value   ; }
      /// Product of Bernstein polynomial and a constant
      Bernstein __mul__     ( const double value ) const { return m_bernstein * value   ; }
      /// Product of Bernstein polynomial and a constant
      Bernstein __rmul__    ( const double value ) const { return m_bernstein * value   ; }
      /// Subtract a constant from Benrstein polynomial
      Bernstein __sub__     ( const double value ) const { return m_bernstein - value   ; }
      /// Constant minus Bernstein polynomial
      Bernstein __rsub__    ( const double value ) const { return value - m_bernstein   ; }
      /// Divide Bernstein polynomial by a constant
      Bernstein __truediv__ ( const double value ) const { return m_bernstein / value   ; }
      Bernstein __div__     ( const double value ) const { return __truediv__ ( value ) ; }
      /// Negate Bernstein polynomial
      Bernstein __neg__     () const { return -m_bernstein ; }
      // ======================================================================
      Bernstein __add__  ( const Bernstein& a ) const { return m_bernstein + a ; }
      Bernstein __radd__ ( const Bernstein& a ) const { return a + m_bernstein ; }
      Bernstein __sub__  ( const Bernstein& a ) const { return m_bernstein - a ; }
      Bernstein __rsub__ ( const Bernstein& a ) const { return a - m_bernstein ; }
      Bernstein __mul__  ( const Bernstein& a ) const { return m_bernstein * a ; }
      Bernstein __rmul__ ( const Bernstein& a ) const { return a * m_bernstein ; }
       // ======================================================================
    public:
      // ======================================================================
      /// get the tag
      std::size_t tag () const { return m_bernstein.tag () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// swap two objects 
      void swap ( Positive& right ) 
      {
        Ostap::Math::swap ( m_bernstein , right.m_bernstein ) ;
        Ostap::Math::swap ( m_sphereA   , right.m_sphereA   ) ;
        Ostap::Math::swap ( m_sphereR   , right.m_sphereR   ) ;
        //
        std::swap ( m_rs  , right.m_rs  ) ;
        std::swap ( m_v1  , right.m_v1  ) ;
        std::swap ( m_v2  , right.m_v2  ) ;
        std::swap ( m_aux , right.m_aux ) ;
      } 
      // ======================================================================
    private : 
      // ======================================================================
      /// update bernstein coefficients
      bool updateBernstein () ;
      // ======================================================================
    protected:
      // ======================================================================
      /// the actual bernstein polynomial
      Ostap::Math::Bernstein m_bernstein ; // the actual bernstein polynomial
      /// (alpha/beta) sphere
      Ostap::Math::NSphere   m_sphereA   ;
      /// sphere of roots 
      Ostap::Math::NSphere   m_sphereR   ;
      // ======================================================================
    private:
      // ======================================================================
      std::vector<long double> m_rs  ;
      std::vector<long double> m_v1  ;
      std::vector<long double> m_v2  ;
      std::vector<long double> m_aux ;
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
    inline Bernstein operator+ ( const Positive&  a , const Bernstein& b ) 
    { return a.bernstein () + b ; }
    inline Bernstein operator- ( const Positive&  a , const Bernstein& b ) 
    { return a.bernstein () - b ; }
    inline Bernstein operator+ ( const Bernstein& b , const Positive&  a ) { return a + b ; }
    inline Bernstein operator- ( const Bernstein& b , const Positive&  a ) 
    { return b - a.bernstein () ; }
    inline Bernstein operator* ( const Positive&  a , const Bernstein& b ) 
    { return a.bernstein () * b ; }
    inline Bernstein operator* ( const Bernstein& b , const Positive&  a ) { return a * b ; }
    // ========================================================================
    /// Swapping function for positive bernstein polynomials 
    inline void swap ( Positive&      a , Positive&      b ) { a.swap ( b ) ; }
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
    class PositiveEven
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      PositiveEven
      ( const unsigned short        N     =  1 ,
        const double                xmin  =  0 ,
        const double                xmax  =  1 ) ;
      // ======================================================================
      /// constructor from N phases
      PositiveEven 
      ( const std::vector<double>&  phases     ,
        const double                xmin  =  0 ,
        const double                xmax  =  1 ) ;
      // =====================================================================
      /** constructor from the sequence of   parameters 
       *  @param begin  start-iterator for sequence of coefficients 
       *  @param end     start-iterator for sequence of coefficients 
       *  @param xmin low-edge 
       *  @param xmax high-edge 
       */
      // template<typename ITERATOR,
      //          typename value_type = typename std::iterator_traits<ITERATOR>::value_type ,
      //          typename = std::enable_if<std::is_convertible<value_type,long double>::value> >
      // PositiveEven 
      // ( ITERATOR     begin    , 
      //   ITERATOR     end      , 
      //   const double xmin = 0 ,
      //   const double xmax = 1 ) 

      //   ...


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
      std::size_t npars () const { return m_positive.npars () ; }
      /// get the parameter value
      double  par       ( const unsigned short k ) const
      { return m_positive.par ( k ) ; }
      double  parameter ( const unsigned short k ) const { return par ( k )  ; }
      /// set k-parameter
      bool setPar       ( const unsigned short k , const double value )        
      { return m_positive.setPar ( k , value ) ? updateBernstein() : false ; }
      /// set k-parameter
      bool setParameter ( const unsigned short k , const double value )
      { return setPar   ( k , value ) ; }
      /// get all parameters (copy by value) 
      std::vector<double> pars  () const { return m_positive.pars () ; }
      // ======================================================================
      template<typename ITERATOR,
               typename value_type = typename std::iterator_traits<ITERATOR>::value_type ,
               typename = std::enable_if<std::is_convertible<value_type,long double>::value> >
      inline bool setPars 
      ( ITERATOR begin , 
        ITERATOR end   ) 
      { return m_positive.setPars ( begin , end ) ? updateBernstein() : false ; }
      // ======================================================================
    public:  // some characteristics
      // ======================================================================
      /// degree
      unsigned short degree      () const { return m_even.degree() ; }
      // ======================================================================
    public:
      // ======================================================================
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
      double integral   () const { return 1 ; }
      /// get the integral between low and high
      double integral   ( const double low , const double high ) const 
      { return m_even.integral ( low , high ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the underlying Bernstein Even polynomial
      const Ostap::Math::BernsteinEven& bernsteinEven () const { return m_even ; }
      /// get the underlying Bernstein Even polynomial
      const Ostap::Math::BernsteinEven& even          () const { return m_even ; }
      /// get the underlying Bernstein polynomial
      const Ostap::Math::Bernstein&     bernstein () const { return m_even.bernstein()  ; }
      /// get the indefinite integral
      Bernstein indefinite_integral ( const double C = 0 ) const
      { return m_even.indefinite_integral ( C )  ; }
      /// get the derivative
      Bernstein derivative          () const
      { return m_even.derivative          ()     ; }
      /// get the derivative at point x 
      double derivative             ( const double x ) const
      { return m_even.derivative          ( x  ) ; }
      // ======================================================================
    public:  /// basic operations  for python
      // ======================================================================
      /// Sum of Bernstein polynomial and a constant
      BernsteinEven __add__     ( const double value ) const { return m_even + value        ; }
      /// Sum of Bernstein polynomial and a constant
      BernsteinEven __radd__    ( const double value ) const { return m_even + value        ; }
      /// Product of Bernstein polynomial and a constant
      BernsteinEven __mul__     ( const double value ) const { return m_even * value        ; }
      /// Product of Bernstein polynomial and a constant
      BernsteinEven __rmul__    ( const double value ) const { return m_even * value        ; }
      /// Subtract a constant from Benrstein polynomial
      BernsteinEven __sub__     ( const double value ) const { return m_even - value        ; }
      /// Constant minus Bernstein polynomial
      BernsteinEven __rsub__    ( const double value ) const { return value  - m_even       ; }
      /// Divide Bernstein polynomial by a constant
      BernsteinEven __truediv__ ( const double value ) const { return m_even / value        ; }
      BernsteinEven __div__     ( const double value ) const { return __truediv__ ( value ) ; }
      // ======================================================================
      Bernstein __add__  ( const Bernstein& a ) const { return m_even + a ; }
      Bernstein __radd__ ( const Bernstein& a ) const { return a + m_even ; }
      Bernstein __sub__  ( const Bernstein& a ) const { return m_even - a ; }
      Bernstein __rsub__ ( const Bernstein& a ) const { return a - m_even ; }
      Bernstein __mul__  ( const Bernstein& a ) const { return m_even * a ; }
      Bernstein __rmul__ ( const Bernstein& a ) const { return a * m_even ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const { return m_even.tag () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// swap two objects 
      void swap ( PositiveEven& right ) 
      {
        Ostap::Math::swap ( m_positive , right.m_positive ) ;
        Ostap::Math::swap ( m_even     , right.m_even     ) ;
      } 
      // ======================================================================
    private: 
      // ======================================================================
      /// update bernstein coefficients
      bool updateBernstein () ;
      // ======================================================================
    protected:
      // ======================================================================
      /// the actual bernstein even polynomial
      Ostap::Math::BernsteinEven m_even     ; // the actual bernstein polynomial
      /// helper polymnomial to set parameters 
      Ostap::Math::Positive      m_positive ; // helper polymnomial to set parameters 
      // ======================================================================
    } ;
    // ========================================================================
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
    inline Bernstein operator+ ( const PositiveEven& a , const Bernstein& b ) 
    { return a.bernstein () + b ; }
    inline Bernstein operator- ( const PositiveEven& a , const Bernstein& b ) 
    { return a.bernstein () - b ; }
    inline Bernstein operator+ ( const Bernstein& b , const PositiveEven&  a ) { return a + b ; }
    inline Bernstein operator- ( const Bernstein& b , const PositiveEven&  a ) 
    { return b - a.bernstein () ; }
    inline Bernstein operator* ( const PositiveEven&  a , const Bernstein& b ) 
    { return a.bernstein () * b ; }
    inline Bernstein operator* ( const Bernstein& b , const PositiveEven&  a ) { return a * b ; }
    // ========================================================================
    /// Swapping function for positive even bernstein polynomials 
    inline void swap ( PositiveEven&  a , PositiveEven&  b ) { a.swap ( b ) ; }
    // ========================================================================
    /** @class Monotonic
     *  The "positive" monotonic polynomial of order N
     * 
     *  Conceptually, the monotonic increasing polynomial is parameterised as 
     *  \f[ I_n (x) \propto \int_0^x P_{n-1}(y) dy + C \f]
     *  similarly monitonically decreasing polynomial is parameterized as 
     *  \f[ D_n (x) \propto \int_x^1 P_{n-1}(y) dy + C \f]
     *  where \f$P_n(x)\f$ is a positive polynomial and 
     *  \f$ C \f$ is an integration constant.  
     *  Taking the normalization constraints, 
     *  the polynomilas  are described as 
     *  \f[ I_n (x) = \cos^2 \phi_0 A \int_0^x P_{n-1}(y|\phi_i) dy                    + \sin^2\phi0 \f] 
     *   and 
     *  \f[ D_n (x) = \cos^2 \phi_0 B \left( 1 - \int_0^x P_{n-1}(y|\phi_i) dy \right) + \sin^2\phi0 \f]
     *  where paramers \f$ A \f$ and \f$ B \f$ are chosen to provide 
     * \f$ A \int_0^1 \left ( \int_0^x P_{n-1}(y) dy     \right) dx \equiv 1 \f$  and 
     * \f$ B \int_0^1 \left ( 1 - \int_0^x P_{n-1}(y) dy \right) dx \equiv 1 \f$, 
     *  (note that \f$ \int_0^1 P_n(y)dy \equiv 1 \f$).
     * 
     *  The (n-1) parameters \f$ \phi_i\f$, \f$ 1 \le i < n \f$, are used to 
     *  parameterise normalized positive polynomial \f$ P_{n-1}(x) \f$
     *  
     *  @see Ostap::Math::Positive 
     *  @see Ostap::Math::Bernstein
     *  @see Vanya BELYAEV Ivan.Belyaev@itep.ru
     */
    class Monotonic
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from the order
       *  @param N    degree of polynomial
       *  @param xmin low-edge 
       *  @param xmax high-edge 
       *  @param increasing  increasing  or decreasing ?
       */
      Monotonic
      ( const unsigned short       N          =    1 ,
        const double               xmin       =    0 ,
        const double               xmax       =    1 ,
        const bool                 increasing = true ) ;
      // ======================================================================
      /** constructor from N parameters/phases 
       *  @param phases  list of parameters/phase 
       *  @param xmin low-edge 
       *  @param xmax high-edge 
       *  @param increasing  increasing  or decreasing ?
       */
      Monotonic
      ( const std::vector<double>& pars              ,
        const double               xmin       =    0 ,
        const double               xmax       =    1 ,
        const bool                 increasing = true ) ;
      // ======================================================================
      /** constructor from the sequence of parameters 
       *  @param begin  start-iterator for sequence of coefficients 
       *  @param end     start-iterator for sequence of coefficients 
       *  @param xmin low-edge 
       *  @param xmax high-edge 
       */
      template<typename ITERATOR,
               typename value_type = typename std::iterator_traits<ITERATOR>::value_type ,
               typename = std::enable_if<std::is_convertible<value_type,long double>::value> >
      Monotonic 
      ( ITERATOR     begin             , 
        ITERATOR     end               , 
        const double xmin       = 0    ,
        const double xmax       = 1    , 
        const bool   increasing = true ) 
        : Monotonic ( std::distance (  begin , end ) , xmin , xmax , increasing ) 
      { setPars  ( begin , end ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double operator () ( const double x ) const { return m_bernstein ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get number of parameters
      std::size_t npars () const { return m_bernstein.degree() ; }
      /// set k-parameter
      bool setPar       ( const unsigned short k , const double value ) 
      { 
        const unsigned short nA = m_sphere  .npars () ;
        const unsigned short nP = m_positive.npars () ;
        const bool update = 
          k < nA      ? m_sphere  .setPar ( k      , value ) :
          k < nA + nP ? m_positive.setPar ( k - nA , value ) : false ;
        return update ? updateBernstein () : false ;
      }
      /// set k-parameter
      bool setParameter ( const unsigned short k , const double value )
      { return setPar   ( k , value ) ; }
      /// get the parameter value
      double  par       ( const unsigned short k ) const
      { 
        const unsigned short nA = m_sphere  .npars () ;
        const unsigned short nP = m_positive.npars () ;
        return 
          k < nA      ? m_sphere  .par ( k ) :
          k < nA + nP ? m_positive.par ( k - nA ) : 0.0 ;
      }
      // ======================================================================
      /// get all parameters (phases on sphere)
      std::vector<double>        pars  () const ;
      /// get bernstein coefficients
      const std::vector<double>& bpars () const { return m_bernstein.pars () ; }
      // ======================================================================
      /** set several/all parameters at once 
       *  @param begin  start itertaor for the sequence of coefficients 
       *  @param end    end   itertaor for the sequence of coefficients 
       *  @return true if at least one parameter is actually changed 
       */
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type ,
                typename = std::enable_if<std::is_convertible<value_type,long double>::value> >
      inline bool setPars 
      ( ITERATOR begin , 
        ITERATOR end   ) 
      {
        const auto NN = std::distance ( begin  , end ) ;
        const unsigned short nS = m_sphere.npars () ;
        const bool updatedS =           m_sphere  .setPars ( begin      , end ) ;
        const bool updatedP = nS < NN ? m_positive.setPars ( begin + nS , end ) : false ;
        return updatedS || updatedP ? updateBernstein ( ) : false ;
      }
      // ======================================================================
      /** set several/all parameters at once 
       *  @param pars (NIPUT) vector of parameters 
       *  @return true if at least one parameter is actually changed 
       */
      inline bool setPars ( const std::vector<double>& pars ) 
      { return setPars ( pars.begin() , pars.end() ) ; }
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
      // ======================================================================
      /// increasing ?
      bool increasing () const { return degree() < 1 ||  m_increasing ; }
      /// decreasing ?
      bool decreasing () const { return degree() < 1 || !m_increasing ; }
      /// monotonic
      bool monotonic  () const { return true  ; }
      //// constant 
      bool constant   () const { return m_bernstein.constant () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the minimal value of function
      double fun_min () const ; // get the minimal value of function
      /// get the maximal value of function
      double fun_max () const ; // get the maximal value of function
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between xmin and xmax
      double integral   () const { return 1 ; } 
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
    public:  /// basic operations  for python
      // ======================================================================
      /// Sum of Bernstein polynomial and a constant
      Bernstein __add__     ( const double value ) const { return m_bernstein + value   ; }
      /// Sum of Bernstein polynomial and a constant
      Bernstein __radd__    ( const double value ) const { return m_bernstein + value   ; }
      /// Product of Bernstein polynomial and a constant
      Bernstein __mul__     ( const double value ) const { return m_bernstein * value   ; }
      /// Product of Bernstein polynomial and a constant
      Bernstein __rmul__    ( const double value ) const { return m_bernstein * value   ; }
      /// Subtract a constant from Benrstein polynomial
      Bernstein __sub__     ( const double value ) const { return m_bernstein - value   ; }
      /// Constant minus Bernstein polynomial
      Bernstein __rsub__    ( const double value ) const { return value - m_bernstein   ; }
      /// Divide Bernstein polynomial by a constant
      Bernstein __truediv__ ( const double value ) const { return m_bernstein / value   ; }
      Bernstein __div__     ( const double value ) const { return __truediv__ ( value ) ; }
      /// Negate Bernstein polynomial
      Bernstein __neg__     () const { return -m_bernstein ; }
      // ======================================================================
      Bernstein __add__  ( const Bernstein& a ) const { return m_bernstein + a ; }
      Bernstein __radd__ ( const Bernstein& a ) const { return a + m_bernstein ; }
      Bernstein __sub__  ( const Bernstein& a ) const { return m_bernstein - a ; }
      Bernstein __rsub__ ( const Bernstein& a ) const { return a - m_bernstein ; }
      Bernstein __mul__  ( const Bernstein& a ) const { return m_bernstein * a ; }
      Bernstein __rmul__ ( const Bernstein& a ) const { return a * m_bernstein ; }
       // ======================================================================
    public:
      // ======================================================================
      /// get the tag
      std::size_t tag () const { return m_bernstein.tag () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// swap two objects 
      void swap ( Monotonic& right ) 
      {
        Ostap::Math::swap ( m_bernstein  , right.m_bernstein  ) ;
        Ostap::Math::swap ( m_positive   , right.m_positive   ) ;
        Ostap::Math::swap ( m_sphere     , right.m_sphere     ) ;
        std::swap         ( m_increasing , right.m_increasing ) ;
        std::swap         ( m_aux        , right.m_aux        ) ;
      } 
      // ======================================================================
    private : 
      // ======================================================================
      /// update bernstein coefficients
      bool updateBernstein () ;
      // ======================================================================
    protected:
      // ======================================================================
      /// the actual bernstein polynomial
      Ostap::Math::Bernstein   m_bernstein  ; // the actual bernstein polynomial
      /// the positive helper polynomial 
      Ostap::Math::Positive    m_positive   ; // the positive helper polynomial
      /// parameters sphere
      Ostap::Math::NSphere     m_sphere     ; // parameters sphere
      /// increasing ?
      bool                     m_increasing ; // increasing ?
      /// auxillary array for computations
      std::vector<long double> m_aux        ; // auxillary array for computations
      // ======================================================================
    } ;
    // ========================================================================
    ///  Monotonic plus      constant
    inline Bernstein operator+( const Monotonic& p , const double v )
    { return p.bernstein() + v ; }
    ///  Monotonic multiply  constant
    inline Bernstein operator*( const Monotonic& p , const double v )
    { return p.bernstein() * v ; }
    ///  Monotonic minus     constant
    inline Bernstein operator-( const Monotonic& p , const double v )
    { return p.bernstein() - v ; }
    ///  Monotonic divide constant
    inline Bernstein operator/( const Monotonic& p , const double v )
    { return p.bernstein() / v ; }
    ///  Constant plus  Monotonic
    inline Bernstein operator+( const double v , const Monotonic& p ) { return p + v  ; }
    ///  Constant times Mononitic
    inline Bernstein operator*( const double v , const Monotonic& p ) { return p * v  ; }
    ///  Constant minus Monotonic 
    inline Bernstein operator-( const double v , const Monotonic& p )
    { return v - p.bernstein() ; }
    // ========================================================================
    inline Bernstein operator+ ( const Monotonic& a , const Bernstein& b ) 
    { return a.bernstein () + b ; }
    inline Bernstein operator- ( const Monotonic& a , const Bernstein& b ) 
    { return a.bernstein () - b ; }
    inline Bernstein operator+ ( const Bernstein& b , const Monotonic& a ) { return a + b ; }
    inline Bernstein operator- ( const Bernstein& b , const Monotonic& a ) 
    { return b - a.bernstein () ; }
    inline Bernstein operator* ( const Monotonic& a , const Bernstein& b ) 
    { return a.bernstein () * b ; }
    inline Bernstein operator* ( const Bernstein& b , const Monotonic& a ) { return a * b ; }
    // ========================================================================
    /// Swapping function for positive monotonic polynomials 
    inline void swap ( Monotonic&     a , Monotonic&     b ) { a.swap ( b ) ; }
    // ========================================================================
    /** @class Convex
     *  The "positive" polynomial of order N with
     *  fixed sign of first and second derivatives
     *
     *  Conceptually, the monotonic increasing convex polynomial is parameterised as 
     *  \f[ C^{(I)}_n (x) = \cos^2 \phi_0 A \int_0^x \int_0^y P_{n-2}(z) dz dy + \sin^2 \phi_0 I_1(x) \f]
     *  where \f$P_n(x)\f$ is a positive polynomial and 
     *  \f$ I_n(x)\f$ is positive monotonically increasing polynomial 
     *
     *  The parameters \f$ \phi_i \f$ where \f$ 0 \le i < N \f$, are 
     *  -  parameter \f$ \phi_0 \f$ is used to parameterize 
     *  -  parameter \f$ \phi_1 \f$ is used to parameterize \f$ I_1 (x) \f$ 
     *  -  parameters \f$ \phi_{2...}\f$ ar eused to parameterize \f$ P_{n-1}(x)\f$ 
     *
     *  @see Ostap::Math::Monotonic 
     *  @see Ostap::Math::Positive
     *  @see Ostap::Math::Bernstein 
     */
    class Convex 
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from the order
       *  @param N  degree of polynomial 
       *  @param xmin  low-edge 
       *  @param xmax  high-edge 
       *  @param increasing increasing or decreasing ? 
       *  @param convex     convex or concave? 
       */
      Convex
      ( const unsigned short       N          =    1 ,
        const double               xmin       =    0 ,
        const double               xmax       =    1 ,
        const bool                 increasing = true ,
        const bool                 convex     = true ) ;
      // ======================================================================
      /** constructor from parameters/phases 
       *  @param pars list of parameters/phases 
       *  @param xmin  low-edge 
       *  @param xmax  high-edge 
       *  @param increasing increasing or decreasing ? 
       *  @param convex     convex or concave? 
       */ 
      Convex
      ( const std::vector<double>& pars              ,
        const double               xmin       =    0 ,
        const double               xmax       =    1 ,
        const bool                 increasing = true ,
        const bool                 convex     = true ) ;
      // ======================================================================
      /** constructor from the sequence of parameters 
       *  @param begin  start-iterator for sequence of coefficients 
       *  @param end     start-iterator for sequence of coefficients 
       *  @param xmin low-edge 
       *  @param xmax high-edge 
       *  @param increasing increasing or decreasing ? 
       *  @param convex     convex or concave? 
       */
      template<typename ITERATOR,
               typename value_type = typename std::iterator_traits<ITERATOR>::value_type ,
               typename = std::enable_if<std::is_convertible<value_type,long double>::value> >
      Convex ( ITERATOR     begin             , 
               ITERATOR     end               , 
               const double xmin       = 0    ,
               const double xmax       = 1    , 
               const bool   increasing = true ,
               const bool   convex     = true )
        : Convex ( std::distance (  begin , end ) , xmin , xmax , increasing , convex ) 
      { setPars  ( begin , end ) ; }
    // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double operator () ( const double x ) const { return m_bernstein ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get number of parameters
      std::size_t npars () const { return m_bernstein.degree() ; }
      /// set k-parameter
      bool setPar       ( const unsigned short k , const double value ) 
      {
        const unsigned short nA = m_sphereA .npars () ;
        const unsigned short nI = m_sphereI .npars () ;
        const unsigned short nP = m_positive.npars () ;
        const bool update = 
          k < nA           ? m_sphereA .setPar ( k           , value ) :
          k < nA + nI      ? m_sphereI .setPar ( k - nA      , value ) : 
          k < nA + nI + nP ? m_positive.setPar ( k - nA - nI , value ) : false ;
        return update ? updateBernstein () : false ;
      }
      /// set k-parameter
      bool setParameter ( const unsigned short k , const double value )
      { return setPar   ( k , value ) ; }
      /// get the parameter value
      double  par       ( const unsigned short k ) const
      {
        const unsigned short nA = m_sphereA .npars () ;
        const unsigned short nI = m_sphereI .npars () ;
        const unsigned short nP = m_positive.npars () ;
        return
          k < nA           ? m_sphereA .par ( k      ) : 
          k < nA + nI      ? m_sphereI .par ( k - nA ) :
          k < nA + nI + nP ? m_positive.par ( k - nA ) : 0.0 ;
      }
      // ================================================================
      /** set several/all parameters at once 
       *  @param begin  start itertaor for the sequence of coefficients 
       *  @param end    end   itertaor for the sequence of coefficients 
       *  @return true if at least one parameter is actually changed 
       */
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type ,
                typename = std::enable_if<std::is_convertible<value_type,long double>::value> >
      inline bool setPars 
      ( ITERATOR begin , 
        ITERATOR end   ) 
      {
        const auto NN = std::distance ( begin  , end ) ;
        const unsigned short nA = m_sphereA .npars () ;
        const unsigned short nI = m_sphereI .npars () ;
        const bool updatedA =                m_sphereA .setPars ( begin           , end ) ;
        const bool updatedI = nA      < NN ? m_sphereI .setPars ( begin + nA      , end ) : false ;
        const bool updatedP = nA + nI < NN ? m_positive.setPars ( begin + nA + nI , end ) : false ;
        return updatedA || updatedI || updatedP ? updateBernstein ( ) : false ;
      }
      // ======================================================================
      /** set several/all parameters at once 
       *  @param pars (NIPUT) vector of parameters 
       *  @return true if at least one parameter is actually changed 
       */
      inline bool setPars ( const std::vector<double>& pars ) 
      { return setPars ( pars.begin() , pars.end() ) ; }
      // ======================================================================
      /// get all parameters (by value) 
      std::vector<double>        pars  () const ;
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
      // ======================================================================
      /// convex     ?
      bool convex     () const { return degree () < 2 ||  m_convex     ; }
      /// convex     ?
      bool concave    () const { return degree () < 2 || !m_convex     ; }
      /// increasing ?
      bool increasing () const { return degree () < 1 ||  m_increasing ; }
      /// decreasing ?
      bool decreasing () const { return degree () < 1 || !m_increasing ; }
      /// monotonic
      bool monotonic  () const { return  true  ; }
      //// constant 
      bool constant   () const { return m_bernstein.constant () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the minimal value of function
      double fun_min () const ; // get the minimal value of function
      /// get the maximal value of function
      double fun_max () const ; // get the maximal value of function
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between xmin and xmax
      double integral   () const { return 1 ; }
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
      /// get the parameter sphere: alpha  & beta 
      const Ostap::Math::NSphere&   asphere   () const { return m_sphereA   ; }
      /// get the parameter sphere: linear integration "constant"
      const Ostap::Math::NSphere&   isphere   () const { return m_sphereI   ; }
      /// get the indefinite integral
      Bernstein indefinite_integral ( const double C = 0 ) const
      { return m_bernstein.indefinite_integral ( C ) ; }
      /// get the derivative
      Bernstein derivative          () const
      { return m_bernstein.derivative          () ; }
      // ======================================================================
    public:  /// basic operations  for python
      // ======================================================================
      /// Sum of Bernstein polynomial and a constant
      Bernstein __add__     ( const double value ) const { return m_bernstein + value   ; }
      /// Sum of Bernstein polynomial and a constant
      Bernstein __radd__    ( const double value ) const { return m_bernstein + value   ; }
      /// Product of Bernstein polynomial and a constant
      Bernstein __mul__     ( const double value ) const { return m_bernstein * value   ; }
      /// Product of Bernstein polynomial and a constant
      Bernstein __rmul__    ( const double value ) const { return m_bernstein * value   ; }
      /// Subtract a constant from Benrstein polynomial
      Bernstein __sub__     ( const double value ) const { return m_bernstein - value   ; }
      /// Constant minus Bernstein polynomial
      Bernstein __rsub__    ( const double value ) const { return value - m_bernstein   ; }
      /// Divide Bernstein polynomial by a constant
      Bernstein __truediv__ ( const double value ) const { return m_bernstein / value   ; }
      Bernstein __div__     ( const double value ) const { return __truediv__ ( value ) ; }
      /// Negate Bernstein polynomial
      Bernstein __neg__     () const { return -m_bernstein ; }
      // ======================================================================
      Bernstein __add__  ( const Bernstein& a ) const { return m_bernstein + a ; }
      Bernstein __radd__ ( const Bernstein& a ) const { return a + m_bernstein ; }
      Bernstein __sub__  ( const Bernstein& a ) const { return m_bernstein - a ; }
      Bernstein __rsub__ ( const Bernstein& a ) const { return a - m_bernstein ; }
      Bernstein __mul__  ( const Bernstein& a ) const { return m_bernstein * a ; }
      Bernstein __rmul__ ( const Bernstein& a ) const { return a * m_bernstein ; }
       // ======================================================================
    public:
      // ======================================================================
      /// get the tag
      std::size_t tag () const { return m_bernstein.tag () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// swap two objects 
      void swap ( Convex& right ) 
      {
        Ostap::Math::swap ( m_bernstein  , right.m_bernstein  ) ;
        Ostap::Math::swap ( m_positive   , right.m_positive   ) ;
        Ostap::Math::swap ( m_sphereA    , right.m_sphereA    ) ;
        Ostap::Math::swap ( m_sphereI    , right.m_sphereI    ) ;
        std::swap         ( m_increasing , right.m_increasing ) ;
        std::swap         ( m_convex     , right.m_convex     ) ;
        std::swap         ( m_aux        , right.m_aux        ) ;
      } 
      // ======================================================================
    private :
      // ======================================================================
      /// update bernstein coefficients
      bool updateBernstein () ;
      // ======================================================================
    protected:
      // ======================================================================
      /// the actual Bernstein polynomial
      Ostap::Math::Bernstein           m_bernstein  ; // the actual bernstein polynomial
      /// helper posivity polymnomial 
      Ostap::Math::Positive            m_positive   ; // helper positive polymnomial
      /// parameters sphere: alpha & beta 
      Ostap::Math::NSphere             m_sphereA    ; // alpha & beta 
      /// parameters sphere: integration linear function 
      Ostap::Math::NSphere             m_sphereI    ;
      /// increasing or decreasing ?
      bool                             m_increasing ; // increasing or decreasing?
      /// convex or concave ?
      bool                             m_convex     ; // convex or concave ?
      /// helper worklspace 
      mutable std::vector<long double> m_aux        ; // helper workspace 
      // ======================================================================
    } ;
    // ========================================================================
    ///  Convex plus      constant
    inline Bernstein operator+( const Convex& p , const double v )
    { return p.bernstein() + v ; }
    ///  Convex multiply  constant
    inline Bernstein operator*( const Convex& p , const double v )
    { return p.bernstein() * v ; }
    ///  Convex minus     constant
    inline Bernstein operator-( const Convex& p , const double v )
    { return p.bernstein() - v ; }
    ///  Convex divide constant
    inline Bernstein operator/( const Convex& p , const double v )
    { return p.bernstein() / v ; }
    ///  Consstant plus  Convex
    inline Bernstein operator+( const double v , const Convex& p ) { return p + v  ; }
    ///  Constant  times Convex
    inline Bernstein operator*( const double v , const Convex& p ) { return p * v  ; }
    ///  Constant  minus Convex
    inline Bernstein operator-( const double v , const Convex& p )
    { return v - p.bernstein() ; }
    // ========================================================================
    inline Bernstein operator+ ( const Convex&  a , const Bernstein& b ) 
    { return a.bernstein () + b ; }
    inline Bernstein operator- ( const Convex&  a , const Bernstein& b ) 
    { return a.bernstein () - b ; }
    inline Bernstein operator+ ( const Bernstein& b , const Convex&  a ) { return a + b ; }
    inline Bernstein operator- ( const Bernstein& b , const Convex&  a ) 
    { return b - a.bernstein () ; }
    inline Bernstein operator* ( const Convex&    a , const Bernstein& b ) 
    { return a.bernstein () * b ; }
    inline Bernstein operator* ( const Bernstein& b , const Convex&    a ) { return a * b ; }
    // ========================================================================
    /// Swapping function for positive monotonic convex/concave polynomials 
    inline void swap ( Convex&        a , Convex&        b ) { a.swap ( b ) ; }
    // ========================================================================
    /** @class ConvexOnly
     *  The "positive" polynomial of order N with
     *  fixed sign the second derivatives
     *  Actually it is a sum of basic bernstein polynomials with
     *  non-negative coefficients
     */
    class ConvexOnly 
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
      bool setPar       ( const unsigned short k , const double value ) 
      { 
        const bool update = m_sphere.setPhase ( k , value ) ;
        return update ? updateBernstein () : false ;
      }
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
      // ======================================================================
      /// convex     ?
      bool convex     () const { return degree () < 2 ||  m_convex     ; }
      /// convex     ?
      bool concave    () const { return degree () < 2 || !m_convex     ; }
      //// constant 
      bool constant   () const { return m_bernstein.constant   () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between xmin and xmax
      double integral   () const { return 1 ; }
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
    public:  /// basic operations  for python
      // ======================================================================
      /// Sum of Bernstein polynomial and a constant
      Bernstein __add__     ( const double value ) const { return m_bernstein + value   ; }
      /// Sum of Bernstein polynomial and a constant
      Bernstein __radd__    ( const double value ) const { return m_bernstein + value   ; }
      /// Product of Bernstein polynomial and a constant
      Bernstein __mul__     ( const double value ) const { return m_bernstein * value   ; }
      /// Product of Bernstein polynomial and a constant
      Bernstein __rmul__    ( const double value ) const { return m_bernstein * value   ; }
      /// Subtract a constant from Benrstein polynomial
      Bernstein __sub__     ( const double value ) const { return m_bernstein - value   ; }
      /// Constant minus Bernstein polynomial
      Bernstein __rsub__    ( const double value ) const { return value - m_bernstein   ; }
      /// Divide Bernstein polynomial by a constant
      Bernstein __truediv__ ( const double value ) const { return m_bernstein / value   ; }
      Bernstein __div__     ( const double value ) const { return __truediv__ ( value ) ; }
      /// Negate Bernstein polynomial
      Bernstein __neg__     () const { return -m_bernstein ; }
      // ======================================================================
      Bernstein __add__  ( const Bernstein& a ) const { return m_bernstein + a ; }
      Bernstein __radd__ ( const Bernstein& a ) const { return a + m_bernstein ; }
      Bernstein __sub__  ( const Bernstein& a ) const { return m_bernstein - a ; }
      Bernstein __rsub__ ( const Bernstein& a ) const { return a - m_bernstein ; }
      Bernstein __mul__  ( const Bernstein& a ) const { return m_bernstein * a ; }
      Bernstein __rmul__ ( const Bernstein& a ) const { return a * m_bernstein ; }
      // ======================================================================
      public:
      // ======================================================================
      /// get the tag
      std::size_t tag () const { return m_bernstein.tag () ; }
      // ======================================================================
     public:
      // ======================================================================
      /// swap two objects 
      void swap ( ConvexOnly& right ) 
      {
        Ostap::Math::swap ( m_bernstein  , right.m_bernstein  ) ;
        Ostap::Math::swap ( m_sphere     , right.m_sphere     ) ;
        std::swap         ( m_convex     , right.m_convex     ) ;
      } 
      // ======================================================================
   private :
      // ======================================================================
      /// update bernstein coefficients
      bool updateBernstein () ;
      // ======================================================================
    private :
      // ======================================================================
      /// the actual bernstein polynomial
      Ostap::Math::Bernstein m_bernstein  ; // the actual bernstein polynomial
      /// parameters sphere
      Ostap::Math::NSphere   m_sphere     ;
      /// convex ?
      bool                   m_convex     ; // convex ?
      // ======================================================================
    } ;
    // ========================================================================
    ///  ConvexOnly plus      constant
    inline Bernstein operator+( const ConvexOnly& p , const double v )
    { return p.bernstein() + v ; }
    ///  ConvexOnly multiply  constant
    inline Bernstein operator*( const ConvexOnly& p , const double v )
    { return p.bernstein() * v ; }
    ///  ConvexOnly minus     constant
    inline Bernstein operator-( const ConvexOnly& p , const double v )
    { return p.bernstein() - v ; }
    ///  ConvexOnly divide constant
    inline Bernstein operator/( const ConvexOnly& p , const double v )
    { return p.bernstein() / v ; }
    ///  Constant plus  ConvexOnly
    inline Bernstein operator+( const double v , const ConvexOnly& p ) { return p + v  ; }
    ///  Constant times ConvexOnly
    inline Bernstein operator*( const double v , const ConvexOnly& p ) { return p * v  ; }
    ///  Constant minus ConvexOnly
    inline Bernstein operator-( const double v , const ConvexOnly& p )
    { return v - p.bernstein() ; }
    // ========================================================================
    inline Bernstein operator+ ( const ConvexOnly&  a , const Bernstein& b ) 
    { return a.bernstein () + b ; }
    inline Bernstein operator- ( const ConvexOnly&  a , const Bernstein& b ) 
    { return a.bernstein () - b ; }
    inline Bernstein operator+ ( const Bernstein& b , const ConvexOnly&  a ) { return a + b ; }
    inline Bernstein operator- ( const Bernstein& b , const ConvexOnly&  a ) 
    { return b - a.bernstein () ; }
    inline Bernstein operator* ( const ConvexOnly& a , const Bernstein&  b ) 
    { return a.bernstein () * b ; }
    inline Bernstein operator* ( const Bernstein&  b , const ConvexOnly& a ) { return a * b ; }
    // ========================================================================
    /// Swapping function for positive convex/concave polynomials 
    inline void swap ( ConvexOnly&    a , ConvexOnly&    b ) { a.swap ( b ) ; }
    // ========================================================================
    namespace Utils 
    {
      // ======================================================================
      /** get "positive-pseudo-roots" for the positive polynomial in 
       *  a form by S.Karlin and L.S. Shapley 
       *  such choice  of roots gives flat polynomial.
       *  The choice of roots is motivated by two identities 
       *  \f[ \begin{array}{l}  
       *       T^2_{\mathrm{n}}(x) + \left(1-x^2\right) U^2_{\mathrm{n-1}}(x) = 1 \\ 
       *     \left(1+x\right)V^2_{\mathrm{n}}(x) + \left(1-x\right) W^2_{\mathrm{n}}(x) = 1
       *     \end{array}\f]
       *  where 
       *   -  \f$ T_{\mathrm{n}} \f$ is chebyshev polynomial of the 1st  kind,
       *   -  \f$ U_{\mathrm{n}} \f$ is chebyshev polynomial of the 2nd  kind,
       *   -  \f$ V_{\mathrm{n}} \f$ is chebyshev polynomial of the 3rd  kind,
       *   -  \f$ W_{\mathrm{n}} \f$ is chebyshev polynomial of the 4th  kind
       *
       *  With such "pseudo-roots" one  has 
       *  - for even N: \f$ \alpha   B_1(s) + (1-alpha) x ( 1 - x ) B_2(s) = 1 \f$
       *  - for odd N : \f$ \alpha x B_1(s) + (1-alpha)   ( 1 - x ) B_2(s) = 1 \f$  
       *
       *  where \f$ B_1(x)\f$ is a normalized polynomial that has 
       *  roots  \f$ r_0, r_0, r_2, r_2, ... \f$, and 
       *  B_2(x) is a normalized polynomial that as root
       *  \f$ r_1, r_1, r_3, r_3, ...\f$
       * 
       *  - The positivity is calculated for \f$ 0 \le x \le 1\f$ interval 
       *  - All  pseudo-rootsbelons to this interval 
       *  - \f$ B_{1,2}(s)\f$ are normalized as \f$ \int_0^1 B_i dx = 1 \f$
       *
       *  @param N  (NIPUT) polynomial degree 
       *  @param pproots (UPDATE) positive pseudo-roots 
       *  @return parameter \f$ \alpha \f$
       */
      double positive_pseudo_roots 
      ( const unsigned short N       , 
        std::vector<double>& pproots ) ;
      // ======================================================================
    }
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_BERNSTEIN1D_H
// ============================================================================
