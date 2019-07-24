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
     *  A special case of BErnstein polynomial with symmetry:
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
       *  the actual degree of polynomial will be 2*N
       *  @param N  parameter that defiend the order of polynomial (2*N)
       *  @param xmin low edge
       *  @param xmax high edge
       */
      BernsteinEven ( const unsigned short N     = 0 ,
                      const double         xmin  = 0 ,
                      const double         xmax  = 1 ) ;
      /** constructor from list of coefficients
       *  @param pars vector of parameters 
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
      bool           zero  () const { return m_bernstein.zero() ; }
      // ======================================================================
      /** set k-parameter
       *  @param k index
       *  @param value new value
       *  @return true if parameter is actually changed
       */
      bool setPar          ( const unsigned short k , const double value ) ;
      // ======================================================================
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
      /// get the tag 
      std::size_t tag () const { return m_bernstein.tag () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// swap two objects 
      void swap ( BernsteinEven& right ) 
      {
        std::swap         ( m_N         , right.m_N         ) ;
        Ostap::Math::swap ( m_bernstein , right.m_bernstein ) ;
      } 
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
    inline Bernstein operator+ ( const BernsteinEven&  a , const Bernstein& b ) 
    { return a.bernstein () + b ; }
    inline Bernstein operator- ( const BernsteinEven&  a , const Bernstein& b ) 
    { return a.bernstein () - b ; }
    inline Bernstein operator+ ( const Bernstein& b , const BernsteinEven&  a ) { return a + b ; }
    inline Bernstein operator- ( const Bernstein& b , const BernsteinEven&  a ) 
    { return b - a.bernstein () ; }
    inline Bernstein operator* ( const BernsteinEven&  a , const Bernstein& b ) 
    { return a.bernstein () * b ; }
    inline Bernstein operator* ( const Bernstein& b , const BernsteinEven&  a ) { return a * b ; }
    // ========================================================================
    /// Swapping function for even bernstein polynomials 
    inline void swap ( BernsteinEven& a , BernsteinEven& b ) { a.swap ( b ) ; }
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
    class Positive 
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
      /// decreasing?
      bool decreasing  () const { return m_bernstein.decreasing () ; }
      /// increasing?
      bool increasing  () const { return m_bernstein.increasing () ; }
      /// monotonic?
      bool monotonic   () const { return m_bernstein.monotonic  () ; }
      /// constant 
      bool constant    () const { return m_bernstein.constant   () ; }
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
      void swap ( Positive& right ) 
      {
        Ostap::Math::swap ( m_bernstein , right.m_bernstein ) ;
        Ostap::Math::swap ( m_sphere    , right.m_sphere    ) ;
      } 
      // ======================================================================
    private : 
      // ======================================================================
      /// update bernstein coefficiencts
      bool updateBernstein () ;
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
      bool setPar       ( const unsigned short k , const double value )        
      {
        const bool update = m_sphere.setPhase ( k , value ) ;
        return update ? updateBernstein () : false ; ;
      }
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
      double integral   () const { return 1 ; }
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
        Ostap::Math::swap ( m_sphere , right.m_sphere ) ;
        Ostap::Math::swap ( m_even   , right.m_even   ) ;
      } 
      // ======================================================================
    private: 
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
     *  Actually it is a sum of basic bernstein polynomials with
     *  non-negative coefficients that form the monotonic sequence 
     */
    class Monotonic
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      Monotonic
      ( const unsigned short       N          =    1 ,
        const double               xmin       =    0 ,
        const double               xmax       =    1 ,
        const bool                 increasing = true ) ;
      // ======================================================================
      /// constructor from N phases
      Monotonic
      ( const std::vector<double>& pars              ,
        const double               xmin       =    0 ,
        const double               xmax       =    1 ,
        const bool                 increasing = true ) ;
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
      /// increasing ?
      bool increasing () const { return degree() < 1 ||  m_increasing ; }
      /// decreasing ?
      bool decreasing () const { return degree() < 1 || !m_increasing ; }
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
        Ostap::Math::swap ( m_sphere     , right.m_sphere     ) ;
        std::swap         ( m_increasing , right.m_increasing ) ;
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
      Ostap::Math::Bernstein m_bernstein  ; // the actual bernstein polynomial
      /// parameters sphere
      Ostap::Math::NSphere   m_sphere     ;
      /// increasing ?
      bool                   m_increasing ; // increasing ?
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
     *  Actually it is a sum of basic bernstein polynomials with
     *  non-negative coefficients
     */
    class Convex 
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
      void swap ( Convex& right ) 
      {
        Ostap::Math::swap ( m_bernstein  , right.m_bernstein  ) ;
        Ostap::Math::swap ( m_sphere     , right.m_sphere     ) ;
        std::swap         ( m_increasing , right.m_increasing ) ;
        std::swap         ( m_convex     , right.m_convex     ) ;
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
      Ostap::Math::Bernstein m_bernstein  ; // the actual bernstein polynomial
      /// parameters sphere
      Ostap::Math::NSphere   m_sphere     ;
      /// increasing ?
      bool                   m_increasing ; // increasing ?
      /// convex ?
      bool                   m_convex     ; // convex ?
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
      /// increasing ?
      bool increasing () const { return m_bernstein.increasing () ; }
      /// decreasing ?
      bool decreasing () const { return m_bernstein.decreasing () ; }
      /// monotonic
      bool monotonic  () const { return m_bernstein.monotonic  () ; }
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
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_BERNSTEIN1D_H
// ============================================================================
