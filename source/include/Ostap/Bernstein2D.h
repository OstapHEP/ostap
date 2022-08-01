#ifndef OSTAP_BERNSTEIN2D_H 
#define OSTAP_BERNSTEIN2D_H 1
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
#include "Ostap/Parameters.h"
#include "Ostap/NSphere.h"
#include "Ostap/Bernstein.h"
// ============================================================================
/** @file Ostap/Bernstein.h
 *  Set of useful math-functions, related to 2D-Bernstein polynomials
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
    class Bernstein2D   ;
    class Bernstein2DSym;    
    // ========================================================================
    /** @class Bernstein2D
     *  The Bernstein's polynomial of order Nx*Ny
     *  \f[  B_{n^x,n^y}(x,y) \equiv 
     *    \sum_{i=0}{n^{x}}
     *    \sum_{j=0}{n^{y}} \alpha_{i,j} B_{n^{x}}^i(x) B_{n^{y}}^j(y) \f]
     *  where 
     *  - \f$ B_n^k(x) \f$ is basic Bernstein polynomial
     *  @see Ostap::Math::Bernstein
     */
    class Bernstein2D : public Ostap::Math::Parameters 
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      Bernstein2D 
      ( const unsigned short       nX    =  1 ,
        const unsigned short       nY    =  1 ,
        const double               xmin  =  0 ,
        const double               xmax  =  1 ,
        const double               ymin  =  0 ,
        const double               ymax  =  1 ) ;
      // ======================================================================
      /** As a product of two 1D-polynomials:
       *  \f[  B_{n^x,n^y}(x,y) \equiv 
       *      B^{n^x}(x)B^{n^y}(y) = 
       *  \left(\sum_{i=0}{n^{x}} \alpha_{i} B_{n^{x}}^i(x)]\right)
       *  \left(\sum_{j=0}{n^{y}} \beta_{j} B_{n^{y}}^j(y)]\right) = 
       *    \sum_{i=0}{n^{x}}
       *    \sum_{j=0}{n^{y}} \alpha_{i}\beta_{j} B_{n^{x}}^i(x) B_{n^{y}}^j(y) \f]
       */          
      Bernstein2D
      ( const Bernstein& bx , 
        const Bernstein& by ) ;
      // ======================================================================
      /// from symmetric variant 
      Bernstein2D 
      ( const Bernstein2DSym& right ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double evaluate    ( const double x , const double y ) const ;
      /// get the value
      double operator () ( const double x , const double y ) const 
      { return evaluate ( x , y ) ; }
      // ======================================================================
    public: // getters & setters 
      // ======================================================================
      using Ostap::Math::Parameters::par    ;
      using Ostap::Math::Parameters::setPar ;
      // ======================================================================      
      /// get (l,m)-parameter
      double  par       
      ( const unsigned short l ,
        const unsigned short m ) const 
      { return l > m_nx ? 0 : m > m_ny ? 0 : par ( index ( l , m ) ) ; }
      // ======================================================================
      bool setPar   
      ( const unsigned short l     ,
        const unsigned short m     ,
        const double         value ) 
        
      { 
        return 
          l > m_nx ? false :
          m > m_ny ? false : setPar ( index ( l , m ) , value ) ; 
      }
      // ======================================================================
    public:
      // ======================================================================
      ///  convert (l,m)-index into single k-index
      unsigned int index 
      ( const unsigned short l , 
        const unsigned short m ) const 
      {
        return
          l > m_nx ? -1 :
          m > m_ny ? -1 :
          1u * l * ( m_ny + 1 ) + m ;  
      }
      // ======================================================================
    public:
      // ======================================================================
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
      double integral  
      ( const double xlow  , 
        const double xhigh ,
        const double ylow  , 
        const double yhigh ) const ;
      // ======================================================================
      /** integral over x-dimension
       *  \f[ \int_{x_{low}}^{x_{high}} \mathcal{B}(x,y) \mathrm{d}x\f]
       *  @param y     variable
       *  @param xlow  low  edge in y
       *  @param xhigh high edge in y
       */
      double integrateX 
      ( const double y     ,
        const double xlow  , 
        const double xhigh ) const ;
      // ======================================================================
      /** integral over y-dimension
       *  \f[ \int_{y_{low}}^{y_{high}} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param x     variable
       *  @param ylow  low  edge in x
       *  @param yhigh high edge in x
       */
      double integrateY 
      ( const double x     ,
        const double ylow  , 
        const double yhigh ) const ;
      // ======================================================================
    public: // special cases
      // ======================================================================
      /** get the integral over 2D-region
       *  \f[  x_min < x < x_max, y_min< y< y_max\f]
       */
      double integral   () const ;
      /** integral over x-dimension
       *  \f[ \int_{x_{min}}^{x_{max}} \mathcal{B}(x,y) \mathrm{d}x\f]
       *  @param y     variable
       */
      double integrateX ( const double y    ) const ;
      /** integral over y-dimension
       *  \f[ \int_{y_{min}}^{y_{max}} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param x     variable
       */
      double integrateY ( const double x    ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// simple  manipulations with polynoms: shift it!
      Bernstein2D& operator += ( const double a ) ;
      /// simple  manipulations with polynoms: shift it!
      Bernstein2D& operator -= ( const double a ) ;
      /// simple  manipulations with polynoms: scale it!
      Bernstein2D& operator *= ( const double a ) ;
      /// simple  manipulations with polynoms: scale it!
      Bernstein2D& operator /= ( const double a ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// negate it!
      Bernstein2D  operator-() const ;
      // ======================================================================
    public: // python! 
      // ======================================================================
      /// Sum of Bernstein polynomial and a constant
      Bernstein2D __add__     ( const double value ) const ;
      /// Sum of Bernstein polynomial and a constant
      Bernstein2D __radd__    ( const double value ) const ;
      /// Product of Bernstein polynomial and a constant
      Bernstein2D __mul__     ( const double value ) const ;
      /// Product of Bernstein polynomial and a constant
      Bernstein2D __rmul__    ( const double value ) const ;
      /// Subtract a constant from Benrstein polynomial
      Bernstein2D __sub__     ( const double value ) const ;
      /// Constant minus Bernstein polynomial
      Bernstein2D __rsub__    ( const double value ) const ;
      /// Divide Benrstein polynomial by a constant
      Bernstein2D __truediv__ ( const double value ) const ;
      Bernstein2D __div__     ( const double value ) const { return __truediv__ ( value ) ;  }
      /// Negate Bernstein polynomial
      Bernstein2D __neg__   () const ;
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
    public:
      // ======================================================================
      /// swap  two 2D-polynomials 
      void swap ( Bernstein2D&  right ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag value 
      std::size_t tag () const ; // get the tag value 
      // ======================================================================
    private: // helper functions to make the calculations
      // ======================================================================
      /// helper function to make calculations
      double calculate
      ( const std::vector<double>& fx , 
        const std::vector<double>& fy ) const ;
      // ======================================================================
    private:
      // ======================================================================
      // polynom order in x-dimension
      unsigned short m_nx ; // polynom order in x-dimension
      // polynom order in y-dimension
      unsigned short m_ny ; // polynom order in y-dimension
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
    ///  Bernstein plus      constant
    inline Bernstein2D operator+( const Bernstein2D& p , const double v )
    { return Bernstein2D ( p ) += v ; } //  Bernstein plus constant
    ///  Bernstein multiply  constant
    inline Bernstein2D operator*( const Bernstein2D& p , const double v )
    { return Bernstein2D ( p ) *= v ; } //  Bernstein plus constant
    ///  Bernstein minus constant
    inline Bernstein2D operator-( const Bernstein2D& p , const double v )
    { return Bernstein2D ( p ) -= v ; } //  Bernstein plus constant
    ///  Bernstein divide constant
    inline Bernstein2D operator/( const Bernstein2D& p , const double v )
    { return Bernstein2D ( p ) /= v ; } //  Bernstein plus constant
    ///  Constant plus  Bernstein
    inline Bernstein2D operator+( const double v , const Bernstein2D& p ) 
    { return p +   v  ; }
    ///  Constant times Bernstein
    inline Bernstein2D operator*( const double v , const Bernstein2D& p ) 
    { return p *   v  ; }
    ///  Constant minus Bernstein
    inline Bernstein2D operator-( const double v , const Bernstein2D& p ) 
    { return v + (-p) ; }
    // ========================================================================
    /// swap  two 2D-polynomials 
    inline void  swap ( Bernstein2D& a  , Bernstein2D& b ) { a.swap ( b ) ; }
    // ========================================================================
    /** @class Positive2D
     *  The "positive" 2D-polynomial of order Nx*Ny
     *  Actually it is a sum of basic bernstein 2D-polynomials with
     *  non-negative coefficients
     *  \f[  P_{n^x,n^y}(x,y) \equiv 
     *    \sum_{i=0}{n^{x}}
     *    \sum_{j=0}{n^{y}} \alpha_{i,j} B_{n^{x}}^i(x) B_{n^{y}}^j(y) \f]
     *  where 
     *  - \f$ B_n^k(x) \f$ is basic Bernstein polynomial
     *  - \f$ \alpha_{ij} \ge 0\f$
     *  - \f$ \sum{i,j} \alpha_{i,j} = 1 \f$ 
     *
     *  Clearly \f$ P_{N^x,N^y}(x,y) \ge 0\f$  
     *  @see Ostap::Math::Bernstein2D
     */
    class Positive2D 
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      Positive2D 
      ( const unsigned short       Nx    =  1 ,
        const unsigned short       Ny    =  1 ,
        const double               xmin  =  0 ,
        const double               xmax  =  1 ,
        const double               ymin  =  0 ,
        const double               ymax  =  1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      ///  get the  value 
      double evaluate    ( const double x , const double y ) const 
      { return m_bernstein ( x , y ) ; }        
      /// get the value
      double operator () ( const double x , const double y ) const
      { return evaluate    ( x , y ) ; }
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
      /// get all parameters/phases 
      const  std::vector<double>& pars() const 
      { return m_sphere.phases() ; }
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
      double integral
      ( const double xlow  , 
        const double xhigh ,
        const double ylow  , 
        const double yhigh ) const ;
      /** integral over x-dimension
       *  \f[ \int_{x_{low}}^{x_{high}} \mathcal{B}(x,y) \mathrm{d}x\f]
       *  @param y     variable
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       */
      double integrateX 
      ( const double y     ,
        const double xlow  ,
        const double xhigh ) const
      { return m_bernstein.integrateX ( y , xlow , xhigh ) ; }
      /** integral over y-dimension
       *  \f[ \int_{y_{low}}^{y_{high}} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param x     variable
       *  @param ylow  low  edge in x
       *  @param yhigh high edge in x
       */
      double integrateY 
      ( const double x     ,
        const double ylow  ,
        const double yhigh ) const
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
       *  @param y     variable
       */
      double integrateX ( const double y    ) const
      { return m_bernstein.integrateX ( y ) ; }
      /** integral over y-dimension
       *  \f[ \int_{y_{min}}^{y_{max}} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param x     variable
       */
      double integrateY ( const double x    ) const
      { return m_bernstein.integrateY ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// swap  two 2D-polynomials 
      void swap ( Positive2D&  right ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const { return m_bernstein.tag () ; } 
      // ======================================================================
    public: // ingeredients
      // ======================================================================
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
    /// p + v 
    inline Bernstein2D operator+( const Positive2D& p , const double      v )
    { return p.bernstein() + v ; } 
    inline Bernstein2D operator*( const Positive2D& p , const double      v )
    { return p.bernstein() * v ; } 
    inline Bernstein2D operator-( const Positive2D& p , const double      v )
    { return p.bernstein() - v ; } 
    inline Bernstein2D operator/( const Positive2D& p , const double      v )
    { return p.bernstein() / v ; } 
    inline Bernstein2D operator+( const double      v , const Positive2D& p )
    { return p +      v ; } 
    inline Bernstein2D operator*( const double      v , const Positive2D& p )
    { return p *      v ; } 
    inline Bernstein2D operator-( const double      v , const Positive2D& p )
    { return v + -1 * p ; } 
    // ========================================================================
    /// swap  two 2D-polynomials 
    inline void  swap ( Positive2D& a  , Positive2D& b ) { a.swap ( b ) ; }
    // ========================================================================
    /** @class Bernstein2DSym
     *  The symmetric Bernstein's polynomial of order N*N.
     *  Symmetric  verison of Ostap::Math::Bernstein2D
     *  
     *  \f[  B_{n}(x,y) \equiv 
     *    \sum_{i=0}{n}
     *    \sum_{j=0}{n} \alpha_{i,j} B_n^i(x) B_n^j(y) \f]
     *  where 
     *  - \f$ B_n^k(x) \f$ is basic Bernstein polynomial
     *  - \f$ \alpha_{ji} = \alpha_{ij}\f$ 
     * 
     *  Clearly  it is symmetric function \f$ B_{n}(y,x)=B_{n}(x,y) \f$
     *  @see Ostap::Math::Bernstein2D
     *  @see Ostap::Math::Bernstein
     */
    class Bernstein2DSym : public Ostap::Math::Parameters 
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      Bernstein2DSym
      ( const unsigned short       n     =  1 ,
        const double               xmin  =  0 ,
        const double               xmax  =  1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double evaluate    ( const double x , const double y ) const ;
      /// get the value
      double operator () ( const double x , const double y ) const 
      { return evaluate ( x , y ) ; }
      // ======================================================================
    public:
      // ======================================================================
      using Ostap::Math::Parameters::par    ;
      using Ostap::Math::Parameters::setPar ;
      // ======================================================================
      /// set (l,m)-parameter
      bool setPar      
      ( const unsigned short l     ,
        const unsigned short m     ,
        const double         value ) 
      {
        return
          l > m_n ? false : 
          m > m_n ? false : setPar ( index ( l , m ) , value ) ;
      }
      /// get (l,m)-parameter
      double  par       
      ( const unsigned short l ,
        const unsigned short m ) const 
      {
        return 
          l > m_n ? 0 : 
          m > m_n ? 0 : par ( index ( l , m ) ) ;
      }
      // ======================================================================
    public:
      // ======================================================================
      ///  convert (l,m)-index into single k-index
      unsigned int index 
      ( const unsigned short l , 
        const unsigned short m ) const 
      {
        return
          m > l   ? index ( m , l )  :
          l > m_n ? -1               :  // NB !!
          1u * l * ( l + 1 ) / 2 + m ;
      }
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
      double integral 
      ( const double xlow  , 
        const double xhigh ,
        const double ylow  , 
        const double yhigh ) const ;
      /** integral over x-dimension
       *  \f[ \int_{x_{low}}^{x_{high}} \mathcal{B}(x,y) \mathrm{d}x\f]
       *  @param y     variable
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       */
      double integrateX
      ( const double y     ,
        const double xlow  , 
        const double xhigh ) const ;
      /** integral over y-dimension
       *  \f[ \int_{y_{low}}^{x_{high}} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param x     variable
       *  @param ylow  low  edge in x
       *  @param yhigh high edge in x
       */
      double integrateY
      ( const double x     ,
        const double ylow  ,
        const double yhigh ) const ;
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
       *  @param y     variable
       */
      double integrateX ( const double y ) const ;
      /** integral over y-dimension
       *  \f[ \int_{y_{min}}^{x_{max}} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param x     variable
       */
      double integrateY ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// simple  manipulations with polynoms: shift it!
      Bernstein2DSym& operator += ( const double a ) ;
      /// simple  manipulations with polynoms: shift it!
      Bernstein2DSym& operator -= ( const double a ) ;
      /// simple  manipulations with polynoms: scale it!
      Bernstein2DSym& operator *= ( const double a ) ;
      /// simple  manipulations with polynoms: scale it!
      Bernstein2DSym& operator /= ( const double a ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// negate it!
      Bernstein2DSym  operator-() const ;
      // ======================================================================
    public:
      // ======================================================================
      /// Sum of Bernstein polynomial and a constant
      Bernstein2DSym __add__     ( const double value ) const ;
      /// Sum of Bernstein polynomial and a constant
      Bernstein2DSym __radd__    ( const double value ) const ;
      /// Product of Bernstein polynomial and a constant
      Bernstein2DSym __mul__     ( const double value ) const ;
      /// Product of Bernstein polynomial and a constant
      Bernstein2DSym __rmul__    ( const double value ) const ;
      /// Subtract a constant from Benrstein polynomial
      Bernstein2DSym __sub__     ( const double value ) const ;
      /// Constant minus Bernstein polynomial
      Bernstein2DSym __rsub__    ( const double value ) const ;
      /// Divide Benrstein polynomial by a constant
      Bernstein2DSym __truediv__ ( const double value ) const ;
      Bernstein2DSym __div__     ( const double value ) const { return __truediv__ ( value ) ; }
      /// Negate Bernstein polynomial
      Bernstein2DSym __neg__   () const ;
      // ======================================================================
    public: // few helper functions to expose internals
      // ======================================================================
      /// evaluate the basic polynomials
      double basic  ( const unsigned short i , const double         x ) const
      { return ( i > m_n || x < m_xmin || x < m_xmax ) ? 0.0 : m_b[i](x) ; }
      /// expose some internals
      const Bernstein& basic ( const unsigned short i ) const { return m_b[i] ; }
      // ======================================================================
    public:
      // ======================================================================
      /// swap  two 2D-polynomials 
      void swap ( Bernstein2DSym&  right ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag value 
      std::size_t tag () const ;
      // ======================================================================
    private: // helper functions to make the calculations
      // ======================================================================
      /// helper function to make calculations
      double calculate 
      ( const std::vector<double>& fx , 
        const std::vector<double>& fy ) const ;
      // ======================================================================
    private:
      // ======================================================================
      // polynom order
      unsigned short m_n  ; // polynom order in x-dimension
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
    ///  Bernstein plus      constant
    inline Bernstein2DSym operator+( const Bernstein2DSym& p , const double v )
    { return Bernstein2DSym ( p ) += v ; } //  Bernstein plus constant
    ///  Bernstein multiply  constant
    inline Bernstein2DSym operator*( const Bernstein2DSym& p , const double v )
    { return Bernstein2DSym ( p ) *= v ; } //  Bernstein plus constant
    ///  Bernstein minus constant
    inline Bernstein2DSym operator-( const Bernstein2DSym& p , const double v )
    { return Bernstein2DSym ( p ) -= v ; } //  Bernstein plus constant
    ///  Bernstein divide constant
    inline Bernstein2DSym operator/( const Bernstein2DSym& p , const double v )
    { return Bernstein2DSym ( p ) /= v ; } //  Bernstein plus constant
    ///  Constant plus  Bernstein
    inline Bernstein2DSym operator+( const double v , const Bernstein2DSym& p ) { return p +   v  ; }
    ///  Constant times Bernstein
    inline Bernstein2DSym operator*( const double v , const Bernstein2DSym& p ) { return p *   v  ; }
    ///  Constant minus Bernstein
    inline Bernstein2DSym operator-( const double v , const Bernstein2DSym& p ) { return v + (-p) ; }
    // ========================================================================
    /// swap  two 2D-polynomials 
    inline void  swap ( Bernstein2DSym& a  , Bernstein2DSym& b ) { a.swap ( b ) ; }
    // ========================================================================
    /** @class Positive2DSym
     *  The "positive" symmetrical polynomial of order Nx*Ny
     *  Actually it is a sum of basic bernstein 2D-polynomials with
     *  non-negative coefficients.
     *  Symmetrized version of Ostap::Math::Positive2D
     *  \f[  P_{n,n}(x,y) \equiv 
     *    \sum_{i=0}^{n}
     *    \sum_{j=0}^{n} \alpha_{i,j} B_n^i(x) B_n^j(y) \f]
     *  where 
     *  - \f$ B_n^k(x) \f$ is basic Bernstein polynomial
     *  - \f$ \alpha_{ij} \ge 0\f$
     *  - \f$ \alpha_{ji} = \alpha_{ij}\f$
     *  - \f$ \sum_{i,j} \alpha_{i,j} = 1 \f$ 
     *
     *  Clearly \f$ P_{n}(x,y) \ge 0\f$  and \f$ P_{n}(y,x)=P_{n}(x,y) \f$
     *  @see Ostap::Math::Positive2D
     *  @see Ostap::Math::Bernstein2DSym
     */
    class Positive2DSym
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      Positive2DSym 
      ( const unsigned short       Nx    =  1 ,
        const double               xmin  =  0 ,
        const double               xmax  =  1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double evaluate    ( const double x , const double y ) const 
      { return m_bernstein ( x ,  y ) ; }
      /// get the value
      double operator () ( const double x , const double y ) const 
      { return evaluate ( x , y ) ; }
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
      /// get all parameters/phases 
      const std::vector<double>& pars() const { return m_sphere.phases() ; }
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
      double integral  
      ( const double xlow  , 
        const double xhigh ,
        const double ylow  , 
        const double yhigh ) const ;
      // ======================================================================
      /** integral over x-dimension
       *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param y     variable
       *  @param xlow  low  edge in y
       *  @param xhigh high edge in y
       */
      double integrateX 
      ( const double y     ,
        const double xlow  ,
        const double xhigh ) const ;
      // ======================================================================
      /** integral over x-dimension
       *  \f[ \int_{x_low}^{x_high} \mathcal{B}(x,y) \mathrm{d}x\f]
       *  @param x     variable
       *  @param ylow  low  edge in x
       *  @param yhigh high edge in x
       */
      double integrateY 
      ( const double x     ,
        const double ylow  , 
        const double yhigh ) const ;
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
       *  @param y     variable
       */
      double integrateX ( const double y ) const ;
      /** integral over y-dimension
       *  \f[ \int_{y_{min}}^{y_{max}} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param x     variable
       */
      double integrateY ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// swap  two 2D-polynomials 
      void swap ( Positive2DSym&  right ) ;
      // ======================================================================
      /// get the tag value 
      std::size_t tag () const { return m_bernstein.tag() ; }
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
    /// swap  two 2D-polynomials 
    inline void  swap ( Positive2DSym& a  , Positive2DSym& b ) { a.swap ( b ) ; }
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_BERNSTEIN2D_H
// ============================================================================
