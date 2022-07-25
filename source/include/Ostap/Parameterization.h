// ============================================================================
#ifndef OSTAP_PARAMETERIZATION_H 
#define OSTAP_PARAMETERIZATION_H 1
// ============================================================================
// Include files
// ============================================================================
// STD/STL
// ============================================================================
#include <vector>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Polynomials.h"
#include "Ostap/Math.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    /// 1D-parameterization as sum of Legendre polynomials 
    class  LegendreSum ; // 1D-parameterization as sum of Legendre polynomials 
    // ========================================================================
    /** @class  LegendreSum2  Ostap/Parameterization.h
     *  2D-parameterization  as sum of Legendre polynomials 
     *  \f$ S(x,y) = \sum_{i0}^{n_x}\sum_{j=0}^{n_y} c_{i,j} P_i(x^\prime)P_j(y^\prime)\f$,
     * where:
     *  - \f$ x^{\prime} = \frac{2x-x_{min}-x_{max}}{x_{max}-x_{min}}\f$,
     *  - \f$ y^{\prime} = \frac{2y-y_{min}-y_{max}}{y_{max}-y_{min}}\f$. 
     *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
     *  @date 2019-06-30
     */
    class LegendreSum2 : public Parameters 
    {
    public:
      // ======================================================================
      /** constructor 
       * \f$ S(x,y) = \sum_{i=0}^{n_x}\sum_{j=0}^{n_y} c_{i,j}P_i(x^\prime)P_j(y^\prime)\f$,
       * where:
       *  - \f$ x^{\prime} = \frac{2x-x_{min}-x_{max}}{x_{max}-x_{min}}\f$,
       *  - \f$ y^{\prime} = \frac{2y-y_{min}-y_{max}}{y_{max}-y_{min}}\f$.
       */
      LegendreSum2 
      ( const unsigned short NX =  0 , 
        const unsigned short NY =  0 , 
        const double   xmin     = -1 , 
        const double   xmax     =  1 , 
        const double   ymin     = -1 , 
        const double   ymax     =  1 ) ;
      // ======================================================================
      /** constructor orm the product of two Legendre sums
       *  \f$ S(x,y) = S_x(x)\times S_y(y) \f$ 
       *  @param sx (INPUT) the first  Legendre sum 
       *  @param sy (INPUT) the second Legendre sum 
       */
      LegendreSum2 
      ( const LegendreSum& sx , 
        const LegendreSum& sy ) ;
      // ======================================================================
    public: 
      // ======================================================================
      /// get the value
      double evaluate    ( const double x , const double y ) const ;
      // ======================================================================
      /// get the value
      double operator () ( const double x , const double y ) const 
      { return 
          x < m_xmin || x > m_xmax ? 0 : 
          y < m_ymin || y > m_ymax ? 0 : evaluate ( x , y ) ; }
      // =====================================================================
    public: 
      // ======================================================================
      /// get lower edge
      double xmin  () const { return m_xmin ; }
      /// get upper edge
      double xmax  () const { return m_xmax ; }
      // ======================================================================
      /// get lower edge
      double ymin  () const { return m_ymin ; }
      /// get upper edge
      double ymax  () const { return m_ymax ; }
      // ======================================================================
    public:
      // ======================================================================
      std::size_t degree_x  () const { return m_NX ; }
      std::size_t degree_y  () const { return m_NY ; }
      std::size_t nx        () const { return m_NX ; }
      std::size_t ny        () const { return m_NY ; }
      // ======================================================================
    public:
      // ======================================================================
      using Parameters::par    ;
      using Parameters::setPar ;
      /// get 2D-parameter 
      inline double par 
      ( const unsigned short ix , 
        const unsigned short iy ) const 
      { return Parameters::par ( index ( ix , iy ) ) ;  }
      // ======================================================================
      /// set 2D-parameter
      bool setPar
      ( const unsigned short ix     ,
        const unsigned short iy     ,
        const double         value  ) 
      { return Parameters::setPar ( index ( ix , iy ) , value ) ; }
      // ======================================================================
    public:
      // ======================================================================
      double x  ( const double tx ) const
      { return  0.5 * ( tx * ( m_xmax - m_xmin ) +   m_xmax + m_xmin ) ; }
      double tx ( const double  x ) const
      { return (  2 *    x   - m_xmax - m_xmin ) / ( m_xmax - m_xmin ) ; }
      // ======================================================================
      double y  ( const double ty ) const
      { return  0.5 * ( ty * ( m_ymax - m_ymin ) +   m_ymax + m_ymin ) ; }
      double ty ( const double  y ) const
      { return (  2 *    y   - m_ymax - m_ymin ) / ( m_ymax - m_ymin ) ; }
      // ======================================================================
    public: 
      // ======================================================================
      /** update  the Legendre expansion by addition of one "event" with 
       *  the given weight
       *  @code
       *  LegendreSum2 sum = ... ;
       *  for ( auto x : .... ) { sum.fill ( x , y ) ; }
       *  @endcode
       *  This is a useful function to make an unbinned parameterization 
       *  of certain distribution and/or efficiency 
       */
      bool fill
      ( const double x          , 
        const double y          , 
        const double weight = 1 ) ;
      // ======================================================================
    public: // several useful operators and operations 
      // ======================================================================
      LegendreSum2  operator+  ( const double b ) const 
      { LegendreSum2 c { *this } ; c += b ; return c ; }
      LegendreSum2  operator-  ( const double b ) const 
      { LegendreSum2 c { *this } ; c -= b ; return c ; }
      LegendreSum2  operator*  ( const double b ) const 
      { LegendreSum2 c { *this } ; c *= b ; return c ; }
      LegendreSum2  operator/  ( const double b ) const 
      { LegendreSum2 c { *this } ; c /= b ; return c ; }
      // ======================================================================
      LegendreSum2& operator+= ( const double value ) { m_pars[0] += value ; return *this ; }
      LegendreSum2& operator-= ( const double value ) { m_pars[0] -= value ; return *this ; }
      LegendreSum2& operator*= ( const double value ) 
      { Ostap::Math::scale ( m_pars ,     value ) ; return *this ; }
      LegendreSum2& operator/= ( const double value ) 
      { Ostap::Math::scale ( m_pars , 1.0/value ) ; return *this ; }
      /// negate it! 
      LegendreSum2  operator-    () const ;
      // ======================================================================
    public: // please python
      // ======================================================================
      LegendreSum2& __iadd__     ( const double value ) { return  (*this) += value ; }
      LegendreSum2& __isub__     ( const double value ) { return  (*this) -= value ; }
      LegendreSum2& __imult__    ( const double value ) { return  (*this) *= value ; }
      LegendreSum2& __idiv__     ( const double value ) { return  (*this) /= value ; }
      LegendreSum2& __itruediv__ ( const double value ) { return  (*this) /= value ; }
      // ======================================================================
      LegendreSum2  __add__      ( const double value ) { return  (*this) +  value ; }
      LegendreSum2  __sub__      ( const double value ) { return  (*this) -  value ; }
      LegendreSum2  __mult__     ( const double value ) { return  (*this) *  value ; }
      LegendreSum2  __div__      ( const double value ) { return  (*this) /  value ; }
      LegendreSum2  __truediv__  ( const double value ) { return  (*this) /  value ; }
      // ======================================================================
      LegendreSum2  __radd__     ( const double value ) { return  (*this) +  value ; }
      LegendreSum2  __rsub__     ( const double value ) { return -(*this) +  value ; }
      LegendreSum2  __rmult__    ( const double value ) { return  (*this) *  value ; }
      // ======================================================================
      /// negate it! 
      LegendreSum2  __neg__      () const { return  -(*this) ; }
      // ======================================================================
    public : // projections/integrations  as functions 
      // ======================================================================
      /** integrate over x dimension 
       *  \f$ f(y) =  \int_{x_{min}}^{x_{max}} F(x,y) {\mathrm{d}} x \f$
       */
      LegendreSum integralX  () const ;
      // ======================================================================
      /** integrate over Y dimension 
       *  \f$ f(x) =  \int_{y_{min}}^{y_{max}} F(x,y) {\mathrm{d}} y \f$
       */
      LegendreSum integralY  () const ;
      // ======================================================================
    public:
      // ======================================================================
      /** integrate over x dimension 
       *  \f$ f(y) =  \int_{x_{low}}^{x_{high}} F(x,y) {\mathrm{d}} x \f$
       */
      LegendreSum integralX  ( const double xlow  , const double xhigh ) const ;
      // ======================================================================
      /** integrate over y dimension 
       *  \f$ f(x) =  \int_{y_{low}}^{y_{high}} F(x,y) {\mathrm{d}} y \f$
       */
      LegendreSum integralY  ( const double ylow  , const double yhigh ) const ;
      // ======================================================================
    public:
      // =======================================================================
      /** get the integral 
       *  \f$ \int_{x_{low}}^{x_{high}}
       *      \int_{y_{low}}^{y_{high}} F(x,y){\mathrm{d}} x {\mathrm{d}} y \f$
       */
      double      integral  
      ( const double xlow  ,   
        const double xhigh , 
        const double ylow  ,   
        const double yhigh ) const ;
      // =====================================================================
      /** get the integral 
       *  \f$ \int_{x_{min}}^{x_{max}}
       *      \int_{y_{min}}^{y_{max}}F(x,y){\mathrm{d}} x {\mathrm{d}} y \f$
       */
      double      integral   ()                     const ;
      // ======================================================================
    public:
      // ======================================================================
      /// "transpose": \f$ S(x,y) \equiv F(y,x) \f$
      LegendreSum2 transpose() const ;
      /// "transpose" : \f$ S(x,y) \equiv F(y,x) \f$
      inline LegendreSum2 T() const { return transpose () ; }
      // ======================================================================
    private: 
      // ======================================================================
      /// 2D-index into 1D-index 
      inline std::size_t index 
      ( const unsigned short ix , 
        const unsigned short iy ) const 
      { return ix * ( m_NY + 1 ) + iy ; }
      // ======================================================================
      double calculate () const
      {
        long double value = 0 ;
        for ( unsigned int ix = 0 ; ix <= m_NX ; ++ix ) 
        { for ( unsigned int iy = 0 ; iy <= m_NY ; ++iy ) 
          { value += m_pars [ index ( ix , iy ) ] * 1.0L
              * m_cache_x [ ix ] * m_cache_y[ iy ] ; } }
        return value ;
      }
      // ======================================================================
    private : 
      // ======================================================================
      /// x-degree 
      std::size_t         m_NX   ; // x-degree
      /// y-degree 
      std::size_t         m_NY   ; //  y-degree
      // ======================================================================
      /// x-min 
      double              m_xmin ; // x-min
      /// x-max 
      double              m_xmax ; // x-max
      // ======================================================================      
      /// y-min 
      double              m_ymin ; // y-min
      /// y-max 
      double              m_ymax ; // y-max
      // ======================================================================
    private: // cache 
      // ======================================================================
      mutable std::vector<double> m_cache_x ;
      mutable std::vector<double> m_cache_y ;      
      // ======================================================================
    };
    // ========================================================================
    inline LegendreSum2 operator+( const double b , const LegendreSum2& a ) { return  a + b ; }    
    inline LegendreSum2 operator-( const double b , const LegendreSum2& a ) { return -a + b ; }
    inline LegendreSum2 operator*( const double b , const LegendreSum2& a ) { return  a * b ; }
    // ========================================================================
    /// Decartes product of two Legendre sums 
    inline LegendreSum2 operator*
    ( const LegendreSum& a , 
      const LegendreSum& b ) { return  LegendreSum2 ( a , b ) ; }
    // ========================================================================
    /** @class  LegendreSum3  Ostap/Parameterization.h
     *  3D-parameterization  as sum of Legendre polynomials 
     *  \f[ S(x,y,z) = \sum_{i0}^{n_x}\sum_{j=0}^{n_y}\sum_{k=0}^{n_z}c_{i,j}
     *      P_i(x^\prime)
     *      P_j(y^\prime)
     *      P_k(z^\prime)\f],
     * where:
     *  - \f$ x^{\prime} = \frac{2x-x_{min}-x_{max}}{x_{max}-x_{min}}\f$,
     *  - \f$ y^{\prime} = \frac{2y-y_{min}-y_{max}}{y_{max}-y_{min}}\f$. 
     *  - \f$ z^{\prime} = \frac{2z-z_{min}-z_{max}}{z_{max}-z_{min}}\f$. 
     *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
     *  @date 2019-06-30
     */
    class LegendreSum3 : public Parameters 
    {
    public:
      // ======================================================================
      /** constructor 
       *  \f[ S(x,y,z) = \sum_{i0}^{n_x}\sum_{j=0}^{n_y}\sum_{k=0}^{n_z}c_{i,j}
       *      P_i(x^\prime)
       *      P_j(y^\prime)
       *      P_k(z^\prime)\f],
       * where:
       *  - \f$ x^{\prime} = \frac{2x-x_{min}-x_{max}}{x_{max}-x_{min}}\f$,
       *  - \f$ y^{\prime} = \frac{2y-y_{min}-y_{max}}{y_{max}-y_{min}}\f$. 
       *  - \f$ z^{\prime} = \frac{2z-z_{min}-z_{max}}{z_{max}-z_{min}}\f$. 
       */
      LegendreSum3 
      ( const unsigned short NX =  0 , 
        const unsigned short NY =  0 , 
        const unsigned short NZ =  0 , 
        const double   xmin     = -1 , 
        const double   xmax     =  1 , 
        const double   ymin     = -1 , 
        const double   ymax     =  1 ,
        const double   zmin     = -1 , 
        const double   zmax     =  1 ) ;
      // ======================================================================
      /** constructor orm the product of two Legendre sums
       *  \f$ S(x,y,z) = S_x(x)\times S_y(y) \times S_z(z) \f$ 
       *  @param sx (INPUT) the first  Legendre sum 
       *  @param sy (INPUT) the second Legendre sum 
       *  @param sz (INPUT) the third Legendre sum 
       */
      LegendreSum3 
      ( const LegendreSum&  sx , 
        const LegendreSum&  sy ,
        const LegendreSum&  sz ) ;
      // ======================================================================
      /** constructor form the product of two Legendre sums
       *  \f$ S(x,y,z) = S_{xy}(x,y) \times S_z(z) \f$ 
       *  @param sxy (INPUT) the first  Legendre sum 
       *  @param sz  (INPUT) the second Legendre sum 
       */
      LegendreSum3 
      ( const LegendreSum2& sxy , 
        const LegendreSum&  sz  ) ;
      // ======================================================================
      /** constructor form the product of two Legendre sums
       *  \f$ S(x,y,z) = S_x(x) \times S_{yz}(y,z)\f$ 
       *  @param sx   (INPUT) the first  Legendre sum 
       *  @param syz  (INPUT) the second Legendre sum 
       */
      LegendreSum3 
      ( const LegendreSum&  sx   , 
        const LegendreSum2& syz  ) ;
      // ======================================================================
    public: 
      // ======================================================================
      /// get the value
      double evaluate 
      ( const double x , 
        const double y , 
        const double z ) const ;
      // ======================================================================
      /// get the value
      double operator () 
      ( const double x , 
        const double y , 
        const double z ) const 
      { return 
          x < m_xmin || x > m_xmax ? 0 : 
          y < m_ymin || y > m_ymax ? 0 : 
          z < m_zmin || z > m_zmax ? 0 : evaluate ( x , y , z ) ; }
      // =====================================================================
    public:
      // ======================================================================
      /// get lower edge
      double xmin  () const { return m_xmin ; }
      /// get upper edge
      double xmax  () const { return m_xmax ; }
      // ======================================================================
      /// get lower edge
      double ymin  () const { return m_ymin ; }
      /// get upper edge
      double ymax  () const { return m_ymax ; }
      // ======================================================================
      /// get lower edge
      double zmin  () const { return m_zmin ; }
      /// get upper edge
      double zmax  () const { return m_zmax ; }
      // ======================================================================
    public:
      // ======================================================================
      std::size_t degree_x  () const { return m_NX ; }
      std::size_t degree_y  () const { return m_NY ; }
      std::size_t degree_z  () const { return m_NZ ; }
      std::size_t nx        () const { return m_NX ; }
      std::size_t ny        () const { return m_NY ; }
      std::size_t nz        () const { return m_NZ ; }
      // ======================================================================
    public:
      // ======================================================================
      double x  ( const double tx ) const
      { return  0.5 * ( tx * ( m_xmax - m_xmin ) +   m_xmax + m_xmin ) ; }
      double tx ( const double  x ) const
      { return (  2 *    x   - m_xmax - m_xmin ) / ( m_xmax - m_xmin ) ; }
      // ======================================================================
      double y  ( const double ty ) const
      { return  0.5 * ( ty * ( m_ymax - m_ymin ) +   m_ymax + m_ymin ) ; }
      double ty ( const double  y ) const
      { return (  2 *    y   - m_ymax - m_ymin ) / ( m_ymax - m_ymin ) ; }
      // ======================================================================
      double z  ( const double tz ) const
      { return  0.5 * ( tz * ( m_zmax - m_zmin ) +   m_zmax + m_zmin ) ; }
      double tz ( const double  z ) const
      { return (  2 *    z   - m_zmax - m_zmin ) / ( m_zmax - m_zmin ) ; }
      // ======================================================================
    public:
      // ======================================================================
      using Parameters::par    ;
      using Parameters::setPar ;
      /// get 3D-parameter 
      inline double par 
      ( const unsigned short ix , 
        const unsigned short iy ,
        const unsigned short iz ) const 
      { return par ( index ( ix , iy , iz ) ) ;  }    
      // ======================================================================
      /// set 3D-parameter
      bool setPar
      ( const unsigned short ix     ,
        const unsigned short iy     ,
        const unsigned short iz     ,
        const double       value  ) 
      { return setPar ( index ( ix , iy , iz ) , value ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /** update  the Legendre expansion by addition of one "event" with 
       *  the given weight
       *  @code
       *  LegendreSum3 sum = ... ;
       *  for ( auto x : .... ) { sum.fill ( x , y , z ) ; }
       *  @endcode
       *  This is a useful function to make an unbinned parameterization 
       *  of certain distribution and/or efficiency 
       */
      bool fill 
      ( const double x          , 
        const double y          , 
        const double z          , 
        const double weight = 1 ) ;
      // ======================================================================
    public: // several useful operators and operations 
      // ======================================================================
      LegendreSum3  operator+  ( const double b ) const 
      { LegendreSum3 c { *this } ; c += b ; return c ; }
      LegendreSum3  operator-  ( const double b ) const 
      { LegendreSum3 c { *this } ; c -= b ; return c ; }
      LegendreSum3  operator*  ( const double b ) const 
      { LegendreSum3 c { *this } ; c *= b ; return c ; }
      LegendreSum3  operator/  ( const double b ) const 
      { LegendreSum3 c { *this } ; c /= b ; return c ; }
      // ======================================================================
      LegendreSum3& operator+= ( const double value ) { m_pars[0] += value ; return *this ; }
      LegendreSum3& operator-= ( const double value ) { m_pars[0] -= value ; return *this ; }
      LegendreSum3& operator*= ( const double value ) 
      { Ostap::Math::scale ( m_pars ,     value ) ; return *this ; }
      LegendreSum3& operator/= ( const double value ) 
      { Ostap::Math::scale ( m_pars , 1.0/value ) ; return *this ; }
      /// negate it! 
      LegendreSum3  operator-    () const ;
      // ======================================================================
    public: // please python
      // ======================================================================
      LegendreSum3& __iadd__     ( const double value ) { return  (*this) += value ; }
      LegendreSum3& __isub__     ( const double value ) { return  (*this) -= value ; }
      LegendreSum3& __imult__    ( const double value ) { return  (*this) *= value ; }
      LegendreSum3& __idiv__     ( const double value ) { return  (*this) /= value ; }
      LegendreSum3& __itruediv__ ( const double value ) { return  (*this) /= value ; }
      // ======================================================================
      LegendreSum3  __add__      ( const double value ) { return  (*this) +  value ; }
      LegendreSum3  __sub__      ( const double value ) { return  (*this) -  value ; }
      LegendreSum3  __mult__     ( const double value ) { return  (*this) *  value ; }
      LegendreSum3  __div__      ( const double value ) { return  (*this) /  value ; }
      LegendreSum3  __truediv__  ( const double value ) { return  (*this) /  value ; }
      // ======================================================================
      LegendreSum3  __radd__     ( const double value ) { return  (*this) +  value ; }
      LegendreSum3  __rsub__     ( const double value ) { return -(*this) +  value ; }
      LegendreSum3  __rmult__    ( const double value ) { return  (*this) *  value ; }
      // ======================================================================
      /// negate it! 
      LegendreSum3  __neg__      () const { return  -(*this) ; }
      // ======================================================================
    public : // projections/integrations
      // ======================================================================
      /** integrate over x dimension 
       *  \f$ f_x(y,z) =  \int_{x_{min}}^{x_{max}} f(x,y,z) {\mathrm{d}} x \f$
       */
      LegendreSum2  integralX  () const ;
      // ======================================================================
      /** integrate over y dimension 
       *  \f$ f_y(x,z) =  \int_{y_{min}}^{y_{max}} f(x,y,z) {\mathrm{d}} y \f$
       */
      LegendreSum2 integralY  () const ;
      // ======================================================================
      /** integrate over z dimension 
       *  \f$ f_z(x,y) =  \int_{z_{min}}^{z_{max}} f(x,y,z) {\mathrm{d}} z \f$
       */
      LegendreSum2 integralZ  () const ;
      // ======================================================================
    public:
      // ======================================================================
      /** integrate over x dimension 
       *  \f$ f_x(y,z) =  \int_{x_{low}}^{x_{high}} f(x,y,z) {\mathrm{d}} x \f$
       */
      LegendreSum2 integralX ( const double xlow , const double xhigh ) const ;
      // ======================================================================
      /** integrate over y dimension 
       *  \f$ f_y(x,z) =  \int_{y_{low}}^{y_{high}} f(x,y,z) {\mathrm{d}} y \f$
       */
      LegendreSum2 integralY ( const double ylow , const double yhigh ) const ;
      // ======================================================================
      /** integrate over z dimension 
       *  \f$ f_z(x,y) =  \int_{z_{low}}^{z_{high}} f(x,y,z) {\mathrm{d}} z \f$
       */
      LegendreSum2 integralZ ( const double zlow , const double zhigh ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /** integral 
       *  \f$ \int_{x_{low}}^{x_{high}}\int_{y_{low}}^{y_{high}}
       *    \int_{z_{low}}^{z_{high}} f(x,y,z) {\mathrm{d}} x {\mathrm{d}} y {\mathrm{d}} z \f$ 
       */
      double integral 
      ( const double xlow , const double xhigh , 
        const double ylow , const double yhigh , 
        const double zlow , const double zhigh ) const ;
      // ======================================================================
      /** integral 
       *  \f$ \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}}
       *    \int_{z_{min}}^{z_{max}} f(x,y,z) {\mathrm{d}} x {\mathrm{d}} y {\mathrm{d}} z \f$ 
       */
      double integral () const ;
      // ======================================================================
    private: 
      // ======================================================================
      /// 3D-index into 1D-index 
      inline std::size_t index 
      ( const unsigned short ix , 
        const unsigned short iy , 
        const unsigned short iz ) const 
      { return ( ix * ( m_NY + 1 ) + iy ) * ( m_NZ + 1 ) + iz ; }
      // ====================================================================== 
      double calculate () const
      {
        long double value = 0 ;
        for ( unsigned short ix = 0 ; ix <= m_NX ; ++ix ) 
        { for ( unsigned short iy = 0 ; iy <= m_NY ; ++iy ) 
          { for ( unsigned short iz = 0 ; iz <= m_NZ ; ++iz ) 
            { value += m_pars [ index ( ix , iy , iz ) ] * 1.0L 
                * m_cache_x[ ix ] * m_cache_y[ iy ] * m_cache_z[ iz ] ; } } }
        return value ;
      }
      // ======================================================================
    private:
      // ======================================================================
      /// x-degree 
      std::size_t         m_NX   ; // x-degree
      /// y-degree 
      std::size_t         m_NY   ; //  y-degree
      /// z-degree 
      std::size_t         m_NZ   ; //  z-degree
      // ======================================================================
      /// x-min 
      double              m_xmin ; // x-min
      /// x-max 
      double              m_xmax ; // x-max
      // ======================================================================      
      /// y-min 
      double              m_ymin ; // y-min
      /// y-max 
      double              m_ymax ; // y-max
      // ======================================================================
      /// z-min 
      double              m_zmin ; // z-min
      /// z-max 
      double              m_zmax ; // z-max
      // ======================================================================
    private: // cache 
      // ======================================================================
      mutable std::vector<double> m_cache_x ;
      mutable std::vector<double> m_cache_y ;      
      mutable std::vector<double> m_cache_z ;      
      // ======================================================================
    } ;
    // ========================================================================
    inline LegendreSum3 operator+( const double b , const LegendreSum3& a ) { return  a + b ; }    
    inline LegendreSum3 operator-( const double b , const LegendreSum3& a ) { return -a + b ; }
    inline LegendreSum3 operator*( const double b , const LegendreSum3& a ) { return  a * b ; }
    /// Decartes product of two Legendre sums 
    inline LegendreSum3 operator*
    ( const LegendreSum2& a , 
      const LegendreSum&  b ) { return  LegendreSum3 ( a , b ) ; }
    /// Decartes product of two Legendre sums 
    inline LegendreSum3 operator*
    ( const LegendreSum&  a , 
      const LegendreSum2& b ) { return  LegendreSum3 ( a , b ) ; }
    // ========================================================================
    /** @class  LegendreSum4  Ostap/Parameterization.h
     *  4D-parameterization  as sum of Legendre polynomials 
     *  \f$ S(x,y,z,u) = \sum_{i0}^{n_x}\sum_{j=0}^{n_y}
     *   \sum_{k=0}^{n_z}c_{i,j}\sum_{m=0}^{n_u}c_{i,j,k,m}
     *      P_i(x^\prime)
     *      P_j(y^\prime)
     *      P_k(z^\prime)
     *      P_m(u^\prime)\f$,
     * where:
     *  - \f$ x^{\prime} = \frac{2x-x_{min}-x_{max}}{x_{max}-x_{min}}\f$,
     *  - \f$ y^{\prime} = \frac{2y-y_{min}-y_{max}}{y_{max}-y_{min}}\f$. 
     *  - \f$ z^{\prime} = \frac{2z-z_{min}-z_{max}}{z_{max}-z_{min}}\f$. 
     *  - \f$ u^{\prime} = \frac{2u-u_{min}-u_{max}}{u_{max}-u_{min}}\f$. 
     *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
     *  @date 2019-06-30
     */
    class LegendreSum4 : public Parameters 
    {
    public:
      // ======================================================================
      /** constructor 
       *  \f$ S(x,y,z,u) = \sum_{i0}^{n_x}\sum_{j=0}^{n_y}
       *             \sum_{k=0}^{n_z}\sum_{m=0}^{n_u}c_{i,j,k,m}
       *      P_i(x^\prime)
       *      P_j(y^\prime)
       *      P_k(z^\prime)
       *      P_k(u^\prime)\f$,
       * where:
       *  - \f$ x^{\prime} = \frac{2x-x_{min}-x_{max}}{x_{max}-x_{min}}\f$,
       *  - \f$ y^{\prime} = \frac{2y-y_{min}-y_{max}}{y_{max}-y_{min}}\f$, 
       *  - \f$ z^{\prime} = \frac{2z-z_{min}-z_{max}}{z_{max}-z_{min}}\f$, 
       *  - \f$ u^{\prime} = \frac{2u-u_{min}-u_{max}}{u_{max}-u_{min}}\f$. 
       */
      LegendreSum4 
      ( const unsigned short NX =  0 , 
        const unsigned short NY =  0 , 
        const unsigned short NZ =  0 , 
        const unsigned short NU =  0 , 
        const double   xmin     = -1 , 
        const double   xmax     =  1 , 
        const double   ymin     = -1 , 
        const double   ymax     =  1 ,
        const double   zmin     = -1 , 
        const double   zmax     =  1 ,
        const double   umin     = -1 , 
        const double   umax     =  1 ) ;
      // ======================================================================
    public: 
      // ======================================================================
      /** constructor orm the product of two Legendre sums
       *  \f$ S(x,y,z) = S_x(x)\times S_y(y) \times S_z(z) \f$ 
       *  @param sx (INPUT) the first  Legendre sum 
       *  @param sy (INPUT) the second Legendre sum 
       *  @param sz (INPUT) the third  Legendre sum 
       *  @param su (INPUT) the fourth Legendre sum 
       */
      LegendreSum4 
      ( const LegendreSum&  sx , 
        const LegendreSum&  sy ,
        const LegendreSum&  sz ,
        const LegendreSum&  su ) ;
      // ======================================================================
    public: 
      // ======================================================================
      /// get the value
      double evaluate
      ( const double x , 
        const double y , 
        const double z ,
        const double u ) const ;
      // ======================================================================
      /// get the value
      double operator ()
      ( const double x , 
        const double y , 
        const double z ,
        const double u ) const 
      { return 
          x < m_xmin || x > m_xmax ? 0 : 
          y < m_ymin || y > m_ymax ? 0 : 
          z < m_zmin || z > m_zmax ? 0 : 
          u < m_umin || u > m_umax ? 0 : evaluate ( x , y , z , u ) ; }
      // =====================================================================
    public:
      // ======================================================================
      /// get lower edge
      double xmin  () const { return m_xmin ; }
      /// get upper edge
      double xmax  () const { return m_xmax ; }
      // ======================================================================
      /// get lower edge
      double ymin  () const { return m_ymin ; }
      /// get upper edge
      double ymax  () const { return m_ymax ; }
      // ======================================================================
      /// get lower edge
      double zmin  () const { return m_zmin ; }
      /// get upper edge
      double zmax  () const { return m_zmax ; }
      // ======================================================================
      /// get lower edge
      double umin  () const { return m_umin ; }
      /// get upper edge
      double umax  () const { return m_umax ; }
      // ======================================================================
    public:
      // ======================================================================
      std::size_t degree_x  () const { return m_NX ; }
      std::size_t degree_y  () const { return m_NY ; }
      std::size_t degree_z  () const { return m_NZ ; }
      std::size_t degree_u  () const { return m_NU ; }
      std::size_t nx        () const { return m_NX ; }
      std::size_t ny        () const { return m_NY ; }
      std::size_t nz        () const { return m_NZ ; }
      std::size_t nu        () const { return m_NU ; }
      // ======================================================================
    public:
      // ======================================================================
      double x  ( const double tx ) const
      { return  0.5 * ( tx * ( m_xmax - m_xmin ) +   m_xmax + m_xmin ) ; }
      double tx ( const double  x ) const
      { return (  2 *    x   - m_xmax - m_xmin ) / ( m_xmax - m_xmin ) ; }
      // ======================================================================
      double y  ( const double ty ) const
      { return  0.5 * ( ty * ( m_ymax - m_ymin ) +   m_ymax + m_ymin ) ; }
      double ty ( const double  y ) const
      { return (  2 *    y   - m_ymax - m_ymin ) / ( m_ymax - m_ymin ) ; }
      // ======================================================================
      double z  ( const double tz ) const
      { return  0.5 * ( tz * ( m_zmax - m_zmin ) +   m_zmax + m_zmin ) ; }
      double tz ( const double  z ) const
      { return (  2 *    z   - m_zmax - m_zmin ) / ( m_zmax - m_zmin ) ; }
      // ======================================================================
      double u  ( const double tu ) const
      { return  0.5 * ( tu * ( m_umax - m_umin ) +   m_umax + m_umin ) ; }
      double tu ( const double  u ) const
      { return (  2 *    u   - m_umax - m_umin ) / ( m_umax - m_umin ) ; }
      // ======================================================================
    public:
      // ======================================================================
      using Parameters::par    ;
      using Parameters::setPar ;
      /// get 4D-parameter 
      inline double par
      ( const unsigned short ix , 
        const unsigned short iy ,
        const unsigned short iz ,
        const unsigned short iu ) const 
      { return par ( index ( ix , iy , iz , iu ) ) ;  }    
      // ======================================================================
      /// set 4D-parameter
      bool setPar
      ( const unsigned short ix     ,
        const unsigned short iy     ,
        const unsigned short iz     ,
        const unsigned short iu     ,
        const double       value  ) 
      { return  setPar ( index ( ix , iy , iz , iu ) , value ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /** update  the Legendre expansion by addition of one "event" with 
       *  the given weight
       *  @code
       *  LegendreSum3 sum = ... ;
       *  for ( auto x : .... ) { sum.fill ( x , y , z , u ) ; }
       *  @endcode
       *  This is a useful function to make an unbinned parameterization 
       *  of certain distribution and/or efficiency 
       */
      bool fill
      ( const double x          , 
        const double y          , 
        const double z          , 
        const double u          , 
        const double weight = 1 ) ;
      // ======================================================================
    public: // several useful operators and operations 
      // ======================================================================
      LegendreSum4  operator+  ( const double b ) const 
      { LegendreSum4 c { *this } ; c += b ; return c ; }
      LegendreSum4  operator-  ( const double b ) const 
      { LegendreSum4 c { *this } ; c -= b ; return c ; }
      LegendreSum4  operator*  ( const double b ) const 
      { LegendreSum4 c { *this } ; c *= b ; return c ; }
      LegendreSum4  operator/  ( const double b ) const 
      { LegendreSum4 c { *this } ; c /= b ; return c ; }
      // ======================================================================
      LegendreSum4& operator+= ( const double value ) { m_pars[0] += value ; return *this ; }
      LegendreSum4& operator-= ( const double value ) { m_pars[0] -= value ; return *this ; }
      LegendreSum4& operator*= ( const double value ) 
      { Ostap::Math::scale ( m_pars ,     value ) ; return *this ; }
      LegendreSum4& operator/= ( const double value ) 
      { Ostap::Math::scale ( m_pars , 1.0/value ) ; return *this ; }
      /// negate it! 
      LegendreSum4  operator-    () const ;
      // ======================================================================
    public: // please python
      // ======================================================================
      LegendreSum4& __iadd__     ( const double value ) { return  (*this) += value ; }
      LegendreSum4& __isub__     ( const double value ) { return  (*this) -= value ; }
      LegendreSum4& __imult__    ( const double value ) { return  (*this) *= value ; }
      LegendreSum4& __idiv__     ( const double value ) { return  (*this) /= value ; }
      LegendreSum4& __itruediv__ ( const double value ) { return  (*this) /= value ; }
      // ======================================================================
      LegendreSum4  __add__      ( const double value ) { return  (*this) +  value ; }
      LegendreSum4  __sub__      ( const double value ) { return  (*this) -  value ; }
      LegendreSum4  __mult__     ( const double value ) { return  (*this) *  value ; }
      LegendreSum4  __div__      ( const double value ) { return  (*this) /  value ; }
      LegendreSum4  __truediv__  ( const double value ) { return  (*this) /  value ; }
      // ======================================================================
      LegendreSum4  __radd__     ( const double value ) { return  (*this) +  value ; }
      LegendreSum4  __rsub__     ( const double value ) { return -(*this) +  value ; }
      LegendreSum4  __rmult__    ( const double value ) { return  (*this) *  value ; }
      // ======================================================================
      /// negate it! 
      LegendreSum4  __neg__      () const { return  -(*this) ; }
      // ======================================================================
    public : // projections/integrations
      // ======================================================================
      /** integrate over x dimension 
       *  \f$ f(y,z,u) =  \int_{x_{min}}^{x_{max}} F(x,y,z,u) {\mathrm{d}} x \f$
       */
      LegendreSum3  integralX  () const ;
      // ======================================================================
      /** integrate over y dimension 
       *  \f$ f(x,z,u) =  \int_{y_{min}}^{y_{max}} F(x,y,z,u) {\mathrm{d}} y \f$
       */
      LegendreSum3  integralY  () const ;
      // ======================================================================
      /** integrate over z dimension 
       *  \f$ f(x,y,u) =  \int_{z_{min}}^{z_{max}} F(x,y,z,u) {\mathrm{d}} z \f$
       */
      LegendreSum3  integralZ  () const ;
      // ======================================================================
      /** integrate over u dimension 
       *  \f$ f(x,y,z) =  \int_{u_{min}}^{u_{max}} F(x,y,z,u) {\mathrm{d}} u \f$
       */
      LegendreSum3  integralU  () const ;
      // ======================================================================
    public:
      // ======================================================================
      /** integrate over x dimension 
       *  \f$ f(y,z,u) =  \int_{x_{low}}^{x_{high}} F(x,y,z,u) {\mathrm{d}} x \f$
       */
      LegendreSum3 integralX ( const double xlow , const double xhigh ) const ;
      // ======================================================================
      /** integrate over y dimension 
       *  \f$ f(x,z,u) =  \int_{y_{low}}^{y_{high}} F(x,y,z,u) {\mathrm{d}} y \f$
       */
      LegendreSum3 integralY ( const double ylow , const double yhigh ) const ;
      // ======================================================================
      /** integrate over z dimension 
       *  \f$ f(x,y,u) =  \int_{z_{low}}^{z_{high}} F(x,y,z,u) {\mathrm{d}} z \f$
       */
      LegendreSum3 integralZ ( const double zlow , const double zhigh ) const ;
      // ======================================================================
      /** integrate over u dimension 
       *  \f$ f(x,y,z) =  \int_{u_{low}}^{u_{high}} F(x,y,z,u) {\mathrm{d}} u \f$
       */
      LegendreSum3 integralU ( const double ulow , const double uhigh ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /** integral 
       *  \f$ \int_{x_{low}}^{x_{high}}\int_{y_{low}}^{y_{high}}
       *      \int_{z_{low}}^{z_{high}}\int_{u_{low}}^{u_{high}}
       *       f(x,y,z,u) {\mathrm{d}} x {\mathrm{d}} y {\mathrm{d}} z {\mathrm{d}} u \f$ 
       */
      double integral
      ( const double xlow , const double xhigh , 
        const double ylow , const double yhigh , 
        const double zlow , const double zhigh ,
        const double ulow , const double uhigh) const ;
      // ======================================================================
      /** integral 
       *  \f$ \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}}
       *      \int_{z_{min}}^{z_{max}}\int_{u_{min}}^{u_{max}} 
       *      f(x,y,z,u) {\mathrm{d}} x {\mathrm{d}} y {\mathrm{d}} z {\mathrm{d}} u \f$ 
       */
      double integral () const ;
      // ======================================================================
    private: 
      // ======================================================================
      /// 4D-index into 1D-index 
      inline unsigned int index 
      ( const unsigned short ix , 
        const unsigned short iy , 
        const unsigned short iz ,
        const unsigned short iu ) const 
      { return ( ( ix * ( m_NY + 1 ) + iy ) * ( m_NZ + 1 ) + iz ) * ( m_NU + 1 ) + iu ; }
      // ====================================================================== 
      double calculate () const
      {
        long double value = 0 ;
        for ( unsigned short ix = 0 ; ix <= m_NX ; ++ix ) 
        { for ( unsigned short iy = 0 ; iy <= m_NY ; ++iy ) 
          { for ( unsigned short iz = 0 ; iz <= m_NZ ; ++iz ) 
            { for ( unsigned short iu = 0 ; iu <= m_NU ; ++iu ) 
              { value += m_pars [ index ( ix , iy , iz , iu ) ] * 1.0L
                  * m_cache_x[ ix ] * m_cache_y[ iy ] 
                  * m_cache_z[ iz ] * m_cache_u[ iu ] ; } } } }
        return value ;
      }
      // ======================================================================
    private:
      // ======================================================================
      /// x-degree 
      std::size_t         m_NX   ; // x-degree
      /// y-degree 
      std::size_t         m_NY   ; //  y-degree
      /// z-degree 
      std::size_t         m_NZ   ; //  z-degree
      /// u-degree 
      std::size_t         m_NU   ; //  u-degree
      // ======================================================================
      /// x-min 
      double              m_xmin ; // x-min
      /// x-max 
      double              m_xmax ; // x-max
      // ======================================================================      
      /// y-min 
      double              m_ymin ; // y-min
      /// y-max 
      double              m_ymax ; // y-max
      // ======================================================================
      /// z-min 
      double              m_zmin ; // z-min
      /// z-max 
      double              m_zmax ; // z-max
      // ======================================================================
      /// u-min 
      double              m_umin ; // u-min
      /// u-max 
      double              m_umax ; // u-max
      // ======================================================================
    private: // cache 
      // ======================================================================
      mutable std::vector<double> m_cache_x ;
      mutable std::vector<double> m_cache_y ;      
      mutable std::vector<double> m_cache_z ;      
      mutable std::vector<double> m_cache_u ;      
      // ======================================================================
    } ;
    // ========================================================================
    inline LegendreSum4 operator+( const double b , const LegendreSum4& a ) { return  a + b ; }    
    inline LegendreSum4 operator-( const double b , const LegendreSum4& a ) { return -a + b ; }
    inline LegendreSum4 operator*( const double b , const LegendreSum4& a ) { return  a * b ; }
    // /// Decartes product of two Legendre sums 
    // inline LegendreSum3 operator*
    // ( const LegendreSum2& a , 
    //   const LegendreSum&  b ) { return  LegendreSum3 ( a , b ) ; }
    // /// Decartes product of two Legendre sums 
    // inline LegendreSum3 operator*
    // ( const LegendreSum&  a , 
    //   const LegendreSum2& b ) { return  LegendreSum3 ( a , b ) ; }
    // ========================================================================

    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
}; //                                                The end of namespace Ostap
// ============================================================================
//                                                                     The END 
// ============================================================================
#endif // OSTAP_PARAMETERIZATION_H
// ============================================================================
