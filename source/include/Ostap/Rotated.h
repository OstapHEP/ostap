// ============================================================================
#ifndef OSTAP_ROTATED_H 
#define OSTAP_ROTATED_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <cmath>
#include <functional>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Workspace.h"
#include "Ostap/Integrator.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    template <class PEAK>
    inline double pivot ( const PEAK& f ) { return f.peak() ; }
    // ========================================================================
    /** @class RotatedProduct
     *  2D-model that represents "rotated" product of two distributions 
     *  \f[ f(x,y) =  F_1\left(x^{\prime}\right)F_2\left(y^{\prime}\right)\f],
     *  where 
     *  - \f$ x^{\prime} = \delta x \cos \phi + \delta y \sin \phi + p_x\f$
     *  - \f$ y^{\prime} = \delta y \cos \phi - \delta x \sin \phi + p_y\f$
     *  - \f$ \delta x   = x - p_x \f$ 
     *  - \f$ \delta y   = y - p_y \f$ 
     *  - \f$(p_x,p_y)\f$ is a pivot point
     */
    template <class SIGNAL1 , class SIGNAL2 = SIGNAL1>
    class RotatedProduct 
    {
    public:
      // ======================================================================
      /// the first signal 
      typedef SIGNAL1 Signal1 ; // the first signal 
      /// the second signal 
      typedef SIGNAL2 Signal2 ; // the second signal 
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from two signals and the rotation phase 
       *  @param  s1 the first  signal
       *  @param  s2 the second signal
       *  @param  phase the  phase 
       */
      RotatedProduct ( const SIGNAL1& s1        , 
                       const SIGNAL2& s2        , 
                       const double   phase = 0 ) 
        : m_signal1  ( s1    ) 
        , m_signal2  ( s2    ) 
        , m_phase    ( phase ) 
      {}
      /// default constructor
      RotatedProduct() = default ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double operator () ( const double x , const double y ) const 
      { return this->evaluate ( x , y ) ; }
      // ======================================================================
      double    evaluate ( const double x , const double y ) const 
      {
        const double px = pivot ( this->m_signal1 ) ;
        const double py = pivot ( this->m_signal2 ) ;
        const double dx = x - px ;
        const double dy = y - py ;
        const double sp = std::sin ( this->m_phase ) ;
        const double cp = std::cos ( this->m_phase ) ;
        const double xp = cp * dx + sp * dy + px ;
        const double yp = cp * dy - sp * dx + py ;
        //
        return this->m_signal1 ( xp ) * this->m_signal2 ( yp ) ;
      }
      // ======================================================================
    public:
      // ======================================================================
      const  Signal1& signal1 () const { return m_signal1 ; }
      const  Signal2& signal2 () const { return m_signal2 ; }
      // ======================================================================
      inline Signal1& signal1 ()       { return m_signal1 ; }
      inline Signal2& signal2 ()       { return m_signal2 ; }
      // ======================================================================
      /// get the phase 
      double      phase () const { return m_phase ; }
      /// set the phase 
      inline bool setPhase ( const double value )
      {
        if      ( value >  M_PI ) { return this->setPhase ( value - M_PI ) ; }
        else if ( value < -M_PI ) { return this->setPhase ( value + M_PI ) ; }
        if ( value == m_phase )   { return false ; }
        m_phase = value ;
        return true ;
      }
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
      double integral ( const double xlow , const double xhigh ,
                        const double ylow , const double yhigh ) const 
      {
        static const Ostap::Math::Integrator i {} ;
        Ostap::Math::Integrator::function2 f2 = std::cref ( *this ) ;
        return i.integrate_with_cache ( this->tag() , f2 , xlow , xhigh , ylow , yhigh , this->m_workspace ) ;
      }
      // ======================================================================
      /** integral over x-dimension
       *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param y     variable
       *  @param xlow  low  edge in y
       *  @param xhigh high edge in y
       */
      double integrateX ( const double y    ,
                          const double xlow , const double xhigh ) const 
      {
        static const Ostap::Math::Integrator i {} ;
        Ostap::Math::Integrator::function2 f2 = std::cref ( *this ) ;
        double r1 = i.integrateX_with_cache ( this->tag() , f2 , y , xlow , xhigh , this->m_workspace ) ;
        // double r2 = i.integrateX ( f2 , y , xlow , xhigh , this->m_workspace ) ;
        // std::cout << "IX 1/2: " <<  r1 << "/" << r2  << " tag =" << this->tag() << " y=" << y << " min/max=" <<  xlow << "/" << xhigh <<  std::endl ;
        // return i.integrateX_with_cache ( this->tag() , f2 , y , xlow , xhigh , this->m_workspace ) ;
        // return i.integrateX ( f2 , y , xlow , xhigh , this->m_workspace ) ;
        return r1 ;
      }
      // ======================================================================
      /** integral over y-dimension
       *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param x     variable
       *  @param ylow  low  edge in x
       *  @param yhigh high edge in x
       */
      double integrateY ( const double x    ,
                          const double ylow , const double yhigh ) const 
      {
        static const Ostap::Math::Integrator i {} ;
        Ostap::Math::Integrator::function2 f2 = std::cref ( *this ) ;
        double r1 = i.integrateY_with_cache ( this->tag() , f2 , x , ylow , yhigh , this->m_workspace ) ;
        // double r2 = i.integrateY ( f2 , x , ylow , yhigh , this->m_workspace ) ;
        // std::cout << "IY 1/2: " <<  r1 << "/" << r2  << " tag =" << this->tag() << " x=" << x << " min/max=" << ylow << "/" << yhigh <<  std::endl ;
        return r1 ;
        // return i.integrateY_with_cache ( this->tag() , f2 , x , ylow , yhigh , this->m_workspace ) ;
        // return i.integrateY ( f2 , x , ylow , yhigh , this->m_workspace ) ;
      }      
      // ======================================================================
      double test ( double x ,  double y ) const 
      {
        Ostap::Math::Integrator::function2 f2 = std::cref ( *this ) ;
        auto fx  = std::bind ( f2 , std::placeholders::_1 , y ) ;
        auto fy  = std::bind ( f2 , x , std::placeholders::_1 ) ; 
        //
        return f2 ( x , y ) ;
      }
      // ======================================================================
    public:
      // ======================================================================
      std::size_t tag ()  const 
      {
        std::size_t seed = std::hash<double>() ( this->m_phase ) ;
        //
        seed ^= this->m_signal1.tag() + 0x9e3779b9 + (seed<<6) + (seed>>2) ; 
        seed ^= this->m_signal2.tag() + 0x9e3779b9 + (seed<<6) + (seed>>2) ;
        //
        return seed ;
      }
      // ======================================================================
    private:
      // ======================================================================
      /// the first signal 
      Signal1 m_signal1 { }     ; // the first signal 
      /// the second signal 
      Signal2 m_signal2 { }     ; // the second signal 
      /// rotation phase 
      double  m_phase   { 0.0 } ; // rotation phase 
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace   ; // integration workspace
      // ======================================================================
    } ;
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                     The  END 
// ============================================================================
#endif // OSTAP_ROTATED_H
// ============================================================================
