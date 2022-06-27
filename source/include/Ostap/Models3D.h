// ============================================================================
#ifndef OSTAP_MODELS3D_H 
#define OSTAP_MODELS3D_H 1
// ============================================================================
// include files
// ============================================================================
// STD & STL
// ============================================================================
#include <cmath>
// ============================================================================
// ROOT 
// ============================================================================
#include  "Math/EulerAngles.h" 
#include  "Math/Rotation3D.h" 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Models.h"
#include "Ostap/Integrator.h"
// ============================================================================
/// forwards declarations 
namespace ROOT { namespace Math { class EulerAngles ; } }
// ============================================================================
/** @file Ostap/Models3D.h
 *  collection of 3D-models
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2022-06-23
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /** @class Gauss3D 
     *  Simple 3D rotated Gaussian function 
     *  @date 2022-06-22
     */
    class Gauss3D 
    {
      // ======================================================================
    public: 
      // ======================================================================
      /** constructor 
       *  @param muX    x-location 
       *  @param muY    y-location 
       *  @param muZ    z-location 
       *  @param sigmaX x-sigma
       *  @param sigmaY y-sigma
       *  @param sigmaZ z-sigma
       *  @param phi   Euler angle phi 
       *  @param theta Euler angle theta 
       *  @param psi   Euler angle psi 
       */
      Gauss3D 
      ( const double muX    = 0 , 
        const double muY    = 0 , 
        const double muZ    = 0 , 
        const double sigmaX = 1 , 
        const double sigmaY = 1 ,
        const double sigmaZ = 1 ,
        const double phi    = 0 ,
        const double theta  = 0 ,
        const double psi    = 0 ) ;
      // =====================================================================
    public:
      // ======================================================================
      /// get the value
      double operator () ( const double x , 
                           const double y , 
                           const double z ) const ;
      // ======================================================================
    public:  // getters 
      // ======================================================================      
      // mux                                    
      double muX    () const { return m_muX    ; }
      // mux                                   
      double muY    () const { return m_muY    ; }
      // mux                                   
      double muZ    () const { return m_muZ    ; }
      // sigmax                                    
      double sigmaX () const { return m_sigmaX ; }
      // sigmay                                   
      double sigmaY () const { return m_sigmaY ; }
      // sigmay                                   
      double sigmaZ () const { return m_sigmaZ ; }
      // Euler phi 
      double phi    () const { return m_phi    ; }
      // Euler theta  
      double theta  () const { return m_theta  ; }
      // Euler psi   
      double psi    () const { return m_psi    ; }
      // ======================================================================
    public:  // setters 
      // ======================================================================      
      bool setMuX      ( const double value  ) ;
      bool setMuY      ( const double value  ) ;
      bool setMuZ      ( const double value  ) ;
      bool setSigmaX   ( const double value  ) ;
      bool setSigmaY   ( const double value  ) ;
      bool setSigmaZ   ( const double value  ) ;
      // ======================================================================
      bool setPhi      ( const double value  ) ;
      bool setTheta    ( const double value  ) ;
      bool setPsi      ( const double value  ) ;
      // ======================================================================      
      bool setEuler    ( const double phi    , 
                         const double theta  , 
                         const double psi    ) ;
      // ======================================================================      
      bool setEuler    ( const ROOT::Math::EulerAngles& angles ) ;
      bool setEuler    ( const ROOT::Math::Rotation3D&  angles ) ;
      // ======================================================================      
    public:
      // ======================================================================      
      /// get 3D rotation 
      const ROOT::Math::Rotation3D& rotation    () const { return m_rotation ; }
      /// get Euler angles 
      ROOT::Math::EulerAngles       eulerAngles () const 
      { return ROOT::Math::EulerAngles ( m_phi , m_theta , m_phi ) ; }
      // ======================================================================      
    public:
      // ======================================================================      
      /// get the integral over the whole 2D-region
      double integral () const ;
      // ======================================================================      
      /** get the integral over 2D-region
       *  \f[ \int_{x_{low}}^{x_{high}}
       *      \int_{y_{low}}^{y_{high}}
       *      \int_{z_{low}}^{z_{high}}
       *        \mathcal{B}(x,y,z) \mathrm{d}x\mathrm{d}y\mathrm{d}z \f]
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       *  @param zlow  low  edge in z
       *  @param zhigh high edge in z
       */
      double integral 
      ( const double xlow  ,
        const double xhigh ,
        const double ylow  , 
        const double yhigh ,
        const double zlow  , 
        const double zhigh ) const ;
      // ======================================================================
      /** integral over x&y-dimensions
       *  \f[ \int_{x_{low}}^{x_{high}}
       *      \int_{y_{low}}^{y_{high}} \mathcal{B}(x,y,z) \mathrm{d}x\mathrm{d}y\f]
       *  @param z     variable
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integrateXY 
      ( const double z    ,                          
        const double xlow , const double xhigh ,
        const double ylow , const double yhigh ) const ;
      // ======================================================================
      /** integral over x&z-dimensions
       *  \f[ \int_{x_{low}}^{x_{high}}
       *      \int_{z_{low}}^{z_{high}} \mathcal{B}(x,y,z) \mathrm{d}x\mathrm{d}z\f]
       *  @param y     variable
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       *  @param zlow  low  edge in y
       *  @param zhigh high edge in y
       */
      double integrateXZ 
      ( const double y    ,                          
        const double xlow , const double xhigh ,
        const double zlow , const double zhigh ) const ;      
      // ======================================================================
      /** integral over y&z-dimensions
       *  \f[ \int_{y_{low}}^{y_{high}}
       *      \int_{z_{low}}^{z_{high}} \mathcal{B}(x,y,z) \mathrm{d}y\mathrm{d}z\f]
       *  @param x     variable
       *  @param ylow  low  edge in x
       *  @param yhigh high edge in x
       *  @param zlow  low  edge in y
       *  @param zhigh high edge in y
       */
      double integrateYZ 
      ( const double x    ,                          
        const double ylow , const double yhigh ,
        const double zlow , const double zhigh ) const ;
      // ======================================================================
      /** integral over x-dimension
       *  \f[ \int_{x_{low}}^{x_{high}} \mathcal{B}(x,y,z) \mathrm{d}x\f]
       *  @param x     variable
       *  @param z     variable
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integrateX
      ( const double y    ,
        const double z    ,                          
        const double xlow , const double xhigh ) const ;
      // ======================================================================
      /** integral over y-dimension
       *  \f[ \int_{y_{low}}^{y_{high}} \mathcal{B}(x,y,z) \mathrm{d}y\f]
       *  @param y     variable
       *  @param z     variable
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       */
      double integrateY 
      ( const double x    ,
        const double z    ,
        const double ylow , const double yhigh ) const ;
      // ======================================================================
      /** integral over z-dimension
       *  \f[ \int_{z_{low}}^{z_{high}} \mathcal{B}(x,y,z) \mathrm{d}z\f]
       *  @param x     variable
       *  @param y     variable
       *  @param zlow  low  edge in z
       *  @param zhigh high edge in z
       */
      double integrateZ 
      ( const double x    ,
        const double y    ,
        const double zlow , const double zhigh ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private: 
      // ======================================================================
      /// mux 
      double m_muX        { 0.0 } ; // mux
      /// muy 
      double m_muY        { 0.0 } ; // muy 
      /// muz 
      double m_muZ        { 0.0 } ; // muz 
      /// sigmax 
      double m_sigmaX     { 1.0 } ; // sigmax 
      /// sigmay 
      double m_sigmaY     { 1.0 } ; // sigmay 
      /// sigmay 
      double m_sigmaZ     { 1.0 } ; // sigmaz 
      /// Euler phi 
      double m_phi        { 0.0 } ; // phi 
      /// Euler theta 
      double m_theta      { 0.0 } ; // theta 
      /// Euler psi
      double m_psi        { 0.0 } ; // phi 
      /// Rotation 
      ROOT::Math::Rotation3D m_rotation {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// workspace for numerical integration
      Ostap::Math::WorkSpace m_workspace {} ;
      // ======================================================================
    } ;
    // ========================================================================
  } //                                         The end of namespace Ostap::Math  
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_MODELS3D_H
// ============================================================================
