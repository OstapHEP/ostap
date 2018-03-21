// ============================================================================
#ifndef OSTAP_POLARIZATION_H 
#define OSTAP_POLARIZATION_H 1
// ============================================================================
// Include files
// ============================================================================
// STD& STD:
// ============================================================================
#include <array>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Vector4DTypes.h"
// ============================================================================
/** @file 
 *  Collection of  functions to deal with polarization  axes
 *  @see  M.Beneke, M.Kramer, M.Vanttiner, Phys.Rev. D57 (1998) 4258
 *  @see https://doi.org/10.1103/PhysRevD.57.4258
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    /** @class Polarization
     *  helper class/namespace to hold the function to deal with polarization 
     *  axes and frames 
     *  @see  M.Beneke, M.Kramer, M.Vanttiner, Phys.Rev. D57 (1998) 4258
     *  @see https://doi.org/10.1103/PhysRevD.57.4258
     */
    class Polarization 
    {
      // =====================================================================
    public:
      // =====================================================================
      /** @enum Frames 
       *  the list of polarization frames
       */
      // =====================================================================
      enum class Frames : unsigned short 
      {
        Recoil            = 0 , //  <--- Helicity
        GottfriedJackson      ,
        Target                ,
        CollinsSoper          ,
      };  
      // =====================================================================
      /** @typedef Frame
       *  three axes for polarization frame 
       */
      typedef std::array<Ostap::LorentzVector,3> Frame ;
      // =====================================================================
      /** @typedef Angles 
       *  \$ (\cos \theta,\phi)\$ structure  
       */
      typedef std::array<double,2> Angles ;
      // =====================================================================
    public:
      // =====================================================================
      /** get three  axes for polarization frame 
       *  @param f  the frame 
       *  @param P     4-momenta of the particle 
       *  @param beam1 4-momenta of the first colliding particle
       *  @param beam2 4-momenta of the second colliding particle
       *  @return three polarization axes: x,y&z
       */
      static Frame frame ( Frames f , 
                           const Ostap::LorentzVector& P     ,
                           const Ostap::LorentzVector& beam1 ,
                           const Ostap::LorentzVector& beam2 ) ;        
      // =====================================================================
      /** get the angles \$ (\cos \theta,\phi)\$ for the particle 
       *  in the defined frame
       *  @param p the particle
       *  @param f the frame 
       *  @return  \$ (\cos \theta,\phi)\$ structure
       */
      static Angles angles ( const Ostap::LorentzVector& p , 
                      const Frame&                f ) ;
      // ====================================================================
      /** get the angles \$ (\cos \theta,\phi)\$ for the particle 
       *  in the rest frame of particle m,  and the beam-momenta p1& p2  
       *  @param p the particle
       *  @param f the frame 
       *  @param m the particle that defined the frame 
       *  @param p1 4-momenta of the first colliding particle
       *  @param p2 4-momenta of the second colliding particle
       *  @return  \$ (\cos \theta,\phi)\$ structure
       */
      static Angles angles ( const Ostap::LorentzVector& p     , 
                             Frames                      f     , 
                             const Ostap::LorentzVector& m     , 
                             const Ostap::LorentzVector& beam1 , 
                             const Ostap::LorentzVector& beam2 );
      // =====================================================================
    } ; //                                The end of Ostap::Math::Polarization
    // =======================================================================
    /** boost Lorentz vector into  rest-frame of another Lorentz vector 
     *  @param what   the vextro to be bosted 
     *  @param frame  the 4-vector of the frame 
     *  @return boosted vector 
     */
      Ostap::LorentzVector boost 
      ( const Ostap::LorentzVector& what  ,
        const Ostap::LorentzVector& frame ) ;  
    // =======================================================================
  } //                                        The end of namespace Ostap::Math
  // =========================================================================
} //                                                The end of namespace Ostap 
// ===========================================================================
//                                                                     The END 
// ===========================================================================
#endif // OSTAP_POLARIZATION_H
// ===========================================================================
