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
#include "Ostap/Utils.h"
// ============================================================================
/** @file Ostap/Polarization.h
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
      typedef std::array<Ostap::LorentzVector,4> Frame ;
      // =====================================================================
      /** @typedef PolVectors
       *  three polarization vectors (-1,0,1)
       */
      typedef std::array<Ostap::ComplexLorentzVector,3> PolVectors ;
      // =====================================================================
      /// helper structure to keep (cos_theta,phi) result
      struct Angles 
      {
        double cos_theta ;
        double phi       ;
      } ;  
      // =====================================================================
      /** @typedef Cosines
       *  The structire to keep the direction cosines
       */
      typedef std::array<double,3> Cosines ;
      // =====================================================================
      /** use Madison convention?
       *  @see "Polarization phenomena in nuclear reactions: 
       *       proceedings of 3rd international symposium on polarization 
       *       phenomena in nuclear reactions", 
       *       Eds. H.H. Barschall and W. Haeberli, 
       *       University of Wisconsin Press, Madison WI U.S.A. 1971. 
       */
      using UseMadisonConvention = Ostap::Utils::tagged_bool<struct UseMadisonConvention_tag>;
      // ===================================================================== 
    public:
      // =====================================================================
      /** get the three axes for polarization frame 
       *  @param f  the frame 
       *  @param P     4-momenta of the particle 
       *  @param beam1 4-momenta of the first colliding particle
       *  @param beam2 4-momenta of the second colliding particle
       *  @param madison use Madison convention?
       *  @return three polarization axes: x,y&z
       */
      static Frame frame
      ( const Frames                f                                    , 
        const Ostap::LorentzVector& P                                    ,
        const Ostap::LorentzVector& beam1                                ,
        const Ostap::LorentzVector& beam2                                , 
        const UseMadisonConvention  madison = UseMadisonConvention{true} ) ;
      // =====================================================================
    public: // polarization vectors for the given frame 
      // =====================================================================
      /** get the polarization vectors from the frame 
       *  @param f the frame 
       *  @return polarization vectors (-1,0,+1)
       */
      static PolVectors vectors ( const Frame& f ) ;
      // =====================================================================
    public: // direction cosines 
      // =====================================================================
      /** get the direction cosines of the particle direction 
       *  in the specified reference frame 
       *  @param p the particle
       *  @param f the frame 
       *  @return  direction  cosines 
       */
      static Cosines cosines 
      ( const Ostap::LorentzVector& p , 
        const Frame&                f ) ;
      // ====================================================================
      /** get the direction cosines of the particle direction 
       *  in the rest frame of particle m,  and the beam-momenta p1& p2  
       *  @param p the particle
       *  @param f the frame 
       *  @param m the particle that defined the frame 
       *  @param p1 4-momenta of the first colliding particle
       *  @param p2 4-momenta of the second colliding particle
       *  @param madison use Madison convention?
       *  @return  \f$ (\cos \theta,\phi)\f$ structure
       */
      static Cosines cosines
      ( const Ostap::LorentzVector& p                                    , 
        const Frames                f                                    , 
        const Ostap::LorentzVector& m                                    , 
        const Ostap::LorentzVector& beam1                                , 
        const Ostap::LorentzVector& beam2                                ,
        const UseMadisonConvention  madison = UseMadisonConvention{true} ) ;
      // =====================================================================
      /** get the angles \f$(\cos\theta,\phi)\f$ for the particle 
       *  in the defined frame
       *  @param p the particle
       *  @param f the frame 
       *  @return  \f$(\cos\theta,\phi)\f$ structure
       */
      static Angles angles
      ( const Ostap::LorentzVector& p , 
        const Frame&                f ) ;
      // ====================================================================
      /** get the angles \f$(\cos \theta,\phi)\f$ for the particle 
       *  in the rest frame of particle m,  and the beam-momenta p1& p2  
       *  @param p the particle
       *  @param f the frame 
       *  @param m the particle that defined the frame 
       *  @param p1 4-momenta of the first colliding particle
       *  @param p2 4-momenta of the second colliding particle
       *  @param madison use Madison convention?
       *  @return  \f$ (\cos \theta,\phi)\f$ structure
       */
      static Angles angles
      ( const Ostap::LorentzVector& p                                    , 
        const Frames               f                                    , 
        const Ostap::LorentzVector& m                                    , 
        const Ostap::LorentzVector& beam1                                , 
        const Ostap::LorentzVector& beam2                                ,
        const UseMadisonConvention  madison = UseMadisonConvention{true} ) ;
      // =====================================================================      
    } ; //                                The end of Ostap::Math::Polarization
    // =======================================================================
  } //                                        The end of namespace Ostap::Math
  // =========================================================================
} //                                                The end of namespace Ostap 
// ===========================================================================
//                                                                     The END 
// ===========================================================================
#endif // OSTAP_POLARIZATION_H
// ===========================================================================
