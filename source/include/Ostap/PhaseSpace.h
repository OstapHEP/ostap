// ============================================================================
#ifndef OSTAP_PHASESPACE_H 
#define OSTAP_PHASESPACE_H 1
// ============================================================================
// Include files
// ============================================================================
// STD &  STL
// ============================================================================
#include <complex>
#include <vector>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Workspace.h"
#include "Ostap/Dalitz.h"
// ============================================================================
/** @file Ostap/PhaseSpace.h
 *  collection of functions related to the phase space calculations 
 *  @author Vanya Belyaev
 */  
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /** @class PhaseSpace2
     *  Function to represent two-body phase space
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class  PhaseSpace2
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from two masses
      PhaseSpace2 ( const double m1 = 0 ,
                    const double m2 = 1 ) ;
      /// destructor
      ~PhaseSpace2 () ;                                         // deststructor
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate 2-body phase space
      double operator () ( const double x ) const ;
      /// integral
      double integral    ( const double xmin , const double xmax ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the momentum at center of mass
      double                q_  ( const double x ) const ;
      /// ditto but as complex
      std::complex<double>  q1_ ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      double m1        () const { return m_m1 ; }
      double m2        () const { return m_m2 ; }
      double lowEdge   () const { return m1 () + m2 () ; }
      double threshold () const { return m1 () + m2 () ; }
      // ======================================================================
    public:
      // ======================================================================
      ///  set the first mass 
      bool setM1  ( const double value ) ; // set the  first mass 
      ///  set the second mass 
      bool setM2  ( const double value ) ; // set the second  mass 
     // ======================================================================
    public:
      // ======================================================================
      // get the tag 
      std::size_t tag() const ;
      // ======================================================================
    private:
      // ======================================================================
      /// the first mass
      double m_m1 ; // the first  mass
      /// the second mass
      double m_m2 ; // the second mass
      // ======================================================================
    public :
      // ======================================================================
      /** calculate the particle momentum in rest frame
       *  \f[ q = \frac{1}{2}\frac{ \lambda^{\frac{1}{2}}
       *        \left( m^2 , m_1^2, m_2^2 \right) }{ m }\f],
       *  @param m the mass
       *  @param m1 the mass of the first particle
       *  @param m2 the mass of the second particle
       *  @return the momentum in rest frame (physical values only)
       */
      static double  q
      ( const double m  ,
        const double m1 ,
        const double m2 ) ;
      // ======================================================================
      /** calculate the particle momentum in rest frame
       *  - real for physical case 
       *  - imaginary for non-physical case (below the threshold)
       *  @param m the mass
       *  @param m1 the mass of the first particle
       *  @param m2 the mass of the second particle
       *  @return the momentum in rest frame  (imaginary for non-physical branch)
       */
      static std::complex<double> q1
      ( const double m  ,
        const double m1 ,
        const double m2 ) ;
      // ======================================================================
      /** calculate the phase space for   m -> m1 + m2
       *  \f[ \Phi = \frac{1}{8\pi} \left( \frac{ \lambda^{\frac{1}{2}}
       *       \left( m^2 , m_1^2, m_2^2 \right) }{ m^2 }\right)^{2L+1}\f],
       *  where \f$\lambda\f$ is a triangle function
       *  @param m the mass
       *  @param m1 the mass of the first particle
       *  @param m2 the mass of the second particle
       *  @param L  the orbital momentum 
       *  @return two-body phase space
       */
      static double phasespace
      ( const double         m      ,
        const double         m1     ,
        const double         m2     ,
        const unsigned short L  = 0 ) ;
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace ;    // integration workspace
      // ======================================================================
    } ;  
    // ========================================================================
    /** @class PhaseSpace3
     *  Function to represent three-body phase space
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class  PhaseSpace3
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from three masses
       *  @param m1 the mass of the first  particle
       *  @param m2 the mass of the second particle
       *  @param m3 the mass of the third  particle
       *  @param l1 the angular momentum between 1st and 2nd particle
       *  @param l2 the angular momentum between the pair and 3rd particle
       */
      PhaseSpace3 ( const double         m1 = 0 ,
                    const double         m2 = 1 ,
                    const double         m3 = 2 ,
                    const unsigned short l1 = 0 ,
                    const unsigned short l2 = 0 ) ;
      /// deststructor
      ~PhaseSpace3 () ;                                         // deststructor
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate 3-body phase space
      double operator () ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      double lowEdge () const { return m_m1 + m_m2 + m_m3 ; }
      // ======================================================================
    public:
      // ======================================================================
      /// helper phase space ("23L")
      double ps2_aux  ( const double m12  ) const ;
      /// get the integral between low and high limits
      double integral ( const double low  ,
                        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag ()  const ;
      // ======================================================================
    private:
      // ======================================================================
      /// the mass of the first particle
      double         m_m1 ; // the mass of the first particle
      /// the mass of the second particle
      double         m_m2 ; // the mass of the second particle
      /// the mass of the third particle
      double         m_m3 ; // the mass of the third  particle
      /// the orbital momentum of the first pair
      unsigned short m_l1 ; // the orbital momentum of the first pair
      /// the orbital momentum between the pair and the third particle
      unsigned short m_l2 ; // the orbital momentum between the pair and the third particle
      // ======================================================================
    private:
      // ======================================================================
      /// the temporary mass
      mutable double m_tmp ; /// the temporary mass
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace  ;    // integration workspace
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace2 ;    // integration workspace
      // ======================================================================
    } ;
    // ========================================================================
    /** @class PhaseSpaceLeft
     *  Function to represent N-body phase space near left-threshold
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class  PhaseSpaceLeft
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from threshold and number of particles and  scale
      PhaseSpaceLeft ( const double                    threshold = 0 ,
                       const unsigned short            num       = 2 , 
                       const double                    scale     = 1 ) ;      
      /// constructor from the list of masses
      PhaseSpaceLeft ( const std::vector<double>&      masses        , 
                       const double                    scale     = 1 ) ;      
      /// special case: true 2-body phasespace 
      PhaseSpaceLeft ( const char*  /* tag */   , 
                       const double m1          , 
                       const double m2          ,
                       const double scale    =  1 ) ;
      /// destructor
      ~PhaseSpaceLeft () ;                                       // destructor
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate N-body phase space near the left threhsold
      double operator   () ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the threshold 
      double threshold  () const { return m_threshold ; }
      /// get the scale 
      double scale      () const { return m_scale     ; }
      // ======================================================================
    public: // integrals
      // ======================================================================
      double integral   ( const double xmin , const double xmax ) const ;
      // ======================================================================
    public:
      // ======================================================================
      bool setThreshold ( const double x ) ;
      bool setScale     ( const double x ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag  
      std::size_t tag   () const ; // get the hash 
      // ======================================================================
    private:
      // ======================================================================
      /// the threshold
      double         m_threshold ; // the threshold
      /// number of particles
      unsigned short m_num       ; // number of particles
      /// the scale  factor 
      double         m_scale     ; // the scale factor
      /// true 2-body phase-space 
      PhaseSpace2    m_ps2       ;
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace {} ;    // integration workspace
      // ======================================================================
    } ;
    // ========================================================================
    /** @class PhaseSpaceRight
     *  simple function to represent L/N-body phase space near right-threshold
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class  PhaseSpaceRight
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from threshold and number of particles
      PhaseSpaceRight ( const double         threshold = 10 ,
                        const unsigned short l         = 2  ,
                        const unsigned short n         = 3  ) ;
      /// deststructor
      ~PhaseSpaceRight () ;                                     // deststructor
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate N/L-body phase space near right  threhsold
      double operator () ( const double x ) const ;
      // ======================================================================
    public: // integrals
      // ======================================================================
      double integral ( const double xmin , const double xmax ) const ;
      // ======================================================================
    public:
      // ======================================================================
      bool setThreshold ( const double x ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag  
      std::size_t tag () const ; // get the hash 
      // ======================================================================
    private:
      // ======================================================================
      /// the threshold
      double         m_threshold ; // the threshold
      /// number of particles
      unsigned short m_N         ; // number of particles
      /// number of particles
      unsigned short m_L         ; // number of particles
      // ======================================================================
    } ;  
    // ========================================================================
    /** @class PhaseSpaceNL
     *  Function epresenting the approximation for
     *  the mass distribution of \f$l\f$-particles from \f$n\f$-body
     *  phase space decay:
     *  \f[ \Phi_{l,n}(x;x_{low},x_{high}) \equiv  
     *       C y^{\frac{3l-5}{2}}\left(1-y\right)^{\frac{3(n-l)-2}{2}}\f]
     *  where 
     *  - \f$ y\equiv \frac{x-x_{low}}{x_{high}-x_{low}}\f$ 
     *  - \f$ C\f$ is a normalization constant, such as 
     *        \f$ \int_{x_{low}}^{x_{high}} \Phi_{l,n}(x) d x = 1\f$ 
     *  - \f$x_{low}= \sum_{i}^{l}m_i\f$ - is a lower threshodl for mass 
     *     of \f$l\f$-particles 
     *  - \f$x_{high}= M - \sum{i=l+1}{n}m_i\f$, is an upper threshold for 
     *     mass of \f$l\f$-particles from \f$n\f$-body decay of the particle
     *     with mass \f$M\f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2011-11-30
     */
    class PhaseSpaceNL
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from thresholds and number of particles
       *  @param low  the low-mass  threshold
       *  @param high the high-mass threshold
       *  @param l    how many particles we consider
       *  @param n    total number of particles ( N>L!)
       */
      PhaseSpaceNL ( const double         low  =  0 ,
                     const double         high = 10 ,
                     const unsigned short l    =  2 ,
                     const unsigned short n    =  3 ) ;
      /// destructor
      ~PhaseSpaceNL () ;                                     // deststructor
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate N/L-body phase space
      double operator () ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      double          lowEdge () const { return m_threshold1 ; }
      double         highEdge () const { return m_threshold2 ; }
      unsigned short       L  () const { return m_L ; }
      unsigned short       N  () const { return m_N ; }
      // ======================================================================
    public:
      // ======================================================================
      /// set the thresholds
      bool setThresholds
      ( const double mn ,
        const double mx ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral
      double integral () const ;
      /// get the integral between low and high limits
      double integral ( const double low  ,
                        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag  
      std::size_t tag () const ; // get the hash 
      // ======================================================================
    private:
      // ======================================================================
      /// the threshold
      double         m_threshold1 ; // the threshold
      double         m_threshold2 ; // the threshold
      /// number of particles
      unsigned short m_N          ; // number of particles
      /// number of particles
      unsigned short m_L          ; // number of particles
      /// normalization
      double         m_norm       ; // normalization
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace ;    // integration workspace
      // ======================================================================
    } ;
    // ========================================================================
    /** @class PSDalitz 
     *  @see Ostap::Kinematics::Dalitz
     */
    class PSDalitz 
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from Dalizt plot          
      PSDalitz ( const Ostap::Kinematics::Dalitz& dalitz ) ;
      /// constructor from all masses  
      PSDalitz ( const double M  = 1 , 
                 const double m1 = 0 , 
                 const double m2 = 0 , 
                 const double m3 = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /** get the value of PDF 
       *  @see Ostap::Kinematics::Dalitz::dRdm12 
       */
      double operator () ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the underlaying Dalitz plot 
      const Ostap::Kinematics::Dalitz& dalitz() const { return m_dalitz ; }
      /// get the overall mass 
      double M  () const { return m_dalitz.M  () ; }
      /// the first mass 
      double m1 () const { return m_dalitz.m1 () ; }
      /// the second mass 
      double m2 () const { return m_dalitz.m2 () ; }
      /// the third mass 
      double m3 () const { return m_dalitz.m3 () ; }
      // ======================================================================
    public:
      // ======================================================================
      double xmin () const { return m1() + m2 () ;  }
      double xmax () const { return M () - m3 () ;  }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral
      double integral () const ;
      /// get the integral between low and high limits
      double integral ( const double low  ,
                        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag  
      std::size_t tag () const ; // get the hash 
      // ======================================================================
    private:
      // ======================================================================
      /// Dalitz plot istself 
      Ostap::Kinematics::Dalitz m_dalitz ; // Dalitz plot istself 
      /// normalization constant
      double                    m_norm   ; // normalization constant
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace ;    // integration workspace
      // ======================================================================
    };
    // ========================================================================
    /** @class PhaseSpace23L
     *  simple function to represent the phase
     *   space of 2 particles from 3-body decays:
     *   \f$ f \propto q^{2\ell+1}p^{2L+1}\f$, where
     *     \f$\ell\f$ is the orbital momentum of the pair of particles,
     *    and \f$L\f$ is the orbital momentum between the pair and
     *    the third particle.
     *   E.g. taking \f$\ell=0, L=1\f$, one can get the S-wave contribution for
     *   \f$\pi^+\pi^-\f$-mass from \f$B^0\rightarrow J/\psi\pi^+\pi^-\f$ decay.
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2012-04-01
     */
    class  PhaseSpace23L
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from four masses and angular momenta
       *  @param m1 the mass of the first  particle
       *  @param m2 the mass of the second particle
       *  @param m3 the mass of the third  particle
       *  @param m  the mass of the mother particle (m>m1+m2+m3)
       *  @param L  the angular momentum between the first pair and
       *  the third particle
       *  @param l  the angular momentum between the first and the second particle
       */
      PhaseSpace23L ( const double         m1 = 0.5 ,
                      const double         m2 = 0.5 ,
                      const double         m3 = 3   ,
                      const double         m  = 5   ,
                      const unsigned short L  = 1   ,
                      const unsigned short l  = 0   ) ;
      /// deststructor
      ~PhaseSpace23L () ;                                     // deststructor
      // ======================================================================
    public:
      // ======================================================================
      /// calculate the phase space
      double operator () ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// calculate the phase space
      double ps23L ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      double m1        () const { return m_m1 ; }
      double m2        () const { return m_m2 ; }
      double m3        () const { return m_m3 ; }
      double m         () const { return m_m  ; }
      unsigned short l () const { return m_l  ; }
      unsigned short L () const { return m_L  ; }
      // ======================================================================
      double lowEdge   () const { return m1 () + m2 () ; }
      double highEdge  () const { return m  () - m3 () ; }
      /// get the momentum of 1st particle in rest frame of (1,2)
      double         q ( const double x ) const ;
      /// get the momentum of 3rd particle in rest frame of mother
      double         p ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral
      double integral () const ;
      /// get the integral between low and high limits
      double integral ( const double low  ,
                        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag  
      std::size_t tag () const ; // get the hash 
      // ======================================================================
    private:
      // ======================================================================
      /// the first mass
      double m_m1 ;            // the first mass
      /// the second mass
      double m_m2 ;            // the second mass
      /// the third  mass
      double m_m3 ;            // the third mass
      /// the mass of mother particle
      double m_m  ;            // the mass of mother particle
      /// the orbital momentum between the first and the second particle
      unsigned short m_l ; // the orbital momentum between the 1st and 2nd
      /// the orbital momentum between the pair an dthe third particle
      unsigned short m_L ; // the orbital momentum between the (12) and 3rd
      // ======================================================================
      /// helper normalization parameter
      double m_norm ; // helper normalization parameter
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace ;    // integration workspace
      // ======================================================================
    } ;
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_PHASESPACE_H
// ============================================================================
