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
      inline double operator () ( const double m ) const 
      { return phasespace ( m , m_m1 , m_m2 ) ; }
      /// integral
      double integral    ( const double xmin , const double xmax ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get a phase space 
      inline double        rho  ( const double m ) const 
      { return phasespace ( m , m_m1 , m_m2 ) ;}
      // ======================================================================
      /** get (a complex) phase space 
       *  real for x >= threshold , imaginary for x< threshold 
       */
      std::complex<double> rho1 ( const double m ) const { return rho1_s (  m * m ) ; }
      // ======================================================================
    public: 
      // ======================================================================
      /// get a phase space as function of s 
      inline double        rho_s  ( const double s ) const 
      { return phasespace_s ( s , m_m1 * m_m1 , m_m2 * m_m2 ) ;}
      // ======================================================================
      /** get (a complex) phase space 
       *  real for x >= threshold , imaginary for x< threshold 
       */
      std::complex<double> rho1_s ( const double s ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the momentum at center of mass
      inline double                q    ( const double m ) const
      { return q    ( m , m_m1 , m_m2 ) ; }
      /// ditto but as complex
      inline std::complex<double>  q1   ( const double m ) const 
      { return q1   ( m , m_m1 , m_m2 ) ; }
      /// get the momentum at given s 
      inline double                q_s  ( const double s ) const 
      { return q_s  ( s , m_m1 * m_m1 , m_m2 * m_m2 ) ; }
      /// ditto but as complex 
      inline std::complex<double>  q1_s ( const double s ) const 
      { return q1_s ( s , m_m1 * m_m1 , m_m2 * m_m2 ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /** get the mass for the given momentum
       *  \f$ m = \sqrt{m_1^2+q^2} + \sqrt{m_2^2+q^2}\f$
       */
      double q2m ( const double q ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// the first mass 
      inline double m1        () const { return m_m1          ; }
      /// the second mass 
      inline double m2        () const { return m_m2          ; }
      /// threshold 
      inline double lowEdge   () const { return m1 () + m2 () ; }
      /// threshold 
      inline double threshold () const { return m1 () + m2 () ; }
      ///  threshols for s 
      inline double s_threshold () const 
      { const double a = threshold() ; return a * a ; }
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
      /** calculate the particle momentum in the rest frame
       *  \f[ q = \frac{1}{2}\frac{ \lambda^{\frac{1}{2}}
       *        \left( m^2 , m_1^2, m_2^2 \right) }{ m }\f],
       *  @param s    the squared mass
       *  @param m2_1 the squared mass of the first particle
       *  @param m2_2 the srqured mass of the second particle
       *  @return the momentum in rest frame (physical values only)
       */
      static double  q_s
      ( const double s    ,
        const double m2_1 ,
        const double m2_2 ) ;
      // ======================================================================
      /** calculate the particle momentum in rest frame
       *  - real for physical case 
       *  - imaginary for non-physical case (below the threshold)
       *  @param s    the squared mass
       *  @param m2_1 the squared  mass of the first particle
       *  @param m2_2 the squared mass of the second particle
       *  @return the momentum in rest frame  (imaginary for non-physical branch)
       */
      static std::complex<double> q1_s
      ( const double s    ,
        const double m2_1 ,
        const double m2_2 ) ;
      // ======================================================================
      /** calculate the phase space for   m -> m1 + m2
       *  \f[ \Phi = \frac{1}{8\pi} \left( \frac{ \lambda^{\frac{1}{2}}
       *       \left( m^2 , m_1^2, m_2^2 \right) }{ m^2 }\right)^{2L+1}\f],
       *  where \f$\lambda\f$ is a triangle function
       *  @param m the mass
       *  @param m1 the mass of the first particle
       *  @param m2 the mass of the second particle
       *  @return two-body phase space
       */
      static double phasespace
      ( const double         m        ,
        const double         m1       ,
        const double         m2       ,  
        const unsigned short L    = 0 ) 
      { return phasespace_s ( m * m , m1 * m1 , m2 * m2 , L ) ; }
      // ======================================================================
      /** calculate the phase space for   m -> m1 + m2
       *  \f[ \Phi = \frac{1}{8\pi} \left( \frac{ \lambda^{\frac{1}{2}}
       *       \left( m^2 , m_1^2, m_2^2 \right) }{ m^2 }\right)^{2L+1}\f],
       *  where \f$\lambda\f$ is a triangle function
       *  @param s    the squared mass
       *  @param m2_1 the squared mass of the first particle
       *  @param m2_1 the squared mass of the second particle
       *  @return two-body phase space
       */
      static double phasespace_s
      ( const double         s        ,
        const double         m2_1     ,
        const double         m2_2     , 
        const unsigned short L    = 0 ) ;
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace ;    // integration workspace
      // ======================================================================
    } ;  
    // ========================================================================
    /** @class sPhaseSpace2
     *  Two-body phase space as function of s 
     *  @see Ostap::Math::PhaseSpace2
     *  @see Ostap::Math::PhaseSpace2::phasespace_s
     *  @see Ostap::Math::PhaseSpace2::phasespace_s
     */
    class sPhaseSpace2 
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from two masses
      sPhaseSpace2 ( const double m1 = 0 , 
                     const double m2 = 1 ) ;
      // ======================================================================
      /// Two-body phase space as function of s
      inline double operator () ( const double s ) const 
      { return PhaseSpace2::phasespace_s ( s , m_m2_1 , m_m2_2 ) ; }
      // ======================================================================
    private: 
      // ======================================================================
      /// the first mass squared
      double m_m2_1 ; // the first mass squared
      /// the second mass squared
      double m_m2_2 ; // the second mass squared
      // ======================================================================
    } ;
    // ========================================================================
    /** @class PhaseSpace3s
     *  Symmetric form of 3-body phase space 
     *  @see Davydychev, A and Delbourgo, R., 
     *       "Three body phase space: Symmetrical treatments",
     *        "{15th Biennial Congress of the Australian Institute of
     *         Physics Sydney, Australia, July 8-11, 2002}",
     *  @see http://arxiv.org/abs/hep-th/0209233
     *
     * three-body phase space, analytic symmetric expression via 
     *  elliptic  integrals 
     *  @see https://indico.cern.ch/event/368497/contributions/1786992/attachments/1134067/1621999/davydychev.PDF
     *  @see http://cds.cern.ch/record/583358/files/0209233.pdf
     *  @see https://www.researchgate.net/publication/2054534_Three-body_phase_space_symmetrical_treatments
     *  @see Ostap::Kinematics::phasespace3     
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2019-10-31
     */
    class PhaseSpace3s
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from three masses
       *  @param m1 the mass of the first  particle
       *  @param m2 the mass of the second particle
       *  @param m3 the mass of the third  particle
       */
      PhaseSpace3s ( const double         m1 = 0 ,
                     const double         m2 = 1 ,
                     const double         m3 = 2 ) ;
      // ======================================================================
    public:
      // ======================================================================      
      /// evaluate 3-body phase space
      double evaluate    ( const double x ) const ;
      /// evaluate 3-body phase space
      double operator () ( const double x ) const { return  evaluate ( x )  ; }
      // ======================================================================
    public:
      // ======================================================================
      double m1 () const { return m_m1 ; }
      double m2 () const { return m_m2 ; }
      double m3 () const { return m_m3 ; }
      // ======================================================================
      double lowEdge () const { return m_m1 + m_m2 + m_m3 ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between low and high limits
      double integral ( const double low  ,
                        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /** three-body phase space, analytic symmetric expression via 
       *  elliptic  integrals 
       *  @see https://indico.cern.ch/event/368497/contributions/1786992/attachments/1134067/1621999/davydychev.PDF
       *  @see http://cds.cern.ch/record/583358/files/0209233.pdf
       *  @see https://www.researchgate.net/publication/2054534_Three-body_phase_space_symmetrical_treatments
       *  @see Ostap::Kinematics::phasespace3
       */
      static double phasespace ( const double x  , 
                                 const double m1 , 
                                 const double m2 , 
                                 const double m3 ) ;
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
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace  {} ; // integration workspace
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
      /** constructor from three masses
       *  @param l1 the angular momentum between 1st and 2nd particle
       *  @param l2 the angular momentum between the pair and 3rd particle
       */
      PhaseSpace3 ( const PhaseSpace3s&  ps3    , 
                    const unsigned short l1 = 0 ,
                    const unsigned short l2 = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /** evaluate 3-body phase space
       *  \f[ R_3 ( M ) = \frac{pi^2}{4M^2}\int_{m2+m3}^{M-m_1} \drac{ds_2}{s_2}
       *   \lambda^{1/2}\left ( s_2 , M^2   , m_1^2\right) 
       *   \lambda^{1/2}\left ( s_2 , m_2^2 , m_3^2\right) 
       *  \f] 
       *  @see E.Byckling, K.Kajantie, "Particle kinematics", John Wiley & Sons,
       *              London, New York, Sydney, Toronto, 1973, Eq. (V.2.17)
       */
      double evaluate    ( const double x ) const ;
      // ======================================================================
      /** evaluate 3-body phase space
       *  \f[ R_3 ( M ) = \frac{pi^2}{4M^2}\int_{m2+m3}^{M-m_1} \drac{ds_2}{s_2}
       *   \lambda^{1/2}\left ( s_2 , M^2   , m_1^2\right) 
       *   \lambda^{1/2}\left ( s_2 , m_2^2 , m_3^2\right) 
       *  \f] 
       *  @see E.Byckling, K.Kajantie, "Particle kinematics", John Wiley & Sons,
       *              London, New York, Sydney, Toronto, 1973, Eq. (V.2.17)
       */
      double operator () ( const double x ) const { return evaluate  ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      double m1 () const { return m_m1 ; }
      double m2 () const { return m_m2 ; }
      double m3 () const { return m_m3 ; }
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
    public:
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
      mutable double m_tmp { 0 } ; /// the temporary mass
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
      PhaseSpaceLeft ( const PhaseSpace2& ps2        , 
                       const double       scale =  1 ) ;
      /// destructor
      ~PhaseSpaceLeft () ;                                       // destructor
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate N-body phase space near the left threshold
      double operator   () ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the threshold 
      double         threshold  () const { return m_threshold ; }
      /// get the number of particles : 0 means true 2-body! 
      unsigned short N          () const { return m_num       ; }
      /// get the scale 
      double         scale      () const { return m_scale     ; }
      // ======================================================================
      const PhaseSpace2& ps2    () const { return m_ps2       ; }
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
      double         m_threshold         ; // the threshold
      /// number of particles
      unsigned short m_num       { 0   } ; // number of particles
      /// the scale  factor 
      double         m_scale     { 1.0 } ; // the scale factor
      /// true 2-body phase-space 
      PhaseSpace2    m_ps2       {}      ;
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
     *  for  \f$ 2\le l < n \$ it is defined as  
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
     * 
     *  It also includes two specific cases 
     *  - \f$ 0 == l \f$, \f$ 1 \le n \f$
     *  \f[ \Phi_{0,n}(x;x_{low},x_{high}) \propto   
     *     \left(1-y\right)^{ \frac{3n-2}{2} } \f]
     *  - \f$ 2 \le \f$, \f$ n = 0  \f$
     *  \f[ \Phi_{l,0}(x;x_{low},x_{high}) \propto   
     *     \left(y\right)^{ \frac{3l-5}{2} } \f]
     *
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
      unsigned short       L  () const { return m_L          ; }
      unsigned short       N  () const { return m_N          ; }
      // ======================================================================
      double xmin () const { return  lowEdge () ; }
      double xmax () const { return highEdge () ; }      
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
      /** Get a full integrated phase space over Dalitz plot 
       *  \f$  R(s) = \int \int R(s_1,s_2) {\mathrm{d}} s_1 {\mathrm{d}} s_2 =
       *  \int _{(m_2+m_3)^2}^{ (\sqrt{s}-m_1)^2}
       *   \frac{{\mathrm{d}} s_2}{s_2}
       *   \lambda^{1/2}(s_2,s,m_1^2)
       *   \lambda^{1/2}(s_2,m_2^2,m_3^2)\f$ 
       */
      double phasespace () const { return 1.0 / m_norm ; }
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
    /** @class M2Q 
     *  \f$ m \rightarrow q \f$ transformation
     *  @see Ostap::Math::PhaseSpace2::q
     *  @see Ostap::Math::PhaseSpace2::q_q
     */
    class M2Q  
    {
    public :
      // ======================================================================
      /// constructor from two masses 
      M2Q ( const double m1 = 0 , const double m2 = 0 ) ;
      /// the only one important method 
      inline double operator () ( const double m ) const 
      { return Ostap::Math::PhaseSpace2::q_s ( m * m , m_m2_1 , m_m2_2 ) ; }
      // ======================================================================
    private:
      // ======================================================================
      /// the first mass squared 
      double m_m2_1 { 0 } ;
      /// the second mass squared 
      double m_m2_2 { 0 } ;
      // ======================================================================
    };
    // ========================================================================
    /** @class S2Q 
     *  \f$ s \rightarrow q \f$ transformation
     *  @see Ostap::Math::PhaseSpace2::q_s
     */
    class S2Q  
    {
    public :
      // ======================================================================
      /// constructor from two masses 
      S2Q ( const double m1 = 0 , const double m2 = 0 ) ;
      /// the only one important method 
      inline double operator () ( const double s ) const 
      { return Ostap::Math::PhaseSpace2::q_s ( s , m_m2_1 , m_m2_2 ) ; }
      // ======================================================================
    private:
      // ======================================================================
      /// the first mass squared 
      double m_m2_1 { 0 } ;
      /// the second mass squared 
      double m_m2_2 { 0 } ;
      // ======================================================================
    };
    // ========================================================================
    /** @class Q2M 
     *  \f$ q \rightarrow m \f$ transformation
     *  @see Ostap::Math::PhaseSpace2::q
     */
    class Q2M  
    {
    public :
      // ======================================================================
      /// constructor from two masses 
      Q2M ( const double m1 = 0 , const double m2 = 0 ) ;
      /// the only one important method 
      inline double operator () ( const double q ) const 
      { 
        const double q2 = q <= 0 ? 0.0 : q * q ;
        return std::sqrt ( m_m2_1 + q2 )  + std::sqrt ( m_m2_2 + q2 ) ;
      }
      // ======================================================================
    private:
      // ======================================================================
      /// the first mass squared 
      double m_m2_1 { 0 } ;
      /// the second mass squared 
      double m_m2_2 { 0 } ;
      // ======================================================================
    };
    // ========================================================================
    /** @class Q2S 
     *  \f$ q \rightarrow s \f$ transformation
     *  @see Ostap::Math::PhaseSpace2::q
     */
    class Q2S  
    {
    public :
      // ======================================================================
      /// constructor from two masses 
      Q2S ( const double m1 = 0 , const double m2 = 0 ) ;
      /// the only one important method 
      inline double operator () ( const double q ) const 
      { 
        const double q2   = q <= 0 ? 0.0 : q * q ;
        const double e2_1 = m_m2_1 + q2 ;
        const double e2_2 = m_m2_2 + q2 ;
        //
        return e2_1 + e2_2 + 2 * std::sqrt ( e2_1 * e2_2 ) ;
      }
      // ======================================================================
    private:
      // ======================================================================
      /// the first mass squared 
      double m_m2_1 { 0 } ;
      /// the second mass squared 
      double m_m2_2 { 0 } ;
      // ======================================================================
    };
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_PHASESPACE_H
// ============================================================================
