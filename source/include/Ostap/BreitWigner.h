// ============================================================================
#ifndef OSTAP_BREITWIGNER_H 
#define OSTAP_BREITWIGNER_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <cmath>
#include <complex>
#include <memory>
#include <functional>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Workspace.h"
#include "Ostap/PhaseSpace.h"
#include "Ostap/Dalitz.h"
#include "Ostap/Models.h"
// ============================================================================
// ROOT
// ============================================================================
#include "RVersion.h"
// ============================================================================
/** @file Ostap/BreitWigner.h
 *
 * Set of useful models for describing signal peaks with the natural width:
 *  - Breit-Wigner
 *  - Flatte 
 *  - LASS  (kappa) 
 *  - Bugg  (sigma-pole)
 *  - Gounaris-Sakurai
 *  - ...
 *
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
    /// base class for formfactors
    class FormFactor ;
    // ========================================================================
    /** @namespace Ostap::Math::FormFactors
     *  form-factor functions for Breit-Wigner
     *   - Blatt-Weisskopf form-factors 
     *   - Jackson's form-factors 
     */
    namespace FormFactors
    {
      // ======================================================================
      /** @typedef rho_fun
       *  the \f$\rho(\omega)\f$ function from Jackson
       *  Arguments
       *    - the        mass
       *    - the pole   mass
       *    - the first  daughter mass
       *    - the second daughter mass
       */
      typedef std::function<double(double,double,double,double)> rho_fun ;
      // ======================================================================
      /** parameterization for \f$\rho(\omega)\f$-function from (A.1)
       *  J.D.Jackson,
       *  "Remarks on the Phenomenological Analysis of Resonances",
       *  In Nuovo Cimento, Vol. XXXIV, N.6
       */
      enum JacksonRho {
        Jackson_0  = 0 , /// \f$\rho(\omega) = 1 \f$
        Jackson_A2     , /// \f$ 1^- \rightarrow 0^- 0^- \f$ , l = 1
        Jackson_A3     , /// \f$          1^- \rightarrow 0^- 1^- \f$ , l = 1
        Jackson_A4     , /// \f$ \frac{3}{2}+ \rightarrow 0^- \frac{1}{2}^+ \f$ , l = 1
        Jackson_A5     , /// \f$ \frac{3}{2}- \rightarrow 0^- \frac{1}{2}^+ \f$ , l = 2
        Jackson_A7       /// recommended for rho0 -> pi+ pi-
      } ;
      // ======================================================================
    } //                          the end of namespace Ostap::Math::FormFactors
    // ========================================================================
    /** @class ChannelBW
     *  Simple definition of the decay channel for the Breit-Wigner function 
     *  It defines three function 
     *  - \f$ N^2_a(s)     \f$  
     *  - \f$ D_a(s)       \f$
     *  - \f$ \varrho_a(s) \f$
     *  
     *  With these factors for the channel \f$ a \f$ the BW-amplitude is 
     *  \f[ \mathcal{A}_a (s) \propto \frac{1} { m_0^2 - s - iD_a(s,m_0) } \f] 
     *  The amplitude can be scaled with \f$ N_a(s,m_0)\f$, if needed.  
     *
     *  The final mass distribution is 
     *  \f[ F(m) \propto \frac{2 m}{\pi} \varrho ( m^2 ) N^2 ( m^2 , m_0 ) \left| \mathcal{A}(m^2) \right|^2 \f$ 
     * 
     *  For many simple cases one has \f$ \varrho (s ) N^2(s,m_0) = D(s,m_0) \$.  
     *
     *  For multi-channel case the amplitude in the channel \f$ a \f$ is 
     *  \f[ \mathcal{A}_a (s) = \frac{1}{ m_0^2 - s - i\sum_b D_b(s,m_0) } \f] 
     *  
     */ 
    class ChannelBW
    {
      // =======================================================================
    public:
      // =======================================================================
      /** constructor from the partial width (or squared coupling constant).
       *  - For above threshold channels \f$ m^2_0 >s_{threshold} \f$ 
       *  it it partial width 
       *  - For sub threshold channels \f$ m^2_0 < s_{threshold} \f$ 
       *  it is proportional to the squared coupling constant 
       */ 
      ChannelBW ( const double gamma  ) ; 
      ChannelBW ( const ChannelBW&    ) = default ;
      ChannelBW (       ChannelBW&&   ) = default ;      
      /// virtual destructor 
      virtual ~ChannelBW () ; // virtual destructor 
      /// clone-method 
      virtual  ChannelBW*           clone () const = 0 ;
      // =======================================================================
    public:  // the main methods 
      // =======================================================================
      /// squared  numerator for the amplitude 
      virtual double               N2
      ( const double s  , 
        const double m0 ) const = 0 ;      
      /// term in the denominator for the amplitide
      virtual std::complex<double> D    
      ( const double s  , 
        const double m0 ) const = 0 ;      
      /** get the phase space factor  \f$ \varrho(s) \f$
       *  optionally normalized at the point \f$ m_n \f$ 
       */
      virtual double rho_s 
      ( const double s  , 
        const double mn ) const = 0 ;
      /// get the opening threshold \f$ s_{threshold} \$ for the channel 
      virtual double s_threshold () const = 0 ;
      // =======================================================================
    public: //  helper methods 
      // =======================================================================
      /// unique tag/label  
      virtual std::size_t tag       () const = 0 ;
      /// describe the channel 
      virtual std::string describe  () const = 0 ;
      // =======================================================================
    public: // interpret it as (a partial) width 
      // =======================================================================
      /// get the partial width for the channel 
      virtual double  gamma0 () const { return m_gamma0 ; } // get the partial width
      /// set the partial width for this channel 
      virtual bool setGamma0 ( const double value ) ;       // set the partial width 
      // =======================================================================
    public: // interpret it as (a squared) coupling constant 
      // =======================================================================
      /// squared coupling constant 
      inline double g2       () const { return gamma0 () ; } // squared coupling constant 
      /// set a squared coupling constant 
      inline bool   setG2    ( const double value ) { return setGamma0 ( value ) ; }
      // =======================================================================
    public:
      // =======================================================================
      /** get the single channel amplitude
       *  \f[ \mathcal{A} =\left( m_0^2 - s - i D ( s , m_0 ) \right)^{-1} \f] 
       *  @param s   \f$ s  \f$ -parameter
       *  @param m0  \f$ m_0\f$-parameter  
       *  @return amplitude 
       */
      std::complex<double> amplitude ( const double s , const double m0 ) const ;
      // =======================================================================
    private : 
      // =======================================================================
      /// the decay width for this channel 
      double m_gamma0 { 0 } ; // the decay width for this channel 
      // =======================================================================      
    } ;
    //==========================================================================
// #if ROOT_VERSION_CODE >= ROOT_VERSION(6,23,1)
//     // =========================================================================
//     /** @class ChannelGeneric
//      *  Generic description of the channel. 
//      *  \f[ \begin{array}{lcl}
//      *   N^2(s,m_0)      & = & m_0 \Gamma_0 f_{N^2}(s,m_0) \\
//      *   D (s,m_0)       & = & m_0 \Gamma_0 f_{D}(s,m_0)   \\
//      *  \varrho (s, m_n) & = & \Theta\left(s-s_{\mathrm{threshold}}\right) f_{\varrho}(s,m_n)
//      *  \end{array}\,,\f]
//      *  where \f$ f_{N^2}\f$,  \f$ f_{D}\f$ and \f$ f_{\varrho}\f$ are provdied externally
//      */
//     class ChannelGeneric : public ChannelBW 
//     {
//     public:
//       // =======================================================================
//       /** @typedef Fun_N2 
//        *  the function type for \f$ N^2(s,m_0)\f$
//        */
//       typedef std::function<double(double,double)>                 Fun_N2  ;
//       // =======================================================================
//       /** @typedef Fun_D 
//        *  the function type for \f$ D(s,m_0)\f$
//        */
//       typedef std::function<std::complex<double>(double,double)>   Fun_D   ;
//       // =======================================================================
//       /** @typedef Fun_DR 
//        *  the (real-valued) type function type for \f$ D(s,m_0)\f$
//        */
//       typedef std::function<double(double,double)>                 Fun_DR  ;
//       // ======================================================================
//       /** @typedef Fun_rho 
//        *  the function type for \f$ \varrho(s,m_n) \f$
//        */
//       typedef std::function<double(double,double)>                 Fun_rho ;
//       // ======================================================================
//     public : 
//       // ======================================================================
//       /// full constructor with all functions specified
//       ChannelGeneric ( const double       gamma                   , 
//                        const Fun_N2&      fN2                     , 
//                        const Fun_D&       fD                      , 
//                        const Fun_rho&     frho                    , 
//                        const double       sthreshold              , 
//                        const std::size_t  tag                     ,
//                        const std::string& description = "Generic" ) ;    
//       // ======================================================================= 
//       /// Constructor  with only real-valued functions 
//       ChannelGeneric ( const double       gamma                   , 
//                        const Fun_N2&      fN2                     , 
//                        Fun_DR             fDr                     , 
//                        const Fun_rho&     frho                    , 
//                        const double       sthreshold              , 
//                        const std::size_t  tag                     ,
//                        const std::string& description = "Generic" ) ;
//       // =======================================================================
//       /// short constructor  with real-valued functions 
//       ChannelGeneric ( const double       gamma                   ,
//                        Fun_DR             fDr                     , 
//                        const double       sthreshold              , 
//                        const std::size_t  tag                     ,
//                        const std::string& description = "Generic" ) ;
//       // =======================================================================
//       ///  copy constructor
//       ChannelGeneric ( const ChannelGeneric& right ) = default ;
//       // =======================================================================
//       /// clone method
//       ChannelGeneric*  clone() const override ; // clone method
//       // =======================================================================
//     public:
//       // =======================================================================
//       template <typename... ARGS>
//       static inline ChannelGeneric
//       create ( const double gamma ,
//                ARGS ...     args  ) { return ChannelGeneric ( gamma , args ... ) ; }
//       // =======================================================================
//     public:
//       // =======================================================================
//       /** squared  numerator for the amplitude 
//        * \f[ N^2(s,m_0) = m_0 \Gamma0 f_{N^2}(s,m_0)\f]
//        */
//       double               N2
//       ( const double s  , 
//         const double m0 ) const override { return m0 * gamma0 () * m_fN2 ( s , m0 ) ; }      
//       // ======================================================================
//       /** term in the denominator for the amplitide
//        * \f[ D (s,m_0) = m_0 \Gamma0 f_{D}(s,m_0)\f]
//        */
//       // ======================================================================
//       std::complex<double> D    
//       ( const double s  , 
//         const double m0 ) const override { return m0 * gamma0 () * m_fD  ( s , m0 ) ; }
//       // ======================================================================
//       /** get the phase space factor  \f$ \varrho(s) \f$
//        *  optionally normalized at the point \f$ m_n \f$ 
//        * \f[ \varrho (s, m_n) = \Theta\left(s-s_{\mathrm{threshold}}\right) f_{\varrho}(s,m_n)\f] 
//        */
//       double rho_s 
//       ( const double s  , 
//         const double mn ) const override 
//       { return s <= m_sthreshold ? 0.0 : m_frho ( s , mn ) ; }
//       /// get the opening threshold \f$ s_{\mathrm{threshold}} \$ for the channel 
//       double s_threshold () const override { return m_sthreshold ; }
//       // =======================================================================
//     public: //  helper methods 
//       // =======================================================================
//       /// unique tag/label  
//       std::size_t tag       () const override ;
//       /// describe the channel 
//       std::string describe  () const override ;
//       // =======================================================================
//     private :
//       // =======================================================================
//       /// function N2 
//       Fun_N2  m_fN2             ; // function N2 
//       /// function D
//       Fun_D   m_fD              ; // function D
//       /// function rho
//       Fun_rho m_frho            ; // function rho
//       /// s-threhold 
//       double  m_sthreshold      ; // s-threhold 
//       /// unique tag 
//       std::size_t m_tag         ; // unique tag 
//       /// description 
//       std::string m_description ; // description 
//       // =======================================================================
//     } ;
//     // =========================================================================
// #endif

    // =========================================================================
    /** @class ChannelCW
     *  Trivial "constant-width" channel 
     *  \f[ \begin{array}{ncl}
     *  N & = & m_0 \Gamma_0 \\ 
     *  D & = & m_0 \Gamma_0 
     *  \end{array} \f]
     *  - the constant width  
     *  - masses of daughter particles
     */ 
    class ChannelCW : public ChannelBW 
    {
    public:
      // ======================================================================
      /** constructor from all parameters and *NO* formfactor
       *  @param gamma the width 
       *  @param m1    the mass of the 1st daughter
       *  @param m2    the mass of the 2nd daughter
       */
      ChannelCW
      ( const double                  gamma = 0.150 , 
        const double                  m1    = 0.139 , 
        const double                  m2    = 0.139 ) ;
      // ======================================================================
      /// clone the channel 
      ChannelCW* clone () const override ;
      // ======================================================================
    public: // trivial accessors  
      // ======================================================================
      /// get the mass of the 1st daughter 
      double            m1         () const { return m_ps2.m1 ()         ; }
      /// get the mass of the 2nd daughter  
      double            m2         () const { return m_ps2.m2 ()         ; }
      /// phase space function 
      const Ostap::Math::PhaseSpace2& ps2 () const { return m_ps2 ; }
      // ======================================================================
    public: // two main methods
      // ======================================================================
      /** the first main method: numerator
       *   \f[ N^2(s,m_0) = m_0\Gamma_0\f]
       */
      double N2 
      ( const double s  , 
        const double m0 ) const override ;
      // ======================================================================
      /** the second main method: term to the denominator 
       *  \f[ D(s,m_0) = m_0\Gamma_0 \f] 
       */
      std::complex<double> D 
      ( const double s  , 
        const double m0 ) const override ;
      // ======================================================================
      /** get the phase space factor  \f$ \varrho \f$, 
       *  (optionally normalized at point \f$ m_n \f$)
       */
      double rho_s 
      ( const double s  , 
        const double mn ) const override ;
      /// get the opening threshold \f$ s_{threshold} \$ for the channel 
      double s_threshold () const override ;
      // =======================================================================
    public:
      // ======================================================================
      // unique tag
      std::size_t tag      () const override ; // unique tag
      // ======================================================================
      /// describe the channel 
      std::string describe () const override ;
      // ======================================================================
    protected : 
      // ======================================================================
      /// two body phase space 
      Ostap::Math::PhaseSpace2                 m_ps2 ; // two body phase space 
      // ======================================================================
    };
    // =========================================================================
    /** @class ChannelWidth
     *  Description of the channel with generic mass-dependent width 
     *  \f[ \begin{array}{lcl}
     *   N^2(s,m_0)      & = & m_0 \Gamma0 \frac{w(s)}{w(m^2_0)} \\
     *   D (s,m_0)       & = & m_0 \Gamma0 \frac{w(s)}{w(m^2_0)} \\
     *  \varrho (s, m_n) & = & \Theta\left(s-s_{threshold}\right)
     *  \end{array}\,,\f]
     */
    class ChannelWidth : public ChannelBW 
    {
    public:
      // =======================================================================
      /** @typedef Width 
       *  the function type for the mass-dependent width
       */
      typedef std::function<double(double)> Width ;
      // =======================================================================
    public : 
      // ======================================================================
      /// full constructor with all functions specified 
      ChannelWidth
      ( const double       gamma                        ,
        Width              width                        ,  
        const double       sthreshold                   , 
        const std::size_t  tag                          ,
        const std::string& description = "ChannelWidth" ) ;
      // =======================================================================
      /// templated constructor
      template <class WIDTH>
      ChannelWidth 
      ( const double       gamma                        ,
        WIDTH              width                        ,  
        const double       sthreshold                   , 
        const std::size_t  tag                          ,
        const std::string& description = "ChannelWidth" ) 
        : ChannelBW     ( gamma ) 
        , m_w           ( width ) 
        , m_sthreshold  ( std::abs ( sthreshold ) ) 
        , m_tag         ( tag   ) 
        , m_description ( description ) 
      {}
      /// clone the channel 
      ChannelWidth* clone () const override ;
      // ======================================================================
    public:
      // =======================================================================
      template <class WIDTH>
      static inline ChannelWidth 
      create
      ( const double       gamma                        , 
        WIDTH              width                        , 
        const double       sthreshold                   , 
        const std::size_t  tag                          ,
        const std::string& description = "ChannelWidth" ) 
      { return ChannelWidth ( gamma , width , sthreshold , tag , description ) ; }
      // ======================================================================
    public:
      // =======================================================================
      /** squared  numerator for the amplitude 
       * \f[ N^2(s,m_0) = m_0 \Gamma0 \frac{w(s)}{w(m^2_0)} \f] 
       */
      double N2
      ( const double s  , 
        const double m0 ) const override 
      { return s <= m_sthreshold ? 0.0 : m0 * gamma0 () * m_w ( s ) / m_w ( m0 * m0 ) ; }      
      // ======================================================================
      /** term in the denominator for the amplitide
       * \f[ D (s,m_0) = m_0 \Gamma0 \frac{w(s)}{w(m^2_0)} \f] 
       */
      // ======================================================================
      std::complex<double> D    
      ( const double s  , 
        const double m0 ) const override 
      { return s <= m_sthreshold ? 0.0 : m0 * gamma0 () * m_w ( s ) / m_w ( m0 * m0 ) ; }      
      // ======================================================================
      /** get the phase space factor  \f$ \varrho(s) \f$
       *  optionally normalized at the point \f$ m_n \f$ 
       * \f[ \varrho (s, m_n) = \Theta\left(s-s_{threshold}\right) \f] 
       */
      double rho_s 
      ( const double s  , 
        const double mn ) const override 
      { return s <= m_sthreshold ? 0.0 : 1 ; }
      /// get the opening threshold \f$ s_{threshold} \$ for the channel 
      double s_threshold () const override { return m_sthreshold ; }
      // =======================================================================
    public: //  helper methods 
      // =======================================================================
      /// unique tag/label  
      std::size_t tag       () const override ;
      /// describe the channel 
      std::string describe  () const override { return m_description ; }
      // =======================================================================
    private :
      // =======================================================================
      /// mass-dependent width 
      Width       m_w           ; // mass-dependent width 
      /// s-threshold 
      double      m_sthreshold  ; // s-threshold 
      /// unique tag/lable 
      std::size_t m_tag         ; // unique tag/lable 
      /// description 
      std::string m_description ;
      // ======================================================================
    } ;   
    // =========================================================================
    /** @class ChannelGamma
     *  Description of the channel with generic mass-dependent width 
     *  \f[ \begin{array}{lcl}
     *   N^2(s,m_0)      & = & m_0 \Gamma0 \gamma(s) } \          \
     *   D (s,m_0)       & = & m_0 \Gamma0 \gamma(s) } \          \
     *  \varrho (s, m_n) & = & \Theta\left(s-s_{threshold}\right)
     *  \end{array}\,,\f]
     */
    class ChannelGamma : public ChannelBW 
    {
    public:
      // =======================================================================
      /** @typedef Width 
       *  the function type for the mass-dependent width
       */
      typedef std::function<double(double)> Width ;
      // =======================================================================
    public : 
      // ======================================================================
      /// full constructor with all functions specified 
      ChannelGamma
      ( const double       gamma                        ,
        Width              width                        ,  
        const double       sthreshold                   , 
        const std::size_t  tag                          ,
        const std::string& description = "ChannelGamma" ) ;
      // =======================================================================
      /// templated constructor
      template <class WIDTH>
      ChannelGamma 
      ( const double       gamma                        ,
        WIDTH              width                        ,  
        const double       sthreshold                   , 
        const std::size_t  tag                          ,
        const std::string& description = "ChannelGamma" ) 
        : ChannelBW     ( gamma ) 
        , m_gamma       ( width ) 
        , m_sthreshold  ( std::abs ( sthreshold ) ) 
        , m_tag         ( tag   ) 
        , m_description ( description ) 
      {}
      /// clone the channel 
      ChannelGamma* clone () const override ;
      // ======================================================================
    public:
      // =======================================================================
      template <class WIDTH>
      static inline ChannelGamma 
      create
      ( const double       gamma                        , 
        WIDTH              width                        , 
        const double       sthreshold                   , 
        const std::size_t  tag                          ,
        const std::string& description = "ChannelGamma" ) 
      { return ChannelGamma ( gamma , width , sthreshold , tag , description ) ; }
      // ======================================================================
    public:
      // =======================================================================
      /** squared  numerator for the amplitude 
       * \f[ N^2(s,m_0) = m_0 \Gamma0 \gamma(s) \f] 
       */
      double N2
      ( const double s  , 
        const double m0 ) const override 
      { return s <= m_sthreshold ? 0.0 : m0 * gamma0 () * m_gamma ( s ) ; }      
      // ======================================================================
      /** term in the denominator for the amplitide
       * \f[ D (s,m_0) = m_0 \Gamma0 \gamma (s) \f] 
       */
      // ======================================================================
      std::complex<double> D    
      ( const double s  , 
        const double m0 ) const override 
      { return s <= m_sthreshold ? 0.0 : m0 * gamma0 () * m_gamma ( s ) ; }      
      // ======================================================================
      /** get the phase space factor  \f$ \varrho(s) \f$
       *  optionally normalized at the point \f$ m_n \f$ 
       * \f[ \varrho (s, m_n) = \Theta\left(s-s_{threshold}\right) \f] 
       */
      double rho_s 
      ( const double s  , 
        const double mn ) const override 
      { return s <= m_sthreshold ? 0.0 : 1 ; }
      /// get the opening threshold \f$ s_{threshold} \$ for the channel 
      double s_threshold () const override { return m_sthreshold ; }
      // =======================================================================
    public: //  helper methods 
      // =======================================================================
      /** get a value of mass-dependent width  (as function of s )
       *  \f$ \gamma(s) \f$ 
       */
      inline double gamma ( const double s ) const { return m_gamma ( s ) ; }
      // =======================================================================
    public: //  helper methods 
      // =======================================================================
      /// unique tag/label  
      std::size_t tag       () const override ;
      /// describe the channel 
      std::string describe  () const override { return m_description ; }
      // =======================================================================
    private :
      // =======================================================================
      /// mass-dependent width 
      Width       m_gamma           ; // mass-dependent width 
      /// s-threshold 
      double      m_sthreshold  ; // s-threshold 
      /// unique tag/lable 
      std::size_t m_tag         ; // unique tag/lable 
      /// description 
      std::string m_description ;
      // ======================================================================
    } ;
    // =========================================================================
    /** @class ChannelQ
     *  Description of the very simple S-wave channel
     *  \f[ \begin{array}{lcl}
     *   N^2(s,m_0)      & = & m_0 \Gamma0 q ( s )  } \\
     *   D (s,m_0)       & = &     \Gamma0 q ( s )  } \\
     *  \varrho (s, m_n) & = & \Theta\left(s-s_{threshold}\right)
     *  \end{array}\,,\f]
     */
    class ChannelQ : public ChannelBW 
    {
    public:
      // =======================================================================
      /** @typedef Width 
       *  the function type for the mass-dependent width
       */
      typedef std::function<double(double)> Width ;
      // =======================================================================
    public : 
      // ======================================================================
      /// full constructor with all functions specified 
      ChannelQ
      ( const double gamma ,
        const double m1    , 
        const double m2    ) ;
      // =======================================================================
      /// clone the channel 
      ChannelQ* clone () const override ;
      // =======================================================================
    public:
      // =======================================================================
      /// get the mass of the 1st daughter 
      double            m1         () const { return m_ps2.m1 ()         ; }
      /// get the mass of the 2nd daughter  
      double            m2         () const { return m_ps2.m2 ()         ; }
      /// phase space function 
      const Ostap::Math::PhaseSpace2& ps2 () const { return m_ps2 ; }
      // =======================================================================
    public: 
      // =======================================================================
      /** squared numerator for the amplitude 
       * \f[ N^2(s,m_0) = m_0 \Gamma0 q(s) \f] 
       */
      double N2
      ( const double s  , 
        const double m0 ) const override 
      { return s <= m_ps2.s_threshold () ? 0.0 : m0 * gamma0 () * m_ps2.q_s ( s ) ; }      
      // ======================================================================
      /** term in the denominator for the amplitide
       * \f[ D (s,m_0) = \Gamma0 q(s) \f] 
       */
      // ======================================================================
      std::complex<double> D    
      ( const double s  , 
        const double m0 ) const override 
      { return s <= m_ps2.s_threshold () ? 0.0 : m0 * gamma0 () * m_ps2.q_s ( s ) ; }      
      // ======================================================================
      /** get the phase space factor  \f$ \varrho(s) \f$
       *  optionally normalized at the point \f$ m_n \f$ 
       * \f[ \varrho (s, m_n) = \Theta\left(s-s_{threshold}\right) \f] 
       */
      double rho_s 
      ( const double s  , 
        const double mn ) const override 
      { return s <= m_ps2.s_threshold () ? 0.0 : 1 ; }      
      /// get the opening threshold \f$ s_{threshold} \$ for the channel 
      double s_threshold () const override { return m_ps2.s_threshold () ; }
      // =======================================================================
    public: //  helper methods 
      // =======================================================================
      /// unique tag/label  
      std::size_t tag       () const override ;
      // =======================================================================
      /// describe the channel 
      std::string describe  () const override ;
      // =======================================================================
    private :
      // =======================================================================
      /// two  body phase space 
      Ostap::Math::PhaseSpace2  m_ps2 ; // two  body phase space 
      // =======================================================================
    } ;

    // =========================================================================
    /** @class Channel
     *  Simple definition for the open decay channel, \f$ m_0 > m_1 + m_2 \f$) 
     *  @see Ostap::Math::ChannelBW 
     *  \f[ \begin{array}{ncl}
     *  N & = &  \phantom{i} g \left( \frac{q}/{q_s} \right)^L \frac{F_L(q)}{F_L(q_0)} \\ 
     *  D & = &  i m_0 \Gamma_0 \frac{\varrho(m)}{\varrho(m_0)}  
     *  \left( \frac{q}{q_s}\right)^{2L} \frac{F^2_L(q)}{F^2_L(q_0)}\,,  
     *  \end{array} \f]
     * where
     *  - \f$ \varrho\f$  is two-body phase space
     *  - \f$ F_Lf\f$     is a phenomenological formfactor, e.g. Blatt-Weisskopf factors 
     *  
     *  - the (partial) width constant 
     *  - masses of daughter particles
     *  - orbital momentum 
     *  - form factor
     */ 
    class Channel : public ChannelCW 
    {
    public:
      // ======================================================================
      /** constructor from all parameters and *NO* formfactor
       *  @param gamma the width 
       *  @param m1    the mass of the 1st daughter
       *  @param m2    the mass of the 2nd daughter
       *  @param L     the orbital momentum
       */
      Channel
      ( const double                  gamma = 0.150 , 
        const double                  m1    = 0.139 , 
        const double                  m2    = 0.139 , 
        const unsigned short          L     = 0     ) ;
      // ======================================================================
      /** constructor from all parameters and Jackson's formfactor 
       *  @param gamma the width 
       *  @param m1    the mass of the 1st daughter
       *  @param m2    the mass of the 2nd daughter
       *  @param L     the orbital momentum
       *  @param r     the Jackson's formfactor 
       */
      Channel
      ( const double                  gamma         , 
        const double                  m1            , 
        const double                  m2            , 
        const unsigned short          L             ,
        const FormFactors::JacksonRho r             ) ;
      // ======================================================================
      /** constructor from all parameters and generic formfactor 
       *  @param gamma the width 
       *  @param m1    the mass of the 1st daughter
       *  @param m2    the mass of the 2nd daughter
       *  @param L     the orbital momentum
       *  @param f     the formfactor 
       */
      Channel
      ( const double                  gamma , 
        const double                  m1    , 
        const double                  m2    , 
        const unsigned short          L     ,
        const FormFactor&             f     ) ;
      /// copy constructor 
      Channel ( const Channel&  right )           ;
      /// move constructor 
      Channel (       Channel&& right ) = default ;
      // ======================================================================
      /// clone the channel 
      Channel* clone () const override ;
      // ======================================================================
    public: // trivial accessors  
      // ======================================================================
      /// get the orbital momentum 
      inline unsigned short    L          () const { return m_L                 ; }        
      /// get the formfactor 
      inline const FormFactor* formfactor () const { return m_formfactor.get () ; }
      // ======================================================================
    public: // two main methods
      // ======================================================================
      /** the first main method: numerator
       *  \f[ N^2(s,m_0) = m_0 \Gamma_0 
       *      \left( \frac{q}{q_0}\right)^{2L}        
       *      \frac{F^2_L(q)}{F^2(q_0)} \f] 
       */
      double N2   
      ( const double s  , 
        const double m0 ) const override ;
      // ======================================================================
      /** the second main method: term to the denominator 
       *  \f[ D( s , m_0 ) =  m_0 \Gamma_0        \
       *     \frac{\varrho(m^2)}{\varrho(m_0^2)} 
       *     \left( \frac{q}{q_0}\right)^{2L}        
       *     \frac{F^2_L(q)}{F^2(q_0)} \f] 
       */
      std::complex<double> D
      ( const double s  , 
        const double m0 ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      // unique tag
      std::size_t tag      () const override ; // unique tag
      // ======================================================================
      /// describe the channel 
      std::string describe () const override ;
      // ======================================================================
    private :
      // ======================================================================
      /// orbital momentum 
      unsigned short                           m_L   { 0 } ; 
      /// the formfactor
      std::unique_ptr<Ostap::Math::FormFactor> m_formfactor { nullptr } ;
      // ======================================================================
    };
    // =========================================================================
    /** @class Channel0
     *  Simple definition for the decay channel, 
     *  that can be also applicable for \f$ m_0 < m_1 + m_2 \f$ 
     *
     *  @see Ostap::Math::ChannelBW 
     *  @see Ostap::Math::Channe 
     * 
     *  \f[ \begin{array}{ncl}
     *  N^2(s,m_0) & = &  \phantom{i} g \left( \frac{q}/{q_s} \right)^L F_L(q) \\ 
     *  D  (s,m_0) & = &  i\Gamma_0 \varrho(m)  
     *  \left( \frac{q}{q_s}\right)^{2L} F^2_L(q)\,,  
     *  \end{array} \f]
     * where
     *  - \f$ \varrho\f$ is two-body phase space (can be complex!)
     *  - \f$ F_Lf\f$    is a phenomenological formfactor, e.g. Blatt-Weisskopf factors 
     *  - \f$ q_s\f$     is a momentum scale \f$q_s>0\f$
     *  
     *  If  \f$q_s\f$ is not specified the formulae are:
     *  \f[ \begin{array}{ncl}
     *  N^2 (s,m_0)& = &  \phantom{i} g \left( q \right)^L F_L(q) \\ 
     *  D (s,m_0) & = &  i\Gamma_0 \varrho(m)  
     *  \left( q\right)^{2L} F^2_L(q)\,, \end{array} \f] 
     *
     *  - the squared coupling constant  
     *  - masses of daughter particles
     *  - orbital momentum 
     *  - form factor
     */ 
    class Channel0 : public Channel
    {
    public:
      // ======================================================================
      /** constructor from all parameters and *NO* formfactor
       *  @param gamma the width 
       *  @param m1    the mass of the 1st daughter
       *  @param m2    the mass of the 2nd daughter
       *  @param L     the oribital momentum
       *  @param qs    the scale momentum 
       */
      Channel0
      ( const double                  gamma = 0.150 , 
        const double                  m1    = 0.139 , 
        const double                  m2    = 0.139 , 
        const unsigned short          L     = 0     , 
        const double                  qs    = 0     ) ;
      // ======================================================================
      /** constructor from all parameters and Jackson's formfactor 
       *  @param gamma the width 
       *  @param m1    the mass of the 1st daughter
       *  @param m2    the mass of the 2nd daughter
       *  @param L     the oribital momentum
       *  @param r     the Jackson's formfactor 
       *  @param qs    the scale momentum 
       */
      Channel0
      ( const double                  gamma         , 
        const double                  m1            , 
        const double                  m2            , 
        const unsigned short          L             ,
        const FormFactors::JacksonRho r             ,
        const double                  qs    = 0     ) ;
      // ======================================================================
      /** constructor from all parameters and generic formfactor 
       *  @param gamma the width 
       *  @param m1    the mass of the 1st daughter
       *  @param m2    the mass of the 2nd daughter
       *  @param L     the oribital momentum
       *  @param f     the formfactor 
       *  @param qs    the scale momentum 
       */
      Channel0
      ( const double                  gamma         , 
        const double                  m1            , 
        const double                  m2            , 
        const unsigned short          L             ,
        const FormFactor&             f             ,
        const double                  qs    = 0     ) ;
      /// copy constructor 
      Channel0 ( const Channel0&  right )           ;
      /// move constructor 
      Channel0 (       Channel0&& right ) = default ;
      /// clone the channel 
      Channel0* clone () const override ;
      // ======================================================================
    public: // trivial accessors  
      // ======================================================================
      /// get the momentum 
      double qs() const { return m_qs ; }
      // ======================================================================
    public: // two main methods
      // ======================================================================
      /// the first main method: numerator
      double               N2
      ( const double s  , 
        const double m0 ) const override ;
      // ======================================================================
      /// the second main method: term to the denominator 
      std::complex<double> D 
      ( const double s  , 
        const double m0 ) const override ;
      // =======================================================================
      /** get the phase space factor  \f$ \varrho(s) \f$
       *  - optionally normalized at the point \f$ m_n \f$ 
       *  - optionally nomalized  at \f$ m_q = m(q_s) \f$
       */
      double rho_s 
      ( const double s  , 
        const double mn ) const override  ;
      // =======================================================================
    public:
      // ======================================================================
      // unique tag
      std::size_t tag() const override ; // unique tag
      // ======================================================================
    public: // print it 
      // ======================================================================
      ///  describe the channel 
      std::string describe () const override ;
      // ======================================================================
    private :
      // ======================================================================
      /// momentum   scale 
      double m_qs { 0 } ; // momentum  scale
      // ======================================================================
    };
    
    // =========================================================================
    /** @class ChannelGRL
     *  Description of the channel 
     *  \f[ \begin{array}{lcl}
     *   N^2(s,m_0)      & = & m_0 \Gamma0 f_{N^2}(s,m_0)\ \
     *   D (s,m_0)       & = & m_0 \Gamma0 ( f_{D}(s) + i ( f_L(s) - f_L(m^2_0 ) ) \\
     *  \varrho (s, m_n) & = & \Theta\left(s-s_{threshold}\right) f_{\varrho}(s,m_n)
     *  \end{array}\,,\f]
     *  where \f$ f_{N^2}]\f$,  \f$ f_{D}]\f$, \f$ f_L \f$  and 
     *  \f$ f_{\varrho}]\f$ are provdied externally
     *  - Integresting special case is when   
     *    \f$ f_L (s,m_0^2) \f$ and \f$ f_D \f$ as real and imaginary parts 
     *   of the amplitude related via the dispersion relation  
     *  with the single subtraction 
     *  \f[ f_L(s) = - \frac{s}{\pi} 
     *  \int \frac{f_D(s^\prime d s^{\prime}}{s^\prime(s^\prime -s ) } \f]
     */
    class ChannelGLR : public ChannelBW 
    {
    public:
      // ======================================================================
      template <class FUNCTION1,
                class FUNCTION2,
                class FUNCTION3,
                class FUNCTION4>
      ChannelGLR
      ( const double       gamma            ,
        FUNCTION1          fN2              ,
        FUNCTION2          fD               ,
        FUNCTION3          fL               ,
        FUNCTION4          fRho             ,
        const double       s0               ,
        const std::string& description = "" , 
        const std::size_t  tag         = 0  )
        : ChannelBW ( gamma )
        , m_fN2         ( fN2         )
        , m_fD          ( fD          )
        , m_fL          ( fL          )
        , m_fRho        ( fRho        )
        , m_sthreshold  ( s0          )
        , m_tag         ( tag         ) 
        , m_description ( description )
      {}
      // ======================================================================
      template <class FUNCTION1,
                class FUNCTION2,
                class FUNCTION3>
      ChannelGLR
      ( const double       gamma            ,
        FUNCTION1          fN2              ,
        FUNCTION2          fD               ,
        FUNCTION3          fL               ,
        const double       s0               ,
        const std::string& description = "" , 
        const std::size_t  tag         = 0  )
        : ChannelBW ( gamma )
        , m_fN2         ( fN2         )
        , m_fD          ( fD          )
        , m_fL          ( fL          )
        , m_fRho        ( [s0] ( const double s ) -> double { return s <= s0 ? 0 : 1 ; } )
        , m_sthreshold  ( s0          )
        , m_tag         ( tag         ) 
        , m_description ( description )
      {}  
      // =======================================================================
      ///  copy constructor
      ChannelGLR ( const ChannelGLR& right ) = default ;
      // =======================================================================
      /// clone method
      ChannelGLR*  clone() const override ; // clone method
      // =======================================================================
    public:
      // =======================================================================
      template <typename... ARGS>
      static inline ChannelGLR
      create
      ( const double gamma ,
        ARGS ...     args  ) { return ChannelGLR ( gamma , args ... ) ; }
      // =======================================================================
    public:
      // =======================================================================
      /** squared  numerator for the amplitude/width  
       * \f[ N^2(s,m_0) = m_0 \Gamma0 f_{N^2}(s)\f]
       */
      double               N2
      ( const double s  , 
        const double m0 ) const override { return m0 * gamma0 () * m_fN2 ( s ) ; }      
      // ======================================================================
      /** term in the denominator for the amplitide
       * \f[ D (s,m_0) = m_0 \Gamma0 ( f_{D}(s) + i( f_L(s) - f_L ( m_0^2 )) \f]
       */
      // ======================================================================
      std::complex<double> D    
      ( const double s  , 
        const double m0 ) const override 
      { return m0 * gamma0 () * std::complex<double> ( m_fD ( s ) ,
                                                       m_fL ( s ) -
                                                       m_fL ( m0 * m0 ) ) ; }    
      // ======================================================================
      /** get the phase space factor  \f$ \varrho(s) \f$
       * \f[ \varrho (s, m_n) = \Theta\left(s-s_{threshold}\right) f_{\varrho}(s)\f] 
       */
      double rho_s 
      ( const double    s  , 
        const double /* mn */ ) const override 
      { return s <= m_sthreshold ? 0.0 : m_fRho ( s ) ; }
      /// get the opening threshold \f$ s_{threshold} \$ for the channel 
      double s_threshold () const override { return m_sthreshold ; }
      // =======================================================================
    public: //  helper methods 
      // =======================================================================
      /// unique tag/label  
      std::size_t tag       () const override ;
      /// describe the channel 
      std::string describe  () const override { return m_description ; }
      // =======================================================================
    private :
      // =======================================================================
      /// function N2
      std::function<double(double)>  m_fN2   ; // function N2
      /// function fD
      std::function<double(double)>  m_fD    ; // function fD 
      /// function fL
      std::function<double(double)>  m_fL    ; // function fL
      /// function fRho
      std::function<double(double)>  m_fRho  ; // function fRho
      /// s-threhold 
      double                         m_sthreshold  ; // s-threhold 
      /// unique tag 
      std::size_t                    m_tag         ; // unique tag 
      /// description 
      std::string                    m_description ; // description 
      // ======================================================================
    } ;
    // =========================================================================
    /** @class ChannelNRL
     *  Description of the channel for non-relativistic Breit-Wigner function 
     *  \f[ \begin{array}{lcl}
     *   N^2(s,m_0)      & = & f_{N^2}(s,m_0)                       \\
     *   D (s,m_0)       & = & \frac{1}{2} \Gamma0 * f_{\Gamma}(s)  \\
     *  \varrho (s, m_n) & = & \Theta\left(s-s_{threshold}\right) f_{\varrho}(s,m_n)
     *  \end{array}\,,\f]
     *  where \f$ f_{N^2}]\f$,  \f$ f_{D}]\f$, \f$ f_L \f$  and 
     *  \f$ f_{\varrho}]\f$ are provdied externally
     *  - Integresting special case is when   
     *    \f$ f_L (s,m_0^2) \f$ and \f$ f_D \f$ as real and imaginary parts 
     *   of the amplitude related via the dispersion relation  
     *  with the single subtraction 
     *  \f[ f_L(s) = - \frac{s}{\pi} 
     *  \int \frac{f_D(s^\prime d s^{\prime}}{s^\prime(s^\prime -s ) } \f]
     */
    class ChannelNRL: public ChannelBW 
    {
    public:
      // ======================================================================
      template <class FUNCTION1,
                class FUNCTION2,
                class FUNCTION3>
      ChannelNRL
      ( const double       gamma            ,
        FUNCTION1          fN2              ,
        FUNCTION2          fGamma           ,
        FUNCTION3          fRho             ,
        const double       s0               ,
        const bool         fake             ,
        const std::string& description = "" , 
        const std::size_t  tag         = 0  )
        : ChannelBW     ( gamma       )
        , m_fN2         ( fN2         )
        , m_fGamma      ( fGamma      )
        , m_fRho        ( fRho        )
        , m_sthreshold  ( s0          )
        , m_fake        ( fake        ) 
        , m_tag         ( tag         ) 
        , m_description ( description )
      {}
      // ======================================================================
      template <class FUNCTION1,
                class FUNCTION2>
      ChannelNRL
      ( const double       gamma            ,
        FUNCTION1          fN2              ,
        FUNCTION2          fGamma           ,
        const double       s0               ,
        const bool         fake             ,
        const std::string& description = "" , 
        const std::size_t  tag         = 0  )
        : ChannelBW     ( gamma       )
        , m_fN2         ( fN2         )
        , m_fGamma      ( fGamma      )
        , m_fRho        ( [] ( const double /* s */ ) -> double { return 1 ; } )
        , m_sthreshold  ( s0          )
        , m_fake        ( fake        ) 
        , m_tag         ( tag         ) 
        , m_description ( description )
      {}
      // ======================================================================
      template <class FUNCTION2>
      ChannelNRL
      ( const double       gamma            ,
        FUNCTION2          fGamma           ,
        const double       s0               ,
        const bool         fake             ,
        const std::string& description = "" , 
        const std::size_t  tag         = 0  )
        : ChannelBW     ( gamma       )
        , m_fN2         ( fGamma      )
        , m_fGamma      ( fGamma      )
        , m_fRho        ( [] ( const double /* s */ ) -> double { return 1 ; } )
        , m_sthreshold  ( s0          )
        , m_fake        ( fake        ) 
        , m_tag         ( tag         ) 
        , m_description ( description )
      {}
      // =======================================================================
      ///  copy constructor
      ChannelNRL ( const ChannelNRL& right ) = default ;
      // =======================================================================
      /// clone method
      ChannelNRL*  clone() const override ; // clone method
      // =======================================================================
    public:
      // =======================================================================
      template <typename... ARGS>
      static inline ChannelNRL
      create
      ( const double gamma ,
        ARGS ...     args  ) { return ChannelNRL ( gamma , args ... ) ; }
      // =======================================================================
    public:
      // =======================================================================
      /** squared  numerator for the amplitude/width  
       * \f[ N^2(s,m_0) = m_0 \Gamma0 f_{N^2}(s)\f]
       */
      double               N2
      ( const double s  , 
        const double m0 ) const override { return m0 * gamma0 () * m_fN2 ( s ) ; }      
      // ======================================================================
      /** term in the denominator for the amplitide
       * \f[ D (s,m_0) = \frac{1}{2}\Gamma_0 f_{\Gamma} (s) + ... \f]
       */
      // ======================================================================
      std::complex<double> D    
      ( const double s  , 
        const double m0 ) const override 
      {
        const double t1 = 0.5 * gamma0() * m_fGamma ( s ) ;
        if ( !m_fake ) { return std::complex<double> ( t1 , 0 ) ; }
        const double m  = std::sqrt ( s ) ;
        const double t2 = ( m - m0 ) - ( m0 * m0 - s ) ;
        return std::complex<double> ( t1 , t2 ) ;
      }    
      // ======================================================================
      /** get the phase space factor  \f$ \varrho(s) \f$
       * \f[ \varrho (s, m_n) = \Theta\left(s-s_{threshold}\right) f_{\varrho}(s)\f] 
       */
      double rho_s 
      ( const double    s  , 
        const double /* mn */ ) const override 
      { return s <= m_sthreshold ? 0.0 : m_fRho ( s ) ; }
      /// get the opening threshold \f$ s_{threshold} \$ for the channel 
      double s_threshold () const override { return m_sthreshold ; }
      // =======================================================================
    public: //  helper methods 
      // =======================================================================
      /// unique tag/label  
      std::size_t tag       () const override ;
      /// describe the channel 
      std::string describe  () const override { return m_description ; }
      // =======================================================================
    private :
      // =======================================================================
      /// function N2
      std::function<double(double)>  m_fN2         ; // function N2
      /// function fGamma
      std::function<double(double)>  m_fGamma      ; // function fGamma 
      /// function fRho
      std::function<double(double)>  m_fRho        ; // function fRho
      /// s-threhold 
      double                         m_sthreshold  ; // s-threhold
      /// add fake term to have real non-relativiic BReit?
      bool                           m_fake        ; // add fake term? 
      /// unique tag 
      std::size_t                    m_tag         ; // unique tag 
      /// description 
      std::string                    m_description ; // description 
      // ======================================================================
    } ;    
    // ========================================================================
    /** @class  ChannelFlatte
     *  Describe Flatte-like channel 
     *  @see Ostap::Math::ChannelBW 
     *  \f[ \begin{array}{ncl}
     *  N^2(s,m_0)& = & m_0 * g * 16\pi \\ 
     *  D  (s,m_0)& = & m_0 * g \frac{2q}{\sqrt{s}}  
     *  \end{array} \f]
     */
    class ChannelFlatte : public ChannelCW
    {
    public :
      // =====================================================================
      /** constructor from all parameters 
       *  @param g     the coupling constant  (dimensionless)
       *  @param m1    the mass of the 1st daughter
       *  @param m2    the mass of the 2nd daughter
       */
      ChannelFlatte 
      ( const double g  = 0.1     , 
        const double m1 = 0.13957 ,   // GeV/c2 
        const double m2 = 0.13957 );  // GeV/c2
      /// clone method 
      ChannelFlatte* clone() const override ;
      // ======================================================================
    public: // define/override base-class abstract methods 
      // ======================================================================
      /** the first main method: numerator
       * \f[ N^2(s,m_0) = m_0 g \f] 
       */
      double               N2
      ( const double s  , 
        const double m0 ) const override ;
      // ======================================================================
      /** the second main method: term to the denominator 
       *  \f[ D(s,m_0) =  m_0 g \frac{2q}{\sqrt{s}} \f], 
       *  @attention it is purely imaginary below threshold!
       */
      std::complex<double> D   
      ( const double s  , 
        const double m0 ) const override ;
      // =======================================================================
    public:
      // ======================================================================
      // unique tag
      std::size_t tag      () const override ; // unique tag
      // ======================================================================
      /// describe the channel 
      std::string describe () const override ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class  ChannelFlatteBugg
     *  Bugg's modification of Flatte channel 
     *  @see D.V. Bugg, "Re-analysis of data on a(0)(1450) and a(0)(980)"
     *           Phys.Rev.D 78 (2008) 074023
     *  @see https://doi.org/10.1103/PhysRevD.78.074023
     *  @see https://arxiv.org/abs/0808.2706
     *
     *  Describe Flatte-like channel  for \f$ f_0(980)\f$ 
     *  @see Ostap::Math::ChannelFlatte 
     *  @see Ostap::Math::ChannelBW 
     *  \f[ \begin{array}{ncl}
     *  N^2(s,m_0)& = & m_0 * g * 16\pi ) \\ 
     *  D  (s,m_0)& = & m_0 * g \varrho  F^2(-alpha  k_{KK}^2)
     *  \end{array} \f]
     *  where 
     *  \f$ \varrho(s) = f_c \frac{2q_c}{s} + f_n \frac{2q_n}{s} \f$ 
     */
    class ChannelFlatteBugg : public ChannelFlatte 
    {
    public :
      // =====================================================================
      /** constructor from all parameters 
       *  @param g        the coupling constant  (dimensionless)
       *  @param mcharged the mass of the (charged) daughters 
       *  @param mneutral the mass of the neutral daughters 
       *  @param mK       the mass of the charged kaon 
       *  @param alpha    formfactor 
       *  @param fc       the first  isospin factor  
       *  @param fn       the second isospin factor  
       */
      ChannelFlatteBugg 
      ( const double g        = 0.1     , 
        const double mcharged = 0.13957 ,   // GeV/c2  
        const double mneutral = 0.13498 ,   // GeV/c2 
        const double mK       = 0.49368 ,   // GeV/c2 
        const double alpha    = 2.0     ,   // GeV^{-2} form-factor 
        const double fc       = 2.0/3.0 ,   // isospin factor  
        const double fn       = 1.0/3.0 ) ; // isospin factor  
      /// clone method 
      ChannelFlatteBugg* clone() const override ;
      // ======================================================================
    public: // define/override base-class abstract methods 
      // ======================================================================
      /** the second main method: term to the denominator 
       *  \f[ D(s,m_0) =  m_0 g \frac{2q}{\sqrt{s}} \f], 
       *  @attention it is purely imaginary below threshold!
       */
      std::complex<double> D   
      ( const double s  , 
        const double m0 ) const override ;
      // =======================================================================
    public:
      // =====================================================================
      /// mass of charged mode 
      double mcharged () const { return m1()        ; }
      /// mass of neutral mode 
      double mneutral () const { return m_ps2n.m1() ; }
      /// mass of kaon 
      double mK       () const { return m_ps2k.m1() ; }
      /// formfactor 
      double alpha    () const { return m_alpha     ; }
      /// isospin factor for charged mode 
      double fc       () const { return m_fc        ; }
      /// isospin factor for neuyral mode 
      double fn       () const { return m_fn        ; }
      // =======================================================================
    public:
      // ======================================================================
      // unique tag
      std::size_t tag      () const override ; // unique tag
      // ======================================================================
      /// describe the channel 
      std::string describe () const override ;
      // ======================================================================
    protected : 
      // ======================================================================
      /// formfactor 
      double                    m_alpha ; // formfactor 
      /// isospin factor for charged mode 
      double                    m_fc    ; // isospin factor for charged mode 
      /// isospin factor for neutral mode 
      double                    m_fn    ; // isospin factor for neutral mode 
      /// two body phase space for neutral mode
      Ostap::Math::PhaseSpace2  m_ps2n  ; // two body phase space 
      /// two body phase space for kaon mode 
      Ostap::Math::PhaseSpace2  m_ps2k  ; // two body phase space 
      // ======================================================================
    } ;

    // ========================================================================
    /** @class BW
     *
     *  J.D.Jackson,
     *  "Remarks on the Phenomenological Analysis of Resonances",
     *  In Nuovo Cimento, Vol. XXXIV, N.6
     *
     *  @see http://www.springerlink.com/content/q773737260425652/
     *  
     *  @see also http://pdg.lbl.gov/2019/reviews/rpp2018-rev-resonances.pdf
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2011-11-30
     */
    class BW
    {
    public:
      // ======================================================================
      /// constructor from all parameters
      BW ( const double     m0        , 
           const ChannelBW& channel   , 
           const double     scale = 1 ) ;     
      ///    copy constructor
      BW ( const BW&  bw ) ;
      /// move construtor
      BW (       BW&&    ) = default ;
      /// virtual destructor 
      virtual ~BW() ;
      /// clone it 
      virtual BW* clone() const = 0 ;
      // ======================================================================
    protected :
      // ======================================================================
      /// default (empty) constructor 
      BW 
      ( const double m0    = 1 , 
        const double scale = 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /** calculate the Breit-Wigner shape
       *  \f$\frac{1}{\pi}\frac{\omega\Gamma(\omega)}
       *   { (\omega_0^2-\omega^2)^2-\omega_0^2\Gamma^2(\omega)-}\f$
       */
      virtual double operator() ( const double m ) const 
      { return breit_wigner ( m ) ; }
      // ======================================================================
    public:  // amplitude 
      // ======================================================================
      /** Get Breit-Wigner amplitude
       *  \f[ A(m) = \frac{1}{m^2_0 - m^2 - \sum_a D_a ( m^2 ) } \f]
       */
      virtual std::complex<double> amplitude ( const double m ) const ;
      // ======================================================================      
      /** Get Breit-Wigner lineshape in channel \f$ a\f$ : 
       *  \f[ F_a(m) = 2m \varrho(s) N^2_a(s,m_0) 
       *   \frac{\Gamma_{tot}}{\Gamma_{0,a}} \left| \mathcal{A} \right|^2 \f] 
       *  @param m the mass point 
       */
      double breit_wigner ( const double m ) const 
      { return m <= threshold() ? 0.0 : breit_wigner ( m , amplitude ( m ) ) ; }
      // ======================================================================
      /** Get Breit-Wigner lineshape in channel \f$ a\f$ : 
       *  \f[ F_a(m) = 2m \varrho(s) N^2_a(s,m_0) 
       *    \frac{\Gamma_{tot}}{\Gamma_{0,a}} \left| \mathcal{A}  \right|^2 \f] 
       *  @param m the mass point 
       *  @param A the amplitide at this point 
       */
      double breit_wigner 
      ( const double                m , 
        const std::complex<double>& A ) const ;
      // ======================================================================
      /// get factor \f$ N^2(s,m_0^2)\f$ from the main channel
      double N2    ( const double s ) const 
      { return s <= s_threshold () ? 0.0 : channel() -> N2 ( s , m0 () ) ; }      
      /// get factor \f$ \varrho(s,m_n^2) \f$ from the main channel 
      double rho_s ( const double s ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// pole position 
      double         m0     () const { return m_m0    ; }
      /// pole position 
      double         mass   () const { return   m0 () ; }
      /// pole position 
      double         peak   () const { return   m0 () ; }
      /** The sum of "gamma" for each channels
       *  - If all channnels are open channels and for 
       *  each channel \f$ \Gamma_i \f$ represents 
       *  the partial width for this channel, 
       *   the result corrrsponds to  a total width 
       *   of the Breit-Wigner
       */
      double         gamma  () const ;
      /// get the scale factor 
      double         scale  () const { return m_scale ; }
      // =====================================================================
    public: // number of channels 
      // ======================================================================
      /** get the decay channel with index \f$ i\f$
       *  - index <code>0</code> corresponds to "the main" dchannel
       *  @param i channel index 
       *  @return the channel with given index, or <code>nullptr</code> otherwise 
       */
      const ChannelBW*    channel     ( const unsigned short i = 0 ) const
      { return i < m_channels.size() ? m_channels[i].get() : nullptr ; }
      // ======================================================================
      /// get number of channels 
      inline unsigned int nChannels   () const { return m_channels.size()  ; }
      // ======================================================================
      /// get the threshold value (cached in constructor)
      double              threshold   () const { return m_threshold ; }
      /// get the threshold value (cached in constructor)
      double              s_threshold () const { return m_threshold * m_threshold ; }
      // ======================================================================
    public:
      // ======================================================================
      /// set pole position 
      bool setM0     ( const double x ) ;
      /// set pole position 
      bool setMass   ( const double x ) { return setM0     ( x ) ; }
      /// set pole position 
      bool setPeak   ( const double x ) { return setM0     ( x ) ; }
      /// set total width at pole 
      bool setGamma  ( const double x ) ;
      /// set scale factor 
      bool setScale ( const double value ) ;
      // ======================================================================
    public: /// integration 
      // ======================================================================
      /// get the integral
      double integral () const ;
      /// get the integral between low and high limits
      double integral 
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the partial gamma for the certain channel 
      double gamma    ( const unsigned short i ) const 
      { return i < nChannels() ? m_channels[i] -> gamma0()            : 0.0   ; }
      /// set the partial gamma for the certain decay
      bool   setGamma ( const unsigned short i , const double value ) 
      { return i < nChannels() ? m_channels[i] -> setGamma0 ( value ) : false ; }
      // ======================================================================
    public:
      // ======================================================================
      /// unique tag/label 
      virtual std::size_t tag() const ;
      // ======================================================================
    protected : 
      // ======================================================================
      /// add one more channel 
      void add ( const ChannelBW& ) ;
      /// add several channels 
      template <typename ...CHANNELS>
      void add
      ( const ChannelBW&    channel  , 
        const CHANNELS& ... channels ) 
      {
        this -> add ( channel     ) ;
        this -> add ( channels... ) ;
      }
      // ======================================================================
    private:
      // ======================================================================
      /// the mass
      double m_m0        { 1 } ; // the mass
      /// the threshold 
      double m_threshold { 0 } ; // the threshold 
      /// additional scale factor 
      double m_scale     { 1 } ; // additional scale factor
      // ======================================================================
    protected:
      // ======================================================================
      /// the channel(s) 
      mutable std::vector<std::unique_ptr<ChannelBW> > m_channels ; // the channel(s)
      // ======================================================================
    private : // integration workspace 
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace m_workspace ;    // integration workspace
      // ======================================================================
    } ;
    // ========================================================================
    /** @class BreitWigner
     *
     *  J.D.Jackson,
     *  "Remarks on the Phenomenological Analysis of Resonances",
     *  In Nuovo Cimento, Vol. XXXIV, N.6
     *  
     *  http://www.springerlink.com/content/q773737260425652/
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2011-11-30
     */
    class  BreitWigner : public BW
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from all parameters
       *  @param m0   the pole position 
       *  @param gam0 the nominal width at the pole 
       *  @param m1   the mass of the 1st daughter particle 
       *  @param m2   the mass of the 2nd daughter particle 
       *  @param L    the orbital momentum 
       */
      BreitWigner
      ( const double         m0     = 0.770 ,
        const double         gam0   = 0.150 ,
        const double         m1     = 0.139 ,
        const double         m2     = 0.139 ,
        const unsigned short L      = 0     , 
        const double         scale  = 1     ) ;
      // ======================================================================
      /** constructor from all parameters
       *  @param m0   the pole position 
       *  @param gam0 the nominal width at the pole 
       *  @param m1   the mass of the 1st daughter particle 
       *  @param m2   the mass of the 2nd daughter particle 
       *  @param L    the orbital momentum 
       *  @param F    the Jackson's formfactor 
       */
      // ======================================================================
      BreitWigner
      ( const double         m0         ,
        const double         gam0       ,
        const double         m1         ,
        const double         m2         ,
        const unsigned short L          ,
        const FormFactors::JacksonRho F , 
        const double         scale  = 1 ) ;
      // ======================================================================
      /** constructor from all parameters
       *  @param m0   the pole position 
       *  @param gam0 the nominal width at the pole 
       *  @param m1   the mass of the 1st daughter particle 
       *  @param m2   the mass of the 2nd daughter particle 
       *  @param L    the orbital momentum 
       *  @param F    the formfactor 
       */
      BreitWigner
      ( const double         m0         ,
        const double         gam0       ,
        const double         m1         ,
        const double         m2         ,
        const unsigned short L          ,
        const FormFactor&    F          ,
        const double         scale  = 1 ) ;
      // ======================================================================
      /// constructor from the channel 
      BreitWigner
      ( const double         m0          ,
        const ChannelBW&     channel     , 
        const double         scale   = 1 ) ;
      /// copy constructor
      BreitWigner ( const BreitWigner&  bw ) ;
      /// move constructor
      BreitWigner (       BreitWigner&& bw ) = default ;
      // ======================================================================
      /// clone it 
      BreitWigner* clone() const override ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Rho0
     *  \f$ \rho^{0} \rightarrow \pi^+ \pi^- \f$
     *  J.D.Jackson,
     *  "Remarks on the Phenomenological Analysis of Resonances",
     *   Nuovo Cimento, Vol. XXXIV, N.6
     *  @see Ostap::Math::BreitWigner::Jackson_A7
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class  Rho0 : public Ostap::Math::BreitWigner
    {
    public:
      // ======================================================================
      // constructor from all parameters
      Rho0
      ( const double m0       = 770   ,     // MeV
        const double gam0     = 150   ,     // MeV
        const double pi_mass  = 139.6 ,     // MeV 
        const double scale    = 1     ) ;
      /// destructor
      virtual ~Rho0 () ;
      // ======================================================================
      /// clone method 
      Rho0* clone()   const override { return new Rho0 ( *this ) ; }
      // ======================================================================
    public:
      // ======================================================================
      double  m1 () const { return m_m1 ; }
      double  m2 () const { return m_m1 ; }      
      // ======================================================================
    private: 
      // ======================================================================
      double m_m1 ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Kstar0
     *  \f$ K^{*0} \rightarrow K^+ \pi^- \f$
     *  J.D.Jackson,
     *  "Remarks on the Phenomenological Analysis of Resonances",
     *  In Nuovo Cimento, Vol. XXXIV, N.6
     *  @see Ostap::Math::BreitWigner::Jackson_A2
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2014-04-27
     */
    class  Kstar0 : public Ostap::Math::BreitWigner
    {
    public:
      // ======================================================================
      /// constructor from all parameters
      Kstar0
      ( const double m0       = 770   ,     // MeV
        const double gam0     = 150   ,     // MeV
        const double k_mass   = 493.7 ,     // MeV
        const double pi_mass  = 139.6 ,     // MeV
        const double scale    = 1     ) ;
      /// destructor
      virtual ~Kstar0 () ;
      // ======================================================================
      /// clone method 
      Kstar0* clone()   const override { return new Kstar0 ( *this ) ; }
      // ======================================================================
    public:
      // ======================================================================
      double  m1 () const { return m_m1 ; }
      double  m2 () const { return m_m2 ; }      
      // ======================================================================
    private: 
      // ======================================================================
      double m_m1 ;
      double m_m2 ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Phi0
     *  \f$ \phi \rightarrow K^+ K^- \f$
     *  J.D.Jackson,
     *  "Remarks on the Phenomenological Analysis of Resonances",
     *  In Nuovo Cimento, Vol. XXXIV, N.6
     *  @see Ostap::Math::BreitWigner::Jackson_A2
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2014-04-27
     */
    class  Phi0 : public Ostap::Math::BreitWigner
    {
    public:
      // ======================================================================
      /// constructor from all parameters
      Phi0
      ( const double m0       = 1019.5 ,     // MeV
        const double gam0     =    4.3 ,     // MeV
        const double k_mass   =  493.7 ,     // MeV
        const double scale    = 1      ) ;
      /// destructor
      virtual ~Phi0 () ;
      // ======================================================================
      /// clone method 
      Phi0* clone()   const override { return new Phi0 ( *this ) ; }
      // ======================================================================
    public:
      // ======================================================================
      double  m1 () const { return m_m1 ; }
      double  m2 () const { return m_m1 ; }      
      // ======================================================================
    private: 
      // ======================================================================
      double m_m1 ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class BreitWignerMC
     *  function to describe Breit-Wigner signal with several channels,
     *  including Flatte's behaviour 
     *  @see http://pdg.lbl.gov/2019/reviews/rpp2018-rev-resonances.pdf
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2018-11-30
     */
    class  BreitWignerMC : public BW 
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from single channel
       *  @param m0   the pole position 
       *  @param c1   the 1st channel 
       */
      BreitWignerMC 
      ( const double     m0 = 0.770       ,
        const ChannelBW& c1 = Channel ( ) ) ;
      // ======================================================================
      /// template constructor with variadic number of channels
      template <typename... CHANNELS>
      BreitWignerMC ( const double       m0       , 
                      const ChannelBW&   c1       , 
                      const CHANNELS&... channels )
        : BreitWignerMC ( m0 , c1 ) 
      { this->add ( channels... ) ; }
      // ======================================================================
      /// copy constructor
      BreitWignerMC ( const BreitWignerMC&  bw ) = default ;
      /// move constructor
      BreitWignerMC (       BreitWignerMC&& bw ) = default ;
      // ======================================================================
      /// clone it 
      BreitWignerMC* clone() const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// add new channel 
      void addChannel ( const ChannelBW& channel ) { this->add ( channel ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// template creator 
      template <typename... CHANNELS>
      static inline BreitWignerMC
      create ( const double       m0       , 
               const ChannelBW&   c1       , 
               const CHANNELS&... channels )
      { return BreitWignerMC ( m0 , c1 , channels... ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Flatte
     *
     *  S.M. Flatte
     *  "Coupled-channel analysis of the \f$\pi\eta\f$ and \f$K\bar{K}\f$
     *  systems near \f$K\bar{K}\f$threshold"
     *  Physics Letters B, Volume 63, Issue 2, 19 July 1976, Pages 224-227
     *
     *  http://www.sciencedirect.com/science/article/pii/0370269376906547
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2011-11-30
     */
    class  Flatte : public BW 
    {
    public:
      // ======================================================================
      /** constructor from all parameters
       *  \f$ f \rightarrow A_1 + A_2\f$
       *  @param m0    the mass
       *  @param m0g1  parameter \f$ m_0\times g_1\f$
       *  @param g2og1 parameter \f$ g2/g_1       \f$
       *  @param mA1   mass of A1
       *  @param mA2   mass of A2
       *  @param mB1   mass of B1
       *  @param mB2   mass of B2
       */
      Flatte 
      ( const double m0    = 980    ,
        const double m0g1  = 165    ,
        const double g2og1 = 4.21   ,   // dimensionless 
        const double mA1   = 139.57 ,
        const double mA2   = 139.57 ,
        const double mB1   = 493.68 ,
        const double mB2   = 493.68 , 
        const double g0    = 0      ,  // the constant width for "other" decays
        const double scale = 1      ) ; 
      /// copy constructor 
      Flatte ( const Flatte&  right ) ;
      /// move constructor 
      Flatte (       Flatte&& right ) = default ;
      // ======================================================================
      /// clone it!
      Flatte* clone() const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// unique tag
      std::size_t tag() const override ; // unique tag
      // ======================================================================
    public:
      // ======================================================================
      /// coupling constant for the main channel 
      double g1    () const { return m_channels [ 0 ] -> gamma0 () ; }
      /// coupling coinstant for the coupled channel 
      double g2    () const { return m_channels [ 1 ] -> gamma0 () ; }
      /// additional constant width for "extra-channels"
      double gam0  () const { return m_channels [ 2 ] -> gamma0 () ; }
      double g0    () const { return gam0 () ; }
      // ======================================================================
    public :  /// derived quantities 
      // ======================================================================
      /// m  * g1 
      double m0g1  () const { return m0 () * g1 () ; }
      /// g2 / g1 
      double g2og1 () const { return g2 () / g1 () ; }
      // ======================================================================
    public :  /// setters 
      // ======================================================================
      bool setG1   ( const double value ) { return m_channels[0] -> setGamma0 ( value ) ; }
      bool setG2   ( const double value ) { return m_channels[1] -> setGamma0 ( value ) ; }
      bool setGam0 ( const double value ) { return m_channels[2] -> setGamma0 ( value ) ; }
      // ======================================================================
    public :
      // ======================================================================
      double mA1 () const { return m_A1 ; }
      double mA2 () const { return m_A2 ; }
      double mB1 () const { return m_B1 ; }
      double mB2 () const { return m_B2 ; }
      // ======================================================================
    private:
      // ======================================================================
      double m_A1 {} ;
      double m_A2 {} ;
      double m_B1 {} ;
      double m_B2 {} ;
      // ======================================================================
    } ;
    // ========================================================================  
    /** @class FlatteBugg 
     *  Bugg's modification of Flatte channel 
     *  @see D.V. Bugg, "Re-analysis of data on a(0)(1450) and a(0)(980)"
     *           Phys.Rev.D 78 (2008) 074023
     *  @see https://doi.org/10.1103/PhysRevD.78.074023
     *  @see https://arxiv.org/abs/0808.2706
     *  Describe Flatte-like channel  for \f$ f_0(980)\f$ 
     *  @see Ostap::Math::ChannelFlatteBugg 
     */
    class  FlatteBugg : public BW 
    {
    public:
      // ======================================================================
      /** constructor from all parameters
       *  \f$ f \rightarrow A_1 + A_2\f$
       *  @param m0      the mass
       *  @param g1      parameter \f$ g_1    \f$
       *  @param g2og1   parameter \f$ g2/g_1 \f$
       *  @param alpha   parameter alpha (formfactor) 
       *  @param mpiplus mass of the charged pion 
       *  @param mpizero mass of the neutral pion 
       *  @param mKplus  mass of the charged kaon 
       *  @param mKzero  mass of the neutral kaon 
       *  @aram  g0      constant with for "other" decays
       */
      FlatteBugg
      ( const double m0       = 0.980   , // GeV/c2
        const double g1       = 0.165   , // GeV/c2 
        const double g2og1    = 4.21    , // dimensionless  
        const double alpha    = 2.0     , // GeV^-2 
        const double mpiplus  = 0.13957 , // pi+ mass in GeV 
        const double mpizero  = 0.13498 , // pi0 mass in GeV 
        const double mKplus   = 0.49368 , // K+ mass in GeV 
        const double mKzero   = 0.49761 , // K0 mass in GeV 
        const double g0       = 0       ,  // the constant width for "other" decays
        const double scale    = 1       ) ; 
      /// copy constructor 
      FlatteBugg ( const FlatteBugg&  right ) ;
      /// move constructor 
      FlatteBugg (       FlatteBugg&& right ) = default ;
      // ======================================================================
      /// clone it!
      FlatteBugg* clone() const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// unique tag
      std::size_t tag() const override ; // unique tag
      // ======================================================================
    public:
      // ======================================================================
      /// coupling constant for the main channel 
      double g1    () const { return m_channels [ 0 ] -> gamma0 () ; }
      /// coupling coinstant for the coupled channel 
      double g2    () const { return m_channels [ 1 ] -> gamma0 () ; }
      /// additional constant width for "extra-channels"
      double gam0  () const { return m_channels [ 2 ] -> gamma0 () ; }
      double g0    () const { return gam0 () ; }
      // ======================================================================
    public :  /// derived quantities 
      // ======================================================================
      /// m  * g1 
      double m0g1  () const { return m0 () * g1 () ; }
      /// g2 / g1 
      double g2og1 () const { return g2 () / g1 () ; }
      // ======================================================================
    public :  /// other quantities 
      // ======================================================================
      double alpha   () const { return m_alpha    ; }
      double mpiplus () const { return m_mpiplus  ; }
      double mpizero () const { return m_mpizero  ; }
      double mKplus  () const { return m_mKplus   ; }
      double mKzero  () const { return m_mKzero   ; }      
      // ======================================================================
    public :  /// setters 
      // ======================================================================
      bool setG1   ( const double value ) { return m_channels[0] -> setGamma0 ( value ) ; }
      bool setG2   ( const double value ) { return m_channels[1] -> setGamma0 ( value ) ; }
      bool setGam0 ( const double value ) { return m_channels[2] -> setGamma0 ( value ) ; }
      // ======================================================================
    private:
      // ======================================================================
      double m_alpha   {} ;
      double m_mpiplus {} ;
      double m_mpizero {} ;
      double m_mKplus  {} ;
      double m_mKzero  {} ;
      // ======================================================================
    } ;

    // ========================================================================  
    /** @namespace Ostap::Math::Jackson
     *   - Jackson's form-factors 
     */
    namespace Jackson
    {
      // ======================================================================
      /** the simplest function: constant
       *  @see Ostap::Math::BreitWigner
       *  @see Ostap::Math::BreitWigner::rho_fun
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2011-11-30
       */      
      double jackson_0 
      ( double /* m  */ ,
        double /* m0 */ ,
        double /* m1 */ ,
        double /* m2 */ ) ;
      // ======================================================================
      /** the simple function for \f$ 1^- \rightarrow 0^- 0^- \f$, l = 1
       *  \f$\rho(\omega)= \omega^{-1} \f$
       *  @see Ostap::Math::BreitWigner
       *  @see Ostap::Math::BreitWigner::rho_fun
       *  @param m the invariant mass
       *  @return the value of rho-function
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2011-11-30
       */      
      double jackson_A2
      ( double    m     ,
        double /* m0 */ ,
        double /* m1 */ ,
        double /* m2 */ ) ;
      // ======================================================================
      /** the simple function for \f$ 1^- \rightarrow 0^- 1^- \f$, l = 1
       *  \f$\rho(\omega)= \omega \f$
       *  @see Ostap::Math::BreitWigner
       *  @see Ostap::Math::BreitWigner::rho_fun
       *  @param m the invariant mass
       *  @return the value of rho-function
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2011-11-30
       */
      double jackson_A3 
      ( double    m     ,
        double /* m0 */ ,
        double /* m1 */ ,
        double /* m2 */ ) ;
      // ======================================================================
      /** the simple function for
       *  \f$ \frac{3}{2}^+ \rightarrow \frac{1}{2}^+ 0^- \f$, l = 1
       *  \f$\rho(\omega)= \frac{ ( \omega + M )^2 - m^2 }{ \omega^2} \f$
       *  @see Ostap::Math::BreitWigner
       *  @see Ostap::Math::BreitWigner::rho_fun
       *  @param m the invariant mass
       *  @param m1 the invariant mass of the first  (spinor) particle
       *  @param m2 the invariant mass of the secodn (scalar) particle
       *  @return the value of rho-function
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2011-11-30
       */
      double jackson_A4 
      ( double    m     ,
        double /* m0 */ ,
        double    m1    ,
        double    m2     ) ;
      // ======================================================================
      /** the simple function for
       *  \f$ \frac{3}{2}^- \rightarrow \frac{1}{2}^+ 0^- \f$, l = 2
       *  \f$\rho(\omega)= \left[ ( \omega + M )^2 - m^2 \right]^{-1} \f$
       *  @see Ostap::Math::BreitWigner
       *  @see Ostap::Math::BreitWigner::rho_fun
       *  @param m the invariant mass
       *  @param m1 the invariant mass of the first  (spinor) particle
       *  @param m2 the invariant mass of the secodn (scalar) particle
       *  @return the value of rho-function
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2011-11-30
       */      
      double jackson_A5
      ( double    m     ,
        double /* m0 */ ,
        double    m1    ,
        double    m2     ) ;
      // ======================================================================
      /** the simple function for \f$\rho^0 \rightarrow \pi^+ \pi^-\f$ and           
       *  \f$ 1- \rightarrow 0^- 0^- \f$, l = 1
       *  \f$ \rho(\omega)= \left[ q_0^2 + q^2 \right]^{-1} \f$
       *  @see Ostap::Math::BreitWigner
       *  @see Ostap::Math::BreitWigner::rho_fun
       *  @param m the invariant mass
       *  @param m the nominam   mass
       *  @return the value of rho-function
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2011-11-30
       */
      double jackson_A7 
      ( double    m     ,
        double    m0    ,
        double    m1    ,
        double    m2    ) ;
      // ======================================================================
    } //                                               end of namespace Jackson
    // ========================================================================
    /** @class FormFactor
     *  abstract class to implement various formfactors
     */
    class FormFactor
    {
    public :
      // ======================================================================
      /** the only important method: the squared ratio 
       *  of formfactors \f$  \frac{F^2(m)}{F^2(m_0)} \$
       */
      virtual double operator()
      ( const double m  , const double m0 ,
        const double m1 , const double m2 ) const  = 0 ;
      // ======================================================================
      /// virtual destructor
      virtual ~FormFactor () ;
      /// clone method ("virtual constructor" )
      virtual FormFactor* clone    () const = 0 ;
      /// describe the formfactor
      virtual std::string describe () const = 0 ;
      /// some unuque tag/label 
      virtual std::size_t tag      () const = 0 ;      
      // ======================================================================
    } ;
    // ========================================================================
    namespace FormFactors
    {
      // ======================================================================
      /** Formfactor for Breit-Wigner amplitude
       *  parameterization for \f$\rho(\omega)\f$-function from (A.1)
       *  J.D.Jackson,
       *  "Remarks on the Phenomenological Analysis of Resonances",
       *  In Nuovo Cimento, Vol. XXXIV, N.6
       */
      class Jackson : public Ostap::Math::FormFactor
      {
      public:
        // ====================================================================
        /// constructor from enum
        Jackson
        ( const Ostap::Math::FormFactors::JacksonRho rho = 
          Ostap::Math::FormFactors::Jackson_0 ) ;
        /// virtual destructor
        virtual ~Jackson  () ;
        /// clone method ("virtual constructor")
        Jackson* clone() const   override;
        // ====================================================================
        /** the only important method the squared ratio 
         *  of formfactors \f$  \frac{F^2(m)}{F^2(m_0)} \$
         */
        double operator() 
        ( const double m  , 
          const double m0 ,
          const double m1 , 
          const double m2 ) const override;
        // ====================================================================
        /// describe the formfactor 
        std::string describe () const override { return m_what ; }
        // ====================================================================
        // unique tag/label
        std::size_t tag      () const override ;
        // ====================================================================
      public:
        // ====================================================================
        /// get the rho-index
        Ostap::Math::FormFactors::JacksonRho rho () const { return m_rho ; }
        // ====================================================================
      private:
        // ====================================================================
        /// the function itself
        Ostap::Math::FormFactors::JacksonRho m_rho  ;
        /// print it 
        std::string                          m_what ;
        // ====================================================================
      } ;
      // ======================================================================
      /** Blatt-Weisskopf formfactor/barrier factor
       *  actually it is "traslation" of
       *  Blatt-Weiskopf barrier factor into in Jackson's "rho"-function
       */
      class BlattWeisskopf : public Ostap::Math::FormFactor
      {
      public:
        // ====================================================================
        /// orbital momentum
        enum Case {
          Zero  = 0 ,
          One   = 1 ,
          Two   = 2 ,
          Three = 3 ,
          Four  = 4 ,
          Five  = 5
        } ;
        // ====================================================================
      public:
        // ====================================================================
        /// constructor from enum and barrier factor
        BlattWeisskopf ( const Case   L , const double b ) ;
        /// default constructor (needed for  serialization)
        BlattWeisskopf () ;
        /// virtual destructor
        virtual ~BlattWeisskopf () ;
        /// clone method ("virtual constructor")
        BlattWeisskopf* clone() const   override;
        // ====================================================================
        /** the only important method the squared ratio 
         *  of formfactors \f$  \frac{F^2(m)}{F^2(m_0)} \$
         */
        double operator() ( const double m  , const double m0 ,
                            const double m1 , const double m2 ) const override;
        // ====================================================================
        /// describe the formfactor 
        std::string describe () const override { return m_what ; }
        // ====================================================================
        // unique tag/label
        std::size_t tag      () const override ;
        // ====================================================================
      public:
        // ====================================================================
        /// orbital momentum 
        Case    L       () const { return m_L ; }
        /// breakup 
        double  breakup () const { return m_b ; }
        // ====================================================================
      protected:
        // ====================================================================
        /// get the sqyaured ratio of squared barrier factor
        double b ( const double z , const double z0  ) const ;
        // ====================================================================
      private:
        // ====================================================================
        /// orbital momentum
        Case   m_L ; // orbital momentum
        /// Break-up 
        double m_b ; // Break-up 
        // ====================================================================
        std::string  m_what ;
        // ====================================================================
      } ;
      // ======================================================================
      /** Generic formfactor for Breit-Wigner amplitude
       */
      class GenericFF : public Ostap::Math::FormFactor
      {
      public:
        // ====================================================================
        /** @typedef formfactor  
         *  the actual type of the squared ratio of formfactors 
         *   \f[ f ( m , m_0 , m_1  , m_2 ) = 
         *   \frac{    F_L^2( q , q_s ) } { F_L^2 ( q_0 , q_s ) } \f] 
         */
        typedef std::function<double(double,double,double,double)> formfactor ;
        // ====================================================================
      public: 
        // ====================================================================
        /** constructor from the generic object, unique tag and description
         *  @param ff  the formfactor 
         *  @param tag the unique tag 
         *  @param description  description 
         */
        template <class FORMFACTOR>
        GenericFF ( FORMFACTOR         ff                        ,
                    const std::size_t  tag                       , 
                    const std::string& description = "GenericFF" )
          : FormFactor() 
          , m_ff          ( ff  ) 
          , m_tag         ( tag ) 
          , m_description ( description )
        {}
        /// copy constrictor 
        GenericFF ( const GenericFF&  f ) = default ;
        /// move constrictor 
        GenericFF (       GenericFF&& f ) = default ;
        // ====================================================================
        /// clone operaion
        GenericFF* clone () const override ; // clone operaion
        // ====================================================================
      public:
        // ====================================================================
        template <class FORMFACTOR>
        static inline GenericFF 
        create ( FORMFACTOR         ff                        ,
                 const std::size_t  tag                       , 
                 const std::string& description = "GenericFF" )
        { return GenericFF ( ff   ,  tag , description ) ; }
        // ====================================================================
      public:
        // ====================================================================
        /** the only important method the squared ratio 
         *  of formfactors \f$  \frac{F^2(m)}{F^2(m_0)} \$
         */
        double operator() 
        ( const double m  , const double m0 ,
          const double m1 , const double m2 ) const override
        { return m_ff ( m , m0 , m1 , m2 ) ; }
        // ====================================================================
      public:
        // ====================================================================
        /// describe the formfactor
        std::string describe () const override { return m_description ; }
        /// some unuque tag/label 
        std::size_t tag      () const override { return m_tag ; }
        // ====================================================================
      private: 
        // ====================================================================
        /// the formfactor object itself 
        formfactor  m_ff          ; // the formfactor object itself 
        /// unique tag/label
        std::size_t m_tag         ; // unique tag/label
        /// description 
        std::string m_description ;
        // ====================================================================
      } ;
      // ======================================================================
      /** @class NoFormFactor
       *  "No-formfactor" 
       */
      class NoFormFactor : public Ostap::Math::FormFactor
      {
        // ====================================================================
      public:
        // ====================================================================
        /// constructor from enum and barrier factor
        NoFormFactor () ;
        /// virtual destructor
        virtual ~NoFormFactor() ;
        /// clone method ("virtual constructor")
        NoFormFactor* clone() const   override;
        // ====================================================================
        /** the only important method the squared ratio 
         *  of formfactors \f$  \frac{F^2(m)}{F^2(m_0)} \$
         */
        double operator() ( const double /* m  */ , 
                            const double /* m0 */ ,
                            const double /* m1 */ , 
                            const double /* m2 */ ) const override { return 1 ; }
        // ====================================================================
        /// describe the formfactor 
        std::string describe () const override ;
        // ====================================================================
        // unique tag/label
        std::size_t tag      () const override ;
        // ====================================================================
      } ;
      // ======================================================================
    } // end of namespace Ostap:Math::FormFactors
    // ========================================================================


    // ========================================================================
    // VARIOUS BEASTS 
    // ========================================================================


    // // ========================================================================
    // /** @class Bugg
    //  *  parametrisation of sigma-pole for
    //  *  two pion mass distribution
    //  *
    //  *  The parameterization of sigma pole by
    //  *  B.S.Zou and D.V.Bugg, Phys.Rev. D48 (1993) R3948.
    //  *
    //  *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
    //  *  @date 2012-04-01
    //  */
    // class  Bugg
    // {
    // public:
    //   // ======================================================================
    //   /** constructor from all masses and angular momenta
    //    *  @param M  mass of sigma (very different from the pole positon!)
    //    *  @param g2 width parameter g2 (4pi width)
    //    *  @param b1 width parameter b1  (2pi coupling)
    //    *  @param b2 width parameter b2  (2pi coupling)
    //    *  @param a  parameter a (the exponential cut-off)
    //    *  @param s1 width parameter s1  (cut-off for 4pi coupling)
    //    *  @param s2 width parameter s2  (cut-off for 4pi coupling)
    //    *  @param m1 the mass of the first  particle
    //    */
    //   Bugg    ( const double         M  = 0.9264        ,  // GeV
    //             const double         g2 = 0.0024        ,  // GeV
    //             const double         b1 = 0.5848        ,  // GeV
    //             const double         b2 = 1.6663        ,  // GeV-1
    //             const double         a  = 1.082         ,  // GeV^2
    //             const double         s1 = 2.8           ,  // GeV^2
    //             const double         s2 = 3.5           ,
    //             const double         m1 =  139.6 / 1000 ) ; // GeV
    //   /// destructor
    //   ~Bugg () ;
    //   // ======================================================================
    // public:
    //   // ======================================================================
    //   /// calculate the Bugg shape
    //   double operator() ( const double x ) const { return pdf ( x ) ; }
    //   /// calculate the Bugg shape      
    //   double pdf        ( const double x ) const ;
    //   // ======================================================================
    // public:
    //   // ======================================================================
    //   /// get the amlitude  (not normalized!)
    //   std::complex<double> amplitude (  const double x ) const ;
    //   /// get the phase space factor (taking into account L)
    //   double phaseSpace ( const double x ) const { return m_ps  ( x ) ; }
    //   // ======================================================================
    // public:
    //   // ======================================================================
    //   // phase space variables
    //   // ======================================================================
    //   double m1        () const { return m_ps.m1 () ; }
    //   double m2        () const { return m_ps.m2 () ; }
    //   // ======================================================================
    //   double lowEdge   () const { return m_ps. lowEdge() ; }
    //   // ======================================================================
    // public:
    //   // ======================================================================
    //   /// get the running width by Bugg
    //   std::complex<double>
    //   gamma ( const double x ) const ; // get the running width by Bugg
    //   // ======================================================================
    // public:
    //   // ======================================================================
    //   /// adler factor
    //   double               adler       ( const double x ) const ; // adler factor
    //   /// ratio of 2pi-phase spaces
    //   double               rho2_ratio  ( const double x ) const ;
    //   /// ratio of 4pi-phase spaces
    //   std::complex<double> rho4_ratio  ( const double x ) const ;
    //   /// b-factor for 2-pi coupling
    //   double b ( const double x ) const { return  b1 () + x * x * b2 () ; }
    //   // ======================================================================
    // private:
    //   // ======================================================================
    //   /// approximation for  4pi-phase space
    //   std::complex<double> rho4        ( const double x ) const ;
    //   // ======================================================================
    // public:
    //   // ======================================================================
    //   // sigma & Bugg variables
    //   // ======================================================================
    //   /// pole position  
    //   double M     () const  { return m_M       ; }
    //   /// m^2
    //   double M2    () const  { return m_M * m_M ; }
    //   /// pole positon  
    //   double m0    () const  { return   M ()    ; }
    //   /// pole positon      
    //   double mass  () const  { return   M ()    ; }
    //   /// pole positon
    //   double peak  () const  { return   M ()    ; }
    //   // ======================================================================
    //   /// g2 
    //   double g2    () const  { return m_g2   ; }
    //   /// b1 
    //   double b1    () const  { return m_b1   ; }
    //   /// b2 
    //   double b2    () const  { return m_b2   ; }
    //   /// s1 
    //   double s1    () const  { return m_s1   ; }
    //   /// s2 
    //   double s2    () const  { return m_s2   ; }
    //   /// a 
    //   double a     () const  { return m_a    ; }
    //   // ======================================================================
    //   /// set pole position
    //   bool setM    ( const double value  ) ;
    //   /// set pole position
    //   bool setM0   ( const double value  ) { return setM ( value )  ; }
    //   /// set pole position
    //   bool setMass ( const double value  ) { return setM ( value )  ; }
    //   /// set pole position
    //   bool setPeak ( const double value  ) { return setM ( value )  ; }
    //   // ======================================================================
    //   /// g2
    //   bool setG2   ( const double value  ) ;
    //   /// b1 
    //   bool setB1   ( const double value  ) ;
    //   /// b2 
    //   bool setB2   ( const double value  ) ;
    //   /// s1 
    //   bool setS1   ( const double value  ) ;
    //   /// s2 
    //   bool setS2   ( const double value  ) ;
    //   /// a
    //   bool setA    ( const double value  ) ;
    //   // ======================================================================
    // public:
    //   // ======================================================================
    //   /// get the integral
    //   /// double integral () const ;
    //   /// get the integral between low and high limits
    //   double integral ( const double low  ,
    //                     const double high ) const ;
    //   // ======================================================================
    // private:
    //   // ======================================================================
    //   // sigma & Bugg variables
    //   // ======================================================================
    //   /// mass of sigma (very different from the pole positon!)
    //   double m_M  ; // mass of sigma (very different from the pole positon!)
    //   /// width parameter g2 (4pi width)
    //   double m_g2 ; // width parameter g2 (4-p width)
    //   /// width parameter b1  (2pi coupling)
    //   double m_b1 ; // width parameter b1  (2pi coupling)
    //   /// width parameter b2  (2pi coupling)
    //   double m_b2 ; // width parameter b2  (2pi coupling)
    //   /// width parameter s1  (cut-off for 4pi coupling)
    //   double m_s1 ; // width parameter b1  (cut-off for 4pi coupling)
    //   /// width parameter s2  (cut-off for 4pi coupling)
    //   double m_s2 ; // width parameter b2  (cut-off for 4pi coupling)
    //   /// parameter a (the exponential cut-off)
    //   double m_a  ; // parameter a (the exponential cut-off)
    //   // ======================================================================
    // private:
    //   // ======================================================================
    //   /// phase space
    //   Ostap::Math::PhaseSpace2   m_ps         ; // phase space
    //   // ======================================================================
    // private:
    //   /// integration workspace
    //   Ostap::Math::WorkSpace     m_workspace  ;    // integration workspace
    //   // ======================================================================
    // } ;
    // // ========================================================================
    // /** @class Swanson
    //  *  Swanson's parameterization of S-wave cusp
    //  *  @see LHCb-PAPER-2016-019 appendix D
    //  *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
    //  *  @date 2016-06-11
    //  */
    // class  Swanson
    // {
    // public:
    //   // ======================================================================
    //   /// constructor from all parameters (numbers are arbitrary...)
    //   Swanson ( const double         m1     = 0.139 ,   // the first  real particle
    //             const double         m2     = 0.139 ,   // the second real particle
    //             const double         m1_0   = 0.135 ,   // the first  particle for cusp
    //             const double         m2_0   = 0.135 ,   // the second particle for cusp
    //             const double         beta_0 = 0.300 ,   // beta_0 parameter
    //             const unsigned short L      = 0     ) ; // orbital momentum for real particles
    //   /// constructor from all parameters
    //   Swanson ( const double         m1             ,   // the first  real particle
    //             const double         m2             ,   // the second real particle
    //             const double         m1_0           ,   // the first  particle for cusp
    //             const double         m2_0           ,   // the second particle for cusp
    //             const double         beta_0         ,   // beta_0 parameter
    //             const unsigned short L              ,   // orbital momentum for real particles
    //             const Ostap::Math::FormFactors::JacksonRho  r ) ; //  formfactor
    //   /// constructor from all parameters
    //   Swanson ( const double         m1             ,   // the first  real particle
    //             const double         m2             ,   // the second real particle
    //             const double         m1_0           ,   // the first particle for cusp
    //             const double         m2_0           ,   // the second particle for cusp
    //             const double         beta_0         ,   // beta_0 parameter
    //             const unsigned short L              ,   // orbital momentum for real particles
    //             const Ostap::Math::FormFactor&    f ) ; // formfactor
    //   /// constructor from all parameters
    //   Swanson ( const BreitWigner&   bw             ,   // breit-wigner
    //             const double         m1_0           ,   // the first  particle for cusp
    //             const double         m2_0           ,   // the second particle for cusp
    //             const double         beta_0         ) ; // beta_0 parameter
    //   /// copy constructor
    //   Swanson ( const Swanson&  sw ) ;
    //   /// destructor
    //   virtual ~Swanson() ;
    //   // ======================================================================
    // public:
    //   // ======================================================================
    //   /// calculate the Swanson shape
    //   double operator () ( const double x ) const { return swanson ( x ) ; }
    //   /// calculate the Swanson shape
    //   double swanson     ( const double x ) const ;
    //   // ======================================================================
    // public:
    //   // ======================================================================
    //   /// calculate complex amplitude
    //   std::complex<double> amplitude ( const double x ) const ;
    //   // ======================================================================
    // public: // getters
    //   // ======================================================================
    //   /// get beta_0 parameter
    //   double  beta0 () const { return m_beta0 ; }  // get beta_0 parameter
    //   /// mass of the first particle
    //   double  m1    () const { return m_m1    ; }  // mass of first particle
    //   /// mass of the second particle
    //   double  m2    () const { return m_m2    ; }  // mass of the second particle
    //   // ======================================================================
    // public: // derived getters
    //   // ======================================================================
    //   double mmin () const { return m_bw.m1() + m_bw.m2() ; }
    //   double cusp () const { return    m_m1   +    m_m2   ; }
    //   // ======================================================================
    // public: // setters
    //   // ======================================================================
    //   /// set new value for beta_0
    //   bool setBeta0  ( const double value ) ;
    //   /// set new value for beta_0
    //   bool setBeta_0 ( const double value ) { return setBeta0 ( value ) ; }
    //   /// set new valeu for m1
    //   bool setM1_0   ( const double value ) ;
    //   /// set new valeu for m2
    //   bool setM2_0   ( const double value ) ;
    //   // ======================================================================
    // public:
    //   // ======================================================================
    //   /// get the integral between low and high limits
    //   virtual double integral  ( const double low  ,
    //                              const double high ) const ;
    //   // ======================================================================
    // private:
    //   // ======================================================================
    //   /// assignement operator is disabled
    //   Swanson& operator=( const Swanson& ) ; // no assignement
    //   // ======================================================================
    // private:
    //   // ======================================================================
    //   /// use Breit-Wigner to keep parameters of real particles
    //   Ostap::Math::BreitWigner   m_bw ;
    //   /// the mass of the first  particle
    //   double            m_m1         ; // the mass of the first  particle
    //   /// the mass of the second particle
    //   double            m_m2         ; // the mass of the second particle
    //   /// beta0 parameter
    //   double            m_beta0      ; // beta0 parameter
    //   // ======================================================================
    // private:
    //   // ======================================================================
    //   /// integration workspace
    //   Ostap::Math::WorkSpace m_workspace ;    // integration workspace
    //   // ======================================================================
    // } ;
    // // ========================================================================


    // ========================================================================
    // "2-from-3" variants 
    // ========================================================================
    
    // ========================================================================
    /** @class Channel23L 
     *  helper class to represent resonances in (12) system from M->1+2+3 decays
     *  where the orbital momentum between (12)  and (2) is known. 
     *(
     *  - \f$ N_a(s)       \f$ delegates to original channel 
     *  - \f$ D_a(s)       \f$ delegates to the original channel 
     *  - \f$ \varrho_a(s) \f$ phase space 23L 
     */
    class Channel23L: public ChannelBW 
    {
    public:
      // ======================================================================
      /// constructor from the channel and phase-space 
      Channel23L
      ( const ChannelBW&                  ch , 
        const Ostap::Math::PhaseSpace23L& ps ) ;
      /// constructor from the channel and dalitz configuration 
      Channel23L
      ( const ChannelBW&                  ch , 
        const Ostap::Kinematics::Dalitz&  dp , 
        const unsigned short              L2 ) ;
      Channel23L
      ( const ChannelCW&                  ch , 
        const double                      m3 ,  
        const double                      M  ,
        const unsigned short              L2 ) ;
      Channel23L ( const Channel23L& right ) ;
      Channel23L* clone () const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// the first main method: numerator
      double               N2 
      ( const double s  , 
        const double m0 ) const override { return m_channel -> N2 ( s , m0 ) ; }
      // ======================================================================
      /// the second main method: term to the denominator 
      std::complex<double> D   
      ( const double s  , 
        const double m0 ) const override { return m_channel -> D ( s , m0 ) ; }
      // ======================================================================
      /** get the phase space factor  \f$ \varrho \f$
       *  optionally normalized at point \f$ m_n \f$
       */
      double rho_s
      ( const double s  , 
        const double mn ) const override ;
      // ======================================================================
      /// get the opening threshold \f$ s_{threshold} \$ for the channel 
      double s_threshold () const override 
      { 
        const double le = m_ps.lowEdge() ;
        return std::max ( m_channel->s_threshold() , le * le ) ;
      }
      // =======================================================================
    public: //  get thephase space
      // ======================================================================= \
      /// get the phase space factors 
      const Ostap::Math::PhaseSpace23L& ps23L () const { return m_ps ; }
      // =======================================================================
    public: //  helper methods 
      // =======================================================================
      /// unique tag/label  
      std::size_t tag       () const override ;
      /// describe the channel 
      std::string describe  () const override ;
      // ======================================================================
    public: // access to the original channel 
      // ======================================================================
      /// get the original channel 
      inline       ChannelBW* channel ()       { return m_channel.get() ; }
      /// get the original channel 
      inline const ChannelBW* channel () const { return m_channel.get() ; }
      // ======================================================================
    public:
      // ======================================================================
      double    gamma0 () const override 
      { return m_channel->gamma0() ; }
      bool   setGamma0 ( const double value ) override
      { return m_channel->setGamma0 ( value ) ; }
      // ======================================================================
    private:
      // ======================================================================
      /// the original channel 
      std::unique_ptr<ChannelBW> m_channel ;
      /// the phase space
      Ostap::Math::PhaseSpace23L m_ps      ;    // the phase space      
      // ======================================================================
    };
    // ========================================================================
    /** @class ChannelNR3
     *  Describe a non-resonant 3-body decay channel \f$ m \rigthtarrow m_1 m_2 m_3 \f$
     *  The basic functions are : 
     *  \f[\begin{array}{lcc} 
     *   N^2(s,m_0)      & = & m_0 \Gamma_0 \frac{\varrho_3(s)}{\varrho_3(m_0^2)} \\ 
     *   D(s,m_0)        & = & m_0 \Gamma_0 \frac{\varrho_3(s)}{\varrho_3(m_0^2)} \\ 
     *   \varrho (s,m_0) & = & 1 
     *   \end{array}\,,\f]
     *  where \f$ \rho_3(s)\f$ is a threebody phase space      
     */
    class ChannelNR3 : public ChannelBW 
    {
    public: 
      // =====================================================================
      /// constructor from (partial) width and three masses 
      ChannelNR3
      ( const double gamma = 1 , 
        const double m1    = 1 , 
        const double m2    = 1 , 
        const double m3    = 1 ) ;
      // =====================================================================
      /// clone method
      ChannelNR3* clone() const override ;  
      // =====================================================================
    public:
      // =====================================================================
      /**  squared  numerator for the amplitude 
       * \f$ N^2(s,m_0) =  m_0 \Gamma_0 \frac{\varrho_3(s)}{\varrho_3(m_0^2)} \f$ 
       */
      double N2
      ( const double s  , 
        const double m0 ) const override ;      
      // ======================================================================
      /** term in the denominator for the amplitide
       *  \f$ D(s,m_0) = m_0 \Gamma_0 \frac{\varrho_3(s)}{\varrho_3(m_0^2)} \f$
       */
      std::complex<double> D    
      ( const double s  , 
        const double m0 ) const override ;
      // ======================================================================
      /** get the phase space factor  \f$ \varrho(s) \f$
       *  optionally normalized at the point \f$ m_n \f$ 
       */
      double rho_s 
      ( const double    s     , 
        const double /* mn */ ) const override 
      { return s <= m_sthreshold ? 0.0 : m_sthreshold ; }
      /// get the opening threshold \f$ s_{threshold} \$ for the channel 
      double s_threshold () const override { return m_sthreshold ; }
      // =====================================================================
    public: //  helper methods 
      // =======================================================================
      /// unique tag/label  
      std::size_t tag       () const override ;
      /// describe the channel 
      std::string describe  () const override ;
      // ======================================================================
    public:
      // =====================================================================
      double m1 () const { return m_m1 ; }
      double m2 () const { return m_m2 ; }
      double m3 () const { return m_m3 ; }      
      // ======================================================================
    private: 
      // =====================================================================
      /// the first  mass 
      double m_m1         ; // the first mass 
      /// the second mass 
      double m_m2         ; // the second mass 
      /// the third  mass 
      double m_m3         ; // the third mass
      /// the s-threshold 
      double m_sthreshold ; // the s-thrrshold 
      // =====================================================================
    } ;
    // ========================================================================
    /** @class GammaBW3 
     *  Running width/phase-space function 
     *  for 3-body decays 
     *  \f[ \Gamma(s) =  \frac{pi^2}{4s} \frac{1}{(2\pi)^5}
     *  \int\int ds_1 ds_2 \frac{1}{2J_i+1}\sum_i\sum_f
     *   \left|\mathcal{A}\left(s,s_1, s_2\right)\right|^2\f] 
     *  @attention note the power of \f$s\f$ in denumerator! 
     */
    class GammaBW3
    {
    public :
      // ======================================================================
      /** @typedef MatrixElement2
       *  Squared module of amplitude, averaged over initial and 
       *  summed over the final spin states
       *  \f[ \frac{1}{2J_i+1}\sum_i \sum_f 
       *   \left|\mathcal{A}\left(s,s_1,s_2\right)\right|^2 \f] 
       */
      typedef std::function<double(double,double,double)> MatrixElement2 ;
      // ======================================================================      
    public :
      // ======================================================================
      /** constructor from the Dalitz configuration and 
       *  the squared matrix element
       *  @param dalitz Dalizt configriation
       *  @param me2 squared matrix element \f$ M^2 (s,s_1,s_2) \equiv \frac{1}{2J+1}\sum_i \left| \mathcal{A} \right| \f$
       */
      GammaBW3 
      ( const Ostap::Kinematics::Dalitz0& dalitz     , 
        MatrixElement2                    me2        , 
        const std::size_t                 tag    = 0 ,
        const unsigned short               n1    = 0 , 
        const unsigned short               n2    = 0 ) ;
      // =====================================================================
      /// templated constructor 
      template <class ME2> 
      GammaBW3 
      ( const  Ostap::Kinematics::Dalitz0& dalitz  , 
        ME2                                me2     , 
        const std::size_t                  tag     , 
        const unsigned short               n1  = 0 , 
        const unsigned short               n2  = 0 )
        : m_me2    ( me2    ) 
        , m_dalitz ( dalitz ) 
        , m_tag    ( tag    )
        , m_n1     ( n1 ) 
        , m_n2     ( n2 )
      {}
      // ======================================================================      
    public: // the main method 
      // ======================================================================      
      /// the main method 
      double operator ()  ( const double s ) const ;
      // ======================================================================      
    public:
      // ======================================================================
      /** templated creator 
       *  (PyROOT does not like templated constructors) 
       */
      template <class ME2> 
      static inline GammaBW3 
      create 
      ( const Ostap::Kinematics::Dalitz0& dalitz     , 
        ME2                               me2        , 
        const std::size_t                 tag    = 0 ,   
        const unsigned short              n1     = 0 , 
        const unsigned short              n2     = 0 )
      { return GammaBW3 ( dalitz , me2 , tag , n1 , n2 ) ; }
      // ======================================================================      
    public: //  accessors 
      // ======================================================================      
      /// matrix element 
      const MatrixElement2&             me2    () const { return m_me2    ; }
      /// Dalitz configurtaion 
      const Ostap::Kinematics::Dalitz0& dalitz () const { return m_dalitz ; }
      // ======================================================================      
      /// s-threshold 
      inline double s_threshold () const { return m_dalitz.s_min () ; }
      /// tag?
      std::size_t  tag          () const { return m_tag ; }
      // ======================================================================
    private:
      // ======================================================================
      /// Matrix element 
      MatrixElement2             m_me2    {   } ; // Matrix element 
      /// Dalitz configuration 
      Ostap::Kinematics::Dalitz0 m_dalitz {   } ; // Dalitz configuration       
      /// unique tag/key/label (if specified)
      std::size_t                m_tag    {   } ;
      /// useful for decays with narrow resonances 
      unsigned short             m_n1     { 0 } ;
      /// useful for decays with narrow resonances 
      unsigned short             m_n2     { 0 } ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class ChannelDalitz
     *  Describe a three-body decays with given matrix element  (squared absolute 
     *  value of the amplitude averaged over the initial spin states 
     *  The basic functions are : 
     *  \f[\begin{array}{lcc} 
     *   N^2(s,m_0)      & = & m_0 \Gamma_0 \frac{\varrho_D(s)}{\varrho_D(m_0^2)} \\ 
     *   D(s,m_0)        & = & m_0 \Gamma_0 \frac{\varrho_D(s)}{\varrho_D(m_0^2)} \\ 
     *   \varrho (s,m_0) & = & 1 
     *   \end{array}\,,\f]
     *  where \f[ 
     *   \varrho_D(s) \equiv s^{-3/2} \int\int ds_1ds_2 M^2(s,s_1,s_2) \f]
     */
    class ChannelDalitz : public ChannelWidth
    {
    public:
      // =====================================================================
      /// constructor from (partial) width
      ChannelDalitz 
      ( const double                                 gamma  , 
        const Ostap::Kinematics::Dalitz0&            dalitz , 
        Ostap::Math::GammaBW3::MatrixElement2        me2    ,
        const std::size_t                            tag    , 
        const std::string& description = "ChannelDalitz"    ) ;      
      // ======================================================================
      /** templated constructor from (partial) width and matrix element
       *  @param gamma ( partial width) 
       *  @param dalitz Dalitz configuration
       *  @param me2   squared matrix element 
       *  @param tag   unique tag/label 
       *  @param description descrfiption 
       */
      template <class ME2> 
      ChannelDalitz 
      ( const double                                 gamma  , 
        const Ostap::Kinematics::Dalitz0&            dalitz , 
        ME2                                          me2    ,
        const std::size_t                            tag    , 
        const std::string& description = "ChannelDalitz"    ) 
        : ChannelWidth ( gamma , GammaBW3 ( dalitz , me2 , tag ) , 
                         dalitz.s_min () , tag , description )
      {}
      // =====================================================================
      /// clone method 
      // =====================================================================
      ChannelDalitz* clone() const override ;
      // =====================================================================
    public:
      // =====================================================================
      /** templated creator from (partial) width and matrix element
       *  *(PyROOT does not like  templated construictors) 
       *  @param gamma ( partial width) 
       *  @param dalitz Dalitz configuration
       *  @param me2   squared matrix element 
       *  @param tag   unique tag/label 
       *  @param description descrfiption 
       */
      template <class ME2> 
      static inline ChannelDalitz
      create 
      ( const double                                 gamma  , 
        const Ostap::Kinematics::Dalitz0&            dalitz , 
        ME2                                          me2    ,
        const std::size_t                            tag    , 
        const std::string& description = "ChannelDalitz"    ) 
      { return ChannelDalitz ( gamma , dalitz , me2 , tag ) ; }
      // ======================================================================
    };
    // ========================================================================
    /** @class  ChannelGS 
     *  Gounaris-Sakurai parameteriztaion of \f$ \rho \rigtharrow \pi+ \pi^- \$f
     *  @see Gounaris, G.J. and Sakurai, J.J.
     *       "Finite width corrections to the vector meson dominance prediction for rho ---> e+ e-}".
     *        Phys. Rev. Lett 21, (1968) 244 
     *        doi = "10.1103/PhysRevLett.21.244"
     *  @see https://doi.org/10.1103/PhysRevLett.21.244 
     *  @see  Lichard, Peter and Vojik, Martin,
     *       "An Alternative parametrization of the pion form-factor and the mass and width of rho(770)}",
     *        hep-ph/0611163, 2006,
     *  @see https://arxiv.org/abs/hep-ph/0611163
     */
    class ChannelGS : public ChannelBW 
    {
    public :
      // ======================================================================
      /// constructor with gamma and pion mass
      ChannelGS
      ( const double gamma = 150 ,
        const double mpi   = 139 ) ;
      /// clone method
      ChannelGS* clone() const override ;
      // ======================================================================
    public:
      // =====================================================================
      /**  squared  numerator for the amplitude 
       * \f$ N^2(s,m_0) =  m_0 \Gamma_0 \frac{\varrho_3(s)}{\varrho_3(m_0^2)} \f$ 
       */
      double N2
      ( const double s  , 
        const double m0 ) const override ;      
      // ======================================================================
      /** term in the denominator for the amplitide
       *  \f$ D(s,m_0) = m_0 \Gamma_0 \frac{\varrho_3(s)}{\varrho_3(m_0^2)} \f$
       */
      std::complex<double> D    
      ( const double s  , 
        const double m0 ) const override ;
      // ======================================================================
      /** get the phase space factor  \f$ \varrho(s) \f$
       *  optionally normalized at the point \f$ m_n \f$ 
       */
      double rho_s 
      ( const double s  , 
        const double mn ) const override ;
      /// get the opening threshold \f$ s_{threshold} \$ for the channel 
      double s_threshold () const override { return m_sthreshold ; }
      // =====================================================================
    public: //  helper methods 
      // =======================================================================
      /// unique tag/label  
      std::size_t tag       () const override ;
      /// describe the channel 
      std::string describe  () const override ;
      // ======================================================================
    public :
      // ======================================================================
      /// h-function
      double h       ( const double s ) const ;
      /// derivative of h-function
      double h_prime ( const double s ) const ;
      /// mpi
      double mpi     () const { return m_mpi ; }
      // ======================================================================
    private:
      // ======================================================================
      /// pion mass 
      double m_mpi        ; // pion mass
      /// s-threshold: \f$ 4m^2_{\pi} \f$ 
      double m_sthreshold ; // s-threshold: \f$ 4m^2_{\pi} \f$ 
      // ======================================================================
    } ;
    // ========================================================================

    // ========================================================================
    /** @class LASS
     *  The LASS parameterization. 
     *  It describes the 0+ component of the Kpi spectrum ("kappa")
     *  It consists of the K*(1430) resonance together with an
     *  effective range non-resonant component
     * 
     * \f[ \begin{array} {rcl}
     *  \mathcal{A}(m) & = & A_{\mathrm{B}}  +  \\ 
     *                       A_{\mathrm{BW}} {\mathrm{e}}^{i\phi} \\
     *  A_{\mathrm{B}}  & = & \sin \delta \mathrm{e}^{i\delta}   \\ 
     *  \cot \delta     & = & \frac{1}{aq} + \frac{1}{2}bq         \\ 
     &  A_{\mathrm{BW}} & = & \frac{m_0\Gamma_1}{ \left(m^2_0 - m^2 right)
     *                        - im_0 (\Gamma_1 + \Gamma_2 ) }      \\ 
     *  \Gamma_i        & = & q_i \Gamma_{R,i}       \\
     *  \phi            & = & 2\delta \end{array} \f] 
     *
     *  @see D. Aston et al.,  "A Study of K- pi+ Scattering in the Reaction 
     *                         K- p ---> K- pi+ n at 11-GeV/c}",
     *                         Nucl .Phys. B, 296 (1988) 493
     *  @see  DOI: 10.1016/0550-3213(88)90028-4
     *  @see https://doi.org/10.1016/0550-3213(88)90028-4  
     *  @see https://lib-extopc.kek.jp/preprints/PDF/1987/8708/8708374.pdf
     * 
     *  @see P. Estabrooks,"Where and what are the scalar mesons?",
     *                      Phys.Rev.D, 19 ()1979) 2678, 
     *                      doi = "10.1103/PhysRevD.19.2678",
     *  @see https://doi.org/10.1103/PhysRevD.19.2678
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2013-10-05
     */
    class LASS : public BW 
    {
    public: 
      // ======================================================================
      /** constructor from all masses and angular momenta
       *  @param m0 the mass of K*(1450) 
       *  @param g0 the width of  K*(1430)
       *  @param m1 the mass of the first  particle (kaon)
       *  @param m2 the mass of the second particle (pion)
       *  @param m3 the mass of the third  particle (eta')
       *  @param a  the LASS parameter
       *  @param b  the LASS parameter
       *  @param e  the LASS parameter
       */
      LASS
      ( const double m0 = 1429    ,   // K*(1450) mass
        const double g0 =  287    ,   // K*(1430) width
        const double m1 =  493.7  ,   // kaon mass 
        const double m2 =  139.6  ,   // pion mass 
        const double m3 =  957.8  ,   // eta' mass 
        const double a  = 4.03e-3 , 
        const double b  = 1.29e-3 ,
        const double e  = 1.00    ) ; // elasticity 
      // ======================================================================
      /// clone method 
      LASS* clone () const  override ;
      // ======================================================================
    public :
      // ======================================================================
      /** LASS amplitude 
       * \f[ \begin{array} {rcl}
       *  \mathcal{A}(m) & = & A_{\mathrm{B}}  +  \\ 
       *                       A_{\mathrm{BW}} {\mathrm{e}}^{i\phi} \   \
       *  A_{\mathrm{B}}  & = & \sin \deeelta \mathrm{e}^{i\delta}   \\ 
       *  \cot \delta     & = & \frac{1}{aq} + \frac{1}{2}bq         \\ 
       &  A_{\mathrm{BW}} & = & \frac{M_R\Gamma_1}{ \left(\M^2_R - M^2 right)
       *                        - iM_R (\Gamma_1 + \Gamma_2 ) }      \\ 
       *  \Gamma_i        & = & q_i \Gamma_{R,i}       \\
       *  \phi            & = & 2\delta  \end{array} \f] 
       */
      std::complex<double>
      amplitude          ( const double m ) const override ;
      /// evaluate LASS function 
      double operator () ( const double m ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// unique label/tag 
      std::size_t tag () const override ;
      // ======================================================================
    public : // getters 
      // ======================================================================     
      /// a-parameter of LASS function
      double a () const { return m_a ; }
      /// b-parameter of LASS function
      double b () const { return m_b ; }
      /// elasticity 
      double e () const { return m_e ; }
      // ======================================================================
    public :
      // ======================================================================
      double m1 () const { return m_ps2.m1 () ; }
      double m2 () const { return m_ps2.m2 () ; }
      double m3 () const { return m_m3        ; }
      // ======================================================================
    public : // setters 
      // ======================================================================
      /// a 
      bool setA ( const double value ) ;
      /// b
      bool setB ( const double value ) ;
      /// elatisicity 
      bool setE ( const double value ) ;
      // ======================================================================
    private:
      // ======================================================================
      /// a-parametter of LASS function 
      double                     m_a ;  // a-parametter of LASS functuion 
      /// b-parametter of LASS functuion 
      double                     m_b ;  // b-parametter of LASS functuion 
      /// elasticity parameter 
      double                     m_e ;  // elasticity parameter
      /// phase space 
      Ostap::Math::PhaseSpace2   m_ps2 ; // phase space
      /// keep m3 for easy access 
      double                     m_m3 ;
      // ======================================================================
    } ;  
    // ========================================================================
    /** @class BWPS
     *  Breit-Wigner function modulated with some phase-space function
     *  - it can approximate the distorted Breit-Wigner shapes 
     *    from multibody decays 
     *
     *  \f[ f(x) \equiv F_{\mathrm{BW}}(x) \Phi_{l,n}(x)  P_k(x) \f]
     *  - \f$ \Phi_{l,n} \f$  is a phase-space function 
     *  - \f$ P_{k} \f$  is a polynomial 
     * 
     *  The function \f$  F_{\mathrm{BW}}(x) \f$ if defined as 
     *  - for <code>use_rho=true</code> and <code>use_N2=true</code> 
     *  \f$ F_{\mathrm{\BW}}(x) \f$ is a Breit-Wigner lineshape 
     *  - <code>use_rho=true</code> and <code>use_N2=false</code> 
     *  \f$ F_{\mathrm{\BW}}(x) = x \left| \mathcal{A}_{\mathrm{BW}}(x)\right|^2
     *   \varrho ( x^2 ) \f$ i
     *  - <code>use_rho=false</code> and <code>use_N2=true</code> 
     *  \f$ F_{\mathrm{\BW}}(x) = x \left| \mathcal{A}_{\mathrm{BW}}(x)\right|^2
     *   N^2_{\mathrm{BW}} ( x ) \f$ i
     *  - <code>use_rho=false</code> and <code>use_N2=false</code> 
     *  \f$ F_{\mathrm{\BW}}(x) = x \left| \mathcal{A}_{\mathrm{BW}}(x)\right|^2 \f$ 
     *  where \f$ \mathcal{A}_{\mathrm{BW}}(x)  \f$ is a complex 
     *  Breit-Wigner amplitude 
     */
    class BWPS
    {
    public:
      // ======================================================================
      /** constructor from Breit-Wigner, Phase-space and flags 
       *  @param bw Breit-Wigner shape 
       *  @param ps phase-space function 
       *  @param use_rho  use rho-function from Breit-Wigner 
       *  @param use_N2   use N2-function from Breit-Wigner 
       */
      BWPS
      ( const Ostap::Math::BW&            bw      , 
        const Ostap::Math::PhaseSpacePol& ps      , 
        const bool use_rho                = true  , 
        const bool use_N2                 = true  ) ;
      // ======================================================================
      /** constructor from Breit-Wigner, phase-space and flags 
       *  @param bw Breit-Wigner shape 
       *  @param ps phase-space function 
       *  @param use_rho  use rho-function from Breit-Wigner 
       *  @param use_N2   use N2-function from Breit-Wigner 
       */
      BWPS
      ( const Ostap::Math::BW&            bw      , 
        const Ostap::Math::PhaseSpaceNL&  ps      , 
        const bool use_rho                = true  , 
        const bool use_N2                 = true  ) ;
      // ======================================================================
      /// copy constructor 
      BWPS ( const BWPS&  ) ;
      /// move constructor 
      BWPS (       BWPS&& ) = default ;
      // =======================================================================
    public: 
      // ======================================================================
      ///  fictive public default constructor (needed for serialization)
      BWPS() {};
      // =======================================================================
    public: 
      // ======================================================================
      /// evaluate the function 
      double evaluate   ( const double x ) const ;
      /// evaluat ethe function 
      double operator() ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      inline const Ostap::Math::BW& breit_wigner() const { return *m_bw.get() ; }
      inline       Ostap::Math::BW& breit_wigner()       { return *m_bw.get() ; }
      inline const PhaseSpacePol&   phase_space () const { return  m_ps       ; }
      inline       PhaseSpacePol&   phase_space ()       { return  m_ps       ; }      
      // ======================================================================
    public:
      // ======================================================================
      /// get an integral 
      double integral () const ;
      /// get an integral 
      double integral ( const double xmin , const double xmax ) const ;
      // ======================================================================
    public:
      // ======================================================================
      double xmin () const { return std::max ( m_ps.xmin() , m_bw->threshold() ) ; }
      double xmax () const { return            m_ps.xmax() ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get number of polynomial parameters
      std::size_t npars () const { return m_ps.npars () ; }      
      /// set k-parameter
      bool setPar       ( const unsigned short k , const double value )
      { return m_ps.setPar ( k , value ) ; }
      /// set k-parameter
      bool setParameter ( const unsigned short k , const double value )
      { return setPar   ( k , value ) ; }
      /// get the parameter value
      double  par       ( const unsigned short k ) const
      { return m_ps.par ( k ) ; }
      /// get the parameter value
      double  parameter ( const unsigned short k ) const
      { return m_ps.par ( k ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the  gamma for the certain channel 
      double gamma    ( const unsigned short i ) const 
      { return m_bw->gamma ( i ) ; }
      /// set the gamma for the certain decay
      bool   setGamma ( const unsigned short i , const double value ) 
      { return m_bw->setGamma ( i , value ) ; }
      /// get the total gamma 
      double    gamma () const { return m_bw->gamma() ; }
      /// set the total gamma 
      bool   setGamma ( const double value ) 
      { return m_bw->setGamma ( value ) ; }
      /// get number of channels 
      unsigned int nChannels() const { return m_bw->nChannels ()  ; }
      /// pole position 
      double       m0     ()   const { return m_bw->m0()     ; }
      /// set pole position 
      bool      setM0     ( const double x ) { return m_bw->setM0 ( x ) ; }
      // get the amplitude 
      std::complex<double> amplitude ( const double m ) const 
      { return m_bw->amplitude  ( m ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// some unique tag 
      std::size_t tag() const ;
      // ======================================================================
    public:
      // ======================================================================
      bool use_rho () const { return m_rho ; }
      bool use_N2  () const { return m_N2  ; }      
      // ======================================================================
    public:
      // ======================================================================
      /// use rho-factor from BreitWigner ? 
      bool                              m_rho { true    } ;
      /// use N2-factor from BreitWigner ? 
      bool                              m_N2  { true    } ;
      /// Breit-wigner 
      std::unique_ptr<Ostap::Math::BW>  m_bw  { nullptr } ;
      /// Phasespace * pol 
      Ostap::Math::PhaseSpacePol        m_ps  {} ;
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace      m_workspace {} ; // integration workspace
      // ======================================================================
    };
    // ========================================================================
    /** @class BW3L
     *  Breit-Wigner function modulated with \f$ p^{2L+1}\f$`factor
     *  - it can approximate the mass distrbition from 3-body decays   
     *    e.g.  \f$ \eta^{\prime)  \rigtharrow \left(\rho^0 
     *               \rigtharrow \pi^+ \pi^-\right)\gamma \f$~decays
     *    or similar  configurations  
     * 
     *  \f[ f(x) \equiv F_{\mathrm{BW}}(x) p(x|M_0,m_3)^{2L+1} \f]
     *  - \f$ p(x|M,m_3) \f$ is a momentumm of the 3rd particle, \f$P_3\f$ 
     *       in the \f$ P \rightarrow \left( P_{\mathrm{BW}} \rightharrow 
     *      P_1 P_2 \right) P_3 \f$ decay chain
     *  - \f$ M \f$ is a (fixed) mass of "mother" particle \f$P\f$
     *  - \f$ m_1\f$ is a (fixed) mass of 1st particle \f$P_1\f$
     *  - \f$ m_2\f$ is a (fixed) mass of 2nd particle \f$P_2\f$
     *  - \f$ m_3\f$ is a (fixed) mass of 3rd particle \f$P_3\f$
     *  - \f$ x \equiv m_{23} \f$ is a mass intermediate Breit-Wigner particle \f$P_{\mathrm{BW}}\f$
     *  - \f$ L \f$  is an orbital momentum between \f$ P_{\mathrm{BW}}\f$ and \f$ P_3\f$
     * 
     *  It is assumed that  \f$ m_1\f$  and \f$ m_2\f$ parameters 
     *  are in agreement with the Breit-Wigner definition 
     */
    class BW3L
    {
    public:
      // ======================================================================
      /** constructor from Breit-Wigner, Phase-space and flags 
       *  @param bw Breit-Wigner shape 
       *  @param M  mass of the "mother" particle 
       *  @param m1 mass of the 1st      particle 
       *  @param m2 mass of the 2nd      particle 
       *  @param m0 mass of the 3rd      particle 
       *  @param L  the orbital momentum between  
       *  system of 1st and 2nd particles and the 3rd particle
       */
      BW3L
      ( const Ostap::Math::BW& bw , 
        const double           M  ,   
        const double           m1 ,         
        const double           m2 ,       
        const double           m3 ,       
        const unsigned short   L  ) ;
      // ======================================================================
      /// copy constructor 
      BW3L ( const BW3L&  ) ;
      /// move constructor 
      BW3L (       BW3L&& ) = default ;
      // =======================================================================
    public: 
      // ======================================================================
      ///  fictive public default constructor (needed for serialization)
      BW3L() {};
      // =======================================================================
    public: 
      // ======================================================================
      /// evaluate the function 
      double evaluate   ( const double x ) const ;
      /// evaluat ethe function 
      double operator() ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      inline const Ostap::Math::BW& breit_wigner() const { return *m_bw.get() ; }
      inline       Ostap::Math::BW& breit_wigner()       { return *m_bw.get() ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get an integral 
      double integral () const ;
      /// get an integral 
      double integral ( const double xmin , const double xmax ) const ;
      // ======================================================================
    public:
      // ======================================================================
      double xmin () const { return std::max ( m_bw->threshold() , m_m1 + m_m2 ) ; }
      double xmax () const { return                                m_M  - m_m3   ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the mass of mother pareticle 
      double M  () const { return m_M ; } // get  mass of mother particle 
      /// get the mass of 1st daughter particle 
      double m1 () const { return m_m1 ; } // get mass of 1st daughter particle 
      /// get the mass of 2nd daughter particle 
      double m2 () const { return m_m2 ; } // get mass of 2nd daughter particle 
      /// get the mass of 3rd daughter particle 
      double m3 () const { return m_m3 ; } // get mass of 3rd daughter particle 
      /// get the orbital momentum between (1,2) and (3) 
      unsigned short L () const { return m_L ; } // get the orbital momentum between (1,2) and (3)
      // ======================================================================
    public:
      // ======================================================================
      /// get the  gamma for the certain channel 
      double gamma    ( const unsigned short i ) const 
      { return m_bw->gamma ( i ) ; }
      /// set the gamma for the certain decay
      bool   setGamma ( const unsigned short i , const double value ) 
      { return m_bw->setGamma ( i , value ) ; }
      /// get the total gamma 
      double    gamma () const { return m_bw->gamma() ; }
      /// set the total gamma 
      bool   setGamma ( const double value ) 
      { return m_bw->setGamma ( value ) ; }
      /// get number of channels 
      unsigned int nChannels() const { return m_bw->nChannels ()  ; }
      /// pole position 
      double       m0     ()   const { return m_bw->m0()     ; }
      /// set pole position 
      bool      setM0     ( const double x ) { return m_bw->setM0 ( x ) ; }
      // get the amplitude 
      std::complex<double> amplitude ( const double m ) const 
      { return m_bw->amplitude  ( m ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// some unique tag 
      std::size_t tag() const ;
      // ======================================================================
    private:
      // ======================================================================
      /// Breit-wigner 
      std::unique_ptr<Ostap::Math::BW>  m_bw  { nullptr } ;
      /// mass of mother particle 
      double          m_M  { 1   } ; // mass of mother particle 
      /// mass of 1st daughter particle 
      double          m_m1 { 0   } ; // mass of 1st daughter particle 
      /// mass of 2nd daughter particle 
      double          m_m2 { 0   } ; // mass of 2nd daughter particle 
      /// mass of 3rd daughter particle 
      double          m_m3 { 0   } ; // mass of 3rd daughter particle 
      /// orbital momentum between (1,2) and (3) 
      unsigned  short m_L  { 0   } ; // orbital momentum between (1,2) and (3) 
      // =======================================================================
      /** momentum of 3rd daughter at
       *   \f$ m^{\ast}_{12} = \frac{1}{2}(m_{12}^{\mathrm{min}} + m_{12}^{\mathrm{max}})\f$ 
       */
      double          m_p0 { 0.5 } ; // momentum of 3rd daughter at some mid-point
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace      m_workspace {} ; // integration workspace
      // ======================================================================
    }; 
    // ========================================================================
    /** @class A1 
     * Squared Breit-Wigner amplitude 
     * @see Ostap::Math::BW 
     */
    class A2 
    {
    public :
      // ======================================================================
      /// constructor from the breit-wigner
      A2 
      ( const BW&    bw           ,
        const double scale = 1.0  ) ;
      // ======================================================================
      /// copy constructor 
      A2
      ( const A2&  bw ) ;
      /// Move constructor 
      A2 (       A2&& bw ) = default ;
      // ======================================================================
    public:
      // ======================================================================
      double operator() ( const double s ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the breint wigner
      const BW& bw    () const { return *m_bw.get() ; }
      double    scale () const { return m_scale      ; }
      // ======================================================================
    public:
      // ======================================================================
      /// unique tag 
      std::size_t tag () const ;
      // ======================================================================
    private: 
      // ======================================================================
      /// Breit-Wigner itself 
      std::unique_ptr<BW>  m_bw    { nullptr } ; // Breit-Wigner 
      // scale factor 
      double               m_scale { 1.0     } ; // scale-factor 
      // ======================================================================
    } ;  
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_BREITWIGNER_H
// ============================================================================
