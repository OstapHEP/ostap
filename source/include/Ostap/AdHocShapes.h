// ============================================================================
#ifndef OSTAP_ADHOCSHAPES_H
#define OSTAP_ADHOCSHAPES_H 1
// ============================================================================
// include files
// ============================================================================
// STD & STL
// ============================================================================
#include <vector>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Bernstein1D.h"
#include "Ostap/Workspace.h"
#include "Ostap/PhaseSpace.h"
#include "Ostap/MoreMath.h"
// ============================================================================
/** @file Ostap/AdHocShapes.h
 *  Set of useful ad-hoc, sometem phsyicmotovated, shapes
 *
 *  - Exponential modulated by polynomial
 *  - Phase space modulated by polynomial 
 *  - (left) Phase space tiemx exponential modulated by polynomial
 *  - Sigmoid/kink function modulated by polynomial
 *  - Difference of two exponents 
 *  - Difference of two exponents modulated by positive polynomial
 *  - Argus & GenArgus
 *  - 
 *
 *  @see Ostap::Math::ExpoPol
 *  @see Ostap::Math::PhaseSpacePol
 *  @see Ostap::Math::PhaseSpaceLeftExpoPol
 *  @see Ostap::Math::Sigmoid
 *  @see Ostap::Math::TwoExpos
 *  @see Ostap::Math::TwoExpoPositive
 *  @see Ostap::Math::Argus 
 *  @see Ostap::Math::GenArgus 
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
    /** @class ExpoPositive
     *  Useful function for parameterizing smooth background:
     *  product of the exponential and positive polinonmial
     *  @see Ostap::Math::Positive
     */
    class  ExpoPositive  final : public Ostap::Math::PolyFactor1D
    {
    public:
      // ======================================================================
      /// constructor from the order
      ExpoPositive 
      ( const unsigned short       N     =  0 ,
        const double               tau   =  0 , // exponent
        const double               xmin  =  0 ,
        const double               xmax  =  1 ) ;
      // ======================================================================
      /// constructor from N phases
      ExpoPositive 
      ( const std::vector<double>& pars       ,
        const double               tau   =  0 , // exponent
        const double               xmin  =  0 ,
        const double               xmax  =  1 ) ;
      // ======================================================================
      /// constructor from polynom and exponential 
      ExpoPositive 
      ( const Ostap::Math::Positive& pol        , 
        const double                 tau   =  0 ) ;// exponent
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      inline double operator () ( const double x ) const { return  evaluate ( x ) ; } 
      /// get the value
      double         evaluate   ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get exponential
      inline double tau    () const { return m_tau ;}
      /// set new value for the exponent
      bool          setTau ( const  double value ) ;
      // ======================================================================
    public: // own parameters: tau
      // ======================================================================
      /// own parameters: tau
      inline std::size_t         npars_own () const { return 1 ; }
      inline std::vector<double> own_pars  () const { return { tau() } ;  }
      inline double own_par   ( const unsigned short k ) const
      {	return 0 == k ? tau () : 0.0 ; }
      inline bool   setOwnPar
      ( const unsigned short k     , 
	const double         value ) 
      { return 0 == k ? setTau ( value ) : false ; }
      /// all parameters 
      inline std::vector<double> all_pars () const
      { return Ostap::Math::Parameters::join ( tau () , pars () ) ;  }      
      // ======================================================================
    public:
      // ======================================================================      
      /// get the value \f$ x_{min}\$ such that  \f$ x_{min} \le p(x) \f$ 
      double min_value () const ;
      /// get the value \f$ x_{max}\$ such that  \f$ x_{max} \ge p(x) \f$ 
      double max_value () const ;
      // =====================================================================
    public:
      // ======================================================================
      double integral () const ;
      double integral
      ( const double low , 
	const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// exponential slope 
      double      m_tau { 0 } ; // exponential slope 
      // ======================================================================
    };
    // ========================================================================
    /** @class PhaseSpaceLeftExpoPol
     *  Function to represent the product of l-body phase space, 
     *  positive polynomial and the exponential function 
     *  \f[ \Phi_{l}^{(N)(x)} \propto
     *      \Phi_{l}(x;x_{low}) \mathrm{e}^{-\left|\tau\right| x } P_{N}(x) \f]
     *  where :
     *  -  \f$  \Phi_{l}(x;x_{low}) \f$  is a phase space of 
     *     l-particles near the threshold 
     *  -  \f$ P_{N}(x) \f$ is a positive polynomial of degree N
     *  
     *  @see Ostap::Math::PhaseSpaceLeft
     *  @see Ostap::Math::Positive
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2018-10-21
     */
    class PhaseSpaceLeftExpoPol final : public Ostap::Math::PolyFactor1D 
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from threshold and number of particles
       *  @param threshold_L the low-mass  threshold
       *  @param l           how many particles we consider
       *  @param N           degree of polynomial
       *  @param tau         the exponent 
       *  @param xhigh       the high edge 
       */
      PhaseSpaceLeftExpoPol 
      ( const double         threshold_L =  0 ,   // low threshold 
	const unsigned short l           =  2 ,   // number of particles 
	const unsigned short N           =  1 ,   // degree of polynomial
	const double         tau         =  0 ,   // the exponent 
	const double         xhigh       =  1 ) ; // high edge 
      // =====================================================================
      /** constructor from threshold and number of particles
       *  @param threshold_L the low-mass  threshold
       *  @param l           how many particles we consider
       *  @param N           degree of polynomial
       *  @param tau         the exponent 
       *  @param xlow        the low  edge 
       *  @param xhigh       the high edge 
       */
      PhaseSpaceLeftExpoPol
      ( const double         threshold_L ,   // low threshold 
	const unsigned short l           ,   // number of particles 
	const unsigned short N           ,   // degree of polynomial
	const double         tau         ,   // the exponent 
	const double         xlow        ,   // low edge 
	const double         xhigh       ) ; // high edge 
      // =====================================================================
      /** constructor from the phase space and polynomial degree
       *  @param ps          phase space factor
       *  @param N           degree of polynomial
       *  @param tau         the exponent 
       *  @param xhigh       the high edge 
       */
      PhaseSpaceLeftExpoPol
      ( const PhaseSpaceLeft& ps        ,
	const unsigned short  N     = 1 ,   // degree of polynomial
	const double          tau   = 0 ,   // the exponent 
	const double          xhigh = 1 ) ; // high edge 
      // =========================================================================
      /** constructor from the phase space and polynomial degree
       *  @param ps          phase space factor
       *  @param N           degree of polynomial
       *  @param tau         the exponent 
       *  @param xlow        the low  edge 
       *  @param xhigh       the high edge 
       */
      PhaseSpaceLeftExpoPol
      ( const PhaseSpaceLeft& ps    ,
	const unsigned short  N     ,   // degree of polynomial
	const double          tau   ,   // the exponent 
	const double          xlow  ,   // low edge 
	const double          xhigh ) ; // high edge
      // ======================================================================
      /** constructor from the phase space and polynomial degree
       *  @param ps          phase space factor
       *  @param N           degree of polynomial
       *  @param tau         the exponent 
       *  @param xlow        the low  edge 
       *  @param xhigh       the high edge 
       */
      PhaseSpaceLeftExpoPol
      ( const PhaseSpaceLeft&        ps  ,   // pjase space 
	const Ostap::Math::Positive& pol ,   // polynomial 
	const double                 tau ) ; // the exponent 
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate modulated phase space
      double evaluate    ( const double x ) const ;
      /// evaluate modulated phase space
      inline double operator () ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public: 
      // ======================================================================
      /// exponential slope 
      inline double tau       () const { return m_tau ; } 
      /// phase space scale 
      inline double scale     () const { return m_phasespace.scale     () ; }
      /// get the threshold  
      inline double threshold () const { return m_phasespace.threshold () ; }
      // ======================================================================
    public: // own parameters: tau & scale  
      // ======================================================================
      /// own parameters: tau, scale
      inline std::size_t         npars_own () const { return 2 ; }
      inline std::vector<double> own_pars  () const { return { tau() , scale () } ;  }
      inline double own_par   ( const unsigned short k ) const
      {	return 0 == k ? tau () : 1 == k ? scale () : 0.0 ; }
      inline bool   setOwnPar
      ( const unsigned short k     , 
	const double         value ) 
      { return
	  0 == k ? setTau   ( value ) :
	  1 == k ? setScale ( value ) : false ; }
      /// all parameters 
      inline std::vector<double> all_pars () const
      { return Ostap::Math::Parameters::join ( tau () , scale () , pars () ) ; }      
      // ======================================================================
    public:
      // ======================================================================
      /// get the phase space 
      inline const Ostap::Math::PhaseSpaceLeft& phasespace () const
      { return m_phasespace ; }
      // ======================================================================
    public:
      // ======================================================================
      /// set the new exponent 
      bool        setTau   ( const double value ) ;
      /// set the   scale  
      inline bool setScale ( const double value ) 
      { return m_phasespace.setScale ( value ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between low and high limits
      double integral
      ( const double low  ,
	const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::PhaseSpaceLeft* operator->() const 
      { return &m_phasespace ; }
      // ======================================================================
    private:
      // ======================================================================
      /// the phase space
      Ostap::Math::PhaseSpaceLeft m_phasespace ; // the phase space
      double                      m_tau        ; // the exponent
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace      m_workspace {} ; // integration workspace
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Sigmoid
     *  (shifted&scaled) Sigmoid/kink function, modulated by the positive polynomial
     *  \f$ f(x) = ( ( 1 - f ) + f_{\sigma} (z) + f  ) \times P_{pos} (x) \f$, where 
     *   - \f$ z = \frac{x-x_0}{\sigma} \f$ 
     *   - \f$ f_{\sigma}(x)   \ge 0 }\f$ is Sigmoid function 
     *   - \f$ P_{pos}(x)      \ge 0 }\f$ is Positive polynomial
     *   - shift \f$ \f = \sin^2 \updelta 0  \f$ 
     * 
     *  All sigmoid fuctions \f$ \sigma(z) \f$ are normalized & scaled such
     *  - \f$ \sigma(-\infty) =0\f$ 
     *  - \f$ \sigma(+\infty) =1\f$ 
     *  - \f$ \sigma^\prime(0)=1\f$ 
     *
     *  @see Ostap::Math::logistic 
     *  @see Ostap::Math::gd  
     *  @see Ostap::Math::smooth_transition 
     *  @see Ostap::Math::smooth_step 
     *
     *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
     *  @date 2015-02-07
     */
    class  Sigmoid  final :  public Ostap::Math::PolyFactor1D 
    {
    public:
      // ============================================================
      /// constructor from polynomial and parameters "alpha" and "x0"
      Sigmoid
      ( const Ostap::Math::Positive&   poly               ,
	const double                   scale = 1          ,
	const double                   x0    = 0          ,
	const double                   delta = 0          , 
	const Ostap::Math::SigmoidType st    = Ostap::Math::SigmoidType::Hyperbolic ) ;
      /// constructor from polynomial and parameter "alpha"
      Sigmoid
      ( const unsigned short           N     = 0          ,
	const double                   xmin  = 0          ,
	const double                   xmax  = 1          ,
	const double                   scale = 1          ,
	const double                   x0    = 0          , 
	const double                   delta = 0          , 
	const Ostap::Math::SigmoidType st    = Ostap::Math::SigmoidType::Hyperbolic ) ;
      /// constructor from polynomial and parameter "alpha"
      Sigmoid
      ( const std::vector<double>&     pars               ,
	const double                   xmin  = 0          ,
	const double                   xmax  = 1          ,
	const double                   scale = 1          ,
	const double                   x0    = 0          ,
	const double                   delta = 0          , 	
	const Ostap::Math::SigmoidType st    = Ostap::Math::SigmoidType::Hyperbolic ) ;
      // ========================================================================
      /// constructor from polynomial and parameters "alpha" and "x0"
      Sigmoid
      ( const std::string&             sigmoid_name       ,
	const Ostap::Math::Positive&   poly               ,
	const double                   scale = 1          ,
	const double                   x0    = 0          ,
	const double                   delta = 0          ) ;
      /// constructor from polynomial and parameter "alpha"
      Sigmoid
      ( const std::string&             sigmoid_name       , 
	const unsigned short           N     = 0          ,
	const double                   xmin  = 0          ,
	const double                   xmax  = 1          ,
	const double                   scale = 1          ,
	const double                   x0    = 0          , 
	const double                   delta = 0          ) ;
      /// constructor from polynomial and parameters "alpha" and "x0"
      Sigmoid
      ( const std::string&             sigmoid_name       , 
	const std::vector<double>&     pars               ,
	const double                   xmin  = 0          ,
	const double                   xmax  = 1          ,
	const double                   scale = 1          ,
	const double                   x0    = 0          ,
	const double                   delta = 0          ) ;      
      // ======================================================================
    public:
      // ======================================================================
      /// get the value of the function: 
      inline double operator () ( const double x ) const
      { 
        const double s2 = m_sin2delta ;
        const double c2 = 1 - s2      ;
        return xmin () <= x && x <= xmax() ?
               m_positive ( x ) * ( c2 * sigmoid ( x ) + s2 ) : 0.0 ; 
      }
      // ======================================================================
    public:
      // ======================================================================
      /** Get the actual sigmoid/kink value
       *  All sigmoids are normalized to have the same slope at the x=x0
       */
      inline double sigmoid ( const double x ) const
      {
	const double z = ( x - m_x0 ) / m_scale ;
	return Ostap::Math::sigmoid ( z , m_type ) ;
      }
      // ======================================================================
    public: // sigmoid getters 
      // ======================================================================
      inline double                   x0           () const { return     m_x0        ; }
      inline double                   scale        () const { return     m_scale     ; }
      inline double                   delta        () const { return     m_delta     ; } 
      inline Ostap::Math::SigmoidType sigmoid_type () const { return     m_type      ; }
      inline double                   sin2delta    () const { return     m_sin2delta ; }
      inline double                   cos2delta    () const { return 1 - m_sin2delta ; }
      /// the name of sigmoid function 
      std::string                     sigmoid_name () const ;
      // ======================================================================
    public: // sigmoid setters  
      // ======================================================================      
      // set new value of x0 
      bool setX0          ( const double value ) ;
      // set new value of scale 
      bool setScale       ( const double value ) ;
      // set new value of delta 
      bool setDelta       ( const double value ) ;
      // ======================================================================
    public: // own parameters: x0 , scale & delta 
      // ======================================================================
      /// own parameters: x0, scale , delta
      inline std::size_t         npars_own () const { return 3 ; }
      inline std::vector<double> own_pars  () const { return { m_x0 , m_scale , m_delta } ; }
      inline double own_par   ( const unsigned short k ) const
      {	return
	  0 == k ? x0    () :
	  1 == k ? scale () :
	  2 == k ? delta () : 0.0 ; }
      inline bool   setOwnPar
      ( const unsigned short k     , 
	const double         value ) 
      { return
	  0 == k ? setX0    ( value ) :
	  1 == k ? setScale ( value ) :
	  2 == k ? setDelta ( value ) : false ; }
      /// all parameters 
      inline std::vector<double> all_pars () const
      { return Ostap::Math::Parameters::join ( x0 () , scale () , delta () , pars () ) ;  }      
      // ======================================================================
    public:       
      // ======================================================================
      /// get the value \f$ x_{min}\$ such that  \f$ x_{min} \le p(x) \f$ 
      double min_value () const ;
      /// get the value \f$ x_{max}\$ such that  \f$ x_{max} \ge p(x) \f$ 
      double max_value () const ;
      // =====================================================================
    public:
      // ======================================================================
      /// get the integral between xmin and xmax
      double integral   () const ;
      /// get the integral between low and high
      double integral  
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// sigmoid location 
      double                   m_x0        { 0          } ; // sigmoid location 
      /// sigmoid scale 
      double                   m_scale     { 1          } ; // sigmoid scale 
      /// sigmoid delta  
      double                   m_delta     { 0          } ; // sigmoid delta 
      /// sigmoid type 
      Ostap::Math::SigmoidType m_type      { Ostap::Math::SigmoidType::Hyperbolic } ;
      /// constant fraction f = sin^2 delta 
      double                   m_sin2delta { 0          } ; // sin^2 delta 
      // ======================================================================
    private:
      // ======================================================================
      /// workspace for integration
      Ostap::Math::WorkSpace m_workspace ;
      // ======================================================================
    };
        // ========================================================================
    /** @class TwoExpos
     *  simple difference of two exponents
     *  \f$ f \propto
     *        \mathrm{e}^{-a_1    x}       -\mathrm{e}^{-a_2 x} =
     *        \mathrm{e}^{-\alpha x}\left(1-\mathrm{e}^{-\delta x}\right) \f$
     *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
     *  @date 2015-02-07
     */
    class  TwoExpos 
    {
    public:
      // ======================================================================
      TwoExpos 
      ( const double alpha = 1 ,
        const double delta = 1 ,
        const double x0    = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      inline double operator() ( const double x ) const { return evaluate ( x ) ; } 
      double        evaluate   ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get alpha
      inline double alpha () const { return m_alpha ; }
      /// get delta
      inline double delta () const { return m_delta ; }
      /// get x0
      inline double x0    () const { return m_x0    ; }
      // ======================================================================
    public:
      // ======================================================================
      /// slope for the first  exponent
      inline double a1         () const { return m_alpha           ; }
      /// slope for the second exponent
      inline double a2         () const { return m_alpha + m_delta ; }
      /// mean-value (for -inf,+inf) interval
      double mean       () const ;
      /// mode
      double mode       () const ;
      /// variance
      double variance   () const ;
      /// dispersion
      double dispersion () const { return variance () ; }
      /// sigma
      double sigma      () const ;
      // get normalization constant
      double norm       () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// slope for the first  exponent
      inline double tau1       () const { return -a1() ; }
      /// slope for the second exponent
      inline double tau2       () const { return -a2() ; }
      // ======================================================================
    public:
      // ======================================================================
      bool setAlpha ( const double value ) ;
      bool setDelta ( const double value ) ;
      bool setX0    ( const double value ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between -inf and +inf
      double integral    () const ;
      /// get the integral between low and high
      double integral 
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================      
      /// get the derivative at given value
      double derivative  ( const double x    ) const ;
      /// get the second at given value
      double derivative2 ( const double x    ) const ;
      /// get the Nth derivative at given value
      double derivative 
      ( const double   x  ,
        const unsigned N  ) const ;
      // ======================================================================
    public: 
      // ======================================================================
      /// get the value \f$ x_{max}\$ such that  \f$ x_{max} \ge p(x) \f$ 
      double max_value () const ;
      // =====================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      double m_alpha ;
      double m_delta ;
      double m_x0    ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class TwoExpoPositive
     *  simple difference of two exponents modulated with positive polynomials
     *  @see TwoExpos
     *  @see Positive
     *  @see ExpoPositive
     *  \f$ f(x) = e_2(x) * p_n(x) \f$, where
     *  \f$ e_2(x) \propto
     *        \mathrm{e}^{-a_1    x}       -\mathrm{e}^{-a_2 x} =
     *        \mathrm{e}^{-\alpha x}\left(1-\mathrm{e}^{-\delta x}\right) \f$
     *  and $p_2(s)$ is positive polynomial function
     *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
     *  @date 2015-03-28
     */
    class  TwoExpoPositive final : public Ostap::Math::PolyFactor1D
    {
    public:
      // ======================================================================
      TwoExpoPositive
        ( const unsigned short N = 1 ,
          const double alpha     = 1 ,
          const double delta     = 1 ,
          const double x0        = 0 ,
          const double xmin      = 0 ,
          const double xmax      = 1 ) ;
      // ======================================================================
      TwoExpoPositive
        ( const std::vector<double>& pars ,
          const double alpha  = 1 ,
          const double delta  = 1 ,
          const double x0     = 0 ,
          const double xmin   = 0 ,
          const double xmax   = 1 ) ;
      // ======================================================================
      TwoExpoPositive
        ( const Positive& poly    ,
          const double alpha  = 1 ,
          const double delta  = 1 ,
          const double x0     = 0 ) ;
      // ======================================================================
      TwoExpoPositive
        ( const Positive& poly   ,
          const TwoExpos& expos  ) ;
      // ======================================================================
      TwoExpoPositive
        ( const TwoExpos& expos  ,
          const Positive& poly   ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double operator() ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get alpha
      inline double alpha () const { return m_2exp.alpha () ; }
      /// get delta
      inline double delta () const { return m_2exp.delta () ; }
      /// get x0
      inline double x0    () const { return m_2exp.x0    () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// slope for the first  exponent
      inline double a1         () const { return m_2exp.a1   () ; }
      /// slope for the second exponent
      inline double a2         () const { return m_2exp.a2   () ; }
      /// slope for the first  exponent
      inline double tau1       () const { return m_2exp.tau1 () ; }
      /// slope for the second exponent
      inline double tau2       () const { return m_2exp.tau2 () ; }
      // ======================================================================
    public:
      // ======================================================================
      inline bool setAlpha ( const double value ) { return m_2exp.setAlpha ( value ) ; }
      inline bool setDelta ( const double value ) { return m_2exp.setDelta ( value ) ; }
      inline bool setX0    ( const double value ) { return m_2exp.setX0    ( value ) ; }
      // ======================================================================
    public: // own parameters: alpga, delta , x0
      // ======================================================================
      /// own parameters: tau, scale
      inline std::size_t         npars_own () const { return 3 ; }
      inline std::vector<double> own_pars  () const { return { alpha () , delta () , x0 () } ;  }
      inline double own_par   ( const unsigned short k ) const
      {	return
	  0 == k ? alpha () :
	  1 == k ? delta () :
	  2 == k ? x0    () : 0.0 ; }
      inline bool   setOwnPar
      ( const unsigned short k     , 
	const double         value )
      { return
	  0 == k ? setAlpha ( value ) :
	  1 == k ? setDelta ( value ) :
	  2 == k ? setX0    ( value ) : false ; }
      /// all parameters       
      inline std::vector<double> all_pars () const
      { return Ostap::Math::Parameters::join ( alpha () , delta () , x0 () , pars () ) ; }      
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between xmin and xmax
      double integral    () const ;
      /// get the integral between low and high
      double integral
      ( const double low  ,
	const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value \f$ x_{min}\$ such that  \f$ x_{min} \le p(x) \f$ 
      double min_value () const ;
      /// get the value \f$ x_{max}\$ such that  \f$ x_{max} \ge p(x) \f$ 
      double max_value () const ;
      // =====================================================================
    public:
      // ======================================================================
      /// get the underlying exponents
      const Ostap::Math::TwoExpos&  twoexpos  () const { return m_2exp      ; }
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag() const ;
      // ======================================================================
    private:
      // ======================================================================
      Ostap::Math::TwoExpos m_2exp     ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class HORNSdini 
     *  \f[ f(x;a,\delta, \phi) = 
     *  \frac{3}{2\delta}\left( z \right)^2
     *  \left( \cos^2( \phi + \frac{\pi}{4}) ( 1 + z ) +
     *         \sin^2( \phi + \frac{\pi}{4}) ( 1 - z ) \right) \f]
     *  where  \f$ z = \frac{ x - ( a - \delta ) } { \delta } \f$ 
     *  for \f$ a \le x \le a + 2\delta\$ and zero otherwise 
     *  
     * The first factor accouns for two-horn parabolic shape, 
     * and the second factor accounts for the linear correction factor 
     * ("efficiency")
     *
     *  - For the actual use it needs to be convoluted with resolution function
     * 
     *  @see https://arxiv.org/abs/2010.08483
     *  @see https://doi.org/10.48550/arXiv.2010.08483
     *  @see Aaij, R. et al., "Measurement of the CKM angle $\gamma$
     *       in $B^\pm\to D K^\pm$ and $B^\pm \to D \pi^\pm$ decays with
     *       $D \to K_\mathrm S^0 h^+ h^-$", JHEP  02 (2021) 169.
     *
     */
    class HORNSdini 
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor fron all parameters
       *  @param a position of th eleft paraboilic horn
       *  @param delta distance fron left to right parabolic horn 
       *  @param phi  correction parameter ("efficiency")
       */
      HORNSdini
      ( const double a     = 0 , 
        const double delta = 1 ,
        const double phi   = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evalaute the function 
      double evaluate ( const double x )  const ;
      /// evaluate the function 
      inline double operator () ( const double x )  const 
      { return evaluate ( x ) ;  }
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// left horn  
      inline double a     () const { return m_a     ; }
      /// right horn
      inline double b     () const { return m_a + 2 * m_delta ; }
      /// delta 
      inline double delta () const { return m_delta ; }
      /// phi 
      inline double phi   () const { return m_phi  ; }
      // ======================================================================
    public :
      // ======================================================================
      inline double xmin () const { return a () ; }
      inline double xmax () const { return b () ; }
      // ======================================================================
    public: // setters
      // ======================================================================
      bool setA      ( const double value ) ;
      bool setDelta  ( const double value ) ;
      bool setPhi    ( const double value ) ;
      // ======================================================================      
    public:
      // ====================================================================== 
      /// get the integral 
      double integral  () const ;
      /// get the integral between low and high
      double integral 
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      double m_a        { 0 } ;
      double m_delta    { 1 } ;
      double m_phi      { 0 } ;
      double m_cos2_phi { 0 } ;
      double m_sin2_phi { 0 } ;
      // ======================================================================      
    } ;
    // ========================================================================
    /** @class HILLdini 
     *  \f[ f(x;a,\delta, \phi) = 
     *  \frac{3}{2\delta}\left( 1 - z \right)^2
     *  \left( \cos^2( \phi + \frac{\pi}{4}) ( 1 + z ) +
     *         \sin^2( \phi + \frac{\pi}{4}) ( 1 - z ) \right) \f]
     *  where  \f$ z = \frac{ x - ( a - \delta ) } { \delta } \f$ 
     *  for \f$ a \le x \le a + 2\delta\$ and zero otherwise 
     *  
     * The first factor account for a parabolic shape, 
     * and the second factor accounts for the linear correction factor 
     * ("efficiency")
     *
     *  - For the actual use it needs to be convoluted with resolution function
     *
     *  @see https://arxiv.org/abs/2010.08483
     *  @see https://doi.org/10.48550/arXiv.2010.08483
     *  @see Aaij, R. et al., "Measurement of the CKM angle $\gamma$
     *       in $B^\pm\to D K^\pm$ and $B^\pm \to D \pi^\pm$ decays with
     *       $D \to K_\mathrm S^0 h^+ h^-$", JHEP  02 (2021) 169.       
     */
    class HILLdini 
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor fron all parameters
       *  @param a position of th eleft paraboilic horn
       *  @param delta distance fron left to right parabolic horn 
       *  @param phi  correction parameter ("efficiency")
       */
      HILLdini
      ( const double a     = 0 , 
        const double delta = 1 ,
        const double phi   = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evalaute the function 
      double evaluate ( const double x )  const ;
      /// evalaute the function 
      inline double operator () ( const double x )  const 
      { return evaluate ( x ) ;  }
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// left horn  
      inline double a     () const { return m_a     ; }
      /// right horn
      inline double b     () const { return m_a + 2 * m_delta ; }
      /// delta 
      inline double delta () const { return m_delta ; }
      /// phi 
      inline double phi   () const { return m_phi  ; }
      // ======================================================================
    public :
      // ======================================================================
      inline double xmin () const { return a () ; }
      inline double xmax () const { return b () ; }
      // ======================================================================
    public: // setters
      // ======================================================================
      bool setA      ( const double value ) ;
      bool setDelta  ( const double value ) ;
      bool setPhi    ( const double value ) ;
      // ======================================================================      
    public:
      // ====================================================================== 
      /// get the integral 
      double integral  () const ;
      /// get the integral between low and high
      double integral 
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      double m_a        { 0 } ;
      double m_delta    { 1 } ;
      double m_phi      { 0 } ;
      double m_cos2_phi { 0 } ;
      double m_sin2_phi { 0 } ;
      // ======================================================================      
    } ;    
    // ========================================================================
    /** @class CutOffGauss 
     *  Useful function for smooth Gaussian cut-off:
     *  \f[ f(x;x_0;\sigma) = \left\{ 
     *    \begin{array}{ll}
     *    1  & \mathrm{for~} x \le x_0  \\
     *    \mathrm{e}^{-\frac{1}{2}\left( \frac{ (x-x_0)^2}{\sigma^2} \right)}
     *       & \mathrm{for~} x >   x_0 
     *    \end{array}\right. \f] 
     */
    class CutOffGauss 
    {
    public:
      // ======================================================================
      /** Constructor from all parameters
       *  @param right dump direction
       *  @param x0    threshold value 
       *  @param sigma sigma  
       */
      CutOffGauss 
      ( const bool   right = true , 
        const double x0    = 0    , 
        const double sigma = 1    ) ;
      // ======================================================================
    public :
      // ======================================================================
      /// the main method 
      double operator() ( const double x ) const ;
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// right ?
      inline bool right   () const { return m_right ; }
      /// x_0 
      inline double x0    () const { return m_x0    ; }
      /// sigma
      inline double sigma () const { return m_sigma ; }      
      // ======================================================================
    public: // setters 
      // ======================================================================
      /// update x_0
      bool setX0    ( const double value ) ;
      /// update sigma 
      bool setSigma ( const double value ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between low and high
      double integral 
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      bool   m_right ;
      double m_x0    ;
      double m_sigma ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class CutOffStudent
     *  Useful function for smooth Student's t=-like (power-law) cut-off:
     *  \f[ f(x;x_0;\sigma) = \left\{ 
     *    \begin{array}{ll}
     *    1  & \mathrm{for~} x \le x_0  \\
     *    \left( \frac{1}{\nu} \left( \frac{(x-x_0)}{\sigma^2} \right)^{ - \frac{\nu+1}{2}} \right) 
     *       & \mathrm{for~} x >   x_0 
     *    \end{array}\right. \f] 
     */
    class CutOffStudent
    {
    public:
      // ======================================================================
      /** Constructor from all parameters
       *  @param right dump direction
       *  @param x0    threshold value 
       *  @param nu    parameter nu 
       *  @param sigma parameter sigma  
       */
      CutOffStudent 
      ( const bool   right = true , 
        const double x0    = 0    , 
        const double n     = 1    ,
        const double sigma = 1    ) ;
      // ======================================================================
    public :
      // ======================================================================
      /// the main method 
      double operator() ( const double x ) const ;
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// right ?
      inline bool right   () const { return m_right ; }
      /// x_0 
      inline double x0    () const { return m_x0    ; }
      /// n 
      inline double nu    () const { return m_nu    ; }      
      /// sigma
      inline double sigma () const { return m_sigma ; }      
      // ======================================================================
    public: // setters 
      // ======================================================================
      /// update x_0
      bool setX0    ( const double value ) ;
      /// update nu 
      bool setNu    ( const double value ) ;
      /// update sigma 
      bool setSigma ( const double value ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between low and high
      double integral
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      bool   m_right ;
      double m_x0    ;
      double m_nu    ;
      double m_sigma ;
      double m_C     ;
      // ======================================================================
    } ;
    // ========================================================================
    // ========================================================================
    /** @class Argus 
     *  Slightly modified version of Argus distribution, with 
     *  support in the interval  \f$ \mu - c \le x \le \mu \f$
     *  @see https://en.wikipedia.org/wiki/ARGUS_distribution
     *  @see ARGUS Collaboration, H. Albrecht et al., 
     *      "Measurement of the polarization in the decay B → J/ψK*". 
     *      Physics Letters B. 340 (3): 217–220.
     *  @see doi:10.1016/0370-2693(94)01302-0.
     *  @see https://doi.org/10.1016%2F0370-2693%2894%2901302-0
     */
    class Argus 
    {
      // ======================================================================
    public :
      // ======================================================================
      // constructor for all elements 
      Argus 
      ( const double mu  = 1 , 
        const double  c  = 1 , 
        const double chi = 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate function 
      double        evaluate   ( const double x ) const ;
      /// get PDF
      inline double pdf        ( const double x ) const { return evaluate ( x ) ; }
      /// get PDF
      inline double operator() ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public:  // gettters 
      // ====================================================================== 
      /// parameter mu 
      inline double mu  () const { return m_mu  ; }
      /// parameter c 
      inline double c   () const { return m_c   ; }
      /// parameter chi 
      inline double chi () const { return m_chi ; }        
      // ======================================================================
    public: // setters 
      // ======================================================================
      /// set mu parameter
      bool setMu  ( const double value ) ;
      /// set c parameter
      bool setC   ( const double value ) ;
      /// set chi parameter
      bool setChi ( const double value ) ;
      // ======================================================================
    public: // properties 
      // ======================================================================
      /// mean of the distribution 
      double mean     () const ;
      /// mode of the distribution 
      double mode     () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// xmin
      inline double xmin () const { return m_mu - m_c ; }
      /// xmax
      inline double xmax () const { return m_mu       ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral 
      double integral   () const ;
      /// get the integral between low and high
      double integral  
      ( const double low  ,
        const double high ) const ;
      /// get CDF 
      double cdf        ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private: // helper function
      // ======================================================================
      /** helper function 
       *  \f$ \Psi ( \chi ) = \Phi(\chi )  - \chi \phi  (\chi ) - \frac{1}{2} \f$ 
       */
      double psi ( const double value ) const ;
      // ======================================================================
   private:
      // ======================================================================
      /// parameter mu 
      double   m_mu   {  1 } ; // parameter mu      
      /// parameter c 
      double   m_c    {  1 } ; // parameter c 
      /// parameter chi 
      double   m_chi  {  1 } ; // parameter chi
      /// normalization 
      double   m_norm { -1 } ; // normalization
      // ======================================================================
    } ;
    // ========================================================================
    /** @class GenArgus 
     *  Slightly modified version of generalized Argus distribution, with 
     *  support in the interval  \f$ \mu - c \le x \le \mu \f$
     *  @see https://en.wikipedia.org/wiki/ARGUS_distribution
     *  @see ARGUS Collab oration, H. Albrecht et al., 
     *      "Measurement of the polarization in the decay B → J/ψK*". 
     *      Physics Letters B. 340 (3): 217–220.
     *  @see doi:10.1016/0370-2693(94)01302-0.
     *  @see https://doi.org/10.1016%2F0370-2693%2894%2901302-0
     */
    class GenArgus 
    {
      // ======================================================================
    public :
      // ======================================================================
      // constructor for all elements 
      GenArgus 
      ( const double mu  = 1   , 
        const double  c  = 1   ,  
        const double chi = 1   ,
        const double dp  = 1.5 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate function 
      double        evaluate   ( const double x ) const ;
      /// get PDF
      inline double pdf        ( const double x ) const { return evaluate ( x ) ; }
      /// get PDF
      inline double operator() ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public:  // gettters 
      // ====================================================================== 
      /// parameter mu 
      inline double mu  () const { return m_mu  ; }
      /// parameter c 
      inline double c   () const { return m_c   ; }
      /// parameter chi 
      inline double chi () const { return m_chi ; }        
      /// parameter dp 
      inline double dp  () const { return m_dp  ; }        
      // ======================================================================
    public:
      // ======================================================================
      /// parameter p 
      inline double p () const { return m_dp - 1 ; }
      // ======================================================================
    public: // setters 
      // ======================================================================
      /// set mu parameter
      bool setMu  ( const double value ) ;
      /// set c parameter
      bool setC   ( const double value ) ;
      /// set chi parameter
      bool setChi ( const double value ) ;
      /// set dp parameter
      bool setDp  ( const double value ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// xmin
      inline double xmin () const { return m_mu - m_c ; }
      /// xmax
      inline double xmax () const { return m_mu       ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral 
      double integral   () const ;
      /// get the integral between low and high
      double integral
      ( const double low  ,
	      const double high ) const ;
      /// get CDF 
      double cdf        ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
   private:
      // ======================================================================
      /// parameter mu 
      double   m_mu   {  1 } ; // parameter mu      
      /// parameter c 
      double   m_c    {  1 } ; // parameter c 
      /// parameter chi 
      double   m_chi  {  1 } ; // parameter chi
      /// parameter dp 
      double   m_dp   {  1 } ; // parameter chi
      /// normalization 
      double   m_norm { -1 } ; // normalization
      // ======================================================================
    } ;


    
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namesapce Ostap
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_ADHOCSHAPES_H
// ============================================================================
