// ============================================================================
// Include files
// ============================================================================
// STD & STL  
// ============================================================================
#include <cmath>
#include <array>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/PhaseSpace.h"
#include "Ostap/BreitWigner.h"
#include "Ostap/Clenshaw.h"
#include "Ostap/MoreMath.h"
#include "Ostap/Kinematics.h"
#include "Ostap/DalitzIntegrator.h"
#include "Ostap/Workspace.h"
#include "Ostap/Models.h"
#include "Ostap/StatusCode.h"
// ============================================================================
// local
// ============================================================================
#include "local_math.h"
#include "local_gsl.h"
#include "Integrator1D.h"
#include "status_codes.h"
// ============================================================================
/** @file 
 *  implementation of useful models for describing signal peaks with the natural width \
 *  - Breit-Wigner
 *  - Flatte 
 *  - LASS  (kappa) 
 *  - Bugg  (sigma-pole)
 *  - Gounaris-Sakurai
 *
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2010-04-19
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  /// imaginary unit 
  const std::complex<double> s_j { 0.0 , 1.0 } ;
  // ==========================================================================
  /// @var s_iPI  
  const double s_iPI = 1.0 / M_PI ;
  // ==========================================================================
  /** @var s_WS 
   *  local integration workspace 
   */
  const Ostap::Math::WorkSpace s_WS ;
  // ==========================================================================
} //                                            The end of  anonymous namespace 
// ============================================================================
// Rho-functions from Jackson
// ============================================================================
/* the simplest function: constant
 *  @see Ostap::Math::BreitWigner
 *  @see Ostap::Math::BreitWigner::rho_fun
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2011-11-30
 */
// ============================================================================
double Ostap::Math::Jackson::jackson_0
( double /* m  */ ,
  double /* m0 */ ,
  double /* m1 */ ,
  double /* m2 */ ) { return 1 ; }
// ============================================================================
/* the simple function for \f$ 1^- \rightarrow 0^- 0^- \f$, l = 1
 *  \f$\rho(\omega)= \omega^{-1}\f$
 *  @see Ostap::Math::BreitWigner
 *  @see Ostap::Math::BreitWigner::rho_fun
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2011-11-30
 */
// ============================================================================
double Ostap::Math::Jackson::jackson_A2
( double    m     ,
  double /* m0 */ ,
  double /* m1 */ ,
  double /* m2 */ ) { return 1./m ; }
// ============================================================================
/*  the simple function for \f$ 1^- \rightarrow 0^- 1^- \f$, l = 1
 *  \f$\rho(\omega)= \omega \f$
 *  @see Ostap::Math::BreitWigner
 *  @see Ostap::Math::BreitWigner::rho_fun
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2011-11-30
 */
// ============================================================================
double Ostap::Math::Jackson::jackson_A3
( double    m     ,
  double /* m0 */ ,
  double /* m1 */ ,
  double /* m2 */ ) { return m ; }
// ============================================================================
/*  the simple function for
 *  \f[ \frac{3}{2}^+ \rightarrow \frac{1}{2}^+ 0^- \f], l = 1
 *  $\rho(\omega)= \frac{ ( \omega + M )^2 - m^2 }{ \omega^2} \f$
 *  @see Ostap::Math::BreitWigner
 *  @see Ostap::Math::BreitWigner::rho_fun
 *  @param m the invariant mass
 *  @param m1 the invariant mass of the first  (spinor) particle
 *  @param m2 the invariant mass of the second (scalar) particle
 *  @return the value of rho-function
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2011-11-30
 */
// ============================================================================
double Ostap::Math::Jackson::jackson_A4
( double    m     ,
  double /* m0 */ ,
  double    m1    ,
  double    m2    )
{
  const double a = m + m1 ;
  //
  return ( a * a  - m2 * m2 ) / ( m * m ) ;
}
// ============================================================================
/*  the simple function for
 *  \f$ \frac{3}{2}^- \rightarrow \frac{1}{2}^+ 0^- \f$, l = 2
 *  $\rho(\omega)= \left[ ( \omega + M )^2 - m^2 \right]^{-1} \f$
 *  @see Ostap::Math::BreitWigner
 *  @see Ostap::Math::BreitWigner::rho_fun
 *  @param m the invariant mass
 *  @param m1 the invariant mass of the first  (spinor) particle
 *  @param m2 the invariant mass of the second (scalar) particle
 *  @return the value of rho-function
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2011-11-30
 */
// ============================================================================
double Ostap::Math::Jackson::jackson_A5
( double    m     ,
  double /* m0 */ ,
  double    m1    ,
  double    m2    )
{
  const double a = m + m1 ;
  //
  return 1 / ( a * a  - m2 * m2 ) ;
}
// ============================================================================
/*  the simple function for \f$\rho^- \rightarrow \pi^+ \pi^-\f$
 *  \f$ 1- \rightarrow 0^- 0^- \f$, l = 1
 *  $\rho(\omega)= \left[ q_0^2 + q^2 \right]^{-1}f$
 *  @see Ostap::Math::BreitWigner
 *  @see Ostap::Math::BreitWigner::rho_fun
 *  @param m  the invariant mass
 *  @param m0 the nominal   mass
 *  @param m1 the invariant mass of the first  particle
 *  @param m2 the invariant mass of the second particle
 *  @return the value of rho-function
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2011-11-30
 */
// ============================================================================
double Ostap::Math::Jackson::jackson_A7
( double    m   ,
  double    m0  ,
  double    m1  ,
  double    m2  )
{
  //
  const double q  = Ostap::Math::PhaseSpace2::q ( m  , m1 , m2 ) ;
  const double q0 = Ostap::Math::PhaseSpace2::q ( m0 , m1 , m2 ) ;
  //
  if ( 0 >= q && 0 >= q0 ) { return 1 ; }
  //
  return 1. / ( q * q + q0 * q0 ) ;
}
// ============================================================================
// FORMFACTORS 
// ============================================================================
// virtual destructor 
// ============================================================================
Ostap::Math::FormFactor::~FormFactor (){}
// ============================================================================
// constructor from enum
// ============================================================================
Ostap::Math::FormFactors::Jackson::Jackson
( const Ostap::Math::FormFactors::JacksonRho rho ) 
  : Ostap::Math::FormFactor()
  , m_rho  ( rho         )
  , m_what ( "Jackson()" ) 
{
  switch ( m_rho )
  {
  case   Ostap::Math::FormFactors::Jackson_0  :
    m_what = "Jackson(Jackson_0)"  ;
    break ;
  case   Ostap::Math::FormFactors::Jackson_A2 :
    m_what = "Jackson(Jackson_A2)" ;
    break ;
  case   Ostap::Math::FormFactors::Jackson_A3 :
    m_what = "Jackson(Jackson_A3)" ;
    break ;
  case   Ostap::Math::FormFactors::Jackson_A4 :
    m_what = "Jackson(Jackson_A4)" ;
    break ;
  case   Ostap::Math::FormFactors::Jackson_A5 :
    m_what = "Jackson(Jackson_A5)" ;
    break ;
  case   Ostap::Math::FormFactors::Jackson_A7 :
    m_what = "Jackson(Jackson_A7)" ;
    break ;
  default : ;
  }
}
// ============================================================================
// unique tag/label
// ============================================================================
std::size_t Ostap::Math::FormFactors::Jackson::tag      () const
{ return std::hash<std::string>{} ( m_what ) ; }
// ============================================================================
// virtual destructor 
// ============================================================================
Ostap::Math::FormFactors::Jackson::~Jackson(){}
// ============================================================================
// clone method ("virtual constructor")
// ============================================================================
Ostap::Math::FormFactors::Jackson* 
Ostap::Math::FormFactors::Jackson:: clone() const 
{ return new Ostap::Math::FormFactors::Jackson ( *this ) ; }
// ============================================================================
// the only important method 
// ============================================================================
double Ostap::Math::FormFactors::Jackson::operator() 
  ( const double m  , const double m0 ,
    const double m1 , const double m2 ) const
{
  return 
    Ostap::Math::FormFactors::Jackson_0 == m_rho ? 
    ( m / m0 ) * Ostap::Math::Jackson::jackson_0   ( m , m0 , m1 , m2 ) : 
    Ostap::Math::FormFactors::Jackson_A2 == m_rho ? 
    ( m / m0 ) * Ostap::Math::Jackson::jackson_A2  ( m , m0 , m1 , m2 ) : 
    Ostap::Math::FormFactors::Jackson_A3 == m_rho ? 
    ( m / m0 ) * Ostap::Math::Jackson::jackson_A3  ( m , m0 , m1 , m2 ) : 
    Ostap::Math::FormFactors::Jackson_A4 == m_rho ? 
    ( m / m0 ) * Ostap::Math::Jackson::jackson_A4  ( m , m0 , m1 , m2 ) : 
    Ostap::Math::FormFactors::Jackson_A5 == m_rho ? 
    ( m / m0 ) * Ostap::Math::Jackson::jackson_A5  ( m , m0 , m1 , m2 ) : 
    Ostap::Math::FormFactors::Jackson_A7 == m_rho ? 
    ( m / m0 ) * Ostap::Math::Jackson::jackson_A7  ( m , m0 , m1 , m2 ) : 
    ( m / m0 ) * Ostap::Math::Jackson::jackson_0   ( m , m0 , m1 , m2 ) ;
}
// ============================================================================
// Blatt-Weisskopf formfactors 
// ============================================================================
namespace
{
  // ==========================================================================
  // Coefficients for Blatt-Weisskopf formfactors 
  // ==========================================================================
  // const std::array<int,1> s_BW_0 { {                               1 } } ;
  // const std::array<int,2> s_BW_1 { {                            1, 1 } } ;
  const std::array<int,3> s_BW_2 { {                        9,  3, 1 } } ;
  const std::array<int,4> s_BW_3 { {                 225,  45,  6, 1 } } ;
  const std::array<int,5> s_BW_4 { {         11025, 1575, 135, 10, 1 } } ;
  const std::array<int,6> s_BW_5 { { 893025, 99225, 6300, 315, 15, 1 } } ;  
  // ==========================================================================
}    
// ============================================================================
// constructor
// ============================================================================
Ostap::Math::FormFactors::BlattWeisskopf::BlattWeisskopf
( const Ostap::Math::FormFactors::BlattWeisskopf::Case L , 
  const double                                         b )
  : Ostap::Math::FormFactor() 
  , m_L ( L ) 
  , m_b ( b )
  , m_what ( "BlattWeisskopf(" + std::to_string( (int) L ) + "," + std::to_string ( b ) + ")" )
{
  Ostap::Assert ( Zero  == L ||
                  One   == L ||
                  Two   == L ||
                  Three == L ||
                  Four  == L ||
                  Five  == L                                 ,
                  "Illegal orbital momentum!"                ,  
                  "Ostap::Math::FormFactors::BlattWeisskopf" ,
                  INVALID_PARS , __FILE__ , __LINE__         ) ; 
}

// ============================================================================
// default constructor (needed for  serialization)
// ============================================================================
Ostap::Math::FormFactors::BlattWeisskopf::BlattWeisskopf()
  : Ostap::Math::FormFactor() 
  , m_L    ( Ostap::Math::FormFactors::BlattWeisskopf::Zero ) 
  , m_b    ( 0.0 )
  , m_what ( "BlattWeisskopf(0,0.0)" )
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Math::FormFactors::BlattWeisskopf::~BlattWeisskopf(){}
// ============================================================================
// clone method ("virtual constructor")
// ============================================================================
Ostap::Math::FormFactors::BlattWeisskopf*
Ostap::Math::FormFactors::BlattWeisskopf::clone() const 
{ return new Ostap::Math::FormFactors::BlattWeisskopf(*this) ; }
// ============================================================================
// get the barrier factor 
// ============================================================================
double Ostap::Math::FormFactors::BlattWeisskopf::b 
( const double z   , 
  const double z0  ) const 
{
  if ( Zero == m_L || s_equal ( z , z0 ) ) { return 1 ; }
  //
  const long double r2 =
    //
    One   == m_L ? ( 1 + z0 ) / ( 1  + z )  :
    //
    Two   == m_L ? 
    Ostap::Math::Clenshaw::monomial_sum ( s_BW_2.rbegin() , s_BW_2.rend  () , z0 ).first /
    Ostap::Math::Clenshaw::monomial_sum ( s_BW_2.rbegin() , s_BW_2.rend  () , z  ).first :
    //
    Three == m_L ? 
    Ostap::Math::Clenshaw::monomial_sum ( s_BW_3.rbegin() , s_BW_3.rend  () , z0 ).first /
    Ostap::Math::Clenshaw::monomial_sum ( s_BW_3.rbegin() , s_BW_3.rend  () , z  ).first :
    //
    Four == m_L ? 
    Ostap::Math::Clenshaw::monomial_sum ( s_BW_4.rbegin() , s_BW_4.rend  () , z0 ).first /
    Ostap::Math::Clenshaw::monomial_sum ( s_BW_4.rbegin() , s_BW_4.rend  () , z  ).first :
    //
    Five == m_L ? 
    Ostap::Math::Clenshaw::monomial_sum ( s_BW_5.rbegin() , s_BW_5.rend  () , z0 ).first /
    Ostap::Math::Clenshaw::monomial_sum ( s_BW_5.rbegin() , s_BW_5.rend  () , z  ).first : 
    //
    1.0L ;
  //
  return r2 ;
}
// ============================================================================
/*  the only important method the squared ratio 
 *  of formfactors \f$  \frac{F^2(m)}{F^2(m_0)} \$
 */
// ============================================================================
double Ostap::Math::FormFactors::BlattWeisskopf::operator() 
  ( const double m  , const double m0 ,
    const double m1 , const double m2 ) const 
{
  //
  if ( s_equal ( m , m0 ) ) { return    1   ; }
  if ( s_zero  ( m_b    ) ) { return m0 / m ; }
  /// get the momenta 
  const double q  = Ostap::Math::PhaseSpace2::q ( m  , m1 , m2 ) ;
  const double q0 = Ostap::Math::PhaseSpace2::q ( m0 , m1 , m2 ) ;
  ///
  const double z  = q  * m_b ;
  const double z0 = q0 * m_b ;
  //
  return b ( z * z , z0 * z0 ) ;
}
// ============================================================================
// unique tag/label
// ============================================================================
std::size_t Ostap::Math::FormFactors::BlattWeisskopf::tag      () const
{
  static const std::string s_name = "BlattWeisskopf";
  return Ostap::Utils::hash_combiner ( s_name , m_what , (int) m_L , m_b ) ; 
}
// ============================================================================
// clone operation 
// ============================================================================
Ostap::Math::FormFactors::GenericFF*
Ostap::Math::FormFactors::GenericFF::clone() const 
{ return new Ostap::Math::FormFactors::GenericFF ( *this ) ; }
// ============================================================================

// ============================================================================
// constructor 
// ============================================================================
Ostap::Math::FormFactors::NoFormFactor::NoFormFactor () 
  : FormFactor() 
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Math::FormFactors::NoFormFactor::~NoFormFactor () {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Math::FormFactors::NoFormFactor*
Ostap::Math::FormFactors::NoFormFactor::clone() const 
{ return new NoFormFactor() ; }
// ===========================================================================
// unique tag/label
// ============================================================================
std::size_t Ostap::Math::FormFactors::NoFormFactor::tag      () const
{ 
  static const std::string s_name = "NoFormFactor";
  return Ostap::Utils::hash_combiner ( s_name ,  describe () ) ; 
}
// ============================================================================
// describe the formfactor
// ============================================================================
std::string Ostap::Math::FormFactors::NoFormFactor::describe () const
{ return "NoFormFactor"; }
// ============================================================================

   

            





// ============================================================================
// Breit-Wigner Channel 
// ============================================================================
// constructor from the partial width 
// ============================================================================
Ostap::Math::ChannelBW::ChannelBW
( const double gamma0 ) 
  : m_gamma0 ( 0 ) { setGamma0  ( gamma0 ) ; }
// ============================================================================
// virtual destructor 
// ============================================================================
Ostap::Math::ChannelBW::~ChannelBW (){}
// ============================================================================
// set new width 
// ============================================================================
bool Ostap::Math::ChannelBW::setGamma0 ( const double value ) 
{
  const double v = std::abs ( value ) ;
  if ( s_equal ( v , m_gamma0 ) ) { return false ; } // RETURN
  m_gamma0 = v ;
  //
  if ( s_zero ( m_gamma0 ) ) { m_gamma0 = 0 ; }
  //
  return true ;
}
// ============================================================================
/*  get the single channel amplitude
 *  \f[ \mathcal{A} =\left( m_0^2 - s - i D (s,m_0) \right)^{-1} \f] 
 *  @param s   \f$ s  \f$ -parameter
 *  @param m0  \f$ m_0\f$-parameter  
 *  @return amplitude 
 */
// ============================================================================
std::complex<double> 
Ostap::Math::ChannelBW::amplitude 
( const double s , 
  const double m0 ) const 
{ 
  const std::complex<double> d = m0 * m0 - s - s_j * D ( s , m0 ) ;
  return 1.0 / d ;
}
// ============================================================================
// unique tag/label
// ============================================================================
std::size_t Ostap::Math::ChannelBW::tag      () const
{ 
  static const std::string s_name = "ChannelBW";
  return Ostap::Utils::hash_combiner ( s_name , m_gamma0 ) ; 
}
// ============================================================================


// ============================================================================
// Breit-Wigner channel with mass-depedent width 
// ============================================================================
// full constructor with all functions specified 
// ============================================================================
Ostap::Math::ChannelWidth::ChannelWidth 
( const double       gamma       ,
  Width              width       ,  
  const double       sthreshold  , 
  const std::size_t  tag         ,
  const std::string& description ) 
  : ChannelBW     ( gamma ) 
  , m_w           ( width ) 
  , m_sthreshold  ( std::abs ( sthreshold ) ) 
  , m_tag         ( tag   ) 
  , m_description ( description ) 
{}
// ============================================================================
// clone method 
// ============================================================================
Ostap::Math::ChannelWidth*
Ostap::Math::ChannelWidth::clone () const
{ return new Ostap::Math::ChannelWidth ( *this ) ; }
// ============================================================================
// unique tag for this lineshape 
// ============================================================================
std::size_t Ostap::Math::ChannelWidth::tag () const
{ 
  static const std::string s_name = "ChannelWidth" ;
  return Ostap::Utils::hash_combiner 
    ( s_name           ,
      ChannelBW::tag() , 
      m_sthreshold     , 
      m_description    , 
      m_tag            ) ; 
}
// ============================================================================


// ============================================================================
// Breit-Wigner channel with mass-depedent width 
// ============================================================================
// full constructor with all functions specified 
// ============================================================================
Ostap::Math::ChannelGamma::ChannelGamma
( const double                     gamma       ,
  Ostap::Math::ChannelGamma::Width width       ,  
  const double                     sthreshold  , 
  const std::size_t                tag         ,
  const std::string&               description ) 
  : ChannelBW     ( gamma ) 
  , m_gamma       ( width ) 
  , m_sthreshold  ( std::abs ( sthreshold ) ) 
  , m_tag         ( tag   ) 
  , m_description ( description ) 
{}
// ============================================================================
// clone method 
// ============================================================================
Ostap::Math::ChannelGamma*
Ostap::Math::ChannelGamma::clone () const
{ return new Ostap::Math::ChannelGamma ( *this ) ; }
// ============================================================================
// unique tag for this lineshape 
// ============================================================================
std::size_t Ostap::Math::ChannelGamma::tag () const
{ 
  static const std::string s_name = "ChannelGamma" ;
  return Ostap::Utils::hash_combiner 
    ( s_name            , 
      ChannelBW::tag () , 
      m_sthreshold      , 
      m_description     , 
      m_tag             ) ; 
}
// ============================================================================






// ============================================================================
// Breit-Wigner constant width channel 
// ============================================================================
/*  constructor from all parameters and *NO* formfactor
 *  @param gamma the width 
 *  @param m1    the mass of the 1st daughter
 *  @param m2    the mass of the 2nd daughter
 */
// ============================================================================
Ostap::Math::ChannelCW::ChannelCW 
( const double                  gamma , 
  const double                  m1    , 
  const double                  m2    ) 
  : ChannelBW ( gamma   ) 
  , m_ps2     ( m1 , m2 ) 
{}
// ============================================================================
// clone the channel 
// ============================================================================
Ostap::Math::ChannelCW*
Ostap::Math::ChannelCW::clone() const 
{ return new Ostap::Math::ChannelCW ( *this ) ; }
// ============================================================================
// get the phase space factor  \f$ \varrho(s) \f$
// ============================================================================
double Ostap::Math::ChannelCW::rho_s 
( const double s  , 
  const double mn ) const
{ 
  const double a = m_ps2.rho_s ( s ) ;
  return mn <= m_ps2.threshold () ? a : a / m_ps2.rho ( mn ) ;
} 
// ============================================================================
// unique tag for this lineshape 
// ============================================================================
std::size_t Ostap::Math::ChannelCW::tag () const
{
  static const std::string s_name = "ChannelCW" ;
  return Ostap::Utils::hash_combiner
    ( s_name            , 
      ChannelBW::tag () , 
      m_ps2. tag     () ) ; 
}
// ============================================================================
// describe the channel 
// ============================================================================
std::string Ostap::Math::ChannelCW::describe() const 
{
  return 
    "ChannelCW(" + std::to_string ( gamma0 () ) + 
    ","          + std::to_string ( m1     () ) + 
    ","          + std::to_string ( m2     () ) + ")" ;
}
// ============================================================================
// the first main method: numerator
// ============================================================================
double Ostap::Math::ChannelCW::N2 
( const double s  , 
  const double m0 ) const 
{
  return s <= m_ps2.s_threshold() ? 0.0 : m0 * gamma0 () ;
}
// ============================================================================
// the second main method: (constant) term to the denominator 
// ============================================================================
std::complex<double> 
Ostap::Math::ChannelCW::D
( const double s  , 
  const double m0 ) const { return m0 * gamma0 () ; }
// ============================================================================
// get the opening threshold for the channel 
// ============================================================================
double Ostap::Math::ChannelCW::s_threshold () const 
{ return m_ps2.s_threshold () ; }


// ============================================================================
// Ostap::Math::Channel
// ============================================================================

// ============================================================================
/*  constructor from all parameters and no formfactor
 *  @param gamma the width 
 *  @param m1    the mass of the 1st daughter
 *  @param m2    the mass of the 2nd daughter
 *  @param L     the oribital momentum
 */
// ============================================================================
Ostap::Math::Channel::Channel
( const double                  gamma , 
  const double                  m1    , 
  const double                  m2    , 
  const unsigned short          L     )
  : ChannelCW    ( gamma , m1 , m2 ) 
  , m_L          ( L        )
  , m_formfactor ( new Ostap::Math::FormFactors::NoFormFactor () )
{}
// ============================================================================
/*  constructor from all parameters and no formfactor
 *  @param gamma the width 
 *  @param m1    the mass of the 1st daughter
 *  @param m2    the mass of the 2nd daughter
 *  @param L     the oribital momentum
 *  @param r     the Jackson's formfactor 
 */
// ============================================================================
Ostap::Math::Channel::Channel
( const double                                gamma , 
  const double                                m1    , 
  const double                                m2    , 
  const unsigned short                        L     ,
  const Ostap::Math::FormFactors::JacksonRho  r     )
  : ChannelCW    ( gamma , m1 , m2 )
  , m_L          ( L )
  , m_formfactor ( new Ostap::Math::FormFactors::Jackson ( r ) )  
{}
// ============================================================================
/*  constructor from all parameters and generic formfactor 
 *  @param gamma the width 
 *  @param m1    the mass of the 1st daughter
 *  @param m2    the mass of the 2nd daughter
 *  @param L     the oribital momentum
 *  @param f     the formfactor 
 *  @param qs    the scale momentum 
 */
// ============================================================================
Ostap::Math::Channel::Channel
( const double                   gamma , 
  const double                   m1    , 
  const double                   m2    , 
  const unsigned short           L     ,
  const Ostap::Math::FormFactor& f     )
  : ChannelCW    ( gamma , m1 ,  m2  )
  , m_L          ( L         )
  , m_formfactor ( f.clone() )  
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Math::Channel::Channel 
( const Ostap::Math::Channel& right ) 
  : ChannelCW    ( right       ) 
  , m_L          ( right.m_L   ) 
  , m_formfactor ( right.m_formfactor ? right.m_formfactor->clone() : nullptr ) 
{}
// ============================================================================
// clone the channel 
// ============================================================================
Ostap::Math::Channel*
Ostap::Math::Channel::clone () const 
{ return new Ostap::Math::Channel(*this) ; }
// ============================================================================
// the first main method: numerator 
// ============================================================================
double Ostap::Math::Channel::N2 
( const double s  , 
  const double m0 ) const 
{
  if ( s <= ps2().s_threshold () ) { return 0 ; }   // RETURN 
  //
  // (1) trivial factor
  double result  = m0 * gamma0 () ;  
  //
  // (2) formfactor 
  const FormFactor* F   = formfactor () ;
  if ( nullptr != F ) 
  {
    const double x = std::sqrt ( s ) ;
    result *= (*F)( x , m0 , m1 () , m2 () ) ; 
  }
  //
  // (3) orbital barrier 
  if ( 0 != m_L ) 
  {
    const double q  = ps2 () . q_s ( s       ) ;
    const double q0 = ps2 () . q_s ( m0 * m0 ) ;
    if (  0 < q0 ) { result *= std::pow ( q / q0 , 2 * m_L ) ; }
  }
  //
  return result ;
}
// ============================================================================
// the second main method: term to the denominator 
// ============================================================================
std::complex<double> 
Ostap::Math::Channel::D 
( const double s  , 
  const double m0 ) const
{
  //
  if ( s <= ps2().s_threshold () ) { return 0 ; }  
  //
  double result = m0 * gamma0 () ;
  //
  // (1) phase space factors 
  const double rho  = ps2 () . rho_s ( s       ) ;
  const double rho0 = ps2 () . rho_s ( m0 * m0 ) ;
  //
  if ( 0 < rho0 ) { result *= rho / rho0 ; }                         // NB!!
  //
  // (2) formfactor 
  const FormFactor* F   = formfactor () ;
  if ( nullptr != F ) 
  {
    const double x = std::sqrt ( s ) ;
    result *= (*F)( x , m0 , m1 () , m2 () ) ; 
  } 
  //
  // (3) orbital barrier 
  if ( 0 != m_L )
  {
    const double q  = ps2() . q_s ( s       ) ;
    const double q0 = ps2() . q_s ( m0 * m0 ) ;
    if ( 0 < q0 ) { result *=  std::pow ( q / q0  , 2 * m_L ) ; }    // NB!
  }
  //
  return result ;
}
// ============================================================================
// unique tag for this lineshape 
// ============================================================================
std::size_t Ostap::Math::Channel::tag() const
{ 
  static const std::string s_name = "Channel" ;
  return Ostap::Utils::hash_combiner
    ( s_name                         , 
      Ostap::Math::ChannelCW::tag()  , 
      m_L                            , 
      m_formfactor  ? m_formfactor->tag() : 0 ) ;
}
// ============================================================================
// describe the channel 
// ============================================================================
std::string Ostap::Math::Channel::describe() const 
{
  return 
    "Channel(" + std::to_string ( gamma0 () ) + 
    ","         + std::to_string ( m1     () ) + 
    ","         + std::to_string ( m2     () ) + 
    ( m_formfactor ? ( "," + m_formfactor->describe() ) : "" )  
    + ")";
}
// ============================================================================


// ============================================================================
// Ostap:::Math::Channel0 
// ============================================================================

// ============================================================================
/*  constructor from all parameters and no formfactor
 *  @param gamma the width 
 *  @param m1    the mass of the 1st daughter
 *  @param m2    the mass of the 2nd daughter
 *  @param L     the oribital momentum
 *  @param qs    the scale momentum 
 */
// ============================================================================
Ostap::Math::Channel0::Channel0
( const double                  gamma , 
  const double                  m1    , 
  const double                  m2    , 
  const unsigned short          L     , 
  const double                  qs    )
  : Channel ( gamma , m1 , m2 , L ) 
  , m_qs    ( s_zero ( qs ) ? 0.0 : std::abs ( qs ) )
{}
// ============================================================================
/*  constructor from all parameters and no formfactor
 *  @param gamma the width 
 *  @param m1    the mass of the 1st daughter
 *  @param m2    the mass of the 2nd daughter
 *  @param L     the oribital momentum
 *  @param r     the Jackson's formfactor 
 *  @param qs    the scale momentum 
 */
// ============================================================================
Ostap::Math::Channel0::Channel0
( const double                                gamma , 
  const double                                m1    , 
  const double                                m2    , 
  const unsigned short                        L     ,
  const Ostap::Math::FormFactors::JacksonRho  r     ,
  const double                                qs    ) 
  : Channel ( gamma , m1 , m2 , L , r )
  , m_qs    ( s_zero ( qs ) ? 0.0 : std::abs ( qs ) )
{}
// ============================================================================
/*  constructor from all parameters and generic formfactor 
 *  @param gamma the width 
 *  @param m1    the mass of the 1st daughter
 *  @param m2    the mass of the 2nd daughter
 *  @param L     the oribital momentum
 *  @param f     the formfactor 
 *  @param qs    the scale momentum 
 */
// ============================================================================
Ostap::Math::Channel0::Channel0
( const double                   gamma , 
  const double                   m1    , 
  const double                   m2    , 
  const unsigned short           L     ,
  const Ostap::Math::FormFactor& f     ,
  const double                   qs    ) 
  : Channel ( gamma , m1 , m2 , L , f ) 
  , m_qs    ( s_zero ( qs ) ? 0.0 : std::abs ( qs ) )
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Math::Channel0::Channel0 
( const Ostap::Math::Channel0&  right ) 
  : Channel ( right      )
  , m_qs    ( right.m_qs )
{}
// ============================================================================
// clone the channel 
// ============================================================================
Ostap::Math::Channel0*
Ostap::Math::Channel0::clone () const 
{ return new Ostap::Math::Channel0(*this) ; }
// ============================================================================
// the first main method: numerator 
// ============================================================================
double Ostap::Math::Channel0::N2 
( const double s  , 
  const double m0 ) const 
{
  //
  // (1) coupling
  double result  = g2 () ;
  //
  // (2) formfactor 
  const FormFactor* F   = formfactor () ;
  if ( nullptr != F ) 
  { 
    const double x = std::sqrt ( s ) ;
    result *= (*F)( x , m0 , m1 () , m2 () ) ; 
  }
  //
  // (3) orbital barrier 
  if ( 0 != L () )
  {
    const double q = std::abs ( m_ps2.q1_s ( s ) ) ;
    result *= 0 < m_qs ? std::pow ( q / m_qs  , 2 * L () ) : std::pow ( q , 2 * L () ) ;
  }
  //
  return result ;
}
// ============================================================================
// the second main method: term to the denominator 
// ============================================================================
std::complex<double> 
Ostap::Math::Channel0::D 
( const double s  , 
  const double m0 ) const
{
  //
  std::complex<double> result = g2 () ;
  //
  // (1) phase space factor 
  result *= m_ps2.rho1_s ( s ) ;
  //
  // (2) formfactor 
  const FormFactor* F   = formfactor () ;
  if ( nullptr != F ) 
  { 
    const double x = std::sqrt ( s ) ;
    result *= (*F)( x , m0 , m1 () , m2 () ) ; 
  }
  //
  // (3) orbital barrier 
  if ( 0 != L () )
  {
    const double q = std::abs ( m_ps2.q1_s ( s ) ) ;
    result *= 0 < m_qs ? std::pow ( q / m_qs  , 2 * L() ) : std::pow ( q , 2 *  L () ) ;
  }
  //
  return result ;
}
// ============================================================================
/*  get the phase space factor  \f$ \varrho(s) \f$
 *  - optionally normalized at the point \f$ m_n \f$ 
 *  - optionally nomalized  at \f$ m_q = m(q_s) \f$
 */
// ============================================================================
double Ostap::Math::Channel0::rho_s 
( const double s  , 
  const double mn ) const
{ 
  const double r = m_ps2.rho_s ( s ) ;
  return 
    mn > ps2().threshold () ? r / ps2 ().rho ( mn                 ) :
    0  < m_qs               ? r / ps2 ().rho ( ps2().q2m ( m_qs ) ) : r ;
}
// ============================================================================
// unique tag for this lineshape 
// ============================================================================
std::size_t Ostap::Math::Channel0::tag() const
{ 
  static const std::string s_name = "Channel0" ;
  return Ostap::Utils::hash_combiner 
    ( s_name , Ostap::Math::Channel::tag() , qs () ) ; 
}
// ============================================================================
// describe the channel 
// ============================================================================
std::string Ostap::Math::Channel0::describe() const 
{
  return 
    "Channel0(" + std::to_string ( gamma0 () ) + 
    ","         + std::to_string ( m1     () ) + 
    ","         + std::to_string ( m2     () ) + 
    ( formfactor () ? ( "," + formfactor()->describe() ) : "" )  + 
    ","         + std::to_string ( qs     () ) + ")";
}
// ============================================================================


// ============================================================================
// constructor
// ============================================================================
Ostap::Math::BW::BW
( const double                  m0      ,
  const Ostap::Math::ChannelBW& channel ,
  const double                  scale   )
  : m_m0         ( std::abs ( m0    )  )
  , m_threshold  ( std::sqrt ( channel.s_threshold  () ) )  //  cache it! 
  , m_scale      ( std::abs ( scale )  )
  , m_channels   ( )    
  , m_workspace  ( 10000 )
{
  add ( channel ) ;
} 
// ============================================================================
// protected default (empty) constructor 
// ============================================================================
Ostap::Math::BW::BW
( const double  m0    , 
  const double  scale )
  : m_m0         ( std::abs ( m0    ) )
  , m_threshold  ( 0 )  //  cache it! 
  , m_scale      ( std::abs ( scale ) ) 
  , m_channels   (   ) 
  , m_workspace  ( 10000 )
{}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Math::BW::BW
( const Ostap::Math::BW& bw ) 
  : m_m0         ( bw.m_m0        ) 
  , m_threshold  ( bw.m_threshold )
  , m_channels   ( ) 
  , m_scale      ( bw.m_scale     )
  , m_workspace  ( 10000 )
{
  for ( const auto& c : bw.m_channels ) { add ( *c ) ; }
}
// ============================================================================
/// destructor 
// ============================================================================
Ostap::Math::BW::~BW(){}
// ============================================================================
//  calculate the Breit-Wigner amplitude
// ============================================================================
std::complex<double>
Ostap::Math::BW::amplitude ( const double x ) const
{
  const double s   = x * x ;
  const double m_0 = m0 () ;
  std::complex<double> d = m_0 * m_0 - s ;
  for ( const auto& c : m_channels ) { d -= s_j * c->D ( s , m_0 ) ; }
  //
  return 1.0 / d ;
}
// ============================================================================
// unique tag for this lineshape 
// ============================================================================
std::size_t Ostap::Math::BW::tag() const
{ 
  static const std::string s_name { "BW" } ;
  std::size_t seed = Ostap::Utils::hash_combiner ( s_name , m0 () , m_threshold , m_scale ) ;
  for ( const auto&c : m_channels ) 
  { seed = Ostap::Utils::hash_combiner ( s_name , seed , c->tag() ) ; }
  return seed ;
}
// ============================================================================
// get factor \f$ \varrho(s,m_n^2) \f$ from the main channel 
// ============================================================================
double Ostap::Math::BW::rho_s ( const double s ) const
{
  if ( s <= s_threshold () ) { return 0 ; }
  /// get the main/first channel 
  const ChannelBW* c    = channel () ;
  /// get rho-factor 
  return c -> rho_s ( s , threshold () < m0() ? m0 () : 1.2 * threshold () ) ;
}
// ============================================================================
/*  Get Breit-Wigner lineshape in channel \f$ a\f$ : 
 *  \f[ F_a(m) = 2m \varrho(s) N_a(s,m_0) 
 *    \frac{\Gamma_{tot}}{\Gamma_{0,a}} \left| \mathcal{A}  \right|^2 \f] 
 *  @param x the mass point 
 *  @param A the amplitide at this point 
 */
// ============================================================================
double Ostap::Math::BW::breit_wigner 
( const double                m ,
  const std::complex<double>& A ) const 
{
  if ( m <= threshold () ) { return 0 ; }
  //
  const double s = m * m ;
  //
  /// get factor \f$ \varrho(s,m_n)\f$ 
  const double rhos = rho_s ( s ) ;
  if ( rhos <= 0 ) { return 0 ; }
  //
  /// get factor \f$ N^2(s,m_0) \f$ 
  const double     n2   = N2 ( s ) ;
  if ( n2   <= 0 ) { return 0 ; }
  //
  // rescale the overall  normalization   facror
  const double gs = 1 < m_channels.size () ? gamma () / channel()->gamma0 () : 1.0 ;
  //
  return 2 * m * n2 * std::norm ( A ) * rhos * gs / M_PI * m_scale ;
  //
}
// ============================================================================
// get the nominal total width at the pole 
// ============================================================================
double Ostap::Math::BW::gamma () const
{
  double gamtot = 0 ;
  for ( const auto& c : m_channels ) { gamtot += c->gamma0 () ; }
  return gamtot ;
}
// ============================================================================
bool Ostap::Math::BW::setM0     ( const double x )
{
  const double v       = std::abs ( x ) ;
  if ( s_equal ( v , m_m0 ) ) { return false ; } // RETURN
  m_m0   = v ;
  return true ;
}
// ============================================================================
bool Ostap::Math::BW::setScale ( const double x )
{
  const double v       = std::abs ( x ) ;
  if ( s_equal ( v , m_scale ) ) { return false ; } // RETURN
  m_scale   = v ;
  return true ;
}
// ============================================================================
bool Ostap::Math::BW::setGamma ( const double value )
{
  const double v = std::abs ( value ) ;
  const double g = gamma  () ;
  if ( s_equal ( v , g  ) ) { return false ; } // RETURN
  const double scale = v / g ;
  for ( auto& c : m_channels ) { c->setGamma0 ( c->gamma0 () *  scale ) ; }
  return true ;
}

// ============================================================================
// get the integral between low and high limits
// ============================================================================
double  Ostap::Math::BW::integral
( const double low  ,
  const double high ) const
{
  if ( s_equal ( low , high ) ) { return                 0.0 ; } // RETURN
  if (           low > high   ) { return - integral ( high ,
                                                      low  ) ; } // RETURN
  //
  const double t = threshold() ;
  if ( t >= high ) { return                    0   ; }
  if ( t >  low  ) { return integral  ( t , high ) ; }
  //
  // split into reasonable sub intervals according to gamma-parameters 
  double g0 = gamma () ;
  if ( 0 < m_m0 ) 
  {
    // check the imaginary part at the pole and use it as "width"
    const double gX = std::abs ( std::imag ( 1.0 / amplitude ( m_m0 ) ) / m_m0 ) ;
    if ( 0 < gX ) { g0 = gX ; }
  }
  //
  if ( 0 < g0 ) 
  {
    //
    for ( unsigned short i = 1 ; i < 6 ; ++i ) 
    {
      const double x1 = m_m0 + i * g0 ;
      if ( low < x1 && x1 < high ) 
      { return integral ( low , x1 ) + integral ( x1 , high ) ; }
      //
      const double x2 = m_m0 - i * g0 ;
      if ( low < x2 && x2 < high ) 
      { return integral ( low , x2 ) + integral ( x2 , high ) ; } 
    }
  }
  //
  if (  low < m_m0 && m_m0 < high ) 
  { return integral ( low , m_m0 ) + integral ( m_m0 , high ) ; }
  //
  // split, if interval too large
  //
  const double width = std::max ( g0 , 0.0 ) ;
  if ( 0 < width && 10 * width < high - low  )
  {
    return
      integral ( low                   , 0.5 *  ( high + low ) ) +
      integral ( 0.5 *  ( high + low ) ,          high         ) ;
  }
  //
  //  where the tails start?
  //
  const double x_low  = m_m0 - 5 * g0 ;
  const double x_high = m_m0 + 6 * g0 ;
  //
  if ( low < x_low  && x_low  < high ) 
  { return integral ( low , x_low  ) + integral ( x_low  , high ) ; }
  if ( low < x_high && x_high < high ) 
  { return integral ( low , x_high ) + integral ( x_high , high ) ; }
  
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<BW> s_integrator {} ;
  static char s_message[] = "Integral(BW)" ;
  //
  const bool in_tail = high <= x_low || x_high <= low ;
  //

  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag () , 
      &F     , 
      low    , high       ,          // low & high edges
      workspace ( m_workspace ) ,    // workspace      
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // absolute precision
      m_workspace.size () ,          // size of workspace
      s_message           ,          // message  
      __FILE__ , __LINE__ ,          // file&line 
      GSL_INTEG_GAUSS61   ) ;        // integration rule 
  //
  return result ;
}
// ============================================================================
// get the integral b
// ============================================================================
double  Ostap::Math::BW::integral () const
{
  //
  // split into reasonable sub intervals
  //
  const double t      = threshold () ;
  const double g0     = gamma     () ;
  //
  const double x1     = std::max ( t , m_m0 - 10 * g0 ) ;
  const double x2     = std::max ( t , m_m0 + 10 * g0 ) ;
  const double x_high = std::max ( x1 , x2 ) ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<BW> s_integrator {} ;
  static char s_message[] = "Integral(BW/tail)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.qagiu_integrate
    ( tag () , 
      &F     , 
      x_high               ,          // low edge
      workspace ( m_workspace ) ,    // workspace
      s_APRECISION_TAIL    ,          // absolute precision
      s_RPRECISION_TAIL    ,          // relative precision
      m_workspace.size ()  ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result + integral ( t  , x_high );
}
// ============================================================================
void Ostap::Math::BW::add 
( const Ostap::Math::ChannelBW& c ) 
{ 
  if ( m_channels.empty() ) { m_threshold = std::sqrt ( c.s_threshold() ) ; }
  m_channels.push_back ( std::unique_ptr<ChannelBW> (  c.clone() ) ) ;
}
// ============================================================================


// ============================================================================
// Breit-Wigner iself 
// ============================================================================
/*  constructor from all parameters
 *  @param m0   the pole position 
 *  @param gam0 the nominal width at the pole 
 *  @param m1   the mass of the 1st daughter particle 
 *  @param m2   the mass of the 2nd daughter particle 
 *  @param L    the orbital momentum 
 */
// ============================================================================
Ostap::Math::BreitWigner::BreitWigner
( const double         m0     ,
  const double         gam0   ,
  const double         m1     ,
  const double         m2     ,
  const unsigned short L      , 
  const double         scale  )
  : BW ( m0 , Ostap::Math::Channel ( gam0 , m1 , m2 , L ) , scale )
{}
// ============================================================================
/*  constructor from all parameters
 *  @param m0   the pole position 
 *  @param gam0 the nominal width at the pole 
 *  @param m1   the mass of the 1st daughter particle 
 *  @param m2   the mass of the 2nd daughter particle 
 *  @param L    the orbital momentum 
 *  @param F    the Jackson's formfactor 
 */
// ============================================================================
Ostap::Math::BreitWigner::BreitWigner
( const double                               m0     , 
  const double                               gam0   ,
  const double                               m1     ,
  const double                               m2     ,
  const unsigned short                       L      ,
  const Ostap::Math::FormFactors::JacksonRho F      , 
  const double                               scale  )
  : BW ( m0 , Ostap::Math::Channel ( gam0 , m1 , m2 , L , F ) , scale )
{}
// ============================================================================
/*  constructor from all parameters
 *  @param m0   the pole position 
 *  @param gam0 the nominal width at the pole 
 *  @param m1   the mass of the 1st daughter particle 
 *  @param m2   the mass of the 2nd daughter particle 
 *  @param L    the orbital momentum 
 *  @param F    the formfactor 
 */
// ============================================================================
Ostap::Math::BreitWigner::BreitWigner
( const double                               m0    , 
  const double                               gam0  ,
  const double                               m1    ,
  const double                               m2    ,
  const unsigned short                       L     ,
  const Ostap::Math::FormFactor&             F     ,  
  const double                               scale )
  : BW ( m0 , Ostap::Math::Channel ( gam0 , m1 , m2 , L , F ) , scale )
{}
// ============================================================================
/*  constructor from all parameters
 *  @param m0   the pole position 
 *  @param c    the channel 
 */
// ============================================================================
Ostap::Math::BreitWigner::BreitWigner
( const double                  m0    , 
  const Ostap::Math::ChannelBW& c     ,  
  const double                  scale )
  : BW ( m0 , c , scale )
{}
// ============================================================================
/// copy constructor 
// ============================================================================
Ostap::Math::BreitWigner::BreitWigner  
( const Ostap::Math::BreitWigner& bw ) 
  : BW ( bw ) 
{}
// ============================================================================

// ============================================================================
// clone it 
// ============================================================================
Ostap::Math::BreitWigner*
Ostap::Math::BreitWigner::clone() const 
{ return new Ostap::Math::BreitWigner ( *this ) ; }
// ============================================================================


// ============================================================================
// multi-channel version of Breit-Wigner
// ===========================================================================
/*  constructor from single channel
 *  @param m0   the pole position 
 *  @param c1   the 1st channel 
 */
// ===========================================================================
Ostap::Math::BreitWignerMC::BreitWignerMC 
( const double                  m0 ,
  const Ostap::Math::ChannelBW& c1 )
  : BW ( m0 , c1 )
{}
// ============================================================================
// clone it 
// ============================================================================
Ostap::Math::BreitWignerMC*
Ostap::Math::BreitWignerMC::clone() const 
{ return new  Ostap::Math::BreitWignerMC (*this )  ;}
// ============================================================================

   





// ============================================================================
// constructor from all parameters
// ============================================================================
Ostap::Math::Rho0::Rho0
( const double m0       ,
  const double gam0     ,
  const double pi_mass  ,
  const double scale    ) 
  : Ostap::Math::BreitWigner ( m0         ,
                               gam0       ,
                               pi_mass    ,
                               pi_mass    ,
                               1          ,
                               Ostap::Math::FormFactors::Jackson_A7 , 
                               scale      )
  , m_m1 ( std::abs ( pi_mass ) )  
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Rho0::~Rho0(){}

// ============================================================================
// constructor from all parameters
// ============================================================================
Ostap::Math::Kstar0::Kstar0
( const double m0       ,
  const double gam0     ,
  const double k_mass   ,
  const double pi_mass  ,
  const double scale    ) 
  : Ostap::Math::BreitWigner ( m0         ,
                               gam0       ,
                               k_mass     ,
                               pi_mass    ,
                               1          ,
                               Ostap::Math::FormFactors::Jackson_A2 , 
                               scale      )
  , m_m1 ( std::abs ( k_mass  ) )  
  , m_m2 ( std::abs ( pi_mass ) )  
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Kstar0::~Kstar0(){}

// ============================================================================
// constructor from all parameters
// ============================================================================
Ostap::Math::Phi0::Phi0
( const double m0       ,
  const double gam0     ,
  const double k_mass   ,
  const double scale    ) 
  : Ostap::Math::BreitWigner ( m0         ,
                               gam0       ,
                               k_mass     ,
                               k_mass     ,
                               1          ,
                               Ostap::Math::FormFactors::Jackson_A2 , 
                               scale      )
  , m_m1 ( std::abs ( k_mass  ) )  
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Phi0::~Phi0(){}


// // calculate the function
// // ============================================================================
// double Ostap::Math::Rho0FromEtaPrime::operator() ( const double x ) const
// {
//   //
//   if ( m_eta_prime <= x ) { return 0 ; }
//   //
//   const double k_gamma = Ostap::Math::PhaseSpace2::q ( m_eta_prime , x , 0 ) ;
//   if ( 0 >= k_gamma     ) { return 0 ; }
//   //
//   const double rho     = breit_wigner ( x ) ;
//   if ( 0 >= rho         ) { return 0 ; }
//   //
//   return rho * Ostap::Math::POW ( 2 * k_gamma / m_eta_prime , 3 ) * 20 ;
//   //
// }
// // ============================================================================

// ============================================================================
// Channel-Flatte
// ============================================================================
/** constructor from all parameters 
 *  @param g2    thew squared coupling constantt  
 *  @param m1    the mass of the 1st daughter
 *  @param m2    the mass of the 2nd daughter
 */
// ============================================================================
Ostap::Math::ChannelFlatte::ChannelFlatte
( const double g2 , 
  const double m1 , 
  const double m2 )
  : ChannelCW ( g2 , m1 , m2 ) 
{}
// ============================================================================
// clone it!
// ============================================================================
Ostap::Math::ChannelFlatte*
Ostap::Math::ChannelFlatte::clone() const 
{ return new Ostap::Math::ChannelFlatte ( *this ) ; }
// ============================================================================
/*  the first main method: numerator
 * \f[ N(s,m_0) = \frac{1}{\pi}m_0 g \f] 
 */
// ============================================================================
double Ostap::Math::ChannelFlatte::N2
( const double /* s */ , 
  const double    m0   ) const { return m0 * g2 () ; }
// ============================================================================
// the second main method: term to the denominator 
// ============================================================================
std::complex<double> 
Ostap::Math::ChannelFlatte::D   
( const double s  , 
  const double m0 ) const 
{
  // return m0 * g2 () * ps2().q1_s ( s )  ;
  const double sqs = std::sqrt ( s ) ;
  return m0 * g2 () * ( 2 / sqs ) * ps2().q1_s ( s ) ;
}
// ============================================================================
// get unique tag/label 
// ============================================================================
std::size_t Ostap::Math::ChannelFlatte::tag() const
{ 
  static const std::string s_name = "ChannelFlatte" ;
  return Ostap::Utils::hash_combiner
    ( s_name    , Ostap::Math::ChannelCW::tag () ) ;
}
// ============================================================================
// describe the channel 
// ============================================================================
std::string Ostap::Math::ChannelFlatte::describe() const 
{
  return 
    "ChannelFlatte(" + std::to_string ( gamma0 () ) + 
    ","          + std::to_string ( m1     () ) + 
    ","          + std::to_string ( m2     () ) + ")" ;
}
// ============================================================================


// ============================================================================
// Channel-Flatte
// ============================================================================
/** constructor from all parameters 
 *  @param g2    thew squared coupling constantt  
 *  @param m1    the mass of the 1st daughter
 *  @param m2    the mass of the 2nd daughter
 */
// ============================================================================
Ostap::Math::ChannelFlatteBugg::ChannelFlatteBugg
( const double g        , 
  const double mcharged ,   // GeV/c2  
  const double mneutral ,   // GeV/c2 
  const double mK       ,   // GeV/c2 
  const double alpha    ,   // GeV^{-2} form-factor 
  const double fc       ,   // isospin factor  
  const double fn       )   // isospin factor  
  : ChannelFlatte ( g , mcharged , mcharged ) 
  , m_alpha    ( std::abs ( alpha    ) )
  , m_fc       ( std::abs ( fc ) ) 
  , m_fn       ( std::abs ( fn ) )
  , m_ps2n     ( mneutral , mneutral ) 
  , m_ps2k     ( mK       , mK       )
{}
// ============================================================================
// clone it!
// ============================================================================
Ostap::Math::ChannelFlatteBugg*
Ostap::Math::ChannelFlatteBugg::clone() const 
{ return new Ostap::Math::ChannelFlatteBugg ( *this ) ; }
// ============================================================================
// the second main method: term to the denominator 
// ============================================================================
std::complex<double> 
Ostap::Math::ChannelFlatteBugg::D   
( const double s  , 
  const double m0 ) const 
{
  const double sqs = std::sqrt ( std::abs ( s ) ) ;
  //
  // formfactor squared 
  const double ff2 = 
    ( m_ps2k.s_threshold() < s && !s_zero ( m_alpha ) ) ? 
    std::exp ( -2 * m_alpha * std::pow ( m_ps2k.q_s ( s ) , 2 ) ) : 1.0 ;  
  //
  return ( m0 * g2 () * 2.0 * ff2 / sqs ) 
    //        charged                    neutral  
    * ( m_fc * ps2().q1_s  ( s ) + m_fn * m_ps2n.q1_s ( s ) ) ;
}
// ============================================================================
// get unique tag/label 
// ============================================================================
std::size_t Ostap::Math::ChannelFlatteBugg::tag() const
{ 
  static const std::string s_name = "ChannelFlatteBugg" ;
  return Ostap::Utils::hash_combiner 
    ( s_name       ,
      Ostap::Math::ChannelFlatte::tag () , 
      m_alpha           , 
      m_fc              , 
      m_fn              ,
      m_ps2n.tag()      ,
      m_ps2k.tag()      ) ;
}
// ============================================================================
// describe the channel 
// ============================================================================
std::string Ostap::Math::ChannelFlatteBugg::describe() const 
{
  return 
    "ChannelFlatteBugg(" + std::to_string ( gamma0 () ) + 
    ","          + std::to_string ( m1        () ) + 
    ","          + std::to_string ( mneutral  () ) + 
    ","          + std::to_string ( m_ps2k.m1 () ) + 
    ","          + std::to_string ( m_alpha      ) + 
    ","          + std::to_string ( m_fc         ) + 
    ","          + std::to_string ( m_fn         ) + ")" ;
}
// ============================================================================


// ============================================================================
//               Flatte
// ============================================================================
/*  constructor  from three parameters
 *  @param m0    the mass
 *  @param m0g1  parameter \f$ m_0\times g_1\f$
 *  @param g2og2 parameter \f$ g2/g_1       \f$
 */
// ============================================================================
Ostap::Math::Flatte::Flatte
( const double m0    ,
  const double m0g1  ,
  const double g2og1 ,
  const double mA1   ,
  const double mA2   ,
  const double mB1   ,
  const double mB2   ,
  const double g0    ,
  const double scale )
  : BW ( m0 , ChannelFlatte ( 1 , mA1 , mA2 ) , scale ) 
  , m_A1 ( std::abs ( mA1 ) ) 
  , m_A2 ( std::abs ( mA2 ) ) 
  , m_B1 ( std::abs ( mB1 ) ) 
  , m_B2 ( std::abs ( mB2 ) )
{
  add ( ChannelFlatte ( 1  , mB1 , mB2 ) ) ;
  add ( ChannelCW     ( g0 , mA1 , mA1 ) ) ;
  //
  const double g1 = std::abs ( m0g1    / m0 ) ;
  const double g2 = std::abs ( g2og1 ) * g1   ;
  //
  Ostap::Assert ( 3 ==   nChannels()            , 
                  "Invalid number of channels!" , 
                  "Ostap:Math::Flatte::Flatte"  ) ;
  //
  m_channels[0] -> setGamma0 ( g1 ) ;
  m_channels[1] -> setGamma0 ( g2 ) ;
  m_channels[2] -> setGamma0 ( g0 ) ;
  //
}
// ============================================================================
// copy contructor
// ============================================================================
Ostap::Math::Flatte::Flatte 
( const Ostap::Math::Flatte& fl ) 
  : BW ( fl ) 
{}
// ============================================================================
Ostap::Math::Flatte*
Ostap::Math::Flatte::clone() const { return new Ostap::Math::Flatte ( *this ) ; }
// ============================================================================
// unique tag
// ============================================================================
std::size_t Ostap::Math::Flatte::tag () const
{ 
  static const std::string s_name = "Flatte" ;
  return Ostap::Utils::hash_combiner 
    ( s_name                  ,
      Ostap::Math::BW::tag () , 
      m_A1                    , 
      m_A2                    , 
      m_B1                    , 
      m_B2                    ) ;
}

// ============================================================================
/*  constructor from all parameters
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
// ============================================================================
Ostap::Math::FlatteBugg::FlatteBugg
( const double m0      , // GeV/c2
  const double g1      , // GeV/c2 
  const double g2og1   , // dimensionless  
  const double alpha   , // GeV^-2 
  const double mpiplus , // pi+ mass in GeV 
  const double mpizero , // pi0 mass in GeV 
  const double mKplus  , // K+ mass in GeV 
  const double mKzero  , // K0 mass in GeV 
  const double g0      , // the constant width for "other" decays
  const double scale   )
  : BW ( m0 , ChannelFlatteBugg ( g1 , mpiplus , mpizero , mKplus , 0.0 , 2.0/3.0 , 1.0/3.0 ) , scale ) 
  , m_alpha   ( std::abs ( alpha ) )
  , m_mpiplus ( std::abs ( mpiplus ) ) 
  , m_mpizero ( std::abs ( mpizero ) ) 
  , m_mKplus  ( std::abs ( mKplus  ) ) 
  , m_mKzero  ( std::abs ( mKzero  ) )
{
  const double g2 = std::abs ( g2og1 ) * std::abs ( g1 ) ;
  //
  add ( ChannelFlatteBugg ( g2 , mKplus  , mKzero , mKplus , alpha , 0.5 , 0.5 ) ) ;
  add ( ChannelCW         ( g0 , mpiplus , mpiplus ) ) ;
  //
  //
  Ostap::Assert ( 3 ==   nChannels()                , 
                  "Invalid number of channels!"     , 
                  "Ostap:Math::Flatte::FlatteBugg"  ) ;
  //
  m_channels[0] -> setGamma0 ( g1 ) ;
  m_channels[1] -> setGamma0 ( g2 ) ;
  m_channels[2] -> setGamma0 ( g0 ) ;
  //
}
// ============================================================================
// copy contructor
// ============================================================================
Ostap::Math::FlatteBugg::FlatteBugg 
( const Ostap::Math::FlatteBugg& fl ) 
  : BW ( fl ) 
{}
// ============================================================================
Ostap::Math::FlatteBugg*
Ostap::Math::FlatteBugg::clone() const 
{ return new Ostap::Math::FlatteBugg ( *this ) ; }
// ============================================================================
// unique tag
// ============================================================================
std::size_t Ostap::Math::FlatteBugg::tag () const
{ 
  static const std::string s_name = "FlatteBugg" ;
  return Ostap::Utils::hash_combiner 
    ( s_name                  , 
      Ostap::Math::BW::tag () ,
      m_alpha                 , 
      m_mpiplus               , 
      m_mpizero               , 
      m_mKplus                , 
      m_mKzero                ) ;
}

// // ============================================================================
// // Bugg
// // ============================================================================
// /*  constructor from all masses and angular momenta
//  *  @param M  mass of sigma (very different from the pole positon!)
//  *  @param g2 width parameter g2 (4pi width)
//  *  @param b1 width parameter b1  (2pi coupling)
//  *  @param b2 width parameter b2  (2pi coupling)
//  *  @param s1 width parameter s1  (cut-off for 4pi coupling)
//  *  @param s2 width parameter s2  (cut-off for 4pi coupling)
//  *  @param a  parameter a (the exponential cut-off)
//  *  @param m1 the mass of the first  particle
//  */
// // ============================================================================
// Ostap::Math::Bugg::Bugg
// ( const double         M  ,
//   const double         g2 ,
//   const double         b1 ,
//   const double         b2 ,
//   const double         a  ,
//   const double         s1 ,
//   const double         s2 ,
//   const double         m1 )
// //
//   : m_M  ( std::abs ( M  ) )
//   , m_g2 ( std::abs ( g2 ) )
//   , m_b1 ( std::abs ( b1 ) )
//   , m_b2 ( std::abs ( b2 ) )
//   , m_s1 ( std::abs ( s1 ) )
//   , m_s2 ( std::abs ( s2 ) )
//   , m_a  ( std::abs ( a  ) )
// // phase space
//   , m_ps ( m1 , m1 )
// //
//   , m_workspace ()
// {}
// // ============================================================================
// // destructor
// // ============================================================================
// Ostap::Math::Bugg::~Bugg(){}
// // ============================================================================
// double Ostap::Math::Bugg::rho2_ratio ( const double x ) const
// {
//   if ( lowEdge() >= x ) { return 0 ; }
//   //
//   return
//     Ostap::Math::PhaseSpace2::phasespace ( x    , m1() , m2 () ) /
//     Ostap::Math::PhaseSpace2::phasespace ( M () , m1() , m2 () ) ;
// }
// // ============================================================================
// std::complex<double>
// Ostap::Math::Bugg::rho4_ratio ( const double x ) const
// {
//   //
//   if ( 2 * m1() >= x ) { return 0 ; }
//   //
//   return rho4 ( x ) / rho4 ( M() ) ;
// }
// // ============================================================================
// std::complex<double>
// Ostap::Math::Bugg::rho4 ( const double x ) const
// {
//   const double s  = x * x ;
//   //
//   const double r2 = 1 - 16 * m1() * m1() / s ;
//   //
//   const double r  =
//     std::sqrt ( std::abs ( r2 ) ) *
//     ( 1 + std::exp ( ( s1 () - s )  / s2 () ) ) ;
//   //
//   return 0 <= r2 ?
//     std::complex<double> ( r , 0 ) :
//     std::complex<double> ( 0 , r ) ;
// }
// // ============================================================================
// // Adler's pole
// // ============================================================================
// double Ostap::Math::Bugg::adler ( const double x ) const
// {
//   if ( lowEdge() >= x ) { return 0 ; }
//   //
//   const double pole = 0.5 * m1 () * m1 ()  ;
//   //
//   return ( x * x - pole ) / ( M2 () - pole ) ;
// }
// // ============================================================================
// // get the running width by Bugg
// // ============================================================================
// std::complex<double>
// Ostap::Math::Bugg::gamma ( const double x ) const
// {
//   //
//   if ( lowEdge() >= x ) { return 0 ; }
//   //
//   const double s = x * x ;
//   //
//   const double g1 =
//     b     ( x ) *
//     adler ( x ) * std::exp ( -1 * ( s - M2() )  / a() ) ;
//   //
//   return g1 * rho2_ratio ( x ) + g2 () * rho4_ratio ( x ) ;
// }
// // ============================================================================
// // get the amlitude  (not normalized!)
// // ============================================================================
// std::complex<double>
// Ostap::Math::Bugg::amplitude (  const double x ) const
// {
//   if ( lowEdge() >= x ) { return 0 ; }
//   //
//   static const std::complex<double> j ( 0 , 1 ) ;
//   //
//   std::complex<double> d = M2() - x * x  - j * M() * gamma ( x ) ;
//   //
//   return 1.0 / d ;
// }
// // ============================================================================
// // evaluate Bugg
// // ============================================================================
// double Ostap::Math::Bugg::pdf ( const double x ) const
// {
//   //
//   if ( lowEdge() >= x ) { return 0 ; }
//   //
//   const double result = phaseSpace  ( x ) ;
//   if ( 0 >= result ) { return 0 ; }
//   //
//   return result * std::norm ( amplitude ( x ) ) ;
// }
// // ============================================================================
// // set the proper parameters
// // ============================================================================
// bool Ostap::Math::Bugg::setM ( const double x )
// {
//   //
//   const double v = std::abs ( x ) ;
//   if ( s_equal ( v , m_M ) ) { return false ; }
//   //
//   m_M = v ;
//   //
//   return true ;
// }
// // ============================================================================
// // set the proper parameters
// // ============================================================================
// bool Ostap::Math::Bugg::setG2 ( const double x )
// {
//   //
//   const double v = std::abs ( x ) ;
//   if ( s_equal ( v , m_g2 ) ) { return false ; }
//   //
//   m_g2 = v ;
//   //
//   return true ;
// }
// // ============================================================================
// // set the proper parameters
// // ============================================================================
// bool Ostap::Math::Bugg::setB1 ( const double x )
// {
//   //
//   const double v = std::abs ( x ) ;
//   if ( s_equal ( v , m_b1 ) ) { return false ; }
//   //
//   m_b1 = v ;
//   //
//   return true ;
// }
// // ============================================================================
// // set the proper parameters
// // ============================================================================
// bool Ostap::Math::Bugg::setB2 ( const double x )
// {
//   //
//   const double v = std::abs ( x ) ;
//   if ( s_equal ( v , m_b2 ) ) { return false ; }
//   //
//   m_b2 = v ;
//   //
//   return true ;
// }
// // ============================================================================
// // set the proper parameters
// // ============================================================================
// bool Ostap::Math::Bugg::setS1 ( const double x )
// {
//   //
//   const double v = std::abs ( x ) ;
//   if ( s_equal ( v , m_s1 ) ) { return false ; }
//   //
//   m_s1 = v ;
//   //
//   return true ;
// }
// // ============================================================================
// // set the proper parameters
// // ============================================================================
// bool Ostap::Math::Bugg::setS2 ( const double x )
// {
//   //
//   const double v = std::abs ( x ) ;
//   if ( s_equal ( v , m_s2 ) ) { return false ; }
//   //
//   m_s2 = v ;
//   //
//   return true ;
// }
// // ============================================================================
// // set the proper parameters
// // ============================================================================
// bool Ostap::Math::Bugg::setA ( const double x )
// {
//   //
//   const double v = std::abs ( x ) ;
//   if ( s_equal ( v , m_a ) ) { return false ; }
//   //
//   m_a = v ;
//   //
//   return true ;
// }
// // ============================================================================
// // get the integral between low and high limits
// // ============================================================================
// double  Ostap::Math::Bugg::integral
// ( const double low  ,
//   const double high ) const
// {
//   if ( s_equal ( low , high ) ) { return                 0.0 ; } // RETURN
//   if (           low > high   ) { return - integral ( high ,
//                                                       low  ) ; } // RETURN
//   //
//   if ( high <= lowEdge  () ) { return 0 ; }
//   //
//   if ( low  <  lowEdge  () )
//   { return integral ( lowEdge() , high        ) ; }
//   //
//   // use GSL to evaluate the integral
//   //
//   static const Ostap::Math::GSL::Integrator1D<Bugg> s_integrator {} ;
//   static char s_message[] = "Integral(Bugg)" ;
//   //
//   const auto F = s_integrator.make_function ( this ) ;
//   int    ierror   =  0 ;
//   double result   =  1 ;
//   double error    = -1 ;
//   std::tie ( ierror , result , error ) = s_integrator.qag_integrate
//     ( &F , 
//       low , high          ,          // low & high edges
//       workspace ( m_workspace ) ,    // workspace      
//       s_PRECISION         ,          // absolute precision
//       s_PRECISION         ,          // relative precision
//       m_workspace.size () ,          // size of workspace
//       s_message           , 
//       __FILE__ , __LINE__ ) ;
//   //
//   return result ;
// }
// // ============================================================================


// // ============================================================================
// // SWANSON CUSP 
// // ============================================================================
// // constructor
// // ============================================================================
// Ostap::Math::Swanson::Swanson
// ( const double         m1     ,   // the first  real particle 
//   const double         m2     ,   // the second real particle                
//   const double         m1_0   ,   // the first  particle for cusp
//   const double         m2_0   ,   // the second particle for cusp 
//   const double         beta_0 ,   // beta_0 parameter
//   const unsigned short L      )   // orbital momentum for real particles 
//   : m_bw ( ( std::abs ( m1 )  + std::abs ( m2 ) ) * 2.1 , // almost arbitrary 
//            ( std::abs ( m1 )  + std::abs ( m2 ) ) * 0.5 , // almost arbitrary  
//            std::abs ( m1 ) , 
//            std::abs ( m2 ) ,  
//            L               ) 
//   , m_m1         ( std::abs (   m1_0 ) )
//   , m_m2         ( std::abs (   m2_0 ) )
//   , m_beta0      ( std::abs ( beta_0 ) )
//     //
//   , m_workspace  ()
//     //
// {}
// // ============================================================================
// // constructor
// // ============================================================================
// Ostap::Math::Swanson::Swanson
// ( const double         m1             ,   // the first  real particle 
//   const double         m2             ,   // the second real particle                
//   const double         m1_0           ,   // the first  particle for cusp
//   const double         m2_0           ,   // the second particle for cusp 
//   const double         beta_0         ,   // beta_0 parameter
//   const unsigned short L              ,   // orbital momentum for real particles 
//   const Ostap::Math::FormFactors::JacksonRho  r )  //  formfactor
//   : m_bw ( ( std::abs ( m1 )  + std::abs ( m2 ) ) * 2.1 , // almost arbitrary 
//            ( std::abs ( m1 )  + std::abs ( m2 ) ) * 0.5 , // almost arbitrary  
//            std::abs ( m1 ) , 
//            std::abs ( m2 ) ,  
//            L  , r          )            
//   , m_m1         ( std::abs (   m1_0 ) )
//   , m_m2         ( std::abs (   m2_0 ) )
//   , m_beta0      ( std::abs ( beta_0 ) )
//     //
//   , m_workspace  ()
//     //
// {}

// // ============================================================================
// // constructor
// // ============================================================================
// Ostap::Math::Swanson::Swanson
// ( const double         m1             ,   // the first  real particle 
//   const double         m2             ,   // the second real particle                
//   const double         m1_0           ,   // the first  particle for cusp
//   const double         m2_0           ,   // the second particle for cusp 
//   const double         beta_0         ,   // beta_0 parameter
//   const unsigned short L              ,   // orbital momentum for real particles 
//   const Ostap::Math::FormFactor&   f  )  //  formfactor
//   : m_bw ( ( std::abs ( m1 )  + std::abs ( m2 ) ) * 2.1 , // almost arbitrary 
//            ( std::abs ( m1 )  + std::abs ( m2 ) ) * 0.5 , // almost arbitrary  
//            std::abs ( m1 ) , 
//            std::abs ( m2 ) ,  
//            L  , f          )            
//   , m_m1         ( std::abs (   m1_0 ) )
//   , m_m2         ( std::abs (   m2_0 ) )
//   , m_beta0      ( std::abs ( beta_0 ) )
//     //
//   , m_workspace  ()
//     //
// {}


// // ============================================================================
// // constructor
// // ============================================================================
// Ostap::Math::Swanson::Swanson
// ( const Ostap::Math::BreitWigner&   bw             ,   // breit-wigner 
//   const double         m1_0   ,   // the first  particle for cusp
//   const double         m2_0   ,   // the second particle for cusp 
//   const double         beta_0 )   // beta_0 parameter
//   : m_bw     ( bw ) 
//   , m_m1     ( std::abs (   m1_0 ) )
//   , m_m2     ( std::abs (   m2_0 ) )
//   , m_beta0  ( std::abs ( beta_0 ) )
//     //
//   , m_workspace  ()
//     //
// {}
// // ============================================================================
// // copy constructor 
// // ============================================================================
// Ostap::Math::Swanson::Swanson
// ( const Ostap::Math::Swanson& sw ) 
//   : m_bw    ( sw.m_bw    )
//   , m_m1    ( sw.m_m1    )
//   , m_m2    ( sw.m_m2    )
//   , m_beta0 ( sw.m_beta0 )
//     //
//   , m_workspace  ()
//     //
// {}
// // ============================================================================
// // destructor
// // ============================================================================
// Ostap::Math::Swanson::~Swanson (){}
// // ============================================================================
// // set the proper parameters
// // ============================================================================
// bool Ostap::Math::Swanson::setM1_0 ( const double x )
// {
//   //
//   const double v = std::abs ( x ) ;
//   if ( s_equal ( v , m_m1 ) ) { return false ; }
//   //
//   m_m1 = v ;
//   //
//   return true ;
// }
// // ============================================================================
// bool Ostap::Math::Swanson::setM2_0 ( const double x )
// {
//   //
//   const double v = std::abs ( x ) ;
//   if ( s_equal ( v , m_m2 ) ) { return false ; }
//   //
//   m_m2 = v ;
//   //
//   return true ;
// }
// // ============================================================================
// bool Ostap::Math::Swanson::setBeta0( const double x )
// {
//   //
//   const double v = std::abs ( x ) ;
//   if ( s_equal ( v , m_beta0 ) ) { return false ; }
//   //
//   m_beta0 = v ;
//   //
//   return true ;
// }
// // ============================================================================
// //  calculate the Swanson amplitude
// // ============================================================================
// std::complex<double>
// Ostap::Math::Swanson::amplitude ( const double x ) const
// {
//   //
//   const double  f = - s_SQRT2PISQUAREDi*m_beta0/(1/m_m1+1/m_m2) ;
//   //
//   const double zf = 4 * m_m1 * m_m2 / ( m_beta0 * m_beta0 * ( m_m1 + m_m2 ) ) ;
//   const double z  = zf * ( m_m1 + m_m2 - x ) ;
//   //
//   // above threshold, Z is negative 
//   std::complex<double> iZ = 
//     0 <= z ? 
//     std::complex<double>(     std::sqrt (            z   ) , 0 ) :
//     std::complex<double>( 0 , std::sqrt ( std::abs ( z ) )     ) ;
//   //
//   return f * 0.5 * s_SQRTPIHALF * ( 1.0 - s_SQRTPI * iZ * Ostap::Math::erfcx ( iZ ) ) ;
// }
// // ============================================================================
// //  calculate the Swanson shape 
// // ============================================================================
// double Ostap::Math::Swanson::swanson ( const double x ) const
// {
//   if ( m_bw.m1() + m_bw.m2() >= x ) { return 0 ; }
//   //
//   const double               g = m_bw.gamma ( x ).real()  ;
//   if ( 0 >= g                     ) { return 0 ; }  
//   //
//   const std::complex<double> a = amplitude ( x ) ;
//   //
//   return 2 * x * std::norm ( a ) * g / m_bw.gamma0() / M_PI ;
// }
// // ============================================================================
// // get the integral between low and high limits
// // ============================================================================
// double  Ostap::Math::Swanson::integral
// ( const double low  ,
//   const double high ) const
// {
//   if ( s_equal ( low , high ) ) { return                 0.0 ; } // RETURN
//   if (           low > high   ) { return - integral ( high ,
//                                                       low  ) ; } // RETURN
//   //
//   const double x_min  = m_bw.channel().m1() + m_bw.channel().m2() ;
//   if ( x_min >= high ) { return                        0   ; }
//   if ( x_min >  low  ) { return integral  ( x_min , high ) ; }
//   //
//   // split into reasonable sub intervals
//   //
//   const double x1   = x_min +  1 * ( m_m1 + m_m2 ) ;
//   const double x2   = x_min +  2 * ( m_m1 + m_m2 ) ;
//   const double x5   = x_min +  5 * ( m_m1 + m_m2 ) ;
//   const double x10  = x_min + 10 * ( m_m1 + m_m2 ) ;
//   //
//   if ( low <  x1 &&  x1 < high ) { return integral ( low ,  x1 ) + integral (  x1 , high ) ; }
//   if ( low <  x2 &&  x2 < high ) { return integral ( low ,  x2 ) + integral (  x2 , high ) ; }
//   if ( low <  x5 &&  x5 < high ) { return integral ( low ,  x5 ) + integral (  x5 , high ) ; }
//   if ( low < x10 && x10 < high ) { return integral ( low , x10 ) + integral ( x10 , high ) ; }
//   //
//   // use GSL to evaluate the integral
//   //
//   static const Ostap::Math::GSL::Integrator1D<Swanson> s_integrator {} ;
//   static char s_message[] = "Integral(Swanson)" ;
//   //
//   const auto F = s_integrator.make_function ( this ) ;
//   int    ierror   =  0 ;
//   double result   =  1 ;
//   double error    = -1 ;
//   std::tie ( ierror , result , error ) = s_integrator.qag_integrate
//     ( &F , 
//       low , high          ,          // low & high edges
//       workspace ( m_workspace ) ,    // workspace      
//       s_PRECISION         ,          // absolute precision
//       ( x10  <= low  ) ? s_PRECISION_TAIL :
//       s_PRECISION         ,          // relative precision
//       m_workspace.size () ,          // size of workspace
//       s_message           , 
//       __FILE__ , __LINE__ ) ;
//   //
//   return result ;
// }
// // ============================================================================





// ============================================================================
// "2-from-3" shapes 
// ============================================================================
Ostap::Math::Channel23L::Channel23L
( const Ostap::Math::ChannelBW&     ch , 
  const Ostap::Math::PhaseSpace23L& ps ) 
  : ChannelBW ( ch ) 
  , m_channel ( ch.clone() ) 
  , m_ps      ( ps ) 
    
{}
// ============================================================================
// constructor from the channel and dalitz configuration 
// ============================================================================
Ostap::Math::Channel23L::Channel23L
( const Ostap::Math::ChannelBW&     ch , 
  const Ostap::Kinematics::Dalitz&  dp , 
  const unsigned short              L2 ) 
  : ChannelBW ( ch ) 
  , m_channel ( ch.clone() ) 
  , m_ps ( dp.m1() , dp.m2() , dp.m3() , dp.M() , L2 , 0 )
{}
// ============================================================================
// constructor from the channel and configuration 
// ============================================================================
Ostap::Math::Channel23L::Channel23L
( const Ostap::Math::ChannelCW&     ch , 
  const double                      m3 ,  
  const double                      M  ,
  const unsigned short              L2 )      
  : ChannelBW ( ch ) 
  , m_ps ( ch.m1() , ch.m2() , m3 , M , L2 , 0 )
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Math::Channel23L::Channel23L
( const Ostap::Math::Channel23L& right ) 
  : ChannelBW ( right      ) 
  , m_channel ( right.m_channel->clone() )
  , m_ps      ( right.m_ps ) 
{}
// ============================================================================
// clone method 
// ============================================================================
Ostap::Math::Channel23L*
Ostap::Math::Channel23L::clone () const 
{ return new Ostap::Math::Channel23L(*this) ; }
// ============================================================================
// unique tag for this lineshape 
// ============================================================================
std::size_t Ostap::Math::Channel23L::tag () const
{ 
  static const std::string s_name = "Channel23L" ;
  return Ostap::Utils::hash_combiner
    ( s_name            ,
      ChannelBW::tag () , 
      m_channel->tag () , 
      m_ps.tag       () ) ; 
}
// ============================================================================
// describe the channel 
// ============================================================================
std::string Ostap::Math::Channel23L::describe() const 
{ return "Channel23L(" + m_channel->describe()  + ")" ; }
// ============================================================================
/*  get the phase space factor  \f$ \varrho \f$
 *  optionally normalized at point \f$ m_n \f$
 */
// ============================================================================
double Ostap::Math::Channel23L::rho_s 
( const double s  , 
  const double mn ) const 
{ 
  const double a = m_ps ( std::sqrt ( s ) ) ;
  return 
    mn <= m_ps.lowEdge  () ? a :
    mn >= m_ps.highEdge () ? a : a / m_ps ( mn )  ;
}

// ============================================================================
// 3-body decays 
// ============================================================================
// constructor from (partial) width and three masses 
// ============================================================================
Ostap::Math::ChannelNR3::ChannelNR3
( const double gamma ,  
  const double m1 , 
  const double m2 , 
  const double m3 ) 
  : ChannelBW ( gamma ) 
  , m_m1 ( s_zero ( m1 ) ? 0.0 : std::abs ( m1 ) ) 
  , m_m2 ( s_zero ( m2 ) ? 0.0 : std::abs ( m2 ) ) 
  , m_m3 ( s_zero ( m3 ) ? 0.0 : std::abs ( m3 ) )
  , m_sthreshold ( 0 ) 
{
  const double m = m_m1 + m_m2 + m_m3 ;
  m_sthreshold = m * m ;
}
// ============================================================================
// clone method
// ============================================================================
Ostap::Math::ChannelNR3*
Ostap::Math::ChannelNR3::clone () const 
{ return new Ostap::Math::ChannelNR3 ( *this ) ; }
// ============================================================================ 
/*  squared  numerator for the amplitude 
 * \f$ N^2(s,m_0) =  m_0 \Gamma_0 \frac{\varrho_3(s)}{\varrho_3(m_0^2)} \f$ 
 */
// ============================================================================ 
double Ostap::Math::ChannelNR3::N2
( const double s  , 
  const double m0 ) const 
{
  return 
    s <= m_sthreshold ? 0.0 : 
    m0 * gamma0 () *
    Ostap::Kinematics::phasespace3 ( std::sqrt ( s ) , m_m1 , m_m2 , m_m3 ) / 
    Ostap::Kinematics::phasespace3 (             m0  , m_m1 , m_m2 , m_m3 ) ;
}
// ============================================================================
/*  term in the denominator for the amplitide
 *  \f$ D(s,m_0) = m_0 \Gamma_0 \frac{\varrho_3(s)}{\varrho_3(m_0^2)} \f$
 */
// ============================================================================
std::complex<double> 
Ostap::Math::ChannelNR3::D    
( const double s  , 
  const double m0 ) const
{
  return 
    s <= m_sthreshold ? 0.0 : 
    m0 * gamma0 () *
    Ostap::Kinematics::phasespace3 ( std::sqrt ( s ) , m_m1 , m_m2 , m_m3 ) / 
    Ostap::Kinematics::phasespace3 (             m0  , m_m1 , m_m2 , m_m3 ) ;
}
// ============================================================================
// unique tag for this lineshape 
// ============================================================================
std::size_t Ostap::Math::ChannelNR3::tag () const
{ 
  static const std::string s_name = "ChannelNR" ;
  return Ostap::Utils::hash_combiner 
    ( s_name            ,
      ChannelBW::tag () , 
      m_m1              , 
      m_m2              , 
      m_m3              ) ; 
}
// ============================================================================
// describe the channel 
// ============================================================================
std::string Ostap::Math::ChannelNR3::describe() const 
{
  return 
    "ChannelNR3(" + std::to_string ( gamma0 () ) + 
    ","           + std::to_string ( m_m1       ) + 
    ","           + std::to_string ( m_m2       ) + 
    ","           + std::to_string ( m_m3       ) + ")" ;
}
// ============================================================================

// ============================================================================
// GammaBW3
// ============================================================================
/* constructor from the Dalitz configuration and 
 *  the squared matrix element
 *  @param dalitz Dalitz configriation
 *  @param me2 squared matrix element \f$ M^2 (s,s_1,s_2) \equiv \frac{1}{2J+1}\sum_i \left| \mathcal{A} \right| \f$
 */
// ============================================================================
Ostap::Math::GammaBW3::GammaBW3 
( const Ostap::Kinematics::Dalitz0&     dalitz , 
  Ostap::Math::GammaBW3::MatrixElement2 me2    , 
  const std::size_t                     tag    , 
  const unsigned short                  n1     , 
  const unsigned short                  n2     )
  : m_me2    ( me2    ) 
  , m_dalitz ( dalitz )
  , m_tag    ( tag    )
  , m_n1     ( n1     )
  , m_n2     ( n2     )
{}
// ============================================================================
// the main method 
// ============================================================================
double Ostap::Math::GammaBW3::GammaBW3::operator() ( const double s ) const 
{
  if ( s <= m_dalitz.s_min () ) { return 0 ; }
  //
  static const double s_CONST = 0.25 * M_PI * M_PI / std::pow ( 2 * M_PI , 5 ) ;
  //
  return s_CONST * Ostap::Math::DalitzIntegrator::integrate_s1s2 
    ( std::cref ( m_me2 ) , s , m_dalitz , m_tag , m_n1 , m_n2 ) / s  ;
}
// ============================================================================
// constructor from (partial) width
// ============================================================================
Ostap::Math::ChannelDalitz::ChannelDalitz 
( const double                                 gamma       , 
  const Ostap::Kinematics::Dalitz0&            dalitz      , 
  Ostap::Math::GammaBW3::MatrixElement2        me2         ,
  const std::size_t                            tag         , 
  const std::string&                           description )
  : ChannelWidth ( gamma                           , 
                   Ostap::Math::GammaBW3 ( dalitz , me2 , tag ) , 
                   dalitz.s_min()                  , 
                   tag                             , 
                   description                     ) 
{}
// ============================================================================
// clone method 
// ============================================================================
Ostap::Math::ChannelDalitz*
Ostap::Math::ChannelDalitz::clone() const 
{ return new Ostap::Math::ChannelDalitz ( *this ) ; }
// ============================================================================


// ============================================================================
// Gounaris-Sakirai
// ============================================================================

// ============================================================================
// constructor with gamma and pion mass
// ============================================================================
Ostap::Math::ChannelGS::ChannelGS 
( const double gamma ,
  const double mpi   ) 
  : ChannelBW ( gamma ) 
  , m_mpi        ( std::abs ( mpi ) ) 
  , m_sthreshold ( 4  * mpi * mpi   ) 
{}
// ============================================================================
// clone method 
// ============================================================================
Ostap::Math::ChannelGS*
Ostap::Math::ChannelGS::clone() const 
{ return new Ostap::Math::ChannelGS (*this ); }
// ============================================================================
// h-function
// ============================================================================
double Ostap::Math::ChannelGS::h ( const double s ) const 
{
  //
  if      ( s_zero ( s ) ) { return s_iPI; }
  else if ( s <= 0       ) 
  {
    const double q2  = m_sthreshold - s ;
    const double sqs = std::sqrt ( -s  ) ;
    const double sq2 = std::sqrt (  q2 ) ;
    return ( sq2 / sqs ) * std::log  ( ( sqs + sq2 ) / ( 2 * m_mpi ) ) * s_iPI ;
  }
  else if ( s <= m_sthreshold ) 
  {
    const double q2 = m_sthreshold - s ;
    const double a  = std::sqrt ( q2 / s ) ;
    return 0.5 * a * ( 1 - 2 * s_iPI * std::atan ( a ) ) ;
  }
  /// physical region 
  const double q2  = s - m_sthreshold ;
  const double sqs = std::sqrt ( s  ) ;
  const double sq2 = std::sqrt ( q2 ) ;
  return ( sq2 / sqs ) * std::log  ( ( sqs + sq2  ) / ( 2 * m_mpi ) ) * s_iPI ;
  //
}
// ============================================================================

// ============================================================================
// h'-function
// ============================================================================
double Ostap::Math::ChannelGS::h_prime ( const double s ) const 
{
  
  if ( s < 0 ) 
  {
    const double q2 = m_sthreshold - s ;
    //
    const double sqs   = std::sqrt ( -s  ) ;
    const double sq2   = std::sqrt (  q2 ) ;
    
    const double t1    = sq2 / sqs ;
    const double et2   = ( sqs + sq2 ) / ( 2 * m_mpi ) ;
    const double t2    = std::log ( et2 ) * s_iPI ;
    //
    const double dt1ds = 0.5 / t1 * ( 1 + q2 / s ) / s ;
    const double dt2ds = - 0.25 * s_iPI / ( et2 * m_mpi ) * ( 1/sqs + 1/sq2 ) ;
    //
    return dt1ds * t2 + t1 * dt2ds ;
  }
  else if ( s < m_sthreshold ) 
  {
    const double q2 = m_sthreshold - s ;
    
    const double a      = q2 / s  ;
    const double sqa    = std::sqrt ( a ) ;
    const double dads   = -1 * ( 1 + q2 / s ) / s ;
    const double dsqads = 0.5 / sqa * dads ;
    
    const double t1     = 0.5 * sqa ;
    const double t2     = ( 1 - 2  * s_iPI * std::atan ( sqa ) ) ;
    //
    const double dt1ds  = 0.5 * dsqads ;
    const double dt2ds  = -2 * s_iPI / ( 1 + a ) * dsqads ;
    //
    return dt1ds * t2 + t1 * dt2ds ;
  }
  /// physical  region 
  const double q2    = s - m_sthreshold ;
  const double sqs   = std::sqrt ( s  ) ;
  const double sq2   = std::sqrt ( q2 ) ;
  //
  const double t1    = sq2 / sqs ;
  const double et2   = ( sqs + sq2  ) / ( 2 * m_mpi ) ;
  const double t2    = std::log  ( et2 ) * s_iPI ;
  //
  const double dt1ds = 0.5 / t1 * ( 1 - q2 / s ) / s ;
  const double dt2ds = 0.25 * s_iPI / ( et2 * m_mpi ) * ( 1/sqs + 1/sq2 ) ;
  //
  return dt1ds * t2 + t1 * dt2ds ;
  //
}
// ============================================================================
/*  term in the denominator for the amplitide
 *  \f$ D(s,m_0) = m_0 \Gamma_0 \gamma(s) + if(s) \f$
 */
// ============================================================================
std::complex<double> Ostap::Math::ChannelGS::D    
( const double s  , 
  const double m0 ) const 
{
  //
  const double q2 = s       - m_sthreshold ;
  const double q0 = m0 * m0 - m_sthreshold ;
  //
  const double grho = s <= m_sthreshold ?  0.0 : 
    gamma0 () * m0 / std::sqrt ( s ) * ( q2 / q0 ) * std::sqrt ( q2 / q0 ) ;
  //
  const double dh = h ( s ) - h ( m0 * m0 ) ;
  //
  const double fs = 
    2 * gamma0 () * m0 * m0 / ( q0 * std::sqrt ( q0 ) ) * 
    ( q2 * dh + q0 *  h_prime ( m0 * m0 ) * ( m0 * m0 - s ) ) ;
  //
  return m0 * grho + s_j * fs ;
}
// ============================================================================
/*   squared  numerator for the amplitude 
 * \f$ N^2(s,m_0) =  m_0 \Gamma_0 \frac{\varrho_3(s)}{\varrho_3(m_0^2)} \f$ 
 */
// ============================================================================
double Ostap::Math::ChannelGS::N2
( const double s  , 
  const double m0 ) const 
{
  if ( s <= m_sthreshold ) { return 0.0 ; }
  //
  const double q2 = s       - m_sthreshold ;
  const double q0 = m0 * m0 - m_sthreshold ;
  //
  return q2 / q0 ;
}
// ======================================================================
/*  get the phase space factor  \f$ \varrho(s) \f$
 *  optionally normalized at the point \f$ m_n \f$ 
 */
// ======================================================================
double Ostap::Math::ChannelGS::rho_s 
( const double s  , 
  const double mn ) const
{ 
  if ( s < m_sthreshold ) { return 0 ; }
  //
  const double mnsq  = mn    * mn    ;
  const double mpisq = m_mpi * m_mpi ;
  return mnsq <= m_sthreshold ? 
    Ostap::Math::PhaseSpace2::phasespace_s ( s    , mpisq , mpisq  ) :
    Ostap::Math::PhaseSpace2::phasespace_s ( s    , mpisq , mpisq  ) /
    Ostap::Math::PhaseSpace2::phasespace_s ( mnsq , mpisq , mpisq  ) ;
}
// ============================================================================
// unique tag for this lineshape 
// ============================================================================
std::size_t Ostap::Math::ChannelGS::tag () const
{ 
  static const std::string s_name = "ChannelGS" ;
  return Ostap::Utils::hash_combiner 
    ( s_name           , 
      ChannelBW::tag() , 
      m_mpi            , 
      m_sthreshold     ) ; 
}
// ============================================================================
// describe the channel 
// ============================================================================
std::string Ostap::Math::ChannelGS::describe() const 
{
  return 
    "ChannelGS(" + std::to_string ( gamma0 () ) + 
    ","          + std::to_string ( m_mpi     ) + ")" ;
}
// ============================================================================


// ============================================================================
/// full constructor 
// ============================================================================
Ostap::Math::ChannelQ::ChannelQ 
( const double gamma ,
  const double m1    , 
  const double m2    ) 
  : ChannelBW ( gamma ) 
  , m_ps2     ( s_zero ( m1 ) || s_zero ( m1 * m1 ) ? 0.0 : std::abs ( m1 ) , 
                s_zero ( m2 ) || s_zero ( m2 * m2 ) ? 0.0 : std::abs ( m2 ) ) 
{}
// ============================================================================
// clone method 
// ============================================================================
Ostap::Math::ChannelQ*
Ostap::Math::ChannelQ::clone() const 
{ return new Ostap::Math::ChannelQ (*this) ; }
// ============================================================================
// unique tag for this lineshape 
// ============================================================================
std::size_t Ostap::Math::ChannelQ::tag () const
{ 
  static const std::string s_name = "ChannelGS" ;
  return Ostap::Utils::hash_combiner 
    ( s_name            ,
      ChannelBW::tag () , 
      m_ps2    . tag () ) ; 
}
// ============================================================================
// describe the channel 
// ============================================================================
std::string Ostap::Math::ChannelQ::describe() const 
{
  return 
    "ChannelQ(" + std::to_string ( gamma0   () ) +
    ","         + std::to_string ( m_ps2.m1 () ) + 
    ","         + std::to_string ( m_ps2.m2 () ) + ")" ;
}
// ============================================================================




// ============================================================================
// clone method
// ============================================================================
Ostap::Math::ChannelGLR*
Ostap::Math::ChannelGLR::clone() const
{ return new ChannelGLR ( *this ) ; }
// ============================================================================
// unique tag for this lineshape 
// ============================================================================
std::size_t Ostap::Math::ChannelGLR::tag () const
{
  static const std::string s_name = "ChannelGLR " ;
  return Ostap::Utils::hash_combiner 
    ( s_name           ,
      ChannelBW::tag() , 
      m_sthreshold     , 
      m_tag            ,
      m_description    ) ;
}
// ============================================================================



// ============================================================================
// clone method
// ============================================================================
Ostap::Math::ChannelNRL*
Ostap::Math::ChannelNRL::clone() const
{ return new Ostap::Math::ChannelNRL ( *this ) ; }

// ============================================================================
// unique tag for this lineshape 
// ============================================================================
std::size_t Ostap::Math::ChannelNRL::tag () const
{
  static const std::string s_name = "ChannelNRL" ;
  return Ostap::Utils::hash_combiner 
    ( s_name            ,
      ChannelBW::tag () ,
      m_sthreshold      , 
      m_fake            , 
      m_tag             ,
      m_description     ) ;
}
// ============================================================================

// ============================================================================
/*  constructor from all masses and angular momenta
 *  @param m0 the mass of K*(1450) 
 *  @param g0 the width of  K*(1430)
 *  @param m1 the mass of the first  particle (kaon)
 *  @param m2 the mass of the second particle (pion)
 *  @param m3 the mass of the third  particle (eta') 
 *  @param a  the LASS parameter a 
 *  @param b  the LASS parameter b 
 *  @param e  the LASS parameter e (elasticity)
 */
// ============================================================================
Ostap::Math::LASS::LASS
( const double m0 ,   // K*(1450) mass
  const double g0 ,   // K*(1430) width
  const double m1 ,   // kaon mass 
  const double m2 ,   // pion mass 
  const double m3 ,   // eta' mass 
  const double a  , 
  const double r  ,
  const double e  )  // elasticity 
  : BW    ( m0 ) 
  , m_a   ( a  ) 
  , m_b   ( a  ) 
  , m_e   ( s_zero ( e ) ? 0.0 : s_equal ( e , 1 ) ? 1 : e ) 
  , m_ps2 ( m1 , m2 )
  , m_m3  ( std::abs ( m3 ) ) 
{
  //
  Ostap::Assert ( 0 <= m_e  && m_e <= 1 , "Invalid elasticity!" , "LASS" ) ;
  const double g1 =       m_e   * std::abs ( g0 ) ;
  const double g2 = ( 1 - m_e ) * std::abs ( g0 ) ;
  //
  add ( Channel ( g1 , m1 , m2 , 0 ) ) ;
  add ( Channel ( g2 , m1 , m3 , 0 ) ) ;
  //
}
// ======================================================================
// clone method
// ======================================================================
Ostap::Math::LASS*
Ostap::Math::LASS::clone() const 
{ return new Ostap::Math::LASS ( *this ) ; }
// ======================================================================
/* LASS amplitude 
 * \f[ \begin{array} {rcl}
 *  \mathcal{A}(m) & = & A_{\mathrm{B}}  +  \\ 
 *                       A_{\mathrm{BW}} {\mathrm{e}}^{i\phi}  \\
 *  A_{\mathrm{B}}  & = & \sin \deeelta \mathrm{e}^{i\delta}   \\ 
 *  \cot \delta     & = & \frac{1}{aq} + \frac{1}{2}bq         \\ 
 &  A_{\mathrm{BW}} & = & \frac{M_R\Gamma_1}{ \left(\M^2_R - M^2 right)
 *                        - iM_R (\Gamma_1 + \Gamma_2 ) }      \\ 
 *  \Gamma_i        & = & q_i \Gamma_{R,i}       \\
 *  \phi            & = & 2\delta  \end{array} \f] 
 */
// ======================================================================
std::complex<double>
Ostap::Math::LASS::amplitude ( const double m ) const
{
  if ( m <= m_ps2.threshold () ) { return 0 ; }
  //
  const std::complex<double> bw = 
    BW::amplitude ( m ) * channel()->N2 ( m * m , m0 () ) ;
  //
  const double q         = m_ps2.q ( m ) ;
  const double cot_delta = 1 / ( q * m_a ) + 0.5  * m_b * q ;
  const double delta     = std::atan ( 1.0 / cot_delta ) ;
  const double sin_delta = std::sin  ( delta   ) ;
  const double cos_delta = cot_delta * sin_delta ;
  const double phi       = 2 * delta             ;
  const double sin_phi   = 2 * sin_delta * cos_delta     ;
  const double cos_phi   = 2 * cos_delta * cos_delta - 1 ;
  //
  const std::complex<double> bg = sin_delta  * 
    std::complex<double>  ( cos_delta , sin_delta ) ;
  //
  return bg + bw * std::complex<double>  ( cos_phi , sin_phi ) ;
}
// ============================================================================
// evaluate LASS function 
// ============================================================================
double Ostap::Math::LASS::operator () ( const double m ) const 
{
  return m <= m_ps2.threshold () ? 0.0 :
    2 * m * m_ps2 ( m ) * std::norm ( amplitude ( m ) ) / M_PI ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::LASS::setA ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_a ) ) { return false ; }
  //
  m_a = v ;
  //
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::LASS::setB ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_b ) ) { return false ; }
  //
  m_b = v ;
  //
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::LASS::setE ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_e ) ) { return false ; }
  //
  m_e = v ;
  const double g  = gamma () ;
  //
  m_channels [ 0 ] -> setGamma0 (       m_e   * g ) ;
  m_channels [ 1 ] -> setGamma0 ( ( 1 - m_e ) * g ) ;
  //
  return true ;
}
// ============================================================================
// unique label/tag 
// ============================================================================
std::size_t Ostap::Math::LASS::tag () const 
{
  static const std::string s_name = "LASS" ;
  return Ostap::Utils::hash_combiner 
    ( s_name , BW::tag () , m_a , m_b , m_e ) ; 
}
// ============================================================================
/*  constructor from Breit-Wigner, Phase-space and flags 
 *  @param bw Breit-Wigner shape 
 *  @param ps phase-space function 
 *  @param use_rho  use rho-fuction from Breit-Wigner 
 *  @param use_N2   use N2-function from Breit-Wigner 
 */
// ============================================================================
Ostap::Math::BWPS::BWPS 
( const Ostap::Math::BW&            bw , 
  const Ostap::Math::PhaseSpacePol& ps , 
  const bool use_rho                   , 
  const bool use_N2                    ) 
  : m_rho ( use_rho    ) 
  , m_N2  ( use_N2     )
  , m_bw  ( bw.clone() ) 
  , m_ps  ( ps         )
  , m_workspace () 
{}
// ============================================================================
/*  constructor from Breit-Wigner, Phase-space and flags 
 *  @param bw Breit-Wigner shape 
 *  @param ps phase-space function 
 *  @param use_rho  use rho-fuction from Breit-Wigner 
 *  @param use_N2   use N2-function from Breit-Wigner 
 */
// ============================================================================
Ostap::Math::BWPS::BWPS 
( const Ostap::Math::BW&            bw , 
  const Ostap::Math::PhaseSpaceNL&  ps , 
  const bool use_rho                   , 
  const bool use_N2                    ) 
  : m_rho ( use_rho    ) 
  , m_N2  ( use_N2     )
  , m_bw  ( bw.clone() ) 
  , m_ps  ( ps , 0     )
  , m_workspace () 
{
  Ostap::Assert ( ( ps.xmin      () <= bw.threshold () ) && 
                  ( bw.threshold () <  ps.xmax      () )  , 
                  "Invalid xmin/xmax/threhsold setting", 
                  "Ostap::Math::BWPS"  ) ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Math::BWPS::BWPS 
( const Ostap::Math::BWPS& right )
  : m_rho ( right.m_rho )   
  , m_N2  ( right.m_N2  )   
  , m_bw  ( right.m_bw->clone () )   
  , m_ps  ( right.m_ps  )
  , m_workspace()
{}
// ============================================================================
// evaluate the function 
// ============================================================================
double Ostap::Math::BWPS::evaluate ( const double x ) const 
{
  //
  if ( x <= m_ps.xmin() || x >= m_ps.xmax() || x <= m_bw->threshold() ) { return 0 ; }
  //
  double result = m_ps ( x ) ;
  // 
  if ( result <= 0 ) { return 0 ; }
  //
  if      (  m_rho &&  m_N2 ) { result *= (*m_bw)( x ) ; }
  else if (  m_rho && !m_N2 ) 
  { result *= 2 * x * m_bw->rho_s ( x * x ) * std::norm ( m_bw->amplitude ( x ) ) ; }
  else if ( !m_rho &&  m_N2 ) 
  { result *= 2 * x * m_bw->N2    ( x * x ) * std::norm ( m_bw->amplitude ( x ) ) ; }
  else  
  { result *= 2 * x                         * std::norm ( m_bw->amplitude ( x ) ) ; }
  //
  return result ;
}
// ============================================================================
// evaluate the integral 
// ============================================================================
double Ostap::Math::BWPS::integral() const 
{
  const double x1 = std::max ( m_ps.xmin() , m_bw->threshold() ) ;
  const double x2 = m_ps.xmax() ;
  if ( x2 <= x1 ) { return 0 ; }
  return integral ( x1 , x2 ) ;
}
// =============================================================================
// evaluate the integral 
// ============================================================================
double Ostap::Math::BWPS::integral
( const double x_low  , 
  const double x_high ) const 
{
  //
  if       ( s_equal ( x_high , x_low ) ) { return 0 ; } 
  else if  (           x_high < x_low   ) { return -1 * integral ( x_high , x_low ) ; }
  //
  const double x1 = std::max ( m_ps.xmin() , m_bw->threshold() ) ;
  const double x2 = m_ps.xmax() ;
  //
  if ( x2 <= x1 ) { return 0 ; }
  //
  const double xlow  =  std::max ( x1 , x_low  ) ;
  const double xhigh =  std::min ( x2 , x_high ) ;
  //
  if ( xhigh <= xlow ) { return 0 ; }
  //
  if ( ( xhigh - xlow ) > 0.2 * ( x2 - x1 ) ) 
  {
    const double xmid = 0.5 * ( xlow + xhigh ) ;
    return integral ( xlow , xmid ) + integral ( xmid , xhigh ) ; 
  }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<BWPS> s_integrator {} ;
  static char s_message [] = "Integral(BWPS)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag  () , 
      &F      , 
      xlow    , xhigh     ,          // low & high edges
      workspace ( m_workspace ) ,    // workspace
      s_APRECISION         ,          // absolute precision
      s_RPRECISION         ,          // relative precision
      m_workspace.size()  ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
  //
}
// ============================================================================
// unique label/tag 
// ============================================================================
std::size_t Ostap::Math::BWPS::tag () const 
{
  static const std::string s_name = "BWPS" ;
  return Ostap::Utils::hash_combiner  
    ( s_name , m_ps .tag () , m_bw->tag () , m_rho , m_N2 ) ; 
}

// ============================================================================
/*  constructor from Breit-Wigner, Phase-space and flags 
 *  @param bw Breit-Wigner shape 
 *  @param M  mass of the "mother" particle 
 *  @param m1 mass of the 1st      particle 
 *  @param m2 mass of the 2nd      particle 
 *  @param m0 mass of the 3rd      particle 
 *  @param L  the orbital momentum between  
 *  system of 1st and 2nd particles and the 3rd particle
 */
// ============================================================================
Ostap::Math::BW3L::BW3L
( const Ostap::Math::BW& bw , 
  const double           M  ,   
  const double           m1 ,         
  const double           m2 ,       
  const double           m3 ,       
  const unsigned short   L  ) 
  : m_bw ( bw.clone() ) 
  , m_M  ( M  ) 
  , m_m1 ( m1 ) 
  , m_m2 ( m2 ) 
  , m_m3 ( m3 ) 
  , m_L  ( L  ) 
  , m_p0 ( 1  )
  , m_workspace () 
{
  Ostap::Assert ( 0 <  M           ,
                  "Invalid mother mass"           , "Ostap:Math::BW3L" ) ;
  Ostap::Assert ( 0 <= m1          , 
                  "Invalid 1st    mass"           , "Ostap:Math::BW3L" ) ;
  Ostap::Assert ( 0 <= m2          ,
                  "Invalid 2nd    mass"           , "Ostap:Math::BW3L" ) ;
  Ostap::Assert ( 0 <= m3          , 
                  "Invalid 3rd    mass"           , "Ostap:Math::BW3L" ) ;
  Ostap::Assert ( m1 + m2 + m3 < M , 
                  "Invalid combination of masses" , "Ostap:Math::BW3L" ) ;
  Ostap::Assert ( m_bw -> threshold () < M - m3 , 
                  "Inconsistent BW-threshold"     , "Ostap:Math::BW3L" ) ;
  // cache:
  const double xmid = 0.5 * ( xmin () + xmax () ) ;
  m_p0 = Ostap::Math::PhaseSpace2::q ( m_M , xmid , m_m3  ) ;
  //
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Math::BW3L::BW3L 
( const Ostap::Math::BW3L& right )
  : m_bw  ( right.m_bw->clone () )   
  , m_M   ( right.m_M  ) 
  , m_m1  ( right.m_m1 ) 
  , m_m2  ( right.m_m2 ) 
  , m_m3  ( right.m_m3 ) 
  , m_L   ( right.m_L  ) 
  , m_p0  ( right.m_p0 )
  , m_workspace()
{}
// ============================================================================
// evaluate the function 
// ============================================================================
double Ostap::Math::BW3L::evaluate   ( const double x ) const 
{
  if ( x <= xmin() || x >= xmax() ) { return 0 ; }           // RETURN
  //
  const double p  = Ostap::Math::PhaseSpace2::q ( m_M , x , m_m3  ) ;
  if ( p <= 0                     ) { return 0 ; }          //  RETURN
  //
  return (*m_bw)(x) * std::pow ( p / m_p0 , 2 * m_L + 1 ) ;
}
// ============================================================================
// evaluate the integral 
// ============================================================================
double Ostap::Math::BW3L::integral() const 
{ return integral ( xmin () , xmax () ) ; }
// =============================================================================
// evaluate the integral 
// ============================================================================
double Ostap::Math::BW3L::integral
( const double x_low  , 
  const double x_high ) const 
{
  if       ( s_equal ( x_high , x_low ) ) { return 0 ; } 
  else if  (           x_high < x_low   ) { return -1 * integral ( x_high , x_low ) ; }
  //
  const double x_mn = xmin () ;
  const double x_mx = xmax () ;
  //
  if  ( x_mx < x_low || x_mn > x_high ) { return 0 ; }
  //
  const double xlow  = std::max ( x_mn , x_low  ) ;
  const double xhigh = std::min ( x_mx , x_high ) ;
  //
  if ( xhigh <= xlow ) { return 0 ; }
  //
  if ( ( xhigh - xlow ) > 0.2 * ( x_mx - x_mn ) ) 
  {
    const double xmid = 0.5 * ( xlow + xhigh ) ;
    return integral ( xlow , xmid ) + integral ( xmid , xhigh ) ; 
  }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<BW3L> s_integrator {} ;
  static char s_message [] = "Integral(BW3P)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag  () , 
      &F      , 
      xlow    , xhigh     ,          // low & high edges
      workspace ( m_workspace ) ,    // workspace
      s_APRECISION         ,          // absolute precision
      s_RPRECISION         ,          // relative precision
      m_workspace.size()  ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
  //
}
// ============================================================================
// unique label/tag 
// ============================================================================
std::size_t Ostap::Math::BW3L::tag () const 
{ 
  static const std::string s_name = "BW3L" ;
  return Ostap::Utils::hash_combiner 
    ( s_name , m_bw -> tag () , m_M , m_m1 , m_m2 , m_m3 , m_L ) ; 
}
// ============================================================================


// ============================================================================
// constructor from breit-wigner
// ============================================================================
Ostap::Math::A2::A2
( const Ostap::Math::BW& bw    , 
  const double           scale ) 
  : m_bw    ( bw.clone() ) 
  , m_scale ( scale      )  
{}
// ============================================================================
Ostap::Math::A2::A2
( const Ostap::Math::A2& bw ) 
  : m_bw   ( bw.m_bw ? bw.m_bw->clone() : nullptr ) 
  , m_scale ( bw.m_scale      )  
{}
// ============================================================================
double Ostap::Math::A2::operator() ( const double s ) const 
{
  if ( !m_bw ) { return 0 ; }
  if ( s < m_bw -> s_threshold () ) { return 0 ; }
  const double m = std::sqrt ( s ) ;
  return m_scale * std::norm ( m_bw->amplitude ( m ) ) ;
}
// ============================================================================
// unique tag 
// ============================================================================
std::size_t Ostap::Math::A2::tag () const 
{
  static const std::string s_name = "A2" ;
  return Ostap::Utils::hash_combiner 
    ( s_name , m_bw ? m_bw -> tag () : 0 , m_scale ) ; 
}
// ============================================================================

// ============================================================================
//                                                                      The END 
// ============================================================================
