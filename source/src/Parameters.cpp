// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================ \
#include <algorithm>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Hash.h"
#include "Ostap/Names.h"
#include "Ostap/Math.h"
#include "Ostap/MoreMath.h"
#include "Ostap/Parameters.h"
// ============================================================================
// Local
// ============================================================================
#include "local_hash.h"
#include "local_math.h"
#include "status_codes.h"
// ============================================================================
/** @file 
 *  Implementation file for class Ostap::Math::Parameters
 *  @see Ostap::Math::Parameters
 *  @author Vanya BELYAEV
 */
// ============================================================================
/*  full constructor
 *  @param value parameter value  
 *  @param name  parameter name  
 *  @param the_class name of the (owner) class 
 */
// ============================================================================
Ostap::Math::Value::Value
( const double       value      ,
  const std::string& value_name ,
  const std::string& class_name )
  : m_value ( value      )
  , m_name  ()
{
  setFullName ( value_name , class_name ) ;
}
// ============================================================================
/** full constructor
 *  @param value parameter value  
 *  @param name  parameter name  
 *  @param the_class the (owner/holder) object
 */
// ============================================================================
Ostap::Math::Value::Value
( const double          value     , 
  const std::string&    name      ,
  const std::type_info& the_class )
  : Value ( value , name , Ostap::class_name ( the_class ) )
{}
// ============================================================================
// set new value for parameter 
// ============================================================================
bool Ostap::Math::Value::setValue ( const double value )
{
  if ( s_equal ( value , m_value ) ) { return false ; }
  m_value = value ;
  return true ;
}
// ============================================================================
/** set the full name of parameter
 *  @param the_class the name of ower/holder class
 *  @parm  the_name  the parameter  name
 */
// ============================================================================
const std::string& Ostap::Math::Value::setFullName
( const std::string& the_class ,
  const std::string& the_name  )
{
  //
  const std::string c1 { Ostap::strip ( the_class ) } ;
  const std::string c2 { Ostap::strip ( the_name  ) } ;
  //
  if       ( !c1.empty() && !c2.empty() ) { m_name = c1 + "::" + c2 ; }
  else if  ( !c1.empty()                ) { m_name = c1 ; }
  else if  ( !c2.empty()                ) { m_name = c2 ; }
  else                                    { m_name = "" ; }
  //
  return m_name  ;
}
// ============================================================================
/*  set the full name of parameter
 *  @param the_class the typeinfo of ower/holder class
 *  @parm  the_name  the parameter  name
 */
// ============================================================================
const std::string&
Ostap::Math::Value::setFullName
( const std::type_info& the_class ,
  const std::string&    the_name  )
{ return setFullName ( Ostap::class_name ( the_class ) , the_name ) ; }
// ============================================================================
// Unique value 
// ============================================================================
std::size_t Ostap::Math::Value::tag () const 
{
  static const std::string s_name { "Value" } ;
  return Ostap::Utils::hash_combiner ( s_name , m_value , m_name )  ;
}
// ============================================================================

// ============================================================================
/*  full constructor
 *  @param value logarithm of parameter
 *  @param name  parameter name  
 *  @param the_class name of the (owner/holder) class 
 */
// ============================================================================
Ostap::Math::LogValue::LogValue
( const double       log_value ,
  const std::string& name      ,
  const std::string& the_class )
  : m_logValue ( 0 )
  , m_value    ( 0 , name , the_class )
{
  static const std::string s_invalid_log { "Invalid log-value!" } ;  
  Ostap::Assert ( s_EXP_UNDERFLOW < log_value && log_value < s_EXP_OVERFLOW ,		  
		  s_invalid_log                              ,
		  m_value.name ()                            , 
		  INVALID_LOGPARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_logValue = log_value ;
  m_value.setValue ( std::exp ( m_logValue ) ) ;
  //
}
// ============================================================================
/*  full constructor
 *  @param value parameter value  
 *  @param name  parameter name  
 *  @param the_class the (owner/holder) object
 */
// ============================================================================
Ostap::Math::LogValue::LogValue
( const double          value     , 
  const std::string&    name      ,
  const std::type_info& the_class )
  : LogValue ( value , name , Ostap::class_name ( the_class ) )
{}
// =====================================================================
// set new value for parameter 
// ============================================================================
bool Ostap::Math::LogValue::setLogValue ( const double log_value )
{
  if ( s_equal ( log_value , m_logValue ) ) { return false ; }
  static const std::string s_invalid_log { "Invalid log-value!" } ;  
  Ostap::Assert ( s_EXP_UNDERFLOW < log_value && log_value < s_EXP_OVERFLOW ,		  
		  s_invalid_log                               ,
		  m_value.name ()                             , 
		  INVALID_LOGPARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_logValue = log_value ;  
  return m_value.setValue ( std::exp ( m_logValue ) ) ;
}
// =====================================================================
// set new value for     parameter 
// =====================================================================
bool Ostap::Math::LogValue::setValue ( const double value  )
{
  if ( s_equal ( m_value.value() , value ) ) { return false ; }
  //
  static const std::string s_invalid_log { "Invalid log/exp-value!" } ;  
  Ostap::Assert ( s_EXP_UNDERFLOW_EXP < value && value < s_EXP_OVERFLOW_EXP ,		  
		  s_invalid_log                               ,
		  m_value.name ()                             , 
		  INVALID_LOGPARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_logValue = std::log   ( value ) ;  
  return m_value.setValue ( value ) ;
}
// ============================================================================
// Unique value 
// ============================================================================
std::size_t Ostap::Math::LogValue::tag () const 
{
  static const std::string s_name { "LogValue" } ;
  return Ostap::Utils::hash_combiner ( s_name , m_value.tag () )  ;
}
// ============================================================================



// ============================================================================
/*  full constructor
 *  @param value parameter value (external)
 *  @param avalue A0value 
 *  @param bvalue B-value 
 *  @param name  parameter name  
 *  @param the_class name of the (owner/holder) class 
 */
// ============================================================================
Ostap::Math::InRange::InRange
( const double       extvalue  , 
  const double       avalue    ,
  const double       bvalue    , 
  const std::string& name      ,
  const std::string& the_class )
  : m_A         ( avalue   ) 
  , m_B         ( bvalue   ) 
  , m_external  ( extvalue )
  , m_value     ( 0 , name , the_class ) 
{
  static const std::string s_invalid_range { "Invalid minmax-range!" } ;  
  Ostap::Assert ( !s_equal ( m_A , m_B ) ,
		  s_invalid_range ,
		  m_value.name () ,		  
		  INVALID_RANGE   , __FILE__ , __LINE__  ) ;
  //
  m_value.setValue ( z ( extvalue ) ) ;
}
// ============================================================================
/*  full constructor
 *  @param value parameter value (external)
 *  @param avalue A0value 
 *  @param bvalue B-value 
 *  @param name  parameter name  
 *  @param the_class name of the (owner/holder) class 
 */
// ============================================================================
Ostap::Math::InRange::InRange
( const double          extvalue  , 
  const double          avalue    ,
  const double          bvalue    , 
  const std::string&    name      ,
  const std::type_info& the_class )
  : InRange ( extvalue ,
	      avalue   ,
	      bvalue   ,
	      name     ,
	      Ostap::class_name ( the_class ) )
{}	      
// ============================================================================
// external -> internal 
// ============================================================================
double Ostap::Math::InRange::z
( const double x ) const
{
  const double s2 = std::sin ( s_pi2 * x ) ; 
  return m_A + ( m_B - m_A ) * s2 * s2 ;
}  
// ============================================================================
// internal -> external 
// ============================================================================
double Ostap::Math::InRange::x
( const double z ) const
{
  const double zs2 = ( z - m_A ) / ( m_B - m_A ) ;
  if      ( s_zero  ( zs2     ) ) { return m_A ; }
  else if ( s_equal ( zs2 , 1 ) ) { return m_B ; }
  //
  static const std::string s_invalid_z { "Invalid internal range!" } ;  
  Ostap::Assert ( 0 <= zs2 && zs2 <= 1 ,
		  s_invalid_z        ,
		  m_value.name ()    ,
		  INVALID_PARAMETER  , __FILE__ , __LINE__  ) ;
  //
  const double zs = std::sqrt ( zs2 ) ;
  return std::asin ( zs ) ;
}
// =====================================================================
// set new value for parameter 
// ============================================================================
bool Ostap::Math::InRange::setExternal ( const double value )
{
  if ( s_equal ( value , m_external ) ) { return false ; }
  m_external = value ;  
  return m_value.setValue ( z ( m_external ) ) ;
}
// =====================================================================
// set new value for     parameter 
// =====================================================================
bool Ostap::Math::InRange::setValue ( const double value  )
{
  if ( s_equal ( m_value.value() , value ) ) { return false ; }
  m_external  = x ( value ) ;  
  return m_value.setValue ( value ) ;
}
// ============================================================================
// Unique value 
// ============================================================================
std::size_t Ostap::Math::InRange::tag () const 
{
  static const std::string s_name { "InRange" } ;
  return Ostap::Utils::hash_combiner ( s_name , m_value.tag () )  ;
}
// ============================================================================

// ============================================================================
/*  full constructor
 *  @param value parameter value  
 *  @param name  parameter name  
 *  @param the_class name of the (owner/holder) class 
 */
// ===========================================================================
Ostap::Math::Scale::Scale
( const double       value     ,
  const std::string& name      ,
  const std::string& the_class ,
  const bool         positive  )
  : m_scale    ( positive ? std::abs ( value ) : value , name , the_class )
  , m_positive ( positive ) 
{
  static const std::string s_invalid_scale { "Invalid scale-value!" } ;  
  Ostap::Assert ( !s_zero ( scale () ) ,		  
		  s_invalid_scale      ,                               
		  m_scale.name ()      , 
		  INVALID_SCALEPARAMETER , __FILE__ , __LINE__ ) ;
  //
}
// ============================================================================
/*  full constructor
 *  @param value parameter value  
 *  @param name  parameter name  
 *  @param the_class name of the (owner/holder) class 
 */
// ============================================================================
Ostap::Math::Scale::Scale
( const double          value     ,
  const std::string&    name      ,
  const std::type_info& the_class , 
  const bool            positive  )
  : Scale ( value , name , Ostap::class_name ( the_class ) , positive )
{}
// ============================================================================
// set new value for scale parameter 
// ============================================================================
bool Ostap::Math::Scale::setValue ( const double value )
{
  const double new_value = m_positive ? std::abs ( value ) : value ;
  if ( s_equal ( new_value , m_scale.value () ) ) { return false ; }
  //
  static const std::string s_invalid_scale { "Invalid scale-value!" } ;
  Ostap::Assert ( !s_zero ( new_value) ,		  
		  s_invalid_scale      ,                               
		  m_scale.name ()      , 
		  INVALID_SCALEPARAMETER , __FILE__ , __LINE__ ) ;
  //
  return m_scale.setValue ( new_value ) ; 
}
// ============================================================================
// Unique value 
// ============================================================================
std::size_t Ostap::Math::Scale::tag () const 
{
  static const std::string s_name { "Scale" } ;
  return Ostap::Utils::hash_combiner ( s_name , m_scale.tag () , m_positive )  ;
}
// ============================================================================

// ============================================================================
/*  @param scale the value of scale parameter 
 *  @param shift the value of scale parameter 
 *  @param scale_name the name of scale parameter
 *  @param shift_name the name of shift parameter
 *  @param the_class name of the (owner/holder) class 
 */
// ============================================================================
Ostap::Math::ShiftAndScale::ShiftAndScale
( const double       scale      ,
  const double       shift      ,
  const std::string& scale_name ,
  const std::string& shift_name ,
  const std::string& the_class  , 
  const bool         positive   )
  : m_scale ( scale , scale_name , the_class , positive )
  , m_shift ( shift , shift_name , the_class )
{}
// ============================================================================
/*  @param scale the value of scale parameter 
 *  @param shift the value of scale parameter 
 *  @param scale_name the name of scale parameter
 *  @param shift_name the name of shift parameter
 *  @param the_class name of the (owner/holder) class 
 */
// ============================================================================
Ostap::Math::ShiftAndScale::ShiftAndScale
( const double          scale      ,
  const double          shift      ,
  const std::string&    scale_name ,
  const std::string&    shift_name ,
  const std::type_info& the_class  , 
  const bool            positive   )
  : ShiftAndScale ( scale      ,
		    shift      ,
		    scale_name ,
		    shift_name , Ostap::class_name ( the_class ) , positive )
{}
// ============================================================================
/*  set the full name of parameter
 *  @param the_class the name of ower/holder class
 *  @parm  scale_name  the parameter  name
 *  @parm  shift_name  the parameter  name
 */
// ============================================================================
void Ostap::Math::ShiftAndScale::setFullName
( const std::string&    the_class   ,
  const std::string&    scale_name  , 
  const std::string&    shift_name  )
{
  m_scale.setFullName ( the_class , scale_name ) ;
  m_shift.setFullName ( the_class , scale_name ) ;
}
// ============================================================================
/* set the full name of parameter
 *  @param the_class the type-info of owner/holder class
 *  @parm  the_name  the parameter  name
 */
// ============================================================================
void Ostap::Math::ShiftAndScale::setFullName
( const std::type_info& the_class  ,
  const std::string&    scale_name , 
  const std::string&    shift_name )
{
  setFullName ( Ostap::class_name ( the_class ) ,
		scale_name ,
		shift_name ) ; 
}
// ============================================================================
// Unique value 
// ============================================================================
std::size_t Ostap::Math::ShiftAndScale::tag () const 
{
  static const std::string s_name { "ShiftAndScale" } ;
  return Ostap::Utils::hash_combiner ( s_name         ,
				       m_shift.tag () , 
				       m_scale.tag () )  ;
}
// ============================================================================


// ============================================================================
/*  @param loga logarithm of a-parameter 
 *  @param logb logarithm of b-parameter 
 *  @param aname the name of a-parameter
 *  @param bname the name of b-parameter
 *  @param the_class name of the (owner/holder) class 
 */
// ===========================================================================
Ostap::Math::AB::AB
( const double       loga      ,
  const double       logb      , 
  const std::string& aname     , 
  const std::string& bname     , 
  const std::string& the_class )
  : m_a ( loga , aname, the_class ) 
  , m_b ( logb , bname, the_class ) 
{}
// ============================================================================
/*  @param loga logarithm of a-parameter 
 *  @param logb logarithm of b-parameter 
 *  @param aname the name of b-parameter
 *  @param bname the name of b-parameter
 *  @param the_class name of the (owner/holder) class 
 */
// ============================================================================
Ostap::Math::AB::AB
( const double          loga      ,
  const double          logb      , 
  const std::string&    aname     , 
  const std::string&    bname     ,
  const std::type_info& the_class )
  : AB ( loga , logb , aname , bname , Ostap::class_name ( the_class ) )
{}
// ============================================================================
/* set the full name of parameter
 *  @param the_class the name of ower/holder class
 *  @parm  pname  the name of a-parameter
 *  @parm  qname  the name of b-parameter
 */
// ============================================================================
void Ostap::Math::AB::setFullName
( const std::string& the_class ,
  const std::string& aname    , 
  const std::string& bname    )
{
  m_a.setFullName ( the_class , aname ) ;
  m_b.setFullName ( the_class , bname ) ; 
}
// ======================================================================			
/*  set the full name of parameter
 *  @param the_class the type-info of owner/holder class
 *  @parm  pname  the name of a-parameter
 *  @parm  qname  the name of b-parameter
 */
// ======================================================================			
void Ostap::Math::AB::setFullName
( const std::type_info& the_class  ,
  const std::string& aname    , 
  const std::string& bname    )
{ return setFullName ( Ostap::class_name ( the_class ) , aname , bname ) ; }
// ============================================================================
// Unique value 
// ============================================================================
std::size_t Ostap::Math::AB::tag () const 
{
  static const std::string s_name { "AB" } ;
  return Ostap::Utils::hash_combiner ( s_name     ,
				       m_a.tag () , 
				       m_b.tag () )  ;
}
// ============================================================================
// helper expression \f$ \log \Beta (a , b ) \f$
// ============================================================================
double Ostap::Math::AB::log_Beta_ab() const
{  return Ostap::Math::lnbeta ( m_a.value () , m_b.value () ) ; }
// ============================================================================
// helper expression \f$ \frac{1}{\Beta (a , b ) } \f$
// ============================================================================
double Ostap::Math::AB::inv_Beta_ab() const
{  return Ostap::Math::ibeta  ( m_a.value () , m_b.value () ) ; }
// ============================================================================

// ============================================================================
/*  @param logp  logarithm of p-parameter 
 *  @param logq  logarithm of q-parameter 
 *  @param pname the name of p-parameter
 *  @param qname the name of q-parameter
 *  @param the_class name of the (owner/holder) class 
 */
// ============================================================================
Ostap::Math::PQ::PQ
( const double       logp      ,
  const double       logq      ,
  const std::string& pname     ,
  const std::string& qname     ,
  const std::string& the_class )
  : m_p ( logp , pname , the_class )
  , m_q ( logq , qname , the_class )
{
  m_log_Beta_pq = Ostap::Math::lnbeta ( m_p.value () , m_q.value () ) ;  
  m_inv_Beta_pq = Ostap::Math::ibeta  ( m_p.value () , m_q.value () ) ;  
}
// ============================================================================
/*  @param logp  logarithm of p-parameter 
 *  @param logq  logarithm of q-parameter 
 *  @param pname the name of p-parameter
 *  @param qname the name of q-parameter
 *  @param the_class name of the (owner/holder) class 
 */
// ============================================================================
Ostap::Math::PQ::PQ
( const double          logp      ,
  const double          logq      ,
  const std::string&    pname     ,
  const std::string&    qname     ,
  const std::type_info& the_class )
  : PQ ( logp , logq , pname , qname , Ostap::class_name ( the_class ) )
{}
// ============================================================================
// update log-P
// ============================================================================
bool Ostap::Math::PQ::setLogP ( const double value  )
{
  if ( !m_p.setLogValue ( value ) ) { return false ; } 
  m_log_Beta_pq = Ostap::Math::lnbeta ( m_p.value () , m_q.value () ) ;
  m_inv_Beta_pq = Ostap::Math::ibeta  ( m_p.value () , m_q.value () ) ;  
  return true ;
} 
// ============================================================================
// update log-Q
// ============================================================================
bool Ostap::Math::PQ::setLogQ ( const double value  )
{
  if ( !m_p.setLogValue ( value ) ) { return false ; } 
  m_log_Beta_pq = Ostap::Math::lnbeta ( m_p.value () , m_q.value () ) ;
  m_inv_Beta_pq = Ostap::Math::ibeta  ( m_p.value () , m_q.value () ) ;  
  return true ;
} 
// ============================================================================
// update P
// ============================================================================
bool Ostap::Math::PQ::setP ( const double value  )
{
  if ( !m_p.setValue ( value ) ) { return false ; } 
  m_log_Beta_pq = Ostap::Math::lnbeta ( m_p.value () , m_q.value () ) ;
  m_inv_Beta_pq = Ostap::Math::ibeta  ( m_p.value () , m_q.value () ) ;  
  return true ;
} 
// ============================================================================
// update Q
// ============================================================================
bool Ostap::Math::PQ::setQ ( const double value  )
{
  if ( !m_q.setValue ( value ) ) { return false ; } 
  m_log_Beta_pq = Ostap::Math::lnbeta ( m_p.value () , m_q.value () ) ;
  m_inv_Beta_pq = Ostap::Math::ibeta  ( m_p.value () , m_q.value () ) ;  
  return true ;
}
// ============================================================================
/*  set the full name of parameter
 *  @param the_class the name of ower/holder class
 *  @parm  pname  the name of p-parameter
 *  @parm  qname  the name of q-parameter
 */
// ============================================================================
void Ostap::Math::PQ::setFullName
( const std::string& the_class ,
  const std::string& pname     , 
  const std::string& qname     )
{
  m_p.setFullName ( the_class , pname ) ;
  m_q.setFullName ( the_class , pname ) ;
}
// ============================================================================
/*  set the full name of parameter
 *  @param the_class the name of ower/holder class
 *  @parm  pname  the name of p-parameter
 *  @parm  qname  the name of q-parameter
 */
// ============================================================================
void Ostap::Math::PQ::setFullName
( const std::type_info& the_class  ,
  const std::string& pname     , 
  const std::string& qname     )
{ return setFullName ( Ostap::class_name ( the_class ) , pname , qname ) ;}
// ============================================================================
// Unique value 
// ============================================================================
std::size_t Ostap::Math::PQ::tag () const 
{
  static const std::string s_name { "PQ" } ;
  return Ostap::Utils::hash_combiner ( s_name     ,
				       m_p.tag () , 
				       m_q.tag () )  ;
}
// ============================================================================











// ============================================================================
// PARAMETERS
// ============================================================================
// constructor from number of parameters 
// ============================================================================
Ostap::Math::Parameters::Parameters 
( const std::size_t np ) 
  : m_pars ( np , 0.0 ) 
{}
// ============================================================================
// constructor from  the list of parameters 
// ============================================================================
Ostap::Math::Parameters::Parameters 
( const std::vector<double>&  pars   ) 
  : m_pars ( pars ) 
{}
// ============================================================================
// constructor from the list of parameters 
// ============================================================================
Ostap::Math::Parameters::Parameters 
(       std::vector<double>&& pars   ) 
  : m_pars ( std::forward<std::vector<double> > ( pars ) ) 
{}
// ============================================================================
// all zero ?
// ============================================================================
bool Ostap::Math::Parameters::zero  () const { return s_vzero ( m_pars ) ; }
// ============================================================================
// set k-parameter
// ============================================================================
bool Ostap::Math::Parameters::_setPar 
( const std::size_t k     , 
  const double      value ,
  const bool        force ) 
{
  if ( m_pars.size() <= k                         ) { return false ; }
  if ( s_equal ( m_pars [ k ] , value ) && !force ) { return false ; }
  m_pars [ k ] = value ;
  return true ;
}
// ============================================================================
// rest all parameters to zero 
// ============================================================================
void Ostap::Math::Parameters::reset () 
{ std::fill ( m_pars.begin() , m_pars.end() ,  0.0 ) ; }
// ============================================================================
// swap two objects 
// ============================================================================
void Ostap::Math::Parameters::swap (  Ostap::Math::Parameters& right ) 
{ std::swap ( m_pars ,  right.m_pars ); }
// ============================================================================
/*  filter out very small terms/parameters 
 *  the term is considered to be very small if
 *   - it is numerically zero: \f$ c_k \approx 0 \f% 
 *   - or if epsilon  > 0: \f$ \left| c_k \right| \le \epsilon \f$ 
 *   - or if scale   != 0: \f$ \left| s \right| + \left| c_k \right| \approx \left| s \right| \f$ 
 *  @param  epsilon  parameter to define "smalness" of terms
 *  @param  scale    parameter to define "smalness" of terms
 *  @return number of nullified terms
 */
// ============================================================================
std::size_t
Ostap::Math::Parameters::remove_noise
( const double epsilon , 
  const double scale   )
{
  std::size_t  num    = 0                  ;
  const bool   eps    = 0 <  epsilon       ;
  const bool   sca    = scale              ;
  const double ascale = std::abs ( scale ) ;
  //
  const std::size_t NN { m_pars.size() } ;
  for ( std::size_t k = 0 ; k < NN ; ++k )
    {
      const double absp = std::abs ( m_pars [ k ] ) ;
      if      ( s_zero ( absp )                           ) { m_pars [ k ] = 0 ; ++num ; }
      else if ( eps && absp <= epsilon                    ) { m_pars [ k ] = 0 ; ++num ; }
      else if ( sca && s_equal ( ascale + absp , ascale ) ) { m_pars [ k ] = 0 ; ++num ; }
    }
  //
  return num ;
}    
// ============================================================================


// ============================================================================
// join two vectors togather 
// ============================================================================
std::vector<double>
Ostap::Math::Parameters::join
( const std::vector<double>& a ,
  const std::vector<double>& b )
{
  //
  if      ( a.empty() ) { return b ; } 
  else if ( b.empty() ) { return a ; } 
  //
  std::vector<double> ab ( a.size () + b.size() ) ;
  //
  std::copy ( a.begin () , a.end () , ab.begin()            ) ;
  std::copy ( b.begin () , b.end () , ab.begin() + a.size() ) ;
  //
  return ab ;
}
// ============================================================================
// join scalar and vector together
// ============================================================================
std::vector<double>
Ostap::Math::Parameters::join
( const double               a ,
  const std::vector<double>& b )
{
  //
  if ( b.empty() ) { return { a } ; }
  //
  std::vector<double> ab ( 1 + b.size() ) ;
  //
  ab[0] = a ;
  std::copy ( b.begin () , b.end () , ab.begin() + 1 ) ;
  //
  return ab ;
}
// ============================================================================
// join 2 scalars and vector together
// ============================================================================
std::vector<double>
Ostap::Math::Parameters::join
( const double               a1 ,
  const double               a2 ,
  const std::vector<double>& b  )
{
  //
  if ( b.empty() ) { return { a1 , a2 } ; }
  //
  std::vector<double> ab ( 2 + b.size() ) ;
  //
  ab [ 0 ] = a1 ;
  ab [ 1 ] = a2  ;
  //
  std::copy ( b.begin () , b.end () , ab.begin() + 2 ) ;
  //
  return ab ;
}
// ============================================================================
// join 3 scalars and vector together
// ============================================================================
std::vector<double>
Ostap::Math::Parameters::join
( const double               a1 ,
  const double               a2 ,
  const double               a3 ,
  const std::vector<double>& b  )
{
  //
  if ( b.empty() ) { return { a1 , a2 , a3 } ; }
  //
  std::vector<double> ab ( 3 + b.size() ) ;
  //
  ab [ 0 ] = a1 ;
  ab [ 1 ] = a2  ;
  ab [ 2 ] = a3  ;
  //
  std::copy ( b.begin () , b.end () , ab.begin() + 3 ) ;
  //
  return ab ;
}
// ============================================================================
// join 4 scalars and vector together
// ============================================================================
std::vector<double>
Ostap::Math::Parameters::join
( const double               a1 ,
  const double               a2 ,
  const double               a3 ,
  const double               a4 ,
  const std::vector<double>& b  )
{
  //
  if ( b.empty() ) { return { a1 , a2 , a3 , a4 } ; }
  //
  std::vector<double> ab ( 4 + b.size() ) ;
  //
  ab [ 0 ] = a1 ;
  ab [ 1 ] = a2  ;
  ab [ 2 ] = a3  ;
  ab [ 3 ] = a4  ;
  //
  std::copy ( b.begin () , b.end () , ab.begin() + 4 ) ;
  //
  return ab ;
}
// ============================================================================
//  join vector & scalar together
// ============================================================================
std::vector<double>
Ostap::Math::Parameters::join
( const std::vector<double>& a , 
  const double               b )
{
  //
  if ( a.empty() ) { return { b } ; }
  //
  std::vector<double> ab ( a.size() + 1 ) ;
  //
  std::copy ( a.begin () , a.end () , ab.begin() ) ;
  ab [ a.size() ] = b ;
  //
  return ab ;
}
// ============================================================================
//  join vector & 2 scalars together
// ============================================================================
std::vector<double>
Ostap::Math::Parameters::join
( const std::vector<double>& a  , 
  const double               b1 , 
  const double               b2 ) 
{
  //
  if ( a.empty() ) { return { b1 , b2 } ; }
  //
  std::vector<double> ab ( a.size() + 2 ) ;
  //
  std::copy ( a.begin () , a.end () , ab.begin() ) ;
  ab [ a.size()     ] = b1 ;
  ab [ a.size() + 1 ] = b2 ;
  //
  return ab ;
}
// ============================================================================
//  join vector & 3 scalars together
// ============================================================================
std::vector<double>
Ostap::Math::Parameters::join
( const std::vector<double>& a  , 
  const double               b1 , 
  const double               b2 , 
  const double               b3 ) 
{
  //
  if ( a.empty() ) { return { b1 , b2 , b3 } ; }
  //
  std::vector<double> ab ( a.size() + 3 ) ;
  //
  std::copy ( a.begin () , a.end () , ab.begin() ) ;
  ab [ a.size()     ] = b1 ;
  ab [ a.size() + 1 ] = b2 ;
  ab [ a.size() + 2 ] = b3 ;
  //
  return ab ;
}
// ============================================================================
//  join vector & 4 scalars together
// ============================================================================
std::vector<double>
Ostap::Math::Parameters::join
( const std::vector<double>& a  , 
  const double               b1 , 
  const double               b2 , 
  const double               b3 , 
  const double               b4 ) 
{
  //
  if ( a.empty() ) { return { b1 , b2 , b3 , b4 } ; }
  // 
  std::vector<double> ab ( a.size() + 4 ) ;
  //
  std::copy ( a.begin () , a.end () , ab.begin() ) ;
  ab [ a.size()     ] = b1 ;
  ab [ a.size() + 1 ] = b2 ;
  ab [ a.size() + 2 ] = b3 ;
  ab [ a.size() + 4 ] = b4 ;
  //
  return ab ;
}




// ============================================================================
//                                                                      The END 
// ============================================================================
