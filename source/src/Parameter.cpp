 // ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================ 
#include <algorithm>
#include <string>
#include <string>
#include <typeinfo>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Hash.h"
#include "Ostap/Names.h"
#include "Ostap/Math.h"
#include "Ostap/MoreMath.h"
#include "Ostap/Parameter.h"
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
namespace
{
  // ==========================================================================
  static const std::string s_invalid_value { "Invalid value!"         } ;
  static const std::string s_invalid_log   { "Invalid log/exp-value!" } ;    
  static const std::string s_invalid_range { "Invalid min/max-range!" } ;  
  static const std::string s_not_same_sign { "Variables must have the same sign!" } ;  
  // ==========================================================================
} // the end of anonymous namespace 
// ===========================================================================
/*  @fn set_par
 *  Helper function to set parameter
 *  @param parameter  (update) parameter
 *  @param value      (input)  new value
 *  @param force      (input)  froce updating
 *  @return true if parameter is modified
 */
// ===========================================================================
bool Ostap::Math::set_par
( double&      parameter ,
  const double value     ,
  const bool   force     )
{
  if ( !force && s_equal ( value , parameter ) ) { return false ; }
  parameter = value ;
  return true ;  
}
// ===========================================================================
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
/*  the sign of the value
 *  @see Ostap::Math::signum
 */
// ============================================================================
std::int8_t Ostap::Math::Value::signum () const
{ return Ostap::Math::signum ( m_value ) ;  }
// ============================================================================
/** set the full name of parameter
 *  @parm  the_name  the parameter  name
 *  @param the_class the name of ower/holder class
 */
// ============================================================================
const std::string&
Ostap::Math::Value::setFullName
( const std::string& the_name  , 
  const std::string& the_class ) 
{
  //
  const std::string c1 { Ostap::strip ( the_name  ) } ;
  const std::string c2 { Ostap::strip ( the_class ) } ;
  //
  if       ( !c1.empty() && !c2.empty() ) { m_name = c2 + "::" + c1 ; }
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
( const std::string&    the_name  , 
  const std::type_info& the_class ) 
{ return setFullName ( the_name , Ostap::class_name ( the_class ) ) ; }
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
 *  @param value of parameter
 *  @param name  parameter name  
 *  @param the_class name of the (owner/holder) class 
 */
// ============================================================================
Ostap::Math::LogValue::LogValue
( const double       value     ,
  const std::string& name      ,
  const std::string& the_class )
  : m_logValue ( 0 )
  , m_value    ( value , name , the_class )
{
  //
  Ostap::Assert ( s_EXP_UNDERFLOW_EXP < value && value < s_EXP_OVERFLOW_EXP ,		  
		  s_invalid_value                            ,
		  m_value.name ()                            , 
		  INVALID_LOGPARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_logValue = std::log   ( value ) ;  
  m_value.setValue ( value ) ;
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
bool Ostap::Math::LogValue::setLogValue
( const double log_value ,   
  const bool   force     ) 
{
  if ( !force && s_equal ( log_value , m_logValue ) ) { return false ; }
  //
  Ostap::Assert ( s_EXP_UNDERFLOW < log_value && log_value < s_EXP_OVERFLOW ,		  
		  s_invalid_log                               ,
		  m_value.name ()                             , 
		  INVALID_LOGPARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_logValue = log_value ;  
  return m_value.setValue ( std::exp ( m_logValue ) , force ) ;
}
// =====================================================================
// set new value for     parameter 
// =====================================================================
bool Ostap::Math::LogValue::setValue
( const double value ,
  const bool   force ) 
{
  if ( !force && s_equal ( m_value.value() , value ) ) { return false ; }
  //
  Ostap::Assert ( s_EXP_UNDERFLOW_EXP < value && value < s_EXP_OVERFLOW_EXP ,		  
		  s_invalid_value                            ,
		  m_value.name ()                            , 
		  INVALID_LOGPARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_logValue = std::log   ( value ) ;  
  return m_value.setValue ( value , force ) ;
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
( const double       value     , 
  const double       avalue    ,
  const double       bvalue    , 
  const std::string& name      ,
  const std::string& the_class )
  : m_min       ( ) 
  , m_max       ( ) 
  , m_external  ( )
  , m_value     ( value , name   , the_class ) 
{
  //
  Ostap::Assert ( !s_equal ( avalue , bvalue ) && std::isfinite ( avalue ) && std::isfinite ( bvalue ) , 
		  s_invalid_range  ,
		  m_value.name ()  ,		  
		  INVALID_RANGE    , __FILE__ , __LINE__  ) ;
  //
  m_min = std::min ( avalue , bvalue ) ;
  m_max = std::max ( avalue , bvalue ) ;
  //
  Ostap::Assert ( ( m_min <= value && value <= m_max )
		  || s_equal ( m_min , value )
		  || s_equal ( m_max , value ) , 
		  s_invalid_value ,
		  m_value.name () ,		  
		  INVALID_RANGE   , __FILE__ , __LINE__  ) ;
  //
  if      ( value != m_min && s_equal ( m_min , value ) ) { m_value.setValue ( m_min , true ) ; }
  else if ( value != m_max && s_equal ( m_max , value ) ) { m_value.setValue ( m_max , true ) ; }
  //
  /// internal -> external 
  m_external = t ( m_value.value () ) ;    
}
// ============================================================================
/*  full constructor
 *  @param value parameter value
 *  @param avalue A0value 
 *  @param bvalue B-value 
 *  @param name  parameter name  
 *  @param the_class name of the (owner/holder) class 
 */
// ============================================================================
Ostap::Math::InRange::InRange
( const double          value     , 
  const double          avalue    ,
  const double          bvalue    , 
  const std::string&    name      ,
  const std::type_info& the_class )
  : InRange ( value  ,
	      avalue ,
	      bvalue ,
	      name   ,
	      Ostap::class_name ( the_class ) )
{}
// ============================================================================
/*  full constructor
 *  @param avalue A-value 
 *  @param bvalue Bvalue 
 *  @param name  parameter name  
 *  @param the_class name of the (owner/holder) class 
 */
// ============================================================================
Ostap::Math::InRange::InRange
( const double       avalue    , 
  const double       bvalue    , 
  const std::string& name      ,
  const std::string& the_class )
  : InRange ( 0.5 * ( avalue + bvalue ) ,
	      avalue    ,
	      bvalue    ,
	      name      ,
	      the_class )
{}
// ============================================================================
/*  full constructor
 *  @param avalue A-value 
 *  @param bvalue B-value 
 *  @param name  parameter name  
 *  @param the_class name of the (owner/holder) class 
 */
// ============================================================================
Ostap::Math::InRange::InRange
( const double          avalue    , 
  const double          bvalue    , 
  const std::string&    name      ,      
  const std::type_info& the_class ) 
  : InRange ( 0.5 * ( avalue + bvalue ) ,
	      avalue    ,
	      bvalue    ,
	      name      ,
	      the_class )
{}
// ============================================================================
// external -> internal 
// ============================================================================
double Ostap::Math::InRange::t
( const double x ) const
{
  //
  const long double dx = 0.5L   * ( m_max - m_min )      ;
  const long double xx = s_pi_2 * ( x     - m_min ) / dx ;
  const long double tt = 1.0L - std::cos ( xx ) ;
  //
  return m_min + dx * tt ;
}  
// ============================================================================
// internal -> external 
// ============================================================================
double Ostap::Math::InRange::x
( const double t ) const
{
  //
  Ostap::Assert ( m_min <= t && t <= m_max , 
		  s_invalid_range          ,
		  m_value.name      ()     ,
		  INVALID_RANGE , __FILE__ , __LINE__ ) ;
  //
  const long double dx = 0.5L * ( m_max - m_min   ) ;
  
  const long double t1 = 1.0L - ( 1.0L * t - m_min ) / dx ;
  const long double t2 = std::acos  ( t1 ) * s_2_pi       ;
  //
  return dx * t2 + m_min ;
}
// =====================================================================
// set new value for parameter 
// ============================================================================
bool Ostap::Math::InRange::setExternal
( const double value ,
  const bool   force )
{
  if ( !force && s_equal ( value , m_external ) ) { return false ; }
  m_external = value ;  
  return m_value.setValue ( t ( m_external ) , force ) ;
}
// =====================================================================
// set new value for     parameter 
// =====================================================================
bool Ostap::Math::InRange::setValue
( const double value ,
  const bool   force )
{
  if ( !m_value.setValue ( value , force ) ) { return false ; }
  m_external  = x ( value ) ;
  return true ;
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
 *  @param value parameter value (external)
 *  @param avalue A0value 
 *  @param bvalue B-value 
 *  @param name  parameter name  
 *  @param the_class name of the (owner/holder) class 
 */
// ============================================================================
Ostap::Math::InRange2::InRange2
( const double       value     , 
  const double       avalue    ,
  const double       bvalue    , 
  const std::string& name      ,
  const std::string& the_class )
  : m_min       ( ) 
  , m_max       ( ) 
  , m_external  ( )
  , m_value     ( value , name   , the_class ) 
{
  //
  Ostap::Assert ( !s_equal ( avalue , bvalue ) && std::isfinite ( avalue ) && std::isfinite ( bvalue ) , 
		  s_invalid_range  ,
		  m_value.name ()  ,		  
		  INVALID_RANGE    , __FILE__ , __LINE__  ) ;
  //
  m_min = std::min ( avalue , bvalue ) ;
  m_max = std::max ( avalue , bvalue ) ;
  //
  Ostap::Assert ( ( m_min < value && value < m_max )
		  || s_equal ( m_min , value )
		  || s_equal ( m_max , value ) , 
		  s_invalid_value ,
		  m_value.name () ,		  
		  INVALID_RANGE   , __FILE__ , __LINE__  ) ;
  //
  /// internal -> external 
  m_external = t ( m_value.value () ) ;    
}
// ============================================================================
/*  full constructor
 *  @param value parameter value
 *  @param avalue A0value 
 *  @param bvalue B-value 
 *  @param name  parameter name  
 *  @param the_class name of the (owner/holder) class 
 */
// ============================================================================
Ostap::Math::InRange2::InRange2
( const double          value     , 
  const double          avalue    ,
  const double          bvalue    , 
  const std::string&    name      ,
  const std::type_info& the_class )
  : InRange2 ( value  ,
	       avalue ,
	       bvalue ,
	       name   ,
	       Ostap::class_name ( the_class ) )
{}
// ============================================================================
/*  full constructor
 *  @param avalue A-value 
 *  @param bvalue Bvalue 
 *  @param name  parameter name  
 *  @param the_class name of the (owner/holder) class 
 */
// ============================================================================
Ostap::Math::InRange2::InRange2
( const double       avalue    , 
  const double       bvalue    , 
  const std::string& name      ,
  const std::string& the_class )
  : InRange2 ( 0.5 * ( avalue + bvalue ) ,
	       avalue    ,
	       bvalue    ,
	       name      ,
	       the_class )
{}
// ============================================================================
/*  full constructor
 *  @param avalue A-value 
 *  @param bvalue B-value 
 *  @param name  parameter name  
 *  @param the_class name of the (owner/holder) class 
 */
// ============================================================================
Ostap::Math::InRange2::InRange2
( const double          avalue    , 
  const double          bvalue    , 
  const std::string&    name      ,      
  const std::type_info& the_class ) 
  : InRange2 ( 0.5 * ( avalue + bvalue ) ,
	       avalue    ,
	       bvalue    ,
	       name      ,
	       the_class )
{}
// ============================================================================
// external -> internal 
// ============================================================================
double Ostap::Math::InRange2::t
( const double x ) const
{
  const long double x0 = 0.5L   * ( m_min + m_max ) ;
  const long double dx = 0.5L   * ( m_max - m_min ) ;
  //
  const long double xx = s_pi_2 * ( x - x0 ) / dx ;  
  const long double tt = 1.0L + s_2_pi * std::atan ( xx ) ;
  //
  return m_min + dx * tt ;
}  
// ============================================================================
// internal -> external 
// ============================================================================
double Ostap::Math::InRange2::x
( const double t ) const
{
  //
  Ostap::Assert ( m_min < t && t < m_max ,
		  s_invalid_range        ,
		  m_value.name      ()   ,
		  INVALID_RANGE , __FILE__ , __LINE__ ) ;
  //
  const long double dx = 0.5L * ( m_max - m_min   ) ;
  const long double x0 = 0.5L * ( m_max + m_min   ) ;
  const long double t1 = ( 1.0L * t - m_min ) / dx  - 1.0L ;
  const long double t2 = std::tan ( s_pi_2 * t1 ) * s_2_pi * dx ;
  //
  return t2 + x0 ;
}
// =====================================================================
// set new value for parameter 
// =====================================================================
bool Ostap::Math::InRange2::setExternal
( const double value ,
  const bool   force )
{
  if ( !force && s_equal ( value , m_external ) ) { return false ; }
  m_external = value ;  
  return m_value.setValue ( t ( m_external ) , force ) ;
}
// =====================================================================
// set new value for     parameter 
// =====================================================================
bool Ostap::Math::InRange2::setValue
( const double value ,
  const bool   force )
{
  if ( !m_value.setValue ( value , force ) ) { return false ; }
  m_external  = x ( value ) ;  
  return true ;
}
// ============================================================================
// Unique value 
// ============================================================================
std::size_t Ostap::Math::InRange2::tag () const 
{
  static const std::string s_name { "InRange2" } ;
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
bool Ostap::Math::Scale::setValue
( const double value ,
  const bool   force )
{
  const double new_value = m_positive ? std::abs ( value ) : value ;
  if ( !force && s_equal ( new_value , m_scale.value () ) ) { return false ; }
  //
  static const std::string s_invalid_scale { "Invalid scale-value!" } ;
  Ostap::Assert ( !s_zero ( new_value) ,		  
		  s_invalid_scale      ,                               
		  m_scale.name ()      , 
		  INVALID_SCALEPARAMETER , __FILE__ , __LINE__ ) ;
  //
  return m_scale.setValue ( new_value , force ) ; 
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
/*  @param p     p-parameter p>0
 *  @param q     q-parameter q>0
 *  @param pname the name of p-parameter
 *  @param qname the name of q-parameter
 *  @param the_class name of the (owner/holder) class 
 */
// ============================================================================
Ostap::Math::PQ::PQ
( const double       p         ,
  const double       q         ,
  const std::string& pname     ,
  const std::string& qname     ,
  const std::string& the_class )
  : m_p ( p , pname , the_class )
  , m_q ( q , qname , the_class )
{
  m_log_Beta = Ostap::Math::lnbeta ( m_p.value () , m_q.value () ) ;  
  m_inv_Beta = Ostap::Math::ibeta  ( m_p.value () , m_q.value () ) ;  
}
// ============================================================================
/*  @param p     p-parameter p>0
 *  @param q     q-parameter q>0
 *  @param pname the name of p-parameter
 *  @param qname the name of q-parameter
 *  @param the_class name of the (owner/holder) class 
 */
// ============================================================================
Ostap::Math::PQ::PQ
( const double          p         ,
  const double          q         ,
  const std::string&    pname     ,
  const std::string&    qname     ,
  const std::type_info& the_class )
  : PQ ( p , q , pname , qname , Ostap::class_name ( the_class ) )
{}
// ============================================================================
// update log-P
// ============================================================================
bool Ostap::Math::PQ::setLogP
( const double value , 
  const bool   force ) 
{
  if ( !m_p.setLogValue ( value , force ) ) { return false ; } 
  m_log_Beta = Ostap::Math::lnbeta ( m_p.value () , m_q.value () ) ;
  m_inv_Beta = Ostap::Math::ibeta  ( m_p.value () , m_q.value () ) ;  
  return true ;
} 
// ============================================================================
// update log-Q
// ============================================================================
bool Ostap::Math::PQ::setLogQ
( const double value , 
  const bool   force )
{
  if ( !m_p.setLogValue ( value , force ) ) { return false ; } 
  m_log_Beta = Ostap::Math::lnbeta ( m_p.value () , m_q.value () ) ;
  m_inv_Beta = Ostap::Math::ibeta  ( m_p.value () , m_q.value () ) ;  
  return true ;
} 
// ============================================================================
// update P
// ============================================================================
bool Ostap::Math::PQ::setP
( const double value , 
  const bool   force ) 
{
  if ( !m_p.setValue ( value , force ) ) { return false ; } 
  m_log_Beta = Ostap::Math::lnbeta ( m_p.value () , m_q.value () ) ;
  m_inv_Beta = Ostap::Math::ibeta  ( m_p.value () , m_q.value () ) ;  
  return true ;
} 
// ============================================================================
// update Q
// ============================================================================
bool Ostap::Math::PQ::setQ
( const double value , 
  const bool   force ) 
{
  if ( !m_q.setValue ( value , force ) ) { return false ; }  
  m_log_Beta = Ostap::Math::lnbeta ( m_p.value () , m_q.value () ) ;
  m_inv_Beta = Ostap::Math::ibeta  ( m_p.value () , m_q.value () ) ;  
  return true ;
}
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
/*  full constructor
 *  @param   var1   the first variable 
 *  @param   var2   the second variable
 *  @param   v1name the name of the first variable 
 *  @param   v2name the name of the second variable
 *  @param the_class name of the (owner/holder) class 
 */       
// ============================================================================
Ostap::Math::AsymVars::AsymVars
( const double       var1      ,
  const double       var2      ,
  const std::string& v1name    ,
  const std::string& v2name    ,
  const std::string& the_class )
  : m_var1  ( var1                  , v1name  , the_class )
  , m_var2  ( var2                  , v2name  , the_class )
  , m_mean  ( 0.5 * ( var1 + var2 ) , "mean"  , the_class , false )
  , m_kappa ( 0.0                   , "kappa" , the_class )
  , m_psi   ( 0.0                   , "psi"   , the_class )
{
  //
  const double v1 = m_var1.value () ;
  const double v2 = m_var2.value () ;
  //
  Ostap::Assert ( Ostap::Math::same_sign ( v1 , v2 ) ,
		  s_not_same_sign   ,
		  m_mean.name ()    ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  const double k  = ( v1 - v2 ) / ( v1 + v2 ) ;
  //  
  m_kappa.setValue ( k ) ;
  //
  const double p  = std::atanh ( k ) ;
  m_psi  .setValue ( p ) ;
  //    
}
// ============================================================================
/*  full constructor
 *  @param   var1   the first variable 
 *  @param   var2   the second variable
 *  @param   v1name the name of the first variable 
 *  @param   v2name the name of the second variable
 *  @param the_class name of the (owner/holder) class 
 */       
// ============================================================================
Ostap::Math::AsymVars::AsymVars
( const double          var1      ,
  const double          var2      ,
  const std::string&    v1name    ,
  const std::string&    v2name    ,
  const std::type_info& the_class )
  : AsymVars ( var1   ,
	       var2   ,
	       v1name ,
	       v2name ,
	       Ostap::class_name ( the_class ) )
{}
// ============================================================================
// set two variables simultaneously 
// ============================================================================
bool Ostap::Math::AsymVars::setVars
( const double var1 ,
  const double var2 )
{
  const bool m1 = m_var1.setValue ( var1 ) ;
  const bool m2 = m_var2.setValue ( var2 ) ;
  if ( !m1 && !m2 ) { return false ; }
  //
  const double v1 = m_var1.value () ;
  const double v2 = m_var2.value () ;
  //
  Ostap::Assert ( Ostap::Math::same_sign ( v1 , v2 ) ,
		  s_not_same_sign   ,
		  m_mean.name ()    ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //    
  const double k  = ( v1 - v2 ) / ( v1 + v2 ) ;
  m_kappa.setValue ( k ) ;
  //
  const double p  = std::atanh ( k ) ;
  m_psi  .setValue ( p ) ;
  //    
  return true ;  
}
// ============================================================================
// set two Mean/Average & kappa
// ============================================================================
bool Ostap::Math::AsymVars::setMeanKappa 
( const double vmean  ,
  const double vkappa )
{
  Ostap::Assert ( -1 < vkappa && vkappa < 1 ,
		  s_invalid_range ,
		  m_kappa.name () ,
		  INVALID_RANGE   , __FILE__ , __LINE__  ) ;
  //
  const bool mm = m_mean .setValue ( vmean  ) ;
  const bool mk = m_kappa.setValue ( vkappa ) ;
  if ( !mm && !mk  ) { return false ; }
  //
  const double mv = m_mean .value () ;
  const double k  = m_kappa.value () ;
  //
  const double v1 = mv * ( 1 + k ) ;
  const double v2 = mv * ( 1 - k ) ;
  //
  const double m1 = m_var1.setValue ( v1 ) ;
  const double m2 = m_var2.setValue ( v2 ) ;
  if ( !m1 && !m2  ) { return false ; }
  //
  const double p  = std::atanh ( k ) ;
  Ostap::Assert ( std::isfinite ( p ) , 
		  s_invalid_range     ,
		  m_psi.name ()       ,
		  INVALID_RANGE       , __FILE__ , __LINE__  ) ;
  //
  if ( !m_psi .setValue ( p ) ) { return false ; }
  //
  return true ; 
}
// ============================================================================
// set two Mean/Average & psi
// ============================================================================
bool Ostap::Math::AsymVars::setMeanPsi  
( const double vmean  ,
  const double vpsi   )
{
  //
  const bool   mm = m_mean .setValue ( vmean ) ;
  const bool   mp = m_psi  .setValue ( vpsi  ) ;
  if ( !mm && !mp  ) { return false ; }
  //
  const double mv = m_mean .value () ;
  const double vp = m_psi  .value () ;
  //
  const double k  = std::tanh ( vp ) ;
  const double mk = m_kappa.setValue ( k ) ;
  //
  const double v1 = mv * ( 1 + k ) ;
  const double v2 = mv * ( 1 - k ) ;
  //
  const double m1 = m_var1.setValue ( v1 ) ;
  const double m2 = m_var2.setValue ( v2 ) ;
  if ( !m1 && !m2 && !mk ) { return false ; }
  //
  return true ;
}
// ============================================================================
// Unique value 
// ============================================================================
std::size_t Ostap::Math::AsymVars::tag () const 
{
  static const std::string s_name { "AsymVars" } ;
  return Ostap::Utils::hash_combiner ( s_name        ,
				       m_var1.tag () ,
				       m_var2.tag () ) ; 
}

// ============================================================================
//                                                                      The END 
// ============================================================================
