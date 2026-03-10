// ============================================================================
#ifndef OSTAP_PARAMETER_H 
#define OSTAP_PARAMETER_H 1
// ============================================================================
// Include files
// ============================================================================
// STD& STL 
// ============================================================================
#include <cstdint>
#include <utility>
#include <vector>
#include <string>
#include <typeinfo>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/SetPar.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    // Scalar and transformed parameters 
    // ========================================================================
    /** @class Value
     *  Trivial (scalar) "Id"-parameter
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
     *  @date 2026-02-21
     */
    class Value
    {
    public :
      // ======================================================================
      /** full constructor
       *  @param value parameter value  
       *  @param name  parameter name  
       *  @param the_class name of the (owner/holder) class 
       */
      Value
      ( const double       value               ,
	const std::string& name      = "value" ,
	const std::string& the_class = ""      ) ;	
      // ======================================================================
      /** full constructor
       *  @param value parameter value  
       *  @param name  parameter name  
       *  @param the_class the (owner/holder) object
       */
      Value
      ( const double          value     , 
	const std::string&    name      ,
	const std::type_info& the_class ) ;
      // ======================================================================
      /** full constructor
       *  @param value parameter value  
       *  @param name  parameter name  
       *  @param the_class the (owner/holder) object
       */
      template <class CLASS>
      Value
      ( const double          value     , 
	const std::string&    name      ,
	const CLASS&          the_class )
	: Value ( value , name , typeid ( the_class ) )
      {}
      // ======================================================================      
    public :
      // ======================================================================
      /// get the value 
      inline double             value () const { return m_value ; }
      /// get the full parameter name  
      inline const std::string& name  () const { return m_name  ; }
      // ======================================================================
      /// the sign of the value
      std::int8_t               signum () const ;
      /// absolute value of parameter
      inline double             abs    () const { return std::abs ( m_value ) ; }
      // ======================================================================      
    public :
      // ======================================================================
      /// implicit conversion to double 
      inline operator double () const { return m_value ; } 
      /// set the value from double
      inline Value& operator=( const double value )
      { setValue ( value ) ; return *this ; }
      // ======================================================================
    public :
      // ======================================================================
      /// set new value for parameter 
      inline bool setValue
      ( const double value         ,
	const bool   force = false )
      { return set_par ( m_value , value , force ) ;  }
      // ======================================================================
    public :
      // ======================================================================
      /** set the full name of parameter
       *  @parm  the_name  the parameter  name
       *  @param the_class the name of ower/holder class
       */
      const std::string& setFullName
      ( const std::string& the_name  , 
	const std::string& the_class ) ;
      // ======================================================================
      /**  set the full name of parameter
       *  @parm  the_name  the parameter  name
       *  @param the_class the typeinfo of ower/holder class
       */
      const std::string& setFullName
      ( const std::string&    the_name  , 
	const std::type_info& the_class ) ;
      // ====================================================================
      /**  set the full name of parameter
       *  @parm  the_name  the parameter  name
       *  @param the_class the typeinfo of ower/holder class
       */
      template <class CLASS>
      inline const std::string& setFullName
      ( const std::string&    the_name  , 
	const CLASS&          the_class )
      { return this->setFullName ( the_name , typeid ( the_class ) ) ; }	      
      // ======================================================================
    public :
      // ======================================================================
      /// unique tag
      std::size_t tag () const ; 
      // ======================================================================			      
    private :
      // ======================================================================
      /// parameter 
      double      m_value { 0 } ; // parameter
      /// full parameter name
      std::string m_name  {}    ; // parameter name
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Transform
     *  Helper base class to represent the tranformed parameter
     *  - the "external" value: \f$ x \in \mathbb{R} \f$ 
     */
    class Transform
    {
      // ======================================================================
    protected : 
      // ======================================================================
      /** full constructor
       *  @param value parameter value  
       *  @param name  parameter name  
       *  @param the_class name of the (owner/holder) class 
       */      
      Transform 
      ( const double       value               ,
	const std::string& name      = "value" ,
	const std::string& the_class = ""      ) ;
      // ======================================================================
      /** full constructor
       *  @param value parameter value  
       *  @param name  parameter name  
       *  @param the_class the (owner/holder) object
       */
      Transform 
      ( const double          value     , 
	const std::string&    name      ,
	const std::type_info& the_class ) ;
      // ======================================================================
      /** full constructor
       *  @param value parameter value  
       *  @param name  parameter name  
       *  @param the_class the (owner/holder) object
       */
      template <class CLASS>
      Transform 
      ( const double          value     , 
	const std::string&    name      ,
	const CLASS&          the_class )
	:  Transform ( value , name , typeid ( the_class ) )
      {}           
      // ======================================================================
    public: // the most important methods from Ostap::Math::Value 
      // ======================================================================
      /// get the value 
      inline double             value     () const { return m_value.value  () ; }
      /// get the full parameter name  
      inline const std::string& name      () const { return m_value.name   () ; }
      /// the sign of the value
      std::int8_t               signum    () const { return m_value.signum () ; } 
      /// absolute value of parameter
      inline double             abs       () const { return m_value.abs    () ; }
      /// implicit conversion to double  
      inline operator           double    () const { return m_value.value  () ; }
      /// get the underlying value
      inline const Value&       the_value () const { return m_value           ; }
      // ======================================================================
    public: 
      // ======================================================================
      /// get external value
      inline double external    () const { return m_external ; }
      // ======================================================================
    public:
      // ======================================================================
      /** set the full name of parameter
       *  @parm  the_name  the parameter  name
       *  @param the_class the name of ower/holder class
       */
      inline const std::string& setFullName
      ( const std::string& the_name  , 
	const std::string& the_class )
      { return m_value.setFullName ( the_name , the_class ) ; } 
      // ======================================================================
      /**  set the full name of parameter
       *  @parm  the_name  the parameter  name
       *  @param the_class the typeinfo of ower/holder class
       */
      inline const std::string& setFullName
      ( const std::string&    the_name  , 
	const std::type_info& the_class ) 
      { return m_value.setFullName ( the_name , the_class ) ; } 
      // ====================================================================
      /**  set the full name of parameter
       *  @parm  the_name  the parameter  name
       *  @param the_class the typeinfo of ower/holder class
       */
      template <class CLASS>
      inline const std::string& setFullName
      ( const std::string&    the_name  , 
	const CLASS&          the_class )
      { return this->setFullName ( the_name , typeid ( the_class ) ) ; }	      
      // ======================================================================
    protected :
      // ======================================================================
      /// the "external" value: \f$ x \in \mathbb{R} \f$ 
      double              m_external { 0 } ; 
      /// the value of parameter 
      Ostap::Math::Value  m_value    { 0 } ;
      // ======================================================================      
    };
    // ========================================================================    
    /** @class LogValue
     *  Trivial (scalar) parameter: \f$ R \rigtharrow (0,+\infty) \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.c
     *  @date 2026-02-21
     */
    class LogValue : public Transform 
    {
    public :
      // ======================================================================
      /** full constructor
       *  @param value of parameter
       *  @param name  parameter name  
       *  @param the_class name of the (owner/holder) class 
       */
      LogValue
      ( const double       value     = 1       ,
	const std::string& name      = "value" ,
	const std::string& the_class = ""      ) ;	
      // ======================================================================
      /** full constructor
       *  @param value parameter value  
       *  @param name  parameter name  
       *  @param the_class the (owner/holder) object
       */
      LogValue
      ( const double          value     , 
	const std::string&    name      ,
	const std::type_info& the_class ) ;
      // ======================================================================
      /** full constructor
       *  @param value parameter value  
       *  @param name  parameter name  
       *  @param the_class the (owner/holder) object
       */
      template <class CLASS>
      LogValue
      ( const double          value     , 
	const std::string&    name      ,
	const CLASS&          the_class )
	: LogValue ( value , name , typeid ( the_class ) )
      {}
      // ======================================================================      
    public : // assignement from double 
      // ======================================================================
      /// set the value from double
      inline LogValue& operator=( const double value )
      { setValue ( value ) ; return *this ; }
      // ======================================================================
    public : 
      // ======================================================================
      /// set new value for  parameter
      bool setValue
      ( const double value         , 
	const bool   force = false ) ;
      // ======================================================================      
    public : // get & set internal/external values 
      // ======================================================================
      /// get log-value 
      inline double logValue   () const { return m_external ; }
      /// set new value for log-parameter 
      bool          setLogValue
      ( const double value         , 
	const bool   force = false ) ; 
      // ======================================================================
    public:
      // ======================================================================      
      /// set external value
      inline bool   setExternal
      ( const double value         , 
	const bool   force = false )
      { return setLogValue ( value , force ) ; }      
      // ======================================================================      
    public: 
      // ======================================================================
      /// unique tag
      std::size_t tag () const ; 
      // ======================================================================			      
    private :
      // ======================================================================
      /// logarithm of value 
      double m_logValue { 0 } ; // logarithm of value 
      /// the value itself 
      Value  m_value    { 1 } ; // the value itself 
      // ======================================================================
    } ;
    // ========================================================================
    /** @class InRange
     *  Trivial parameter with the limits \f$ \mathbb{R} \rightarrow \left[A,B\right]\f$
     *
     *  The actual transformation is
     *  \f[ p ( x ) = \frac{ B - A}{2} \left( 1 - \cos \frac{\pi\left(x-A){B-A}\right) + A \f]      
     *  for \f$ A < B \f$ 
     *  
     *  Such that
     *
     *   -  \f$  A \le p(x) \le B \f$ 
     *   -  \f$  p(A) = A \f$ 
     *   -  \f$  p(B) = B \f$
     *   -  \f$  p( \frac{A+B}{2} ) = \frac{A+B}{2} \f$
     *
     */
    class InRange : public Transform  
    {
    public :
      // =====================================================================
      /** full constructor
       *  @param value parameter value
       *  @param avalue A-value 
       *  @param bvalue B-value 
       *  @param name  parameter name  
       *  @param the_class name of the (owner/holder) class 
       */
      InRange
      ( const double       value     = 0 ,
        const double       avalue    = 0 , 
        const double       bvalue    = 1 , 
	const std::string& name      = "value" ,
	const std::string& the_class = ""      ) ;      
      // =====================================================================
      /** full constructor
       *  @param value parameter
       *  @param avalue A-value 
       *  @param bvalue B-value 
       *  @param name  parameter name  
       *  @param the_class name of the (owner/holder) class 
       */
      InRange
      ( const double          value     ,
        const double          avalue    , 
        const double          bvalue    , 
	const std::string&    name      ,      
	const std::type_info& the_class ) ;
      // =====================================================================
      /** full constructor
       *  @param value parameter value
       *  @param avalue A-value 
       *  @param bvalue B-value 
       *  @param name  parameter name  
       *  @param the_class name of the (owner/holder) class 
       */
      template <class CLASS>
      InRange      
      ( const double         value     ,
        const double         avalue    , 
        const double         bvalue    , 
	const std::string&   name      ,      
	const CLASS&         the_class )
	: InRange ( value  ,
		    avalue ,
		    bvalue ,
		    name   ,
		    typeid ( the_class ) )
      {}
      // =====================================================================
      /** full constructor
       *  @param avalue A-value 
       *  @param bvalue Bvalue 
       *  @param name  parameter name  
       *  @param the_class name of the (owner/holder) class 
       */
      InRange
      ( const double       avalue         , 
        const double       bvalue         , 
	const std::string& name           ,
	const std::string& the_class = "" ) ;      
      // =====================================================================
      /** full constructor
       *  @param avalue A-value 
       *  @param bvalue B-value 
       *  @param name  parameter name  
       *  @param the_class name of the (owner/holder) class 
       */
      InRange
      ( const double          avalue    , 
        const double          bvalue    , 
	const std::string&    name      ,      
	const std::type_info& the_class ) ;
      // =====================================================================
      /** full constructor
       *  @param avalue A-value 
       *  @param bvalue B-value 
       *  @param name  parameter name  
       *  @param the_class name of the (owner/holder) class 
       */
      template <class CLASS>
      InRange      
      ( const double         avalue    , 
        const double         bvalue    , 
	const std::string&   name      ,      
	const CLASS&         the_class )
	: InRange ( avalue ,
		    bvalue ,
		    name   ,
		    typeid ( the_class ) )
      {}      
      // =====================================================================
    public : 
      // ======================================================================
      /// set the value from double
      inline InRange& operator=( const double value )
      { setValue ( value ) ; return *this ; }
      // ======================================================================
    public : // get & set internal/external values
      // ======================================================================
      /// set new value for external -parameter 
      bool          setExternal
      ( const double value         ,
	const bool   force = false ) ;
      /// set new value for     parameter 
      bool          setValue
      ( const double value         ,
	const bool   force = false ) ;
      // ======================================================================
    public: 
      // ======================================================================
      /// minimal value 
      inline double vmin () const  { return m_min ; } 
      /// maximal  value 
      inline double vmax () const  { return m_max ; } 
      // ======================================================================      
    public: 
      // ======================================================================
      /// unique tag
      std::size_t tag () const ; 
      // ======================================================================
    protected : // variable transformation 
      // =====================================================================
      /// external -> internal 
      double t ( const double x ) const ;
      /// internal -> external 
      double x ( const double t ) const ;
      // ======================================================================
    public:
      // ======================================================================
      using Ostap::Math::Transform::value ;
      /// get the value as function external parameter
      inline double value       ( const double p ) const { return t ( p ) ; }
      /// get the value as function external parameter 
      inline double operator () ( const double p ) const { return t ( p ) ; }
      // ======================================================================
    private:
      // =====================================================================
      /// A-limit
      double m_min       {  0   } ; // A-limit
      /// B-limit
      double m_max       {  1   } ; // B-limit
      // =====================================================================
    } ;
    // ========================================================================
    /** @class InRange2
     *  Trivial parameter with the limits \f$ R \rightarrow (A,B) \f$
     *
     *  The actual transformation is 
     *  \f[ p ( x ) = \frac{B - A}{2} * \left ( 1 + \frac{2}{\pi} \atan k (x-x_0) \right) + A \f]
     *  for \f$ A < B \f$,
     *
     *  where :
     *  - \f$ x_0 = \frac{A+B}{2} \f$ 
     *  - \f$ \k = \frac{\pi}{B-A} f$ 
     *  
     *  Such that
     *
     *   -  \f$  A <  p ( y )  < B \f$
     *   -  \f$  p(-\infty)                = A \f$ 
     *   -  \f$  p(+\infty)                = B \f$     
     *   -  \f$  p( \frac{A+B}{2} )        = \frac{A+B}{2}\f$
     *   -  \f$  p^\prime( \frac{A+B}{2} ) = 1 \f$ 
     *
     */
    class InRange2 : public Transform 
    {
    public :
      // =====================================================================
      /** full constructor
       *  @param value parameter value
       *  @param avalue A-value 
       *  @param bvalue B-value 
       *  @param name  parameter name  
       *  @param the_class name of the (owner/holder) class 
       */
      InRange2
      ( const double       value     = 0 ,
        const double       avalue    = 0 , 
        const double       bvalue    = 1 , 
	const std::string& name      = "value" ,
	const std::string& the_class = ""      ) ;      
      // =====================================================================
      /** full constructor
       *  @param value parameter
       *  @param avalue A-value 
       *  @param bvalue B-value 
       *  @param name  parameter name  
       *  @param the_class name of the (owner/holder) class 
       */
      InRange2
      ( const double          value     ,
        const double          avalue    , 
        const double          bvalue    , 
	const std::string&    name      ,      
	const std::type_info& the_class ) ;
      // =====================================================================
      /** full constructor
       *  @param value parameter value
       *  @param avalue A-value 
       *  @param bvalue B-value 
       *  @param name  parameter name  
       *  @param the_class name of the (owner/holder) class 
       */
      template <class CLASS>
      InRange2      
      ( const double         value     ,
        const double         avalue    , 
        const double         bvalue    , 
	const std::string&   name      ,      
	const CLASS&         the_class )
	: InRange2 ( value  ,
		    avalue ,
		    bvalue ,
		    name   ,
		    typeid ( the_class ) )
      {}
      // =====================================================================
      /** full constructor
       *  @param avalue A-value 
       *  @param bvalue B-value 
       *  @param name  parameter name  
       *  @param the_class name of the (owner/holder) class 
       */
      InRange2
      ( const double       avalue         , 
        const double       bvalue         , 
	const std::string& name           ,
	const std::string& the_class = "" ) ;      
      // =====================================================================
      /** full constructor
       *  @param avalue A-value 
       *  @param bvalue B-value 
       *  @param name  parameter name  
       *  @param the_class name of the (owner/holder) class 
       */
      InRange2
      ( const double          avalue    , 
        const double          bvalue    , 
	const std::string&    name      ,      
	const std::type_info& the_class ) ;
      // =====================================================================
      /** full constructor
       *  @param avalue A-value 
       *  @param bvalue B-value 
       *  @param name  parameter name  
       *  @param the_class name of the (owner/holder) class 
       */
      template <class CLASS>
      InRange2      
      ( const double         avalue    , 
        const double         bvalue    , 
	const std::string&   name      ,      
	const CLASS&         the_class )
	: InRange2 ( avalue ,
		     bvalue ,
		     name   ,
		     typeid ( the_class ) )
      {}      
      // =====================================================================
    public : // conversion 
      // ======================================================================
      /// set the value from double
      inline InRange2& operator=( const double value )
      { setValue ( value ) ; return *this ; }
      // ======================================================================
    public : // get & set internal/external values
      // ======================================================================
      /// set new value for external -parameter 
      bool          setExternal
      ( const double value         ,
	const bool   force = false ) ;
      /// set new value for     parameter 
      bool          setValue
      ( const double value         ,
	const bool   force = false ) ;
      // ======================================================================
    public: 
      // ======================================================================
      /// minimal value 
      inline double vmin () const  { return m_min ; } 
      /// maximal  value 
      inline double vmax () const  { return m_max ; } 
      // ======================================================================      
    public: 
      // ======================================================================
      /// unique tag
      std::size_t tag () const ; 
      // ======================================================================
    protected : // variable transformation
      // =====================================================================
      /// external -> internal 
      double t ( const double x ) const ;
      /// internal -> external 
      double x ( const double t ) const ;
      // ======================================================================
    public:
      // ======================================================================
      using Ostap::Math::Transform::value ;
      /// get the value as function external parameter
      inline double value       ( const double p ) const { return t ( p ) ; }
      /// get the value as function external parameter 
      inline double operator () ( const double p ) const { return t ( p ) ; }
      // ======================================================================
    private:
      // =====================================================================
      /// A-limit
      double m_min       {  0   } ; // A-limit
      /// B-limit
      double m_max       {  1   } ; // B-limit
      // =====================================================================
    } ;    
    // ========================================================================
    /** @class Scale 
     *  Trivial scalar parameter scale != 0 
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.c
     *  @date 2026-02-21
     */
    class Scale : public Transform 
    {
    public :
      // ======================================================================
      /** full constructor
       *  @param value the value of scale parameter 
       *  @param name  the name of scale parameer 
       *  @param the_class name of the (owner/holder) class
       *  @param positive force the scale to be positive 
       */
      Scale
      ( const double       value               ,
	const std::string& name      = "scale" ,
	const std::string& the_class = ""      ,
	const bool         positive  = true    ) ;
      // ======================================================================
      /** full constructor
       *  @param value the value of scale parameter 
       *  @param name  the name of scale parameer 
       *  @param the_class the (owner/holder) object
       *  @param positive force the scale to be positive 
       */
      Scale
      ( const double          value     , 
	const std::string&    name      ,
	const std::type_info& the_class , 
	const bool            positive  = true  ) ;
      // ======================================================================
      /** full constructor
       *  @param value the value of scale parameter 
       *  @param name  the name of scale parameer 
       *  @param the_class the (owner/holder) object
       *  @param positive force the scale to be positive 
       */
      template <class CLASS>
      Scale
      ( const double          value     , 
	const std::string&    name      ,
	const CLASS&          the_class , 
	const bool            positive  = true  ) 
	: Scale ( value , name , typeid ( the_class ) , positive )
      {}
      // ======================================================================      
    public :
      // ======================================================================
      /// get the value 
      inline double             scale  () const { return m_value.value  () ; } 
      // ======================================================================
    public :
      // ======================================================================
      /// set the value from double
      inline Scale& operator=( const double value )
      { setValue ( value ) ; return *this ; }
      /// forced to be positive ?
      inline bool               positive () const { return m_positive       ; } 
      // ======================================================================
    public :
      // ======================================================================
      /// set new value for parameter 
      bool setValue
      ( const double value         ,
	const bool   force = false ) ;       
      /// set new value for exetrnal parameter 
      bool setExternal 
      ( const double value         ,
	const bool   force = false ) ;       
      // ======================================================================
      /// get the external parameter 
      inline double inv_value () const { return m_external ; }
      /// get the external parameter 
      inline double inv_scale () const { return m_external ; }
      /// set new value for exetrnal parameter 
      inline bool setInvScale 
      ( const double value         ,
	const bool   force = false )
      { return setExternal ( value , force ) ; }       
      /// set new value for exetrnal parameter 
      inline bool setInvValue 
      ( const double value         ,
	const bool   force = false )
      { return setExternal ( value , force ) ; }       
      // ======================================================================
    public: 
      // ======================================================================
      /// unique tag
      std::size_t tag () const ; 
      // ======================================================================			      
    private : 
      // ======================================================================
      /// forced to be positive ? 
      bool m_positive { true } ; // forced to be positive ? 
      // ======================================================================
    } ;    
    // ========================================================================
    /** @class ShiftAndScale
     *  Helper class to keep two parameters 
     *  - shift
     *  - scale 
     */
    class ShiftAndScale
    {
    public :
      // =====================================================================
      /** @param scale the value of scale parameter 
       *  @param shift the value of shift  parameter 
       *  @param scale_name the name of scale parameter
       *  @param shift_name the name of shift parameter
       *  @param the_class name of the (owner/holder) class 
       */
      ShiftAndScale
      ( const double       scale      = 1        ,
	const double       shift      = 0        ,
	const std::string& scale_name = "scale"  ,
	const std::string& shift_name = "shift"  ,
	const std::string& the_class  = ""       ,
	const bool         positive   = true     ) ;      
      // =====================================================================
      /** @param scale the value of scale parameter 
       *  @param shift the value of shift  parameter 
       *  @param scale_name the name of scale parameter
       *  @param shift_name the name of shift parameter
       *  @param the_class name of the (owner/holder) class 
       */
      ShiftAndScale
      ( const double          scale      ,
	const double          shift      ,
	const std::string&    scale_name ,
	const std::string&    shift_name ,
	const std::type_info& the_class  , 	       	
	const bool            positive   = true     ) ;      
      // =====================================================================
      /** @param scale the value of scale parameter 
       *  @param shift the value of shift  parameter 
       *  @param scale_name the name of scale parameter
       *  @param shift_name the name of shift parameter
       *  @param the_class name of the (owner/holder) class 
       */
      template <class CLASS>
      ShiftAndScale
      ( const double          scale      ,
	const double          shift      ,
	const std::string&    scale_name ,
	const std::string&    shift_name ,
	const CLASS&          the_class  , 
	const bool            positive = true )
	: ShiftAndScale ( scale      ,
			  shift      ,
			  scale_name ,
			  shift_name ,
			  typeid ( the_class ) , 
			  positive   )
      {}
      // =====================================================================
    public :
      // =====================================================================
      // get scale 
      inline double             scale      () const { return m_scale.value () ; }
      // get hift 
      inline double             shift      () const { return m_shift.value () ; }
      /// get scale name  
      inline const std::string& scale_name () const { return m_scale.name ()  ; }
      /// get shift name  
      inline const std::string& shift_name () const { return m_shift.name ()  ; }      
      // =====================================================================
      /** the sign of the value
       *  @see Ostap::Math::signum
       *  @see Ostap::Math::Value::signum
       */
      inline std::int8_t        scale_sign () const { return m_scale.signum () ; }
      /** the sign of the value
       *  @see Ostap::Math::signum
       *  @see Ostap::Math::Value::signum
       */      
      inline std::int8_t        sign_scale () const { return m_scale.signum () ; }
      /// absolute value of the scale
      inline double             scale_abs  () const { return m_scale.abs    () ; }
      /// absoluet value of the scale
      inline double             abs_scale  () const { return m_scale.abs    () ; }
      // =====================================================================
    public : 
      // =====================================================================
      /// scale variable 
      inline const Scale&       scale_var  () const { return m_scale ; }
      /// shift variable 
      inline const Value&       shift_var  () const { return m_shift ; }
      // =====================================================================      
    public: 
      // =====================================================================
      /// set new scale value
      inline bool setScale
      ( const double value         , 
	const bool   force = false )
      { return m_scale.setValue ( value , force ) ; }
      /// set new shift value
      inline bool setShift
      ( const double value         , 
	const bool   force = false )
      { return m_shift.setValue ( value , force ) ; }
      /// set both scale & shift 
      inline bool setScaleShift
      ( const double scale ,
	const double shift , 
	const bool   force = false )
      {
	const bool m1 = m_scale.setValue ( scale , force ) ;
	const bool m2 = m_shift.setValue ( shift , force ) ;
	return m1 && m2 ;
      }
      // =====================================================================					 
    public: 
      // ======================================================================
      /// x -> t transformation
      inline double t ( const double x ) const
      { return ( x - m_shift.value () ) / m_scale.value () ;  } 
      // ======================================================================			      
      /// t -> x transfrormation
      inline double x ( const double t ) const
      { return   t * m_scale.value ()   + m_shift.value () ;  } 
      // ======================================================================			      
    protected :
      // ======================================================================
      /// unique tag
      std::size_t tag () const ; 
      // ======================================================================			      
    protected :
      // =====================================================================
      /// scale parameter 
      Scale m_scale { 1 } ; // scale parameter 
      /// shift parameter 
      Value m_shift { 0 } ; // shift parameter 
      // =====================================================================
    } ;
    // ========================================================================
    /** @class PQ
     *  keep two positive variables P & Q
     *  - variables positive and defied via their logarithms
     *  - cache \f$  ln \Beta ( p , q )         \f$
     *  - cache \f$  \frac{1}{ \Beta ( p , q )} \f$
     *  @see Ostap::Math::AB 
     *  @see Ostap::Math::ibeta 
     *  @see Ostap::Math::lnbeta
     */
    class PQ
    {
      // ======================================================================
    public :
      // ======================================================================
      /** @param p     p-parameter 
       *  @param q     q-parameter 
       *  @param pname the name of p-parameter
       *  @param qname the name of q-parameter
       *  @param the_class name of the (owner/holder) class 
       */
      PQ
      ( const double       p         =  1  ,
	const double       q         =  1  ,
	const std::string& pname     = "p" ,
	const std::string& qname     = "q" ,
	const std::string& the_class = ""  ) ;
      // ======================================================================
      /** @param p     p-parameter 
       *  @param q     q-parameter 
       *  @param pname the name of p-parameter
       *  @param qname the name of q-parameter
       *  @param the_class name of the (owner/holder) class
       */
      PQ      
      ( const double          p         ,
	const double          q         ,
	const std::string&    pname     ,
	const std::string&    qname     ,
	const std::type_info& the_class ) ;
      // ======================================================================      
      /** @param p     p-parameter 
       *  @param q     q-parameter 
       *  @param pname the name of p-parameter
       *  @param qname the name of q-parameter
       *  @param the_class name of the (owner/holder) class
       */
      template <class CLASS>
      PQ
      ( const double          p         ,
	const double          q         ,
	const std::string&    pname     ,
	const std::string&    qname     ,
	const CLASS&          the_class )
	: PQ ( p , q , pname , qname , typeid ( the_class ) )
      {}      
      // =====================================================================
    public :
      // =====================================================================
      // get p 
      inline double             p     () const { return m_p.value () ; }
      // get q 
      inline double             q     () const { return m_q.value () ; }
      /// get log-P
      inline double             logP  () const { return m_p.logValue () ; } 
      /// get log-Q
      inline double             logQ  () const { return m_q.logValue () ; } 
      /// get p-name  
      inline const std::string& pname () const { return m_p.name () ; }
      /// get q-name  
      inline const std::string& qname () const { return m_q.name () ; }
      // =====================================================================
    public : 
      // =====================================================================
      /// p-parameter 
      inline const LogValue&    pvar  () const { return m_p ; }
      /// q-parameter 
      inline const LogValue&    qvar  () const { return m_q ; }
      // =====================================================================      
    public:   /// setters 
      // =====================================================================
      bool setP    ( const double value , const bool force = false ) ;
      bool setQ    ( const double value , const bool force = false ) ;
      //
      bool setLogP ( const double value , const bool force = false ) ; 
      bool setLogQ ( const double value , const bool force = false ) ; 
      /// set two parameters at once
      inline bool setPQ
      ( const double pvalue ,
	const double qvalue ,
	const bool   force  = false )
      {
	const bool mp = setP ( pvalue , force ) ;
	const bool mq = setQ ( qvalue , force ) ;
	return mp && mq ;
      }
      /// set two parameters at once
      inline bool setLogPQ
      ( const double pvalue ,
	const double qvalue , 
	const bool   force  = false )
      {
	const bool mp = setLogP ( pvalue , force ) ;
	const bool mq = setLogQ ( qvalue , force ) ;
	return mp && mq ;
      } 
      // =====================================================================      
    public: 
      // ======================================================================
      /// unique tag
      std::size_t tag () const ; 
      // ======================================================================			      
    protected : 
      // ======================================================================
      /// get cached value of ln B(p,q)
      inline double log_Beta () const { return m_log_Beta ; } 
      /// get cached value of 1/B(p,q)
      inline double inv_Beta () const { return m_inv_Beta ; } 
      // ======================================================================			      
    private :
      // ======================================================================
      /// parameter P 
      LogValue m_p      { 1 } ; // parameter P 
      /// parameter Q 
      LogValue m_q      { 1 } ; // parameter Q
      /// cached value of \f$ \log \Beta (p,q \f$ 
      double m_log_Beta { 0 } ; // log Beta (p,q) 
      /// cached value of \f$ \frac{1} { \Beta (p,q } \f$ 
      double m_inv_Beta { 0 } ; // inv Beta (p,q) 
      // ======================================================================      
    } ;    
    // ========================================================================
    /** @class AsymVars
     *  Hold two same-sign varables and their asymemtries
     *  - var1
     *  - var2
     *  - average :   var1 + var2 
     *  - kappa   : ( var1 - var2 ) / ( var1 + var2 )
     *  - psi     : kappa = tanh ( psi )
     */
    class AsymVars
    {
    public:
      // ======================================================================
      /** full constructor
       *  @param   var1   the first variable 
       *  @param   var2   the second variable
       *  @param   v1name the name of the first variable 
       *  @param   v2name the name of the second variable
       *  @param the_class name of the (owner/holder) class 
       */       
      AsymVars
      ( const double          var1       = 1      ,
	const double          var2       = 1      ,
	const std::string&    v1name     = "var1" ,
	const std::string&    v2name     = "var2" ,
	const std::string&    the_class  = ""     ) ;
      // =====================================================================
      /** full constructor
       *  @param   var1   the first variable 
       *  @param   var2   the second variable
       *  @param   v1name the name of the first variable 
       *  @param   v2name the name of the second variable
       *  @param   the_class name of the (owner/holder) class 
       */       
      AsymVars
      ( const double          var1      ,
	const double          var2      ,
	const std::string&    v1name    ,
	const std::string&    v2name    ,
      	const std::type_info& the_class ) ;
      // =====================================================================
      /** full constructor
       *  @param   var1   the first variable 
       *  @param   var2   the second variable
       *  @param   v1name the name of the first variable 
       *  @param   v2name the name of the second variable
       *  @param   the_class name of the (owner/holder) class 
       */
      template <class CLASS>
      AsymVars
      ( const double          var1      ,
	const double          var2      ,
	const std::string&    v1name    ,
	const std::string&    v2name    ,
      	const CLASS&          the_class )
	: AsymVars ( var1 , var2 , v1name , v2name , typeid ( the_class ) )
      {}
      // =====================================================================
    public:
      // ======================================================================
      inline double var1    () const { return m_var1  .value () ; } 
      inline double var2    () const { return m_var2  .value () ; }
      inline double mean    () const { return m_mean  .value () ; }
      inline double kappa   () const { return m_kappa .value () ; }
      inline double psi     () const { return m_psi   .value () ; }      
      // =====================================================================
      inline const std::string& v1name     () const { return m_var1 .name () ; } 
      inline const std::string& v2name     () const { return m_var2 .name () ; } 
      inline const std::string& mean_name  () const { return m_mean .name () ; } 
      inline const std::string& kappa_name () const { return m_kappa.name () ; } 
      inline const std::string& psi_name   () const { return m_psi  .name () ; } 
      // =====================================================================
    public:
      // =====================================================================
      /// set two variables simultaneously 
      bool setVars
      ( const double v1 ,
	const double v2 ) ;
      /// set two Mean/Average & kappa 
      bool setMeanKappa 
      ( const double vmean  ,
	const double vkappa ) ;
      /// set two Mean/Average & psi 
      bool setMeanPsi  
      ( const double vmean  ,
	const double vpsi   ) ;      
      // =====================================================================
    public:
      // =====================================================================
      inline const Ostap::Math::Value& var_1     () const { return m_var1  ; }
      inline const Ostap::Math::Value& var_2     () const { return m_var2  ; }
      inline const Ostap::Math::Scale& var_mean  () const { return m_mean  ; }
      inline const Ostap::Math::Value& var_kappa () const { return m_kappa ; }
      inline const Ostap::Math::Value& var_psi   () const { return m_psi   ; }
      // =====================================================================      
    public :
      // ======================================================================
      /// unique tag
      std::size_t tag () const ; 
      // ======================================================================			      
    private :
      // ======================================================================
      /// the first  variable 
      Ostap::Math::Value m_var1  { 1 } ; // the first variable
      /// the second variable       
      Ostap::Math::Value m_var2  { 1 } ; // the second variable 
      /// mean value:
      Ostap::Math::Scale m_mean  { 1 } ; // average 
      /// asymmetry: kappa 
      Ostap::Math::Value m_kappa { 0 } ; // asymmetry kappa 
      /// asymmetry: psi
      Ostap::Math::Value m_psi   { 0 } ; // asymmetry psi 
      // =====================================================================
    } ;
    // ========================================================================
  } //                                             end of namespace Ostap::Math 
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_PARAMETER_H
// ============================================================================
