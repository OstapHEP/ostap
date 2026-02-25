// ============================================================================
#ifndef OSTAP_PARAMETERS_H 
#define OSTAP_PARAMETERS_H 1
// ============================================================================
// Include files
// ============================================================================
// STD& STL 
// ============================================================================
#include <initializer_list>
#include <vector>
#include <algorithm>
#include <typeinfo>
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
     *  Trivial (scalar) Id-parameter
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.c
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
      bool setValue
      ( const double value         ,
	const bool   force = false ) ; 
      // ======================================================================
    public :
      // ======================================================================
      /** set the full name of parameter
       *  @param the_class the name of owner/holder class
       *  @parm  the_name  the parameter  name
       */
      const std::string&
      setFullName
      ( const std::string&    the_class ,
	const std::string&    the_name  ) ;
      // ======================================================================			
      /** set the full name of parameter
       *  @param the_class the type-info of owner/holder class
       *  @parm  the_name  the parameter  name
       */
      const std::string&
      setFullName
      ( const std::type_info& the_class ,
	const std::string&    the_name  ) ;
      // ======================================================================			
      template <class CLASS>
      inline const std::string&
      setFullName
      ( const CLASS&          the_class , 
	const std::string&    the_name  )
      { return setFullName ( typeid ( the_class ) , the_name ) ;  }
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
    /** @class LogValue
     *  Trivial (scalar) parameter \f$ R \rigtharrow (0,+\infty) \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.c
     *  @date 2026-02-21
     */
    class LogValue
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
    public :
      // ======================================================================
      /// get the value 
      inline double             value      () const { return m_value.value () ; }
      /// get the full parameter name  
      inline const std::string& name       () const { return m_value.name  () ; }      
      // ======================================================================
    public : // conversion 
      // ======================================================================
      /// implicit conversion to double 
      inline operator double  () const { return m_value.value ()  ; } 
      /// set the value from double
      inline LogValue& operator=( const double value )
      { setValue ( value ) ; return *this ; }
      // ======================================================================
    public : // get & set log-value 
      // ======================================================================
      /// get log-value 
      inline double logValue   () const { return m_logValue ; }
      /// set new value for log-parameter 
      bool          setLogValue
      ( const double value         , 
	const bool   force = false ) ; 
      /// set new value for     parameter 
      bool          setValue
      ( const double value         , 
	const bool   force = false ) ; 
      // ======================================================================
    public :
      // ======================================================================
      /** set the full name of parameter
       *  @param the_class the name of owner/holder class
       *  @parm  the_name  the parameter  name
       */
      inline const std::string&
      setFullName
      ( const std::string&    the_class ,
	const std::string&    the_name  )
      { return m_value.setFullName ( the_class , the_name ) ; } 
      // ======================================================================			
      /** set the full name of parameter
       *  @param the_class the type-info of owner/holder class
       *  @parm  the_name  the parameter  name
       */
      inline const std::string&
      setFullName
      ( const std::type_info& the_class ,
	const std::string&    the_name  ) 
      { return m_value.setFullName ( the_class , the_name ) ; } 
      // ======================================================================			
      template <class CLASS>
      inline const std::string&
      setFullName
      ( const CLASS&          the_class , 
	const std::string&    the_name  )
      { return this -> setFullName ( typeid ( the_class ) , the_name ) ;  }
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
      /// the value 
      Value  m_value    { 1 } ; // the value 
      // ======================================================================
    } ;
    // ========================================================================
    /** @class InRange
     *  Trivial parameter with limit \f$ R \rightarrow [A,B] \f$
     *
     *  The actual transformation is 
     *  \f[ p = ( B - A ) * \sin^2 \frac{\pi x}{2} + A \f]
     *
     *  such that
     *   -  \f$ p(0)= A \f$ 
     *   -  \f$ p(1)= B \f$
     *   -  \f$   \min {A,B} \le p \le \max { A, B} 
     */
    class InRange
    {
    public :
      // =====================================================================
      /** full constructor
       *  @param value parameter value
       *  @param avalue A-value 
       *  @param bvalue Bvalue 
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
      /// get the value 
      inline double             value      () const { return m_value.value () ; }
      /// get the full parameter name  
      inline const std::string& name       () const { return m_value.name  () ; }
      // ======================================================================
    public : // conversion 
      // ======================================================================
      /// implicit conversion to double 
      inline operator double  () const { return m_value.value ()  ; }
      /// set the value from double
      inline InRange& operator=( const double value )
      { setValue ( value ) ; return *this ; }
      // ======================================================================
    public : // get & set log-value 
      // ======================================================================
      /// get external value  
      inline double external    () const { return m_external  ; }
      /// set new value for external -parameter 
      bool          setExternal
      ( const double value         ,
	const bool   force = false ) ;
      /// set new value for     parameter 
      bool          setValue
      ( const double value         ,
	const bool   force = false ) ;
      // ======================================================================
    public :
      // ======================================================================
      /** set the full name of parameter
       *  @param the_class the name of owner/holder class
       *  @parm  the_name  the parameter  name
       */
      inline const std::string&
      setFullName
      ( const std::string&    the_class ,
	const std::string&    the_name  )
      { return m_value.setFullName ( the_class , the_name ) ; } 
      // ======================================================================			
      /** set the full name of parameter
       *  @param the_class the type-info of owner/holder class
       *  @parm  the_name  the parameter  name
       */
      inline const std::string&
      setFullName
      ( const std::type_info& the_class ,
	const std::string&    the_name  ) 
      { return m_value.setFullName ( the_class , the_name ) ; } 
      // ======================================================================			
      template <class CLASS>
      inline const  std::string&
      setFullName
      ( const CLASS&          the_class , 
	const std::string&    the_name  )
      { return this->setFullName ( typeid ( the_class ) , the_name ) ;  }
      // ======================================================================			      
    public: 
      // ======================================================================
      /// minimal value 
      inline double vmin () const  { return std::min ( m_A , m_B ) ; } 
      /// maximal  value 
      inline double vmax () const  { return std::max ( m_A , m_B ) ; } 
      // ======================================================================      
    public : // variable transformation 
      // =====================================================================
      /// external -> internal 
      double t ( const double x ) const ;
      /// interbnal-> external 
      double x ( const double t ) const ;
      // =====================================================================
    public: 
      // ======================================================================
      /// unique tag
      std::size_t tag () const ; 
      // ======================================================================
    private:
      // =====================================================================
      /// A-limit
      double m_A         { 0 } ; // A-limit
      /// B-limit
      double m_B         { 1 } ; // B-limit
      /// external value 
      double m_external  { 0 }  ; // external value 
      /// the value ;
      Value m_value             ; // the value
      // =====================================================================
    } ;
    // ========================================================================
    /** @class Scale 
     *  Trivial scalar parameter scale != 0 
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.c
     *  @date 2026-02-21
     */
    class Scale
    {
    public :
      // ======================================================================
      /** full constructor
       *  @param value the value of scale parameter 
       *  @param name  the name of scale parameer 
       *  @param the_class name of the (owner/holder) class 
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
      inline double             scale    () const { return m_scale.value  () ; }
      /// get the value 
      inline double             value    () const { return m_scale.value  () ; }
      /// get the full parameter name  
      inline const std::string& name     () const { return m_scale.name   () ; }
      // ======================================================================
    public :
      // ======================================================================
      /// implicit conversion to double 
      inline operator double () const { return m_scale.value () ; } 
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
      // ======================================================================
    public :
      // ======================================================================
      /** set the full name of parameter
       *  @param the_class the name of owner/holder class
       *  @parm  the_name  the parameter  name
       */
      inline const std::string&
      setFullName
      ( const std::string&    the_class ,
	const std::string&    the_name  )
      { return m_scale.setFullName ( the_class , the_name ) ; } 
      // ======================================================================			
      /** set the full name of parameter
       *  @param the_class the type-info of owner/holder class
       *  @parm  the_name  the parameter  name
       */
      inline const std::string&
      setFullName
      ( const std::type_info& the_class ,
	const std::string&    the_name  ) 
      { return m_scale.setFullName ( the_class , the_name ) ; } 
      // ======================================================================			
      template <class CLASS>
      inline const  std::string&
      setFullName
      ( const CLASS&          the_class , 
	const std::string&    the_name  )
      { return this->setFullName ( typeid ( the_class ) , the_name ) ;  }
      // ======================================================================			      
    public: 
      // ======================================================================
      /// unique tag
      std::size_t tag () const ; 
      // ======================================================================			      
    private :
      // ======================================================================
      /// parameter 
      Ostap::Math::Value m_scale    { 1    } ; // parameter
      /// forced to be positive ? 
      bool               m_positive { true } ; // forced to be positive ? 
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
    public : 
      // =====================================================================
      inline const Scale&       scale_var  () const { return m_scale ; }
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
    public :
      // ======================================================================
      /** set the full name of parameter
       *  @param the_class the name of owner/holder class
       *  @parm  scale_name  the scale parameter  name
       *  @parm  shift_name  the hift parameter  name
       */
      void setFullName
      ( const std::string&    the_class   ,
	const std::string&    scale_name  , 
	const std::string&    shift_name  ) ; 
      // ======================================================================			
      /** set the full name of parameter
       *  @param the_class the type-info of owner/holder class
       *  @parm  scale_name  the scale parameter  name
       *  @parm  shift_name  the hift parameter  name
       */
      void setFullName
      ( const std::type_info& the_class  ,
	const std::string&    scale_name , 
	const std::string&    shift_name ) ;
      // ======================================================================			
      template <class CLASS>
      inline void setFullName
      ( const CLASS&          the_class  , 
	const std::string&    scale_name , 
	const std::string&    shift_name ) 
      { this->setFullName ( typeid ( the_class ) , scale_name , shift_name ) ;  }
      // ======================================================================			      
    public: 
      // ======================================================================
      /// x -> t transfrormation
      inline double t ( const double x ) const
      { return ( x - m_shift.value () ) / m_scale.value () ;  } 
      // ======================================================================			      
      /// t -> x transfrormation
      inline double x ( const double t ) const
      { return   t * m_scale.value ()   + m_shift.value () ;  } 
      // ======================================================================			      
    public: 
      // ======================================================================
      /// unique tag
      std::size_t tag () const ; 
      // ======================================================================			      
    private :
      // =====================================================================
      /// scale parameter 
      Scale m_scale { 1 } ; // scale parameter 
      /// shift parameter 
      Value m_shift { 0 } ; // shift parameter 
      // =====================================================================
    } ;
    // ========================================================================
    /** @typedef ScaleAndShift
     *  anothe name for ShiftandScale
     *  @see ShiftAndScale
     */
    typedef ShiftAndScale ScaleAndShift ; 
    // ========================================================================
    /** @class AB
     *  Two positive variables A & B
     *  - variables are positive and defined via their logarithms
     */
    class AB
    {
      // ======================================================================
    public :
      // ======================================================================
      /** @param a a-parameter 
       *  @param b b-parameter 
       *  @param aname the name of a-parameter
       *  @param bname the name of b-parameter
       *  @param the_class name of the (owner/holder) class 
       */
      AB
      ( const double       a        =  1  ,
	const double       b        =  1  ,
	const std::string& aname    = "a" ,
	const std::string& bname    = "b" ,
	const std::string& the_class = ""  ) ;
      // ====================================================================== 
      /** @param a a-parameter 
       *  @param b b-parameter 
       *  @param aname the name of a-parameter
       *  @param bname the name of b-parameter
       *  @param the_class name of the (owner/holder) class 
       */
      AB
      ( const double          a         ,
	const double          b         ,
	const std::string&    aname     ,
	const std::string&    bname     ,
	const std::type_info& the_class ) ;
      // ======================================================================
      /** @param a a-parameter 
       *  @param b b-parameter 
       *  @param aname the name of a-parameter
       *  @param bname the name of b-parameter
       *  @param the_class name of the (owner/holder) class 
       */
      template <class CLASS>
      AB
      ( const double          a         ,
	const double          b         ,
	const std::string&    aname     ,
	const std::string&    bname     ,
	const CLASS&          the_class )
	: AB ( a , b , aname , bname , typeid ( the_class ) ) 
      {}      
      // =====================================================================
    public :
      // =====================================================================
      // get a 
      inline double             a     () const { return m_a.value    () ; }
      // get b 
      inline double             b     () const { return m_b.value    () ; }
      /// get a-name  
      inline const std::string& aname () const { return m_a.name     () ; }
      /// get b-name  
      inline const std::string& bname () const { return m_b.name     () ; }
      /// get log-a
      inline double             logA  () const { return m_a.logValue () ; } 
      /// get log-b
      inline double             logB  () const { return m_b.logValue () ; } 
      // =====================================================================
    public : 
      // =====================================================================
      /// a-parameter 
      inline const LogValue&    avar  () const { return m_a ; }
      /// b-parameter 
      inline const LogValue&    bvar  () const { return m_b ; }
      // =====================================================================      
    public: // setters 
      // =====================================================================
      inline bool setLogA ( const double value , const bool force = false )
      { return m_a.setLogValue ( value , force ) ; } 
      inline bool setLogB ( const double value , const bool force = false )
      { return m_b.setLogValue ( value , force ) ; } 
      inline bool setA    ( const double value , const bool force = false )
      { return m_a.setValue    ( value , force ) ; } 
      inline bool setB    ( const double value , const bool force = false )
      { return m_b.setValue    ( value , force ) ; }
      /// set two parameters at once
      inline bool setLogAB
      ( const double avalue ,
	const double bvalue ,
	const bool   force  = false )
      {
	const bool ma = m_a.setLogValue ( avalue , force ) ;
	const bool mb = m_b.setLogValue ( bvalue , force ) ;
	return ma && mb ;
      }
      /// set two parameters at once
      inline bool setAB
      ( const double avalue ,
	const double bvalue ,
	const bool   force  = false )
      {
	const bool ma = m_a.setValue ( avalue , force ) ;
	const bool mb = m_b.setValue ( bvalue , force ) ;
	return ma && mb ;
      }	
      // =====================================================================      
    public :
      // ======================================================================
      /** set the full name of parameter
       *  @param the_class the name of owner/holder class
       *  @parm  pname  the name of a-parameter
       *  @parm  qname  the name of b-parameter
       */
      void setFullName
      ( const std::string& the_class ,
	const std::string& aname     , 
	const std::string& bname     ) ;
      // ======================================================================			
      /** set the full name of parameter
       *  @param the_class the type-info of owner/holder class
       *  @parm  pname  the name of a-parameter
       *  @parm  qname  the name of b-parameter
       */
      void setFullName
      ( const std::type_info& the_class  ,
	const std::string&    aname     , 
	const std::string&    bname     ) ;
      // ======================================================================			
      /** set the full name of parameter
       *  @param the_class the type-info of owner/holder class
       *  @parm  pname  the name of a-parameter
       *  @parm  qname  the name of b-parameter
       */
      template <class CLASS>
      inline void setFullName
      ( const CLASS&          the_class  , 
	const std::string&    aname , 
	const std::string&    bname ) 
      { this->setFullName ( typeid ( the_class ) , aname , bname ) ;  }
      // ======================================================================			      
    public: // helper  expressions 
      // ======================================================================
      /// helper expression \f$ \log \Beta (a , b ) \f$
      double log_Beta_ab () const ; 
      /// helper expression \f$ \frac{1}{\Beta (a , b ) {\f$
      double inv_Beta_ab () const ; 
      // ======================================================================			      
    public: 
      // ======================================================================
      /// unique tag
      std::size_t tag () const ; 
      // ======================================================================			      
    private :
      // ======================================================================
      /// parameter A
      LogValue m_a { 0 } ; // parameter A 
      /// parameter B 
      LogValue m_b { 0 } ; // parameter B
      // ======================================================================      
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
      bool setLogP ( const double value , const bool force = false ) ; 
      bool setLogQ ( const double value , const bool force = false ) ; 
      bool setP    ( const double value , const bool force = false ) ;
      bool setQ    ( const double value , const bool force = false ) ; 
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
      // =====================================================================      
    public :
      // ======================================================================
      /** set the full name of parameter
       *  @param the_class the name of owner/holder class
       *  @parm  pname  the name of p-parameter
       *  @parm  qname  the name of q-parameter
       */
      void setFullName
      ( const std::string& the_class ,
	const std::string& pname     , 
	const std::string& qname     ) ;
      // ======================================================================			
      /** set the full name of parameter
       *  @param the_class the type-info of owner/holder class
       *  @parm  pname  the name of p-parameter
       *  @parm  qname  the name of q-parameter
       */
      void setFullName
      ( const std::type_info& the_class  ,
	const std::string&    pname     , 
	const std::string&    qname     ) ;
      // ======================================================================			
      /** set the full name of parameter
       *  @param the_class the type-info of owner/holder class
       *  @parm  pname  the name of p-parameter
       *  @parm  qname  the name of q-parameter
       */
      template <class CLASS>
      inline void setFullName
      ( const CLASS&          the_class  , 
	const std::string&    pname , 
	const std::string&    qname ) 
      { this->setFullName ( typeid ( the_class ) , qname , qname ) ;  }
      // ======================================================================			      
    public: 
      // ======================================================================
      /// get cached value of ln B(p,q)
      inline double log_Beta_pq () const { return m_log_Beta_pq ; } 
      /// get cached value of 1/B(p,q)
      inline double inv_Beta_pq () const { return m_inv_Beta_pq ; } 
      // ======================================================================			      
    public: 
      // ======================================================================
      /// unique tag
      std::size_t tag () const ; 
      // ======================================================================			      
    private :
      // ======================================================================
      /// parameter P 
      LogValue m_p         { 1 } ; // parameter P 
      /// parameter Q 
      LogValue m_q         { 1 } ; // parameter Q
      /// cached value of \f$ \log \Beta (p,q \f$ 
      double m_log_Beta_pq { 0 } ; // log Beta (p,q) 
      /// cached value of \f$ \frac{1} { \Beta (p,q } \f$ 
      double m_inv_Beta_pq { 0 } ; // inv Beta (p,q) 
      // ======================================================================      
    } ;    
    // ========================================================================


    // ========================================================================
    // "Vector" parameters : polynomials, etc..
    // ========================================================================
    /** @class Parameters Ostap/Parameters.h
     *  Holder for parameters 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-02-10
     */
    class Parameters 
    {
    public:
      // =======================================================================
      /// actual type of parameters 
      typedef std::vector<double>        PARAMETERS ;
      typedef PARAMETERS::value_type     data_type  ;
      /// iterator type 
      typedef PARAMETERS::const_iterator iterator   ;
      // =======================================================================
    public:
      // =======================================================================
      /// constructor from number of parameters 
      Parameters ( const std::size_t           np = 1 ) ; 
      /// constructor from  the list of parameters 
      Parameters ( const PARAMETERS&  pars   ) ;
      /// constructor from  the list of parameters 
      Parameters (       PARAMETERS&& pars   ) ;
      /// templated constructor from the sequence of parameters 
      template <typename ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,data_type>::value,bool>::type = true >           
      Parameters
      ( ITERATOR begin , 
        ITERATOR end   )
        : m_pars ( begin , end )
      {}
      // ======================================================================
    public:
      // ======================================================================
      /// number of parameters 
      inline std::size_t npars  () const { return m_pars.size()     ; }
      // ======================================================================
      /// get the parameter value
      inline double  par          ( const std::size_t k ) const
      { return ( k < m_pars.size() ) ? m_pars [ k ] : 0.0 ; }
      // ======================================================================
      /// get all parameters:
      const PARAMETERS& pars () const { return m_pars ; }
      // ======================================================================
      /** set k-parameter
       *  @param k index
       *  @param value new value 
       *  @return true if parameter is actually changed 
       */
      inline bool setPar          
      ( const std::size_t  k             , 
        const double       value         ,
        const bool         force = false ) 
      { return k < m_pars.size() ?  _setPar ( k , value , force ) : false ; }
      // ======================================================================
      /** set several/all parameters at once 
       *  @param begin  start iterator for the sequence of coefficients 
       *  @param end    end   iterator for the sequence of coefficients 
       *  @return true if at least one parameter is actually changed 
       */
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type ,                
                typename std::enable_if<std::is_convertible<value_type,long double>::value,bool>::type = true >           
      inline bool setPars 
      ( ITERATOR   begin         , 
        ITERATOR   end           ,
        const bool force = false ) 
      {
        bool update = false ;
        const unsigned int   N = m_pars.size()  ;
        for ( unsigned short k ; k < N && begin != end ;  ++k, ++begin ) 
          { update = this->_setPar ( k  , *begin , force ) ? true : update ; }
        return update ;
      }
      // ======================================================================
      /** set several/all parameters at once 
       *  @param pars (NIPUT) vector of parameters 
       *  @return true if at least one parameter is actually changed 
       */
      inline bool setPars
      ( const PARAMETERS& pars          ,
        const bool        force = false ) 
      { return setPars ( pars.begin() , pars.end() , force ) ; }
      // ======================================================================
      /// all parameters are zero ?
      bool           zero   () const ;
      // ======================================================================
    public: // Reset 
      // ======================================================================
      /// reset all parameters to zero 
      void reset() ;
      /// reset all parameters to zero 
      void Reset() { reset () ; }
      // ======================================================================
    public:
      // ======================================================================
      /** Filter out very small terms/parameters 
       *  the term is considered to be very small if
       *   - it is numerically zero: \f$ c_k \approx 0 \f% 
       *   - or if epsilon  > 0: \f$ \left| c_k \right| \le \epsilon \f$ 
       *   - or if scale   != 0: \f$ \left| s \right| + \left| c_k \right| \approx \left| s \right| \f$ 
       *  @param  epsilon  parameter to define "smalness" of terms
       *  @param  scale    parameter to define "smalness" of terms
       *  @return number of nullified terms
       */
      std::size_t remove_noise 
      ( const double epsilon = 0 , 
        const double scale   = 0 ) ;
      // ======================================================================
    public: // simple  manipulations with parameters 
      // ======================================================================
      /// simple  manipulations with parameters: scale it! 
      /// Parameters& operator *= ( const double a ) ;     // scale it! 
      /// simple  manipulations with parameters: scale it! 
      /// Parameters& operator /= ( const double a ) ;     // scale it! 
      // ======================================================================
    public: // expose constant iterators 
      // ======================================================================
      /// begin iterator
      iterator begin () const { return m_pars.begin () ; }
      /// end   iterator
      iterator end   () const { return m_pars.end   () ; }
      // ======================================================================
    protected:
      // ======================================================================
      /// swap two parameter sets 
      void swap ( Parameters& right ) ;
      // ======================================================================
    private:
      // ======================================================================
      /** set k-parameter
       *  @param k index
       *  @param value new value 
       *  @return true if parameter is actually changed 
       */
    bool _setPar 
    ( const std::size_t k             , 
      const double      value         ,
      const bool        force = false ) ;
      // ======================================================================
    public : // few helper static fuctions 
      // ======================================================================
      /// join two vectors togather 
      static std::vector<double>
      join
      ( const std::vector<double>& a ,
	const std::vector<double>& b ) ;
      /// join scalar & vector together 
      static std::vector<double>
      join
      ( const double               a ,
	const std::vector<double>& b ) ;
      /// join 2 scalars & vector together 
      static std::vector<double>
      join
      ( const double               a1 ,
	const double               a2 ,
	const std::vector<double>& b  ) ;
      /// join 3 scalars & vector together 
      static std::vector<double>
      join
      ( const double               a1 ,
	const double               a2 ,
	const double               a3 ,
	const std::vector<double>& b  ) ;
      /// join 4 scalars & vector together 
      static std::vector<double>
      join
      ( const double               a1 ,
	const double               a2 ,
	const double               a3 ,
	const double               a4 ,
	const std::vector<double>& b  ) ;
      /// join vector & scalar together 
      static std::vector<double>
      join
      ( const std::vector<double>& a , 
	const double               b ) ;
      /// join vector & 2 scalars together 
      static std::vector<double>
      join
      ( const std::vector<double>& a  , 
	const double               b1 , 
	const double               b2 ) ;
      /// join vector & 3 scalars together 
      static std::vector<double>
      join
      ( const std::vector<double>& a  , 
	const double               b1 , 
	const double               b2 , 
	const double               b3 ) ;
      /// join vector & 4 scalars together 
      static std::vector<double>
      join
      ( const std::vector<double>& a  , 
	const double               b1 , 
	const double               b2 , 
	const double               b3 , 
	const double               b4 ) ;
      // ======================================================================
    protected :
      // ======================================================================
      /// parameters 
      PARAMETERS m_pars ; //  vector of parameters 
      // ======================================================================
    } ;
    // ========================================================================
  } //                                             end of namespace Ostap::Math 
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_PARAMETERS_H
// ============================================================================
