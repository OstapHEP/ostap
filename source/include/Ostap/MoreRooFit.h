// ============================================================================
#ifndef OSTAP_MOREROOFIT_H 
#define OSTAP_MOREROOFIT_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <functional>
// ============================================================================
// ROOT/RooFit 
// ============================================================================
#include "RVersion.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooRealConstant.h"
#include "RooRealProxy.h"
#include "RooAbsPdf.h"
#include "RooGlobalFunc.h"
// ============================================================================
namespace Ostap 
{
  // ==========================================================================
  /** @namespace Ostap::MoreRooFit   Ostap/MoreRooFit.h
   *  Collection of small additios to RooFit 
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
   *  @date 2019-11-21
   */
  namespace MoreRooFit 
  {
    // ========================================================================
    /** @class Addition
     *  A simple modification of class RooAddition
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */ 
    class Addition : public RooAddition
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::Addition , 1 ) ;  // sum of RooAbsReal objects
      // ========================================================================
    public:
      // ========================================================================
      /// construct x + y 
      Addition
      ( const std::string& name  ,  
        const std::string& title ,
        RooAbsReal&        x     ,
        RooAbsReal&        y     ) ;
      // ======================================================================
      /// constructor with two variables 
      Addition
      ( RooAbsReal&        x           , 
        RooAbsReal&        y           ,
        const std::string& name  = ""  , 
        const std::string& title = ""  ) 
        : Addition ( name , title , x , y )
      {}
      // ======================================================================
      /// constructor with two variables 
      Addition
      ( RooAbsReal&        x           , 
        const double       y           ,
        const std::string& name  = ""  , 
        const std::string& title = ""  ) 
        : Addition ( name , title , x , RooFit::RooConst ( y ) )
      {}
      /// constructor with two variables 
      Addition ( const double       x           ,
                 RooAbsReal&        y           , 
                 const std::string& name  = ""  , 
                 const std::string& title = ""  ) 
        : Addition ( name , title , RooFit::RooConst ( x ) , y )
      {}
      /// copy 
      Addition ( const Addition&    right       , 
                 const char*        newname = 0 ) ;
      /// destructor
      virtual ~Addition () ;
      /// clone 
      Addition* clone ( const char* newname ) const override ;
      /// fake defautl constructor (needed for serisalization) 
      Addition () = default ;
      // ======================================================================
    public:
      // ======================================================================
      inline const RooAbsReal& x  () const { return m_x .arg() ; }
      inline const RooAbsReal& y  () const { return m_y .arg() ; }    
      // ======================================================================
    private:
      // ======================================================================
      /// variable
      RooRealProxy m_x  {} ; // variable
      RooRealProxy m_y  {} ; // variable
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Addition2
     *  A simple modification of class RooAddition
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */ 
    class Addition2 : public RooAddition
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::Addition2 , 1 ) ;  // sum of RooAbsReal objects
      // ========================================================================
    public:
      // ========================================================================
      /// construct c1*a + c2*b 
      Addition2 
      ( const std::string& name  ,  
        const std::string& title ,
        RooAbsReal&        a     ,
        RooAbsReal&        b     ,
        RooAbsReal&        c1    ,
        RooAbsReal&        c2    ) ;
      // ======================================================================
      /// construct c1*a + c2*b 
      Addition2
      ( const std::string& name   ,   
        const std::string& title  ,
        RooAbsReal&        a      ,
        RooAbsReal&        b      ,
        const double       c1 = 1 , 
        const double       c2 = 1 )
        : Addition2 ( name , title , a , b, 
                      RooFit::RooConst ( c1 ) , 
                      RooFit::RooConst ( c2 ) )
      {}               
      /// copy 
      Addition2 ( const Addition2&    right       , 
                  const char*        newname = 0 ) ;
      /// destructor
      virtual ~Addition2 () ;
      /// clone 
      Addition2* clone ( const char* newname ) const override ;
      /// fake defautl constructor (needed for serisalization) 
      Addition2 () = default ;
      // ======================================================================
    public:
      // ======================================================================
      inline const RooAbsReal& x  () const { return m_x  .arg() ; }
      inline const RooAbsReal& y  () const { return m_y  .arg() ; }    
      inline const RooAbsReal& c1 () const { return m_c1 .arg() ; }
      inline const RooAbsReal& c2 () const { return m_c2 .arg() ; }    
      // ======================================================================
    private:
      // ======================================================================
      /// variable
      RooRealProxy m_x  {} ; // variable
      RooRealProxy m_y  {} ; // variable
      RooRealProxy m_c1 {} ; // variable
      RooRealProxy m_c2 {} ; // variable
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Product
     *  A simple modification of class RooProduct
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */ 
    class Product : public RooProduct
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::Product , 1 ) ;  // sum of RooAbsReal objects
      // ========================================================================
    public:
      // ========================================================================
      /// constructor with two variables 
      Product
      ( const std::string& name  , 
        const std::string& title , 
        RooAbsReal&        a     , 
        RooAbsReal&        b     ) ;
      /// constructor with two variables 
      Product
      ( RooAbsReal&        a           , 
        RooAbsReal&        b           ,
        const std::string& name  = ""  , 
        const std::string& title = ""  ) 
        : Product ( name , title , a , b )
      {}
      /// constructor with two variables 
      Product
      ( RooAbsReal&        a           , 
        const double       b           ,
        const std::string& name  = ""  , 
        const std::string& title = ""  ) 
        : Product ( name , title , a , RooFit::RooConst ( b )  )
      {}
      /// constructor with two variables 
      Product 
      ( const double       a           ,
        RooAbsReal&        b           , 
        const std::string& name  = ""  , 
        const std::string& title = ""  ) 
        : Product ( name , title , RooFit::RooConst ( a ) , b )
      {}
      /// copy 
      Product
      ( const Product&    right       , 
        const char*        newname = 0 ) ;
      /// destructor 
      virtual ~Product () ;
      /// clone 
      Product* clone ( const char* newname ) const override ;
      /// fake defautl constructor (needed for serisalization) 
      Product () = default ;
      // ======================================================================
    public:
      // ======================================================================
      inline const RooAbsReal& x  () const { return m_x .arg() ; }
      inline const RooAbsReal& y  () const { return m_y .arg() ; }    
      // ======================================================================
    private: 
      // ======================================================================
      /// variable
      RooRealProxy m_x  {} ; // variable
      RooRealProxy m_y  {} ; // variable
      // ======================================================================
    }; // 
    // ========================================================================
    /** @class Subtraction 
     *  A simple modification of class RooAddition
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */ 
    class Subtraction : public Addition2
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::Subtraction , 1 ) ;  // Difference of RooAbsReal objects
      // ========================================================================
    public:
      // ========================================================================
      /// constructor with list of variables 
      Subtraction
      ( const std::string& name  , 
        const std::string& title , 
        RooAbsReal&        a     , 
        RooAbsReal&        b     ) ;
      /// constructor with two variables 
      Subtraction
      ( RooAbsReal&        a           , 
        RooAbsReal&        b           ,
        const std::string& name  = ""  , 
        const std::string& title = ""  ) 
        : Subtraction ( name , title , a , b )
      {}
      /// constructor with two variables 
      Subtraction 
      ( RooAbsReal&        a           , 
        const double       b           ,
        const std::string& name  = ""  , 
        const std::string& title = ""  ) 
        : Subtraction ( name , title , a , RooFit::RooConst ( b ) )
      {}
      /// constructor with two variables 
      Subtraction 
      ( const double       a           ,
        RooAbsReal&        b           , 
        const std::string& name  = ""  , 
        const std::string& title = ""  ) 
        : Subtraction ( name , title , RooFit::RooConst ( a ) , b )
      {}
      /// copy 
      Subtraction
      ( const Subtraction& right       , 
        const char*        newname = 0 ) ;
      /// destructor
      virtual ~Subtraction() ;
      /// clone 
      Subtraction* clone ( const char* newname ) const override ;
      /// fake default constructor (needed for serisalization) 
      Subtraction () = default ;
      // ======================================================================
    }; // 
    // ========================================================================
    /** @class Constant 
     *  simple extension for the constant varibale 
     *  - add (fuictive) dependencies
     *  @see RooConstVar 
     */
    class Constant : public RooAbsReal
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::Constant , 1 ) ;  
      // ========================================================================
    public:
      // =======================================================================
      /// constructor with no variables 
      Constant 
      ( const std::string& name      , 
        const std::string& title     , 
        const double       value     ) ;
      /// constructor with one variable 
      Constant
      ( const std::string& name      , 
        const std::string& title     , 
        const double       value     , 
        RooAbsReal&        x         ) ;
      /// constructor with two variables 
      Constant
      ( const std::string& name      , 
        const std::string& title     , 
        const double       value     , 
        RooAbsReal&        x         ,
        RooAbsReal&        y         ) ;
      /// constructor with three variables 
      Constant
      ( const std::string& name      , 
        const std::string& title     , 
        const double       value     ,
        RooAbsReal&        x         ,
        RooAbsReal&        y         ,
        RooAbsReal&        z         ) ;
      /// constructor with list of variables 
      Constant
      ( const std::string& name      , 
        const std::string& title     , 
        const double       value     , 
        const RooArgList&  v         ) ;
      /// copy 
      Constant
      ( const Constant&    right       , 
        const char*        newname = 0 ) ;
      /// clone 
      Constant* clone ( const char* newname ) const override ;
      /// fake default constructor (needed for serisalization)
      Constant () = default ;
      // ======================================================================
    public: // redefined methods 
      // ======================================================================
      /// the actual evaluation of the result 
      Double_t evaluate () const override ;
      // ======================================================================
      /// the actual evaluation of the result
      double   getValV  ( const RooArgSet* a = nullptr ) const override ;
      // ======================================================================
      // srite to the stream 
      void writeToStream(std::ostream& os, bool compact) const override ;
      // ======================================================================
#if ROOT_VERSION(6,24,0)<=ROOT_VERSION_CODE
      // ======================================================================
      RooSpan<const double> 
      getValues ( RooBatchCompute::RunContext& evalData , 
                  const RooArgSet*             aset     ) const override ;
      // ======================================================================
#endif 
      // ======================================================================
   public:
      // ======================================================================
      /// variables 
      inline const RooArgList& vlst () const { return m_vlst ; }
      // ======================================================================
    private:
      // ======================================================================
      /// constant 
      double       m_value { 0 } ; // constant 
      /// variables 
      RooListProxy m_vlst  {}    ; // variables 
      // ======================================================================
    } ;  
    // ========================================================================
    /** @class OneVar
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2020-03-17
     */ 
    class OneVar : public RooAbsReal 
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::OneVar , 1 ) ;  
      // ========================================================================
    public:
      // ========================================================================
      /// constructor with a variable
      OneVar
      ( const std::string& name  , 
        const std::string& title , 
        RooAbsReal&        x     ) ;
      /// copy 
      OneVar
      ( const OneVar&   right       , 
        const char* newname = 0 ) ;
      /// destructor 
      virtual ~OneVar () ;
      /// fake default constructor (needed for serisalization)
      OneVar () = default ;
      // ======================================================================
    public:
      // ======================================================================
      inline const RooAbsReal& x () const { return m_x.arg() ; }
      // ======================================================================
    protected :
      // ======================================================================
      /// variable
      RooRealProxy m_x {} ; // variable
      // ======================================================================
    }; //    
    // ========================================================================
    /** @class FunOneVar 
     *  Class for tranformation of variables 
     */
    class FunOneVar : public OneVar  
    {
      // ======================================================================
      ClassDefOverride ( Ostap::MoreRooFit::FunOneVar , 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /** create the tranformation variable
       *  \f$ z = f(x) \f$ 
       *  @param fun the transformation function 
       *  @param x the first variable 
       *  @param name the name of variable 
       *  @param title the title  of variable 
       */
      template <class FUNCTION> 
      FunOneVar ( const std::string& name  , 
                  const std::string& title , 
                  FUNCTION           fun   , 
                  RooAbsReal&        x     )
        : OneVar ( name , title , x )
        , m_fun  ( fun ) 
      {}
      // =====================================================================
      /** create the tranformation variable
       *  \f$ z = f(x) \f$ 
       *  @param fun the transformation function 
       *  @param x the first variable 
       *  @param name the name of variable 
       *  @param title the title  of variable 
       */
      template <class FUNCTION> 
      FunOneVar ( FUNCTION    fun               , 
                  RooAbsReal& x                 ,                   
                  const std::string& name  = "" , 
                  const std::string& title = "" )
        : FunOneVar ( name ,   title , fun , x )
      {}
      /// copy  
      FunOneVar ( const FunOneVar& right        , 
                  const char*      newname  = 0 ) ;
      // ======================================================================
      /// destructor 
      virtual  ~FunOneVar () ;
      // ======================================================================
      /// clone 
      FunOneVar* clone ( const char* newname ) const override ;
      // ======================================================================
    public :  
      // ======================================================================     
      /** create the tranformation variable
       *  \f$ z = f(x) \f$ 
       *  @param fun the transformation function 
       *  @param x the first variable 
       *  @param name the name of variable 
       *  @param title the title  of variable 
       */
      template <class FUNCTION>
      static inline FunOneVar 
      create ( FUNCTION    fun               , 
               RooAbsReal& x                 ,
               const std::string& name  = "" , 
               const std::string& title = "" )
      { return FunOneVar ( fun , x , name , title ) ; }
      // ======================================================================
    public:
      // ======================================================================
      // fictive default constructor 
      FunOneVar () {}      
      // ======================================================================
    public :  
      // ======================================================================     
      // the actual evaluation of the result 
      Double_t evaluate () const override ;
      // ======================================================================
    private:
      // ======================================================================
      /// actual transformation 
      std::function<double(double)> m_fun {} ; // actual transformation 
      // ======================================================================
    } ;  
    // ========================================================================
    /** @class TwoVars
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2020-03-17
     */ 
    class TwoVars : public OneVar
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::TwoVars , 1 ) ;  
      // ========================================================================
    public:
      // ========================================================================
      /// constructor with a variable
      TwoVars 
      ( const std::string& name  , 
        const std::string& title , 
        RooAbsReal&        x     ,
        RooAbsReal&        y     ) ;
      TwoVars
      ( RooAbsReal&        x     ,
        RooAbsReal&        y     ,
        const std::string& name  , 
        const std::string& title ) 
        : TwoVars ( name , title , x , y ) 
      {}
      /// copy 
      TwoVars 
      ( const TwoVars& right       , 
        const char*    newname = 0 ) ;
      /// destructor 
      virtual ~TwoVars () ;
      /// fake default constructor (needed for serisalization)
      TwoVars () = default ;
      // ======================================================================
    public:
      // ======================================================================
      inline const RooAbsReal& y () const { return m_y.arg () ; }
      // ======================================================================
    protected :
      // ======================================================================
      /// variable
      RooRealProxy m_y {} ; // variable
      // ======================================================================
    }; 
    // ========================================================================
    /** @class TwoVarsFun
     *  Helper class to create function \f$ f(a,b) = f(a*b)\f$
     */
    class FunTwoVars : public TwoVars 
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::FunTwoVars , 1 ) ;  
      // ========================================================================
    public: 
      // ======================================================================
      /** create the tranformation variable as function from two variables 
       *  \f$ z= f(x,y) \f$ 
       *  @param name the name of variable 
       *  @param title the title  of variable 
       *  @param fun the transformation function 
       *  @param x the first variable 
       *  @param y the second variable 
       */
      template <class FUNCTION>
      FunTwoVars ( const std::string& name  , 
                   const std::string& title , 
                   FUNCTION           fun   , 
                   RooAbsReal&        x     , 
                   RooAbsReal&        y     ) 
        : TwoVars ( name , title , x , y ) 
        , m_fun2  ( fun )
      {}  
      // ======================================================================
      /** create the tranformation variable as function from two variables 
       *  \f$ z= f(x,y) \f$ 
       *  @param fun the transformation function 
       *  @param x the first variable 
       *  @param y the second variable 
       *  @param name the name of variable 
       *  @param title the title  of variable 
       */
      template <class FUNCTION>
      FunTwoVars ( FUNCTION fun                    , 
                   RooAbsReal&         x           , 
                   RooAbsReal&         y           ,
                   const std::string&  name  = ""  , 
                   const std::string&  title = ""  ) 
        : FunTwoVars ( name , title , fun , x , y )
      {}
      // ======================================================================
      /** create the tranformation variable as function from two variables 
       *  \f$ z= f(x,y) \f$ 
       *  @param fun the transformation function 
       *  @param x the first variable 
       *  @param y the second variable 
       *  @param name the name of variable 
       *  @param title the title  of variable 
       */
      template <class FUNCTION>
      FunTwoVars ( FUNCTION fun                    , 
                   const double        x           , 
                   RooAbsReal&         y           ,
                   const std::string&  name  = ""  , 
                   const std::string&  title = ""  ) 
        : FunTwoVars ( name , title , fun , RooFit::RooConst ( x ) , y )
      {}
      // ======================================================================
      /** create the tranformation variable as function from two variables 
       *  \f$ z= f(x,y) \f$ 
       *  @param fun the transformation function 
       *  @param x the first variable 
       *  @param y the second variable 
       *  @param name the name of variable 
       *  @param title the title  of variable 
       */
      template <class FUNCTION>
      FunTwoVars ( FUNCTION fun                    , 
                   RooAbsReal&         a           ,
                   const double        y           , 
                   const std::string&  name  = ""  , 
                   const std::string&  title = ""  ) 
        : FunTwoVars ( name , title , fun , x , RooFit::RooConst ( y ) )
      {}
      // ======================================================================
      /// copy 
      FunTwoVars ( const FunTwoVars& right       , 
                   const char*       newname = 0 ) ;
      /// destructor 
      virtual ~FunTwoVars () ;
      /// clone 
      FunTwoVars* clone ( const char* newname ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /** create the tranformation variable as function from two variables 
       *  \f$ z = f(x,y) \f$ 
       *  @param fun the transformation function 
       *  @param x the first variable 
       *  @param y the second variable 
       *  @param name the name of variable 
       *  @param title the title  of variable 
       */
      template <class FUNCTION>
      static inline FunTwoVars 
      create ( FUNCTION fun                    , 
               RooAbsReal&         x           , 
               RooAbsReal&         y           ,
               const std::string&  name  = ""  , 
               const std::string&  title = ""  ) 
      { return FunTwoVars ( fun , x , y , name , title ) ; }
      // ======================================================================
    public:
      // ======================================================================
      // fictive default constructor 
      FunTwoVars() {}      
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;
      // ======================================================================
    private:
      // =====================================================================
      /// the function 
      std::function<double(double,double)> m_fun2 {} ; // the function 
      // =====================================================================      
    } ;
    // ========================================================================
    /** @class Division
     *  Evaluate \f$ \frac{a}{b}\f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Division final : public TwoVars 
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::Division , 2 ) ;  // ratio of RooAbsReal objects
      // ========================================================================
    public:
      // ======================================================================
      Division 
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        x     , 
          RooAbsReal&        y     ) ;
      /// constructor with two variables 
      Division 
        ( RooAbsReal&         x           , 
          RooAbsReal&         y           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Division ( name , title , x , y )
      {}
      /// constructor with two variables 
      Division 
        ( const double        x           , 
          RooAbsReal&         y           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Division ( name , title , RooFit::RooConst ( x ) , y )
      {}
      /// constructor with two variables 
      Division
        ( RooAbsReal&         x           ,
          const double        y           , 
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Division ( name , title , x , RooFit::RooConst ( y ) )
      {}
      // ======================================================================
      /// (fictive) default constructor 
      Division  (){} ; //  =  default ;
      // ======================================================================
      // copy 
      Division 
        ( const Division& right , 
          const char* newname = 0 ) 
        : TwoVars ( right , newname ) 
      {}
      // ======================================================================
      Division* clone ( const char* newname ) const override 
      { return new Division ( *this , newname ) ; }
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Combination
     *  evaluate \f$ f(x,y) = \alpha  x ( \beta + \gamma y )\f$ 
     */
    class Combination final : public TwoVars
    {
      // ======================================================================
      ClassDefOverride(Ostap::MoreRooFit::Combination , 1 ) ;  // combination
      // ======================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Combination 
        ( const std::string& name        , 
          const std::string& title       , 
          RooAbsReal&        x           , 
          RooAbsReal&        y           ,
          const double       alpha = 1   , 
          const double       beta  = 1   , 
          const double       gamma = 1   ) ;
      /// constructor with two variables 
      Combination
        ( RooAbsReal&         x           , 
          RooAbsReal&         y           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  , 
          const double        alpha = 1   , 
          const double        beta  = 1   , 
          const double        gamma = 1   ) 
        : Combination ( name , title , x , y , alpha , beta , gamma )
      {}
      /// constructor with two variables 
      Combination
        ( const double        x           , 
          RooAbsReal&         y           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  , 
          const double        alpha = 1   , 
          const double        beta  = 1   , 
          const double        gamma = 1   )
        : Combination ( name , title , RooFit::RooConst ( x ) , y , alpha , beta , gamma )
      {}
      /// constructor with two variables 
      Combination
        ( RooAbsReal&         x           ,
          const double        y           , 
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ,
          const double        alpha = 1   , 
          const double        beta  = 1   , 
          const double        gamma = 1   )
        : Combination ( name , title , x , RooFit::RooConst ( y ) , alpha , beta , gamma )
      {}
      // ======================================================================
      /// default constructor 
      Combination  () =  default ;
      // ======================================================================
      // copy 
      Combination
        ( const Combination& right , 
          const char* newname = 0 ) 
        : TwoVars ( right , newname ) 
        , m_alpha ( right.m_alpha ) 
        , m_beta  ( right.m_beta  ) 
        , m_gamma ( right.m_gamma ) 
      {}
      // ======================================================================
      Combination* clone ( const char* newname ) const override 
      { return new Combination  ( *this , newname) ; }
      // ======================================================================
    public:
      // ======================================================================
      double alpha () const { return m_alpha ; }
      double beta  () const { return m_beta  ; }
      double gamma () const { return m_gamma ; }      
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;
      // ======================================================================
    private:
      // ======================================================================
      /// alpha 
      double m_alpha  ; // alpha 
      /// beta 
      double m_beta   ; // beta 
      /// gamma 
      double m_gamma  ; // gamma 
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Fraction
     *  Evaluate \f$ \frac{a}{a+b}\f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Fraction final : public TwoVars
    {
      // ======================================================================
      ClassDefOverride(Ostap::MoreRooFit::Fraction , 2 ) ;  // Fraction 
      // ======================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Fraction
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        x     , 
          RooAbsReal&        y     ) ;
      /// constructor with two variables 
      Fraction
        ( RooAbsReal&         x           , 
          RooAbsReal&         y           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Fraction ( name , title , x , y )
      {}
      /// constructor with two variables 
      Fraction
        ( const double        x           , 
          RooAbsReal&         y           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Fraction ( name , title , RooFit::RooConst ( x ) , y )
      {}
      /// constructor with two variables 
      Fraction 
        ( RooAbsReal&         x           ,
          const double        y           , 
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Fraction ( name , title , x , RooFit::RooConst ( y ) )
      {}
      // ======================================================================
      /// default constructor 
      Fraction  () =  default ;
      // ======================================================================
      // copy 
      Fraction ( const Fraction& right , const char* newname = 0 ) 
        : TwoVars ( right , newname ) 
      {}
      // ======================================================================
      Fraction* clone ( const char* newname ) const override 
      { return new Fraction ( *this , newname) ; }
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Asymmetry
     *  Evaluate \f$ c*\frac{a-b}{a+b}\f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Asymmetry final : public TwoVars 
    {
      // ======================================================================
      ClassDefOverride(Ostap::MoreRooFit::Asymmetry , 3 ) ;  // Relative difference 
      // ======================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Asymmetry
        ( const std::string& name         , 
          const std::string& title        , 
          RooAbsReal&        a            , 
          RooAbsReal&        b            , 
          const double       scale  = 1   ) ;
      /// constructor with two variables 
      Asymmetry
        ( RooAbsReal&         a           , 
          RooAbsReal&         b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  , 
          const double        scale = 1   ) 
      : Asymmetry ( name , title , a , b , scale )
      {}
      /// constructor with two variables 
      Asymmetry
        ( RooAbsReal&         a           , 
          const double        b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ,
          const double        scale = 1   ) 
        : Asymmetry ( name , title , a , RooFit::RooConst ( b ) , scale )
      {}
      /// constructor with two variables 
      Asymmetry
        ( const double        a           ,
          RooAbsReal&         b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  , 
          const double        scale = 1   ) 
        : Asymmetry ( name , title , RooFit::RooConst ( a ) , b , scale )
      {}
      // ======================================================================
      /// default constructor 
      Asymmetry () =  default ;
      // ======================================================================
      // copy 
      Asymmetry ( const Asymmetry& right , const char* newname = 0 ) 
        : TwoVars ( right , newname ) 
        , m_scale ( right.m_scale   )
      {}
      // ======================================================================
      Asymmetry* clone ( const char* newname ) const override 
      { return new Asymmetry ( *this , newname ) ; }
      // ======================================================================
      // destructor 
      virtual ~Asymmetry () {};
      // ======================================================================
    public:
      // ======================================================================
      double scale () const { return m_scale ; }
      // ======================================================================
    protected:
      // ======================================================================
      /// the actual evaluation of the result 
      Double_t evaluate () const override ; 
      // ======================================================================
    private:
      // ======================================================================
      /// scaling constant 
      double    m_scale ; // scaling constant 
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Power
     *  Evaluate \f$ a^b \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Power final : public TwoVars 
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::Power , 2 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Power  
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     , 
          RooAbsReal&        b     ) ;
      /// constructor with two variables 
      Power
        ( RooAbsReal&         a           , 
          RooAbsReal&         b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Power ( name , title , a , b )
      {}
      /// constructor with two variables 
      Power 
        ( RooAbsReal&         a           , 
          const double        b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Power ( name , title , a , RooFit::RooConst ( b ) ) 
      {}
      /// constructor with two variables 
      Power 
        ( const double        a           ,
          RooAbsReal&         b           , 
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Power ( name , title , RooFit::RooConst ( a ) , b ) 
      {}
      // ======================================================================
      /// default constructor 
      Power () =  default ;
      // ======================================================================
      // copy 
      Power 
        ( const Power& right , 
          const char* newname = 0 ) 
        : TwoVars ( right , newname ) 
      {}
      // ======================================================================
      Power* clone ( const char* newname ) const override 
      { return new Power ( *this , newname ) ; }
      // ======================================================================
      // destructor 
      virtual ~Power() {} ;
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Abs
     *  Evaluate \f$ abs(ab) \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Abs final : public TwoVars 
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::Abs , 2 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Abs 
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     , 
          RooAbsReal&        b     ) ;
      /// constructor with two variables 
      Abs 
        ( RooAbsReal&         a           , 
          RooAbsReal&         b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Abs ( name , title , a , b )
      {}
      /// constructor with two variables 
      Abs 
        ( RooAbsReal&         a           , 
          const double        b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Abs ( name , title , a , RooFit::RooConst ( b ) )
      {}
      /// constructor with two variables 
      Abs
        ( const double        a           ,
          RooAbsReal&         b           , 
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Abs ( name , title , RooFit::RooConst ( a ) , b )
      {}
      /// constructor with one variable
      Abs
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     ) 
        : Abs ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
      /// fake defautl constructor (needed for serisalization) 
      Abs () = default ;
      // ======================================================================
      // copy 
      Abs ( const Abs& right , const char* newname = 0 ) 
        : TwoVars ( right , newname ) 
      {}
      // ======================================================================
      Abs* clone ( const char* newname ) const override 
      { return new Abs ( *this , newname) ; }
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Exp
     *  Evaluate \f$ exp(ab) \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Exp final : public TwoVars 
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::Exp , 2 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Exp 
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     , 
          RooAbsReal&        b     ) ;
      /// constructor with two variables 
      Exp 
        ( RooAbsReal&         a           , 
          RooAbsReal&         b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Exp ( name , title , a , b )
      {}
      /// constructor with two variables 
      Exp
        ( RooAbsReal&         a           , 
          const double        b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Exp ( name , title , a , RooFit::RooConst ( b )  )
      {}
      /// constructor with two variables 
      Exp
        ( const double        a           ,
          RooAbsReal&         b           , 
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Exp ( name , title , RooFit::RooConst ( a ) , b  )
      {}
      /// constructor with one variable
      Exp 
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     ) 
        : Exp ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
      /// fake defautl constructor (needed for serisalization) 
      Exp () = default ;
      // ======================================================================
      // copy 
      Exp ( const Exp& right , const char* newname = 0 ) 
        : TwoVars ( right , newname ) 
      {}
      // ======================================================================
      Exp* clone ( const char* newname ) const override 
      { return new Exp ( *this , newname) ; }
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Log
     *  Evaluate \f$ \log ab  \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Log final : public TwoVars
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::Log , 2 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Log 
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     , 
          RooAbsReal&        b     ) ;
      /// constructor with two variables
      Log 
        ( RooAbsReal&         a           , 
          RooAbsReal&         b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Log  ( name , title , a , b )
      {}
      /// constructor with two variables
      Log 
        ( RooAbsReal&         a           , 
          const double        b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Log  ( name , title , a , RooFit::RooConst ( b ) )
      {}
      /// constructor with two variables
      Log 
        ( const double        a           ,
          RooAbsReal&         b           , 
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Log  ( name , title , RooFit::RooConst ( a ) , b )
      {}
      /// constructor with one variable
      Log 
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     ) 
        : Log ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
      /// fake defautl constructor (needed for serisalization) 
      Log () = default ;
      // ======================================================================
      // copy 
      Log 
        ( const Log& right , 
          const char* newname = 0 ) 
        : TwoVars ( right , newname ) 
      {}
      // ======================================================================
      Log* clone ( const char* newname ) const override 
      { return new Log ( *this , newname) ; }
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Log10
     *  Evaluate \f$ \log10 ab  \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Log10 final : public TwoVars 
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::Log10 , 2 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Log10
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     , 
          RooAbsReal&        b     ) ;
      /// constructor with two variables
      Log10
        ( RooAbsReal&         a           , 
          RooAbsReal&         b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Log10  ( name , title , a , b )
      {}
      /// constructor with two variables
      Log10
        ( RooAbsReal&         a           , 
          const double        b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Log10  ( name , title , a , RooFit::RooConst ( b ) )
      {}
      /// constructor with two variables
      Log10
        ( const double        a           ,
          RooAbsReal&         b           , 
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Log10  ( name , title , RooFit::RooConst ( a ) , b )
      {}
      /// constructor with one variable
      Log10
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     ) 
        : Log10 ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
      /// fake defautl constructor (needed for serisalization) 
      Log10 () = default ;
      // ======================================================================
      // copy 
      Log10 
        ( const Log10& right , 
          const char* newname = 0 ) 
        : TwoVars ( right , newname ) 
      {}
      // ======================================================================
      Log10* clone ( const char* newname ) const override 
      { return new Log10 ( *this , newname) ; }
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Erf
     *  Evaluate \f$ erf(ab) \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Erf final : public TwoVars
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::Erf , 2 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Erf 
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     , 
          RooAbsReal&        b     ) ;
      /// constructor with two variables 
      Erf 
        ( RooAbsReal&         a           , 
          RooAbsReal&         b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Erf ( name , title , a , b )
      {}
      /// constructor with two variables 
      Erf
        ( RooAbsReal&         a           , 
          const double        b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Erf ( name , title , a , RooFit::RooConst ( b ) )
      {}
      /// constructor with two variables 
      Erf
        ( const double        a           ,
          RooAbsReal&         b           , 
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Erf ( name , title , RooFit::RooConst ( a ) , b )
      {}
      /// constructor with one variable
      Erf 
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     ) 
        : Erf ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
      /// fake defautl constructor (needed for serisalization) 
      Erf () = default ;
      // ======================================================================
      // copy 
      Erf 
        ( const Erf& right , 
          const char* newname = 0 ) 
        : TwoVars ( right , newname ) 
      {}
      // ======================================================================
      Erf* clone ( const char* newname ) const override 
      { return new Erf ( *this , newname) ; }
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Erfc
     *  Evaluate \f$ erfc(ab) \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Erfc final : public TwoVars
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::Erfc , 2 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Erfc 
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     , 
          RooAbsReal&        b     ) ;
      /// constructor with two variables 
      Erfc 
        ( RooAbsReal&         a           , 
          RooAbsReal&         b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Erfc ( name , title , a , b )
      {}
      /// constructor with two variables 
      Erfc
        ( RooAbsReal&         a           , 
          const double        b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Erfc ( name , title , a , RooFit::RooConst ( b ) )
      {}
      /// constructor with two variables 
      Erfc 
        ( const double        a           ,
          RooAbsReal&         b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Erfc ( name , title , RooFit::RooConst ( a ) , b )
      {}
      /// constructor with one variable
      Erfc 
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     ) 
        : Erfc ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
      /// fake defautl constructor (needed for serisalization) 
      Erfc () = default ;
      // ======================================================================
      // copy 
      Erfc 
        ( const Erfc& right , 
          const char* newname = 0 ) 
        : TwoVars ( right , newname ) 
      {}
      // ======================================================================
      Erfc* clone ( const char* newname ) const override 
      { return new Erfc ( *this , newname) ; }
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Gamma
     *  Evaluate \f$ \Gamma(ab) \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Gamma final : public TwoVars 
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::Gamma , 2 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Gamma 
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     , 
          RooAbsReal&        b     ) ;
      /// constructor with two variables 
      Gamma
        ( RooAbsReal&         a           , 
          RooAbsReal&         b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Gamma ( name , title , a , b )
      {}
      /// constructor with two variables 
      Gamma
        ( RooAbsReal&         a           , 
          const double        b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Gamma ( name , title , a , RooFit::RooConst ( b ) )
      {}
      /// constructor with two variables 
      Gamma
        ( const double        a           ,
          RooAbsReal&         b           , 
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Gamma ( name , title , RooFit::RooConst ( a ) , b )
      {}
      /// constructor with one variable
      Gamma 
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     ) 
        : Gamma ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
      /// fake defautl constructor (needed for serisalization) 
      Gamma () = default ;
      // ======================================================================
      // copy 
      Gamma 
        ( const Gamma& right , 
          const char* newname = 0 ) 
        : TwoVars ( right , newname ) 
      {}
      // ======================================================================
      Gamma* clone ( const char* newname ) const override 
      { return new Gamma ( *this , newname) ; }
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;
      // ======================================================================
    }; //
    // ========================================================================
    /** @class LGamma
     *  Evaluate \f$ \log\Gamma(ab) \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class LGamma final : public TwoVars
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::LGamma , 2 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      LGamma  
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     , 
          RooAbsReal&        b     ) ;
      /// constructor with two variables 
      LGamma
        ( RooAbsReal&         a           , 
          RooAbsReal&         b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : LGamma ( name , title , a , b )
      {}
      /// constructor with two variables 
      LGamma 
        ( RooAbsReal&         a           , 
          const double        b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : LGamma ( name , title , a , RooFit::RooConst ( b ) )
      {}
      /// constructor with two variables 
      LGamma
        ( const double        a           ,
          RooAbsReal&         b           , 
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : LGamma ( name , title , RooFit::RooConst ( a ) , b )
      {}
      /// constructor with one variable
      LGamma
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     ) 
        : LGamma ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
      /// fake defautl constructor (needed for serisalization) 
      LGamma () = default ;
      // ======================================================================
      // copy 
      LGamma ( const LGamma& right , const char* newname = 0 ) 
        : TwoVars ( right , newname ) 
      {}
      // ======================================================================
      LGamma* clone ( const char* newname ) const override 
      { return new LGamma ( *this , newname ) ; }
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;
      // ======================================================================
    }; //
    // ========================================================================
    /** @class IGamma
     *  Evaluate \f$ \frac{1}{\Gamma(ab)} \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class IGamma final : public TwoVars 
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::IGamma , 2 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      IGamma 
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     , 
          RooAbsReal&        b     ) ;
      /// constructor with two variables 
      IGamma 
        ( RooAbsReal&         a           , 
          RooAbsReal&         b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : IGamma ( name , title , a , b )
      {}
      /// constructor with two variables 
      IGamma 
        ( RooAbsReal&         a           , 
          const double        b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : IGamma ( name , title , a , RooFit::RooConst ( b ) )
      {}
      /// constructor with two variables 
      IGamma
        ( const double        a           ,
          RooAbsReal&         b           , 
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : IGamma ( name , title , RooFit::RooConst ( a ) , b )
      {}
      /// constructor with one variable
      IGamma
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     ) 
        : IGamma ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
      /// fake defautl constructor (needed for serisalization) 
      IGamma () = default ;
      // ======================================================================
      // copy 
      IGamma 
        ( const IGamma& right , 
          const char* newname = 0 ) 
        : TwoVars ( right , newname ) 
      {}
      // ======================================================================
      IGamma* clone ( const char* newname ) const override 
      { return new IGamma ( *this , newname ) ; }
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Sin
     *  Evaluate \f$ \sin ab  \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Sin final : public TwoVars 
    {
      // ======================================================================
      ClassDefOverride(Ostap::MoreRooFit::Sin , 2 ) ;  // power function
      // ======================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Sin 
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     , 
          RooAbsReal&        b     ) ;
      /// constructor with two variables
      Sin 
        ( RooAbsReal&         a           , 
          RooAbsReal&         b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Sin ( name , title , a , b )
      {}
      /// constructor with two variables
      Sin 
        ( RooAbsReal&         a           , 
          const double        b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Sin ( name , title , a , RooFit::RooConst ( b ) )
      {}
      /// constructor with two variables
      Sin 
        ( const double        a           ,
          RooAbsReal&         b           , 
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Sin ( name , title , RooFit::RooConst ( a ) , b )
      {}
      /// constructor with one variable
      Sin 
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     ) 
        : Sin ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
      /// fake defautl constructor (needed for serisalization) 
      Sin () = default ;
      // ======================================================================
      // copy 
      Sin
        ( const Sin& right , 
          const char* newname = 0 ) 
        : TwoVars ( right , newname ) 
      {}
      // ======================================================================
      Sin* clone ( const char* newname ) const override 
      { return new Sin ( *this , newname ) ; }
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Cos
     *  Evaluate \f$ \cos ab  \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Cos final : public TwoVars 
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::Cos , 2 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Cos 
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     , 
          RooAbsReal&        b     ) ;
      /// constructor with two variables
      Cos 
        ( RooAbsReal&         a           , 
          RooAbsReal&         b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Cos ( name , title , a , b )
      {}
      /// constructor with two variables
      Cos
        ( RooAbsReal&         a           , 
          const double        b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Cos ( name , title , a , RooFit::RooConst ( b ) )
      {}
      /// constructor with two variables
      Cos 
        ( const double        a           ,
          RooAbsReal&         b           , 
          const std::string&  name  = ""  ,            
          const std::string&  title = ""  ) 
        : Cos ( name , title , RooFit::RooConst ( a ) , b )
      {}
      /// constructor with one variable
      Cos 
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     ) 
        : Cos ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
      /// fake defautl constructor (needed for serisalization) 
      Cos () = default ;
      // ======================================================================
      // copy 
      Cos ( const Cos& right , const char* newname = 0 ) 
        : TwoVars ( right , newname ) 
      {}
      // ======================================================================
      Cos* clone ( const char* newname ) const override 
      { return new Cos ( *this , newname ) ; }
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;
      // ======================================================================
    };
    // =========================================================================
    /** @class Tan
     *  Evaluate \f$ \tan ab  \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Tan final : public TwoVars
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::Tan , 2 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Tan 
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     , 
          RooAbsReal&        b     ) ;
      /// constructor with two variables
      Tan 
        ( RooAbsReal&         a           , 
          RooAbsReal&         b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Tan ( name , title , a , b )
      {}
      /// constructor with two variables
      Tan 
        ( RooAbsReal&         a           , 
          const double        b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Tan ( name , title , a , RooFit::RooConst ( b ) )
      {}
      /// constructor with two variables
      Tan
        ( const double        a           ,
          RooAbsReal&         b           , 
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Tan ( name , title , RooFit::RooConst ( a ) , b )
      {}
      /// constructor with one variable
      Tan
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     ) 
        : Tan ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
      /// fake defautl constructor (needed for serisalization) 
      Tan () = default ;
      // ======================================================================
      // copy 
      Tan 
        ( const Tan& right , 
          const char* newname = 0 ) 
        : TwoVars ( right , newname ) 
      {}
      // ======================================================================
      Tan* clone ( const char* newname ) const override 
      { return new Tan ( *this , newname ) ; }
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Sinh
     *  Evaluate \f$ \sinh ab  \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Sinh final : public TwoVars 
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::Sinh , 2 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Sinh 
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     , 
          RooAbsReal&        b     ) ;
      /// constructor with two variables
      Sinh 
        ( RooAbsReal&         a           , 
          RooAbsReal&         b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Sinh ( name , title , a , b )
      {}
      /// constructor with two variables
      Sinh 
        ( RooAbsReal&         a           , 
          const double        b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Sinh ( name , title , a , RooFit::RooConst ( b ) ) 
      {}
      /// constructor with two variables
      Sinh 
        ( const double        a           ,
          RooAbsReal&         b           , 
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Sinh ( name , title , RooFit::RooConst ( a ) , b ) 
      {}
      /// constructor with one variable
      Sinh 
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     ) 
        : Sinh ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
      /// fake defautl constructor (needed for serisalization) 
      Sinh () = default ;
      // ======================================================================
      // copy 
      Sinh 
        ( const Sinh& right , 
          const char* newname = 0 ) 
        : TwoVars ( right , newname ) 
      {}
      // ======================================================================
      Sinh* clone ( const char* newname ) const override 
      { return new Sinh ( *this , newname ) ; }
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Cosh
     *  Evaluate \f$ \cosh ab  \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Cosh final : public TwoVars 
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::Cosh , 2 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Cosh  
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     , 
          RooAbsReal&        b     ) ;
      /// constructor with two variables
      Cosh 
        ( RooAbsReal&         a           , 
          RooAbsReal&         b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Cosh ( name , title , a , b )
      {}
      /// constructor with two variables
      Cosh
        ( RooAbsReal&         a           , 
          const double        b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Cosh ( name , title , a , RooFit::RooConst ( b ) ) 
      {}
      /// constructor with two variables
      Cosh 
        ( const double        a           ,
          RooAbsReal&         b           , 
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Cosh ( name , title , RooFit::RooConst ( a ) , b ) 
      {}
      /// constructor with one variable
      Cosh
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     ) 
        : Cosh ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
      /// fake defautl constructor (needed for serisalization)
      Cosh () = default ;
      // ======================================================================
      // copy 
      Cosh 
        ( const Cosh& right , 
          const char* newname = 0 ) 
        : TwoVars ( right , newname ) 
      {}
      // ======================================================================
      Cosh* clone ( const char* newname ) const override 
      { return new Cosh ( *this , newname ) ; }
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;
      // ======================================================================
    };
    // ========================================================================
    /** @class Tanh
     *  Evaluate \f$ \tanh ab  \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Tanh final : public TwoVars
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::Tanh , 2 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Tanh
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     , 
          RooAbsReal&        b     ) ;
      /// constructor with two variables
      Tanh 
        ( RooAbsReal&         a           , 
          RooAbsReal&         b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Tanh ( name , title , a , b )
      {}
      /// constructor with two variables
      Tanh 
        ( RooAbsReal&         a           , 
          const double        b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Tanh ( name , title , a , RooFit::RooConst ( b ) )
      {}
      /// constructor with two variables
      Tanh 
        ( const double        a           ,
          RooAbsReal&         b           , 
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Tanh ( name , title , RooFit::RooConst ( a ) , b )
      {}
      /// constructor with one variable
      Tanh 
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     ) 
        : Tanh ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
      /// fake defautl constructor (needed for serisalization)
      Tanh () = default ;
      // ======================================================================
      // copy 
      Tanh ( const Tanh& right , const char* newname = 0 ) 
        : TwoVars ( right , newname ) 
      {}
      // ======================================================================
      Tanh* clone ( const char* newname ) const override 
      { return new Tanh ( *this , newname ) ; }
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Sech
     *  Evaluate \f$ \sech ab = \frac{1}{ \cosh a b }   \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Sech final : public TwoVars
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::Sech , 2 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Sech 
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     , 
          RooAbsReal&        b     ) ;
      /// constructor with two variables
      Sech 
        ( RooAbsReal&         a           , 
          RooAbsReal&         b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Sech ( name , title , a , b )
      {}
      /// constructor with two variables
      Sech 
        ( RooAbsReal&         a           , 
          const double        b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Sech ( name , title , a , RooFit::RooConst ( b ) ) 
      {}
      /// constructor with two variables
      Sech 
        ( const double        a           ,
          RooAbsReal&         b           , 
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Sech ( name , title , RooFit::RooConst ( a ) , b ) 
      {}
      /// constructor with one variable
      Sech
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     ) 
        : Sech ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
      /// fake defautl constructor (needed for serisalization)
      Sech () = default ;
      // ======================================================================
      // copy 
      Sech 
        ( const Sech& right , 
          const char* newname = 0 ) 
        : TwoVars ( right , newname ) 
      {}
      // ======================================================================
      Sech* clone ( const char* newname ) const override 
      { return new Sech ( *this , newname ) ; }
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Atan2
     *  Evaluate \f$ atan2 ( a , b)  \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Atan2 final : public TwoVars 
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::Atan2 , 2 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Atan2 
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     , 
          RooAbsReal&        b     ) ;
      /// constructor with two variables
      Atan2 
        ( RooAbsReal&         a           , 
          RooAbsReal&         b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Atan2 ( name , title , a , b )
      {}
      /// constructor with two variables
      Atan2
        ( RooAbsReal&         a           , 
          const double        b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Atan2 ( name , title , a , RooFit::RooConst ( b ) )
      {}
      /// constructor with two variables
      Atan2
        ( const double        a           ,
          RooAbsReal&         b           , 
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : Atan2 ( name , title , RooFit::RooConst ( a ) , b )
      {}
      /// constructor with one variable
      Atan2 
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     ) 
        : Atan2 ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
      /// fake defautl constructor (needed for serisalization)
      Atan2 () = default ;
      // ======================================================================
      // copy 
      Atan2 
        ( const Atan2& right , 
          const char* newname = 0 ) 
        : TwoVars ( right , newname ) 
      {}
      // ======================================================================
      Atan2* clone ( const char* newname ) const override 
      { return new Atan2 ( *this , newname ) ; }
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;
      // ======================================================================
    }; //
    // ========================================================================    
    /** @class MaxV
     *  Evaluate \f$ max ( a , b ) \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class MaxV final : public TwoVars 
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::MaxV , 2 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      MaxV 
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     , 
          RooAbsReal&        b     ) ;
      /// constructor with two variables
      MaxV 
        ( RooAbsReal&         a           , 
          RooAbsReal&         b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : MaxV ( name , title , a , b )
      {}
      /// constructor with two variables
      MaxV 
        ( RooAbsReal&         a           , 
          const double        b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : MaxV ( name , title , a , RooFit::RooConst ( b ) )
      {}
      /// constructor with two variables
      MaxV
        ( const double        a           ,
          RooAbsReal&         b           , 
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : MaxV ( name , title , RooFit::RooConst ( a ) , b )
      {}
      // ======================================================================
      /// fake defautl constructor (needed for serisalization)
      MaxV () = default ;
      // ======================================================================
      // copy 
      MaxV 
        ( const MaxV& right , 
          const char* newname = 0 ) 
        : TwoVars ( right , newname ) 
      {}
      // ======================================================================
      MaxV* clone ( const char* newname ) const override 
      { return new MaxV ( *this , newname ) ; }
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;
      // ======================================================================
    }; //
    // ========================================================================    
    /** @class MinV
     *  Evaluate \f$ min ( a , b ) \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class MinV final : public TwoVars 
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::MinV , 2 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      MinV 
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     , 
          RooAbsReal&        b     ) ;
      /// constructor with two variables
      MinV 
        ( RooAbsReal&         a           , 
          RooAbsReal&         b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : MinV ( name , title , a , b )
      {}
      /// constructor with two variables
      MinV 
        ( RooAbsReal&         a           , 
          const double        b           ,
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : MinV ( name , title , a , RooFit::RooConst ( b ) )
      {}
      /// constructor with two variables
      MinV
        ( const double        a           ,
          RooAbsReal&         b           , 
          const std::string&  name  = ""  , 
          const std::string&  title = ""  ) 
        : MinV ( name , title , RooFit::RooConst ( a ) , b )
      {}
      // ======================================================================
      /// fake defautl constructor (needed for serisalization)
      MinV () = default ;
      // ======================================================================
      // copy 
      MinV 
        ( const MinV& right , 
          const char* newname = 0 ) 
        : TwoVars ( right , newname ) 
      {}
      // ======================================================================
      MinV* clone ( const char* newname ) const override 
      { return new MinV ( *this , newname ) ; }
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Id
     *  Trivial variable: identical transformation 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2020-03-17
     */ 
    class Id final : public OneVar 
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::Id , 1 ) ;  // sum of RooAbsReal objects
      // ========================================================================
    public:
      // ========================================================================
      Id () = default ;
      /// constructor with the name, title and variable 
      Id 
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        a     ) ;
      /// constructor with the name and variable 
      Id 
        ( const std::string& name  , 
          RooAbsReal&        a     ) ;
      /// copy 
      Id 
        ( const Id&   right       , 
          const char* newname = 0 ) ;
      /// destructor 
      virtual ~Id () ;
      /// clone 
      Id* clone ( const char* newname ) const override ;
      // ======================================================================
    public: // delegation 
      // ======================================================================
      Double_t    analyticalIntegral
      ( Int_t            code            ,
        const char*      range = nullptr ) const override ;
      Double_t    analyticalIntegralWN
      ( Int_t            code            ,
        const RooArgSet* normset         ,
        const char*      range = nullptr ) const override ;
      //
      Int_t    getAnalyticalIntegral
      ( RooArgSet&       allVars         ,
        RooArgSet&       analVars        ,
        const char*      range = nullptr ) const override ;
      //
      Int_t    getAnalyticalIntegralWN
      ( RooArgSet&       allVars         ,
        RooArgSet&       analVars        ,
        const RooArgSet* normset         ,
        const char*      range = nullptr ) const override ;
      // ======================================================================
    public: 
      // ======================================================================
      /// the actual evaluation of the result
      double   getValV  ( const RooArgSet* aset = nullptr ) const override 
      { return m_x.arg ().getValV ( aset ) ; }
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override { return m_x ; }
      // ======================================================================
    }; //
    // ========================================================================    
    /** @class AddDeps
     *  helper class to add a fictive dependency of one varibale on others 
     */
    class AddDeps final : public OneVar 
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::AddDeps , 1 ) ;  
      // ========================================================================
    public:
      // ========================================================================
      /// constructor with a variable
      AddDeps
        ( const std::string& name  , 
          const std::string& title , 
          RooAbsReal&        var   , 
          const RooArgList&  deps  ) ; // list of fake dependants 
      /// copy 
      AddDeps
        ( const AddDeps&     right       , 
          const char* newname = 0 ) ;
      /// fake default constructor (needed for serisalization)
      AddDeps () = default ;
      /// destructor 
      virtual ~AddDeps () ;
      /// clone 
      AddDeps* clone ( const char* newname ) const override ;
      // ======================================================================    
    public:
      // ======================================================================
      /// the actual evaluation of the result 
      Double_t evaluate () const override { return m_x ; }
      /// the actual evaluation of the result
      double   getValV  ( const RooArgSet* aset = nullptr ) const override 
      { return m_x.arg ().getValV ( aset ) ; }
      // ======================================================================
#if ROOT_VERSION(6,24,0)<=ROOT_VERSION_CODE
      // ======================================================================
      RooSpan<const double>
      getValues ( RooBatchCompute::RunContext& evalData , 
                  const RooArgSet*             aset     ) const override 
      {  return m_x.arg ().getValues ( evalData , aset ) ; }
      // ======================================================================
#endif 
      // ======================================================================
    public:
      // ======================================================================
      /// get the list 
      const RooArgList& vlst() const { return m_vlst ; }
      // ======================================================================
    private:
      // ======================================================================
      /// variables 
      RooListProxy m_vlst  {}    ; // variables 
      // ======================================================================
    } ;
    // ========================================================================    
    /** @class ProductPdf
     *  Oversimplified product of two PDF
     *  - It is useful to bypass some "features" of RooFit
     *  - it can be rather inefficient 
     *  @attention all native roofit optimisations are not applied 
     *  for this case!      
     */
    class ProductPdf : public RooAbsPdf
    {
      // ========================================================================
      ClassDefOverride(Ostap::MoreRooFit::ProductPdf , 1 ) ;  // sum of RooAbsReal objects
      // ========================================================================
    public:
      // ======================================================================== 
      /** constructor from name, title and two pdfs
       *  @param name  name 
       *  @param title name 
       *  @param pdf1 the first pdf 
       *  @param pdf2 the second pdf 
       */
      ProductPdf 
      ( const char* name  , 
        const char* title , 
        RooAbsPdf&  pdf1  , 
        RooAbsPdf&  pdf2  ) ;
      /// copy constructor 
      ProductPdf
      ( const ProductPdf& right    , 
        const char*       name = 0 ) ;
      /// destructor 
      virtual ~ProductPdf() ;
      /// clone 
      ProductPdf* clone ( const char* newname ) const override ;
      // ========================================================================
      /// fake defautl constructor (needed for serisalization)
      ProductPdf () = default ;
      // ======================================================================
    protected:
      // ========================================================================
      /// the main method 
      Double_t evaluate () const override ;
      // ========================================================================
    public:
      // ========================================================================
      inline const RooAbsReal& x () const { return m_pdf1.arg() ; }
      inline const RooAbsReal& y () const { return m_pdf2.arg() ; }    
      // ========================================================================
    protected:
      // ========================================================================
      /// pdf1
      RooRealProxy m_pdf1 ; // the first pdf 
      /// pdf2
      RooRealProxy m_pdf2 ; // the second pdf 
      // ======================================================================
    } ;
    // ========================================================================    
    /** @class WrapPdf 
     *  backport of <code>RooWrapperPdf</code> from new ROOT 
     * @see RooWrapperPdf 
     */
    class WrapPdf : public RooAbsPdf 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::MoreRooFit::WrapPdf,1)
      // ======================================================================
      /// constructor 
      WrapPdf 
      ( const char *name  , 
        const char *title , 
        RooAbsReal& func  ) ;
      /// copy contructor 
      WrapPdf 
      ( const WrapPdf& other    , 
        const char*    name = 0 ) ;
      /// cloning method 
      WrapPdf* clone( const char* newname )  const override ;
      // ======================================================================
    public :
      // ======================================================================
      /// constructor needed fo rserialization 
      WrapPdf () {} ;
      // ======================================================================
    public: 
      // ====================================================================== 
      // Analytical Integration handling
      bool forceAnalyticalInt
      ( const RooAbsArg& dep ) const override ;
      Int_t getAnalyticalIntegralWN
      ( RooArgSet&       allVars       , 
        RooArgSet&       analVars      , 
        const RooArgSet* normSet       ,
        const char*      rangeName = 0 ) const override ;
      // ======================================================================
      Int_t getAnalyticalIntegral
      ( RooArgSet&       allVars       , 
        RooArgSet&       numVars       ,
        const char*      rangeName = 0 ) const override ;
      // ======================================================================
      // Hints for optimized brute-force sampling
      Int_t getMaxVal
      ( const RooArgSet& vars     ) const override ;
      // ======================================================================
      double maxVal
      ( Int_t            code     ) const override ;
      // ======================================================================
      Int_t minTrialSamples
      ( const RooArgSet& arGenObs ) const override ;
      // Plotting and binning hints
      bool isBinnedDistribution 
      ( const RooArgSet& obs     ) const override ;
      // ======================================================================
      std::list<double>* binBoundaries
      ( RooAbsRealLValue& obs , 
        double            xlo , 
        double            xhi ) const override ;
      // ======================================================================
      std::list<double>* plotSamplingHint
      ( RooAbsRealLValue& obs , 
        double            xlo , 
        double            xhi ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      // evaluate PDF 
      Double_t evaluate ()  const override { return m_func ; }
      // ======================================================================
    private:
      // ======================================================================
      /// the function 
      RooRealProxy  m_func {} ;
      // ======================================================================
    }; //                            The end of clas Ostap::MoreRooFit::WrapPdf 
    // ========================================================================    
  } //                                   The end of namespace Ostap::MoreRooFit  
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_MOREROOFIT_H
// ============================================================================
