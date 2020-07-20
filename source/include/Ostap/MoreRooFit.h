// ============================================================================
#ifndef OSTAP_MOREROOFIT_H 
#define OSTAP_MOREROOFIT_H 1
// ============================================================================
// Include files
// ============================================================================
// ROOT/RooFit 
// ============================================================================
#include "RVersion.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooRealConstant.h"
#include "RooRealProxy.h"
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
      ClassDef(Ostap::MoreRooFit::Addition , 1 ) ;  // sum of RooAbsReal objects
      // ========================================================================
    public:
      // ========================================================================
      Addition () = default ;
      /// constructor with two variables 
      Addition ( const std::string& name  , 
                 const std::string& title , 
                 RooAbsReal&        a     , 
                 RooAbsReal&        b     ) ;
      /// constructor with two variables 
      Addition ( RooAbsReal&        a           , 
                 RooAbsReal&        b           ,
                 const std::string& name  = ""  , 
                 const std::string& title = ""  ) 
        : Addition ( name , title , a , b )
      {}
      /// copy 
      Addition ( const Addition&    right       , 
                 const char*        newname = 0 ) ;
      /// destructor
      virtual ~Addition () ;
      /// clone 
      Addition* clone ( const char* newname ) const override ;
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
      ClassDef(Ostap::MoreRooFit::Product , 1 ) ;  // sum of RooAbsReal objects
      // ========================================================================
    public:
      // ========================================================================
      Product () = default ;
      /// constructor with two variables 
      Product ( const std::string& name  , 
                const std::string& title , 
                RooAbsReal&        a     , 
                RooAbsReal&        b     ) ;
      /// constructor with two variables 
      Product ( RooAbsReal&        a           , 
                RooAbsReal&        b           ,
                const std::string& name  = ""  , 
                const std::string& title = ""  ) 
        : Product ( name , title , a , b )
      {}
      /// copy 
      Product ( const Product&    right       , 
                const char*        newname = 0 ) ;
      /// destructor 
      virtual ~Product () ;
      /// clone 
      Product* clone ( const char* newname ) const override ;
      // ======================================================================
    }; // 
    // ========================================================================
    /** @class Subtraction 
     *  A simple modification of class RooAddition
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */ 
    class Subtraction : public Addition
    {
      // ========================================================================
      ClassDef(Ostap::MoreRooFit::Subtraction , 1 ) ;  // Difference of RooAbsReal objects
      // ========================================================================
    public:
      // ========================================================================
      Subtraction () = default ;
      /// constructor with list of variables 
      Subtraction ( const std::string& name  , 
                    const std::string& title , 
                    RooAbsReal&        a     , 
                    RooAbsReal&        b     ) ;
      /// constructor with two variables 
      Subtraction ( RooAbsReal&        a           , 
                    RooAbsReal&        b           ,
                    const std::string& name  = ""  , 
                    const std::string& title = ""  ) 
        : Subtraction ( name , title , a , b )
      {}
      /// copy 
      Subtraction ( const Subtraction& right       , 
                    const char*        newname = 0 ) ;
      /// destructir
      virtual ~Subtraction() ;
      /// clone 
      Subtraction* clone ( const char* newname ) const override ;
      /// integrals 
      Double_t analyticalIntegral ( Int_t code ,
                                    const char* rangeName = 0 ) const override ;
      Double_t evaluate           () const override ;    
      // ======================================================================
    }; // 
    // ========================================================================
    /** @class OneVar
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2020-03-17
     */ 
    class OneVar : public RooAbsReal 
    {
      // ========================================================================
      ClassDef(Ostap::MoreRooFit::OneVar , 1 ) ;  
      // ========================================================================
    public:
      // ========================================================================
      OneVar () = default ;
      /// constructor with a variable
      OneVar ( const std::string& name  , 
               const std::string& title , 
               RooAbsReal&        x     ) ;
      /// copy 
      OneVar ( const OneVar&   right       , 
               const char* newname = 0 ) ;
      /// destructor 
      virtual ~OneVar () ;
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
      ClassDef ( Ostap::MoreRooFit::FunOneVar , 1 ) ;
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
      ClassDef(Ostap::MoreRooFit::TwoVars , 1 ) ;  
      // ========================================================================
    public:
      // ========================================================================
      TwoVars () = default ;
      /// constructor with a variable
      TwoVars ( const std::string& name  , 
                const std::string& title , 
                RooAbsReal&        x     ,
                RooAbsReal&        y     ) ;
      /// copy 
      TwoVars ( const TwoVars& right       , 
                const char*    newname = 0 ) ;
      /// destructor 
      virtual ~TwoVars () ;
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
      ClassDef(Ostap::MoreRooFit::FunTwoVars , 1 ) ;  
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
      /// copy 
      FunTwoVars ( const FunTwoVars& right       , 
                   const char*       newname = 0 ) ;
      /// destructor 
      virtual ~FunTwoVars () ;
      /// clone 
      FunTwoVars* clone  ( const char* newname ) const override ;
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
    class Division final : public FunTwoVars 
    {
      // ========================================================================
      ClassDef(Ostap::MoreRooFit::Division , 1 ) ;  // ratio of RooAbsReal objects
      // ========================================================================
    public:
      // ======================================================================
      Division  ( const std::string& name  , 
                  const std::string& title , 
                  RooAbsReal&        x     , 
                  RooAbsReal&        y     ) ;
      /// constructor with two variables 
      Division ( RooAbsReal&         x           , 
                 RooAbsReal&         y           ,
                 const std::string&  name  = ""  , 
                 const std::string&  title = ""  ) 
        : Division ( name , title , x , y )
      {}
      // ======================================================================
      /// default constructor 
      Division  () =  default ;
      // ======================================================================
      // copy 
      Division 
        ( const Division& right       ,
          const char*     newname = 0 ) 
        : FunTwoVars ( right , newname ) 
      {}
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Fraction
     *  Evaluate \f$ \frac{a}{a+b}\f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Fraction final : public FunTwoVars
    {
      // ======================================================================
      ClassDef(Ostap::MoreRooFit::Fraction , 1 ) ;  // Fraction 
      // ======================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Fraction  ( const std::string& name  , 
                  const std::string& title , 
                  RooAbsReal&        x     , 
                  RooAbsReal&        y     ) ;
      /// constructor with two variables 
      Fraction ( RooAbsReal&         x           , 
                 RooAbsReal&         y           ,
                 const std::string&  name  = ""  , 
                 const std::string&  title = ""  ) 
        : Fraction ( name , title , x , y )
      {}
      // ======================================================================
      /// default constructor 
      Fraction  () =  default ;
      // ======================================================================
      // copy 
      Fraction ( const Fraction& right       , 
                 const char*     newname = 0 ) 
        : FunTwoVars ( right , newname ) 
      {}
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Asymmetry
     *  Evaluate \f$ \frac{a-b}{a+b}\f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Asymmetry: public FunTwoVars 
    {
      // ======================================================================
      ClassDef(Ostap::MoreRooFit::Asymmetry , 1 ) ;  // Relative difference 
      // ======================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Asymmetry  ( const std::string& name  , 
                   const std::string& title , 
                   RooAbsReal&        a     , 
                   RooAbsReal&        b     ) ;
      /// constructor with two variables 
      Asymmetry ( RooAbsReal&         a           , 
                  RooAbsReal&         b           ,
                  const std::string&  name  = ""  , 
                  const std::string&  title = ""  ) 
        : Asymmetry ( name , title , a , b )
      {}
      // ======================================================================
      /// default constructor 
      Asymmetry () =  default ;
      // ======================================================================
      // copy 
      Asymmetry ( const Asymmetry&    right       , 
                  const char*         newname = 0 ) 
        : FunTwoVars ( right , newname ) 
      {}
      // ======================================================================
      // destructor 
      virtual ~Asymmetry () {};
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Power
     *  Evaluate \f$ a^b \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Power: public FunTwoVars 
    {
      // ========================================================================
      ClassDef(Ostap::MoreRooFit::Power , 1 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Power  ( const std::string& name  , 
               const std::string& title , 
               RooAbsReal&        a     , 
               RooAbsReal&        b     ) ;
      /// constructor with two variables 
      Power ( RooAbsReal&         a           , 
              RooAbsReal&         b           ,
              const std::string&  name  = ""  , 
              const std::string&  title = ""  ) 
        : Power ( name , title , a , b )
      {}
      // ======================================================================
      /// default constructor 
      Power () =  default ;
      // ======================================================================
      // copy 
      Power ( const Power&   right       , 
              const char*    newname = 0 ) 
        : FunTwoVars ( right , newname ) 
      {}
      // ======================================================================
      // destructor 
      virtual ~Power() {} ;
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Abs
     *  Evaluate \f$ abs(ab) \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Abs final : public FunTwoVars 
    {
      // ========================================================================
      ClassDef(Ostap::MoreRooFit::Abs , 1 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Abs  ( const std::string& name  , 
             const std::string& title , 
             RooAbsReal&        a     , 
             RooAbsReal&        b     ) ;
      /// constructor with two variables 
      Abs ( RooAbsReal&         a           , 
            RooAbsReal&         b           ,
            const std::string&  name  = ""  , 
            const std::string&  title = ""  ) 
        : Abs ( name , title , a , b )
      {}
      /// constructor with one variable
      Abs  ( const std::string& name  , 
             const std::string& title , 
             RooAbsReal&        a     ) 
        : Abs ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
      /// default constructor 
      Abs () =  default ;
      // ======================================================================
      // copy 
      Abs ( const Abs&     right       , 
            const char*    newname = 0 ) 
        : FunTwoVars ( right , newname ) 
      {}
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Exp
     *  Evaluate \f$ exp(ab) \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Exp final : public FunTwoVars 
    {
      // ========================================================================
      ClassDef(Ostap::MoreRooFit::Exp , 1 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Exp  ( const std::string& name  , 
             const std::string& title , 
             RooAbsReal&        a     , 
             RooAbsReal&        b     ) ;
      /// constructor with two variables 
      Exp ( RooAbsReal&         a           , 
            RooAbsReal&         b           ,
            const std::string&  name  = ""  , 
            const std::string&  title = ""  ) 
        : Exp ( name , title , a , b )
      {}
      /// constructor with one variable
      Exp  ( const std::string& name  , 
             const std::string& title , 
             RooAbsReal&        a     ) 
        : Exp ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Log
     *  Evaluate \f$ \log ab  \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Log final : public FunTwoVars
    {
      // ========================================================================
      ClassDef(Ostap::MoreRooFit::Log , 1 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Log  ( const std::string& name  , 
             const std::string& title , 
             RooAbsReal&        a     , 
             RooAbsReal&        b     ) ;
      /// constructor with two variables
      Log ( RooAbsReal&         a           , 
            RooAbsReal&         b           ,
            const std::string&  name  = ""  , 
            const std::string&  title = ""  ) 
        : Log  ( name , title , a , b )
      {}
      /// constructor with one variable
      Log  ( const std::string& name  , 
             const std::string& title , 
             RooAbsReal&        a     ) 
        : Log ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Log10
     *  Evaluate \f$ \log10 ab  \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Log10 final : public FunTwoVars 
    {
      // ========================================================================
      ClassDef(Ostap::MoreRooFit::Log10 , 1 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Log10  ( const std::string& name  , 
               const std::string& title , 
               RooAbsReal&        a     , 
               RooAbsReal&        b     ) ;
      /// constructor with two variables
      Log10 ( RooAbsReal&         a           , 
              RooAbsReal&         b           ,
              const std::string&  name  = ""  , 
              const std::string&  title = ""  ) 
        : Log10  ( name , title , a , b )
      {}
      /// constructor with one variable
      Log10  ( const std::string& name  , 
               const std::string& title , 
               RooAbsReal&        a     ) 
        : Log10 ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Erf
     *  Evaluate \f$ erf(ab) \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Erf final : public FunTwoVars
    {
      // ========================================================================
      ClassDef(Ostap::MoreRooFit::Erf , 1 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Erf  ( const std::string& name  , 
             const std::string& title , 
             RooAbsReal&        a     , 
             RooAbsReal&        b     ) ;
      /// constructor with two variables 
      Erf ( RooAbsReal&         a           , 
            RooAbsReal&         b           ,
            const std::string&  name  = ""  , 
            const std::string&  title = ""  ) 
        : Erf ( name , title , a , b )
      {}
      /// constructor with one variable
      Erf  ( const std::string& name  , 
             const std::string& title , 
             RooAbsReal&        a     ) 
        : Erf ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Gamma
     *  Evaluate \f$ \Gamma(ab) \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Gamma: public FunTwoVars 
    {
      // ========================================================================
      ClassDef(Ostap::MoreRooFit::Gamma , 1 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Gamma  ( const std::string& name  , 
             const std::string& title , 
             RooAbsReal&        a     , 
             RooAbsReal&        b     ) ;
      /// constructor with two variables 
      Gamma ( RooAbsReal&         a           , 
              RooAbsReal&         b           ,
              const std::string&  name  = ""  , 
              const std::string&  title = ""  ) 
        : Gamma ( name , title , a , b )
      {}
      /// constructor with one variable
      Gamma  ( const std::string& name  , 
               const std::string& title , 
               RooAbsReal&        a     ) 
        : Gamma ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
    }; //
    // ========================================================================
    /** @class LGamma
     *  Evaluate \f$ \log\Gamma(ab) \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class LGamma: public FunTwoVars
    {
      // ========================================================================
      ClassDef(Ostap::MoreRooFit::LGamma , 1 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      LGamma  ( const std::string& name  , 
                const std::string& title , 
                RooAbsReal&        a     , 
                RooAbsReal&        b     ) ;
      /// constructor with two variables 
      LGamma ( RooAbsReal&         a           , 
               RooAbsReal&         b           ,
               const std::string&  name  = ""  , 
               const std::string&  title = ""  ) 
        : LGamma ( name , title , a , b )
      {}
      /// constructor with one variable
      LGamma  ( const std::string& name  , 
                const std::string& title , 
                RooAbsReal&        a     ) 
        : LGamma ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
    }; //
    // ========================================================================
    /** @class IGamma
     *  Evaluate \f$ \frac{1}{\Gamma(ab)} \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class IGamma: public FunTwoVars 
    {
      // ========================================================================
      ClassDef(Ostap::MoreRooFit::IGamma , 1 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      IGamma  ( const std::string& name  , 
                const std::string& title , 
                RooAbsReal&        a     , 
                RooAbsReal&        b     ) ;
      /// constructor with two variables 
      IGamma ( RooAbsReal&         a           , 
               RooAbsReal&         b           ,
               const std::string&  name  = ""  , 
               const std::string&  title = ""  ) 
        : IGamma ( name , title , a , b )
      {}
      /// constructor with one variable
      IGamma  ( const std::string& name  , 
                const std::string& title , 
                RooAbsReal&        a     ) 
        : IGamma ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Sin
     *  Evaluate \f$ \sin ab  \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Sin: public FunTwoVars 
    {
      // ========================================================================
      ClassDef(Ostap::MoreRooFit::Sin , 1 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Sin  ( const std::string& name  , 
             const std::string& title , 
             RooAbsReal&        a     , 
             RooAbsReal&        b     ) ;
      /// constructor with two variables
      Sin ( RooAbsReal&         a           , 
            RooAbsReal&         b           ,
            const std::string&  name  = ""  , 
            const std::string&  title = ""  ) 
        : Sin ( name , title , a , b )
      {}
      /// constructor with one variable
      Sin ( const std::string& name  , 
            const std::string& title , 
            RooAbsReal&        a     ) 
        : Sin ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Cos
     *  Evaluate \f$ \cos ab  \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Cos: public FunTwoVars 
    {
      // ========================================================================
      ClassDef(Ostap::MoreRooFit::Cos , 1 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Cos  ( const std::string& name  , 
             const std::string& title , 
             RooAbsReal&        a     , 
             RooAbsReal&        b     ) ;
      /// constructor with two variables
      Cos ( RooAbsReal&         a           , 
            RooAbsReal&         b           ,
            const std::string&  name  = ""  , 
            const std::string&  title = ""  ) 
        : Cos ( name , title , a , b )
      {}
      /// constructor with one variable
      Cos ( const std::string& name  , 
            const std::string& title , 
            RooAbsReal&        a     ) 
        : Cos ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
    };
    // =========================================================================
    /** @class Tan
     *  Evaluate \f$ \tan ab  \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Tan final : public FunTwoVars
    {
      // ========================================================================
      ClassDef(Ostap::MoreRooFit::Tan , 1 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Tan  ( const std::string& name  , 
             const std::string& title , 
             RooAbsReal&        a     , 
             RooAbsReal&        b     ) ;
      /// constructor with two variables
      Tan ( RooAbsReal&         a           , 
            RooAbsReal&         b           ,
            const std::string&  name  = ""  , 
            const std::string&  title = ""  ) 
        : Tan ( name , title , a , b )
      {}
      /// constructor with one variable
      Tan ( const std::string& name  , 
            const std::string& title , 
            RooAbsReal&        a     ) 
        : Tan ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Tanh
     *  Evaluate \f$ \tanh ab  \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Tanh: public FunTwoVars
    {
      // ========================================================================
      ClassDef(Ostap::MoreRooFit::Tanh , 1 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Tanh  ( const std::string& name  , 
              const std::string& title , 
              RooAbsReal&        a     , 
              RooAbsReal&        b     ) ;
      /// constructor with two variables
      Tanh ( RooAbsReal&         a           , 
             RooAbsReal&         b           ,
             const std::string&  name  = ""  , 
             const std::string&  title = ""  ) 
        : Tanh ( name , title , a , b )
      {}
      /// constructor with one variable
      Tanh ( const std::string& name  , 
             const std::string& title , 
             RooAbsReal&        a     ) 
        : Tanh ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Atan2
     *  Evaluate \f$ atan2 ( a , b)  \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Atan2: public FunTwoVars 
    {
      // ========================================================================
      ClassDef(Ostap::MoreRooFit::Atan2 , 1 ) ;  // power function
      // ========================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Atan2  ( const std::string& name  , 
              const std::string& title , 
              RooAbsReal&        a     , 
              RooAbsReal&        b     ) ;
      /// constructor with two variables
      Atan2 ( RooAbsReal&         a           , 
              RooAbsReal&         b           ,
              const std::string&  name  = ""  , 
              const std::string&  title = ""  ) 
        : Atan2 ( name , title , a , b )
      {}
      /// constructor with one variable
      Atan2 ( const std::string& name  , 
              const std::string& title , 
              RooAbsReal&        a     ) 
        : Atan2 ( name , title , a , RooRealConstant::value ( 1.0 ) ) 
      {}
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Id
     *  Trivial variable: identical trnaformation 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2020-03-17
     */ 
    class Id final : public OneVar 
    {
      // ========================================================================
      ClassDef(Ostap::MoreRooFit::Id , 1 ) ;  // sum of RooAbsReal objects
      // ========================================================================
    public:
      // ========================================================================
      Id () = default ;
      /// constructor with a variable
      Id ( const std::string& name  , 
           const std::string& title , 
           RooAbsReal&        a     ) ;
      /// copy 
      Id ( const Id&   right       , 
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
    protected:
      // ======================================================================
// #if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
//       // ======================================================================
//       RooSpan<double> evaluateBatch 
//         ( std::size_t begin     , 
//           std::size_t batchSize ) const override
//       { return m_x.arg().evaluateBatch ( begin , batchSize  ) ; }
//       // ======================================================================
// #endif
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;    
      // ======================================================================
    }; //    
    // ========================================================================    
  } //                                   The end of namespace Ostap::MoreRooFit  
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_MOREROOFIT_H
// ============================================================================
