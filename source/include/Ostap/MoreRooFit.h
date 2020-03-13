// ============================================================================
#ifndef OSTAP_MOREROOFIT_H 
#define OSTAP_MOREROOFIT_H 1
// ============================================================================
// Include files
// ============================================================================
// ROOT/RooFit 
// ============================================================================
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
    // ==========================================================================
    /** @class Division
     *  Evaluate \f$ \frac{a}{b}\f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Division : public RooAbsReal
    {
      // ========================================================================
      ClassDef(Ostap::MoreRooFit::Division , 1 ) ;  // ratio of RooAbsReal objects
      // ========================================================================
    public:
      // ======================================================================
      Division  ( const std::string& name  , 
                  const std::string& title , 
                  RooAbsReal&        a     , 
                  RooAbsReal&        b     ) ;
      /// constructor with two variables 
      Division ( RooAbsReal&         a           , 
                 RooAbsReal&         b           ,
                 const std::string&  name  = ""  , 
                 const std::string&  title = ""  ) 
        : Division ( name , title , a , b )
      {}
      // ======================================================================
      /// default constructor 
      Division  () =  default ;
      // ======================================================================
      // copy 
      Division ( const Division& right       , 
                 const char*     newname = 0 ) ;
      // ======================================================================
      // destructor 
      virtual ~Division() ;
      // ======================================================================
      // clone method 
      Division* clone ( const char* newname ) const override ;
      // ====================================================================== 
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;    
      // ======================================================================
    protected :
      // ======================================================================
      /// a-variable
      RooRealProxy m_A {} ; // a-variable
      /// b-variable
      RooRealProxy m_B {} ; /// b-variable
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Fraction
     *  Evaluate \f$ \frac{a}{a+b}\f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Fraction : public Division
    {
      // ======================================================================
      ClassDef(Ostap::MoreRooFit::Fraction , 1 ) ;  // Fraction 
      // ======================================================================
    public:
      // ======================================================================
      /// constructor with two variables 
      Fraction  ( const std::string& name  , 
                  const std::string& title , 
                  RooAbsReal&        a     , 
                  RooAbsReal&        b     ) ;
      /// constructor with two variables 
      Fraction ( RooAbsReal&         a           , 
                 RooAbsReal&         b           ,
                 const std::string&  name  = ""  , 
                 const std::string&  title = ""  ) 
        : Fraction ( name , title , a , b )
      {}
      // ======================================================================
      /// default constructor 
      Fraction  () =  default ;
      // ======================================================================
      // copy 
      Fraction ( const Fraction& right       , 
                 const char*     newname = 0 ) ;
      // ======================================================================
      // destructor 
      virtual ~Fraction() ;
      // ======================================================================
      // clone method 
      Fraction* clone ( const char* newname ) const override ;
      // 
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;    
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Asymmetry
     *  Evaluate \f$ \frac{a-b}{a+b}\f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Asymmetry: public Division
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
                  const char*         newname = 0 ) ;
      // ======================================================================
      // destructor 
      virtual ~Asymmetry () ;
      // ======================================================================
      // clone method 
      Asymmetry * clone ( const char* newname ) const override ;
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;    
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Power
     *  Evaluate \f$ a^b \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Power: public Division 
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
              const char*    newname = 0 ) ;
      // ======================================================================
      // destructor 
      virtual ~Power() ;
      // ======================================================================
      // clone method 
      Power* clone ( const char* newname ) const override ;
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
    class Exp: public Division 
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
      /// default constructor 
      Exp () =  default ;
      // ======================================================================
      // copy 
      Exp ( const Exp&     right       , 
            const char*    newname = 0 ) 
        : Division ( right  , newname ) 
      {}
      // ======================================================================
      // destructor 
      virtual ~Exp () {} ;
      // ======================================================================
      // clone method 
      Exp* clone ( const char* newname ) const override ;
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
    class Log: public Division 
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
      /// default constructor 
      Log () =  default ;
      // ======================================================================
      // copy 
      Log ( const Log&     right       , 
            const char*    newname = 0 ) 
        : Division ( right  , newname ) 
      {}
      // ======================================================================
      // destructor 
      virtual ~Log () {} ;
      // ======================================================================
      // clone method 
      Log* clone ( const char* newname ) const override ;
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
    class Erf: public Division 
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
      /// default constructor 
      Erf () =  default ;
      // ======================================================================
      // copy 
      Erf ( const Erf&     right       , 
            const char*    newname = 0 ) 
        : Division ( right  , newname ) 
      {}
      // ======================================================================
      // destructor 
      virtual ~Erf () {} ;
      // ======================================================================
      // clone method 
      Erf* clone ( const char* newname ) const override ;
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
    class Gamma: public Division 
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
      /// default constructor 
      Gamma () =  default ;
      // ======================================================================
      // copy 
      Gamma ( const Gamma&   right       , 
              const char*    newname = 0 ) 
        : Division ( right  , newname ) 
      {}
      // ======================================================================
      // destructor 
      virtual ~Gamma() {} ;
      // ======================================================================
      // clone method 
      Gamma* clone ( const char* newname ) const override ;
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
    class LGamma: public Division 
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
      /// default constructor 
      LGamma () =  default ;
      // ======================================================================
      // copy 
      LGamma ( const LGamma&  right       , 
               const char*    newname = 0 ) 
        : Division ( right  , newname ) 
      {}
      // ======================================================================
      // destructor 
      virtual ~LGamma() {} ;
      // ======================================================================
      // clone method 
      LGamma* clone ( const char* newname ) const override ;
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
    class IGamma: public Division 
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
      /// default constructor 
      IGamma () =  default ;
      // ======================================================================
      // copy 
      IGamma ( const IGamma&  right       , 
               const char*    newname = 0 ) 
        : Division ( right  , newname ) 
      {}
      // ======================================================================
      // destructor 
      virtual ~IGamma() {} ;
      // ======================================================================
      // clone method 
      IGamma* clone ( const char* newname ) const override ;
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
    class Sin: public Division 
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
      /// default constructor 
      Sin  () =  default ;
      // ======================================================================
      // copy 
      Sin  ( const Sin&     right       , 
             const char*    newname = 0 ) 
        : Division ( right  , newname ) 
      {}
      // ======================================================================
      // destructor 
      virtual ~Sin () {} ;
      // ======================================================================
      // clone method 
      Sin* clone ( const char* newname ) const override ;
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
    class Cos: public Division 
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
      /// default constructor 
      Cos  () =  default ;
      // ======================================================================
      // copy 
      Cos  ( const Cos&     right       , 
             const char*    newname = 0 ) 
        : Division ( right  , newname ) 
      {}
      // ======================================================================
      // destructor 
      virtual ~Cos () {} ;
      // ======================================================================
      // clone method 
      Cos* clone ( const char* newname ) const override ;
      // ====================================================================== 
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;    
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Tan
     *  Evaluate \f$ \tan ab  \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Tan: public Division 
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
      /// default constructor 
      Tan  () =  default ;
      // ======================================================================
      // copy 
      Tan  ( const Tan&     right       , 
             const char*    newname = 0 ) 
        : Division ( right  , newname ) 
      {}
      // ======================================================================
      // destructor 
      virtual ~Tan () {} ;
      // ======================================================================
      // clone method 
      Tan* clone ( const char* newname ) const override ;
      // ====================================================================== 
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;    
      // ======================================================================
    }; //
    // ========================================================================
    /** @class Tanh
     *  Evaluate \f$ \tanh ab  \f$
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class Tanh: public Division 
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
      /// default constructor 
      Tanh () =  default ;
      // ======================================================================
      // copy 
      Tanh ( const Tanh&    right             , 
             const char*    newname = nullptr )
        : Division ( right  , newname ) 
      {}
      // ======================================================================
      // destructor 
      virtual ~Tanh () {} ;
      // ======================================================================
      // clone method 
      Tanh* clone ( const char* newname ) const override ;
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
    class Atan2: public Division 
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
      /// default constructor 
      Atan2 () =  default ;
      // ======================================================================
      // copy 
      Atan2 ( const Atan2&   right             , 
              const char*    newname = nullptr )
        : Division ( right  , newname ) 
      {}
      // ======================================================================
      // destructor 
      virtual ~Atan2 () {} ;
      // ======================================================================
      // clone method 
      Atan2* clone ( const char* newname ) const override ;
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
