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
      // ======================================================================
      /// default constructor 
      Exp () =  default ;
      // ======================================================================
      // copy 
      Exp ( const Exp&     right       , 
            const char*    newname = 0 ) ;
      // ======================================================================
      // destructor 
      virtual ~Exp () ;
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
    /** @class ScaleAndShift 
     *  \f[ f = a +  b c \f]
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
     *  @date 2019-11-21
     */
    class ScaleAndShift : public RooAbsReal 
    {
      // ======================================================================
      ClassDef(Ostap::MoreRooFit::ScaleAndShift, 1 ) ;  // Fraction 
      // ======================================================================
    public:
      // ======================================================================
      /// constructor with three variables 
      ScaleAndShift  
      ( const std::string& name  , 
        const std::string& title , 
        RooAbsReal&    a     , 
        RooAbsReal&    b     , 
        RooAbsReal&    c     ) ;
      // ======================================================================   
      /// constructor with three variables 
      ScaleAndShift  
      ( RooAbsReal&        a          , 
        RooAbsReal&        b          , 
        RooAbsReal&        c          ,
        const std::string& name  = "" , 
        const std::string& title = "" ) 
        : ScaleAndShift ( name ,  title , a , b , c ) 
      {} 
      // ======================================================================
      /// default constructor 
      ScaleAndShift   () =  default ;
      // ======================================================================
      // copy 
      ScaleAndShift  ( const ScaleAndShift& right       , 
                       const char*          newname = 0 ) ;
      // ======================================================================
      // destructor 
      virtual ~ScaleAndShift () ;
      // ======================================================================
      // clone method 
      ScaleAndShift * clone ( const char* newname ) const override ;
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;    
      // ======================================================================
    private:
      // ======================================================================
      RooRealProxy m_a ; // a 
      RooRealProxy m_b ; // b 
      RooRealProxy m_c ; // c 
      // ======================================================================
    };
    // ========================================================================
  } //                                   The end of namespace Ostap::MoreRooFit  
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_MOREROOFIT_H
// ============================================================================
