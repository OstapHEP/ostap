// ============================================================================
#ifndef OSTAP_PYPDF_H 
#define OSTAP_PYPDF_H 1
// ============================================================================
// Include files
// ============================================================================
// RooFit 
// ============================================================================
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooAbsReal.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/OstapPyROOT.h"
#include "Ostap/PyCallable.h"
// ============================================================================
namespace Ostap 
{
  // ==========================================================================
  namespace Models 
  {
    // ========================================================================
    /** @class PyPdf PyPdf.h Ostap/PyPdf.h
     *  Helper intermediate base class to implement "purely-python" RooAbsPdf
     *  @see ostap.fitting.pypdf.PyPDF
     *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
     *  @date   2018-06-06
     */
    class PyPdf : public RooAbsPdf
    {
    public: 
      // ======================================================================
      ClassDefOverride(Ostap::Models::PyPdf, 2 ) ;
      // ======================================================================
    public:
      // ======================================================================
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT 
      // ======================================================================
      /** Standard constructor
       *  @param self python partner for this instance 
       *  @param name      the name of PDF 
       *  @param title     the title  of PDF 
       *  @param variables all variables 
       */
      PyPdf ( PyObject*         self      , 
              const char*       name      , 
              const char*       title     ,
              const RooArgList& variables );
      // ======================================================================
      /** Standard constructor
       *  @param self python partner for this instance 
       *  @param name      the name of PDF 
       *  @param variables all variables 
       *  @param title     the title  of PDF 
       */
      PyPdf ( PyObject*          self      , 
              const std::string& name      , 
              const RooArgList&  variables ,
              const std::string& title = "" ) 
        : PyPdf ( self         , 
                  name.c_str() , 
                  title.empty() ?  name.c_str () : title.c_str() , 
                  variables ) 
      {}
      // ======================================================================
      /** Standard constructor
       *  @param self python partner for this instance 
       *  @param variables all variables 
       *  @param name      the name of PDF 
       *  @param title     the title  of PDF 
       */
      PyPdf ( PyObject*          self      , 
              const RooArgList&  variables ,
              const std::string& name      , 
              const std::string& title = "" ) 
        : PyPdf ( self , name , variables , title ) 
      {}
      // ======================================================================
#else 
      // ======================================================================
      /** Standard constructor
       *  @param name      the name of PDF 
       *  @param title     the title  of PDF 
       *  @param variables all variables 
       */
      PyPdf ( const char*       name      , 
              const char*       title     ,
              const RooArgList& variables );
      // ======================================================================
      /** Standard constructor
       *  @param name      the name of PDF 
       *  @param variables all variables 
       *  @param title     the title  of PDF 
       */
      PyPdf ( const std::string& name       , 
              const RooArgList&  variables  , 
              const std::string& title = "" ) 
        : PyPdf ( name.c_str() , 
                  title.empty() ? name.c_str() : title.c_str() , 
                  variables ) 
      {} 
      // ======================================================================
      /** Standard constructor
       *  @param variables all variables 
       *  @param name      the name of PDF 
       *  @param title     the title  of PDF 
       */
      PyPdf ( const RooArgList&  variables  , 
              const std::string& name       ,
              const std::string& title = "" ) 
        : PyPdf ( name , variables , title ) 
      {}
      // ======================================================================
#endif
      // ======================================================================
      /// copy  constructor 
      PyPdf ( const PyPdf& right , const char* name = nullptr ) ;
      /// virtual destructor 
      virtual ~PyPdf() ;
      /// clone method 
      PyPdf* clone ( const char* name ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      // default constructor (needed for serialization)
      PyPdf() {} // default constructor (needed for serialization)
      // ======================================================================
    public:
      // ======================================================================
      ///  get all variables in a form of the list 
      const RooArgList&   proxylist () const { return m_varlist ; }
      ///  get all variables in a form of the list      
      /// const RooListProxy& variables () const { return m_varlist ; }
      const RooArgList&   variables () const { return m_varlist ; }
      ///  get all variables in a form of the list 
      // const RooListProxy& varlist   () const { return m_varlist ; }
      const RooArgList&   varlist   () const { return m_varlist ; }
      // ======================================================================
    public:
      // ======================================================================
      std::string myid() const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get a variable with index 
      double variable ( const unsigned short index ) const ;
      /// get a variable with name 
      double variable ( const          char*  name ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    public: // analytical integrals 
      // ======================================================================
      Int_t    getAnalyticalIntegral
      ( RooArgSet&     allVars      ,
        RooArgSet&     analVars     ,
        const char*    rangeName    ) const override;
      Double_t analyticalIntegral
      ( Int_t          code         ,
        const char*    rangeName    ) const override;
      // ======================================================================
    public: /// helper methods  for implementation of getAnalyticalIntegral 
      // ======================================================================
      inline const RooArgSet*   allDeps () const { return  m_allDeps  ; } 
      inline       RooArgSet*  analDeps () const { return m_analDeps  ; }
      inline const char*      rangeName () const { return m_rangeName ; }
      inline Int_t              intCode () const { return m_intCode   ; }
      /// move the function from protected to public integrface 
      Bool_t matchArgs ( const RooArgSet& refVars ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// helper function to be redefined in python  
      virtual int    get_analytical_integral () const ;
      /// helper function to be redefined in python  
      virtual double     analytical_integral () const ;
      // ======================================================================
    private:
      // ======================================================================  
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT
      // ======================================================================
      /// python's  "self"
      PyObject*    m_self      { nullptr } ; // python's  "self"
      // ======================================================================
#endif  
      // ======================================================================
      /// all variables as list of variables 
      RooListProxy m_varlist {} ; // all variables as list of variables 
      // ======================================================================  
    private: // helper fields for implementation of getAnalyticalIntegral 
      // ======================================================================  
      mutable const RooArgSet*  m_allDeps   { nullptr } ;
      mutable       RooArgSet*  m_analDeps  { nullptr } ;
      mutable const char*       m_rangeName { nullptr } ;
      mutable Int_t             m_intCode   { 0       } ; 
      // ======================================================================  
    } ;
    // ========================================================================
    /** @class PyPdf2 Ostap/PyPdf.h
     *  ``Light'' version of PyPdf
     *  @see Ostap::Models::PyPDF
     *  @see ostap.fitting.pypdf.PyPDF2
     *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
     *  @date   2018-06-06
     */
    class PyPdf2 : public RooAbsPdf
    {
    public: 
      // ======================================================================
      ClassDefOverride(Ostap::Models::PyPdf2, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /** Standard constructor
       *  @param name      the name of PDF 
       *  @param title     the title  of PDF 
       *  @param function  callable function
       *  @param variables all variables 
       */
      PyPdf2 ( const char*       name      , 
               const char*       title     ,
               PyObject*         function  , 
               const RooArgList& variables );
      // =======================================================================
      /** Standard constructor
       *  @param name      the name of PDF 
       *  @param function  callable function
       *  @param variables all variables 
       *  @param title     the title  of PDF 
       */
      PyPdf2 ( const std::string& name       , 
               PyObject*          function   ,
               const RooArgList&  variables  ,
               const std::string& title = "" )
        : PyPdf2 ( name.c_str() , 
                   title.empty() ? name.c_str() : title.c_str() , 
                   function     ,  
                   variables    )
      {}      
      // =======================================================================
      /** Standard constructor
       *  @param name      the name of PDF 
       *  @param function  callable function
       *  @param variables all variables 
       *  @param title     the title  of PDF 
       */
      PyPdf2 ( const std::string& name       , 
               const RooArgList&  variables  ,
               PyObject*          function   ,
               const std::string& title = "" )
        : PyPdf2 ( name , function , variables , title ) 
      {}
      /// copy  constructor 
      PyPdf2 ( const PyPdf2& right , const char* name = nullptr ) ;
      /// virtual destructor 
      virtual ~PyPdf2() ;
      /// clone method 
      PyPdf2* clone ( const char* name ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      // default constructor (needed for serialization)
      PyPdf2() {} // default constructor (needed for serialization)
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    private:
      // ======================================================================
      // python partner
      PyObject*    m_function  { nullptr } ; // python partner
      PyObject*    m_arguments { nullptr } ; // argument cache
      /// all variables as list of variables 
      RooListProxy m_varlist   {} ; // all variables as list of variables 
      // ======================================================================  
    } ;
    // ========================================================================
  } //                                   The end of the namespace Ostap::Models
  // ==========================================================================
} //                                             The end of the namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_PYPDF_H
// ============================================================================
