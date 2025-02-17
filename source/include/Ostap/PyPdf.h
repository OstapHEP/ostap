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
      ClassDefOverride(Ostap::Models::PyPdf, 5 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /** Standard constructor
       *  @param name      the name of PDF 
       *  @param title     the title  of PDF 
       *  @param variables all variables 
       */
      PyPdf
      ( const std::string& name        , 
        const std::string& title       ,
        const RooArgList&  variables   ) ;
      // =======================================================================
      /** Standard constructor
       *  @param name      the name of PDF 
       *  @param title     the title  of PDF 
       *  @param observables observables 
       *  @param parameters  parameters 
       */
      PyPdf
      ( const std::string& name        ,   
        const std::string& title       ,
        const RooArgList&  observables , 
        const RooArgList&  parameters  ) ;      
      // ======================================================================
      /// copy  constructor 
      PyPdf
      ( const PyPdf& right           ,
        const char*  name  = nullptr ) ;
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
      const RooArgList& varlist () const { return m_varlist  ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get a value for a variable with index 
      double value ( const unsigned short index ) const ;
      /// get a value for a variable with name 
      double value ( const          char*  name ) const ;
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
      inline const RooArgSet*   allDeps () const { return m_allDeps   ; } 
      inline       RooArgSet*  analDeps () const { return m_analDeps  ; }
      inline const char*      rangeName () const { return m_rangeName ; }
      inline Int_t              intCode () const { return m_intCode   ; }
      /// move the function from protected to public integrface 
      Bool_t match_args ( const RooArgSet& vars ) const ;
      /// move the function from protected to public integrface 
      Bool_t match_arg  ( const RooAbsArg& Var  ) const ;
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
    /** @class PyPdfLite Ostap/PyPdf.h
     *  `Light' version of PyPdf
     *  @see Ostap::Models::PyPDF
     *  @see ostap.fitting.pypdf.PyPDF2
     *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
     *  @date   2018-06-06
     */
    class PyPdfLite : public RooAbsPdf
    {
    public: 
      // ======================================================================
      ClassDefOverride(Ostap::Models::PyPdfLite, 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /** Standard constructor
       *  @param name      the name of PDF 
       *  @param title     the title  of PDF 
       *  @param function  callable function
       *  @param variables all variables 
       */
      PyPdfLite
      ( const char*       name      , 
        const char*       title     ,
        PyObject*         function  , 
        const RooArgList& variables ) ;
      // =======================================================================
      /// copy  constructor 
      PyPdfLite
      ( const PyPdfLite& right          ,
        const char*       name = nullptr ) ;
      /// virtual destructor 
      virtual ~PyPdfLite () ;
      /// clone method 
      PyPdfLite* clone ( const char* name ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      // default constructor (needed for serialization)
      PyPdfLite() {} // default constructor (needed for serialization)
      // ======================================================================
    public:
      // ======================================================================
      ///  get all variables in a form of the list
      const RooArgList&   varlist   () const { return m_varlist  ; }
      /** get the underlyaing function 
       *  @attention referenc odut is incremented!
       */
      const PyObject*     function  () const ; 
      // ======================================================================
      std::size_t         numrefs   () const ;  
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    private:
      // ======================================================================
      // python partner
      PyObject*    m_function  { nullptr } ; //! python partner
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
