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
#include "RooSetProxy.h"
#include "RooAbsReal.h"
// ============================================================================
// ROOT
// ============================================================================
#include "TPySelector.h"
// ============================================================================
namespace Ostap 
{
  // ==========================================================================
  namespace Models 
  {
    // ========================================================================
    /** @class PyPdf PyPdf.h Ostap/PyPdf.h
     *  Helper base class to implement "purely-python"
     *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
     *  @date   2018-06-06
     */
    class PyPdf : public RooAbsPdf
    {
    public: 
      // ======================================================================
      ClassDef(Ostap::Models::PyPdf, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /** Standard constructor
       *  @param self python partner for this instance 
       *  @param name      the name of PDF 
       *  @param title     the title  of PDF 
       *  @param variables all variables 
       */
      PyPdf ( PyObject*               self      , 
              const char*             name      , 
              const char*             title     ,
              const RooAbsCollection& variables );
      /// copy  constructor 
      PyPdf ( const PyPdf& right , const  char* name = 0 ) ;
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
      const RooArgList& varlist () const { return m_varlist ; }
      ///  get all varianles in a form of the set  
      const RooArgSet&  varset  () const { return m_varset  ; }      
      // ======================================================================
    public:
      // ======================================================================
      // the actual evaluation of function
      Double_t evaluate() const override;
      // ======================================================================
    private:
      // ======================================================================  
      // python partner
      PyObject*    m_self    { nullptr  } ; // python partner 
      /// all variables as list of variables 
      RooListProxy m_varlist {} ; // all variables as list of variables 
      /// all variables as set  of variables 
      RooSetProxy  m_varset  {} ; // all variables as set  of variables 
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
