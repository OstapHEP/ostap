// ============================================================================
#ifndef OSTAP_PYVAR_H 
#define OSTAP_PYVAR_H 1
// ============================================================================
// Include files
// ============================================================================
// Python
// ============================================================================
#include "Python.h"
// ============================================================================
// ROOT/RooFit
// ============================================================================
#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooListProxy.h"
// ============================================================================
// Ostap
// ============================================================================
#include "PyCallable.h"
// ============================================================================
namespace Ostap 
{
  // ==========================================================================
  namespace Functions 
  {
    // ========================================================================
    /** @class PyVar Ostap/PyVar.h
     *  The attempt to create python analogue of RooFormulaVar
     *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
     *  @date   2019-04-28
     */
    class PyVar : public RooAbsReal
    {
      // ======================================================================
    public:
      // ======================================================================
      ClassDef ( Ostap::Functions::PyVar , 1 ) ;
      // ======================================================================
    public: 
      // ======================================================================
      /** Standard constructor
       *  @param self python object 
       *  @param name      the obkect name 
       *  @param title     the objkect title
       *  @param variables the list of variables 
       */
      PyVar (  PyObject*         self      , 
               const char*       name      , 
               const char*       title     ,
               const RooArgList& variables ) ;
      /// Copy constructor
      PyVar ( const PyVar& right , const char* name = 0 ) ;
      /// virtual destructor
      virtual ~PyVar() ;
      // ======================================================================
    public:
      // ======================================================================
      // default constructor (needed for serialization)
      PyVar () {} // default constructor (needed for serialization)
      // ======================================================================
    public:
      // ======================================================================
      /// clone method 
      PyVar* clone ( const char* name = 0 ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// the actual evaluation of the function 
      Double_t evaluate() const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// get a variable with index 
      double variable ( const unsigned short index ) const ;
      /// get a variable with name 
      double variable ( const          char*  name ) const ;
      // =====================================================================
      /// get all   parameters in a form of list
      const RooListProxy& variables () const { return m_variables    ; }
      const RooListProxy& params    () const { return   variables () ; }
      /// get all   parameters in a form of list
      const RooListProxy& varlist   () const { return   variables () ; }
      // ======================================================================
    private:
      // ======================================================================
      /// python's  "self"
      PyObject*    m_self      { nullptr } ; // python's  "self"
      /// the list of variables/parameters 
      RooListProxy m_variables {} ; // the list of variables/parameters 
      // ======================================================================      
    };
    // ========================================================================
    /** @class PyVar2 Ostap/PyVar.h
     *  ``Light'' version of Ostap::Functions::PyVar
     *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
     *  @date   2019-04-28
     */
    class PyVar2 : public RooAbsReal
    {
      // ======================================================================
    public:
      // ======================================================================
      ClassDef ( Ostap::Functions::PyVar2 , 1 ) ;
      // ======================================================================
    public: 
      // ======================================================================
      /** Standard constructor
       *  @param name      the object name 
       *  @param title     the object title
       *  @param function  python object 
       *  @param variables the list of variables 
       */
      PyVar2 ( const char*       name      , 
               const char*       title     ,
               PyObject*         function  , 
               const RooArgList& variables ) ;
      /// Copy constructor
      PyVar2 ( const PyVar2& right , const char* name = 0 ) ;
      /// virtual destructor
      virtual ~PyVar2() ;
      // ======================================================================
    public:
      // ======================================================================
      // default constructor (needed for serialization)
      PyVar2() {} // default constructor (needed for serialization)
      // ======================================================================
     public:
      // ======================================================================
      /// clone method 
      PyVar2* clone ( const char* name = 0 ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// the actual evaluation of the function 
      Double_t evaluate() const override ;
      // ======================================================================
    private:
      // ======================================================================
      /// python's  function
      PyObject*    m_function  { nullptr } ; // python's  "self"
      PyObject*    m_arguments { nullptr } ; // arguments cache 
      /// the list of variables/parameters 
      RooListProxy m_variables {} ; // the list of variables/parameters 
      // ======================================================================      
    };
    // ========================================================================
  } //                                    The end of namespace Ostap::Functions
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_PYVAR_H
// ============================================================================
