// ============================================================================
#ifndef OSTAP_PYVAR_H 
#define OSTAP_PYVAR_H 1
// ============================================================================
// Include files
// ============================================================================
// ROOT/RooFit
// ============================================================================
#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooListProxy.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/OstapPyROOT.h"
#include "Ostap/PyCallable.h"
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
      ClassDefOverride ( Ostap::Functions::PyVar , 2 ) ;
      // ======================================================================
    public: 
      // ======================================================================
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT 
      // ======================================================================
      /** Standard constructor
       *  @param self python object 
       *  @param name      the object name 
       *  @param title     the object title
       *  @param variables the list of variables 
       */
      PyVar ( PyObject*         self      , 
              const char*       name      , 
              const char*       title     ,
              const RooArgList& variables ) ;
      // ======================================================================
      /** Standard constructor
       *  @param self python object 
       *  @param variables the list of variables 
       *  @param name      the object name 
       *  @param title     the objkect title
       */
      PyVar ( PyObject*         self        , 
              const RooArgList& variables   , 
              const std::string& name       , 
              const std::string& title = "" )
        : PyVar ( self , name.c_str() , 
                  title.empty () ? name.c_str() : title.c_str() , 
                  variables ) 
      {}
      // ======================================================================
      /** Standard constructor
       *  @param self python object 
       *  @param variables the list of variables 
       *  @param name      the object name 
       *  @param title     the objkect title
       */
      PyVar ( PyObject*         self        , 
              const std::string& name       , 
              const RooArgList& variables   , 
              const std::string& title = "" )
        : PyVar ( self , variables , name , title ) 
      {}
      // ======================================================================
#else 
      // ======================================================================
      /** Standard constructor
       *  @param name      the object name 
       *  @param title     the object title
       *  @param variables the list of variables 
       */
      PyVar ( const char*       name      , 
              const char*       title     ,
              const RooArgList& variables ) ;
      // ======================================================================
      /** Standard constructor
       *  @param name      the object name 
       *  @param variables the list of variables 
       *  @param title     the object title
       */
      PyVar ( const std::string& name       , 
              const RooArgList&  variables  , 
              const std::string& title = "" )
        : PyVar ( name.c_str() , 
                  title.empty() ? name.c_str() : title.c_str() , 
                  variables  ) 
      {}
      // =======================================================================
      /** Standard constructor
       *  @param variables the list of variables 
       *  @param name      the object name 
       *  @param title     the object title
       */
      PyVar ( const RooArgList&  variables  , 
              const std::string& name       , 
              const std::string& title = "" ) 
        : PyVar  ( name , variables , title ) 
      {}
      // ======================================================================
#endif 
      // ======================================================================
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
      const RooArgList& variables () const { return m_variables    ; }
      const RooArgList& params    () const { return   variables () ; }
      /// get all   parameters in a form of list
      const RooArgList& varlist   () const { return   variables () ; }
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
      ClassDefOverride ( Ostap::Functions::PyVar2 , 1 ) ;
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
      // ========================================================================
      /** Standard constructor
       *  @param name      the object name 
       *  @param function  python object 
       *  @param variables the list of variables 
       *  @param title     the object title
       */
      PyVar2 ( const std::string& name       , 
               PyObject*          function   , 
               const RooArgList&  variables  , 
               const std::string& title = "" )
        : PyVar2 ( name.c_str() , 
                   title.empty() ? name.c_str() : title.c_str() , 
                   function     , 
                   variables    ) 
      {}
      // ============================================================================
      /** Standard constructor
       *  @param name      the object name 
       *  @param variables the list of variables 
       *  @param function  python object 
       *  @param title     the object title
       */
      PyVar2 ( const std::string& name       , 
               const RooArgList&  variables  , 
               PyObject*          function   , 
               const std::string& title = "" )
        : PyVar2 ( name , function , variables , title ) 
      {}
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
