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
      ClassDefOverride ( Ostap::Functions::PyVar , 5 ) ;
      // ======================================================================
    public: 
      // ======================================================================
      /** Standard constructor
       *  @param name      the object name 
       *  @param title     the object title
       *  @param variables the list of variables 
       */
      PyVar
      ( const char*       name      , 
        const char*       title     ,
        const RooArgList& variables ) ;
      // ======================================================================
      /** Standard constructor
       *  @param name      the object name 
       *  @param variables the list of variables 
       *  @param title     the object title
       */
      PyVar
      ( const std::string& name       , 
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
      PyVar
      ( const RooArgList&  variables  , 
        const std::string& name       , 
        const std::string& title = "" ) 
        : PyVar  ( name , variables , title ) 
      {}
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
      double value ( const unsigned short index ) const ;
      /// get a variable with name 
      double value ( const          char*  name ) const ;
      // =====================================================================
      /// get all   parameters in a form of list
      const RooArgList&   varlist   () const { return m_varlist ; }
      // ======================================================================
      /// get a number of variables 
      std::size_t         nvars      () const { return m_varlist.size() ; }
      //// get values as a vector 
      std::vector<double> get_values () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// the list of variables/parameters 
      RooListProxy m_varlist {} ; // the list of variables/parameters 
      // ======================================================================      
    };
    // ========================================================================
    /** @class PyVarLite Ostap/PyVar.h
     *  ``Light'' version of Ostap::Functions::PyVar
     *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
     *  @date   2019-04-28
     */
    class PyVarLite : public RooAbsReal
    {
      // ======================================================================
    public:
      // ======================================================================
      ClassDefOverride ( Ostap::Functions::PyVarLite , 5 ) ;
      // ======================================================================
    public: 
      // ======================================================================
      /** Standard constructor
       *  @param name      the object name 
       *  @param title     the object title
       *  @param function  python object 
       *  @param variables the list of variables 
       */
      PyVarLite
      ( const char*       name      , 
        const char*       title     ,
        PyObject*         function  ,
        const RooArgList& variables ) ;
      // ========================================================================
      /// Copy constructor
      PyVarLite ( const PyVarLite& right , const char* name = 0 ) ;
      /// virtual destructor
      virtual ~PyVarLite() ;
      // ======================================================================
    public:
      // ======================================================================
      // default constructor (needed for serialization)
      PyVarLite() {} // default constructor (needed for serialization)
      // ======================================================================
    public:
      // ======================================================================
      /// clone method 
      PyVarLite* clone ( const char* name = 0 ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      ///  get all variables in a form of the list
      const RooArgList&   varlist    () const { return m_varlist  ; }
      /// get a number of variables 
      std::size_t         nvars      () const { return m_varlist.size() ; }
      //// get values as a vector 
      std::vector<double> get_values () const ;
      // ======================================================================
    public:
      // ======================================================================
      /** get the underlyaing function 
       *  @attention referenc odut is incremented!
       */
      const PyObject*     function  () const ; 
      // ======================================================================
      /// #of referenced for internap PyObject 
      std::size_t         numrefs   () const ;  
      // ======================================================================
    public:
      // ======================================================================
      /// the actual evaluation of the function 
      Double_t evaluate() const override ;
      // ======================================================================
    private:
      // ======================================================================
      /// python's  function
      PyObject*    m_function  { nullptr } ; //! python partner 
      /// the list of variables/parameters 
      RooListProxy m_varlist {} ; // the list of variables/parameters 
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
