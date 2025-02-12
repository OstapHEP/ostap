#ifndef OSTAP_PYFUNCS_H 
#define OSTAP_PYFUNCS_H 1
// ============================================================================
// Include files
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/OstapPyROOT.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/IFuncs.h"
// ============================================================================
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT 
struct  _object ;
typedef _object PyObject ;
#endif 
// ============================================================================
namespace Ostap 
{
  // ==========================================================================
  /** @namespace Ostap::Functions Ostap/PyFuncs.h
   *  collection of special functions/functors 
   */
  namespace  Functions 
  {
    // ========================================================================
    /** @class PyFuncTree  Ostap/PyFuncs.h
     *  Helper class to implemet "TTree-functions" in python
     *  @see Ostap::IFuncTree
     */
    class PyFuncTree : public Ostap::IFuncTree 
    {
    public :
      // ======================================================================
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT // ==========================
      // ======================================================================
      /** constructor
       *  @param self python object
       *  @param tree pointer to the tree
       */
      PyFuncTree ( PyObject*    self = 0 , 
                   const TTree* tree = 0 );
      /// copy 
      PyFuncTree ( const PyFuncTree& right ) ;
      // ======================================================================
#else // ======================================================================
      // ======================================================================
      /** constructor
       *  @param self python object
       *  @param tree pointer to the tree
       */
      PyFuncTree ( const TTree* tree ) ;
      /// copy 
      PyFuncTree ( const PyFuncTree& right ) = default ; 
      /// default constructor 
      PyFuncTree () : PyFuncTree ( nullptr ) {}
      // ======================================================================
#endif // =====================================================================
      // ======================================================================
      /// destructor 
      virtual ~PyFuncTree() ;
      // ======================================================================
    public:  
      // ======================================================================
      PyFuncTree* clone ( const char* name = "" ) const override ; 
      // ======================================================================
    public:
      // ======================================================================
      /// the basic 
      double operator() ( const TTree* tree = 0 ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// function that needs to be redefiend in python 
      virtual double evaluate () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the pointer to TTree
      const TTree* tree () const { return m_tree ; }
      // ======================================================================
    private:
      // ======================================================================
      /// potentially cached pointer to the tree 
      mutable const TTree*    m_tree { nullptr } ;
      // ======================================================================
    private :
      // ======================================================================
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT // ==========================
      // ======================================================================
      // self reference for python instance 
      PyObject*               m_self { nullptr } ; // self reference for python instance
      // ======================================================================
#endif // =====================================================================
      // ======================================================================
    } ;
    // ========================================================================
    /** @class PyFuncData
     *  Helper class to implemet "RooAbsData-functions" in python
     *  @see Ostap::IFuncData
     */
     class PyFuncData :  public Ostap::IFuncData 
    {
    public :
      // ======================================================================
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT 
      // ======================================================================
      /** constructor
       *  @param self python object
       *  @param data pointer to the data
       */
      PyFuncData ( PyObject*         self = 0 , 
                   const RooAbsData* data = 0 );
      /// copy 
      PyFuncData ( const PyFuncData& right ) ;
      // ======================================================================
#else 
      // ======================================================================
      /** constructor
       *  @param self python object
       *  @param data pointer to the data
       */
      PyFuncData ( const RooAbsData* data ) ;
      PyFuncData () : PyFuncData ( nullptr ) {} ;
      /// copy 
      PyFuncData ( const PyFuncData& right ) = default ;
      // ======================================================================
#endif 
      // ======================================================================
      /// destructor 
      virtual ~PyFuncData () ;
      // ======================================================================
    public:  
      // ======================================================================
      PyFuncData* clone ( const char* name = "" ) const override ; 
      // ======================================================================
    public:
      // ======================================================================
      /// the basic 
      double operator() ( const RooAbsData* data = 0 ) const override ;
      // ======================================================================
      /// function that needs to be redefiend in python 
      virtual double evaluate () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the pointer to TTree
      const RooAbsData* data () const { return m_data ; }
      // ======================================================================
    private:
      // ======================================================================
      /// potentially cached pointer to the tree 
      mutable const RooAbsData* m_data { nullptr } ;
      // ======================================================================
    private :
      // ======================================================================
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT 
      // self reference for python instance 
      PyObject* m_self { nullptr } ; ; // self-reference for python instance 
#endif 
      // ======================================================================
     } ;
    // ========================================================================
  } //                                The END of the namespace Ostap::Functions
  // ==========================================================================
} //                                             The END of the namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_PYFUNCS_H
// ============================================================================
