#ifndef OSTAP_PYFUNCS_H 
#define OSTAP_PYFUNCS_H 1
// ============================================================================
// Include files
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/IFuncs.h"
// ============================================================================
// Python
// ============================================================================
#include "Python.h"
// ============================================================================
namespace Ostap 
{
  // ==========================================================================
  namespace  Functions 
  {
    // ========================================================================
    /** @class PyFuncTree 
     *  Helper class to implemet "TTree-functions" in python
     *  @see Ostap::IFuncTree
     */
     class PyFuncTree :  public Ostap::IFuncTree 
    {
    public :
      // ======================================================================
      /** constructor
       *  @param self python object
       *  @param tree pointer to the tree
       */
      PyFuncTree ( PyObject*    self     , 
                   const TTree* tree = 0 );
      // ======================================================================
      /// destructor 
      virtual ~PyFuncTree() ;
      // ======================================================================
      /// the basic 
      double evaluate ( const TTree* tree = 0 ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the pointer to TTree
      const TTree* tree () const { return m_tree ; }
      // ======================================================================
    private:
      // ======================================================================
      /// potentially cached pointer to the tree 
      mutable const TTree*    m_tree ;
      // ======================================================================
    private :
      // ======================================================================
      // self reference for python instance 
      PyObject*         m_self ; // self reference for python instance 
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
      /** constructor
       *  @param self python object
       *  @param data pointer to the data
       */
      PyFuncData ( PyObject*         self     , 
                   const RooAbsData* data = 0 );
      // ======================================================================
      /// destructor 
      virtual ~PyFuncData () ;
      // ======================================================================
      /// the basic 
      double evaluate ( const RooAbsData* data = 0 ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the pointer to TTree
      const RooAbsData* data () const { return m_data ; }
      // ======================================================================
    private:
      // ======================================================================
      /// potentially cached pointer to the tree 
      mutable const RooAbsData* m_data ;
      // ======================================================================
    private :
      // ======================================================================
      // self reference for python instance 
      PyObject*           m_self ; // self-reference for python instance 
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
