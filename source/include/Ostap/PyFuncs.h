#ifndef OSTAP_PYFUNCS_H 
#define OSTAP_PYFUNCS_H 1
// ============================================================================
// Include files
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/IFuncs.h"
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
      PyFuncData ( const RooAbsData* data ) ;
      PyFuncData () : PyFuncData ( nullptr ) {} ;
      /// copy 
      PyFuncData ( const PyFuncData& right ) = default ;
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
