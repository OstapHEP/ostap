// ============================================================================
#ifndef OSTAP_IFUNCS_H 
#define OSTAP_IFUNCS_H 1
// ============================================================================
// Include files
// ============================================================================
// Forward declarations
// ============================================================================
class TTree       ; // From ROOT 
class RooAbsData  ; // From RooFit
// ============================================================================
namespace Ostap 
{
  // ==========================================================================
  /** @class IFuncTree Ostap/IFuncs.h
   *  Helper abstract class to evaluate certain TTree-functions 
   */
  class IFuncTree 
  {
  public :
    // ========================================================================
    /// evaluate the function from TTree 
    virtual double     operator () ( const TTree* tree ) const = 0 ;
    /// virtual destructor 
    virtual ~IFuncTree  () ;
    // ========================================================================
  };
  // ==========================================================================
  /** @class IFuncData Ostap/IFuncs.h
   *  Helper abstract class to evaluate certain RooAbsData-functions 
   */
  class IFuncData 
  {
  public :
    // ========================================================================
    /// evaluate the function from TTree 
    virtual double operator ()  ( const RooAbsData* tree ) const = 0 ;
    /// virtual destructor 
    virtual ~IFuncData () ;
    // ========================================================================
  };
  // ==========================================================================
} //                                                The END of  namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_IFUNCS_H
// ============================================================================
