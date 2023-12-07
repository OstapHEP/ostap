// ============================================================================
#ifndef OSTAP_GETWEIGHT_H 
#define OSTAP_GETWEIGHT_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <utility>
#include <string>
// ============================================================================
// Forward declarations 
// ============================================================================
class RooAbsReal ;
class RooAbsData ;
class RooDataSet ;
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Utils 
  {
    // ========================================================================
    /** get the name of the weight variable from data set (if and when possible)
     *  @param data (INPUT) dataset
     *  @return the name of weigth variable, if and when possible 
     */
    std::string       getWeight ( const RooAbsData* data ) ;
    // ========================================================================
    /** get the weight variable from data set (if and when possible)
     *  @param data (INPUT) dataset
     *  @return the weigth variable, if and when possible 
     */
    const RooAbsReal* getWeightVar ( const RooAbsData* data ) ;
    // ========================================================================
    /** weight errors are stored for the weighted data set?
     *  @param data (INPUT) dataset
     *  @return true if weight errors are stored
     *  The function checks the <code>StoreError</code> and 
     *   <code>StoreAsymError</code> attributes for the weight variable 
     *  @see RooAbsData
     *  @see RooAbsArg::getAttribute
     */
    bool              storeError ( const RooAbsData* data ) ;
    // ========================================================================
    /** weigth errors are stored for the weighted data set?
     *  @param data (INPUT) dataset
     *  @return true if weight errors are stored
     *  The function checks the <code>StoreAsymError</code> attribute for the weight variable 
     *  @see RooAbsData
     *  @see RooAbsArg::getAttribute
     */
    bool              storeAsymError ( const RooAbsData* data ) ;
    // ========================================================================
    /** make un unweighted dataset from weighted one 
     *  @param weighted_data   (INPUT) input dataset 
     #  @param weigth_var      (INPUT) provde (new) name of weight variable 
     *  @return unweighted dataset + name of weight variable, if and when possible
     */
    std::pair<RooDataSet*,std::string> 
    unweight 
    ( const RooAbsData*  weighted_data      ,
      const std::string& weight_var    = "" ) ;
  // ========================================================================
  } //                                        The end of namespace Ostap::Utils
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_GETWEIGHT_H
// ============================================================================
