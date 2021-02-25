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
    std::string getWeight ( const RooAbsData* data ) ;
    // ========================================================================
    /** make un unweighted dataset from weighted one 
     *  @param weighted_data   (INPUT) input dataset 
     *  @return unweighted dataset, if and when possible
     */
    std::pair<RooDataSet*,std::string> 
    unweight ( const RooAbsData*  weighted_data      ,
               std::string        weight_var    = "" ) ;
    // ========================================================================
  } //                                        The end of namespace Ostap::Utils
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_GETWEIGHT_H
// ============================================================================
