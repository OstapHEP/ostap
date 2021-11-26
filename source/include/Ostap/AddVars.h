// ============================================================================
#ifndef OSTAP_ADDVARS_H 
#define OSTAP_ADDVARS_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <string>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/IFuncs.h"
// ============================================================================
// Forward declarations 
// ============================================================================
class TTree      ; // from ROOT 
class TH1        ; // from ROOT 
class TH2        ; // from ROOT  
class TH3        ; // from ROOT 
class RooDataSet ; // from RooFit 
class RooAbsReal ; // from RooFit 
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Functions
  {
    // ========================================================================
    /** add new variable to dataset
     *  @param  dataset input    dataset
     *  @param  name    variable name 
     *  @param  func    rule to  calculate new variable
     *  @return the added variable 
     */
    const RooAbsReal* add_var 
    ( RooDataSet&             dataset , 
      const std::string&      name    , 
      const Ostap::IFuncData& func    ) ;
    // ========================================================================    
    /** add new variable to dataset
     *  @param  dataset input dataset
     *  @param  name    variable name 
     *  @param  formula formula 
     *  @return the added variable 
     */
    const RooAbsReal* add_var 
    ( RooDataSet&             dataset , 
      const std::string&      name    , 
      const std::string&      formula ) ;
    // ========================================================================
    /** add new variable to dataset, sampled from 1D-histogram
     *  @param  dataset input    dataset
     *  @param  name    variable name 
     *  @param  histo   histogram to be sampled 
     *  @return the added variable 
     *  @see TH1::GetRandom 
     */
    const RooAbsReal* add_var 
    ( RooDataSet&             dataset , 
      const std::string&      name    , 
      const TH1&              histo   ) ;
    // ========================================================================
    /** add new variables to dataset, sampled from 2D-histogram
     *  @param  dataset input    dataset
     *  @param  namex   variable name 
     *  @param  namey   variable name 
     *  @param  histo   histogram to be sampled 
     *  @return the added variable 
     *  @see TH2::GetRandom2
     */
    const RooAbsReal* add_var 
    ( RooDataSet&             dataset , 
      const std::string&      namex   , 
      const std::string&      namey   , 
      const TH2&              histo   ) ;
    // ========================================================================
    /** add new variables to dataset, sampled from 3D-histogram
     *  @param  dataset input    dataset
     *  @param  namex   variable name 
     *  @param  namey   variable name 
     *  @param  namez   variable name 
     *  @param  histo   histogram to be sampled 
     *  @return the added variable 
     *  @see TH3::GetRandom3
     */
    const RooAbsReal* add_var 
    ( RooDataSet&             dataset , 
      const std::string&      namex   , 
      const std::string&      namey   , 
      const std::string&      namez   , 
      const TH3&              histo   ) ;
    // ========================================================================
  } //                                    The end of namespace Ostap::Functions 
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_ADDVARS_H
// ============================================================================
