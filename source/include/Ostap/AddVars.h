// ============================================================================
#ifndef OSTAP_ADDVARS_H 
#define OSTAP_ADDVARS_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <string>
#include <functional>
// ============================================================================\
// Ostap
// ============================================================================
#include "Ostap/ProgressConf.h"
// ============================================================================
// Forward declarations 
// ============================================================================
// ROOT 
// ============================================================================
class TTree      ; // from ROOT 
class TH1        ; // from ROOT 
class TH2        ; // from ROOT  
class TH3        ; // from ROOT 
class RooDataSet ; // from RooFit 
class RooAbsReal ; // from RooFit 
// ============================================================================
// Ostap 
// ============================================================================
namespace Ostap { class IFuncData ; }
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  class IFuncTree ;
  class IFuncData ;
  // =========================================================================
  namespace Utils 
  {
    // ========================================================================
    /// progress bar configrutaion 
    class ProgressConf ; // progress bar configuration 
    // ========================================================================
  }
  // ==========================================================================
  /** @class AddVar 
   *  Helper class to add  the variables into RooDataSet 
   *  @see RooDataSet 
   */
  class AddVars 
  {
    // ========================================================================
  public :
    // ========================================================================
    /// constructor with the progress flag
    AddVars ( const Ostap::Utils::ProgressConf& progress = false) ;
    // ========================================================================    
  public:
    // ========================================================================
    /** add new variable to dataset
     *  @param  dataset input    dataset
     *  @param  name    variable name 
     *  @param  func    rule to  calculate new variable
     *  @param progress configuration of the progress bar
     *  @return the added variable 
     */
    const RooAbsReal* add_var 
    ( RooDataSet&                       dataset  , 
      const std::string&                name     , 
      const Ostap::IFuncData&           func     ) ; 
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
    // Generic 1D function
    // ========================================================================
    /** add new variable to dataset, calculated from generic function 
     *  @param  dataset input    dataset
     *  @param  vname   variable name 
     *  @param  xname   variable name 
     *  @param  fun     the function 
     */
    const RooAbsReal* add_var 
    ( RooDataSet&                   dataset , 
      const std::string&            vname   , 
      const std::string&            xname   , 
      std::function<double(double)> fun     ) ;
    // ========================================================================
    /** add new variable to dataset, calculated from generic function 
     *  @param  dataset input    dataset
     *  @param  vname   variable name 
     *  @param  xname   variable name 
     *  @param  fun     the function 
     */
    template <class FUNCTION>
    inline const RooAbsReal* add_var 
    ( RooDataSet&                   dataset , 
      const std::string&            vname   , 
      FUNCTION                      fun     , 
      const std::string&            xname   ) 
    { return add_var ( dataset , vname , xname , std::cref ( fun ) ) ; }
    // ========================================================================
    // Generic 2D function
    // ========================================================================
    /** add new variable to dataset, calculated from generic function 
     *  @param  dataset input    dataset
     *  @param  vname   variable name 
     *  @param  xname   variable name 
     *  @param  yname   variable name 
     *  @param  fun     the function 
     */
    const RooAbsReal* add_var 
    ( RooDataSet&                          dataset , 
      const std::string&                   vname   , 
      const std::string&                   xname   , 
      const std::string&                   yname   , 
      std::function<double(double,double)> fun     ) ;
    // ========================================================================
    /** add new variable to dataset, calculated from generic function 
     *  @param  dataset input    dataset
     *  @param  vname   variable name 
     *  @param  xname   variable name 
     *  @param  yname   variable name 
     *  @param  fun     the function 
     */
    template <class FUNCTION>
    inline const RooAbsReal* add_var 
    ( RooDataSet&                   dataset , 
      const std::string&            vname   , 
      FUNCTION                      fun     , 
      const std::string&            xname   ,
      const std::string&            yname   ) 
    { return add_var ( dataset , vname , xname , yname , std::cref ( fun ) ) ; }
    // ========================================================================
    // Generic 3d function
    // ========================================================================
    /** add new variable to dataset, calculated from generic function 
     *  @param  dataset input    dataset
     *  @param  vname   variable name 
     *  @param  xname   variable name 
     *  @param  yname   variable name 
     *  @param  zname   variable name 
     *  @param  fun     the function 
     */
    const RooAbsReal* add_var 
    ( RooDataSet&                                 dataset , 
      const std::string&                          vname   , 
      const std::string&                          xname   , 
      const std::string&                          yname   , 
      const std::string&                          zname   , 
      std::function<double(double,double,double)> fun     ) ;
    // ========================================================================
    /** add new variable to dataset, calculated from generic function 
     *  @param  dataset input    dataset
     *  @param  vname   variable name 
     *  @param  xname   variable name 
     *  @param  yname   variable name 
     *  @param  zname   variable name 
     *  @param  fun     the function 
     */
    template <class FUNCTION>
    inline const RooAbsReal* add_var 
    ( RooDataSet&                   dataset , 
      const std::string&            vname   , 
      FUNCTION                      fun     , 
      const std::string&            xname   ,
      const std::string&            yname   ,
      const std::string&            zname   ) 
    { return add_var ( dataset , vname , xname , yname , zname , std::cref ( fun ) ) ; }
    // ========================================================================
  public: 
    // ========================================================================
    /// congfiguration of the progress bar 
    inline const Ostap::Utils::ProgressConf& progress () const
    { return m_progress ; }
    // ========================================================================
  private:
    // ========================================================================
    /// configuration of progress bar 
    Ostap::Utils::ProgressConf m_progress {} ;     // progess bar configuration
    // ========================================================================    
  } ; //                                        The end of class Ostap::AddVar 
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_ADDVARS_H
// ============================================================================
