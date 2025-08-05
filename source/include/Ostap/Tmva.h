// ============================================================================
#ifndef OSTAP_TMVA_H 
#define OSTAP_TMVA_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL 
// ============================================================================
#include <map> 
#include <vector> 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Types.h"
#include "Ostap/StatusCode.h"
#include "Ostap/ProgressConf.h"
// ============================================================================
// Forward declarations
// ============================================================================
class TTree       ; // from ROOT 
class RooDataSet  ; // from RooFit 
class RooAbsReal  ; // from RooFit 
class RooCategory ; // from RooFit 
// ============================================================================
namespace Ostap 
{
  // ==========================================================================
  /** @namespace Ostap::Tmva Ostap/Tmva.h 
   *  collection of helper functions to deal with TMVA
   */
  namespace Tmva
  {
    // ========================================================================
    /** Disable scatter pltof form TMVA 
     *  Unfortunately thee recommended action in PyROOT has no effect
     *
     * <code>
     * <PlotVariables> Will not produce scatter plots ==> 
     * : |  The number of 5 input variables and 0 target values would require 10 two-dimensional
     * : |  histograms, which would occupy the computer's memory. Note that this
     * : |  suppression does not have any consequences for your analysis, other
     * : |  than not disposing of these scatter plots. You can modify the maximum
     * : |  number of input variables allowed to generate scatter plots in your
     * : |  script via the command line:
     * : |  "(TMVA::gConfig().GetVariablePlotting()).fMaxNumOfAllowedVariablesForScatterPlots = <some int>;"
     * </code>      
     *
     */
    void disable_scatter_plots () ;
    // ========================================================================
  }
  // ==========================================================================
  /** @class AddTMVA
   *  Helper class to ads the TTMA/Chopping response to TTree or RooDataSet 
   */
  class AddTMVA
  {
    // ======================================================================
  public:
    // ======================================================================
    typedef Ostap::Dict<std::string>                                  MAP   ;
    typedef std::vector<MAP>                                          MAPS  ;
    typedef Ostap::Names                                              NAMES ; 
    // ======================================================================
  public:
    // ======================================================================
    /// constructor with progress bar configuratios
    AddTMVA ( const Ostap::Utils::ProgressConf& progress = false ) ;
    // ======================================================================
  public: // TMVA response --> RooDataSet 
    // ======================================================================
    /** Add TMVA response to dataset 
     *  The  function add variables  <code>prefix+method+suffix</code> that 
     *  are the responses of TMVA. 
     *  - TMVA  configuration for all methods 
     *  is read from the trained (xml) weight-files 
     *  @param data         (UPDATE) dataset 
     *  @param inputs       (INPUT) map  { varname : formula     }  
     *  @param weight_files (INPUT) map  { method  : weight_file }  
     *  @param prefix       (INPUT) the prefix for added variables 
     *  @param suffix       (INPUT) the suffix for added variables 
     *  @param aux          (INPUT) obligatory for the cuts method
     *                              where it represents the efficiency cutoff
     */ 
    Ostap::StatusCode addResponse
    ( RooDataSet*        data                    ,
      const MAP&         inputs                  , 
      const MAP&         weight_files            ,
      const NAMES&       spectators    = NAMES() , 
      const std::string& options       = ""      ,
      const std::string& prefix        = ""      , 
      const std::string& suffix        = ""      , 
      const double       aux           =0.9      ) const ;
    // =========================================================================
  public: // TMVA response --> TTree 
    // ========================================================================
    /** Add TMVA response to TTree
     *  The  function add branches <code>prefix+method+suffix</code> that 
     *  are the responses of TMVA. 
     *  - TMVA  configuration for all methods 
     *  is read from the trained (xml) weight-files 
     *  @param tree         (UPDATE) input TTree
     *  @param inputs       (INPUT) map  { varname : formula     }  
     *  @param weight_files (INPUT) map  { method  : weight_file }  
     *  @param prefix       (INPUT) the prefix for added variables 
     *  @param suffix       (INPUT) the suffix for added variables 
     *  @param aux          (INPUT) obligatory for the cuts method
     *                              where it represents the efficiency cutoff
     */ 
    Ostap::StatusCode addResponse
    ( TTree*             tree                    ,
      const MAP&         inputs                  , 
      const MAP&         weight_files            ,
      const NAMES&       spectators    = NAMES() , 
      const std::string& options       = ""      ,
      const std::string& prefix        = ""      , 
      const std::string& suffix        = ""      , 
      const double       aux           = 0.9     ) const ;
    // ========================================================================
  public : // Chopping response --> RooDataSet 
    // ========================================================================
    /** Add TMVA/Chopping response to dataset 
     *  The  function add branches <code>prefix+method+suffix</code> that 
     *  are the responses of TMVA. 
     *  - TMVA  configuration for all methods 
     *  is read from the trained (xml) weight-files 
     *  @param data         (UPDATE) input dataset 
     *  @param chopping     (INPUT) chopping varibale 
     *  @param chopping     (INPUT) chopping category 
     *  @param N            (INPUT) number of categories 
     *  @param inputs       (INPUT) map  { varname : formula     }  
     *  @param weight_files (INPUT) list [ map  { method  : weight_file } ] 
     *  @param prefix       (INPUT) the prefix for added variables 
     *  @param suffix       (INPUT) the suffix for added variables 
     *  @param aux          (INPUT) obligatory for the cuts method
     *                              where it represents the efficiency cutoff
     */ 
    Ostap::StatusCode addChoppingResponse 
    ( RooDataSet*          data                    ,
      RooAbsReal&          chopping                , // category function 
      RooCategory&         category                , // category variable 
      const unsigned short N                       , // number of categories 
      const MAP&           inputs                  , // mapping of input variables 
      const MAPS&          weight_files            ,
      const NAMES&         spectators    = NAMES() , 
      const std::string&   options       = ""      ,
      const std::string&   prefix        = ""      , 
      const std::string&   suffix        = ""      ,
      const double         aux           = 0.9     ) const ;
    // ========================================================================
  public : // Chopping response --> TTree
    // ========================================================================
    /** Add TMVA/Chopping response to TTree
     *  The  function add branches <code>prefix+method+suffix</code> that 
     *  are the responses of TMVA. 
     *  - TMVA  configuration for all methods 
     *  is read from the trained (xml) weight-files 
     *  @param tree         (UPDATE) input tree 
     *  @param chopping     (INPUT) chopping variable/expression
     *  @param chopping     (INPUT) chopping category 
     *  @param N            (INPUT) number of categories 
     *  @param inputs       (INPUT) map  { varname : formula     }  
     *  @param weight_files (INPUT) list[ map  { method  : weight_file } ; 
     *  @param prefix       (INPUT) the prefix for added variables 
     *  @param suffix       (INPUT) the suffix for added variables 
     *  @param aux          (INPUT) obligatory for the cuts method
     *                              where it represents the efficiency cutoff
     */ 
    Ostap::StatusCode addChoppingResponse 
    ( TTree*               tree                    ,
      const std::string&   chopping                , // category function 
      const std::string&   category_name           , // category variable 
      const unsigned short N                       , // number of categories 
      const MAP&           inputs                  , // mapping of input variables 
      const MAPS&          weight_files            ,
      const NAMES&         spectators    = NAMES() ,       
      const std::string&   options       = ""      ,
      const std::string&   prefix        = ""      , 
      const std::string&   suffix        = ""      ,
      const double         aux           = 0.9     ) const ;
    // ========================================================================  
  public: 
    // ========================================================================
    /// congfiguration of the progress bar 
    inline const Ostap::Utils::ProgressConf& progress () const
    { return m_progress ; }
    // ========================================================================
  private :
    // ========================================================================
    /// congfiguration of the progress bar 
    Ostap::Utils::ProgressConf m_progress { false } ; 
    // ========================================================================
  } ; //                                    The END of class Ostap::AddResponse 
  // ==========================================================================
} //                                                 The END of namespace Ostap
// ============================================================================
#endif // OSTAP_TMVA_H
// ============================================================================
//                                                                      The END
// ============================================================================
