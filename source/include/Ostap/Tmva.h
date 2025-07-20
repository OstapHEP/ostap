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
  /** @namespace Ostap::TMVA Ostap/Tmva.h 
   *  collection of helper functions to deal with TMVA
   */
  namespace TMVA 
  {
    // ========================================================================
    enum {
      InvalidInputVariables = 201 ,
      InvalidWeightFiles          ,
      InvalidChoppingWeightFiles  ,
      InvalidBookTMVA       = 301 , 
      InvalidDataSet              ,
      InvalidFormula              ,
      InvalidChoppingFormula      ,
      InvalidChoppingCategory     ,
      InvalidVariable             ,
      InvalidBranch               , 
      InvalidTree                 , 
      InvalidEntry          = 401 ,
    } ;  
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
    typedef Ostap::Dict<std::string>                                   MAP  ;
    typedef std::vector<MAP>                                           MAPS ;
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
    ( RooDataSet*        data          ,
      const MAP&         inputs        , 
      const MAP&         weight_files  ,
      const std::string& options = ""  ,
      const std::string& prefix  = ""  , 
      const std::string& suffix  = ""  , 
      const double       aux     = 0.9 ) const ;
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
    ( TTree*             tree          ,
      const MAP&         inputs        , 
      const MAP&         weight_files  ,
      const std::string& options = ""  ,
      const std::string& prefix  = ""  , 
      const std::string& suffix  = ""  , 
      const double       aux     = 0.9 ) const ;
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
    ( RooDataSet*          data                   ,
      RooAbsReal&          chopping               , // category function 
      RooCategory&         category               , // category variable 
      const unsigned short N                      , // number of categories 
      const MAP&           inputs                 , // mapping of input variables 
      const MAPS&          weight_files           ,
      const std::string&   options  = ""          ,
      const std::string&   prefix   = ""          , 
      const std::string&   suffix   = ""          ,
      const double         aux      = 0.9         ) const ;
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
    ( TTree*               tree                   ,
      const std::string&   chopping               , // category function 
      const std::string&   category_name          , // category variable 
      const unsigned short N                      , // number of categories 
      const MAP&           inputs                 , // mapping of input variables 
      const MAPS&          weight_files           ,
      const std::string&   options  = ""          ,
      const std::string&   prefix   = ""          , 
      const std::string&   suffix   = ""          ,
      const double         aux      = 0.9         ) const ;
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
