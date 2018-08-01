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
#include "Ostap/StatusCode.h"
#include "Ostap/IFuncs.h"
// ============================================================================
// Forward declarations
// ============================================================================
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
      InvalidEntry          = 401 ,
    } ;  
    // ========================================================================
    typedef std::map   <std::string,std::string>             MAP   ;
    typedef std::vector<MAP>                                 MAPS  ;
    typedef std::vector<std::pair<std::string,std::string> > PAIRS ;
    // ========================================================================
    /** Add TMVA response to dataset 
     *  The  function add variables  "prefix+methos+suffix" that 
     *  are the responses of TMVA. TMNVA  configurtaion for methods 
     *  is read from the trained (xml) weight-files 
     *  @param data         (UPDATE) dataset 
     *  @param inputs       (INPUT) map  { varname : formula     }  
     *  @param weight_files (INPUT) map  { method  : weight_file }  
     *  @param prefix       (INPUT) the prefix for added variables 
     *  @param suffix       (INPUT) the suffix for added variables 
     *  @param aux          (INPUT) obligatory for the cuts method,
     *                              where it represents the efficiency cutoff
     */ 
    Ostap::StatusCode addResponse
    ( RooDataSet&        data         ,
      const MAP&         inputs       , 
      const MAP&         weight_files ,
      const std::string& prefix = ""  , 
      const std::string& suffix = ""  , 
      const double       aux    = 0.9 ) ;
    // ========================================================================
    /** Add TMVA response to dataset 
     *  The  function add variables  "prefix+methos+suffix" that 
     *  are the responses of TMVA. TMNVA  configurtaion for methods 
     *  is read from the trained (xml) weight-files 
     *  @param data         (UPDATE) dataset 
     *  @param inputs       (INPUT) [ (varnname,formula)   , ...  ] 
     *  @param weight_files (INPUT) [ (method,weight_file) , ...  ]
     *  @param prefix       (INPUT) the prefix for added variables 
     *  @param suffix       (INPUT) the suffix for added variables 
     *  @param aux          (INPUT) obligatory for the cuts method,
     *                              where it represents the efficiency cutoff
     */ 
    Ostap::StatusCode addResponse
    ( RooDataSet&        data         ,
      const PAIRS&       inputs       , 
      const PAIRS&       weight_files , 
      const std::string& prefix = ""  , 
      const std::string& suffix = ""  ,
      const double       aux    = 0.9 ) ;
    // ========================================================================
    /** Add TMVA response to dataset 
     *  The  function add variables  "prefix+methos+suffix" that 
     *  are the responses of TMVA. TMNVA  configurtaion for methods 
     *  is read from the trained (xml) weight-files 
     *  @param data         (UPDATE) dataset 
     *  @param inputs       (INPUT) [ (varnname,formula)   , ...  ] 
     *  @param weight_files (INPUT) map  { method  : weight_file }  
     *  @param prefix       (INPUT) the prefix for added variables 
     *  @param suffix       (INPUT) the suffix for added variables 
     *  @param aux          (INPUT) obligatory for the cuts method,
     *                              where it represents the efficiency cutoff
     */
    Ostap::StatusCode addResponse
    ( RooDataSet&        data         ,
      const PAIRS&       inputs       , 
      const MAP&         weight_files ,
      const std::string& prefix = ""  , 
      const std::string& suffix = ""  ,
      const double       aux    = 0.9 ) ;
    // ========================================================================
    /** Add TMVA response to dataset 
     *  The  function add variables  "prefix+methos+suffix" that 
     *  are the responses of TMVA. TMNVA  configurtaion for methods 
     *  is read from the trained (xml) weight-files 
     *  @param data         (UPDATE) dataset 
     *  @param inputs       (INPUT) map  { varname : formula     }  
     *  @param weight_files (INPUT) [ (method,weight_file) , ...  ]
     *  @param prefix       (INPUT) the prefix for added variables 
     *  @param suffix       (INPUT) the suffix for added variables 
     *  @param aux          (INPUT) obligatory for the cuts method,
     *                              where it represents the efficiency cutoff
     */
    Ostap::StatusCode addResponse
    ( RooDataSet&        data         ,
      const MAP&         inputs       , 
      const PAIRS&       weight_files ,
      const std::string& prefix = ""  , 
      const std::string& suffix = ""  ,
      const double       aux    = 0.9 ) ;
    // ========================================================================
    // Chopping 
    // ========================================================================
    Ostap::StatusCode addChoppingResponse 
    ( RooDataSet&          data                   ,
      RooAbsReal&          chopping               , // category function 
      RooCategory&         category               , // category variable 
      const unsigned short N                      , // number of categories 
      const MAP&           inputs                 , // mapping of input variables 
      const MAPS&          weight_files           ,
      const std::string&   prefix   = ""          , 
      const std::string&   suffix   = ""          ,
      const double         aux      = 0.9         ) ;
    // ========================================================================
  } //                                         The END of namespace Ostap::TMVA 
  // ==========================================================================
} //                                                 The END of namespace Ostap
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_TMVA_H
