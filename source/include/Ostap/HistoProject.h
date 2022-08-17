// ===========================================================================
#ifndef OSTAP_HISTOPROJECT_H 
#define OSTAP_HISTOPROJECT_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <limits>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatusCode.h"
#include "Ostap/DataFrame.h"
// ============================================================================
// Forward declarations 
// =============================================================================
class TH1       ;     // ROOT 
class TH2       ;     // ROOT 
class TH3       ;     // ROOT 
// =============================================================================
class RooAbsData ; // RooFit 
class RooAbsReal ; // RooFit 
// =============================================================================
namespace Ostap
{  
  // ==========================================================================
  /** @class HistoProject Ostap/HistoProject.h
   *  Helper class to project Dataset/DataFrame to histogram 
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date   2015-10-08
   */
  class HistoProject 
  {
  public:
    // ========================================================================
    /** make a projection of RooDataSet into the histogram 
     *  @param data  (INPUT)  input data 
     *  @param histo (UPDATE) histogram 
     *  @param expression (INPUT) expression
     *  @param selection  (INPUT) selection criteria/weight 
     *  @param first (INPUT) the first event to process 
     *  @param last  (INPUT) the last event to process 
     */
    static Ostap::StatusCode project
    ( const RooAbsData*   data            , 
      TH1*                histo           ,
      const std::string&  expression      ,
      const std::string&  selection  = "" ,
      const unsigned long first      = 0                                         ,
      const unsigned long last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into the histogram 
     *  @param data  (INPUT)  input data 
     *  @param histo (UPDATE) histogram 
     *  @param xexpression (INPUT) expression for x-axis 
     *  @param yexpression (INPUT) expression for y-axis 
     *  @param selection  (INPUT) selection criteria/weight 
     *  @param first (INPUT) the first event to process 
     *  @param last  (INPUT) the last event to process 
     */
    static Ostap::StatusCode project2
    ( const RooAbsData*   data            , 
      TH2*                histo           ,
      const std::string&  xexpression     ,
      const std::string&  yexpression     ,
      const std::string&  selection  = "" ,
      const unsigned long first      = 0                                         ,
      const unsigned long last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into the histogram 
     *  @param data  (INPUT)  input data 
     *  @param histo (UPDATE) histogram 
     *  @param xexpression (INPUT) expression for x-axis 
     *  @param yexpression (INPUT) expression for y-axis 
     *  @param zexpression (INPUT) expression for z-axis 
     *  @param selection  (INPUT) selection criteria/weight 
     *  @param first (INPUT) the first event to process 
     *  @param last  (INPUT) the last event to process 
     */
    static Ostap::StatusCode project3
    ( const RooAbsData*   data            , 
      TH3*                histo           ,
      const std::string&  xexpression     ,
      const std::string&  yexpression     ,
      const std::string&  zexpression     ,
      const std::string&  selection  = "" ,
      const unsigned long first      = 0                                         ,
      const unsigned long last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
  public: // RooDataSet 
    // ========================================================================
    /** make a projection of RooDataSet into the histogram 
     *  @param data  (INPUT)  input data 
     *  @param histo (UPDATE) histogram 
     *  @param expression (INPUT) expression
     *  @param selection  (INPUT) selection criteria/weight 
     *  @param first (INPUT) the first event to process 
     *  @param last  (INPUT) the last event to process 
     */
    static Ostap::StatusCode project
    ( const RooAbsData*   data            , 
      TH1*                histo           ,
      const RooAbsReal&   expression      ,
      const RooAbsReal*   selection  = 0  ,
      const unsigned long first      = 0                                         ,
      const unsigned long last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into the histogram 
     *  @param data  (INPUT)  input data 
     *  @param histo (UPDATE) histogram 
     *  @param xexpression (INPUT) expression for x-axis 
     *  @param yexpression (INPUT) expression for y-axis 
     *  @param selection  (INPUT) selection criteria/weight 
     *  @param first (INPUT) the first event to process 
     *  @param last  (INPUT) the last event to process 
     */
    static Ostap::StatusCode project2
    ( const RooAbsData*   data            , 
      TH2*                histo           ,
      const RooAbsReal&   xexpression     ,
      const RooAbsReal&   yexpression     ,
      const RooAbsReal*   selection  = 0  ,
      const unsigned long first      = 0                                         ,
      const unsigned long last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into the histogram 
     *  @param data  (INPUT)  input data 
     *  @param histo (UPDATE) histogram 
     *  @param xexpression (INPUT) expression for x-axis 
     *  @param yexpression (INPUT) expression for y-axis 
     *  @param zexpression (INPUT) expression for z-axis 
     *  @param selection  (INPUT) selection criteria/weight 
     *  @param first (INPUT) the first event to process 
     *  @param last  (INPUT) the last event to process 
     */
    static Ostap::StatusCode project3
    ( const RooAbsData*   data            , 
      TH3*                histo           ,
      const RooAbsReal&   xexpression     ,
      const RooAbsReal&   yexpression     ,
      const RooAbsReal&   zexpression     ,
      const RooAbsReal*   selection  = 0  ,
      const unsigned long first      = 0                                         ,
      const unsigned long last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
  public:  //   DataFrame 
    // ========================================================================
    /** make a projection of DataFrame into the histogram 
     *  @param data  (INPUT)  input data 
     *  @param histo (UPDATE) histogram 
     *  @param expression (INPUT) expression
     *  @param selection  (INPUT) selection criteria/weight 
     */
    static Ostap::StatusCode project
    ( FrameNode           data            , 
      TH1*                histo           ,
      const std::string&  expression      ,
      const std::string&  selection  = "" ) ;
    // ========================================================================
    /** make a projection of RooDataSet into the histogram 
     *  @param data  (INPUT)  input data 
     *  @param histo (UPDATE) histogram 
     *  @param xexpression (INPUT) expression for x-axis 
     *  @param yexpression (INPUT) expression for y-axis 
     *  @param selection  (INPUT) selection criteria/weight 
     */
    static Ostap::StatusCode project2
    ( FrameNode           data            , 
      TH2*                histo           ,
      const std::string&  xexpression     ,
      const std::string&  yexpression     ,
      const std::string&  selection  = "" );
    // ========================================================================
    /** make a projection of RooDataSet into the histogram 
     *  @param data        (INPUT)  input data 
     *  @param histo       (UPDATE) histogram 
     *  @param xexpression (INPUT)  expression for x-axis 
     *  @param yexpression (INPUT)  expression for y-axis 
     *  @param zexpression (INPUT)  expression for z-axis 
     *  @param selection   (INPUT)  selection criteria/weight 
     */
    static Ostap::StatusCode project3
    ( FrameNode           data            , 
      TH3*                histo           ,
      const std::string&  xexpression     ,
      const std::string&  yexpression     ,
      const std::string&  zexpression     ,
      const std::string&  selection  = "" ) ;
    // ========================================================================
  } ;
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_HPROJECT_H
// ============================================================================

