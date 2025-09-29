// ===========================================================================
#ifndef OSTAP_PROJECT_H 
#define OSTAP_PROJECT_H 1
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
#include "Ostap/StatVar.h"
// ============================================================================
// Forward declarations 
// =============================================================================
// ROOT&RooFit 
// =============================================================================
class TH1        ;  // ROOT 
class TH2        ;  // ROOT 
class TH3        ;  // ROOT 
class TTree      ;  // ROOT
class TProfile   ;  // ROOT 
class TProfile2D ;  // ROOT 
class RooAbsData ;  // ROOT/RooFit 
class RooAbsReal ;  // ROOT/RooFit 
// =============================================================================
// Ostap
// =============================================================================
namespace Ostap
{
  // ===========================================================================
  namespace Math
  {
    // =========================================================================
    class ChebyshevSum ; 
    class LegendreSum  ;
    class LegendreSum2 ;
    class LegendreSum3 ;
    class LegendreSum4 ;
    class Bernstein    ;
    class Bernstein2D  ;
    class Bernstein3D  ;
    class ECDF         ;
    class WECDF        ;    
    // =========================================================================
  }
  // ===========================================================================
}
// =============================================================================
/** @file Ostap/HistoProject 
 *  Collecton of useful counters for "projection" of 
 *  data sources (TTree/RooAbsData) into data accumulators:
 * - counters
 * - histograms 
 * - paramterizations 
 * @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 * @data 2025-06-19
 */
namespace Ostap
{  
  // ===========================================================================
  /* @class Project
   * helper class to "project" TTree/RooabsData into 
   * - histograms
   * - polynoimals, 
   * - ...
   */
  class Project : public StatVar
  {
    // ========================================================================
  public :
    // ========================================================================
    /// construtctor with the progress flag
    Project ( const Ostap::Utils::ProgressConf& progress = false) ;
    // ========================================================================    
  public: // 1D histo 
    // ========================================================================
    /** Project data in the 1D-ihstoram
     *  @param data       (INPUT)  input data 
     *  @param histo      (UPDATE) the !D-histogram 
     *  @param expression (INPUT)  expression 
     *  @param selection  (INPUT)  selection/weight 
     *  @param first      (INOPUT) the first event to process (inclusive) 
     *  @param last       (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     *  @see   TH1::Fill
     */
    Ostap::StatusCode project1
    ( TTree*                  data            ,
      TH1*                    histo           ,
      const std::string&      expression      ,
      const std::string&      selection  = "" ,
      const Ostap::EventIndex first      = Ostap::FirstEvent ,
      const Ostap::EventIndex last       = Ostap::LastEvent  ) const ;
    // ==========================================================================
    /** Project data in the 1D-ihstoram
     *  @param data       (INPUT)  input data 
     *  @param histo      (UPDATE) the !D-histogram 
     *  @param expression (INPUT)  expression 
     *  @param selection  (INPUT)  selection/weight 
     *  @param cut_range  (INPUT)  cut-range 
     *  @param first      (INOPUT) the first event to process (inclusive) 
     *  @param last       (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     *  @see   TH1::Fill
     */
    Ostap::StatusCode project1
    ( const RooAbsData*       data            ,
      TH1*                    histo           ,
      const std::string&      expression      ,
      const std::string&      selection  = "" ,
      const std::string&      cut_range  = "" , 
      const Ostap::EventIndex first      = Ostap::FirstEvent ,
      const Ostap::EventIndex last       = Ostap::LastEvent  ) const ;
    // ========================================================================    
  public: // 2D histo 
    // ========================================================================    
    /** Project data in the 2D-ihstoram
     *  @param data        (INPUT)  input data 
     *  @param histo       (UPDATE) the !D-histogram 
     *  @param expression1 (INPUT)  expression1 
     *  @param expression2 (INPUT)  expression2 
     *  @param selection   (INPUT)  selection/weight  
     *  @param first       (INOPUT) the first event to process (inclusive) 
     *  @param last        (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     *  @see   TH1::Fill
     */
    Ostap::StatusCode project2
    ( TTree*                  data            ,
      TH2*                    histo           ,
      const std::string&      expression1     ,
      const std::string&      expression2     ,
      const std::string&      selection  = "" ,
      const Ostap::EventIndex first      = Ostap::FirstEvent ,
      const Ostap::EventIndex last       = Ostap::LastEvent  ) const ;      
    // =======================================================================
    /** Project data in the 2D-ihstoram
     *  @param data        (INPUT)  input data 
     *  @param histo       (UPDATE) the !D-histogram 
     *  @param expression1 (INPUT)  expression1 
     *  @param expression2 (INPUT)  expression2 
     *  @param selection   (INPUT)  selection/weight  
     *  @param first       (INOPUT) the first event to process (inclusive) 
     *  @param last        (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     *  @see   TH1::Fill
     */
    Ostap::StatusCode project2
    ( const RooAbsData*       data            ,
      TH2*                    histo           ,
      const std::string&      expression1     ,
      const std::string&      expression2     ,
      const std::string&      selection  = "" ,
      const std::string&      cut_range  = "" ,  
      const Ostap::EventIndex first      = Ostap::FirstEvent ,
      const Ostap::EventIndex last       = Ostap::LastEvent  ) const ;    
    // =======================================================================
  public: // 1D profile    
    // =======================================================================
    /** Project data in the 1D-profile 
     *  @param data        (INPUT)  input data 
     *  @param histo       (UPDATE) the 1D-profile 
     *  @param expression1 (INPUT)  expression1 
     *  @param expression2 (INPUT)  expression2 
     *  @param selection   (INPUT)  selection/weight  
     *  @param first       (INOPUT) the first event to process (inclusive) 
     *  @param last        (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     *  @see   TH1::Fill
     */
    Ostap::StatusCode project2
    ( TTree*                  data            ,
      TProfile*               histo           ,
      const std::string&      expression1     ,
      const std::string&      expression2     ,
      const std::string&      selection  = "" ,
      const Ostap::EventIndex first      = Ostap::FirstEvent ,
      const Ostap::EventIndex last       = Ostap::LastEvent  ) const ;
    // =======================================================================
    /** Project data in the 1D-historam
     *  @param data        (INPUT)  input data 
     *  @param histo       (UPDATE) the 1D-profile 
     *  @param expression1 (INPUT)  expression1 
     *  @param expression2 (INPUT)  expression2 
     *  @param selection   (INPUT)  selection/weight  
     *  @param first       (INOPUT) the first event to process (inclusive) 
     *  @param last        (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     *  @see   TH1::Fill
     */
    Ostap::StatusCode project2
    ( const RooAbsData*       data            ,
      TProfile*               histo           ,
      const std::string&      expression1     ,
      const std::string&      expression2     ,
      const std::string&      selection  = "" ,
      const std::string&      cut_range  = "" ,  
      const Ostap::EventIndex first      = Ostap::FirstEvent ,
      const Ostap::EventIndex last       = Ostap::LastEvent  ) const ;        
    // =======================================================================
  public: // 3D histo 
    // ========================================================================    
    /** Project data in the 3D-ihstoram
     *  @param data        (INPUT)  input data 
     *  @param histo       (UPDATE) the !D-histogram 
     *  @param expression1 (INPUT)  expression1 
     *  @param expression2 (INPUT)  expression2 
     *  @param expression3 (INPUT)  expression3
     *  @param selection   (INPUT)  selection/weight  
     *  @param first       (INOPUT) the first event to process (inclusive) 
     *  @param last        (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     *  @see   TH1::Fill
     */
    Ostap::StatusCode project3
    ( TTree*                  data            ,
      TH3*                    histo           ,
      const std::string&      expression1     ,
      const std::string&      expression2     ,
      const std::string&      expression3     ,
      const std::string&      selection  = "" ,
      const Ostap::EventIndex first      = Ostap::FirstEvent ,
      const Ostap::EventIndex last       = Ostap::LastEvent  ) const ;      
    // =======================================================================
    /** Project data in the 2D-ihstoram
     *  @param data        (INPUT)  input data 
     *  @param histo       (UPDATE) the !D-histogram 
     *  @param expression1 (INPUT)  expression1 
     *  @param expression2 (INPUT)  expression2 
     *  @param expression3 (INPUT)  expression3 
     *  @param selection   (INPUT)  selection/weight  
     *  @param first       (INOPUT) the first event to process (inclusive) 
     *  @param last        (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     *  @see   TH1::Fill
     */
    Ostap::StatusCode project3
    ( const RooAbsData*       data            ,
      TH3*                    histo           ,
      const std::string&      expression1     ,
      const std::string&      expression2     ,
      const std::string&      expression3     ,
      const std::string&      selection  = "" ,
      const std::string&      cut_range  = "" ,  
      const Ostap::EventIndex first      = Ostap::FirstEvent ,
      const Ostap::EventIndex last       = Ostap::LastEvent  ) const ;    
    // =======================================================================
  public: // 2D profile 
    // ========================================================================    
    /** Project data in the 2D-profile 
     *  @param data        (INPUT)  input data 
     *  @param histo       (UPDATE) the 2D-profile 
     *  @param expression1 (INPUT)  expression1 
     *  @param expression2 (INPUT)  expression2 
     *  @param expression3 (INPUT)  expression3
     *  @param selection   (INPUT)  selection/weight  
     *  @param first       (INOPUT) the first event to process (inclusive) 
     *  @param last        (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     *  @see   TH1::Fill
     */
    Ostap::StatusCode project3
    ( TTree*                  data            ,
      TProfile2D*             histo           ,
      const std::string&      expression1     ,
      const std::string&      expression2     ,
      const std::string&      expression3     ,
      const std::string&      selection  = "" ,
      const Ostap::EventIndex first      = Ostap::FirstEvent ,
      const Ostap::EventIndex last       = Ostap::LastEvent  ) const ;      
    // =======================================================================
    /** Project data in the 2D-profile 
     *  @param data        (INPUT)  input data 
     *  @param histo       (UPDATE) the 2D-profile 
     *  @param expression1 (INPUT)  expression1 
     *  @param expression2 (INPUT)  expression2 
     *  @param expression3 (INPUT)  expression3 
     *  @param selection   (INPUT)  selection/weight  
     *  @param first       (INOPUT) the first event to process (inclusive) 
     *  @param last        (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     *  @see   TH1::Fill
     */
    Ostap::StatusCode project3
    ( const RooAbsData*       data            ,
      TProfile2D*             histo           ,
      const std::string&      expression1     ,
      const std::string&      expression2     ,
      const std::string&      expression3     ,
      const std::string&      selection  = "" ,
      const std::string&      cut_range  = "" ,  
      const Ostap::EventIndex first      = Ostap::FirstEvent ,
      const Ostap::EventIndex last       = Ostap::LastEvent  ) const ;    
    // =======================================================================
  public : // ECDF & WECDD 
    // ======================================================================
    /** Project data in ECDS/WECDF 
     *  @param data       (INPUT)  input data 
     *  @param ecdf       (UPDATE) empirical cumulatiove distribution 
     *  @param expression (INPUT)  expression 
     *  @param selection  (INPUT)  selection (as boolean)
     *  @param first      (INOPUT) the first event to process (inclusive) 
     *  @param last       (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     */
    Ostap::StatusCode project1
    ( TTree*                    data            ,
      Ostap::Math::ECDF&        ecdf            ,
      const std::string&        expression      ,
      const std::string&        selection  = "" ,
      const Ostap::EventIndex   first      = Ostap::FirstEvent ,
      const Ostap::EventIndex   last       = Ostap::LastEvent  , 
      const Ostap::DataType     xmin       = Ostap::MinValue   ,
      const Ostap::DataType     xmax       = Ostap::MaxValue   ) const ;
    // =========================================================================
    /** Project data in ECDS/WECDF 
     *  @param data       (INPUT)  input data 
     *  @param ecdf       (UPDATE) empirical cumulatiove distribution 
     *  @param expression (INPUT)  expression 
     *  @param selection  (INPUT)  selection (as weight) 
     *  @param first      (INOPUT) the first event to process (inclusive) 
     *  @param last       (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     */
    Ostap::StatusCode project1
    ( TTree*                    data            ,
      Ostap::Math::WECDF&       ecdf            ,
      const std::string&        expression      ,
      const std::string&        selection  = "" ,
      const Ostap::EventIndex   first      = Ostap::FirstEvent ,
      const Ostap::EventIndex   last       = Ostap::LastEvent  , 
      const Ostap::DataType     xmin       = Ostap::MinValue   ,
      const Ostap::DataType     xmax       = Ostap::MaxValue   ) const ;
    // =========================================================================
    /** Project data in ECDS/WECDF 
     *  @param data       (INPUT)  input data 
     *  @param ecdf       (UPDATE) empirical cumulatiove distribution 
     *  @param expression (INPUT)  expression 
     *  @param selection  (INPUT)  selection (as weight) 
     *  @param first      (INOPUT) the first event to process (inclusive) 
     *  @param last       (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     */
    Ostap::StatusCode project1
    ( const RooAbsData*         data            ,
      Ostap::Math::WECDF&       ecdf            ,
      const std::string&        expression      ,
      const std::string&        selection  = "" ,
      const std::string&        cut_range  = "" , 
      const Ostap::EventIndex   first      = Ostap::FirstEvent ,
      const Ostap::EventIndex   last       = Ostap::LastEvent  , 
      const Ostap::DataType     xmin       = Ostap::MinValue   ,
      const Ostap::DataType     xmax       = Ostap::MaxValue   ) const ;
    // =======================================================================
  public : // 1D-Chebyshev polynomial: on-flight parameterisation 
    // ======================================================================
    /** Project data in 1D-polyoominal: on-flight parameterisation 
     *  @param data       (INPUT)  input data 
     *  @param poly       (UPDATE) the polynoimal 
     *  @param expression (INPUT)  expression 
     *  @param selection  (INPUT)  selection (as boolean)
     *  @param first      (INOPUT) the first event to process (inclusive) 
     *  @param last       (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     */
    Ostap::StatusCode project1
    ( TTree*                     data            ,
      Ostap::Math::ChebyshevSum& poly             ,
      const std::string&         expression      ,
      const std::string&         selection  = "" ,
      const Ostap::EventIndex    first      = Ostap::FirstEvent ,
      const Ostap::EventIndex    last       = Ostap::LastEvent  ) const ; 
    // =========================================================================
    /** Project data in 1D-polyoominal: on-flight parameterisation 
     *  @param data       (INPUT)  input data 
     *  @param poly       (UPDATE) the polynoimal 
     *  @param expression (INPUT)  expression 
     *  @param selection  (INPUT)  selection (as boolean)
     *  @param first      (INOPUT) the first event to process (inclusive) 
     *  @param last       (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     */
    Ostap::StatusCode project1
    ( const RooAbsData*          data            ,
      Ostap::Math::ChebyshevSum& poly            ,
      const std::string&         expression      ,
      const std::string&         selection  = "" ,
      const std::string&         cut_range  = "" , 
      const Ostap::EventIndex    first      = Ostap::FirstEvent ,
      const Ostap::EventIndex    last       = Ostap::LastEvent  ) const ; 
    // =========================================================================
  public : // 1D-Legendre polynomial: on-flight parameterisation 
    // ======================================================================
    /** Project data in 1D-polyoominal: on-flight parameterisation 
     *  @param data       (INPUT)  input data 
     *  @param poly       (UPDATE) the polynoimal 
     *  @param expression (INPUT)  expression 
     *  @param selection  (INPUT)  selection (as boolean)
     *  @param first      (INOPUT) the first event to process (inclusive) 
     *  @param last       (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     */
    Ostap::StatusCode project1
    ( TTree*                     data            ,
      Ostap::Math::LegendreSum&  poly             ,
      const std::string&         expression      ,
      const std::string&         selection  = "" ,
      const Ostap::EventIndex    first      = Ostap::FirstEvent ,
      const Ostap::EventIndex    last       = Ostap::LastEvent  ) const ; 
    // =========================================================================
    /** Project data in 1D-polyoominal: on-flight parameterisation 
     *  @param data       (INPUT)  input data 
     *  @param poly       (UPDATE) the polynoimal 
     *  @param expression (INPUT)  expression 
     *  @param selection  (INPUT)  selection (as boolean)
     *  @param first      (INOPUT) the first event to process (inclusive) 
     *  @param last       (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     */
    Ostap::StatusCode project1
    ( const RooAbsData*          data            ,
      Ostap::Math::LegendreSum&  poly            ,
      const std::string&         expression      ,
      const std::string&         selection  = "" ,
      const std::string&         cut_range  = "" , 
      const Ostap::EventIndex    first      = Ostap::FirstEvent ,
      const Ostap::EventIndex    last       = Ostap::LastEvent  ) const ; 
    // =========================================================================
  public : // 1D-Bernstein polynomial: on-flight parameterisation 
    // ======================================================================
    /** Project data in 1D-polyoominal: on-flight parameterisation 
     *  @param data       (INPUT)  input data 
     *  @param poly       (UPDATE) the polynoimal 
     *  @param expression (INPUT)  expression 
     *  @param selection  (INPUT)  selection (as boolean)
     *  @param first      (INOPUT) the first event to process (inclusive) 
     *  @param last       (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     */
    Ostap::StatusCode project1
    ( TTree*                     data            ,
      Ostap::Math::Bernstein&    poly             ,
      const std::string&         expression      ,
      const std::string&         selection  = "" ,
      const Ostap::EventIndex    first      = Ostap::FirstEvent ,
      const Ostap::EventIndex    last       = Ostap::LastEvent  ) const ; 
    // =========================================================================
    /** Project data in 1D-polyoominal: on-flight parameterisation 
     *  @param data       (INPUT)  input data 
     *  @param poly       (UPDATE) the polynoimal 
     *  @param expression (INPUT)  expression 
     *  @param selection  (INPUT)  selection (as boolean)
     *  @param first      (INOPUT) the first event to process (inclusive) 
     *  @param last       (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     */
    Ostap::StatusCode project1
    ( const RooAbsData*          data            ,
      Ostap::Math::Bernstein&    poly            ,
      const std::string&         expression      ,
      const std::string&         selection  = "" ,
      const std::string&         cut_range  = "" , 
      const Ostap::EventIndex    first      = Ostap::FirstEvent ,
      const Ostap::EventIndex    last       = Ostap::LastEvent  ) const ; 
    // =========================================================================
  public : // 2D-Bernstein polynomial: on-flight parameterisation 
    // ======================================================================
    /** Project data in 2D-polyoominal: on-flight parameterisation 
     *  @param data        (INPUT)  input data 
     *  @param poly        (UPDATE) the polynoimal 
     *  @param expression1 (INPUT)  expression1 
     *  @param expression2 (INPUT)  expression2 
     *  @param selection   (INPUT)  selection (as boolean)
     *  @param first       (INOPUT) the first event to process (inclusive) 
     *  @param last        (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     */
    Ostap::StatusCode project2
    ( TTree*                     data            ,
      Ostap::Math::Bernstein2D&  poly            ,
      const std::string&         expression1     ,
      const std::string&         expression2     ,
      const std::string&         selection  = "" ,
      const Ostap::EventIndex    first      = Ostap::FirstEvent ,
      const Ostap::EventIndex    last       = Ostap::LastEvent  ) const ; 
    // =========================================================================
    /** Project data in 2D-polyoominal: on-flight parameterisation 
     *  @param data        (INPUT)  input data 
     *  @param poly        (UPDATE) the polynoimal 
     *  @param expression1 (INPUT)  expression1 
     *  @param expression2 (INPUT)  expression2 
     *  @param selection   (INPUT)  selection (as boolean)
     *  @param first       (INOPUT) the first event to process (inclusive) 
     *  @param last        (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     */
    Ostap::StatusCode project2
    ( const RooAbsData*          data            ,
      Ostap::Math::Bernstein2D&  poly            ,
      const std::string&         expression1     ,
      const std::string&         expression2     ,
      const std::string&         selection  = "" ,
      const std::string&         cut_range  = "" , 
      const Ostap::EventIndex    first      = Ostap::FirstEvent ,
      const Ostap::EventIndex    last       = Ostap::LastEvent  ) const ; 
    // =========================================================================
  public : // 2D-Legendre polynomial: on-flight parameterisation 
    // ======================================================================
    /** Project data in 2D-polyoominal: on-flight parameterisation 
     *  @param data        (INPUT)  input data 
     *  @param poly        (UPDATE) the polynoimal 
     *  @param expression1 (INPUT)  expression1 
     *  @param expression2 (INPUT)  expression2 
     *  @param selection   (INPUT)  selection (as boolean)
     *  @param first       (INOPUT) the first event to process (inclusive) 
     *  @param last        (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     */
    Ostap::StatusCode project2
    ( TTree*                     data            ,
      Ostap::Math::LegendreSum2& poly            ,
      const std::string&         expression1     ,
      const std::string&         expression2     ,
      const std::string&         selection  = "" ,
      const Ostap::EventIndex    first      = Ostap::FirstEvent ,
      const Ostap::EventIndex    last       = Ostap::LastEvent  ) const ; 
    // =========================================================================
    /** Project data in 2D-polyoominal: on-flight parameterisation 
     *  @param data        (INPUT)  input data 
     *  @param poly        (UPDATE) the polynoimal 
     *  @param expression1 (INPUT)  expression1 
     *  @param expression2 (INPUT)  expression2 
     *  @param selection   (INPUT)  selection (as boolean)
     *  @param first       (INOPUT) the first event to process (inclusive) 
     *  @param last        (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     */
    Ostap::StatusCode project2
    ( const RooAbsData*          data            ,
      Ostap::Math::LegendreSum2& poly            ,
      const std::string&         expression1     ,
      const std::string&         expression2     ,
      const std::string&         selection  = "" ,
      const std::string&         cut_range  = "" , 
      const Ostap::EventIndex    first      = Ostap::FirstEvent ,
      const Ostap::EventIndex    last       = Ostap::LastEvent  ) const ; 
    // =========================================================================
  public : // 3D-Bernstrein polynomial: on-flight parameterisation 
    // ======================================================================
    /** Project data in 3D-polyoominal: on-flight parameterisation 
     *  @param data        (INPUT)  input data 
     *  @param poly        (UPDATE) the polynoimal 
     *  @param expression1 (INPUT)  expression1 
     *  @param expression2 (INPUT)  expression2 
     *  @param expression3 (INPUT)  expression2 
     *  @param selection   (INPUT)  selection (as boolean)
     *  @param first       (INOPUT) the first event to process (inclusive) 
     *  @param last        (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     */
    Ostap::StatusCode project3
    ( TTree*                     data            ,
      Ostap::Math::Bernstein3D&  poly            ,
      const std::string&         expression1     ,
      const std::string&         expression2     ,
      const std::string&         expression3     ,
      const std::string&         selection  = "" ,
      const Ostap::EventIndex    first      = Ostap::FirstEvent ,
      const Ostap::EventIndex    last       = Ostap::LastEvent  ) const ; 
    // =========================================================================
    /** Project data in 3D-polyoominal: on-flight parameterisation 
     *  @param data        (INPUT)  input data 
     *  @param poly        (UPDATE) the polynoimal 
     *  @param expression1 (INPUT)  expression1 
     *  @param expression2 (INPUT)  expression2 
     *  @param selection   (INPUT)  selection (as boolean)
     *  @param first       (INOPUT) the first event to process (inclusive) 
     *  @param last        (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     */
    Ostap::StatusCode project3
    ( const RooAbsData*          data            ,
      Ostap::Math::Bernstein3D&  poly            ,
      const std::string&         expression1     ,
      const std::string&         expression2     ,
      const std::string&         expression3     ,
      const std::string&         selection  = "" ,
      const std::string&         cut_range  = "" , 
      const Ostap::EventIndex    first      = Ostap::FirstEvent ,
      const Ostap::EventIndex    last       = Ostap::LastEvent  ) const ; 
    // =========================================================================
  public : // 3D-Legendre polynomial: on-flight parameterisation 
    // ======================================================================
    /** Project data in 3D-polyoominal: on-flight parameterisation 
     *  @param data        (INPUT)  input data 
     *  @param poly        (UPDATE) the polynoimal 
     *  @param expression1 (INPUT)  expression1 
     *  @param expression2 (INPUT)  expression2 
     *  @param expression3 (INPUT)  expression2 
     *  @param selection   (INPUT)  selection (as boolean)
     *  @param first       (INOPUT) the first event to process (inclusive) 
     *  @param last        (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     */
    Ostap::StatusCode project3
    ( TTree*                     data            ,
      Ostap::Math::LegendreSum3& poly            ,
      const std::string&         expression1     ,
      const std::string&         expression2     ,
      const std::string&         expression3     ,
      const std::string&         selection  = "" ,
      const Ostap::EventIndex    first      = Ostap::FirstEvent ,
      const Ostap::EventIndex    last       = Ostap::LastEvent  ) const ; 
    // =========================================================================
    /** Project data in 3D-polyoominal: on-flight parameterisation 
     *  @param data        (INPUT)  input data 
     *  @param poly        (UPDATE) the polynoimal 
     *  @param expression1 (INPUT)  expression1 
     *  @param expression2 (INPUT)  expression2 
     *  @param selection   (INPUT)  selection (as boolean)
     *  @param first       (INOPUT) the first event to process (inclusive) 
     *  @param last        (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     */
    Ostap::StatusCode project3
    ( const RooAbsData*          data            ,
      Ostap::Math::LegendreSum3& poly            ,
      const std::string&         expression1     ,
      const std::string&         expression2     ,
      const std::string&         expression3     ,
      const std::string&         selection  = "" ,
      const std::string&         cut_range  = "" , 
      const Ostap::EventIndex    first      = Ostap::FirstEvent ,
      const Ostap::EventIndex    last       = Ostap::LastEvent  ) const ; 
    // =========================================================================
  public : // 3D-Legendre polynomial: on-flight parameterisation 
    // ======================================================================
    /** Project data in 4D-polyoominal: on-flight parameterisation 
     *  @param data        (INPUT)  input data 
     *  @param poly        (UPDATE) the polynoimal 
     *  @param expression1 (INPUT)  expression1 
     *  @param expression2 (INPUT)  expression2 
     *  @param expression3 (INPUT)  expression3 
     *  @param expression4 (INPUT)  expression4 
     *  @param selection   (INPUT)  selection (as boolean)
     *  @param first       (INOPUT) the first event to process (inclusive) 
     *  @param last        (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     */
    Ostap::StatusCode project4
    ( TTree*                     data            ,
      Ostap::Math::LegendreSum4& poly            ,
      const std::string&         expression1     ,
      const std::string&         expression2     ,
      const std::string&         expression3     ,
      const std::string&         expression4     ,
      const std::string&         selection  = "" ,
      const Ostap::EventIndex    first      = Ostap::FirstEvent ,
      const Ostap::EventIndex    last       = Ostap::LastEvent  ) const ; 
    // =========================================================================
    /** Project data in 4D-polyoominal: on-flight parameterisation 
     *  @param data        (INPUT)  input data 
     *  @param poly        (UPDATE) the polynoimal 
     *  @param expression1 (INPUT)  expression1 
     *  @param expression2 (INPUT)  expression2 
     *  @param selection   (INPUT)  selection (as boolean)
     *  @param first       (INOPUT) the first event to process (inclusive) 
     *  @param last        (INOPUT) the last  event to process (exclusive)
     *  @return Status Code 
     */
    Ostap::StatusCode project4
    ( const RooAbsData*          data            ,
      Ostap::Math::LegendreSum4& poly            ,
      const std::string&         expression1     ,
      const std::string&         expression2     ,
      const std::string&         expression3     ,
      const std::string&         expression4     ,
      const std::string&         selection  = "" ,
      const std::string&         cut_range  = "" , 
      const Ostap::EventIndex    first      = Ostap::FirstEvent ,
      const Ostap::EventIndex    last       = Ostap::LastEvent  ) const ; 
    // =========================================================================
  } ; //                                         The end f class Ostap::Project     
  // ==========================================================================
} // The end of namespace Ostap 
// =============================================================================
#endif // OSTAP_PROJECT_H
// =============================================================================
//                                                                       The END 
// =============================================================================
