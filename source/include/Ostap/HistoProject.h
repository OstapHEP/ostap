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
class TH1        ;  // ROOT 
class TH2        ;  // ROOT 
class TH3        ;  // ROOT 
class TProfile   ;  // ROOT 
class TProfile2D ;  // ROOT 
class TTree      ;  // ROOT 
// =============================================================================
class RooAbsData ; // RooFit 
class RooAbsReal ; // RooFit 
// =============================================================================
namespace Ostap { namespace Math { class ChebyshevSum ; } }
namespace Ostap { namespace Math { class LegendreSum  ; } }
namespace Ostap { namespace Math { class LegendreSum2 ; } }
namespace Ostap { namespace Math { class LegendreSum3 ; } }
namespace Ostap { namespace Math { class LegendreSum4 ; } }
namespace Ostap { namespace Math { class Bernstein    ; } }
namespace Ostap { namespace Math { class Bernstein2D  ; } }
namespace Ostap { namespace Math { class Bernstein3D  ; } }
// =============================================================================
namespace Ostap
{  
  // ==========================================================================
  namespace Utils 
  {
    // ========================================================================
    /// forward declaration 
    class ProgressConf ; /// forward declaration 
    // ========================================================================
  }    
  // ==========================================================================
  /** @class HistoProject Ostap/HistoProject.h
   *  Helper class to project Dataset/DataFrame/TTree to histogram 
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date   2015-10-08
   */
  class HistoProject 
  {
  public:
    // ========================================================================
    /** make a projection of RooDataSet into the histogram 
     *  @param data       (INPUT)  input data 
     *  @param histo      (UPDATE) histogram 
     *  @param expression (INPUT) expression
     *  @param selection  (INPUT) selection criteria/weight 
     *  @param first      (INPUT) the first event to process 
     *  @param last       (INPUT) the last event to process 
     */
    static Ostap::StatusCode project
    ( const RooAbsData*   data                 , 
      TH1*                histo                ,
      const std::string&  expression           ,
      const std::string&  selection  = ""      ,
      const char*         range      = nullptr , 
      const unsigned long first      = 0                                         ,
      const unsigned long last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================  
    /** make a projection of RooDataSet into the histogram 
     *  @param data       (INPUT)  input data
     *  @param progress   (INPUT)  configuration of progres bar 
     *  @param histo      (UPDATE) histogram 
     *  @param expression (INPUT) expression
     *  @param selection  (INPUT) selection criteria/weight 
     *  @param range      (INPUT) cut-range 
     *  @param first      (INPUT) the first event to process 
     *  @param last       (INPUT) the last event to process 
     */
    static Ostap::StatusCode project
    ( const RooAbsData*                 data                 , 
      const Ostap::Utils::ProgressConf& progress             ,
      TH1*                              histo                ,
      const std::string&                expression           ,
      const std::string&                selection  = ""      ,
      const char*                       range      = nullptr , 
      const unsigned long               first      = 0                                         ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into the histogram 
     *  @param data        (INPUT)  input data 
     *  @param histo       (UPDATE) histogram 
     *  @param xexpression (INPUT) expression for x-axis 
     *  @param yexpression (INPUT) expression for y-axis 
     *  @param selection   (INPUT) selection criteria/weight 
     *  @param first       (INPUT) the first event to process 
     *  @param last        (INPUT) the last event to process 
     */
    static Ostap::StatusCode project2
    ( const RooAbsData*   data                 , 
      TH2*                histo                ,
      const std::string&  xexpression          ,
      const std::string&  yexpression          ,
      const std::string&  selection  = ""      ,
      const char*         range      = nullptr , 
      const unsigned long first      = 0                                         ,
      const unsigned long last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into the histogram 
     *  @param data        (INPUT)  input data 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param histo       (UPDATE) histogram 
     *  @param xexpression (INPUT) expression for x-axis 
     *  @param yexpression (INPUT) expression for y-axis 
     *  @param selection   (INPUT) selection criteria/weight 
     *  @param first       (INPUT) the first event to process 
     *  @param last        (INPUT) the last event to process 
     */
    static Ostap::StatusCode project2
    ( const RooAbsData*                 data                 ,  
      const Ostap::Utils::ProgressConf& progress             ,
      TH2*                              histo                ,
      const std::string&                xexpression          ,
      const std::string&                yexpression          ,
      const std::string&                selection  = ""      ,
      const char*                       range      = nullptr , 
      const unsigned long               first      = 0                                         ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
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
    ( const RooAbsData*   data                 , 
      TH3*                histo                ,
      const std::string&  xexpression          ,
      const std::string&  yexpression          ,
      const std::string&  zexpression          ,
      const std::string&  selection  = ""      ,
      const char*         range      = nullptr , 
      const unsigned long first      = 0                                         ,
      const unsigned long last       = std::numeric_limits<unsigned long>::max() ) ;

    // ========================================================================
    /** make a projection of RooDataSet into the histogram 
     *  @param data        (INPUT)  input data 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param histo       (UPDATE) histogram 
     *  @param xexpression (INPUT) expression for x-axis 
     *  @param yexpression (INPUT) expression for y-axis 
     *  @param zexpression (INPUT) expression for z-axis 
     *  @param selection   (INPUT) selection criteria/weight 
     *  @param first       (INPUT) the first event to process 
     *  @param last        (INPUT) the last event to process 
     */
    static Ostap::StatusCode project3
    ( const RooAbsData*   data                 , 
      const Ostap::Utils::ProgressConf& progress        ,
      TH3*                histo                ,
      const std::string&  xexpression          ,
      const std::string&  yexpression          ,
      const std::string&  zexpression          ,
      const std::string&  selection  = ""      ,
      const char*         range      = nullptr , 
      const unsigned long first      = 0                                         ,
      const unsigned long last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
  public: // RooDataSet 
    // ========================================================================
    /** make a projection of RooDataSet into the histogram 
     *  @param data       (INPUT)  input data 
     *  @param histo      (UPDATE) histogram 
     *  @param expression (INPUT) expression
     *  @param selection  (INPUT) selection criteria/weight 
     *  @param first      (INPUT) the first event to process 
     *  @param last       (INPUT) the last event to process 
     */
    static Ostap::StatusCode project
    ( const RooAbsData*   data                 , 
      TH1*                histo                ,
      const RooAbsReal&   expression           ,
      const RooAbsReal*   selection  = nullptr ,
      const char*         range      = nullptr , 
      const unsigned long first      = 0                                         ,
      const unsigned long last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into the histogram 
     *  @param data        (INPUT)  input data 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param histo       (UPDATE) histogram 
     *  @param expression  (INPUT)  expression
     *  @param selection   (INPUT)  selection criteria/weight 
     *  @param first       (INPUT)  the first event to process 
     *  @param last        (INPUT)  the last event to process 
     */
    static Ostap::StatusCode project
    ( const RooAbsData*                 data                 , 
      const Ostap::Utils::ProgressConf& progress             ,
      TH1*                              histo                ,
      const RooAbsReal&                 expression           ,
      const RooAbsReal*                 selection  = nullptr ,
      const char*                       range      = nullptr , 
      const unsigned long               first      = 0                                         ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into the histogram 
     *  @param data        (INPUT)  input data 
     *  @param histo       (UPDATE) histogram 
     *  @param xexpression (INPUT) expression for x-axis 
     *  @param yexpression (INPUT) expression for y-axis 
     *  @param selection   (INPUT) selection criteria/weight 
     *  @param first       (INPUT) the first event to process 
     *  @param last        (INPUT) the last event to process 
     */
    static Ostap::StatusCode project2
    ( const RooAbsData*   data                 , 
      TH2*                histo                ,
      const RooAbsReal&   xexpression          ,
      const RooAbsReal&   yexpression          ,
      const RooAbsReal*   selection  = nullptr ,
      const char*         range      = nullptr , 
      const unsigned long first      = 0                                         ,
      const unsigned long last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into the histogram 
     *  @param data        (INPUT)  input data 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param histo       (UPDATE) histogram 
     *  @param xexpression (INPUT) expression for x-axis 
     *  @param yexpression (INPUT) expression for y-axis 
     *  @param selection   (INPUT) selection criteria/weight 
     *  @param first       (INPUT) the first event to process 
     *  @param last        (INPUT) the last event to process 
     */
    static Ostap::StatusCode project2
    ( const RooAbsData*                 data                 , 
      const Ostap::Utils::ProgressConf& progress             ,
      TH2*                              histo                ,
      const RooAbsReal&                 xexpression          ,
      const RooAbsReal&                 yexpression          ,
      const RooAbsReal*                 selection  = nullptr ,
      const char*                       range      = nullptr ,       
      const unsigned long               first      = 0                                         ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into the histogram 
     *  @param data        (INPUT)  input data 
     *  @param histo       (UPDATE) histogram 
     *  @param xexpression (INPUT) expression for x-axis 
     *  @param yexpression (INPUT) expression for y-axis 
     *  @param zexpression (INPUT) expression for z-axis 
     *  @param selection   (INPUT) selection criteria/weight 
     *  @param first       (INPUT) the first event to process 
     *  @param last        (INPUT) the last event to process 
     */
    static Ostap::StatusCode project3
    ( const RooAbsData*   data            , 
      TH3*                histo           ,
      const RooAbsReal&   xexpression     ,
      const RooAbsReal&   yexpression     ,
      const RooAbsReal&   zexpression     ,
      const RooAbsReal*   selection  = nullptr ,
      const char*         range      = nullptr ,       
      const unsigned long first      = 0                                         ,
      const unsigned long last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into the histogram 
     *  @param data        (INPUT)  input data 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param histo       (UPDATE) histogram 
     *  @param xexpression (INPUT)  expression for x-axis 
     *  @param yexpression (INPUT)  expression for y-axis 
     *  @param zexpression (INPUT)  expression for z-axis 
     *  @param selection   (INPUT)  selection criteria/weight 
     *  @param first       (INPUT)  the first event to process 
     *  @param last        (INPUT)  the last event to process 
     */
    static Ostap::StatusCode project3
    ( const RooAbsData*   data            , 
      const Ostap::Utils::ProgressConf& progress        ,
      TH3*                histo                ,
      const RooAbsReal&   xexpression          ,
      const RooAbsReal&   yexpression          ,
      const RooAbsReal&   zexpression          ,
      const RooAbsReal*   selection  = nullptr ,
      const char*         range      = nullptr ,       
      const unsigned long first      = 0                                         ,
      const unsigned long last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
  public:  //   TTree
    // ========================================================================
    /** make a projection of TTree into the histogram 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param data        (INPUT)  input data 
     *  @param histo       (UPDATE) histogram 
     *  @param expression  (INPUT)  expression
     *  @param selection   (INPUT)  selection criteria/weight 
     *  @param first       (INPUT)  the first event to process 
     *  @param last        (INPUT)  the last event to process 
     */
    static Ostap::StatusCode project
    ( TTree*              data            , 
      TH1*                histo           ,
      const std::string&  expression      ,
      const std::string&  selection  = "" ,
      const unsigned long first      = 0  ,
      const unsigned long last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of TTree into the histogram 
     *  @param data        (INPUT)  input data 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param histo       (UPDATE) histogram 
     *  @param expression  (INPUT)  expression
     *  @param selection   (INPUT)  selection criteria/weight 
     *  @param first       (INPUT)  the first event to process 
     *  @param last        (INPUT)  the last event to process 
     */
    static Ostap::StatusCode project
    ( TTree*                            data            , 
      const Ostap::Utils::ProgressConf& progress        ,
      TH1*                              histo           ,
      const std::string&                expression      ,
      const std::string&                selection  = "" ,
      const unsigned long               first      = 0  ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of TTree into the histogram 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param data        (INPUT)  input data 
     *  @param histo       (UPDATE) histogram 
     *  @param xexpression (INPUT)  x-expression
     *  @param yexpression (INPUT)  y-expression
     *  @param selection   (INPUT)  selection criteria/weight 
     *  @param first       (INPUT)  the first event to process 
     *  @param last        (INPUT)  the last event to process 
     */
    static Ostap::StatusCode project2
    ( TTree*              data            , 
      TH2*                histo           ,
      const std::string&  xexpression     ,
      const std::string&  yexpression     ,
      const std::string&  selection  = "" ,
      const unsigned long first      = 0  ,
      const unsigned long last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of TTree into the histogram 
     *  @param data        (INPUT)  input data 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param histo       (UPDATE) histogram 
     *  @param xexpression (INPUT)  x-expression
     *  @param yexpression (INPUT)  y-expression
     *  @param selection   (INPUT)  selection criteria/weight 
     *  @param first       (INPUT)  the first event to process 
     *  @param last        (INPUT)  the last event to process 
     */
    static Ostap::StatusCode project2
    ( TTree*                            data            , 
      const Ostap::Utils::ProgressConf& progress        ,
      TH2*                              histo           ,
      const std::string&                xexpression     ,
      const std::string&                yexpression     ,
      const std::string&                selection  = "" ,
      const unsigned long               first      = 0  ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of TTree into the histogram 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param data        (INPUT)  input data 
     *  @param histo       (UPDATE) histogram 
     *  @param xexpression (INPUT)  x-expression
     *  @param yexpression (INPUT)  y-expression
     *  @param zexpression (INPUT)  z-expression
     *  @param selection   (INPUT)  selection criteria/weight 
     *  @param first       (INPUT)  the first event to process 
     *  @param last        (INPUT)  the last event to process 
     */
    static Ostap::StatusCode project3
    ( TTree*              data            , 
      TH3*                histo           ,
      const std::string&  xexpression     ,
      const std::string&  yexpression     ,
      const std::string&  zexpression     ,
      const std::string&  selection  = "" ,
      const unsigned long first      = 0  ,
      const unsigned long last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of TTree into the histogram 
     *  @param data        (INPUT)  input data 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param histo       (UPDATE) histogram 
     *  @param xexpression (INPUT)  x-expression
     *  @param yexpression (INPUT)  y-expression
     *  @param zexpression (INPUT)  z-expression
     *  @param selection   (INPUT)  selection criteria/weight 
     *  @param first       (INPUT)  the first event to process 
     *  @param last        (INPUT)  the last event to process 
     */
    static Ostap::StatusCode project3
    ( TTree*                            data            , 
      const Ostap::Utils::ProgressConf& progress        ,
      TH3*                              histo           ,
      const std::string&                xexpression     ,
      const std::string&                yexpression     ,
      const std::string&                zexpression     ,
      const std::string&                selection  = "" ,
      const unsigned long               first      = 0  ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
  public:
    // ========================================================================
    /** make a projection of RooDataSet into LegendreSum object 
     *  @param data       (INPUT)  input data 
     *  @param object     (UPDATE) the object 
     *  @param expression (INPUT) expression
     *  @param selection  (INPUT) selection criteria/weight 
     *  @param first      (INPUT) the first event to process 
     *  @param last       (INPUT) the last event to process 
     */
    static Ostap::StatusCode project
    ( const RooAbsData*         data                 , 
      Ostap::Math::LegendreSum& object               ,
      const RooAbsReal&         expression           ,
      const RooAbsReal*         selection  = nullptr ,
      const char*               range      = nullptr , 
      const unsigned long       first      = 0       ,
      const unsigned long       last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into LegendreSum object 
     *  @param data       (INPUT)  input data 
     *  @param progress   (INPUT)  configuration of progres bar 
     *  @param object     (UPDATE) the object 
     *  @param expression (INPUT) expression
     *  @param selection  (INPUT) selection criteria/weight 
     *  @param first      (INPUT) the first event to process 
     *  @param last       (INPUT) the last event to process 
     */
    static Ostap::StatusCode project
    ( const RooAbsData*                 data                 , 
      const Ostap::Utils::ProgressConf& progress             ,
      Ostap::Math::LegendreSum&         object               ,
      const RooAbsReal&                 expression           ,
      const RooAbsReal*                 selection  = nullptr ,
      const char*                       range      = nullptr , 
      const unsigned long               first      = 0       ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into ChebyshevSum object 
     *  @param data       (INPUT)  input data 
     *  @param object     (UPDATE) the object 
     *  @param expression (INPUT) expression
     *  @param selection  (INPUT) selection criteria/weight 
     *  @param first      (INPUT) the first event to process 
     *  @param last       (INPUT) the last event to process 
     */
    static Ostap::StatusCode project
    ( const RooAbsData*          data                 , 
      Ostap::Math::ChebyshevSum& object               ,
      const RooAbsReal&          expression           ,
      const RooAbsReal*          selection  = nullptr ,
      const char*                range      = nullptr , 
      const unsigned long        first      = 0       ,
      const unsigned long        last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into ChebyshevSum object 
     *  @param data       (INPUT)  input data 
     *  @param progress   (INPUT)  configuration of progres bar 
     *  @param object     (UPDATE) the object 
     *  @param expression (INPUT) expression
     *  @param selection  (INPUT) selection criteria/weight 
     *  @param first      (INPUT) the first event to process 
     *  @param last       (INPUT) the last event to process 
     */
    static Ostap::StatusCode project
    ( const RooAbsData*                 data                 , 
      const Ostap::Utils::ProgressConf& progress             ,
      Ostap::Math::ChebyshevSum&        object               ,
      const RooAbsReal&                 expression           ,
      const RooAbsReal*                 selection  = nullptr ,
      const char*                       range      = nullptr , 
      const unsigned long               first      = 0       ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into Bernstein object 
     *  @param data       (INPUT)  input data 
     *  @param object     (UPDATE) the object 
     *  @param expression (INPUT) expression
     *  @param selection  (INPUT) selection criteria/weight 
     *  @param first      (INPUT) the first event to process 
     *  @param last       (INPUT) the last event to process 
     */
    static Ostap::StatusCode project
    ( const RooAbsData*          data                 , 
      Ostap::Math::Bernstein&    object               ,
      const RooAbsReal&          expression           ,
      const RooAbsReal*          selection  = nullptr ,
      const char*                range      = nullptr , 
      const unsigned long        first      = 0       ,
      const unsigned long        last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into Bernstein object 
     *  @param data       (INPUT)  input data 
     *  @param progress   (INPUT)  configuration of progres bar 
     *  @param object     (UPDATE) the object 
     *  @param expression (INPUT) expression
     *  @param selection  (INPUT) selection criteria/weight 
     *  @param first      (INPUT) the first event to process 
     *  @param last       (INPUT) the last event to process 
     */
    static Ostap::StatusCode project
    ( const RooAbsData*                 data                 , 
      const Ostap::Utils::ProgressConf& progress             ,
      Ostap::Math::Bernstein&           object               ,
      const RooAbsReal&                 expression           ,
      const RooAbsReal*                 selection  = nullptr ,
      const char*                       range      = nullptr , 
      const unsigned long               first      = 0       ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into LegendreSum2 object 
     *  @param data        (INPUT)  input data 
     *  @param object      (UPDATE) the object 
     *  @param xexpression (INPUT) x-expression
     *  @param yexpression (INPUT) y-expression
     *  @param selection   (INPUT) selection criteria/weight 
     *  @param first       (INPUT) the first event to process 
     *  @param last        (INPUT) the last event to process 
     */
    static Ostap::StatusCode project2
    ( const RooAbsData*          data                 , 
      Ostap::Math::LegendreSum2& object               ,
      const RooAbsReal&          xexpression          ,
      const RooAbsReal&          yexpression          ,
      const RooAbsReal*          selection  = nullptr ,
      const char*                range      = nullptr , 
      const unsigned long        first      = 0       ,
      const unsigned long        last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into LegendreSum2 object 
     *  @param data        (INPUT)  input data 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param object      (UPDATE) the object 
     *  @param xexpression (INPUT) x-expression
     *  @param yexpression (INPUT) y-expression
     *  @param selection   (INPUT) selection criteria/weight 
     *  @param range       (INPUT) selection range  
     *  @param first       (INPUT) the first event to process 
     *  @param last        (INPUT) the last event to process 
     */
    static Ostap::StatusCode project2
    ( const RooAbsData*                 data                 , 
      const Ostap::Utils::ProgressConf& progress             ,
      Ostap::Math::LegendreSum2&        object               ,
      const RooAbsReal&                 xexpression          ,
      const RooAbsReal&                 yexpression          ,
      const RooAbsReal*                 selection  = nullptr ,
      const char*                       range      = nullptr , 
      const unsigned long               first      = 0       ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into Bernstein2D object 
     *  @param data        (INPUT)  input data 
     *  @param object      (UPDATE) the object 
     *  @param xexpression (INPUT) x-expression
     *  @param yexpression (INPUT) y-expression
     *  @param selection   (INPUT) selection criteria/weight 
     *  @param first       (INPUT) the first event to process 
     *  @param last        (INPUT) the last event to process 
     */
    static Ostap::StatusCode project2
    ( const RooAbsData*          data                 , 
      Ostap::Math::Bernstein2D&  object               ,
      const RooAbsReal&          xexpression          ,
      const RooAbsReal&          yexpression          ,
      const RooAbsReal*          selection  = nullptr ,
      const char*                range      = nullptr , 
      const unsigned long        first      = 0       ,
      const unsigned long        last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into Bernstein2D object 
     *  @param data        (INPUT)  input data 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param object      (UPDATE) the object 
     *  @param xexpression (INPUT) x-expression
     *  @param yexpression (INPUT) y-expression
     *  @param selection   (INPUT) selection criteria/weight 
     *  @param range       (INPUT) selection range  
     *  @param first       (INPUT) the first event to process 
     *  @param last        (INPUT) the last event to process 
     */
    static Ostap::StatusCode project2
    ( const RooAbsData*                 data                 , 
      const Ostap::Utils::ProgressConf& progress             ,
      Ostap::Math::Bernstein2D&         object               ,
      const RooAbsReal&                 xexpression          ,
      const RooAbsReal&                 yexpression          ,
      const RooAbsReal*                 selection  = nullptr ,
      const char*                       range      = nullptr , 
      const unsigned long               first      = 0       ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into LegendreSum3 object 
     *  @param data        (INPUT)  input data 
     *  @param object      (UPDATE) the object 
     *  @param xexpression (INPUT) x-expression
     *  @param yexpression (INPUT) y-expression
     *  @param zexpression (INPUT) z-expression
     *  @param selection   (INPUT) selection criteria/weight 
     *  @param first       (INPUT) the first event to process 
     *  @param last        (INPUT) the last event to process 
     */
    static Ostap::StatusCode project3
    ( const RooAbsData*          data                 , 
      Ostap::Math::LegendreSum3& object               ,
      const RooAbsReal&          xexpression          ,
      const RooAbsReal&          yexpression          ,
      const RooAbsReal&          zexpression          ,
      const RooAbsReal*          selection  = nullptr ,
      const char*                range      = nullptr , 
      const unsigned long        first      = 0       ,
      const unsigned long        last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into LegendreSum3 object 
     *  @param data        (INPUT)  input data 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param object      (UPDATE) the object 
     *  @param xexpression (INPUT) x-expression
     *  @param yexpression (INPUT) y-expression
     *  @param zexpression (INPUT) z-expression
     *  @param selection   (INPUT) selection criteria/weight 
     *  @param range       (INPUT) selection range  
     *  @param first       (INPUT) the first event to process 
     *  @param last        (INPUT) the last event to process 
     */
    static Ostap::StatusCode project3
    ( const RooAbsData*                 data                 , 
      const Ostap::Utils::ProgressConf& progress             ,
      Ostap::Math::LegendreSum3&        object               ,
      const RooAbsReal&                 xexpression          ,
      const RooAbsReal&                 yexpression          ,
      const RooAbsReal&                 zexpression          ,
      const RooAbsReal*                 selection  = nullptr ,
      const char*                       range      = nullptr , 
      const unsigned long               first      = 0       ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into Bernstein3D object 
     *  @param data        (INPUT)  input data 
     *  @param object      (UPDATE) the object 
     *  @param xexpression (INPUT) x-expression
     *  @param yexpression (INPUT) y-expression
     *  @param zexpression (INPUT) z-expression
     *  @param selection   (INPUT) selection criteria/weight 
     *  @param first       (INPUT) the first event to process 
     *  @param last        (INPUT) the last event to process 
     */
    static Ostap::StatusCode project3
    ( const RooAbsData*          data                 , 
      Ostap::Math::Bernstein3D&  object               ,
      const RooAbsReal&          xexpression          ,
      const RooAbsReal&          yexpression          ,
      const RooAbsReal&          zexpression          ,
      const RooAbsReal*          selection  = nullptr ,
      const char*                range      = nullptr , 
      const unsigned long        first      = 0       ,
      const unsigned long        last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into Bernstein3D object 
     *  @param data        (INPUT)  input data 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param object      (UPDATE) the object 
     *  @param xexpression (INPUT) x-expression
     *  @param yexpression (INPUT) y-expression
     *  @param zexpression (INPUT) z-expression
     *  @param selection   (INPUT) selection criteria/weight 
     *  @param range       (INPUT) selection range  
     *  @param first       (INPUT) the first event to process 
     *  @param last        (INPUT) the last event to process 
     */
    static Ostap::StatusCode project3
    ( const RooAbsData*                 data                 , 
      const Ostap::Utils::ProgressConf& progress             ,
      Ostap::Math::Bernstein3D&         object               ,
      const RooAbsReal&                 xexpression          ,
      const RooAbsReal&                 yexpression          ,
      const RooAbsReal&                 zexpression          ,
      const RooAbsReal*                 selection  = nullptr ,
      const char*                       range      = nullptr , 
      const unsigned long               first      = 0       ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into LegendreSum4 object 
     *  @param data        (INPUT)  input data 
     *  @param object      (UPDATE) the object 
     *  @param xexpression (INPUT) x-expression
     *  @param yexpression (INPUT) y-expression
     *  @param zexpression (INPUT) z-expression
     *  @param uexpression (INPUT) u-expression
     *  @param selection   (INPUT) selection criteria/weight 
     *  @param first       (INPUT) the first event to process 
     *  @param last        (INPUT) the last event to process 
     */
    static Ostap::StatusCode project4
    ( const RooAbsData*          data                 , 
      Ostap::Math::LegendreSum4& object               ,
      const RooAbsReal&          xexpression          ,
      const RooAbsReal&          yexpression          ,
      const RooAbsReal&          zexpression          ,
      const RooAbsReal&          uexpression          ,
      const RooAbsReal*          selection  = nullptr ,
      const char*                range      = nullptr , 
      const unsigned long        first      = 0       ,
      const unsigned long        last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into LegendreSum4 object 
     *  @param data        (INPUT)  input data 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param object      (UPDATE) the object 
     *  @param xexpression (INPUT) x-expression
     *  @param yexpression (INPUT) y-expression
     *  @param zexpression (INPUT) z-expression
     *  @param uexpression (INPUT) u-expression
     *  @param selection   (INPUT) selection criteria/weight 
     *  @param range       (INPUT) selection range  
     *  @param first       (INPUT) the first event to process 
     *  @param last        (INPUT) the last event to process 
     */
    static Ostap::StatusCode project4
    ( const RooAbsData*                 data                 , 
      const Ostap::Utils::ProgressConf& progress             ,
      Ostap::Math::LegendreSum4&        object               ,
      const RooAbsReal&                 xexpression          ,
      const RooAbsReal&                 yexpression          ,
      const RooAbsReal&                 zexpression          ,
      const RooAbsReal&                 uexpression          ,
      const RooAbsReal*                 selection  = nullptr ,
      const char*                       range      = nullptr , 
      const unsigned long               first      = 0       ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
  public:
    // ========================================================================
    /** make a projection of RooDataSet into LegendreSum object 
     *  @param data       (INPUT)  input data 
     *  @param object     (UPDATE) the object 
     *  @param expression (INPUT) expression
     *  @param selection  (INPUT) selection criteria/weight 
     *  @param first      (INPUT) the first event to process 
     *  @param last       (INPUT) the last event to process 
     */
    static Ostap::StatusCode project
    ( const RooAbsData*         data                 , 
      Ostap::Math::LegendreSum& object               ,
      const std::string&        expression           , 
      const std::string&        selection  = ""      , 
      const char*               range      = nullptr , 
      const unsigned long       first      = 0       ,
      const unsigned long       last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into LegendreSum object 
     *  @param data       (INPUT)  input data 
     *  @param progress   (INPUT)  configuration of progres bar 
     *  @param object     (UPDATE) the object 
     *  @param expression (INPUT) expression
     *  @param selection  (INPUT) selection criteria/weight 
     *  @param first      (INPUT) the first event to process 
     *  @param last       (INPUT) the last event to process 
     */
    static Ostap::StatusCode project
    ( const RooAbsData*                 data                 , 
      const Ostap::Utils::ProgressConf& progress             ,
      Ostap::Math::LegendreSum&         object               ,
      const std::string&                expression           , 
      const std::string&                selection  = ""      , 
      const char*                       range      = nullptr , 
      const unsigned long               first      = 0       ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into ChebyshevSum object 
     *  @param data       (INPUT)  input data 
     *  @param object     (UPDATE) the object 
     *  @param expression (INPUT) expression
     *  @param selection  (INPUT) selection criteria/weight 
     *  @param first      (INPUT) the first event to process 
     *  @param last       (INPUT) the last event to process 
     */
    static Ostap::StatusCode project
    ( const RooAbsData*          data                 , 
      Ostap::Math::ChebyshevSum& object               ,
      const std::string&         expression           , 
      const std::string&         selection  = ""      ,
      const char*                range      = nullptr , 
      const unsigned long        first      = 0       ,
      const unsigned long        last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into ChebyshevSum object 
     *  @param data       (INPUT)  input data 
     *  @param progress   (INPUT)  configuration of progres bar 
     *  @param object     (UPDATE) the object 
     *  @param expression (INPUT) expression
     *  @param selection  (INPUT) selection criteria/weight 
     *  @param first      (INPUT) the first event to process 
     *  @param last       (INPUT) the last event to process 
     */
    static Ostap::StatusCode project
    ( const RooAbsData*                 data                 , 
      const Ostap::Utils::ProgressConf& progress             ,
      Ostap::Math::ChebyshevSum&        object               ,
      const std::string&                expression           , 
      const std::string&                selection  = ""      , 
      const char*                       range      = nullptr , 
      const unsigned long               first      = 0       ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into Bernstein object 
     *  @param data       (INPUT)  input data 
     *  @param object     (UPDATE) the object 
     *  @param expression (INPUT) expression
     *  @param selection  (INPUT) selection criteria/weight 
     *  @param first      (INPUT) the first event to process 
     *  @param last       (INPUT) the last event to process 
     */
    static Ostap::StatusCode project
    ( const RooAbsData*         data                 , 
      Ostap::Math::Bernstein&   object               ,
      const std::string&        expression           , 
      const std::string&        selection  = ""      , 
      const char*               range      = nullptr , 
      const unsigned long       first      = 0       ,
      const unsigned long       last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into Bernstein  object 
     *  @param data       (INPUT)  input data 
     *  @param progress   (INPUT)  configuration of progres bar 
     *  @param object     (UPDATE) the object 
     *  @param expression (INPUT) expression
     *  @param selection  (INPUT) selection criteria/weight 
     *  @param first      (INPUT) the first event to process 
     *  @param last       (INPUT) the last event to process 
     */
    static Ostap::StatusCode project
    ( const RooAbsData*                 data                 , 
      const Ostap::Utils::ProgressConf& progress             ,
      Ostap::Math::Bernstein&           object               ,
      const std::string&                expression           , 
      const std::string&                selection  = ""      , 
      const char*                       range      = nullptr , 
      const unsigned long               first      = 0       ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;    
    // ========================================================================
    /** make a projection of RooDataSet into LegendreSum2 object 
     *  @param data        (INPUT)  input data 
     *  @param object      (UPDATE) the object 
     *  @param xexpression (INPUT) x-expression
     *  @param yexpression (INPUT) y-expression
     *  @param selection   (INPUT) selection criteria/weight 
     *  @param first       (INPUT) the first event to process 
     *  @param last        (INPUT) the last event to process 
     */
    static Ostap::StatusCode project2
    ( const RooAbsData*          data                 , 
      Ostap::Math::LegendreSum2& object               ,
      const std::string&         xexpression           , 
      const std::string&         yexpression           , 
      const std::string&         selection  = ""      , 
      const char*                range      = nullptr , 
      const unsigned long        first      = 0       ,
      const unsigned long        last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into LegendreSum2 object 
     *  @param data       (INPUT)  input data 
     *  @param progress   (INPUT)  configuration of progres bar 
     *  @param object     (UPDATE) the object 
     *  @param expression (INPUT) expression
     *  @param selection  (INPUT) selection criteria/weight 
     *  @param first      (INPUT) the first event to process 
     *  @param last       (INPUT) the last event to process 
     */
    static Ostap::StatusCode project2
    ( const RooAbsData*                 data                 , 
      const Ostap::Utils::ProgressConf& progress             ,
      Ostap::Math::LegendreSum2&        object               ,
      const std::string&                xexpression          , 
      const std::string&                yexpression          , 
      const std::string&                selection  = ""      , 
      const char*                       range      = nullptr , 
      const unsigned long               first      = 0       ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into Bernstein2D object 
     *  @param data        (INPUT)  input data 
     *  @param object      (UPDATE) the object 
     *  @param xexpression (INPUT) x-expression
     *  @param yexpression (INPUT) y-expression
     *  @param selection   (INPUT) selection criteria/weight 
     *  @param first       (INPUT) the first event to process 
     *  @param last        (INPUT) the last event to process 
     */
    static Ostap::StatusCode project2
    ( const RooAbsData*          data                 , 
      Ostap::Math::Bernstein2D&  object               ,
      const std::string&         xexpression          , 
      const std::string&         yexpression          , 
      const std::string&         selection  = ""      , 
      const char*                range      = nullptr , 
      const unsigned long        first      = 0       ,
      const unsigned long        last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into Bernstein2D object 
     *  @param data       (INPUT)  input data 
     *  @param progress   (INPUT)  configuration of progres bar 
     *  @param object     (UPDATE) the object 
     *  @param expression (INPUT) expression
     *  @param selection  (INPUT) selection criteria/weight 
     *  @param first      (INPUT) the first event to process 
     *  @param last       (INPUT) the last event to process 
     */
    static Ostap::StatusCode project2
    ( const RooAbsData*                 data                 , 
      const Ostap::Utils::ProgressConf& progress             ,
      Ostap::Math::Bernstein2D&         object               ,
      const std::string&                xexpression          , 
      const std::string&                yexpression          , 
      const std::string&                selection  = ""      , 
      const char*                       range      = nullptr , 
      const unsigned long               first      = 0       ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================    
    /** make a projection of RooDataSet into LegendreSum3 object 
     *  @param data        (INPUT)  input data 
     *  @param object      (UPDATE) the object 
     *  @param xexpression (INPUT) x-expression
     *  @param yexpression (INPUT) y-expression
     *  @param zexpression (INPUT) z-expression
     *  @param selection   (INPUT) selection criteria/weight 
     *  @param first       (INPUT) the first event to process 
     *  @param last        (INPUT) the last event to process 
     */
    static Ostap::StatusCode project3
    ( const RooAbsData*          data                 , 
      Ostap::Math::LegendreSum3& object               ,
      const std::string&         xexpression          , 
      const std::string&         yexpression          , 
      const std::string&         zexpression          , 
      const std::string&         selection  = ""      , 
      const char*                range      = nullptr , 
      const unsigned long        first      = 0       ,
      const unsigned long        last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into LegendreSum3 object 
     *  @param data       (INPUT)  input data 
     *  @param progress   (INPUT)  configuration of progres bar 
     *  @param object     (UPDATE) the object 
     *  @param xexpression (INPUT) x-expression
     *  @param yexpression (INPUT) y-expression
     *  @param zexpression (INPUT) z-expression
     *  @param selection  (INPUT) selection criteria/weight 
     *  @param first      (INPUT) the first event to process 
     *  @param last       (INPUT) the last event to process 
     */
    static Ostap::StatusCode project3
    ( const RooAbsData*                 data                 , 
      const Ostap::Utils::ProgressConf& progress             ,
      Ostap::Math::LegendreSum3&        object               ,
      const std::string&                xexpression          , 
      const std::string&                yexpression          , 
      const std::string&                zexpression          , 
      const std::string&                selection  = ""      , 
      const char*                       range      = nullptr , 
      const unsigned long               first      = 0       ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into Bernstein3D object 
     *  @param data        (INPUT)  input data 
     *  @param object      (UPDATE) the object 
     *  @param xexpression (INPUT) x-expression
     *  @param yexpression (INPUT) y-expression
     *  @param zexpression (INPUT) z-expression
     *  @param selection   (INPUT) selection criteria/weight 
     *  @param first       (INPUT) the first event to process 
     *  @param last        (INPUT) the last event to process 
     */
    static Ostap::StatusCode project3
    ( const RooAbsData*          data                 , 
      Ostap::Math::Bernstein3D&  object               ,
      const std::string&         xexpression          , 
      const std::string&         yexpression          , 
      const std::string&         zexpression          , 
      const std::string&         selection  = ""      ,
      const char*                range      = nullptr ,       
      const unsigned long        first      = 0       ,
      const unsigned long        last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into Bernstein3D object 
     *  @param data       (INPUT)  input data 
     *  @param progress   (INPUT)  configuration of progres bar 
     *  @param object     (UPDATE) the object 
     *  @param xexpression (INPUT) x-expression
     *  @param yexpression (INPUT) y-expression
     *  @param zexpression (INPUT) z-expression
     *  @param selection  (INPUT) selection criteria/weight 
     *  @param first      (INPUT) the first event to process 
     *  @param last       (INPUT) the last event to process 
     */
    static Ostap::StatusCode project3
    ( const RooAbsData*                 data                 , 
      const Ostap::Utils::ProgressConf& progress             ,
      Ostap::Math::Bernstein3D&         object               ,
      const std::string&                xexpression          , 
      const std::string&                yexpression          , 
      const std::string&                zexpression          , 
      const std::string&                selection  = ""      ,
      const char*                       range      = nullptr , 
      const unsigned long               first      = 0       ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================    
    /** make a projection of RooDataSet into LegendreSum4 object 
     *  @param data        (INPUT)  input data 
     *  @param object      (UPDATE) the object 
     *  @param xexpression (INPUT) x-expression
     *  @param yexpression (INPUT) y-expression
     *  @param zexpression (INPUT) z-expression
     *  @param uexpression (INPUT) u-expression
     *  @param selection   (INPUT) selection criteria/weight 
     *  @param first       (INPUT) the first event to process 
     *  @param last        (INPUT) the last event to process 
     */
    static Ostap::StatusCode project4
    ( const RooAbsData*          data                 , 
      Ostap::Math::LegendreSum4& object               ,
      const std::string&         xexpression          , 
      const std::string&         yexpression          , 
      const std::string&         zexpression          , 
      const std::string&         uexpression          , 
      const std::string&         selection  = ""      , 
      const char*                range      = nullptr , 
      const unsigned long        first      = 0       ,
      const unsigned long        last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of RooDataSet into LegendreSum4 object 
     *  @param data       (INPUT)  input data 
     *  @param progress   (INPUT)  configuration of progres bar 
     *  @param object     (UPDATE) the object 
     *  @param xexpression (INPUT) x-expression
     *  @param yexpression (INPUT) y-expression
     *  @param zexpression (INPUT) z-expression
     *  @param selection  (INPUT) selection criteria/weight 
     *  @param first      (INPUT) the first event to process 
     *  @param last       (INPUT) the last event to process 
     */
    static Ostap::StatusCode project4
    ( const RooAbsData*                 data                 , 
      const Ostap::Utils::ProgressConf& progress             ,
      Ostap::Math::LegendreSum4&        object               ,
      const std::string&                xexpression          , 
      const std::string&                yexpression          , 
      const std::string&                zexpression          , 
      const std::string&                uexpression          , 
      const std::string&                selection  = ""      , 
      const char*                       range      = nullptr , 
      const unsigned long               first      = 0       ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
  public:  //   DataFrame 
    // ========================================================================
    /** make a projection of DataFrame into the histogram 
     *  @param data        (INPUT)  input data 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param histo       (UPDATE) histogram 
     *  @param expression  (INPUT) expression
     *  @param selection   (INPUT) selection criteria/weight 
     */
    static Ostap::StatusCode project
    ( FrameNode           data            , 
      TH1*                histo           ,
      const std::string&  expression      ,
      const std::string&  selection  = "" ) ;
    // ========================================================================
    /** make a projection of DataFrame into the histogram 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param data        (INPUT)  input data 
     *  @param histo       (UPDATE) histogram 
     *  @param expression  (INPUT) expression
     *  @param selection   (INPUT) selection criteria/weight 
     */
    static Ostap::StatusCode project
    ( FrameNode                         data            , 
      const Ostap::Utils::ProgressConf& progress        ,
      TH1*                              histo           ,
      const std::string&                expression      ,
      const std::string&                selection  = "" ) ;
    // ========================================================================
    /** make a projection of RooDataSet into the histogram 
     *  @param data        (INPUT)  input data 
     *  @param histo       (UPDATE) histogram 
     *  @param xexpression (INPUT) expression for x-axis 
     *  @param yexpression (INPUT) expression for y-axis 
     *  @param selection   (INPUT) selection criteria/weight 
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
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param histo       (UPDATE) histogram 
     *  @param xexpression (INPUT) expression for x-axis 
     *  @param yexpression (INPUT) expression for y-axis 
     *  @param selection   (INPUT) selection criteria/weight 
     */
    static Ostap::StatusCode project2
    ( FrameNode                         data            , 
      const Ostap::Utils::ProgressConf& progress        ,
      TH2*                              histo           ,
      const std::string&                xexpression     ,
      const std::string&                yexpression     ,
      const std::string&                selection  = "" );
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
    /** make a projection of RooDataSet into the histogram 
     *  @param data        (INPUT)  input data 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param histo       (UPDATE) histogram 
     *  @param xexpression (INPUT)  expression for x-axis 
     *  @param yexpression (INPUT)  expression for y-axis 
     *  @param zexpression (INPUT)  expression for z-axis 
     *  @param selection   (INPUT)  selection criteria/weight 
     */
    static Ostap::StatusCode project3
    ( FrameNode                         data            , 
      const Ostap::Utils::ProgressConf& progress        ,
      TH3*                              histo           ,
      const std::string&                xexpression     ,
      const std::string&                yexpression     ,
      const std::string&                zexpression     ,
      const std::string&                selection  = "" ) ;
    // ========================================================================
  public: 
    // ========================================================================
    // TTeee -> 1D non-histograms 
    // ========================================================================
    /** make a projection of TTree into the object 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param data        (INPUT)  input data 
     *  @param sum         (UPDATE) sum  
     *  @param expression  (INPUT)  expression
     *  @param selection   (INPUT)  selection criteria/weight 
     *  @param first       (INPUT)  the first event to process 
     *  @param last        (INPUT)  the last event to process 
     */
    static Ostap::StatusCode project
    ( TTree*                     data            , 
      Ostap::Math::LegendreSum&  sum             ,
      const std::string&         expression      ,
      const std::string&         selection  = "" ,
      const unsigned long        first      = 0  ,
      const unsigned long        last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of TTree into the object 
     *  @param data        (INPUT)  input data 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param sum         (UPDATE) sum 
     *  @param expression  (INPUT)  expression
     *  @param selection   (INPUT)  selection criteria/weight 
     *  @param first       (INPUT)  the first event to process 
     *  @param last        (INPUT)  the last event to process 
     */
    static Ostap::StatusCode project
    ( TTree*                            data            , 
      const Ostap::Utils::ProgressConf& progress        ,
      Ostap::Math::LegendreSum&         sum             ,
      const std::string&                expression      ,
      const std::string&                selection  = "" ,
      const unsigned long               first      = 0  ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of TTree into the object 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param data        (INPUT)  input data 
     *  @param sum         (UPDATE) sum  
     *  @param expression  (INPUT)  expression
     *  @param selection   (INPUT)  selection criteria/weight 
     *  @param first       (INPUT)  the first event to process 
     *  @param last        (INPUT)  the last event to process 
     */
    static Ostap::StatusCode project
    ( TTree*                     data            , 
      Ostap::Math::ChebyshevSum& sum             ,
      const std::string&         expression      ,
      const std::string&         selection  = "" ,
      const unsigned long        first      = 0  ,
      const unsigned long        last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of TTree into the object 
     *  @param data        (INPUT)  input data 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param sum         (UPDATE) sum 
     *  @param expression  (INPUT)  expression
     *  @param selection   (INPUT)  selection criteria/weight 
     *  @param first       (INPUT)  the first event to process 
     *  @param last        (INPUT)  the last event to process 
     */
    static Ostap::StatusCode project
    ( TTree*                            data            , 
      const Ostap::Utils::ProgressConf& progress        ,
      Ostap::Math::ChebyshevSum&        sum             ,
      const std::string&                expression      ,
      const std::string&                selection  = "" ,
      const unsigned long               first      = 0  ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of TTree into the object 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param data        (INPUT)  input data 
     *  @param sum         (UPDATE) sum  
     *  @param expression  (INPUT)  expression
     *  @param selection   (INPUT)  selection criteria/weight 
     *  @param first       (INPUT)  the first event to process 
     *  @param last        (INPUT)  the last event to process 
     */
    static Ostap::StatusCode project
    ( TTree*                     data            , 
      Ostap::Math::Bernstein&    sum             ,
      const std::string&         expression      ,
      const std::string&         selection  = "" ,
      const unsigned long        first      = 0  ,
      const unsigned long        last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of TTree into the object 
     *  @param data        (INPUT)  input data 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param sum         (UPDATE) sum 
     *  @param expression  (INPUT)  expression
     *  @param selection   (INPUT)  selection criteria/weight 
     *  @param first       (INPUT)  the first event to process 
     *  @param last        (INPUT)  the last event to process 
     */
    static Ostap::StatusCode project
    ( TTree*                            data            , 
      const Ostap::Utils::ProgressConf& progress        ,
      Ostap::Math::Bernstein&           sum             ,
      const std::string&                expression      ,
      const std::string&                selection  = "" ,
      const unsigned long               first      = 0  ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
  public: 
    // ========================================================================
    // TTeee -> 2D non-histograms 
    // ========================================================================
    /** make a projection of TTree into the object
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param data        (INPUT)  input data 
     *  @param sum         (UPDATE) sum 
     *  @param xexpression (INPUT)  x-expression
     *  @param yexpression (INPUT)  y-expression
     *  @param selection   (INPUT)  selection criteria/weight 
     *  @param first       (INPUT)  the first event to process 
     *  @param last        (INPUT)  the last event to process 
     */
    static Ostap::StatusCode project2
    ( TTree*                     data            , 
      Ostap::Math::LegendreSum2& sum             ,
      const std::string&         xexpression     ,
      const std::string&         yexpression     ,
      const std::string&         selection  = "" ,
      const unsigned long        first      = 0  ,
      const unsigned long        last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of TTree into the object 
     *  @param data        (INPUT)  input data 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param sum         (UPDATE) sum 
     *  @param xexpression (INPUT)  x-expression
     *  @param yexpression (INPUT)  y-expression
     *  @param selection   (INPUT)  selection criteria/weight 
     *  @param first       (INPUT)  the first event to process 
     *  @param last        (INPUT)  the last event to process 
     */
    static Ostap::StatusCode project2
    ( TTree*                            data            , 
      const Ostap::Utils::ProgressConf& progress        ,
      Ostap::Math::LegendreSum2&        sum             ,
      const std::string&                xexpression     ,
      const std::string&                yexpression     ,
      const std::string&                selection  = "" ,
      const unsigned long               first      = 0  ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of TTree into the object
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param data        (INPUT)  input data 
     *  @param sum         (UPDATE) sum 
     *  @param xexpression (INPUT)  x-expression
     *  @param yexpression (INPUT)  y-expression
     *  @param selection   (INPUT)  selection criteria/weight 
     *  @param first       (INPUT)  the first event to process 
     *  @param last        (INPUT)  the last event to process 
     */
    static Ostap::StatusCode project2
    ( TTree*                     data            , 
      Ostap::Math::Bernstein2D&  sum             ,
      const std::string&         xexpression     ,
      const std::string&         yexpression     ,
      const std::string&         selection  = "" ,
      const unsigned long        first      = 0  ,
      const unsigned long        last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of TTree into the object 
     *  @param data        (INPUT)  input data 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param sum         (UPDATE) sum 
     *  @param xexpression (INPUT)  x-expression
     *  @param yexpression (INPUT)  y-expression
     *  @param selection   (INPUT)  selection criteria/weight 
     *  @param first       (INPUT)  the first event to process 
     *  @param last        (INPUT)  the last event to process 
     */
    static Ostap::StatusCode project2
    ( TTree*                            data            , 
      const Ostap::Utils::ProgressConf& progress        ,
      Ostap::Math::Bernstein2D&         sum             ,
      const std::string&                xexpression     ,
      const std::string&                yexpression     ,
      const std::string&                selection  = "" ,
      const unsigned long               first      = 0  ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
  public: 
    // ========================================================================
    // TTeee -> 3D non-histograms 
    // ========================================================================
    /** make a projection of TTree into the object
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param data        (INPUT)  input data 
     *  @param sum         (UPDATE) sum 
     *  @param xexpression (INPUT)  x-expression
     *  @param yexpression (INPUT)  y-expression
     *  @param zexpression (INPUT)  z-expression
     *  @param selection   (INPUT)  selection criteria/weight 
     *  @param first       (INPUT)  the first event to process 
     *  @param last        (INPUT)  the last event to process 
     */
    static Ostap::StatusCode project3
    ( TTree*                     data            , 
      Ostap::Math::LegendreSum3& sum             ,
      const std::string&         xexpression     ,
      const std::string&         yexpression     ,
      const std::string&         zexpression     ,
      const std::string&         selection  = "" ,
      const unsigned long        first      = 0  ,
      const unsigned long        last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of TTree into the object 
     *  @param data        (INPUT)  input data 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param sum         (UPDATE) sum 
     *  @param xexpression (INPUT)  x-expression
     *  @param yexpression (INPUT)  y-expression
     *  @param zexpression (INPUT)  z-expression
     *  @param selection   (INPUT)  selection criteria/weight 
     *  @param first       (INPUT)  the first event to process 
     *  @param last        (INPUT)  the last event to process 
     */
    static Ostap::StatusCode project3
    ( TTree*                            data            , 
      const Ostap::Utils::ProgressConf& progress        ,
      Ostap::Math::LegendreSum3&        sum             ,
      const std::string&                xexpression     ,
      const std::string&                yexpression     ,
      const std::string&                zexpression     ,
      const std::string&                selection  = "" ,
      const unsigned long               first      = 0  ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of TTree into the object
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param data        (INPUT)  input data 
     *  @param sum         (UPDATE) sum 
     *  @param xexpression (INPUT)  x-expression
     *  @param yexpression (INPUT)  y-expression
     *  @param zexpression (INPUT)  z-expression
     *  @param selection   (INPUT)  selection criteria/weight 
     *  @param first       (INPUT)  the first event to process 
     *  @param last        (INPUT)  the last event to process 
     */
    static Ostap::StatusCode project3
    ( TTree*                     data            , 
      Ostap::Math::Bernstein3D&  sum             ,
      const std::string&         xexpression     ,
      const std::string&         yexpression     ,
      const std::string&         zexpression     ,
      const std::string&         selection  = "" ,
      const unsigned long        first      = 0  ,
      const unsigned long        last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of TTree into the object 
     *  @param data        (INPUT)  input data 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param sum         (UPDATE) sum 
     *  @param xexpression (INPUT)  x-expression
     *  @param yexpression (INPUT)  y-expression
     *  @param zexpression (INPUT)  z-expression
     *  @param selection   (INPUT)  selection criteria/weight 
     *  @param first       (INPUT)  the first event to process 
     *  @param last        (INPUT)  the last event to process 
     */
    static Ostap::StatusCode project3
    ( TTree*                            data            , 
      const Ostap::Utils::ProgressConf& progress        ,
      Ostap::Math::Bernstein3D&         sum             ,
      const std::string&                xexpression     ,
      const std::string&                yexpression     ,
      const std::string&                zexpression     ,
      const std::string&                selection  = "" ,
      const unsigned long               first      = 0  ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
  public: 
    // ========================================================================
    // TTeee -> 4D non-histograms 
    // ========================================================================
    /** make a projection of TTree into the object
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param data        (INPUT)  input data 
     *  @param sum         (UPDATE) sum 
     *  @param xexpression (INPUT)  x-expression
     *  @param yexpression (INPUT)  y-expression
     *  @param zexpression (INPUT)  z-expression
     *  @param zexpression (INPUT)  z-expression
     *  @param selection   (INPUT)  selection criteria/weight 
     *  @param first       (INPUT)  the first event to process 
     *  @param last        (INPUT)  the last event to process 
     */
    static Ostap::StatusCode project4
    ( TTree*                     data            , 
      Ostap::Math::LegendreSum4& sum             ,
      const std::string&         xexpression     ,
      const std::string&         yexpression     ,
      const std::string&         zexpression     ,
      const std::string&         uexpression     ,
      const std::string&         selection  = "" ,
      const unsigned long        first      = 0  ,
      const unsigned long        last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** make a projection of TTree into the object 
     *  @param data        (INPUT)  input data 
     *  @param progress    (INPUT)  configuration of progres bar 
     *  @param sum         (UPDATE) sum 
     *  @param xexpression (INPUT)  x-expression
     *  @param yexpression (INPUT)  y-expression
     *  @param zexpression (INPUT)  z-expression
     *  @param selection   (INPUT)  selection criteria/weight 
     *  @param first       (INPUT)  the first event to process 
     *  @param last        (INPUT)  the last event to process 
     */
    static Ostap::StatusCode project4
    ( TTree*                            data            , 
      const Ostap::Utils::ProgressConf& progress        ,
      Ostap::Math::LegendreSum4&        sum             ,
      const std::string&                xexpression     ,
      const std::string&                yexpression     ,
      const std::string&                zexpression     ,
      const std::string&                uexpression     ,
      const std::string&                selection  = "" ,
      const unsigned long               first      = 0  ,
      const unsigned long               last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
  } ;
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_HPROJECT_H
// ============================================================================

