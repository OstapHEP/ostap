// ============================================================================
#ifndef OSTAP_STATVAR_H
#define OSTAP_STATVAR_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <limits>
// ============================================================================
// Forward declarations 
// =============================================================================
class TTree      ; // ROOT 
class TChain     ; // ROOT 
class RooAbsData ; // RooFit
// =============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatEntity.h"
#include "Ostap/WStatEntity.h"
#include "Ostap/ValueWithError.h"
#include "Ostap/Covariance.h"
#include "Ostap/SymmetricMatrixTypes.h"
#include "Ostap/DataFrame.h"
#include "Ostap/StatusCode.h"
// ============================================================================
// forward declarations 
// ============================================================================
template <class SCALAR> class TMatrixTSym ; // ROOT 
// ============================================================================
namespace Ostap { namespace Math { class  Statistic  ; } }
namespace Ostap { namespace Math { class WStatistic  ; } }
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  /** @class StatVar Ostap/StatVar.h
   *  Helper class to get statistical 
   *  infomation  about the variable/expression 
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date   2013-10-13
   */
  class StatVar 
  {
  public:
    // ========================================================================
    /// the actual type for statistic 
    typedef Ostap::WStatEntity        Statistic   ;
    /// the actual type for vector of statistic 
    typedef std::vector<Statistic>    Statistics  ;
    /// variable names 
    typedef std::vector<std::string>  Names       ;
    /// covariance
    typedef Ostap::Math::Covariance   Covariance  ; 
    /// (weighted) covariance
    typedef Ostap::Math::WCovariance  WCovariance ; 
    // ========================================================================
  public:
    // ========================================================================
    /** @struct Interval 
     *  the actual type for interval
     */
    struct Interval
    {
      Interval
      (  const double l = 0 ,
	 const double h = 0 )
        : low  ( std::min ( l , h ) )
        , high ( std::max ( l , h ) )
      {}
      // ======================================================================
      /// low edge of the interval 
      double low  { 0 } ; // low edge of the interval 
      /// high edge of the interval 
      double high { 0 } ; // high edge of the interval 
      // ======================================================================
    };
    // ========================================================================
    /** @struct Quantile 
     *  The actual type for quantile and statitsics 
     */    
    struct Quantile 
    {
      Quantile
      ( const double        q = 0 ,
	const unsigned long n = 0 )
        : quantile ( q ) 
        , nevents  ( n ) 
      {}
      // ======================================================================
      /// quantile value  
      double        quantile { 0 } ; // quantile value  
      /// number of events used for estimation
      unsigned long nevents  { 0 } ; // number of events used for estimation
      // ======================================================================
    };
    // ========================================================================
    /** @struct Quantiles 
     *  The actual type for quantiles and statitsics 
     */    
    struct Quantiles 
    {
      // =======================================================================
      Quantiles
      ( const std::vector<double>& q = std::vector<double>() , 
	const unsigned long        n = 0  )
        : quantiles ( q ) 
        , nevents   ( n ) 
      {}
      // ======================================================================
      /// quantile values  
      std::vector<double> quantiles {} ; // quantile values  
      /// number of events used for estimation
      unsigned long       nevents  { 0 } ; // number of events used for estimation
      // ======================================================================
    };
    // ========================================================================
    /** @struct QInterval
     *  The actual type for interval with statitsics 
     */    
    struct QInterval    
    {
      // =====================================================================
      QInterval
      ( const Interval&     i = Interval () ,
	const unsigned long n = 0           )
        : interval ( i ) 
        , nevents  ( n ) 
      {}
      // ======================================================================
      /// the interval
      Interval      interval {}    ; // the interval
      /// number of events used for estimation
      unsigned long nevents  { 0 } ; // number of events used for estimation
      // ======================================================================
    };
    // ========================================================================
  public: 
    // ========================================================================
    static_assert ( std::numeric_limits<unsigned long>::is_specialized   ,
                    "Numeric_limist<unsigned long> are not specialized!" ) ;
    static constexpr unsigned long LAST { std::numeric_limits<unsigned long>::max() } ;
    // ========================================================================
  public:
    // ========================================================================
    /** build statistic for the <code>expression</code>
     *  @param tree       (INPUT) the tree 
     *  @param expression (INPUT) the expression
     *  @param cuts       (INPUT) the selection criteria 
     *  @param first      (INPUT) the first entry 
     *  @param last       (INPUT) the last entry      *
     *
     *  @code
     *  tree = ... 
     *  stat = tree.statVar( 'S_sw' ) 
     *  @endcode 
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2013-10-13
     */
    static Ostap::StatEntity statVar
    ( TTree*              tree              , 
      const std::string&  expression        , 
      const unsigned long first      = 0    ,
      const unsigned long last       = LAST ) ;
    // ========================================================================
    /** build statistic for the <code>expression</code>
     *  @param tree       (INPUT) the tree 
     *  @param expression (INPUT) the expression
     *  @param cuts       (INPUT) the selection criteria 
     *  @param first      (INPUT) the first entry 
     *  @param last       (INPUT) the last entry      *
     *
     *  @code
     *  tree = ... 
     *  stat = tree.statVar( 'S_sw' ,'pt>1000') 
     *  @endcode 
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2013-10-13
     */
    static Statistic statVar
    ( TTree*              tree              , 
      const std::string&  expression        , 
      const std::string&  cuts              ,
      const unsigned long first      = 0    ,
      const unsigned long last       = LAST ) ;
    // ========================================================================
  public:
    // ========================================================================
    /** check if there is at least one entry that satisfies criteria 
     *  @param tree       (INPUT) the tree 
     *  @param cuts       (INPUT) criteria 
     *  @param first      (INPUT) the first entry 
     *  @param last       (INPUT) the last entry
     *  @return true if there exist at leats one entry 
     */
    static bool hasEntry 
    ( TTree*              tree              , 
      const std::string&  cuts              ,   
      const unsigned long first      = 0    ,
      const unsigned long last       = LAST ) ;      
    // ========================================================================
    /** check if there is at least one entry that satisfies criteria 
     *  @param data       (INPUT) data 
     *  @param cuts       (INPUT) criteria 
     *  @param first      (INPUT) the first entry 
     *  @param last       (INPUT) the last entry
     *  @return true if there exist at leats one entry 
     */
    static bool hasEntry 
    ( const RooAbsData*   data               , 
      const std::string&  cuts               , 
      const unsigned long first      = 0     ,
      const unsigned long last       = LAST  ) ;
    // ========================================================================
    /** check if there is at least one entry that satisfies criteria 
     *  @param data       (INPUT) data 
     *  @param cuts       (INPUT) criteria 
     *  @param cut_range  (INPUT) cut range 
     *  @param first      (INPUT) the first entry 
     *  @param last       (INPUT) the last entry
     *  @return true if there exist at leats one entry 
     */
    static bool hasEntry 
    ( const RooAbsData*   data               , 
      const std::string&  cuts               , 
      const std::string&  cut_range          , 
      const unsigned long first      = 0     ,
      const unsigned long last       = LAST  ) ;
    // ========================================================================
  public:
    // ========================================================================
    /** build statistic for the <code>expressions</code>
     *  @param tree        (INPUT)  the tree 
     *  @param result      (UPDATE) the output statistics for specified expressions 
     *  @param expressions (INPUT)  the list of  expressions
     *  @param first       (INPUT)  the first entry to process 
     *  @param last        (INPUT)  the last entry to process (not including!)
     *  @return number of processed entries 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2018-11-04
     */
    static unsigned long statVars
    ( TTree*                          tree              , 
      std::vector<Statistic>&         result            , 
      const Names&                    expressions       ,
      const unsigned long             first      = 0    ,
      const unsigned long             last       = LAST ) ;
    // ========================================================================
    /** build statistic for the <code>expressions</code>
     *  @param tree        (INPUT)  the tree 
     *  @param result      (UPDATE) the output statistics for specified expressions 
     *  @param expressions (INPUT)  the list of  expressions
     *  @param cuts        (INPUT)  the selection criteria 
     *  @param first       (INPUT)  the first entry to process 
     *  @param last        (INPUT)  the last entry to process (not including!)
     *  @return number of processed entries 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2018-11-04
     */
    static unsigned long statVars
    ( TTree*                          tree               ,       
      std::vector<Statistic>&         result             , 
      const Names&                    expressions        ,
      const std::string&              cuts               ,
      const unsigned long             first       = 0    ,
      const unsigned long             last        = LAST ) ;
    // ========================================================================    
  public:
    // ========================================================================
    /** build statistic for the <code>expression</code>
     *  @param data       (INPUT) the data 
     *  @param expression (INPUT) the expression
     *  @param first      (INPUT) the frist event to process
     *  @param last       (INPUT) process till "last"-event
     *
     *  @code
     *  data = ... 
     *  stat = data.statVar( 'S_sw' , 'pt>10') 
     *  @endcode 
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2015-02-15
     */
    static Statistic statVar 
    ( const RooAbsData*   data               , 
      const std::string&  expression         , 
      const unsigned long first      = 0     ,
      const unsigned long last       = LAST  ) 
    { return statVar( data , expression , std::string() , "" , first , last ) ; }
    // ========================================================================
    /** build statistic for the <code>expression</code>
     *  @param data       (INPUT) the data 
     *  @param expression (INPUT) the expression
     *  @param cuts       (INPUT) the selection
     *  @param first      (INPUT) the frist event to process
     *  @param last       (INPUT) process till "last"-event
     *
     *  @code
     *  data = ... 
     *  stat = data.statVar( 'S_sw' , 'pt>10') 
     *  @endcode 
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2015-02-15
     */
    static Statistic statVar 
    ( const RooAbsData*   data               , 
      const std::string&  expression         , 
      const std::string&  cuts               , 
      const unsigned long first      = 0     ,
      const unsigned long last       = LAST  ) 
    { return statVar( data , expression , cuts , "" , first , last ) ; }
    // ========================================================================
    /** build statistic for the <code>expression</code>
     *  @param data       (INPUT) the data 
     *  @param expression (INPUT) the expression
     *  @param cuts       (INPUT) the selection
     *  @param cutrange   (INPUT) the cut-range
     *  @param first      (INPUT) the frist event to process
     *  @param last       (INPUT) process till "last"-event
     *
     *  @code
     *  data = ... 
     *  cut  = TCut ( ... ) 
     *  stat = data.statVar( 'S_sw' , cut ) 
     *  @endcode 
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2015-02-15
     */
    static Statistic statVar 
    ( const RooAbsData*   data              , 
      const std::string&  expression        , 
      const std::string&  cuts              , 
      const std::string&  cut_range         ,
      const unsigned long first      = 0    ,
      const unsigned long last       = LAST ) ;
    // ========================================================================
  public:
    // ========================================================================
    /** build statistic for the <code>expressions</code>
     *  @param data        (INPUT)  input data 
     *  @param result      (UPDATE) the output statistics for specified expressions 
     *  @param expressions (INPUT)  the list of  expressions
     *  @param first       (INPUT)  the first entry to process 
     *  @param last        (INPUT)  the last entry to process (not including!)
     *  @return number of processed entries 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2018-11-04
     */
    static unsigned long statVars
    ( const RooAbsData*               data              , 
      std::vector<Statistic>&         result            , 
      const Names&                    expressions       ,
      const unsigned long             first      = 0    ,
      const unsigned long             last       = LAST ) 
    { return statVars ( data , result , expressions , std::string() , std::string () , first , last ) ; }
    // ========================================================================
    /** build statistic for the <code>expressions</code>
     *  @param data        (INPUT)  input data 
     *  @param result      (UPDATE) the output statistics for specified expressions 
     *  @param expressions (INPUT)  the list of  expressions
     *  @param cuts        (INPUT)  the selection  
     *  @param first       (INPUT)  the first entry to process 
     *  @param last        (INPUT)  the last entry to process (not including!)
     *  @return number of processed entries 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2021-06-04
     */
    static unsigned long statVars
    ( const RooAbsData*               data              , 
      std::vector<Statistic>&         result            , 
      const Names&                    expressions       ,
      const std::string&              cuts              ,
      const unsigned long             first      = 0    ,
      const unsigned long             last       = LAST ) 
    { return statVars ( data , result , expressions , cuts  , std::string() , first , last ) ; }
    // ========================================================================
    /** build statistic for the <code>expressions</code>
     *  @param data        (INPUT)  input data 
     *  @param result      (UPDATE) the output statistics for specified expressions 
     *  @param expressions (INPUT)  the list of  expressions
     *  @param cuts        (INPUT)  the selection 
     *  @param cut_range   (INPUT)  cut range  
     *  @param first       (INPUT)  the first entry to process 
     *  @param last        (INPUT)  the last entry to process (not including!)
     *  @return number of processed entries 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2021-06-04
     */
    static unsigned long statVars
    ( const RooAbsData*               data              , 
      std::vector<Statistic>&         result            , 
      const std::vector<std::string>& expressions       ,
      const std::string&              cuts              ,
      const std::string&              cut_range         ,
      const unsigned long             first      = 0    ,
      const unsigned long             last       = LAST ) ;
    // ========================================================================
  public:
    // ========================================================================
    /** build statistic for the <code>expression</code>
     *  @param data       (INPUT) the data 
     *  @param expression (INPUT) the expression
     *  @param cuts       (INPUT) the selection/weight
     *
     *  @code
     *  data = ... 
     *  stat = data.statVar( 'S_sw' , 'pt>10') 
     *  @endcode 
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2018-06-18
     */
    static Statistic statVar 
    ( FrameNode           frame           , 
      const std::string&  expression      , 
      const std::string&  cuts       = "" ) ;
    // ========================================================================    
  public:
    // ========================================================================
    /** get the number of equivalent entries 
     *  \f$ n_{eff} \equiv = \frac{ (\sum w)^2}{ \sum w^2} \f$
     *  @param tree  (INPUT) the tree 
     *  @param cuts  (INPUT) selection  criteria 
     *  @param  first  (INPUT) the first  event to process 
     *  @param  last   (INPUT) the last event to  process
     *  @return number of equivalent entries 
     */
    static double nEff 
    ( TTree&               tree          , 
      const std::string&   cuts   = ""   , 
      const unsigned long  first  = 0    ,
      const unsigned long  last   = LAST ) ;
    // ========================================================================
    /** get the number of equivalent entries 
     *  \f$ n_{eff} \equiv = \frac{ (\sum w)^2}{ \sum w^2} \f$
     *  @param tree      (INPUT) the tree 
     *  @param cuts      (INPUT) selection crietria  
     *  @param first     (INPUT) the first  event to process 
     *  @param last      (INPUT) the last event to  process
     *  @param cut_range (INPUT) cut range 
     *  @return number of equivalent entries 
     */
    static double nEff 
    ( const RooAbsData&    tree             , 
      const std::string&   cuts      = ""   , 
      const std::string&   cut_range = ""   , 
      const unsigned long  first     = 0    ,
      const unsigned long  last      = LAST ) ;
    // ========================================================================
    /** get the number of equivalent entries 
     *  \f$ n_{eff} \equiv = \frac{ (\sum w)^2}{ \sum w^2} \f$
     *  @param tree      (INPUT) the frame 
     *  @param cuts      (INPUT) selection crietria/weight  
     *  @param first     (INPUT) the first  event to process 
     *  @param last      (INPUT) the last event to  process
     *  @return number of equivalent entries 
     */
    static double nEff 
    ( FrameNode            frame            ,
      const std::string&   cuts      = ""   ) ;
    // ========================================================================
  public: 
    // ========================================================================
    /** calculate the covariance of several expressions 
     *  @param tree  (INPUT)  the inpout tree 
     *  @param vars  (INPUT)  expressions 
     *  @param cuts  (INPUT)  the selection criteria 
     *  @param stats (UPDATE) the statistics 
     *  @param cov2  (UPDATE) the covariance matrix 
     *  @return number of processed events 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2023-02-28
     */
    static unsigned long statCov 
    ( TTree*                          tree         ,
      const std::vector<std::string>& vars         , 
      const std::string&              cuts         ,
      std::vector<Statistic>&         stats        ,  
      TMatrixTSym<double>&            cov2         , 
      const unsigned long             first = 0    ,
      const unsigned long             last  = LAST ) ;
    // ========================================================================
 public:  // new methods 
    // ========================================================================
    /** calculate the covariance of two expressions 
     *  @param tree  (INPUT)  the inpout tree 
     *  @param exp1  (INPUT)  the first  expression
     *  @param exp2  (INPUT)  the second expression
     *  @return Covariance 
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
     *  @date   2024-07-22-
     */
    static Covariance statCov 
    ( TTree*               tree          , 
      const std::string&   exp1          , 
      const std::string&   exp2          , 
      const unsigned long  first  = 0    ,
      const unsigned long  last   = LAST ) ;
    // ========================================================================
    /** calculate the covariance of two expressions 
     *  @param tree  (INPUT)  the inpout tree 
     *  @param exp1  (INPUT)  the first  expression
     *  @param exp2  (INPUT)  the second expression
     *  @param cuts  (INPUT)  cuts/weigths expression 
     *  @return Covariance 
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
     *  @date   2024-07-22
     */
    static WCovariance statCov 
    ( TTree*               tree          , 
      const std::string&   exp1          , 
      const std::string&   exp2          ,
      const std::string&   cuts          ,      
      const unsigned long  first  = 0    ,
      const unsigned long  last   = LAST ) ;
    // ========================================================================
  public: 
    // ========================================================================
    /** calculate the covariance of several expressions 
     *  @param tree  (INPUT)  the inpout tree 
     *  @param vars  (INPUT)  expressions 
     *  @param stats (UPDATE) the statistics 
     *  @param cov2  (UPDATE) the covariance matrix 
     *  @return number of processed events 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2023-02-28
     */
    static unsigned long statCov 
    ( TTree*                          tree         ,
      const std::vector<std::string>& vars         , 
      std::vector<Statistic>&         stats        ,  
      TMatrixTSym<double>&            cov2         , 
      const unsigned long             first = 0    ,
      const unsigned long             last  = LAST ) ;
    // ========================================================================
  public:
    // ========================================================================
    /** calculate the covariance of several expressions 
     *  @param tree      (INPUT)  the inpout tree 
     *  @param vars      (INPUT)  expressions 
     *  @param cuts      (INPUT)  the selection criteria 
     *  @param stats     (UPDATE) the statistics 
     *  @param cov2      (UPDATE) the covariance matrix 
     *  @param cut_range (INPUT)  range  
     *  @return number of processed events 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2023-02-28
     */
    static unsigned long statCov 
    ( const RooAbsData*    tree             , 
      const std::vector<std::string>& vars           ,  
      const std::string&              cuts           ,
      std::vector<Statistic>&         stats          ,  
      TMatrixTSym<double>&            cov2           , 
      const std::string&              cut_range = "" ,
      const unsigned long             first = 0      ,
      const unsigned long             last  = LAST   ) ;
    // ========================================================================
    /** calculate the covariance of several expressions 
     *  @param tree      (INPUT)  the inpout tree 
     *  @param vars      (INPUT)  expressions 
     *  @param stats     (UPDATE) the statistics 
     *  @param cov2      (UPDATE) the covariance matrix 
     *  @param cut_range (INPUT)  range  
     *  @return number of processed events 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2023-02-28
     */
    static unsigned long statCov 
    ( const RooAbsData*               tree           , 
      const std::vector<std::string>& vars           ,  
      std::vector<Statistic>&         stats          ,  
      TMatrixTSym<double>&            cov2           , 
      const std::string&              cut_range = "" ,
      const unsigned long             first = 0      ,
      const unsigned long             last  = LAST   ) ;
    // ========================================================================
  public:
    // ========================================================================
    /** calculate the covariance of two expressions 
     *  @param tree  (INPUT)  the input  tree 
     *  @param exp1  (INPUT)  the first  expresiion
     *  @param exp2  (INPUT)  the second expresiion
     *  @param cuts  (INPUT)  selection/weight expression 
     *  @return number of processed events 
     *  
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2014-03-27
     */
    static WCovariance statCov
    ( const RooAbsData*    tree             , 
      const std::string&   exp1             , 
      const std::string&   exp2             , 
      const std::string&   cuts             , 
      const std::string&   cut_range = ""   ,
      const unsigned long  first     = 0    ,
      const unsigned long  last      = LAST ) ;
    // ========================================================================
  public:
    // ========================================================================
    /** calculate the covariance of two expressions 
     *  @param frame (INPUT)  data frame 
     *  @param exp1  (INPUT)  the first  expresiion
     *  @param exp2  (INPUT)  the second expresiion
     *  @return Covariancce
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2024-07-22
     */
    static Covariance statCov
    ( FrameNode            data , 
      const std::string&   exp1 , 
      const std::string&   exp2 ) ; 
    // ========================================================================
    /** calculate the covariance of two expressions 
     *  @param frame (INPUT)  data frame 
     *  @param exp1  (INPUT)  the first  expresiion
     *  @param exp2  (INPUT)  the second expresiion
     *  @param cuts  (INPUT)  selection/weight expresiion
     *  @return Covariancce
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2024-07-22
     */
    static WCovariance statCov
    ( FrameNode            data , 
      const std::string&   exp1 , 
      const std::string&   exp2 , 
      const std::string&   cuts ) ; 
    // ========================================================================
  public:
    // ========================================================================
    /** calculate the moment of order "order" relative to the center "center"
     *  @param  tree   (INPUT) input tree 
     *  @param  expr   (INPUT) expression  (must  be valid TFormula!)
     *  @param  order  (INPUT) the order 
     *  @param  center (INPUT) the center 
     *  @param  cuts   (INPUT) cuts 
     *  @param  first  (INPUT) the first  event to process 
     *  @param  last   (INPUT) the last event to  process
     *  @return the moment 
     */
    static double get_moment 
    ( TTree&               tree          ,  
      const unsigned short order         , 
      const std::string&   expr          , 
      const double         center = 0    ,
      const std::string&   cuts   = ""   , 
      const unsigned long  first  = 0    ,
      const unsigned long  last   = LAST ) ;
    // ========================================================================    
    /** calculate the moment of order "order" relative to the center "center"
     *  @param  data      (INPUT) input data
     *  @param  expr      (INPUT) expression  (must  be valid TFormula!)
     *  @param  order     (INPUT) the order 
     *  @param  center    (INPUT) the center 
     *  @param  cuts      (INPUT) cuts 
     *  @param  first     (INPUT) the first  event to process 
     *  @param  last      (INPUT) the last event to  process
     *  @param  cut_range (INPUT) cut range 
     *  @return the moment 
     */
    static double get_moment 
    ( const RooAbsData&    data             ,  
      const unsigned short order            , 
      const std::string&   expr             , 
      const double         center    = 0    ,
      const std::string&   cuts      = ""   , 
      const std::string&   cut_range = ""   , 
      const unsigned long  first     = 0    ,
      const unsigned long  last      = LAST ) ;
    // ========================================================================    
  public:
    // ========================================================================    
    /** calculate the moment of order "order" relative to the center "center"
     *  @param  frame  (INPUT) input frame 
     *  @param  expr   (INPUT) expression 
     *  @param  order  (INPUT) the order 
     *  @param  center (INPUT) the center 
     *  @param  cuts   (INPUT) the selection criteria/weight 
     *  @return the moment 
     */
    static double get_moment 
    ( FrameNode            frame       ,  
      const unsigned short order       , 
      const std::string&   expr        , 
      const double         center      ,    
      const std::string&   cuts   = "" ) ;
    // ========================================================================
  public:
    // ========================================================================    
    /** calculate the moment of order "order"
     *  @param  tree  (INPUT) input tree 
     *  @param  order (INPUT) the order 
     *  @param  expr  (INPUT) expression  (must  be valid TFormula!)
     *  @param  cuts  (INPUT) cuts 
     *  @param  first (INPUT) the first  event to process 
     *  @param  last  (INPUT) the last event to  process
     *  @return the moment 
     */ 
    static Ostap::Math::ValueWithError moment
    ( TTree&               tree         ,  
      const unsigned short order        ,
      const std::string&   expr         , 
      const std::string&   cuts  = ""   , 
      const unsigned long  first = 0    ,
      const unsigned long  last  = LAST ) ;
    // ========================================================================
    /** calculate the moment of order "order"
     *  @param  data  (INPUT) input data
     *  @param  order (INPUT) the order 
     *  @param  expr  (INPUT) expression  (must  be valid TFormula!)
     *  @param  cuts  (INPUT) cuts 
     *  @param  first (INPUT) the first  event to process 
     *  @param  last  (INPUT) the last event to  process
     *  @param  cut_range (INPUT) cut range 
     *  @return the moment 
     */ 
    static Ostap::Math::ValueWithError moment
    ( const  RooAbsData&   data             ,  
      const unsigned short order            ,
      const std::string&   expr             , 
      const std::string&   cuts      = ""   , 
      const std::string&   cut_range = ""   , 
      const unsigned long  first     = 0    ,
      const unsigned long  last      = LAST ) ;
    // ========================================================================
  public:
    // ========================================================================
    /** calculate the moment of order "order"
     *  @param  tree  (INPUT) input tree 
     *  @param  order (INPUT) the order 
     *  @param  expr  (INPUT) expression  (must  be valid TFormula!)
     *  @param  cuts  (INPUT) the selection criteria/weight
     *  @return the moment 
     */ 
    static Ostap::Math::ValueWithError moment
    ( FrameNode            frame      ,  
      const unsigned short order      ,
      const std::string&   expr       , 
      const std::string&   cuts  = "" ) ;
    // ========================================================================
  public:
    // ========================================================================
    /** calculate the central moment of order "order"
     *  @param  tree  (INPUT) input tree 
     *  @param  order (INPUT) the order 
     *  @param  expr  (INPUT) expression  (must  be valid TFormula!)
     *  @param  cuts  (INPUT) cuts 
     *  @param  first (INPUT) the first  event to process 
     *  @param  last  (INPUT) the last event to  process
     *  @return the moment 
     */ 
    static Ostap::Math::ValueWithError central_moment 
    ( TTree&               tree         ,  
      const unsigned short order        ,
      const std::string&   expr         , 
      const std::string&   cuts  = ""   , 
      const unsigned long  first = 0    ,
      const unsigned long  last  = LAST ) ;
    // ========================================================================
    /** calculate the central moment of order "order"
     *  @param  data  (INPUT) input data
     *  @param  order (INPUT) the order 
     *  @param  expr  (INPUT) expression  (must  be valid TFormula!)
     *  @param  cuts  (INPUT) cuts 
     *  @param  first (INPUT) the first  event to process 
     *  @param  last  (INPUT) the last event to  process
     *  @param  cut_range (INPUT) cut range 
     *  @return the moment 
     */ 
    static Ostap::Math::ValueWithError central_moment 
    ( const RooAbsData&    data             ,  
      const unsigned short order            ,
      const std::string&   expr             , 
      const std::string&   cuts      = ""   , 
      const std::string&   cut_range = ""   , 
      const unsigned long  first     = 0    ,
      const unsigned long  last      = LAST ) ;
    // ========================================================================    
  public:
    // ========================================================================    
    /** calculate the central moment of order "order"
     *  @param  tree  (INPUT) input tree 
     *  @param  order (INPUT) the order 
     *  @param  expr  (INPUT) expression 
     *  @param  cuts  (INPUT) the selection criteria/weight
     *  @return the moment 
     */ 
    static Ostap::Math::ValueWithError central_moment 
    ( FrameNode            frame    ,  
      const unsigned short order    ,
      const std::string&   expr     , 
      const std::string&   cuts     )  ;
    // ========================================================================    
  public:
    // ========================================================================    
    /** calculate the skewness of the  distribution
     *  @param  tree  (INPUT) input tree 
     *  @param  expr  (INPUT) expression  (must  be valid TFormula!)
     *  @param  cuts  (INPUT) cuts 
     *  @param  first (INPUT) the first  event to process 
     *  @param  last  (INPUT) the last event to  process
     *  @return the skewness 
     */
    static Ostap::Math::ValueWithError skewness 
    ( TTree&               tree         ,  
      const std::string&   expr         , 
      const std::string&   cuts  = ""   , 
      const unsigned long  first = 0    ,
      const unsigned long  last  = LAST ) ;
    // ========================================================================    
    /** calculate the skewness of the  distribution
     *  @param  data      (INPUT) input data
     *  @param  expr      (INPUT) expression  (must  be valid TFormula!)
     *  @param  cuts      (INPUT) cuts 
     *  @param  first     (INPUT) the first  event to process 
     *  @param  last      (INPUT) the last event to  process
     *  @param  cut_range (INPUT) cut range 
     *  @return the skewness 
     */
    static Ostap::Math::ValueWithError skewness 
    ( const RooAbsData&    data             ,  
      const std::string&   expr             , 
      const std::string&   cuts      = ""   , 
      const std::string&   cut_range = ""   , 
      const unsigned long  first     = 0    ,
      const unsigned long  last      = LAST ) ;
    // ========================================================================    
  public:
    // ========================================================================    
    /** calculate the skewness of the  distribution
     *  @param  frame (INPUT) frame
     *  @param  expr  (INPUT) expression
     *  @param  cuts  (INPUT) the selection criteria/weight 
     *  @return the skewness 
     */
    static Ostap::Math::ValueWithError skewness 
    ( FrameNode            frame        ,  
      const std::string&   expr         , 
      const std::string&   cuts  = ""   ) ;
    // ========================================================================    
  public:
    // ========================================================================    
    /** calculate the (excess) kurtosis of the  distribution
     *  @param  tree  (INPUT) input tree 
     *  @param  expr  (INPUT) expression  (must  be valid TFormula!)
     *  @param  cuts  (INPUT) cuts 
     *  @param  first (INPUT) the first  event to process 
     *  @param  last  (INPUT) the last event to  process
     *  @return the (excess) kurtosis
     */
    static Ostap::Math::ValueWithError kurtosis
    ( TTree&               tree         ,  
      const std::string&   expr         , 
      const std::string&   cuts  = ""   , 
      const unsigned long  first = 0    ,
      const unsigned long  last  = LAST ) ;
    // ========================================================================    
    /** calculate the (excess) kurtosis of the  distribution
     *  @param  data  (INPUT) input data
     *  @param  expr  (INPUT) expression  (must  be valid TFormula!)
     *  @param  cuts  (INPUT) cuts 
     *  @param  first (INPUT) the first  event to process 
     *  @param  last  (INPUT) the last event to  process
     *  @param  cut_range (INPUT) cut range 
     *  @return the (excess) kurtosis
     */
    static Ostap::Math::ValueWithError kurtosis
    ( const RooAbsData&    data             ,  
      const std::string&   expr             , 
      const std::string&   cuts      = ""   , 
      const std::string&   cut_range = ""   , 
      const unsigned long  first     = 0    ,
      const unsigned long  last      = LAST ) ;
    // ========================================================================    
  public:
    // ========================================================================    
    /** calculate the (excess) kurtosis of the  distribution
     *  @param  frame (INPUT) input frame 
     *  @param  expr  (INPUT) expression  (must  be valid TFormula!)
     *  @param  cuts  (INPUT) the selection criteria/weight 
     *  @param  first (INPUT) the first  event to process 
     *  @param  last  (INPUT) the last event to  process
     *  @return the (excess) kurtosis
     */
    static Ostap::Math::ValueWithError kurtosis
    ( FrameNode            frame        ,  
      const std::string&   expr         , 
      const std::string&   cuts  = ""   );
    // ========================================================================    
  public:
    // ========================================================================    
    /**  get quantile of the distribution  
     *   @param tree  (INPUT) the input tree 
     *   @param q     (INPUT) quantile value   0 < q < 1  
     *   @param expr  (INPUT) the expression 
     *   @param cuts  (INPUT) selection cuts 
     *   @param  first (INPUT) the first  event to process 
     *   @param  last  (INPUT) the last event to  process
     *   @return the quantile value 
     */
    static Quantile quantile
    ( TTree&              tree         ,
      const double        q            , //  0<q<1 
      const std::string&  expr         , 
      const std::string&  cuts  = ""   , 
      const unsigned long first = 0    ,
      const unsigned long last  = LAST ) ;
    // ========================================================================    
    /**  get (approximate) quantile of the distribution  using P^2 algortihm
     *   @param tree  (INPUT) the input tree 
     *   @param q     (INPUT) quantile value   0 < q < 1  
     *   @param expr  (INPUT) the expression 
     *   @param cuts  (INPUT) selection cuts 
     *   @param  first (INPUT) the first  event to process 
     *   @param  last  (INPUT) the last event to  process
     *   @return the quantile value 
     */
    static Quantile p2quantile
    ( TTree&              tree         ,
      const double        q            , //  0<q<1 
      const std::string&  expr         , 
      const std::string&  cuts  = ""   , 
      const unsigned long first = 0    ,
      const unsigned long last  = LAST ) ;
    // ========================================================================    
  public:
    // ========================================================================
    /**  get quantile of the distribution  
     *   @param data   (INPUT) the input data
     *   @param q      (INPUT) quantile value   0 < q < 1  
     *   @param expr   (INPUT) the expression 
     *   @param cuts   (INPUT) selection cuts 
     *   @param  first (INPUT) the first  event to process 
     *   @param  last  (INPUT) the last event to  process
     *   @param  cut_range (INPUT) cut range 
     *   @return the quantile value 
     */
    static Quantile quantile
    ( const RooAbsData&   data             ,
      const double        q                , //  0<q<1 
      const std::string&  expr             , 
      const std::string&  cuts      = ""   , 
      const std::string&  cut_range = ""   , 
      const unsigned long first     = 0    ,
      const unsigned long last      = LAST ) ;
    // ========================================================================    
    /**  get (approximate) quantile of the distribution  using p^2 algorithm
     *   @param data   (INPUT) the input data
     *   @param q      (INPUT) quantile value   0 < q < 1  
     *   @param expr   (INPUT) the expression 
     *   @param cuts   (INPUT) selection cuts 
     *   @param  first (INPUT) the first  event to process 
     *   @param  last  (INPUT) the last event to  process
     *   @param  cut_range (INPUT) cut range 
     *   @return the quantile value 
     */
    static Quantile p2quantile
    ( const RooAbsData&   data             ,
      const double        q                , //  0<q<1 
      const std::string&  expr             , 
      const std::string&  cuts      = ""   , 
      const std::string&  cut_range = ""   , 
      const unsigned long first     = 0    ,
      const unsigned long last      = LAST ) ;
    // ========================================================================    
  public:
    // ========================================================================    
    /**  get quantile of the distribution  
     *   @param frame  (INPUT) the input data
     *   @param q      (INPUT) quantile value   0 < q < 1  
     *   @param expr   (INPUT) the expression 
     *   @param cuts   (INPUT) selection criteria
     *   @return the quantile value 
     */
    static Quantile quantile
    ( FrameNode           frame            ,
      const double        q                , //  0<q<1 
      const std::string&  expr             , 
      const std::string&  cuts      = ""   ) ;
    // ========================================================================    
    /**  get approximate quantile of the distribution usnig P^2 algorithm 
     *   @param frame  (INPUT) the input data
     *   @param q      (INPUT) quantile value   0 < q < 1  
     *   @param expr   (INPUT) the expression 
     *   @param cuts   (INPUT) selection criteria
     *   @return the quantile value 
     */
    static Quantile p2quantile
    ( FrameNode           frame            ,
      const double        q                , //  0<q<1 
      const std::string&  expr             , 
      const std::string&  cuts      = ""   ) ;
    // ========================================================================    
  public:
    // ========================================================================    
    /**  get quantiles of the distribution  
     *   @param tree  (INPUT) the input tree 
     *   @param q     (INPUT) quantile value   0 < q < 1  
     *   @param expr  (INPUT) the expression 
     *   @param cuts  (INPUT) selection cuts 
     *   @param  first (INPUT) the first  event to process 
     *   @param  last  (INPUT) the last event to  process
     *   @return the quantile value 
     */
    static Quantiles quantiles
    ( TTree&                     tree             ,
      const std::vector<double>& quantiles        , 
      const std::string&         expr             , 
      const std::string&         cuts      = ""   , 
      const unsigned long        first     = 0    ,
      const unsigned long        last      = LAST ) ;
    // ========================================================================
    /**  get quantiles of the distribution  
     *   @param data  (INPUT) the input data
     *   @param q     (INPUT) quantile value   0 < q < 1  
     *   @param expr  (INPUT) the expression 
     *   @param cuts  (INPUT) selection cuts 
     *   @param  first (INPUT) the first  event to process 
     *   @param  last  (INPUT) the last event to  process
     *   @return the quantile value 
     */
    static Quantiles quantiles
    ( const RooAbsData&          data              ,
      const std::vector<double>& quantiles         , 
      const std::string&         expr              , 
      const std::string&         cuts       = ""   , 
      const std::string&         cut_range  = ""   , 
      const unsigned long        first      = 0    ,
      const unsigned long        last       = LAST ) ;
    // ========================================================================
  public:
    // ========================================================================    
    /**  get (approximate) quantiles of the distribution  using P^2 algorithm
     *   @param tree  (INPUT) the input tree 
     *   @param q     (INPUT) quantile value   0 < q < 1  
     *   @param expr  (INPUT) the expression 
     *   @param cuts  (INPUT) selection cuts 
     *   @param  first (INPUT) the first  event to process 
     *   @param  last  (INPUT) the last event to  process
     *   @return the quantile value 
     */
    static Quantiles p2quantiles
    ( TTree&                     tree             ,
      const std::vector<double>& quantiles        , 
      const std::string&         expr             , 
      const std::string&         cuts      = ""   , 
      const unsigned long        first     = 0    ,
      const unsigned long        last      = LAST ) ;
    // ========================================================================
    /**  get (approximate) quantiles of the distribution using P^2 algorithm  
     *   @param data  (INPUT) the input data
     *   @param q     (INPUT) quantile value   0 < q < 1  
     *   @param expr  (INPUT) the expression 
     *   @param cuts  (INPUT) selection cuts 
     *   @param  first (INPUT) the first  event to process 
     *   @param  last  (INPUT) the last event to  process
     *   @return the quantile value 
     */
    static Quantiles p2quantiles
    ( const RooAbsData&          data              ,
      const std::vector<double>& quantiles         , 
      const std::string&         expr              , 
      const std::string&         cuts       = ""   , 
      const std::string&         cut_range  = ""   , 
      const unsigned long        first      = 0    ,
      const unsigned long        last       = LAST ) ;
    // ========================================================================
  public:
    // ========================================================================    
    /**  get quantiles of the distribution  
     *   @param frame (INPUT) the input frame
     *   @param q     (INPUT) quantile value   0 < q < 1  
     *   @param expr  (INPUT) the expression 
     *   @param cuts  (INPUT) selection cuts 
     *   @return the quantile value 
     */
    static Quantiles quantiles
    ( FrameNode                  frame            ,
      const std::vector<double>& quantiles        , 
      const std::string&         expr             , 
      const std::string&         cuts      = ""   ) ;
    // ========================================================================
  public:
    // ========================================================================    
    /**  get approximate  quantiles of the distribution  using P^2 algorithm
     *   @param frame (INPUT) the input frame
     *   @param q     (INPUT) quantile value   0 < q < 1  
     *   @param expr  (INPUT) the expression 
     *   @param cuts  (INPUT) selection cuts 
     *   @return the quantile value 
     */
    static Quantiles p2quantiles
    ( FrameNode                  frame            ,
      const std::vector<double>& quantiles        , 
      const std::string&         expr             , 
      const std::string&         cuts      = ""   ) ;
    // ========================================================================
  public:
    // ========================================================================    
    /**  get the interval of the distribution  
     *   @param tree  (INPUT) the input tree 
     *   @param q1    (INPUT) quantile value   0 < q1 < 1  
     *   @param q2    (INPUT) quantile value   0 < q2 < 1  
     *   @param expr  (INPUT) the expression 
     *   @param cuts  (INPUT) selection cuts 
     *   @param  first (INPUT) the first  event to process 
     *   @param  last  (INPUT) the last event to  process
     *   @return the quantile value 
     *   @code
     *   Tree& tree = ... ;
     *   /// get 90% interval:
     *   Interval ab = interval ( tree , 0.05 , 0.95 , 'mass' , 'pt>3' ) ;
     *   @endcode 
     */
    static QInterval interval 
    ( TTree&              tree         ,
      const double        q1           , //  0<q1<1 
      const double        q2           , //  0<q2<1 
      const std::string&  expr         , 
      const std::string&  cuts  = ""   , 
      const unsigned long first = 0    ,
      const unsigned long last  = LAST ) ;
    // ========================================================================    
    /**  get the interval of the distribution  
     *   @param data  (INPUT) the input data
     *   @param q1    (INPUT) quantile value   0 < q1 < 1  
     *   @param q2    (INPUT) quantile value   0 < q2 < 1  
     *   @param expr  (INPUT) the expression 
     *   @param cuts  (INPUT) selection cuts 
     *   @param  first (INPUT) the first  event to process 
     *   @param  last  (INPUT) the last event to  process
     *   @return the quantile value 
     *   @code
     *   const RooAbsData& data = ... ;
     *   /// get 90% interval:
     *   Interval ab = interval ( data , 0.05 , 0.95 , 'mass' , 'pt>3' ) ;
     *   @endcode 
     */
    static QInterval interval 
    ( const RooAbsData&   data             ,
      const double        q1               , //  0<q1<1 
      const double        q2               , //  0<q2<1 
      const std::string&  expr             , 
      const std::string&  cuts      = ""   , 
      const std::string&  cut_range = ""   , 
      const unsigned long first     = 0    ,
      const unsigned long last      = LAST ) ;
    // ========================================================================    
  public:
    // ========================================================================    
    /**  get the (approximate) interval of the distribution using P^2 algotirhm 
     *   @param tree  (INPUT) the input tree 
     *   @param q1    (INPUT) quantile value   0 < q1 < 1  
     *   @param q2    (INPUT) quantile value   0 < q2 < 1  
     *   @param expr  (INPUT) the expression 
     *   @param cuts  (INPUT) selection cuts 
     *   @param  first (INPUT) the first  event to process 
     *   @param  last  (INPUT) the last event to  process
     *   @return the quantile value 
     *   @code
     *   Tree& tree = ... ;
     *   /// get 90% interval:
     *   Interval ab = p2interval ( tree , 0.05 , 0.95 , 'mass' , 'pt>3' ) ;
     *   @endcode 
     */
    static QInterval p2interval 
    ( TTree&              tree         ,
      const double        q1           , //  0<q1<1 
      const double        q2           , //  0<q2<1 
      const std::string&  expr         , 
      const std::string&  cuts  = ""   , 
      const unsigned long first = 0    ,
      const unsigned long last  = LAST ) ;
    // ========================================================================    
    /**  get the approximate  interval of the distribution using P^2 algorithm 
     *   @param data  (INPUT) the input data
     *   @param q1    (INPUT) quantile value   0 < q1 < 1  
     *   @param q2    (INPUT) quantile value   0 < q2 < 1  
     *   @param expr  (INPUT) the expression 
     *   @param cuts  (INPUT) selection cuts 
     *   @param  first (INPUT) the first  event to process 
     *   @param  last  (INPUT) the last event to  process
     *   @return the quantile value 
     *   @code
     *   const RooAbsData& data = ... ;
     *   /// get 90% interval:
     *   Interval ab = p2interval ( data , 0.05 , 0.95 , 'mass' , 'pt>3' ) ;
     *   @endcode 
     */
    static QInterval p2interval 
    ( const RooAbsData&   data             ,
      const double        q1               , //  0<q1<1 
      const double        q2               , //  0<q2<1 
      const std::string&  expr             , 
      const std::string&  cuts      = ""   , 
      const std::string&  cut_range = ""   , 
      const unsigned long first     = 0    ,
      const unsigned long last      = LAST ) ;
    // ========================================================================    
  public:
    // ========================================================================    
    /**  get the interval of the distribution  
     *   @param tree  (INPUT) the input tree 
     *   @param q1    (INPUT) quantile value   0 < q1 < 1  
     *   @param q2    (INPUT) quantile value   0 < q2 < 1  
     *   @param expr  (INPUT) the expression 
     *   @param cuts  (INPUT) selection cuts 
     *   @return the quantile value 
     *   @code
     *   FRAME& frame = ... ;
     *   /// get 90% interval:
     *   Interval ab = interval ( frame , 0.05 , 0.95 , 'mass' , 'pt>3' ) ;
     *   @endcode 
     */
    static QInterval interval 
    ( FrameNode           frame        ,
      const double        q1           , //  0<q1<1 
      const double        q2           , //  0<q2<1 
      const std::string&  expr         , 
      const std::string&  cuts  = ""   ) ;
    // ========================================================================    
  public:
    // ========================================================================    
    /**  get the approximate  interval of the distribution  
     *   @param tree  (INPUT) the input tree 
     *   @param q1    (INPUT) quantile value   0 < q1 < 1  
     *   @param q2    (INPUT) quantile value   0 < q2 < 1  
     *   @param expr  (INPUT) the expression 
     *   @param cuts  (INPUT) selection cuts 
     *   @return the quantile value 
     *   @code
     *   FRAME& frame = ... ;
     *   /// get 90% interval:
     *   Interval ab = p2interval ( frame , 0.05 , 0.95 , 'mass' , 'pt>3' ) ;
     *   @endcode 
     */
    static QInterval p2interval 
    ( FrameNode           frame        ,
      const double        q1           , //  0<q1<1 
      const double        q2           , //  0<q2<1 
      const std::string&  expr         , 
      const std::string&  cuts  = ""   ) ;
    // ========================================================================    
  public:  
    // ========================================================================
    /** get the moment as Ostap::Math::Moment_<N>
     *  @see Ostap::Math::Moment_
     *  @see Ostap::Math::Moment
     */
    static Ostap::StatusCode the_moment
    ( TTree*                  tree          , 
      Ostap::Math::Statistic& moment        , 
      const std::string&      expression    , 
      const unsigned long     first  = 0    ,
      const unsigned long     last   = LAST ) ;
    // ========================================================================    
    /** get the moment as Ostap::Math::WMoment_<N>
     *  @see Ostap::Math::WMoment_
     *  @see Ostap::Math::WMoment
     */
    static Ostap::StatusCode the_moment
    ( TTree*                   tree          , 
      Ostap::Math::WStatistic& moment        , 
      const std::string&       expression    , 
      const std::string&       selection     , 
      const unsigned long      first  = 0    ,
      const unsigned long      last   = LAST ) ;
    // ========================================================================    
    /** get the moment as Ostap::Math::WMoment_<N>
     *  @see Ostap::Math::WMoment_
     *  @see Ostap::Math::WMoment
     */
    static Ostap::StatusCode the_moment
    ( const RooAbsData*        data , 
      Ostap::Math::WStatistic& moment              , 
      const std::string&       expression          , 
      const std::string&       selection = ""      , 
      const std::string&       cutrange  = ""      ,      
      const unsigned long      first     = 0       ,
      const unsigned long      last      = LAST    ) ;
    // ========================================================================
    /** get the moment as Ostap::Math::WMoment_<N>
     *  @see Ostap::Math::WMoment_
     *  @see Ostap::Math::WMoment
     */
    static Ostap::StatusCode the_moment
    ( const RooAbsData*        data                , 
      Ostap::Math::WStatistic& moment              , 
      const std::string&       expression          , 
      const std::string&       selection           , 
      const unsigned long      first               ,
      const unsigned long      last      = LAST    ) ;
    // ========================================================================
    /** get the moment as Ostap::Math::WMoment_<N>
     *  @see Ostap::Math::WMoment_
     *  @see Ostap::Math::WMoment
     */
    static Ostap::StatusCode the_moment
    ( const RooAbsData*        data                , 
      Ostap::Math::WStatistic& moment              , 
      const std::string&       expression          , 
      const unsigned long      first               ,
      const unsigned long      last      = LAST    ) ;
    // ========================================================================
  public:
    // ========================================================================
    /// data column 
    typedef std::vector<double> Column ;    
    /// data table 
    typedef std::vector<Column> Table  ;
    // ========================================================================    
  public:
    // ========================================================================    
    /** get variables from dataset in form of the table 
     *  @param data input dataset
     *  @param vars list of variables
     *  @param cuts selection criteria 
     *  @param table output table
     *  @param weights column of weigths (empty for non-weighted data) 
     *  @param first first entry 
     *  @param last  last entry 
     */
    static unsigned long 
    get_table 
    ( const RooAbsData*                  data         , 
      const Names&                       vars         , 
      const std::string&                 cuts         , 
      Table&                             table        ,
      Column&                            weights      , 
      const unsigned long                first = 0    ,
      const unsigned long                last  = LAST ) ;
    // ========================================================================
    /** get variables from dataset in form of the table 
     *  @param data input dataset
     *  @param vars list of variables
     *  @param table output table
     *  @param weights column of weigths (empty for non-weighted data) 
     *  @param first first entry 
     *  @param last  last entry 
     */
    static unsigned long 
    get_table 
    ( const RooAbsData*                  data         , 
      const Names&                       vars         , 
      Table&                             table        ,
      Column&                            weights      , 
      const unsigned long                first = 0    ,
      const unsigned long                last  = LAST ) ;    
    // ========================================================================
  public:
    // ========================================================================    
    /** get variables from dataset in form of the table 
     *  @param data input dataset
     *  @param vars list of variables
     *  @param cuts selection criteria 
     *  @param cutrange cut range 
     *  @param table output table
     *  @param weights column of weigths (empty for non-weighted data) 
     *  @param first first entry 
     *  @param last  last entry 
     */
    static unsigned long 
    get_table
    ( const RooAbsData*                  data         , 
      const Names&                       vars         , 
      const std::string&                 cuts         , 
      Table&                             table        ,
      Column&                            weights      , 
      const std::string&                 cutrange     ,
      const unsigned long                first = 0    ,
      const unsigned long                last  = LAST ) ;
    // ========================================================================
  } ;  
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_STATVAR_H
// ============================================================================
