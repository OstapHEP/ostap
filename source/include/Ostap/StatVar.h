// $Id$ 
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
class TCut       ; // ROOT 
class RooAbsData ; // RooFit 
// =============================================================================
// Ostap
// ============================================================================
#include "Ostap/WStatEntity.h"
#include "Ostap/SymmetricMatrixTypes.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  /** @class StatVar Ostap/StatVar.h
   *  Helper class to get statistical 
   *  infomation  about the variable/expression 
   *
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date   2013-10-13
   */
  class StatVar 
  {
  public:
    // ========================================================================
    /// the actual type for styatistic 
    typedef Ostap::WStatEntity  Statistic ;
    // ========================================================================
  public: 
    // ========================================================================
    /** build statistic for the <code>expression</code>
     *  @param tree (INPUT) the tree 
     *  @param expression (INPUT) the expression
     *
     *  @code
     *  tree = ... 
     *  stat = tree.statVar ( 'S_sw' ) 
     *  @endcode 
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2013-10-13
     */
    static Statistic statVar
      ( TTree*              tree                                                   , 
        const std::string&  expression                                             ,
        const unsigned long first      = 0                                         ,
        const unsigned long last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** build statistic for the <code>expression</code>
     *  @param tree       (INPUT) the tree 
     *  @param expression (INPUT) the expression
     *  @param cuts       (INPUT) the selection criteria 
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
      ( TTree*              tree        , 
        const std::string&  expression  , 
        const std::string&  cuts        ,
        const unsigned long first      = 0                                         ,
        const unsigned long last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** build statistic for the <code>expression</code>
     *  @param tree       (INPUT) the tree 
     *  @param expression (INPUT) the expression
     *  @param cuts       (INPUT) the selection criteria 
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
      ( TTree*              tree        , 
        const std::string&  expression  , 
        const TCut&         cuts        ,
        const unsigned long first      = 0                                         ,
        const unsigned long entries    = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** calculate the covariance of two expressions 
     *  @param tree  (INPUT)  the inpout tree 
     *  @param exp1  (INPUT)  the first  expresiion
     *  @param exp2  (INPUT)  the second expresiion
     *  @param stat1 (UPDATE) the statistic for the first  expression
     *  @param stat2 (UPDATE) the statistic for the second expression
     *  @param cov2  (UPDATE) the covariance matrix 
     *  @return number of processed events 
     *  
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2014-03-27
     */
    static unsigned long statCov 
      ( TTree*               tree   , 
        const std::string&   exp1   , 
        const std::string&   exp2   , 
        Statistic&           stat1  ,  
        Statistic&           stat2  ,  
        Ostap::SymMatrix2x2& cov2   , 
        const unsigned long  first   = 0 ,
        const unsigned long  last    = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** calculate the covariance of two expressions 
     *  @param tree  (INPUT)  the inpout tree 
     *  @param exp1  (INPUT)  the first  expresiion
     *  @param exp2  (INPUT)  the second expresiion
     *  @param cuts  (INPUT)  the selection criteria 
     *  @param stat1 (UPDATE) the statistic for the first  expression
     *  @param stat2 (UPDATE) the statistic for the second expression
     *  @param cov2  (UPDATE) the covariance matrix 
     *  @return number of processed events 
     *  
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2014-03-27
     */
    static unsigned long statCov 
      ( TTree*               tree   ,
        const std::string&   exp1   , 
        const std::string&   exp2   , 
        const std::string&   cuts   ,
        Statistic&           stat1  ,  
        Statistic&           stat2  ,  
        Ostap::SymMatrix2x2& cov2   , 
        const unsigned long  first   = 0                                         ,
        const unsigned long  last    = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** calculate the covariance of two expressions 
     *  @param tree  (INPUT)  the inpout tree 
     *  @param exp1  (INPUT)  the first  expresiion
     *  @param exp2  (INPUT)  the second expresiion
     *  @param cuts  (INPUT)  the selection criteria 
     *  @param stat1 (UPDATE) the statistic for the first  expression
     *  @param stat2 (UPDATE) the statistic for the second expression
     *  @param cov2  (UPDATE) the covariance matrix 
     *  @return number of processed events 
     *  
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2014-03-27
     */
    static unsigned long statCov 
      ( TTree*               tree   ,
        const std::string&   exp1   , 
        const std::string&   exp2   , 
        const TCut&          cuts   ,
        Statistic&           stat1  ,  
        Statistic&           stat2  ,  
        Ostap::SymMatrix2x2& cov2   , 
        const unsigned long  first   = 0                                           ,
        const unsigned long  last    = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
  public: // the same but with RooFit 
    // ========================================================================
    /** build statistic for the <code>expression</code>
     *  @param data (INPUT) the data 
     *  @param expression (INPUT) the expression
     *
     *  @code
     *  data = ... 
     *  stat = data.statVar( 'S_sw' ) 
     *  @endcode 
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2015-02-15
     */
    static Statistic statVar 
      ( const RooAbsData*   data        ,   
        const std::string&  expression  , 
        const unsigned long first   = 0                                           ,
        const unsigned long last    = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** build statistic for the <code>expression</code>
     *  @param data       (INPUT) the data 
     *  @param expression (INPUT) the expression
     *  @param cuts       (INPUT) the selection
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
      ( const RooAbsData*   data        , 
        const std::string&  expression  , 
        const std::string&  cuts        , 
        const unsigned long first   = 0                                         ,
        const unsigned long last    = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** build statistic for the <code>expression</code>
     *  @param data       (INPUT) the data 
     *  @param expression (INPUT) the expression
     *  @param cuts       (INPUT) the selection
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
      ( const RooAbsData*   data        , 
        const std::string&  expression  , 
        const TCut&         cuts        , 
        const unsigned long first   = 0                                         ,
        const unsigned long last    = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** calculate the covariance of two expressions 
     *  @param tree  (INPUT)  the input  tree 
     *  @param exp1  (INPUT)  the first  expresiion
     *  @param exp2  (INPUT)  the second expresiion
     *  @param stat1 (UPDATE) the statistic for the first  expression
     *  @param stat2 (UPDATE) the statistic for the second expression
     *  @param cov2  (UPDATE) the covariance matrix 
     *  @return number of processed events 
     *  
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2014-03-27
     */
    static unsigned long statCov
      ( const RooAbsData*    tree   , 
        const std::string&   exp1   , 
        const std::string&   exp2   , 
        Statistic&           stat1  ,  
        Statistic&           stat2  ,  
        Ostap::SymMatrix2x2& cov2   , 
        const unsigned long  first   = 0 ,
        const unsigned long  last    = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** calculate the covariance of two expressions 
     *  @param tree  (INPUT)  the input  tree 
     *  @param exp1  (INPUT)  the first  expresiion
     *  @param exp2  (INPUT)  the second expresiion
     *  @param cuts  (INPUT)  the selection 
     *  @param stat1 (UPDATE) the statistic for the first  expression
     *  @param stat2 (UPDATE) the statistic for the second expression
     *  @param cov2  (UPDATE) the covariance matrix 
     *  @return number of processed events 
     *  
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2014-03-27
     */
    static unsigned long statCov
      ( const RooAbsData*    tree   , 
        const std::string&   exp1   , 
        const std::string&   exp2   , 
        const std::string&   cuts   , 
        Statistic&           stat1  ,  
        Statistic&           stat2  ,  
        Ostap::SymMatrix2x2& cov2   , 
        const unsigned long  first   = 0 ,
        const unsigned long  last    = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
  public: // generic cases, require some decoration/postprocessing  in python
    // ========================================================================
    /** calculate the covariances for generic case 
     *  @param tree  (INPUT)  the inpout tree 
     *  @param vars  (INPUT)  the list of variables 
     *  @param cuts  (INPUT)  the selection criteria/weights 
     *  @param stats (UPDATE) list of statistics 
     *  @param covs  (UPDATE) elements of "covariance matrix" 
     *  @return number of processed events 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2017-02-17
     */
    static unsigned long _statCov 
    ( TTree*               tree   , 
      const std::vector<std::string>& vars      , 
      const std::string&              cuts      ,
      std::vector<Statistic>&         stats     ,  
      std::vector<double>&            covs      ,
      const unsigned long  first   = 0          ,
      const unsigned long  last    = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================    
    /** calculate the covariances for generic case 
     *  @param tree  (INPUT)  the inpout tree 
     *  @param vars  (INPUT)  the list of variables 
     *  @param cuts  (INPUT)  the selection criteria/weights 
     *  @param stats (UPDATE) list of statistics 
     *  @param covs  (UPDATE) elements of "covariance matrix" 
     *  @return number of processed events 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2017-02-17
     */
    static unsigned long _statCov 
    ( TTree*               tree   , 
      const std::vector<std::string>& vars      , 
      const TCut&                     cuts      ,
      std::vector<Statistic>&         stats     ,  
      std::vector<double>&            covs      ,
      const unsigned long  first   = 0          ,
      const unsigned long  last    = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** calculate the covariances for generic case 
     *  @param tree  (INPUT)  the inpout tree 
     *  @param vars  (INPUT)  the list of variables 
     *  @param stats (UPDATE) list of statistics 
     *  @param covs  (UPDATE) elements of "covariance matrix" 
     *  @return number of processed events 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2017-02-17
     */
    static unsigned long _statCov 
    ( TTree*               tree   , 
      const std::vector<std::string>& vars      , 
      std::vector<Statistic>&         stats     ,  
      std::vector<double>&            covs      ,
      const unsigned long  first   = 0          ,
      const unsigned long  last    = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
  public: // the same but with RooFit 
    // ========================================================================
    /** calculate the covariances for generic case 
     *  @param data  (INPUT)  the inpout dataset 
     *  @param vars  (INPUT)  the list of variables 
     *  @param cuts  (INPUT)  the selection criteria/weights 
     *  @param stats (UPDATE) list of statistics 
     *  @param covs  (UPDATE) elements of "covariance matrix" 
     *  @return number of processed events 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2017-02-17
     */
    static unsigned long _statCov 
    ( const RooAbsData*               data      , 
      const std::vector<std::string>& vars      , 
      const std::string&              cuts      ,
      std::vector<Statistic>&         stats     ,  
      std::vector<double>&            covs      ,
      const unsigned long  first   = 0          ,
      const unsigned long  last    = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================    
    /** calculate the covariances for generic case 
     *  @param data  (INPUT)  the inpout dataset 
     *  @param vars  (INPUT)  the list of variables 
     *  @param stats (UPDATE) list of statistics 
     *  @param covs  (UPDATE) elements of "covariance matrix" 
     *  @return number of processed events 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2017-02-17
     */
    static unsigned long _statCov 
    ( const RooAbsData*               data      , 
      const std::vector<std::string>& vars      , 
      std::vector<Statistic>&         stats     ,  
      std::vector<double>&            covs      ,
      const unsigned long  first   = 0          ,
      const unsigned long  last    = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
  };
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_STATVAR_H
// ============================================================================
