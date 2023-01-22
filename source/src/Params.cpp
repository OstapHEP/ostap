// ============================================================================
// Include files 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Params.h"
#include "Ostap/Formula.h"
#include "Ostap/Notifier.h"
#include "Ostap/ProgressBar.h"
#include "Ostap/Parameterization.h"
#include "Ostap/Polynomials.h"
#include "Ostap/Bernstein.h"
#include "Ostap/Bernstein2D.h"
#include "Ostap/Bernstein3D.h"
// ============================================================================
// ROOT
// ============================================================================
#include "TTree.h"
// ============================================================================
// Local 
// ============================================================================
#include "Exception.h"
// ============================================================================
/** @file 
 *  Implementation file for header Ostap/Params.h
 *  @date 2019-07-03 
 *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
 */
// ============================================================================
namespace
{
  // ==========================================================================  
  template <class OBJECT> 
  unsigned long 
  _param1_
  ( TTree*                            data       , 
    const Ostap::Utils::ProgressConf& progress   ,
    OBJECT&                           obj        ,
    const std::string&                expression ,
    const std::string&                selection  ,
    const unsigned long               first      ,
    const unsigned long               last       ,
    const double                      xmin       , 
    const double                      xmax       )
  {
    if ( 0 == data         ) { return 0 ; }
    //
    const unsigned long nEntries = std::min ( last , (unsigned long) data->GetEntries() ) ;
    if ( nEntries <= first  ) { return 0 ; }
    //
    Ostap::Formula xvar ( expression , data ) ;
    Ostap::Assert ( !(!xvar) && xvar.ok ()    , 
                    "Invalid expression"      ,   
                    "Ostap::DataParam::parameterize"  , 310 ) ;
    //
    std::unique_ptr<Ostap::Formula> cut { nullptr } ;
    if  ( !selection.empty() ) 
    { 
      cut = std::make_unique<Ostap::Formula>( selection , data ) ; 
      Ostap::Assert ( !(!cut) && cut->ok () , 
                      "Invalid selection"   ,
                      "Ostap::DataParam::parameterize" , 315 ) ;
    }
    //
    unsigned long filled = 0 ;
    //
    Ostap::Utils::Notifier notify ( data , &xvar , cut.get() ) ;
    std::vector<double> results {} ;
    /// make an explicit loop 
    Ostap::Utils::ProgressBar bar (  nEntries - first , progress ) ;
    for ( unsigned long entry = first ; entry < nEntries ; ++entry , ++bar )
    {
      long ievent = data->GetEntryNumber ( entry ) ;
      Ostap::Assert ( 0 <= ievent , 
                      "Errot in TTree::GetEventNumber" , 
                      "Ostap::DataParam::parameterize" , 316 ) ;
      //
      ievent      = data->LoadTree ( ievent ) ;      
      Ostap::Assert ( 0 <= ievent , 
                      "Error in TTree::LoadTree " , 
                      "Ostap::DataParam::parameterize" , 317 ) ;
      //
      const double weight = cut ? cut -> evaluate() : 1.0 ;
      if ( !weight    ) { continue ;}
      //
      xvar.evaluate ( results ) ;
      for ( double x : results ) 
      { 
        // check the x-range  & fill 
        if ( xmin <= x && x <= xmax && obj.Fill ( x , weight ) ) { ++filled ; }
      }
    }
    return filled ;
  }
  // ==========================================================================
  template <class OBJECT>
  unsigned long 
  _param2_
  ( TTree*                            data        , 
    const Ostap::Utils::ProgressConf& progress    ,
    OBJECT&                           obj         ,
    const std::string&                xexpression ,
    const std::string&                yexpression ,
    const std::string&                selection   ,
    const unsigned long               first       ,
    const unsigned long               last        ,
    const double                      xmin        , 
    const double                      xmax        ,
    const double                      ymin        , 
    const double                      ymax        )
  {
    if ( 0 == data          ) { return 0 ; }
    //
    const unsigned long nEntries = std::min ( last , (unsigned long) data->GetEntries() ) ;
    if ( nEntries <= first  ) { return 0 ; }
    //
    Ostap::Formula xvar ( xexpression , data ) ;
    Ostap::Assert ( !(!xvar) && xvar.ok () , 
                    "Invalid xexpression"  , 
                    "Ostap::DataParam::parameterize" , 310 ) ;
    //
    Ostap::Formula yvar ( yexpression , data ) ;
    Ostap::Assert ( !(!yvar) && yvar.ok () , 
                    "Invalid yexpression"  ,
                    "Ostap::DataParam::parameterize" , 311 ) ;
    //
    std::unique_ptr<Ostap::Formula> cut { nullptr } ;
    if  ( !selection.empty() ) 
    { 
      cut = std::make_unique<Ostap::Formula>( selection , data ) ; 
      Ostap::Assert ( !(!cut) && cut->ok () , 
                      "Invalid selection"   ,
                      "Ostap::DataParam::parameterize" , 315 ) ;
    }
    //
    unsigned long filled = 0 ;
    // 
    Ostap::Utils::Notifier notify ( data , &xvar , &yvar , cut.get() ) ;
    std::vector<double> xresults {} ;
    std::vector<double> yresults {} ;
    /// make an explicit  loop 
    Ostap::Utils::ProgressBar bar (  nEntries - first , progress ) ;
    for ( unsigned long entry = first ; entry < nEntries ; ++entry , ++bar )
    {
      //
      long ievent = data->GetEntryNumber ( entry ) ;
      Ostap::Assert ( 0 <= ievent , 
                      "Errot in TTree::GetEventNumber" , 
                      "Ostap::DataParam::parameterize" , 316 ) ;
      //
      ievent      = data->LoadTree ( ievent ) ;      
      Ostap::Assert ( 0 <= ievent , 
                      "Error in TTree::LoadTree " , 
                      "Ostap::DataParam::parameterize" , 317 ) ;
      //
      const double weight = cut ? cut -> evaluate() : 1.0 ;
      if ( !weight    ) { continue ;}
      //
      xvar.evaluate ( xresults ) ;
      yvar.evaluate ( yresults ) ;
      for ( double x : xresults )
      { 
        if ( xmax < x || x< xmin ) { continue ; }
        for ( double y : yresults )
        { 
          if ( ymin <= y && y <= ymax && obj.fill ( x , y , weight ) ) { ++filled ; }
        }
      }
    }
    return filled ;
  }
  // ==========================================================================
  template <class OBJECT>
  Ostap::StatusCode 
  _param3_
  ( TTree*                            data        , 
    const Ostap::Utils::ProgressConf& progress    ,
    OBJECT&                           obj         ,
    const std::string&                xexpression ,
    const std::string&                yexpression ,
    const std::string&                zexpression ,
    const std::string&                selection   ,
    const unsigned long               first       ,
    const unsigned long               last        , 
    const double                      xmin        , 
    const double                      xmax        ,
    const double                      ymin        , 
    const double                      ymax        ,
    const double                      zmin        , 
    const double                      zmax        )
  {
    if ( 0 == data          ) { return 0 ; }
    //
    const unsigned long nEntries = std::min ( last , (unsigned long) data->GetEntries() ) ;
    if ( nEntries <= first  ) { return 0 ; }
    //
    Ostap::Formula xvar ( xexpression , data ) ;
    Ostap::Assert ( !(!xvar) && xvar.ok () , 
                    "Invalid xexpression"  , 
                    "Ostap::DataParam::parameterize" , 310 ) ;
    //
    Ostap::Formula yvar ( yexpression , data ) ;
    Ostap::Assert ( !(!yvar) && yvar.ok () , 
                    "Invalid yexpression"  , 
                    "Ostap::DataParam::parameterize" , 311 ) ;
    //
    Ostap::Formula zvar ( zexpression , data ) ;
    Ostap::Assert ( !(!zvar) && zvar.ok () , 
                    "Invalid zexpression"  , 
                    "Ostap::DataParam::parameterize" , 312 ) ;
    //
    std::unique_ptr<Ostap::Formula> cut { nullptr } ;
    if  ( !selection.empty() ) 
    { 
      cut = std::make_unique<Ostap::Formula>( selection , data ) ; 
      Ostap::Assert ( !(!cut) && cut->ok () , 
                      "Invalid selection"  ,
                      "Ostap::DataParam::parameterize" , 315 ) ;
    }
    //
    unsigned long filled = 0 ;
    //
    Ostap::Utils::Notifier notify ( data , &xvar , &yvar , &zvar , cut.get() ) ;
    std::vector<double> xresults {} ;
    std::vector<double> yresults {} ;
    std::vector<double> zresults {} ;
    /// make an explicit  loop 
    Ostap::Utils::ProgressBar bar (  nEntries - first , progress ) ;
    for ( unsigned long entry = first ; entry < nEntries ; ++entry , ++bar )
    {      //
      long ievent = data->GetEntryNumber ( entry ) ;
      Ostap::Assert ( 0 <= ievent , 
                      "Errot in TTree::GetEventNumber" , 
                      "Ostap::DataParam::parameterize" , 316 ) ;
      //
      ievent      = data->LoadTree ( ievent ) ;      
      Ostap::Assert ( 0 <= ievent , 
                      "Error in TTree::LoadTree " , 
                      "Ostap::DataParam::parameterize" , 317 ) ;
      //
      const double weight = cut ? cut -> evaluate() : 1.0 ;
      if ( !weight    ) { continue ;}
      //
      xvar.evaluate ( xresults ) ;
      yvar.evaluate ( yresults ) ;
      zvar.evaluate ( zresults ) ;
      for ( double x : xresults )
      { 
        if ( xmax < x || x < xmin ) { continue ; }
        for ( double y : yresults )
        { 
          if ( xmax < x || y < ymin ) { continue ; }
          for ( double z : zresults )
          { 
            if ( zmin <= z && z <= zmax && obj.fill ( x , y , z , weight ) ) { ++filled ; }
          } 
        } 
      }
    }
    return filled ;
  }
  // ==========================================================================
  template <class OBJECT>
  Ostap::StatusCode 
  _param4_
  ( TTree*                            data        , 
    const Ostap::Utils::ProgressConf& progress    ,
    OBJECT&                           obj         ,
    const std::string&                xexpression ,
    const std::string&                yexpression ,
    const std::string&                zexpression ,
    const std::string&                uexpression ,
    const std::string&                selection   ,
    const unsigned long               first       ,
    const unsigned long               last        ,
    const double                      xmin        , 
    const double                      xmax        ,
    const double                      ymin        , 
    const double                      ymax        ,
    const double                      zmin        , 
    const double                      zmax        ,
    const double                      umin        , 
    const double                      umax        )
  {
    if ( 0 == data          ) { return 0 ; }
    //
    const unsigned long nEntries = std::min ( last , (unsigned long) data->GetEntries() ) ;
    if ( nEntries <= first  ) { return 0 ; }
    //
    Ostap::Formula xvar ( xexpression , data ) ;
    Ostap::Assert ( !(!xvar) && xvar.ok () , 
                    "Invalid xexpression"  , 
                    "Ostap::DataParam::parameterize" , 310 ) ;
    //
    Ostap::Formula yvar ( yexpression , data ) ;
    Ostap::Assert ( !(!yvar) && yvar.ok () , 
                    "Invalid yexpression"  ,
                    "Ostap::DataParam::parameterize" , 311 ) ;
    //
    Ostap::Formula zvar ( zexpression , data ) ;
    Ostap::Assert ( !(!zvar) && zvar.ok () , 
                    "Invalid zexpression"  ,
                    "Ostap::DataParam::parameterize" , 312 ) ;
    //
    Ostap::Formula uvar ( uexpression , data ) ;
    Ostap::Assert ( !(!uvar) && uvar.ok () , 
                    "Invalid uexpression"  , 
                    "Ostap::DataParam::parameterize" , 313 ) ;
    //
    std::unique_ptr<Ostap::Formula> cut { nullptr } ;
    if  ( !selection.empty() ) 
    { 
      cut = std::make_unique<Ostap::Formula>( selection , data ) ; 
      Ostap::Assert ( !(!cut) && cut->ok () , 
                      "Invalid selection"   ,
                      "Ostap::DataParam::parameterize" , 315 ) ;
    }
    //
    unsigned long filled = 0 ;
    //
    Ostap::Utils::Notifier notify ( data , &xvar , &yvar , &zvar , &uvar , cut.get() ) ;
    std::vector<double> xresults {} ;
    std::vector<double> yresults {} ;
    std::vector<double> zresults {} ;
    std::vector<double> uresults {} ;
    /// make an explicit  loop 
    Ostap::Utils::ProgressBar bar (  nEntries - first , progress ) ;
    for ( unsigned long entry = first ; entry < nEntries ; ++entry , ++bar )
    {
      //
      long ievent = data->GetEntryNumber ( entry ) ;
      Ostap::Assert ( 0 <= ievent , 
                      "Errot in TTree::GetEventNumber" , 
                      "Ostap::DataParam::parameterize" , 316 ) ;
      //
      ievent      = data->LoadTree ( ievent ) ;      
      Ostap::Assert ( 0 <= ievent , 
                      "Error in TTree::LoadTree " , 
                      "Ostap::DataParam::parameterize" , 317 ) ;
      //
      const double weight = cut ? cut -> evaluate() : 1.0 ;
      if ( !weight    ) { continue ;}
      //
      xvar.evaluate ( xresults ) ;
      yvar.evaluate ( yresults ) ;
      zvar.evaluate ( zresults ) ;
      uvar.evaluate ( uresults ) ;
      for ( double x : xresults )
      { 
        if ( xmax < x || x < xmin ) { continue ; }
        for ( double y : yresults )
        { 
          if ( xmax < x || y < ymin ) { continue ; }
          for ( double z : zresults )
          { 
            if ( zmax < z || z < zmin ) { continue ; }
            for ( double u : uresults )
            { 
              if ( umin <= u && u <= umax && obj.fill ( x , y , z , u , weight ) ) { ++filled ; }
            } 
          }   
        } 
      }
    }
    return filled ;
  }
  // ==========================================================================
}    
// ============================================================================
// 1D stuff 
// ============================================================================
/** fill Legendre sum with data from the Tree 
 *  @see Ostap::Math::LegendreSum 
 *  @see Ostap::Math::LegendreSum::fill
 *  @param tree       (INPUT)  the input tree 
 *  @param sum        (UPDATE) the parameterization object 
 *  @param expression (INPUT)  expression to be parameterized
 *  @param first      (INPUT)  the first event in Tree 
 *  @param last       (INPUT)  the last  event in Tree 
 *  @code
 *  Tree*  tree = ...
 *  LegendreSum s ( 5 , -1 , 1 ) ;
 *  DataParam::parameterize ( tree , s , "x" ) ;
 *  @endcode
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-07-3
 */ 
// ============================================================================
unsigned long 
Ostap::DataParam::parameterize 
( TTree*                    tree       , 
  Ostap::Math::LegendreSum& sum        , 
  const std::string&        expression , 
  const unsigned long       first      ,
  const unsigned long       last       ) 
{ return parameterize ( tree , sum , expression , "" ,  first , last ) ; }
// ============================================================================
/*  fill Legendre sum with data from the Tree 
 *  @see Ostap::Math::LegendreSum 
 *  @see Ostap::Math::LegendreSum::fill
 *  @param tree       (INPUT)  the input tree 
 *  @param sum        (UPDATE) the parameterization object 
 *  @param expression (INPUT)  expression to be parameterized
 *  @param seelction  (INPUT)  selection/weight to be used 
 *  @param first      (INPUT)  the first event in Tree 
 *  @param last       (INPUT)  the last  event in Tree 
 *  @return  sum of weigths  used in parameterization
 *  @code
 *  Tree*  tree = ...
 *  LegendreSum s ( 5 , -1 , 1 ) ;
 *  DataParam::parameterize ( tree , s , "x"  , "y>10" ) ;
 *  @endcode
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-07-3
 */ 
// ============================================================================
unsigned long 
Ostap::DataParam::parameterize
( TTree*                    tree       , 
  Ostap::Math::LegendreSum& sum        , 
  const std::string&        expression , 
  const std::string&        selection  , 
  const unsigned long       first      ,
  const unsigned long       last       ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  /// reset the sum 
  sum *= 0 ;
  return _param1_ ( tree        , 
                    progress    , 
                    sum         , 
                    expression  , 
                    selection   , 
                    first       , 
                    last        , 
                    sum.xmin()  ,
                    sum.xmax()  ) ;
}
// ============================================================================
/*  fill Chebyshev sum with data from the Tree 
 *  @see Ostap::Math::ChebyshevSum 
 *  @see Ostap::Math::ChebyshevSum::fill
 *  @param tree       (INPUT)  the input tree 
 *  @param sum        (UPDATE) the parameterization object 
 *  @param expression (INPUT)  expression to be parameterized
 *  @param first      (INPUT)  the first event in Tree 
 *  @param last       (INPUT)  the last  event in Tree 
 *  @code
 *  Tree*  tree = ...
 *  ChebyshevSum s ( 5 , -1 , 1 ) ;
 *  DataParam::parameterize ( tree , s , "x" ) ;
 *  @endcode
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-07-3
 */ 
// ============================================================================
unsigned long 
Ostap::DataParam::parameterize
( TTree*                     tree       , 
  Ostap::Math::ChebyshevSum& sum        , 
  const std::string&        expression , 
  const unsigned long       first      ,
  const unsigned long       last       ) 
{ return parameterize ( tree , sum , expression , "" ,  first , last ) ; }
// ========================================================================
/* fill Chebyshev sum with data from the Tree 
 *  @see Ostap::Math::ChebyshevSum 
 *  @see Ostap::Math::ChebyshevSum::fill
 *  @param tree       (INPUT)  the input tree 
 *  @param sum        (UPDATE) the parameterization object 
 *  @param expression (INPUT)  expression to be parameterized
 *  @param seelction  (INPUT)  selection/weight to be used 
 *  @param first      (INPUT)  the first event in Tree 
 *  @param last       (INPUT)  the last  event in Tree 
 *  @code
 *  Tree*  tree = ...
 *  ChebyshevSum s ( 5 , -1 , 1 ) ;
 *  DataParam::parameterize ( tree , s , "x"  , "y>10" ) ;
 *  @endcode
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-07-3
 */
// ============================================================================
unsigned long 
Ostap::DataParam::parameterize
( TTree*                     tree       , 
  Ostap::Math::ChebyshevSum& sum        , 
  const std::string&         expression , 
  const std::string&         selection  , 
  const unsigned long        first      ,
  const unsigned long        last       ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  /// reset the sum 
  sum *= 0 ;
  return _param1_ ( tree        , 
                    progress    , 
                    sum         , 
                    expression  , 
                    selection   , 
                    first       , 
                    last        , 
                    sum.xmin()  ,
                    sum.xmax()  ) ;
}
// ============================================================================
/*  fill Bernstein sum with data from the Tree 
 *  @see Ostap::Math::Bernstein 
 *  @see Ostap::Math::Bernstein::fill
 *  @param tree       (INPUT)  the input tree 
 *  @param sum        (UPDATE) the parameterization object 
 *  @param expression (INPUT)  expression to be parameterized
 *  @param first      (INPUT)  the first event in Tree 
 *  @param last       (INPUT)  the last  event in Tree 
 *  @code
 *  Tree*  tree = ...
 *  Bernstein s ( 5 , -1 , 1 ) ;
 *  DataParam::parameterize ( tree , s , "x" ) ;
 *  @endcode
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-07-3
 */ 
// ============================================================================
unsigned long 
Ostap::DataParam::parameterize
( TTree*                     tree       , 
  Ostap::Math::Bernstein&    sum        , 
  const std::string&        expression , 
  const unsigned long       first      ,
  const unsigned long       last       ) 
{ return parameterize ( tree , sum , expression , "" ,  first , last ) ; }
// ========================================================================
/*  fill Bernstein sum with data from the Tree 
 *  @see Ostap::Math::Bernstein 
 *  @see Ostap::Math::Bernstein::fill
 *  @param tree       (INPUT)  the input tree 
 *  @param sum        (UPDATE) the parameterization object 
 *  @param expression (INPUT)  expression to be parameterized
 *  @param seelction  (INPUT)  selection/weight to be used 
 *  @param first      (INPUT)  the first event in Tree 
 *  @param last       (INPUT)  the last  event in Tree 
 *  @code
 *  Tree*  tree = ...
 *  Bernstein s ( 5 , -1 , 1 ) ;
 *  DataParam::parameterize ( tree , s , "x"  , "y>10" ) ;
 *  @endcode
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-07-3
 */ 
// ============================================================================
unsigned long 
Ostap::DataParam::parameterize
( TTree*                     tree       , 
  Ostap::Math::Bernstein&    sum        , 
  const std::string&         expression , 
  const std::string&         selection  , 
  const unsigned long        first      ,
  const unsigned long        last       ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  /// reset the sum 
  sum *= 0 ;
  return _param1_ ( tree        , 
                    progress    , 
                    sum         , 
                    expression  , 
                    selection   , 
                    first       , 
                    last        , 
                    sum.xmin()  ,
                    sum.xmax()  ) ;
}
// ============================================================================
// 2D 
// ============================================================================
/*  fill Legendre sum with data from the Tree 
 *  @see Ostap::Math::LegendreSum2 
 *  @see Ostap::Math::LegendreSum2::fill
 *  @param tree        (INPUT)  the input tree 
 *  @param sum         (UPDATE) the parameterization object 
 *  @param xexpression (INPUT)  x-expression to be parameterized
 *  @param yexpression (INPUT)  y-expression to be parameterized
 *  @param first       (INPUT)  the first event in Tree 
 *  @param last        (INPUT)  the last  event in Tree 
 *  @code
 *  Tree*  tree = ...
 *  LegendreSum2 s ( 5 , 3 , -1 , 1 , -2 , 2 ) ;
 *  DataParam::parameterize ( tree , s , "x" , "y/z") ;
 *  @endcode
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-07-3
 */ 
// ============================================================================
unsigned long 
Ostap::DataParam::parameterize
( TTree*                     tree        , 
  Ostap::Math::LegendreSum2& sum         ,
  const std::string&         xexpression , 
  const std::string&         yexpression , 
  const unsigned long        first       ,
  const unsigned long        last        ) 
{ return parameterize ( tree , sum , xexpression , yexpression , "" ,  first , last ) ; }
// ========================================================================
/*  fill Legendre sum with data from the Tree 
 *  @see Ostap::Math::LegendreSum2 
 *  @see Ostap::Math::LegendreSum2::fill
 *  @param tree        (INPUT)  the input tree 
 *  @param sum         (UPDATE) the parameterization object 
 *  @param xexpression (INPUT)  x-expression to be parameterized
 *  @param yexpression (INPUT)  y-expression to be parameterized
 *  @param selection   (INPUT)  selection/weight to be used 
 *  @param first       (INPUT)  the first event in Tree 
 *  @param last        (INPUT)  the last  event in Tree 
 *  @code
 *  Tree*  tree = ...
 *  LegendreSum2 s ( 5 , 2 ,  -1 , 1 , -4  , -5 ) ;
 *  DataParam::parameterize ( tree , s , "x"  , "z" , "y>10" ) ;
 *  @endcode
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-07-3
 */ 
// ============================================================================
unsigned long 
Ostap::DataParam::parameterize
( TTree*                     tree        , 
  Ostap::Math::LegendreSum2& sum         , 
  const std::string&         xexpression , 
  const std::string&         yexpression , 
  const std::string&         selection   , 
  const unsigned long        first       ,
  const unsigned long        last        ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  /// reset the sum 
  sum *= 0 ;
  return _param2_ ( tree        , 
                    progress    , 
                    sum         , 
                    xexpression , 
                    yexpression , 
                    selection   , 
                    first       , 
                    last        , 
                    sum.xmin()  ,
                    sum.xmax()  ,
                    sum.ymin()  ,
                    sum.ymax()  ) ;
}
// ============================================================================
/** fill Bernstein with data from the Tree 
 *  @see Ostap::Math::Bernstein2D
 *  @see Ostap::Math::Bernstein2D::fill
 *  @param tree        (INPUT)  the input tree 
 *  @param sum         (UPDATE) the parameterization object 
 *  @param xexpression (INPUT)  x-expression to be parameterized
 *  @param yexpression (INPUT)  y-expression to be parameterized
 *  @param first       (INPUT)  the first event in Tree 
 *  @param last        (INPUT)  the last  event in Tree 
 *  @code
 *  Tree*  tree = ...
 *  Bernsteinn2D s ( 5 , 3 , -1 , 1 , -2 , 2 ) ;
 *  DataParam::parameterize ( tree , s , "x" , "y/z") ;
 *  @endcode
 *  @attention it is less CPU efficient than Ostap::Math::LegendreSum2 
 *  @see Ostap::Math::LegendreSum2 
 *  @see Ostap::Math::LegendreSum2::fill
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2022-12-26
 */
// ============================================================================
unsigned long 
Ostap::DataParam::parameterize
( TTree*                     tree        , 
  Ostap::Math::Bernstein2D&  sum         , 
  const std::string&         xexpression , 
  const std::string&         yexpression , 
  const unsigned long        first       ,
  const unsigned long        last        ) 
{ return parameterize ( tree , sum , xexpression , yexpression , "" ,  first , last ) ; }
// ============================================================================
/*  fill Bernstein with data from the Tree 
 *  @see Ostap::Math::Bernstein2D
 *  @see Ostap::Math::Bernstein2D::fill
 *  @param tree        (INPUT)  the input tree 
 *  @param sum         (UPDATE) the parameterization object 
 *  @param xexpression (INPUT)  x-expression to be parameterized
 *  @param yexpression (INPUT)  y-expression to be parameterized
 *  @param selection   (INPUT)  selection/weight to be used 
 *  @param first       (INPUT)  the first event in Tree 
 *  @param last        (INPUT)  the last  event in Tree 
 *  @code
 *  Tree*  tree = ...
 *  Bernstein2D s ( 5 , 2 ,  -1 , 1 , -4  , -5 ) ;
 *  DataParam::parameterize ( tree , s , "x"  , "z" , "y>10" ) ;
 *  @endcode
 *  @attention it is less CPU efficient than Ostap::Math::LegendreSum2 
 *  @see Ostap::Math::LegendreSum2 
 *  @see Ostap::Math::LegendreSum2::fill
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2022-12-26
 */ 
// ============================================================================
unsigned long 
Ostap::DataParam::parameterize
( TTree*                     tree        , 
  Ostap::Math::Bernstein2D&  sum         , 
  const std::string&         xexpression , 
  const std::string&         yexpression , 
  const std::string&         selection   , 
  const unsigned long        first       ,
  const unsigned long        last        ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  /// reset the sum 
  sum *= 0 ;
  return _param2_ ( tree        , 
                    progress    , 
                    sum         , 
                    xexpression , 
                    yexpression , 
                    selection   , 
                    first       , 
                    last        , 
                    sum.xmin()  ,
                    sum.xmax()  ,
                    sum.ymin()  ,
                    sum.ymax()  ) ;
}
// ============================================================================
// 3D 
// ============================================================================`
/*  fill Legendre sum with data from the Tree 
 *  @see Ostap::Math::LegendreSum3 
 *  @see Ostap::Math::LegendreSum3::fill
 *  @param tree        (INPUT)  the input tree 
 *  @param sum         (UPDATE) the parameterization object 
 *  @param xexpression (INPUT)  x-expression to be parameterized
 *  @param yexpression (INPUT)  y-expression to be parameterized
 *  @param zexpression (INPUT)  z-expression to be parameterized
 *  @param first       (INPUT)  the first event in Tree 
 *  @param last        (INPUT)  the last  event in Tree 
 *  @code
 *  Tree*  tree = ...
 *  LegendreSum3 s ( 5 , 3 , 2 , -1 , 1 , -2 , 2 , 0 , 4 ) ;
 *  DataParam::parameterize ( tree , s , "x" , "y" , "y/z") ;
 *  @endcode
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-07-3
 */ 
// ============================================================================
unsigned long 
Ostap::DataParam::parameterize
( TTree*                     tree        , 
  Ostap::Math::LegendreSum3& sum         , 
  const std::string&         xexpression , 
  const std::string&         yexpression , 
  const std::string&         zexpression , 
  const unsigned long        first       ,
  const unsigned long        last        ) 
{ return parameterize ( tree , sum , xexpression , yexpression , zexpression , "" ,  first , last ) ; }
// ========================================================================
/*  fill Legendre sum with data from the Tree 
 *  @see Ostap::Math::LegendreSum3 
 *  @see Ostap::Math::LegendreSum3::fill
 *  @param tree        (INPUT)  the input tree 
 *  @param sum         (UPDATE) the parameterization object 
 *  @param xexpression (INPUT)  x-expression to be parameterized
 *  @param yexpression (INPUT)  y-expression to be parameterized
 *  @param zexpression (INPUT)  z-expression to be parameterized
 *  @param selection   (INPUT)  selection/weight to be used 
 *  @param first       (INPUT)  the first event in Tree 
 *  @param last        (INPUT)  the last  event in Tree 
 *  @code
 *  Tree*  tree = ...
 *  LegendreSum3 s ( 5 , 2 , 4 ,  -1 , 1 , -4  , -5 , 0 , 3 ) ;
 *  DataParam::parameterize ( tree , s , "x"  , "y" , "z" , "t>10" ) ;
 *  @endcode
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-07-3
 */ 
// ============================================================================
unsigned long 
Ostap::DataParam::parameterize
( TTree*                     tree        , 
  Ostap::Math::LegendreSum3& sum         , 
  const std::string&         xexpression , 
  const std::string&         yexpression , 
  const std::string&         zexpression , 
  const std::string&         selection   , 
  const unsigned long        first       ,
  const unsigned long        last        ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  /// reset the sum 
  sum *= 0 ;
  return _param3_ ( tree        , 
                    progress    , 
                    sum         , 
                    xexpression , 
                    yexpression , 
                    zexpression , 
                    selection   , 
                    first       , 
                    last        , 
                    sum.xmin()  ,
                    sum.xmax()  ,
                    sum.ymin()  ,
                    sum.ymax()  ,
                    sum.zmin()  ,
                    sum.zmax()  ) ;
}
// ============================================================================
/** fill Bernstein with data from the Tree 
 *  @see Ostap::Math::Bernstein3D
 *  @see Ostap::Math::Bernstein3D::fill
 *  @param tree        (INPUT)  the input tree 
 *  @param sum         (UPDATE) the parameterization object 
 *  @param xexpression (INPUT)  x-expression to be parameterized
 *  @param yexpression (INPUT)  y-expression to be parameterized
 *  @param zexpression (INPUT)  z-expression to be parameterized
 *  @param first       (INPUT)  the first event in Tree 
 *  @param last        (INPUT)  the last  event in Tree 
 *  @code
 *  Tree*  tree = ...
 *  Bernstein3D s ( 5 , 3 , 2 , -1 , 1 , -2 , 2 , 0 , 4 ) ;
 *  DataParam::parameterize ( tree , s , "x" , "y" , "y/z") ;
 *  @endcode
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-07-3
 */ 
// ============================================================================
unsigned long 
Ostap::DataParam::parameterize
( TTree*                     tree        , 
  Ostap::Math::Bernstein3D&  sum         , 
  const std::string&         xexpression , 
  const std::string&         yexpression , 
  const std::string&         zexpression , 
  const unsigned long        first       ,
  const unsigned long        last        ) 
{ return parameterize ( tree , sum , xexpression , yexpression , zexpression , "" ,  first , last ) ; }
// ========================================================================
/** fill Bernstein with data from the Tree 
 *  @see Ostap::Math::Bernstein3D
 *  @see Ostap::Math::Bernstein3D::fill
 *  @param tree        (INPUT)  the input tree 
 *  @param sum         (UPDATE) the parameterization object 
 *  @param xexpression (INPUT)  x-expression to be parameterized
 *  @param yexpression (INPUT)  y-expression to be parameterized
 *  @param zexpression (INPUT)  z-expression to be parameterized
 *  @param selection   (INPUT)  selection/weight to be used 
 *  @param first       (INPUT)  the first event in Tree 
 *  @param last        (INPUT)  the last  event in Tree 
 *  @code
 *  Tree*  tree = ...
 *  Bernstein3D s ( 5 , 3 , 2 , -1 , 1 , -2 , 2 , 0 , 4 ) ;
 *  DataParam::parameterize ( tree , s , "x"  , "y" , "z" , "t>10" ) ;
 *  @endcode
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-07-3
 */ 
// ============================================================================
unsigned long 
Ostap::DataParam::parameterize
( TTree*                     tree        , 
  Ostap::Math::Bernstein3D&  sum         , 
  const std::string&         xexpression , 
  const std::string&         yexpression , 
  const std::string&         zexpression , 
  const std::string&         selection   , 
  const unsigned long        first       ,
  const unsigned long        last        ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  /// reset the sum 
  sum *= 0 ;
  return _param3_ ( tree        , 
                    progress    , 
                    sum         , 
                    xexpression , 
                    yexpression , 
                    zexpression , 
                    selection   , 
                    first       , 
                    last        , 
                    sum.xmin()  ,
                    sum.xmax()  ,
                    sum.ymin()  ,
                    sum.ymax()  ,
                    sum.zmin()  ,
                    sum.zmax()  ) ;
}
// ============================================================================
// 4D 
// ============================================================================
/*  fill Legendre sum with data from the Tree 
 *  @see Ostap::Math::LegendreSum4 
 *  @see Ostap::Math::LegendreSum4::fill
 *  @param tree        (INPUT)  the input tree 
 *  @param sum         (UPDATE) the parameterization object 
 *  @param xexpression (INPUT)  x-expression to be parameterized
 *  @param yexpression (INPUT)  y-expression to be parameterized
 *  @param zexpression (INPUT)  z-expression to be parameterized
 *  @param uexpression (INPUT)  u-expression to be parameterized
 *  @param first       (INPUT)  the first event in Tree 
 *  @param last        (INPUT)  the last  event in Tree 
 *  @code
 *  Tree*  tree = ...
 *  LegendreSum4 s ( 5 , 3 , 2 , 2 , -1 , 1 , -2 , 2 , 0 , 4 , 0 , 1 ) ;
 *  DataParam::parameterize ( tree , s , "x" , "y" , "z" , "t" ) ;
 *  @endcode
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-07-3
 */ 
// ============================================================================
unsigned long 
Ostap::DataParam::parameterize
( TTree*                     tree        , 
  Ostap::Math::LegendreSum4& sum         ,
  const std::string&         xexpression , 
  const std::string&         yexpression , 
  const std::string&         zexpression , 
  const std::string&         uexpression , 
  const unsigned long        first       ,
  const unsigned long        last        ) 
{ return parameterize ( tree , sum , xexpression , yexpression , zexpression , "" ,  first , last ) ; }
// ========================================================================
/*  fill Legendre sum with data from the Tree 
 *  @see Ostap::Math::LegendreSum4 
 *  @see Ostap::Math::LegendreSum4::fill
 *  @param tree        (INPUT)  the input tree 
 *  @param sum         (UPDATE) the parameterization object 
 *  @param xexpression (INPUT)  x-expression to be parameterized
 *  @param yexpression (INPUT)  y-expression to be parameterized
 *  @param zexpression (INPUT)  z-expression to be parameterized
 *  @param uexpression (INPUT)  u-expression to be parameterized
 *  @param selection   (INPUT)  selection/weight to be used 
 *  @param first       (INPUT)  the first event in Tree 
 *  @param last        (INPUT)  the last  event in Tree 
 *  @code
 *  Tree*  tree = ...
 *  LegendreSum4 s ( 5 , 2 , 4 ,  -1 , 1 , -4  , -5 , 0 , 3 ) ;
 *  DataParam::parameterize ( tree , s , "x"  , "y" , "z" , "t" , "q>10" ) ;
 *  @endcode
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-07-3
 */ 
// ============================================================================
unsigned long 
Ostap::DataParam::parameterize
( TTree*                     tree        , 
  Ostap::Math::LegendreSum4& sum         , 
  const std::string&         xexpression , 
  const std::string&         yexpression , 
  const std::string&         zexpression , 
  const std::string&         uexpression , 
  const std::string&         selection   , 
  const unsigned long        first       ,
  const unsigned long        last        ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  /// reset the sum 
  sum *= 0 ;
  return _param3_ ( tree        , 
                    progress    , 
                    sum         , 
                    xexpression , 
                    yexpression , 
                    zexpression , 
                    selection   , 
                    first       , 
                    last        , 
                    sum.xmin()  ,
                    sum.xmax()  ,
                    sum.ymin()  ,
                    sum.ymax()  ,
                    sum.zmin()  ,
                    sum.zmax()  ) ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================
