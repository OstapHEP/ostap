// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <cstring>
// ============================================================================
// ROOT
// ============================================================================
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
// ============================================================================
// Ostap
// ===========================================================================
#include "Ostap/StatVar.h"
#include "Ostap/Project.h"
#include "Ostap/Exception.h"
#include "Ostap/ECDF.h"
#include "Ostap/Polynomials.h"
#include "Ostap/Bernstein.h"
#include "Ostap/Bernstein2D.h"
#include "Ostap/Bernstein3D.h"
#include "Ostap/Parameterization.h"
// ============================================================================
// Local 
// ============================================================================
#include "local_utils.h"
#include "status_codes.h"
#include "hstats.h"
// ============================================================================
/** @file
 *  Implementation file for class Ostap::HistoProject
 *  @see Ostap::HistoProject
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2015-10-08
 */
// ============================================================================
// construtctor with the progress flag
// ============================================================================
Ostap::Project::Project					
( const Ostap::Utils::ProgressConf& progress )
  : Ostap::StatVar ( progress )
{}
// ============================================================================

// ============================================================================
// !D-histos 
// ============================================================================
/*  Project data in the 1D-ihstogram
 *  @param data       (INPUT)  input data 
 *  @param histo      (UPDATE) the !D-histogram 
 *  @param expression (INPUT)  expression 
 *  @param selection  (INPUT)  selection 
 *  @param first      (INOPUT) the first event to process (inclusive) 
 *  @param last       (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 *  @see   TH1::Fill
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project1
( TTree*                  data        ,
  TH1*                    histo       ,
  const std::string&      expression  ,
  const std::string&      selection   ,
  const Ostap::EventIndex first       ,
  const Ostap::EventIndex last        ) const 
{
  if ( !histo || 1 != histo->GetDimension () ) { return INVALID_TH1 ; }
  //
  Ostap::Utils::H1 h1 ( histo )    ;
  //
  const TAxis* xaxis = histo->GetXaxis () ;
  if ( !xaxis ) { return INVALID_XAXIS ; }
  //
  return get_stat ( data                ,
		    h1                  ,
		    expression          ,
		    selection           ,
		    first               ,
		    last                ,
		    xaxis -> GetXmin () ,
		    xaxis -> GetXmax () ) ;
}
// ============================================================================
/*  Project data in the 1D-ihstogram
 *  @param data       (INPUT)  input data 
 *  @param histo      (UPDATE) the !D-histogram 
 *  @param expression (INPUT)  expression 
 *  @param selection  (INPUT)  selection 
 *  @param first      (INOPUT) the first event to process (inclusive) 
 *  @param last       (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 *  @see   TH1::Fill
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project1
( const RooAbsData*       data        ,
  TH1*                    histo       ,
  const std::string&      expression  ,
  const std::string&      selection   ,
  const std::string&      cut_range   ,  
  const Ostap::EventIndex first       ,
  const Ostap::EventIndex last        ) const 
{
  if ( !histo || 1 != histo->GetDimension () ) { return INVALID_TH1 ; }
  //
  Ostap::Utils::H1 h1 ( histo )    ;
  //
  const TAxis* xaxis = histo->GetXaxis () ;
  if ( !xaxis ) { return INVALID_XAXIS ; }
  //
  return get_stat ( data                ,
		    h1                  ,
		    expression          ,
		    selection           ,
		    cut_range           , 
		    first               ,
		    last                ,
		    xaxis -> GetXmin () ,
		    xaxis -> GetXmax () ) ;
}
// ============================================================================

// ============================================================================
// 2D-histos 
// ============================================================================
/*  Project data in the 2D-ihstoram
 *  @param data        (INPUT)  input data 
 *  @param histo       (UPDATE) the !D-histogram 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param selection   (INPUT)  selection/weight  
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 *  @see   TH1::Fill
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project2
( TTree*                  data        ,
  TH2*                    histo       ,
  const std::string&      expression1 ,
  const std::string&      expression2 ,
  const std::string&      selection   ,
  const Ostap::EventIndex first       ,
  const Ostap::EventIndex last        ) const
{
  if ( !histo || 2 != histo->GetDimension () ) { return INVALID_TH2 ; }
  //
  Ostap::Utils::H2 h2 ( histo )    ;
  //
  const TAxis* xaxis = histo->GetXaxis () ;
  if ( !xaxis ) { return INVALID_XAXIS ; }
  //
  const TAxis* yaxis = histo->GetYaxis () ;
  if ( !yaxis ) { return INVALID_YAXIS ; }
  //
  return get_stat ( data                ,
		    h2                  ,
		    expression1         ,
		    expression2         ,
		    selection           ,
		    first               ,
		    last                ,
		    xaxis -> GetXmin () ,
		    xaxis -> GetXmax () ,
		    yaxis -> GetXmin () ,
		    yaxis -> GetXmax () ) ;
}
// ============================================================================
/*  Project data in the 2D-ihstoram
 *  @param data        (INPUT)  input data 
 *  @param histo       (UPDATE) the !D-histogram 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param selection   (INPUT)  selection/weight  
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 *  @see   TH1::Fill
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project2
( const RooAbsData*       data        ,
  TH2*                    histo       ,
  const std::string&      expression1 ,
  const std::string&      expression2 ,
  const std::string&      selection   ,
  const std::string&      cut_range   ,
  const Ostap::EventIndex first       ,
  const Ostap::EventIndex last        ) const
{
  if ( !histo || 2 != histo->GetDimension () ) { return INVALID_TH2 ; }
  //
  Ostap::Utils::H2 h2 ( histo )    ;
  //
  const TAxis* xaxis = histo->GetXaxis () ;
  if ( !xaxis ) { return INVALID_XAXIS ; }
  //
  const TAxis* yaxis = histo->GetYaxis () ;
  if ( !yaxis ) { return INVALID_YAXIS ; }
  //
  return get_stat ( data                ,
		    h2                  ,
		    expression1         ,
		    expression2         ,
		    selection           ,
		    cut_range           , 
		    first               ,
		    last                ,
		    xaxis -> GetXmin () ,
		    xaxis -> GetXmax () ,
		    yaxis -> GetXmin () ,
		    yaxis -> GetXmax () ) ;
}
// ============================================================================


// ============================================================================
// 3D-histos 
// ============================================================================
/* Project data in the 3D-ihstoram
 *  @param data        (INPUT)  input data 
 *  @param histo       (UPDATE) the !D-histogram 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param expression2 (INPUT)  expression3 
 *  @param selection   (INPUT)  selection/weight  
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 *  @see   TH1::Fill
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project3
( TTree*                  data        ,
  TH3*                    histo       ,
  const std::string&      expression1 ,
  const std::string&      expression2 ,
  const std::string&      expression3 ,
  const std::string&      selection   ,
  const Ostap::EventIndex first       ,
  const Ostap::EventIndex last        ) const
{
  if ( !histo || 3 != histo->GetDimension () ) { return INVALID_TH3 ; }
  //
  Ostap::Utils::H3 h3 ( histo )    ;
  //
  const TAxis* xaxis = histo->GetXaxis () ;
  if ( !xaxis ) { return INVALID_XAXIS ; }
  //
  const TAxis* yaxis = histo->GetYaxis () ;
  if ( !yaxis ) { return INVALID_YAXIS ; }
  //
  const TAxis* zaxis = histo->GetZaxis () ;
  if ( !zaxis ) { return INVALID_ZAXIS ; }
  //
  return get_stat ( data                ,
		    h3                  ,
		    expression1         ,
		    expression2         ,
		    expression3         ,
		    selection           ,
		    first               ,
		    last                ,
		    xaxis -> GetXmin () ,
		    xaxis -> GetXmax () ,
		    yaxis -> GetXmin () ,
		    yaxis -> GetXmax () ,
		    zaxis -> GetXmin () ,
		    zaxis -> GetXmax () ) ;
} 
// ============================================================================
/* Project data in the 2D-ihstoram
 *  @param data        (INPUT)  input data 
 *  @param histo       (UPDATE) the !D-histogram 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param expression2 (INPUT)  expression3 
 *  @param selection   (INPUT)  selection/weight  
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 *  @see   TH1::Fill
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project3
( const RooAbsData*       data        ,
  TH3*                    histo       ,
  const std::string&      expression1 ,
  const std::string&      expression2 ,
  const std::string&      expression3 ,
  const std::string&      selection   ,
  const std::string&      cut_range   , 
  const Ostap::EventIndex first       ,
  const Ostap::EventIndex last        ) const
{
  if ( !histo || 3 != histo->GetDimension () ) { return INVALID_TH3 ; }
  //
  Ostap::Utils::H3 h3 ( histo )    ;
  //
  const TAxis* xaxis = histo->GetXaxis () ;
  if ( !xaxis ) { return INVALID_XAXIS ; }
  //
  const TAxis* yaxis = histo->GetYaxis () ;
  if ( !yaxis ) { return INVALID_YAXIS ; }
  //
  const TAxis* zaxis = histo->GetZaxis () ;
  if ( !zaxis ) { return INVALID_ZAXIS ; }
  //
  return get_stat ( data                ,
		    h3                  ,
		    expression1         ,
		    expression2         ,
		    expression3         ,
		    selection           ,
		    cut_range           , 
		    first               ,
		    last                ,
		    xaxis -> GetXmin () ,
		    xaxis -> GetXmax () ,
		    yaxis -> GetXmin () ,
		    yaxis -> GetXmax () ,
		    zaxis -> GetXmin () ,
		    zaxis -> GetXmax () ) ;
}
// ============================================================================


// ============================================================================
// ECDF & WECDF 
// ============================================================================
/*  Project data in ECDS/WECDF 
 *  @param data       (INPUT)  input data 
 *  @param poly       (UPDATE) the polynoimal 
 *  @param expression (INPUT)  expression 
 *  @param selection  (INPUT)  selection (as boolean)
 *  @param first      (INOPUT) the first event to process (inclusive) 
 *  @param last       (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project1
( TTree*                    data       ,
  Ostap::Math::ECDF&        ecdf       ,
  const std::string&        expression ,
  const std::string&        selection  ,
  const Ostap::EventIndex   first      ,
  const Ostap::EventIndex   last       , 
  const Ostap::DataType     xmin       , 
  const Ostap::DataType     xmax       ) const
{
  return get_stat ( data                ,
		    ecdf                ,
		    expression          ,
		    selection           ,
		    first               ,
		    last                ,
		    xmin                ,
		    xmax                ) ;
}
// ============================================================================
/*  Project data in ECDS/WECDF 
 *  @param data       (INPUT)  input data 
 *  @param poly       (UPDATE) the polynoimal 
 *  @param expression (INPUT)  expression 
 *  @param selection  (INPUT)  selection (as weighjt )
 *  @param first      (INOPUT) the first event to process (inclusive) 
 *  @param last       (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project1
( TTree*                    data       ,
  Ostap::Math::WECDF&       ecdf       ,
  const std::string&        expression ,
  const std::string&        selection  ,
  const Ostap::EventIndex   first      ,
  const Ostap::EventIndex   last       , 
  const Ostap::DataType     xmin       , 
  const Ostap::DataType     xmax       ) const
{
  return get_stat ( data                ,
		    ecdf                ,
		    expression          ,
		    selection           ,
		    first               ,
		    last                ,
		    xmin                ,
		    xmax                ) ;
}
// ============================================================================
/*  Project data in ECDS/WECDF 
 *  @param data       (INPUT)  input data 
 *  @param poly       (UPDATE) the polynoimal 
 *  @param expression (INPUT)  expression 
 *  @param selection  (INPUT)  selection (as weight) 
 *  @param first      (INOPUT) the first event to process (inclusive) 
 *  @param last       (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project1
( const RooAbsData*         data       ,
  Ostap::Math::WECDF&       ecdf       ,
  const std::string&        expression ,
  const std::string&        selection  ,
  const std::string&        cut_range  , 
  const Ostap::EventIndex   first      , 
  const Ostap::EventIndex   last       , 
  const Ostap::DataType     xmin       ,
  const Ostap::DataType     xmax       ) const 
{
  return get_stat ( data                ,
		    ecdf                ,
		    expression          ,
		    selection           ,
		    cut_range           , 
		    first               ,
		    last                ,
		    xmin                ,
		    xmax                ) ;
}
// ============================================================================


// ============================================================================
// 1D-chebyshev polynoials: on-flight parameterization
// ============================================================================

// ============================================================================
/*  Project data in 1D-polyom=ominal: on-flight parameterisation 
 *  @param data       (INPUT)  input data 
 *  @param poly       (UPDATE) the polynoimal 
 *  @param expression (INPUT)  expression 
 *  @param selection  (INPUT)  selection (as boolean)
 *  @param first      (INOPUT) the first event to process (inclusive) 
 *  @param last       (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project1
( TTree*                     data       ,
  Ostap::Math::ChebyshevSum& poly       ,
  const std::string&         expression ,
  const std::string&         selection  ,
  const Ostap::EventIndex    first      ,
  const Ostap::EventIndex    last       ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression          ,
		    selection           ,
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        ) ;
}
// =========================================================================
/*  Project data in 1D-polyom=ominal: on-flight parameterisation 
 *  @param data       (INPUT)  input data 
 *  @param poly       (UPDATE) the polynoimal 
 *  @param expression (INPUT)  expression 
 *  @param selection  (INPUT)  selection (as boolean)
 *  @param first      (INOPUT) the first event to process (inclusive) 
 *  @param last       (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project1
( const RooAbsData*          data       ,
  Ostap::Math::ChebyshevSum& poly       ,
  const std::string&         expression ,
  const std::string&         selection  ,
  const std::string&         cut_range  , 
  const Ostap::EventIndex    first      ,
  const Ostap::EventIndex    last       ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression          ,
		    selection           ,
		    cut_range           , 
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        ) ;
}
// =========================================================================
  

// ============================================================================
// 1D-Legendre polynoials: on-flight parameterization
// ============================================================================

// ============================================================================
/*  Project data in 1D-polyom=ominal: on-flight parameterisation 
 *  @param data       (INPUT)  input data 
 *  @param poly       (UPDATE) the polynoimal 
 *  @param expression (INPUT)  expression 
 *  @param selection  (INPUT)  selection (as boolean)
 *  @param first      (INOPUT) the first event to process (inclusive) 
 *  @param last       (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project1
( TTree*                     data       ,
  Ostap::Math::LegendreSum&  poly       ,
  const std::string&         expression ,
  const std::string&         selection  ,
  const Ostap::EventIndex    first      ,
  const Ostap::EventIndex    last       ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression          ,
		    selection           ,
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        ) ;
}
// =========================================================================
/*  Project data in 1D-polyom=ominal: on-flight parameterisation 
 *  @param data       (INPUT)  input data 
 *  @param poly       (UPDATE) the polynoimal 
 *  @param expression (INPUT)  expression 
 *  @param selection  (INPUT)  selection (as boolean)
 *  @param first      (INOPUT) the first event to process (inclusive) 
 *  @param last       (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project1
( const RooAbsData*          data       ,
  Ostap::Math::LegendreSum&  poly       ,
  const std::string&         expression ,
  const std::string&         selection  ,
  const std::string&         cut_range  , 
  const Ostap::EventIndex    first      ,
  const Ostap::EventIndex    last       ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression          ,
		    selection           ,
		    cut_range           , 
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        ) ;
}
// ============================================================================

// ============================================================================
// 1D-Legendre polynoials: on-flight parameterization
// ============================================================================

// ============================================================================
/*  Project data in 1D-polyom=ominal: on-flight parameterisation 
 *  @param data       (INPUT)  input data 
 *  @param poly       (UPDATE) the polynoimal 
 *  @param expression (INPUT)  expression 
 *  @param selection  (INPUT)  selection (as boolean)
 *  @param first      (INOPUT) the first event to process (inclusive) 
 *  @param last       (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project1
( TTree*                     data       ,
  Ostap::Math::Bernstein&    poly       ,
  const std::string&         expression ,
  const std::string&         selection  ,
  const Ostap::EventIndex    first      ,
  const Ostap::EventIndex    last       ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression          ,
		    selection           ,
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        ) ;
}
// =========================================================================
/*  Project data in 1D-polyom=ominal: on-flight parameterisation 
 *  @param data       (INPUT)  input data 
 *  @param poly       (UPDATE) the polynoimal 
 *  @param expression (INPUT)  expression 
 *  @param selection  (INPUT)  selection (as boolean)
 *  @param first      (INOPUT) the first event to process (inclusive) 
 *  @param last       (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project1
( const RooAbsData*          data       ,
  Ostap::Math::Bernstein&    poly       ,
  const std::string&         expression ,
  const std::string&         selection  ,
  const std::string&         cut_range  , 
  const Ostap::EventIndex    first      ,
  const Ostap::EventIndex    last       ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression          ,
		    selection           ,
		    cut_range           , 
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        ) ;
}
// =========================================================================


// =========================================================================
// 2D-Bernstein polynomial: on-flight parameterisation
// =========================================================================
/*  Project data in 2D-polyoominal: on-flight parameterisation 
 *  @param data        (INPUT)  input data 
 *  @param poly        (UPDATE) the polynoimal 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param selection   (INPUT)  selection (as boolean)
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// =========================================================================
Ostap::StatusCode
Ostap::Project::project2
( TTree*                     data        ,
  Ostap::Math::Bernstein2D&  poly        ,
  const std::string&         expression1 ,
  const std::string&         expression2 ,
  const std::string&         selection   ,
  const Ostap::EventIndex    first       ,
  const Ostap::EventIndex    last        ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression1         ,
		    expression2         ,
		    selection           ,
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        , 
		    poly.ymin ()        ,
		    poly.ymax ()        ) ;
}
// =========================================================================
/*  Project data in 2D-polyoominal: on-flight parameterisation 
 *  @param data        (INPUT)  input data 
 *  @param poly        (UPDATE) the polynoimal 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param selection   (INPUT)  selection (as boolean)
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// =========================================================================
Ostap::StatusCode
Ostap::Project::project2
( const RooAbsData*          data        ,
  Ostap::Math::Bernstein2D&  poly        ,
  const std::string&         expression1 ,
  const std::string&         expression2 ,
  const std::string&         selection   ,
  const std::string&         cut_range   , 
  const Ostap::EventIndex    first       ,
  const Ostap::EventIndex    last        ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression1         ,
		    expression2         ,
		    selection           ,
		    cut_range           , 
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        , 
		    poly.ymin ()        ,
		    poly.ymax ()        ) ;
}
// =========================================================================


// =========================================================================
// 2D-Legendre polynomial: on-flight parameterisation
// =========================================================================
/*  Project data in 2D-polyoominal: on-flight parameterisation 
 *  @param data        (INPUT)  input data 
 *  @param poly        (UPDATE) the polynoimal 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param selection   (INPUT)  selection (as boolean)
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// =========================================================================
Ostap::StatusCode
Ostap::Project::project2
( TTree*                     data        ,
  Ostap::Math::LegendreSum2& poly        ,
  const std::string&         expression1 ,
  const std::string&         expression2 ,
  const std::string&         selection   ,
  const Ostap::EventIndex    first       ,
  const Ostap::EventIndex    last        ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression1         ,
		    expression2         ,
		    selection           ,
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        , 
		    poly.ymin ()        ,
		    poly.ymax ()        ) ;
}
// =========================================================================
/*  Project data in 2D-polyoominal: on-flight parameterisation 
 *  @param data        (INPUT)  input data 
 *  @param poly        (UPDATE) the polynoimal 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param selection   (INPUT)  selection (as boolean)
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// =========================================================================
Ostap::StatusCode
Ostap::Project::project2
( const RooAbsData*          data        ,
  Ostap::Math::LegendreSum2& poly        ,
  const std::string&         expression1 ,
  const std::string&         expression2 ,
  const std::string&         selection   ,
  const std::string&         cut_range   , 
  const Ostap::EventIndex    first       ,
  const Ostap::EventIndex    last        ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression1         ,
		    expression2         ,
		    selection           ,
		    cut_range           , 
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        , 
		    poly.ymin ()        ,
		    poly.ymax ()        ) ;
}
// ============================================================================


// ============================================================================
// 3D-Bernstein polynomials 
// ============================================================================
/** Project data in 2D-polyoominal: on-flight parameterisation 
 *  @param data        (INPUT)  input data 
 *  @param poly        (UPDATE) the polynoimal 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param expression3 (INPUT)  expression2 
 *  @param selection   (INPUT)  selection (as boolean)
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project3
( TTree*                     data        ,
  Ostap::Math::Bernstein3D&  poly        ,
  const std::string&         expression1 ,
  const std::string&         expression2 ,
  const std::string&         expression3 ,
  const std::string&         selection   ,
  const Ostap::EventIndex    first       ,
  const Ostap::EventIndex    last        ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression1         ,
		    expression2         ,
		    expression3         ,
		    selection           ,
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        , 
		    poly.ymin ()        ,
		    poly.ymax ()        , 
		    poly.zmin ()        ,
		    poly.zmax ()        ) ;
}
// =========================================================================
/*  Project data in 2D-polyoominal: on-flight parameterisation 
 *  @param data        (INPUT)  input data 
 *  @param poly        (UPDATE) the polynoimal 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param selection   (INPUT)  selection (as boolean)
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// =========================================================================
Ostap::StatusCode
Ostap::Project::project3
( const RooAbsData*          data        ,
  Ostap::Math::Bernstein3D&  poly        ,
  const std::string&         expression1 ,
  const std::string&         expression2 ,
  const std::string&         expression3 ,
  const std::string&         selection   ,
  const std::string&         cut_range   , 
  const Ostap::EventIndex    first       ,
  const Ostap::EventIndex    last        ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression1         ,
		    expression2         ,
		    expression3         ,
		    selection           ,
		    cut_range           , 
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        , 
		    poly.ymin ()        ,
		    poly.ymax ()        ,
		    poly.zmin ()        ,
		    poly.zmax ()        ) ;
}
// ============================================================================

// ============================================================================
// 3D-Legendre polynomials 
// ============================================================================
/** Project data in 2D-polyoominal: on-flight parameterisation 
 *  @param data        (INPUT)  input data 
 *  @param poly        (UPDATE) the polynoimal 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param expression3 (INPUT)  expression2 
 *  @param selection   (INPUT)  selection (as boolean)
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project3
( TTree*                     data        ,
  Ostap::Math::LegendreSum3& poly        ,
  const std::string&         expression1 ,
  const std::string&         expression2 ,
  const std::string&         expression3 ,
  const std::string&         selection   ,
  const Ostap::EventIndex    first       ,
  const Ostap::EventIndex    last        ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression1         ,
		    expression2         ,
		    expression3         ,
		    selection           ,
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        , 
		    poly.ymin ()        ,
		    poly.ymax ()        , 
		    poly.zmin ()        ,
		    poly.zmax ()        ) ;
}
// =========================================================================
/*  Project data in 2D-polyoominal: on-flight parameterisation 
 *  @param data        (INPUT)  input data 
 *  @param poly        (UPDATE) the polynoimal 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param selection   (INPUT)  selection (as boolean)
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// =========================================================================
Ostap::StatusCode
Ostap::Project::project3
( const RooAbsData*          data        ,
  Ostap::Math::LegendreSum3& poly        ,
  const std::string&         expression1 ,
  const std::string&         expression2 ,
  const std::string&         expression3 ,
  const std::string&         selection   ,
  const std::string&         cut_range   , 
  const Ostap::EventIndex    first       ,
  const Ostap::EventIndex    last        ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression1         ,
		    expression2         ,
		    expression3         ,
		    selection           ,
		    cut_range           , 
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        , 
		    poly.ymin ()        ,
		    poly.ymax ()        ,
		    poly.zmin ()        ,
		    poly.zmax ()        ) ;
}
// ============================================================================


// ============================================================================
// 4D-Legendre polynomials 
// ============================================================================
/** Project data in 2D-polyoominal: on-flight parameterisation 
 *  @param data        (INPUT)  input data 
 *  @param poly        (UPDATE) the polynoimal 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param expression3 (INPUT)  expression2 
 *  @param selection   (INPUT)  selection (as boolean)
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project4
( TTree*                     data        ,
  Ostap::Math::LegendreSum4& poly        ,
  const std::string&         expression1 ,
  const std::string&         expression2 ,
  const std::string&         expression3 ,
  const std::string&         expression4 ,
  const std::string&         selection   ,
  const Ostap::EventIndex    first       ,
  const Ostap::EventIndex    last        ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression1         ,
		    expression2         ,
		    expression3         ,
		    expression4         ,
		    selection           ,
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        , 
		    poly.ymin ()        ,
		    poly.ymax ()        , 
		    poly.zmin ()        ,
		    poly.zmax ()        , 
		    poly.umin ()        ,
		    poly.umax ()        ) ;
}
// =========================================================================
/*  Project data in 2D-polyoominal: on-flight parameterisation 
 *  @param data        (INPUT)  input data 
 *  @param poly        (UPDATE) the polynoimal 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param selection   (INPUT)  selection (as boolean)
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// =========================================================================
Ostap::StatusCode
Ostap::Project::project4
( const RooAbsData*          data        ,
  Ostap::Math::LegendreSum4& poly        ,
  const std::string&         expression1 ,
  const std::string&         expression2 ,
  const std::string&         expression3 ,
  const std::string&         expression4 ,
  const std::string&         selection   ,
  const std::string&         cut_range   , 
  const Ostap::EventIndex    first       ,
  const Ostap::EventIndex    last        ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression1         ,
		    expression2         ,
		    expression3         ,
		    expression4         ,
		    selection           ,
		    cut_range           , 
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        , 
		    poly.ymin ()        ,
		    poly.ymax ()        ,
		    poly.zmin ()        ,
		    poly.zmax ()        , 
		    poly.umin ()        ,
		    poly.umax ()        ) ;
}
// ============================================================================




// ============================================================================
//                                                                      The END 
// ============================================================================
