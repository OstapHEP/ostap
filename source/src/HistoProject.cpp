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
// Local: 
// ============================================================================
#include "Ostap/Statistic.h"
#include "Ostap/HistoProject.h"
#include "Ostap/Exception.h"
// ============================================================================
#include "Exception.h"
#include "local_utils.h"
#include "status_codes.h"
// ============================================================================
/** @file
 *  Implementation file for class Ostap::HistoProject
 *  @see Ostap::HistoProject
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2015-10-08
 */
// ============================================================================
// constructor 
// ============================================================================
Ostap::Utils::TH1_Statistic::TH1_Statistic
( TH1* histo )
  : m_histo ( histo )
{
  Ostap::Assert ( histo  && 1 == histo->GetDimension () ,
		  "Inbvald TH1 pointer"                 ,
		  "Ostap:L=:Utils::TH1_Statistic"       , 
		  INVALID_TH1  , __FILE__ , __LINE__    ) ;
}
// ============================================================================
// reset the histogram
// ============================================================================
void Ostap::Utils::TH1_Statistic::reset  ()
{ 
  ESentry senty {} ;
  if ( m_histo ) { m_histo -> Reset () ; m_histo -> Sumw2 () ; }
}
// ============================================================================
// update the content 
// ============================================================================
void Ostap::Utils::TH1_Statistic::update
( const double x )
{ if ( m_histo && std::isfinite ( x ) ) { m_histo->Fill ( x ) ; } }
// ============================================================================
// constructor 
// ============================================================================
Ostap::Utils::TH1_WStatistic::TH1_WStatistic
( TH1* histo )
  : m_histo ( histo )
{
  Ostap::Assert ( histo  && 1 == histo->GetDimension () ,
		  "Inbvald TH1 pointer"                 ,
		  "Ostap:L=:Utils::TH1_WStatistic"      , 
		  INVALID_TH1  , __FILE__ , __LINE__    ) ;
}
// ============================================================================
// reset the histogram
// ============================================================================
void Ostap::Utils::TH1_WStatistic::reset  ()
{
  ESentry senty {} ;
  if ( m_histo ) { m_histo -> Reset () ; m_histo -> Sumw2 () ; }
}
// ============================================================================
// update the content 
// ============================================================================
void Ostap::Utils::TH1_WStatistic::update
( const double x ,
  const double w )
{ if ( m_histo && w
       && std::isfinite ( x )
       && std::isfinite ( w ) ) { m_histo->Fill ( x , w ) ; } }
// ============================================================================
// constructor 
// ============================================================================
Ostap::Utils::TH2_Statistic::TH2_Statistic
( TH2* histo )
  : m_histo ( histo )
{
  Ostap::Assert ( histo  && 2 == histo->GetDimension () ,
		  "Inbvald TH2 pointer"                 ,
		  "Ostap:L=:Utils::TH2_Statistic"       , 
		  INVALID_TH2  , __FILE__ , __LINE__    ) ;
}
// ============================================================================
// reset the histogram
// ============================================================================
void Ostap::Utils::TH2_Statistic::reset  ()
{ 
  ESentry senty {} ;
  if ( m_histo ) { m_histo -> Reset () ; m_histo -> Sumw2 () ; }
}
// ============================================================================
// update the content 
// ============================================================================
void Ostap::Utils::TH2_Statistic::update
( const double x , 
  const double y ) 
{ if ( m_histo
       && std::isfinite ( x )
       && std::isfinite ( y ) ) { m_histo->Fill ( x , y ) ; } }
// ============================================================================
// constructor 
// ============================================================================
Ostap::Utils::TH2_WStatistic::TH2_WStatistic
( TH2* histo )
  : m_histo ( histo )
{
  Ostap::Assert ( histo  && 2 == histo->GetDimension () ,
		  "Invald TH2 pointer"                  ,
		  "Ostap:L=:Utils::TH2_WStatistic"      , 
		  INVALID_TH2  , __FILE__ , __LINE__    ) ;
}
// ============================================================================
// reset the histogram
// ============================================================================
void Ostap::Utils::TH2_WStatistic::reset  ()
{
  ESentry senty {} ;
  if ( m_histo ) { m_histo -> Reset () ; m_histo -> Sumw2 () ; }
}
// ============================================================================
// update the content 
// ============================================================================
void Ostap::Utils::TH2_WStatistic::update
( const double x ,
  const double y , 
  const double w )
{
  if ( m_histo && w
       && std::isfinite ( x )
       && std::isfinite ( y )
       && std::isfinite ( w ) ) { m_histo->Fill ( x , y , w ) ; }
}
// ============================================================================
// constructor 
// ============================================================================
Ostap::Utils::TH3_Statistic::TH3_Statistic
( TH3* histo )
  : m_histo ( histo )
{
  Ostap::Assert ( histo  && 3 == histo->GetDimension () ,
		  "Inbvald TH3 pointer"                 ,
		  "Ostap:L=:Utils::TH3_Statistic"       , 
		  INVALID_TH3  , __FILE__ , __LINE__    ) ;
}
// ============================================================================
// reset the histogram
// ============================================================================
void Ostap::Utils::TH3_Statistic::reset  ()
{ 
  ESentry senty {} ;
  if ( m_histo ) { m_histo -> Reset () ; m_histo -> Sumw2 () ; }
}
// ============================================================================
// update the content 
// ============================================================================
void Ostap::Utils::TH3_Statistic::update
( const double x , 
  const double y , 
  const double z )
{
  if ( m_histo
       && std::isfinite ( x )
       && std::isfinite ( y )
       && std::isfinite ( z ) ) { m_histo->Fill ( x , y , z ) ; } }
// ============================================================================
// constructor 
// ============================================================================
Ostap::Utils::TH3_WStatistic::TH3_WStatistic
( TH3* histo )
  : m_histo ( histo )
{
  Ostap::Assert ( histo  && 3 == histo->GetDimension () ,
		  "Invald TH3 pointer"                  ,
		  "Ostap:L=:Utils::TH3_WStatistic"      , 
		  INVALID_TH3  , __FILE__ , __LINE__    ) ;
}
// ============================================================================
// reset the histo
// ============================================================================
void Ostap::Utils::TH3_WStatistic::reset  ()
{
  ESentry senty {} ;
  if ( m_histo ) { m_histo -> Reset () ; m_histo -> Sumw2 () ; }
}
// ============================================================================
// update the content 
// ============================================================================
void Ostap::Utils::TH3_WStatistic::update
( const double x ,
  const double y , 
  const double z , 
  const double w )
{
  if ( m_histo && w
       && std::isfinite ( x )
       && std::isfinite ( y )
       && std::isfinite ( z )
       && std::isfinite ( w ) ) { m_histo->Fill ( x , y , z , w ) ; } 
}
// ============================================================================

// ============================================================================
//                                                                      The END 
// ============================================================================
