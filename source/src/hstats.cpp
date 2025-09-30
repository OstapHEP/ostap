// =============================================================================
// Include files 
// =============================================================================
// STD& STL 
// =============================================================================
#include <limits>
#include <cmath>
// =============================================================================
// ROOT 
// =============================================================================
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
// =============================================================================
// Ostap
// =============================================================================
#include "Ostap/StatusCode.h"
#include "Ostap/Statistic.h"
// =============================================================================
// Local
// =============================================================================
#include "hstats.h"
#include "local_utils.h"
#include "status_codes.h"
// =============================================================================
Ostap::Utils::H1::H1 ( TH1* histo )
  : m_histo ( histo )
{
  Ostap::Assert ( m_histo && 1 == m_histo->GetDimension () , 
		  "Invalid TH1"                              ,
		  "Ostap::Utils::H1"                         ,
		  INVALID_TH2 , __FILE__ , __LINE__          ) ;
}
// ==============================================================================
Ostap::Utils::H2::H2 ( TH2* histo )
: m_histo ( histo )
{
  Ostap::Assert ( m_histo && 2 == m_histo->GetDimension () , 
		  "Invalid TH2"                              ,
		  "Ostap::Utils::H2"                         ,
		  INVALID_TH2 , __FILE__ , __LINE__          ) ;
}
// ==============================================================================
Ostap::Utils::H3::H3 ( TH3* histo )
:  m_histo ( histo )
{
  Ostap::Assert ( m_histo && 3 == m_histo->GetDimension () , 
		  "Invalid TH3"                              ,
		  "Ostap::Utils::H3"                         ,
		  INVALID_TH3 , __FILE__ , __LINE__          ) ;
}
// =============================================================================
Ostap::Utils::P1::P1 ( TProfile* histo )
:  m_histo ( histo )
{
  Ostap::Assert ( m_histo && 1 == m_histo->GetDimension () , 
		  "Invalid TProfile"                       ,
		  "Ostap::Utils::P1"                       ,
		  INVALID_TPROFILE , __FILE__ , __LINE__   ) ;
}
// =============================================================================
Ostap::Utils::P2::P2 ( TProfile2D* histo )
  :  m_histo ( histo )
{
  Ostap::Assert ( m_histo && 2 == m_histo->GetDimension ()   , 
		  "Invalid TProfile2D"                       ,
		  "Ostap::Utils::P2"                         ,
		  INVALID_TPROFILE2D , __FILE__ , __LINE__   ) ;
}
// =============================================================================
Ostap::Utils::P3::P3 ( TProfile3D* histo )
  :  m_histo ( histo )
{
  Ostap::Assert ( m_histo && 3 == m_histo->GetDimension ()   , 
		  "Invalid TProfile3D"                       ,
		  "Ostap::Utils::P3"                         ,
		  INVALID_TPROFILE3D , __FILE__ , __LINE__   ) ;
}
// =============================================================================


// =============================================================================
// reset 
// =============================================================================

// =============================================================================
void Ostap::Utils::H1::reset()
{
  if ( m_histo )
    {
      ESentry  sentry {} ; 
      m_histo->Reset  () ;
      if ( !m_histo->GetSumw2 () ) { m_histo->Sumw2() ; }
    }
}
// =============================================================================
void Ostap::Utils::H2::reset()
{
  if ( m_histo )
    {
      ESentry  sentry {} ; 
      m_histo->Reset  () ;
      if ( !m_histo->GetSumw2 () ) { m_histo->Sumw2() ; }
    }
}
// =============================================================================
void Ostap::Utils::H3::reset()
{
  if ( m_histo )
    {
      ESentry  sentry {} ; 
      m_histo->Reset  () ;
      if ( !m_histo->GetSumw2 () ) { m_histo->Sumw2() ; }
    }
}
// =============================================================================
void Ostap::Utils::P1::reset()
{
  if ( m_histo )
    {
      ESentry  sentry {} ; 
      m_histo->Reset  () ;
      if ( !m_histo->GetSumw2 () ) { m_histo->Sumw2() ; }
    }
}
// =============================================================================
void Ostap::Utils::P2::reset()
{
  if ( m_histo )
    {
      ESentry  sentry {} ; 
      m_histo->Reset  () ;
      if ( !m_histo->GetSumw2 () ) { m_histo->Sumw2() ; }
    }
}
// =============================================================================
void Ostap::Utils::P3::reset()
{
  if ( m_histo )
    {
      ESentry  sentry {} ; 
      m_histo->Reset  () ;
      if ( !m_histo->GetSumw2 () ) { m_histo->Sumw2() ; }
    }
}
// =============================================================================

// =============================================================================
// update 
// =============================================================================

// =============================================================================
void Ostap::Utils::H1::update
( const double x ,
  const double w )
{
  if ( m_histo && w && std::isfinite ( w )
       && std::isfinite ( x ) )
    { m_histo->Fill ( x , w ) ; }
}
// =============================================================================
void Ostap::Utils::H2::update
( const double x ,
  const double y ,
  const double w )
{
  if ( m_histo && w && std::isfinite ( w )
       && std::isfinite ( x ) 
       && std::isfinite ( y ) )
    { m_histo->Fill ( x , y ,  w ) ; }
}
// =============================================================================
void Ostap::Utils::H3::update
( const double x ,
  const double y ,
  const double z ,
  const double w )
{
  if ( m_histo && w && std::isfinite ( w )
       && std::isfinite ( x ) 
       && std::isfinite ( y ) 
       && std::isfinite ( z ) )
    { m_histo->Fill ( x , y , z , w ) ; }
}
// =============================================================================
void Ostap::Utils::P1::update
( const double x ,
  const double y ,
  const double w )
{
  if ( m_histo && w && std::isfinite ( w )
       && std::isfinite ( x ) 
       && std::isfinite ( y ) )
    { m_histo->Fill ( x , y ,  w ) ; }
}
// =============================================================================
void Ostap::Utils::P2::update
( const double x ,
  const double y ,
  const double z ,
  const double w )
{
  if ( m_histo && w && std::isfinite ( w )
       && std::isfinite ( x ) 
       && std::isfinite ( y ) 
       && std::isfinite ( z ) )
    { m_histo->Fill ( x , y , z , w ) ; }
}
// =============================================================================
void Ostap::Utils::P3::update
( const double x ,
  const double y ,
  const double z ,
  const double t ,
  const double w )
{
  if ( m_histo && w && std::isfinite ( w )
       && std::isfinite ( x ) 
       && std::isfinite ( y ) 
       && std::isfinite ( z ) 
       && std::isfinite ( t ) )
    { m_histo->Fill ( x , y , z , t , w ) ; }
}
// =============================================================================
//                                                                       The END 
// =============================================================================
