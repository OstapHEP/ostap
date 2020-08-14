// ============================================================================
// Include files 
// ============================================================================
//  STD&STL
// ============================================================================
#include <utility>
#include <algorithm>
#include <limits>
// ============================================================================
// GSL
// ============================================================================
#include "gsl/gsl_rstat.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/P2Quantile.h"
// ============================================================================
// local
// ============================================================================
#include "Exception.h"
// ============================================================================
/** @file  
 *  Implementation file for class Ostap::Math::P2Quantile
 *  @date 2020-05-21 
 *  @author Vanya BELYAEV Ivan/Belyaev@itep.ru
 */
// ============================================================================
namespace
{
  const    std::string s_INVALID  = "P2Quantile: invalid quantile!" ;
  constexpr double     s_INFINITY = -1 * std::numeric_limits<double>::max ()  ;
}
// ============================================================================
/* Standard constructor for quantile
 *  @param p quatile \f$ 0 < p > 1 \f$
 */
// ============================================================================
Ostap::Math::GSL::P2Quantile::P2Quantile
( const double p )
  : m_ws   ( nullptr )
  , m_p    ( p )  
{
  Ostap::Assert ( 0 < p && p < 1 , s_INVALID ) ;
}
// ======================================================================
// copy contructor 
// ======================================================================
Ostap::Math::GSL::P2Quantile::P2Quantile
( const Ostap::Math::GSL::P2Quantile& right )
  : m_ws   ( nullptr   )
  , m_p    ( right.m_p )  
{
  if ( right.m_ws )
  {
    m_ws = gsl_rstat_quantile_alloc ( m_p ) ;
    //
    std::copy ( right.m_ws->q    , right.m_ws->q    + 5 , m_ws->q    ) ;
    std::copy ( right.m_ws->npos , right.m_ws->npos + 5 , m_ws->npos ) ;
    std::copy ( right.m_ws->np   , right.m_ws->np   + 5 , m_ws->np   ) ;
    std::copy ( right.m_ws->dnp  , right.m_ws->dnp  + 5 , m_ws->dnp  ) ;
    m_ws -> n = right.m_ws-> n ;
  }
}
// ============================================================================
// move constructor
// ============================================================================
Ostap::Math::GSL::P2Quantile::P2Quantile
( Ostap::Math::GSL::P2Quantile&& right )
  : m_ws   ( nullptr   )
  , m_p    ( right.m_p )  
{
  std::swap ( m_ws  , right.m_ws  ) ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::GSL::P2Quantile::~P2Quantile()
{ if  ( m_ws ) { gsl_rstat_quantile_free ( m_ws ) ; m_ws = nullptr ; } }
// ============================================================================
// copy assignement
// ============================================================================
Ostap::Math::GSL::P2Quantile&
Ostap::Math::GSL::P2Quantile::operator=
( const Ostap::Math::GSL::P2Quantile&  right )
{
  if ( &right == this  ) { return *this ; }
  if ( m_ws != nullptr ) { gsl_rstat_quantile_free ( m_ws ) ; m_ws = nullptr ; }
  m_p = right.m_p ;
  //
  if ( right.m_ws )
  {
    m_ws = gsl_rstat_quantile_alloc ( m_p ) ;
    //
    std::copy ( right.m_ws->q    , right.m_ws->q    + 5 , m_ws->q    ) ;
    std::copy ( right.m_ws->npos , right.m_ws->npos + 5 , m_ws->npos ) ;
    std::copy ( right.m_ws->np   , right.m_ws->np   + 5 , m_ws->np   ) ;
    std::copy ( right.m_ws->dnp  , right.m_ws->dnp  + 5 , m_ws->dnp  ) ;
    //
    m_ws -> n = right.m_ws-> n ;
  }
  return *this ;
}
// ============================================================================
// move  assignement
// ============================================================================
Ostap::Math::GSL::P2Quantile&
Ostap::Math::GSL::P2Quantile::operator=
(       Ostap::Math::GSL::P2Quantile&& right )
{
  if ( &right == this  ) { return *this ; }
  std::swap ( m_ws , right.m_ws ) ;
  std::swap ( m_p  , right.m_p  ) ;
  return *this ;
}
// ============================================================================
// get the quantile value
// ============================================================================
double Ostap::Math::GSL::P2Quantile::value() const
{ return m_ws == nullptr ? s_INFINITY : gsl_rstat_quantile_get ( m_ws ) ; }
// ============================================================================
//                                                                      The END  
// ============================================================================
