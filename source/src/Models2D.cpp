// ============================================================================
// Include files 
// ============================================================================
// STD & STL
// ============================================================================
#include <cmath>
#include <map>
#include <limits>
#include <complex>
#include <algorithm>
#include <numeric>
#include <array>
#include <tuple>
#include <functional>
// ============================================================================
// ROOT
// ============================================================================
#include "TMath.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/MoreMath.h"
#include "Ostap/Power.h"
#include "Ostap/Clenshaw.h"
#include "Ostap/Models2D.h"
// ============================================================================
//  Local 
// ============================================================================
#include "Exception.h"
#include "local_math.h"
#include "local_hash.h"
#include "Integrator1D.h"
#include "Integrator2D.h"
// ============================================================================
/** @file 
 *  Implementation file for classes  from  file Ostap/Models2D.h
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  /** @var s_CACHESIZE  
   *  the  size of cache maps 
   */
  const unsigned int s_CACHESIZE = 1000 ;
  // ==========================================================================
  // get the helper integral 
  // ==========================================================================
  class PSBERN
  {
  public:
    // ========================================================================
    PSBERN ( const Ostap::Math::PhaseSpaceNL* ps        , 
             const Ostap::Math::Bernstein*    bp        )
      : m_ps   ( ps   ) 
      , m_bp   ( bp   ) 
    {}
    PSBERN () =  delete ;
    //
    double operator() ( const double x ) const 
    { return (*m_ps)( x ) * (*m_bp)( x ) ; }
    // ==========================================================================
  private:
    // ========================================================================== 
    const Ostap::Math::PhaseSpaceNL*   m_ps   { nullptr } ; // phase space factor 
    const Ostap::Math::Bernstein*      m_bp   { nullptr } ; // bernstein polinomial 
    // ========================================================================== 
  };
  // ==========================================================================
  /// make 1D integration of the product of the phase space and the bernstein polynomial
  double _integral_
  ( const Ostap::Math::PhaseSpaceNL&    ps        , 
    const Ostap::Math::Bernstein&       bp        ,
    const double                        low       , 
    const double                        high      ,
    const Ostap::Math::WorkSpace&       work      ) 
  {
    //
    if      ( ps.highEdge() <= bp.xmin() || ps. lowEdge() >= bp.xmax() ) { return 0 ; }
    //
    if      ( s_equal ( low , high ) ) { return 0 ; }
    else if ( bp.zero()              ) { return 0 ; }
    else if ( low > high  ) { return _integral_ ( ps , bp , high , low , work ) ; }
    //
    if      ( high <= ps.lowEdge () || high <= bp.xmin () ) { return 0 ; }
    else if ( low  >= ps.highEdge() || low  >= bp.xmax () ) { return 0 ; }
    //
    const double xlow  = std::max ( std::max ( ps. lowEdge() , bp.xmin() ) , low  ) ;
    const double xhigh = std::min ( std::min ( ps.highEdge() , bp.xmax() ) , high ) ;
    //
    if ( xlow >= xhigh   ) { return 0 ; }
    //
    if ( 1 == bp.npars() ) { return bp.par(0) * ps.integral ( xlow , xhigh ) ; }
    //
    // construct the hash 
    const std::size_t tag = std::hash_combine ( bp.tag () , ps.tag () ) ;
    //
    /// integrator for class PSBERN 
    static const Ostap::Math::GSL::Integrator1D<PSBERN> s_integrator ;
    static const char message[] = "Integral(PS*Pol)" ;
    //
    const PSBERN ps_bern { &ps , &bp } ;
    const auto F     = s_integrator.make_function ( &ps_bern ) ;
    //
    // use GSL to evaluate the integral 
    //
    int    ierror    =  0   ;
    double result    =  1.0 ;
    double error     = -1.0 ;
    std::tie ( ierror , result , error ) = s_integrator.gaq_integrate_with_cache  
      ( tag                  , 
        &F                   ,   // the function
        xlow   , xhigh       ,   // low & high edges
        workspace ( work )   ,   // workspace
        s_PRECISION          ,   // absolute precision
        s_PRECISION          ,   // relative precision
        work.size ()         ,   // maximum number of subintervals
        message              ,   // message 
        __FILE__  , __LINE__ ) ; // filename & line number 
    //
    return result ;
  }
  // ==========================================================================
}
// ===========================================================================
// constructor from the order
// ===========================================================================
Ostap::Math::PS2DPol::PS2DPol
( const Ostap::Math::PhaseSpaceNL&   psx  , 
  const Ostap::Math::PhaseSpaceNL&   psy  , 
  const unsigned short               Nx   ,
  const unsigned short               Ny   )
  : m_positive ( Nx, Ny , 
                 psx.lowEdge() , psx.highEdge() , 
                 psy.lowEdge() , psy.highEdge() )
  , m_psx  ( psx   ) 
  , m_psy  ( psy   )
{}
// ===========================================================================
// constructor from the order
// ===========================================================================
Ostap::Math::PS2DPol::PS2DPol
( const Ostap::Math::PhaseSpaceNL&   psx  , 
  const Ostap::Math::PhaseSpaceNL&   psy  , 
  const unsigned short               Nx   ,
  const unsigned short               Ny   ,
  const double                       xmin , 
  const double                       xmax , 
  const double                       ymin , 
  const double                       ymax )
  : m_positive ( Nx , 
                 Ny , 
                 std::max ( psx. lowEdge () , std::min ( xmin , xmax ) ) , 
                 std::min ( psx.highEdge () , std::max ( xmin , xmax ) ) , 
                 std::max ( psy. lowEdge () , std::min ( ymin , ymax ) ) , 
                 std::min ( psy.highEdge () , std::max ( ymin , ymax ) ) ) 
  , m_psx  ( psx   ) 
  , m_psy  ( psy   )
{}
// ===========================================================================
// get the value
// ===========================================================================
double Ostap::Math::PS2DPol::operator () 
  ( const double x , 
    const double y ) const 
{
  //
  if      ( x < m_psx. lowEdge() || x < m_positive.xmin () ) { return 0 ; }
  else if ( x > m_psx.highEdge() || x > m_positive.xmax () ) { return 0 ; }
  else if ( y < m_psy. lowEdge() || y < m_positive.ymin () ) { return 0 ; }
  else if ( y > m_psy.highEdge() || y > m_positive.ymax () ) { return 0 ; }
  //
  return m_positive ( x , y ) * m_psx ( x ) * m_psy ( y ) ;
}
// ============================================================================
// helper function to make calculations
// ============================================================================
double Ostap::Math::PS2DPol::calculate
( const std::vector<double>& fx , 
  const std::vector<double>& fy ) const 
{
  double       result = 0 ;
  const Ostap::Math::Bernstein2D& b2d = m_positive.bernstein() ;
  for  ( unsigned short ix = 0 ; ix <= nX()  ; ++ix )
  {
    for  ( unsigned short iy = 0 ; iy <= nY() ; ++iy )
    { result += b2d.par ( ix , iy ) * fx[ix] * fy[iy] ; } 
  }
  //
  const double scalex = ( nX () + 1 ) / ( xmax () - xmin () ) ;
  const double scaley = ( nY () + 1 ) / ( ymax () - ymin () ) ;
  //
  return result * scalex * scaley ;
}
// ============================================================================
double Ostap::Math::PS2DPol::integral 
( const double xlow , const double xhigh , 
  const double ylow , const double yhigh ) const 
{
  //
  if      ( s_equal ( xlow , xhigh ) || s_equal ( ylow  , yhigh ) ) { return 0 ; }
  //
  else if ( xlow  > xhigh ) { return -1*integral ( xhigh , xlow  , ylow  , yhigh ) ; }
  else if ( ylow  > yhigh ) { return -1*integral ( xlow  , xhigh , yhigh , ylow  ) ; }
  //
  if      ( xhigh <  m_positive.xmin () || xhigh <  m_psx. lowEdge () ) { return 0 ; }
  else if ( xlow  >  m_positive.xmax () || xlow  >  m_psx.highEdge () ) { return 0 ; }
  else if ( yhigh <  m_positive.ymin () || yhigh <  m_psy. lowEdge () ) { return 0 ; }
  else if ( ylow  >  m_positive.ymax () || ylow  >  m_psy.highEdge () ) { return 0 ; }
  //
  const double  x_low  = std::max ( std::max ( m_psx. lowEdge() , m_positive.xmin() ) , xlow  ) ;
  const double  x_high = std::min ( std::min ( m_psx.highEdge() , m_positive.xmax() ) , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  const double  y_low  = std::max ( std::max ( m_psy. lowEdge() , m_positive.ymin() ) , ylow  ) ;
  const double  y_high = std::min ( std::min ( m_psy.highEdge() , m_positive.ymax() ) , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  const unsigned short nx  = m_positive.nX() ;
  const unsigned short ny  = m_positive.nY() ;
  //
  const Bernstein2D&   b2d = m_positive.bernstein() ;
  //
  std::vector<double> fy ( ny + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= ny ; ++i ) 
  { fy[i] = _integral_ ( m_psy , b2d.basicY ( i ) , y_low , y_high , m_workspace ) ; }
  //
  std::vector<double> fx ( nx + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nx ; ++i ) 
  { fx[i] = _integral_ ( m_psx , b2d.basicX ( i ) , x_low , x_high , m_workspace ) ; }
  //
  return calculate  ( fx  , fy ) ;
}
// ============================================================================
double Ostap::Math::PS2DPol::integrateY 
( const double x    , 
  const double ylow , const double yhigh ) const 
{
  //
  if      ( s_equal ( ylow  , yhigh ) ) { return 0 ; }
  else if ( ylow  > yhigh ) { return -1*integrateY ( x , yhigh , ylow  ) ; }
  if      ( x     <  m_positive.xmin () || x     <  m_psx. lowEdge () ) { return 0 ; }
  else if ( x     >  m_positive.xmax () || x     >  m_psx.highEdge () ) { return 0 ; }
  else if ( yhigh <  m_positive.ymin () || yhigh <  m_psy. lowEdge () ) { return 0 ; }
  else if ( ylow  >  m_positive.ymax () || ylow  >  m_psy.highEdge () ) { return 0 ; }
  //
  const double  y_low  = std::max ( std::max ( m_psy. lowEdge() , m_positive.ymin() ) , ylow  ) ;
  const double  y_high = std::min ( std::min ( m_psy.highEdge() , m_positive.ymax() ) , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  const unsigned short nx  = m_positive.nX() ;
  const unsigned short ny  = m_positive.nY() ;
  //
  const Bernstein2D&   b2d = m_positive.bernstein() ;
  //
  std::vector<double> fy ( ny + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= ny ; ++i ) 
  { fy[i] = _integral_ ( m_psy , b2d.basicY ( i ) , y_low , y_high , m_workspace ) ; }
  //
  std::vector<double> fx ( nx + 1 , 0 ) ;
  const double psx = m_psx ( x ) ;
  for  ( unsigned short i = 0 ; i <= nx ; ++i ) 
  { fx[i] = psx * b2d.basicX ( i ) ( x ) ; }
  //
  return calculate  ( fx  , fy )  ;
}
// ============================================================================
double Ostap::Math::PS2DPol::integrateX 
( const double y    , 
  const double xlow , const double xhigh ) const 
{
  //
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; }
  //
  else if ( xlow  > xhigh ) { return -1*integrateX ( y , xhigh , xlow ) ; }
  if      ( xhigh <  m_positive.xmin () || xhigh <  m_psx. lowEdge () ) { return 0 ; }
  else if ( xlow  >  m_positive.xmax () || xlow  >  m_psx.highEdge () ) { return 0 ; }
  else if ( y     <  m_positive.ymin () || y     <  m_psy. lowEdge () ) { return 0 ; }
  else if ( y     >  m_positive.ymax () || y     >  m_psy.highEdge () ) { return 0 ; }
  //
  const double  x_low  = std::max ( std::max ( m_psx. lowEdge() , m_positive.xmin() ) , xlow  ) ;
  const double  x_high = std::min ( std::min ( m_psx.highEdge() , m_positive.xmax() ) , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  const unsigned short nx  = m_positive.nX() ;
  const unsigned short ny  = m_positive.nY() ;
  //
  const Bernstein2D&   b2d = m_positive.bernstein() ;
  //
  std::vector<double> fy ( ny + 1 , 0 ) ;
  const double psy = m_psy ( y ) ;
  for ( unsigned short i = 0 ; i <= ny ; ++i ) 
  { fy[i] = psy * b2d.basicY ( i ) ( y ) ; }
  //
  std::vector<double> fx ( nx + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nx ; ++i ) 
  { fx[i] = _integral_ ( m_psx , b2d.basicX ( i ) , x_low , x_high , m_workspace ) ; }
  //
  return calculate  ( fx  , fy )  ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::PS2DPol::tag () const 
{ return std::hash_combine ( m_positive.tag () , m_psx.tag () , m_psy.tag () ) ; }
// ============================================================================









// ===========================================================================
// constructor from the order
// ===========================================================================
Ostap::Math::PS2DPolSym::PS2DPolSym
( const Ostap::Math::PhaseSpaceNL&   ps   , 
  const unsigned short               N    ) 
  : m_positive ( N, ps.lowEdge() , ps.highEdge() )
  , m_ps ( ps ) 
{}
// ===========================================================================
// constructor from the order
// ===========================================================================
Ostap::Math::PS2DPolSym::PS2DPolSym
( const Ostap::Math::PhaseSpaceNL&   ps   , 
  const unsigned short               N    ,
  const double                       xmin , 
  const double                       xmax )
  : m_positive ( N  , 
                 std::max ( ps. lowEdge () , std::min ( xmin , xmax ) ) , 
                 std::min ( ps.highEdge () , std::max ( xmin , xmax ) ) ) 
  , m_ps   ( ps   ) 
{}
// ===========================================================================
// get the value
// ===========================================================================
double Ostap::Math::PS2DPolSym::operator () 
  ( const double x , 
    const double y ) const 
{
  //
  if      ( x < m_ps. lowEdge() || x < m_positive.xmin () ) { return 0 ; }
  else if ( x > m_ps.highEdge() || x > m_positive.xmax () ) { return 0 ; }
  else if ( y < m_ps. lowEdge() || y < m_positive.ymin () ) { return 0 ; }
  else if ( y > m_ps.highEdge() || y > m_positive.ymax () ) { return 0 ; }
  //
  return m_positive ( x , y ) * m_ps ( x ) * m_ps ( y ) ;
}
// ============================================================================
// helper function to make calculations
// ============================================================================
double Ostap::Math::PS2DPolSym::calculate
( const std::vector<double>& fx , 
  const std::vector<double>& fy ) const 
{
  double       result = 0 ;
  const Ostap::Math::Bernstein2DSym& b2d = m_positive.bernstein() ;
  for  ( unsigned short ix = 0 ; ix <= nX()  ; ++ix )
  {
    result   += b2d.par ( ix , ix ) * fx[ix] * fy[ix] ;
    for  ( unsigned short iy = 0 ; iy < ix ; ++iy )
    { result += b2d.par ( ix , iy ) * ( fx[ix] * fy[iy] + fx[iy] * fy[ix] ) ; } 
  }
  //
  const double scalex = ( nX () + 1 ) / ( xmax () - xmin () ) ;
  const double scaley = scalex ;
  //
  return result * scalex * scaley ;
}
// ============================================================================
double Ostap::Math::PS2DPolSym::integral 
( const double xlow , const double xhigh , 
  const double ylow , const double yhigh ) const 
{
  //
  if      ( s_equal ( xlow , xhigh ) || s_equal ( ylow  , yhigh ) ) { return 0 ; }
  //
  else if ( xlow  > xhigh ) { return -1*integral ( xhigh , xlow  , ylow  , yhigh ) ; }
  else if ( ylow  > yhigh ) { return -1*integral ( xlow  , xhigh , yhigh , ylow  ) ; }
  //
  if      ( xhigh <  m_positive.xmin () || xhigh <  m_ps. lowEdge () ) { return 0 ; }
  else if ( xlow  >  m_positive.xmax () || xlow  >  m_ps.highEdge () ) { return 0 ; }
  else if ( yhigh <  m_positive.ymin () || yhigh <  m_ps. lowEdge () ) { return 0 ; }
  else if ( ylow  >  m_positive.ymax () || ylow  >  m_ps.highEdge () ) { return 0 ; }
  //
  const double  x_low  = std::max ( std::max ( m_ps. lowEdge() , m_positive.xmin() ) , xlow  ) ;
  const double  x_high = std::min ( std::min ( m_ps.highEdge() , m_positive.xmax() ) , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  const double  y_low  = std::max ( std::max ( m_ps. lowEdge() , m_positive.ymin() ) , ylow  ) ;
  const double  y_high = std::min ( std::min ( m_ps.highEdge() , m_positive.ymax() ) , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //  
  const unsigned short n = m_positive.n () ;
  //
  const Ostap::Math::Bernstein2DSym& b2d = m_positive.bernstein() ;
  //
  std::vector<double> fy ( n + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= n ; ++i ) 
  { fy[i] = _integral_ ( m_ps , b2d.basic( i )  , y_low , y_high , m_workspace ) ; }
  //
  std::vector<double> fx ( n + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= n ; ++i ) 
  { fx[i] = _integral_ ( m_ps , b2d.basic( i )  , x_low , x_high , m_workspace ) ; }
  //
  return calculate ( fx  , fy ) ;
}
// ============================================================================
double Ostap::Math::PS2DPolSym::integrateY
( const double x    , 
  const double ylow , const double yhigh ) const 
{
  //
  if      ( s_equal ( ylow  , yhigh ) ) { return 0 ; }
  else if ( ylow  > yhigh ) { return -1*integrateY ( x , yhigh , ylow  ) ; }
  //
  if      ( x     <  m_positive.xmin () || x     <  m_ps. lowEdge () ) { return 0 ; }
  else if ( x     >  m_positive.xmax () || x     >  m_ps.highEdge () ) { return 0 ; }
  else if ( yhigh <  m_positive.ymin () || yhigh <  m_ps. lowEdge () ) { return 0 ; }
  else if ( ylow  >  m_positive.ymax () || ylow  >  m_ps.highEdge () ) { return 0 ; }
  //
  const double  y_low  = std::max ( std::max ( m_ps. lowEdge() , m_positive.ymin() ) , ylow  ) ;
  const double  y_high = std::min ( std::min ( m_ps.highEdge() , m_positive.ymax() ) , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  const unsigned short n = m_positive.n () ;
  //
  const Ostap::Math::Bernstein2DSym& b2d = m_positive.bernstein() ;
  //
  std::vector<double> fy ( n + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= n ; ++i ) 
  { fy[i] = _integral_ ( m_ps , b2d.basic( i )  , y_low , y_high , m_workspace ) ; }
  //
  std::vector<double> fx ( n + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= n ; ++i ) 
  { fx[i] = m_ps ( x ) * b2d.basic( i ) ( x )  ; }
  //
  return calculate ( fx  , fy ) ;
}
// ============================================================================
double Ostap::Math::PS2DPolSym::integrateX
( const double y                         , 
  const double xlow , const double xhigh ) const 
{ return integrateY ( y , xlow , xhigh ) ; }
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::PS2DPolSym::tag () const 
{ return std::hash_combine ( m_positive.tag () , m_ps.tag () ) ; }
// ============================================================================


// ===========================================================================
// constructor from the order
// ===========================================================================
Ostap::Math::PS2DPol2::PS2DPol2
( const Ostap::Math::PhaseSpaceNL&   psx  , 
  const Ostap::Math::PhaseSpaceNL&   psy  , 
  const double                       mmax ,
  const unsigned short               Nx   ,
  const unsigned short               Ny   )
  : m_positive ( Nx, Ny , 
                 psx.lowEdge () , psx.highEdge () , 
                 psy.lowEdge () , psy.highEdge () )
  , m_psx     ( psx                                        )
  , m_psy     ( psy                                        )
  , m_mmax    ( psx.lowEdge  () + psy.lowEdge  () < mmax ? 
                mmax : 
                psx.highEdge () + psy.highEdge () )
  , m_psx_aux ( psx.lowEdge  () , psx.highEdge () , psx.L () , psx.N () -  1 ) 
  , m_psy_aux ( psy.lowEdge  () , psy.highEdge () , psy.L () , psy.N () -  1 ) 
{}
// ===========================================================================
// constructor from the order
// ===========================================================================
Ostap::Math::PS2DPol2::PS2DPol2
( const Ostap::Math::PhaseSpaceNL&   psx  , 
  const Ostap::Math::PhaseSpaceNL&   psy  , 
  const double                       mmax ,
  const unsigned short               Nx   ,
  const unsigned short               Ny   ,
  const double                       xmin , 
  const double                       xmax , 
  const double                       ymin , 
  const double                       ymax )
  : m_positive ( Nx , 
                 Ny , 
                 std::max ( psx. lowEdge () , std::min ( xmin , xmax ) ) , 
                 std::min ( psx.highEdge () , std::max ( xmin , xmax ) ) , 
                 std::max ( psy. lowEdge () , std::min ( ymin , ymax ) ) , 
                 std::min ( psy.highEdge () , std::max ( ymin , ymax ) ) ) 
  , m_psx  ( psx   ) 
  , m_psy  ( psy   )
  , m_mmax ( psx.lowEdge() + psy.lowEdge() < mmax ? mmax : psx.highEdge() + psy.highEdge() )
  , m_psx_aux ( psx.lowEdge() , psx.highEdge() , psx.L () , psx.N () - 1 ) 
  , m_psy_aux ( psy.lowEdge() , psy.highEdge() , psy.L () , psy.N () - 1 ) 
{}
// ===========================================================================
// get the value
// ===========================================================================
double Ostap::Math::PS2DPol2::operator () 
  ( const double x , 
    const double y ) const 
{
  //
  if      ( x < m_psx. lowEdge() || x < m_positive.xmin () ) { return 0 ; }
  else if ( x > m_psx.highEdge() || x > m_positive.xmax () ) { return 0 ; }
  else if ( y < m_psy. lowEdge() || y < m_positive.ymin () ) { return 0 ; }
  else if ( y > m_psy.highEdge() || y > m_positive.ymax () ) { return 0 ; }
  //
  if ( x + y > m_mmax ) { return 0 ; }
  //
  m_psx_aux.setThresholds ( m_psx.lowEdge() , m_mmax - y ) ;
  m_psy_aux.setThresholds ( m_psy.lowEdge() , m_mmax - x ) ;
  //
  return m_positive ( x , y ) * 
    0.5 * ( m_psx ( x ) * m_psy_aux ( y ) + m_psy ( y ) * m_psx_aux ( x ) )  ;  
}
// ============================================================================
double Ostap::Math::PS2DPol2::integral 
( const double xlow , const double xhigh , 
  const double ylow , const double yhigh ) const 
{
  //
  if      ( s_equal ( xlow , xhigh ) || s_equal ( ylow  , yhigh ) ) { return 0 ; }
  //
  else if ( xlow  > xhigh ) { return -1*integral ( xhigh , xlow  , ylow  , yhigh ) ; }
  else if ( ylow  > yhigh ) { return -1*integral ( xlow  , xhigh , yhigh , ylow  ) ; }
  //
  if      ( xhigh <  m_positive.xmin () || xhigh <  m_psx. lowEdge () ) { return 0 ; }
  else if ( xlow  >  m_positive.xmax () || xlow  >  m_psx.highEdge () ) { return 0 ; }
  else if ( yhigh <  m_positive.ymin () || yhigh <  m_psy. lowEdge () ) { return 0 ; }
  else if ( ylow  >  m_positive.ymax () || ylow  >  m_psy.highEdge () ) { return 0 ; }
  //
  const double  x_low  = std::max ( std::max ( m_psx. lowEdge() , m_positive.xmin() ) , xlow  ) ;
  const double  x_high = std::min ( std::min ( m_psx.highEdge() , m_positive.xmax() ) , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  const double  y_low  = std::max ( std::max ( m_psy. lowEdge() , m_positive.ymin() ) , ylow  ) ;
  const double  y_high = std::min ( std::min ( m_psy.highEdge() , m_positive.ymax() ) , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  if ( x_low + y_low >= m_mmax ) { return 0 ; }
  //
  // use cubature   
  static const Ostap::Math::GSL::Integrator2D<PS2DPol2> s_cubature{} ;
  static const char s_message[] = "Integral(PS2DPol2)" ;
  const auto F = s_cubature.make_function ( this , x_low , x_high , y_low , y_high ) ;
  //
  int    ierror  =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  std::tie ( ierror , result , error ) = s_cubature.cubature_with_cache 
    ( tag () , &F , 20000 , s_PRECISION , s_PRECISION , s_message , __FILE__ , __LINE__ ) ;
  // ==========================================================================
  return result ;
}
// ============================================================================
double Ostap::Math::PS2DPol2::integrateY 
( const double x    , 
  const double ylow , const double yhigh ) const 
{
  //
  if      ( s_equal ( ylow  , yhigh ) ) { return 0 ; }
  else if ( ylow  > yhigh ) { return -1*integrateY ( x , yhigh , ylow  ) ; }
  if      ( x     <  m_positive.xmin () || x     <  m_psx. lowEdge () ) { return 0 ; }
  else if ( x     >  m_positive.xmax () || x     >  m_psx.highEdge () ) { return 0 ; }
  else if ( yhigh <  m_positive.ymin () || yhigh <  m_psy. lowEdge () ) { return 0 ; }
  else if ( ylow  >  m_positive.ymax () || ylow  >  m_psy.highEdge () ) { return 0 ; }
  //
  const double  y_low  = std::max ( std::max ( m_psy. lowEdge() , m_positive.ymin() ) , ylow  ) ;
  const double  y_high = std::min ( std::min ( m_psy.highEdge() , m_positive.ymax() ) , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  typedef Ostap::Math::IntegrateY<PS2DPol2> IY ;
  static const Ostap::Math::GSL::Integrator1D<IY> s_integrator ;
  static const char message[] = "IntegrateY(PS2DPol2)" ;
  //
  const IY fy { this , x } ;
  const auto F = s_integrator.make_function ( &fy ) ;
  //  
  int    ierror    =  0   ;
  double result    =  1.0 ;
  double error     = -1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate_with_cache  
    (std::hash_combine ( tag() , 'X' , x )                   , 
      &F                        ,   // the function
      y_low   , y_high          ,   // low & high edges
      workspace ( m_workspace ) ,   // workspace
      s_PRECISION               ,   // absolute precision
      s_PRECISION               ,   // relative precision
      m_workspace.size()                    ,   // maximum number of subintervals
      message                   ,   // message 
      __FILE__  , __LINE__      ) ; // filename & line number 
  //
  return result ;
}
// ============================================================================
double Ostap::Math::PS2DPol2::integrateX 
( const double y    , 
  const double xlow , const double xhigh ) const 
{
  //
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; }
  //
  else if ( xlow  > xhigh ) { return -1*integrateX ( y , xhigh , xlow ) ; }
  if      ( xhigh <  m_positive.xmin () || xhigh <  m_psx. lowEdge () ) { return 0 ; }
  else if ( xlow  >  m_positive.xmax () || xlow  >  m_psx.highEdge () ) { return 0 ; }
  else if ( y     <  m_positive.ymin () || y     <  m_psy. lowEdge () ) { return 0 ; }
  else if ( y     >  m_positive.ymax () || y     >  m_psy.highEdge () ) { return 0 ; }
  //
  const double  x_low  = std::max ( std::max ( m_psx. lowEdge() , m_positive.xmin() ) , xlow  ) ;
  const double  x_high = std::min ( std::min ( m_psx.highEdge() , m_positive.xmax() ) , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  //
  typedef Ostap::Math::IntegrateX<PS2DPol2> IX ;
  static const Ostap::Math::GSL::Integrator1D<IX> s_integrator ;
  static const char message[] = "IntegrateX(PS2DPol2)" ;
  //
  const IX fx { this , y } ;
  const auto F = s_integrator.make_function ( &fx ) ;
  //  
  int    ierror    =  0   ;
  double result    =  1.0 ;
  double error     = -1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate_with_cache  
    (std::hash_combine ( tag() , 'Y' , y )                   , 
      &F                        ,   // the function
      x_low   , x_high          ,   // low & high edges
      workspace ( m_workspace ) ,   // workspace
      s_PRECISION               ,   // absolute precision
      s_PRECISION               ,   // relative precision
      m_workspace.size()                    ,   // maximum number of subintervals
      message                   ,   // message 
      __FILE__  , __LINE__      ) ; // filename & line number 
  //
  return result ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::PS2DPol2::tag () const 
{ return std::hash_combine ( m_positive. tag () , 
                             m_psx     . tag () , 
                             m_psy     . tag () , m_mmax ) ; }
// ============================================================================

// ===========================================================================
// constructor from the order
// ===========================================================================
Ostap::Math::PS2DPol2Sym::PS2DPol2Sym
( const Ostap::Math::PhaseSpaceNL&   ps   , 
  const double                       mmax ,  
  const unsigned short               N    ) 
  : m_positive ( N, ps.lowEdge() , ps.highEdge() )
  , m_ps      ( ps ) 
  , m_mmax    ( 2 * ps.lowEdge() < mmax ? mmax : 2 * ps.highEdge()   )
  , m_psx_aux ( ps.lowEdge() , ps.highEdge() , ps.L () , ps.N () - 1 ) 
  , m_psy_aux ( ps.lowEdge() , ps.highEdge() , ps.L () , ps.N () - 1 ) 
{}
// ===========================================================================
// constructor from the order
// ===========================================================================
Ostap::Math::PS2DPol2Sym::PS2DPol2Sym
( const Ostap::Math::PhaseSpaceNL&   ps   , 
  const double                       mmax ,  
  const unsigned short               N    ,
  const double                       xmin , 
  const double                       xmax )
  : m_positive ( N  , 
                 std::max ( ps. lowEdge () , std::min ( xmin , xmax ) ) , 
                 std::min ( ps.highEdge () , std::max ( xmin , xmax ) ) ) 
  , m_ps      ( ps   ) 
  , m_mmax    ( 2 * ps.lowEdge() < mmax ? mmax : 2 * ps.highEdge()   )
  , m_psx_aux ( ps.lowEdge() , ps.highEdge() , ps.L () , ps.N () - 1 ) 
  , m_psy_aux ( ps.lowEdge() , ps.highEdge() , ps.L () , ps.N () - 1 ) 
{}
// ===========================================================================
// get the value
// ===========================================================================
double Ostap::Math::PS2DPol2Sym::operator () 
  ( const double x , 
    const double y ) const 
{
  //
  if      ( x < m_ps. lowEdge() || x < m_positive.xmin () ) { return 0 ; }
  else if ( x > m_ps.highEdge() || x > m_positive.xmax () ) { return 0 ; }
  else if ( y < m_ps. lowEdge() || y < m_positive.ymin () ) { return 0 ; }
  else if ( y > m_ps.highEdge() || y > m_positive.ymax () ) { return 0 ; }
  //
  if ( x + y > m_mmax ) { return 0 ; }
  //
  m_psx_aux.setThresholds ( m_ps.lowEdge() , m_mmax - y ) ;
  m_psy_aux.setThresholds ( m_ps.lowEdge() , m_mmax - x ) ;
  //
  return m_positive ( x , y ) * 0.5 * 
    ( m_ps ( y ) * m_psx_aux ( x ) + m_ps ( x ) * m_psy_aux ( y ) ) ;
}
// ============================================================================
double Ostap::Math::PS2DPol2Sym::integral 
( const double xlow , const double xhigh , 
  const double ylow , const double yhigh ) const 
{
  //
  if      ( s_equal ( xlow , xhigh ) || s_equal ( ylow  , yhigh ) ) { return 0 ; }
  //
  else if ( xlow  > xhigh ) { return -1*integral ( xhigh , xlow  , ylow  , yhigh ) ; }
  else if ( ylow  > yhigh ) { return -1*integral ( xlow  , xhigh , yhigh , ylow  ) ; }
  //
  if      ( xhigh <  m_positive.xmin () || xhigh <  m_ps. lowEdge () ) { return 0 ; }
  else if ( xlow  >  m_positive.xmax () || xlow  >  m_ps.highEdge () ) { return 0 ; }
  else if ( yhigh <  m_positive.ymin () || yhigh <  m_ps. lowEdge () ) { return 0 ; }
  else if ( ylow  >  m_positive.ymax () || ylow  >  m_ps.highEdge () ) { return 0 ; }
  //
  const double  x_low  = std::max ( std::max ( m_ps. lowEdge() , m_positive.xmin() ) , xlow  ) ;
  const double  x_high = std::min ( std::min ( m_ps.highEdge() , m_positive.xmax() ) , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  const double  y_low  = std::max ( std::max ( m_ps. lowEdge() , m_positive.ymin() ) , ylow  ) ;
  const double  y_high = std::min ( std::min ( m_ps.highEdge() , m_positive.ymax() ) , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  if ( x_low + y_low >= m_mmax ) { return 0 ; }
  //
  // use cubature   
  static const Ostap::Math::GSL::Integrator2D<PS2DPol2Sym> s_cubature{} ;
  static const char s_message[] = "Integral(PS2DPol2Sym)" ;
  const auto F = s_cubature.make_function ( this , x_low , x_high , y_low , y_high ) ;
  //
  int    ierror  =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  std::tie ( ierror , result , error ) = s_cubature.cubature_with_cache 
    ( tag () , &F , 20000 , s_PRECISION , s_PRECISION , s_message , __FILE__ , __LINE__ ) ;
  // ==========================================================================
  return  result ;
}
// ============================================================================
double Ostap::Math::PS2DPol2Sym::integrateY
( const double x    , 
  const double ylow , const double yhigh ) const 
{
  //
  if      ( s_equal ( ylow  , yhigh ) ) { return 0 ; }
  else if ( ylow  > yhigh ) { return -1*integrateY ( x , yhigh , ylow  ) ; }
  //
  if      ( x     <  m_positive.xmin () || x     <  m_ps. lowEdge () ) { return 0 ; }
  else if ( x     >  m_positive.xmax () || x     >  m_ps.highEdge () ) { return 0 ; }
  else if ( yhigh <  m_positive.ymin () || yhigh <  m_ps. lowEdge () ) { return 0 ; }
  else if ( ylow  >  m_positive.ymax () || ylow  >  m_ps.highEdge () ) { return 0 ; }
  //
  const double  y_low  = std::max ( std::max ( m_ps. lowEdge() , m_positive.ymin() ) , ylow  ) ;
  const double  y_high = std::min ( std::min ( m_ps.highEdge() , m_positive.ymax() ) , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  typedef Ostap::Math::IntegrateY<PS2DPol2Sym> IY ;
  static const Ostap::Math::GSL::Integrator1D<IY> s_integrator ;
  static const char message[] = "IntegrateY(PS2DPol2Sym)" ;
  //
  const IY fy { this , x } ;
  const auto F = s_integrator.make_function ( &fy ) ;
  //  
  int    ierror    =  0   ;
  double result    =  1.0 ;
  double error     = -1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate_with_cache  
    ( std::hash_combine ( tag() , 'X' , x )                   , 
      &F                        ,   // the function
      y_low   , y_high          ,   // low & high edges
      workspace ( m_workspace ) ,   // workspace
      s_PRECISION               ,   // absolute precision
      s_PRECISION               ,   // relative precision
      m_workspace.size()                    ,   // maximum number of subintervals
      message                   ,   // message 
      __FILE__  , __LINE__      ) ; // filename & line number 
  // ==========================================================================
  return result ;
}
// ============================================================================
double Ostap::Math::PS2DPol2Sym::integrateX
( const double y                         , 
  const double xlow , const double xhigh ) const 
{ return integrateY ( y , xlow , xhigh ) ; }
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::PS2DPol2Sym::tag () const 
{ return std::hash_combine ( m_positive.tag()  , m_ps.tag () , m_mmax ) ; }
// ============================================================================









// ===========================================================================
// constructor from the order
// ===========================================================================
Ostap::Math::PS2DPol3::PS2DPol3
( const Ostap::Math::PhaseSpaceNL&   psx  , 
  const Ostap::Math::PhaseSpaceNL&   psy  , 
  const double                       mmax ,
  const unsigned short               Nx   ,
  const unsigned short               Ny   )
  : m_psx     ( psx  , Nx )
  , m_psy     ( psy  , Ny )
  , m_mmax    ( psx.lowEdge() + psy.lowEdge() < mmax ? mmax : psx.highEdge() + psy.highEdge() )
  , m_psx_aux ( psx.lowEdge() , psx.highEdge() , psx.L () , psx.N () - 1 ) 
  , m_psy_aux ( psy.lowEdge() , psy.highEdge() , psy.L () , psy.N () - 1 ) 
{}
// ===========================================================================
// constructor from the order
// ===========================================================================
Ostap::Math::PS2DPol3::PS2DPol3
( const Ostap::Math::PhaseSpacePol&  psx  , 
  const Ostap::Math::PhaseSpacePol&  psy  , 
  const double                       mmax )
  : m_psx     ( psx )
  , m_psy     ( psy )
  , m_mmax    ( psx->lowEdge () + psy->lowEdge () < mmax ? mmax : 
                psx->highEdge() + psy->highEdge() )
  , m_psx_aux ( psx->lowEdge ()     , 
                psx->highEdge()     , 
                psx->L       ()     , 
                psx->N       () - 1 ) 
  , m_psy_aux ( psy->lowEdge ()     , 
                psy->highEdge()     , 
                psy->L       ()     , 
                psy->N       () - 1 ) 
{}
// ===========================================================================
// constructor from the order
// ===========================================================================
Ostap::Math::PS2DPol3::PS2DPol3
( const Ostap::Math::PhaseSpaceNL&   psx  , 
  const Ostap::Math::PhaseSpaceNL&   psy  , 
  const double                       mmax ,
  const unsigned short               Nx   ,
  const unsigned short               Ny   ,
  const double                       xmin , 
  const double                       xmax , 
  const double                       ymin , 
  const double                       ymax )
  : m_psx  ( psx   , Nx ,
             std::max ( psx. lowEdge () , std::min ( xmin , xmax ) ) , 
             std::min ( psx.highEdge () , std::max ( xmin , xmax ) ) )
  , m_psy  ( psy   , Ny , 
             std::max ( psy. lowEdge () , std::min ( ymin , ymax ) ) , 
             std::min ( psy.highEdge () , std::max ( ymin , ymax ) ) )
  , m_mmax ( psx.lowEdge() + psy.lowEdge() < mmax ? mmax : psx.highEdge() + psy.highEdge() )
  , m_psx_aux ( psx.lowEdge() , psx.highEdge() , psx.L () , psx.N () - 1 ) 
  , m_psy_aux ( psy.lowEdge() , psy.highEdge() , psy.L () , psy.N () - 1 ) 
{}
// ===========================================================================
// get parameters/phases 
std::vector<double> Ostap::Math::PS2DPol3::pars() const 
{
  std::vector<double> result ( npars() , 0.0 ) ;
  const unsigned short nPx = m_psx.npars() ;
  const unsigned short nPy = m_psy.npars() ;
  //
  for ( unsigned short ix = 0 ; ix < nPx ; ++ix )
  { result [       ix ] = m_psx.par ( ix ) ; }
  for ( unsigned short iy = 0 ; iy < nPy ; ++iy )
  { result [ nPx + iy ] = m_psy.par ( iy ) ; }
  //
  // const std::vector<double>& xps = xpars () ;
  // const std::vector<double>& yps = ypars () ;
  // std::copy ( xps.begin() , xps.end() , result.begin()              ) ;
  // std::copy ( yps.begin() , yps.end() , result.begin() + xps.size() ) ;
  return result ;
}
// ===========================================================================
// get the value
// ===========================================================================
double Ostap::Math::PS2DPol3::operator () 
  ( const double x , 
    const double y ) const 
{
  //
  if      ( x < xmin () || x > xmax() ) {  return 0 ; }
  else if ( y < ymin () || y > ymax() ) {  return 0 ; }
  //
  if ( x + y > m_mmax ) { return 0 ; }
  //
  m_psx_aux.setThresholds ( m_psx->lowEdge() , m_mmax - y ) ;
  m_psy_aux.setThresholds ( m_psy->lowEdge() , m_mmax - x ) ;
  //
  return 0.5 * ( m_psx ( x ) * m_psy_aux ( y ) + m_psy ( y ) * m_psx_aux ( x ) )  ;  
}
// ============================================================================
double Ostap::Math::PS2DPol3::integral 
( const double xlow , const double xhigh , 
  const double ylow , const double yhigh ) const 
{
  //
  if      ( s_equal ( xlow , xhigh ) || s_equal ( ylow  , yhigh ) ) { return 0 ; }
  //
  else if ( xlow  > xhigh ) { return -1*integral ( xhigh , xlow  , ylow  , yhigh ) ; }
  else if ( ylow  > yhigh ) { return -1*integral ( xlow  , xhigh , yhigh , ylow  ) ; }
  //
  if      ( xhigh < xmin () || xlow > xmax() ) { return 0 ; }
  else if ( yhigh < ymin () || ylow > ymax() ) { return 0 ; }
  //
  const double  x_low  = std::max ( xmin () , xlow  ) ;
  const double  x_high = std::min ( xmax () , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  const double  y_low  = std::max ( ymin () , ylow  ) ;
  const double  y_high = std::min ( ymax () , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  if ( x_low + y_low >= m_mmax ) { return 0 ; }
  //
  // use cubature   
  static const Ostap::Math::GSL::Integrator2D<PS2DPol3> s_cubature{} ;
  static const char s_message[] = "Integral(PS2DPol3)" ;
  const auto F = s_cubature.make_function ( this , x_low , x_high , y_low , y_high ) ;
  //
  int    ierror  =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  std::tie ( ierror , result , error ) = s_cubature.cubature_with_cache  
    ( tag () , &F , 20000 , s_PRECISION , s_PRECISION , s_message , __FILE__ , __LINE__ ) ;
  // ==========================================================================
  return result ;
}
// ============================================================================
double Ostap::Math::PS2DPol3::integrateY 
( const double x    , 
  const double ylow , const double yhigh ) const 
{
  //
  if      ( s_equal ( ylow  , yhigh ) ) { return 0 ; }
  else if ( ylow  > yhigh )             { return -1*integrateY ( x , yhigh , ylow  ) ; }
  //
  if      ( x     <  xmin () || x    > xmax() ) { return 0 ; }
  else if ( yhigh <  ymin () || ylow > ymax() ) { return 0 ; }
  //
  const double  y_low  = std::max ( ymin () , ylow  ) ;
  const double  y_high = std::min ( ymax () , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  typedef Ostap::Math::IntegrateY<PS2DPol3> IY ;
  static const Ostap::Math::GSL::Integrator1D<IY> s_integrator ;
  static const char message[] = "IntegrateY(PS2DPol3)" ;
  //
  const IY fy { this , x } ;
  const auto F = s_integrator.make_function ( &fy ) ;
  //  
  int    ierror    =  0   ;
  double result    =  1.0 ;
  double error     = -1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate_with_cache   
    ( tag  ()                   , 
      &F                        ,   // the function
      y_low   , y_high          ,   // low & high edges
      workspace ( m_workspace ) ,   // workspace
      s_PRECISION               ,   // absolute precision
      s_PRECISION               ,   // relative precision
      m_workspace.size()                    ,   // maximum number of subintervals
      message                   ,   // message 
      __FILE__  , __LINE__      ) ; // filename & line number 
  // ==========================================================================
  return result ;
}
// ============================================================================
double Ostap::Math::PS2DPol3::integrateX 
( const double y    , 
  const double xlow , const double xhigh ) const 
{
  //
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; }
  else if ( xlow  > xhigh ) { return -1*integrateX ( y , xhigh , xlow ) ; }
  // 
  if      ( y     <  ymin () || y    > ymax () ) { return 0 ; }
  if      ( xhigh <  xmin () || xlow > xmax () ) { return 0 ; }
  //
  const double  x_low  = std::max ( xmin () , xlow  ) ;
  const double  x_high = std::min ( xmax () , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  typedef Ostap::Math::IntegrateX<PS2DPol3> IX ;
  static const Ostap::Math::GSL::Integrator1D<IX> s_integrator ;
  static const char message[] = "IntegrateX(PS2DPol3)" ;
  //
  const IX fx { this , y } ;
  const auto F = s_integrator.make_function ( &fx ) ;
  //  
  int    ierror    =  0   ;
  double result    =  1.0 ;
  double error     = -1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate_with_cache  
    (std::hash_combine ( tag() , 'Y' , y )                   , 
      &F                        ,   // the function
      x_low   , x_high          ,   // low & high edges
      workspace ( m_workspace ) ,   // workspace
      s_PRECISION               ,   // absolute precision
      s_PRECISION               ,   // relative precision
      m_workspace.size()                    ,   // maximum number of subintervals
      message                   ,   // message 
      __FILE__  , __LINE__      ) ; // filename & line number 
  // ==========================================================================
  return result ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::PS2DPol3::tag () const 
{ return std::hash_combine ( m_psx.tag () , m_psy.tag () , m_mmax ) ; }
// ============================================================================


// ===========================================================================
// constructor from the order
// ===========================================================================
Ostap::Math::PS2DPol3Sym::PS2DPol3Sym
( const Ostap::Math::PhaseSpaceNL&   ps   , 
  const double                       mmax ,  
  const unsigned short               N    ) 
  : m_ps      ( ps ,  N ) 
  , m_mmax    ( 2 * ps.lowEdge() < mmax ? mmax : 2 * ps.highEdge()   )
  , m_psx_aux ( ps.lowEdge() , ps.highEdge() , ps.L () , ps.N () - 1 ) 
  , m_psy_aux ( ps.lowEdge() , ps.highEdge() , ps.L () , ps.N () - 1 ) 
{}
// ===========================================================================
Ostap::Math::PS2DPol3Sym::PS2DPol3Sym
( const Ostap::Math::PhaseSpacePol&  ps   ,
  const double                       mmax ) 
 : m_ps      ( ps ) 
 , m_mmax    ( 2 * ps->lowEdge() < mmax ? mmax : 2 * ps->highEdge()   )
 , m_psx_aux ( ps->lowEdge() , ps->highEdge() , ps->L () , ps->N () - 1 ) 
 , m_psy_aux ( ps->lowEdge() , ps->highEdge() , ps->L () , ps->N () - 1 ) 
{}
// ===========================================================================
// constructor from the order
// ===========================================================================
Ostap::Math::PS2DPol3Sym::PS2DPol3Sym
( const Ostap::Math::PhaseSpaceNL&   ps   , 
  const double                       mmax ,  
  const unsigned short               N    ,
  const double                       xmin , 
  const double                       xmax )
  : m_ps      ( ps , N ,  xmin , xmax ) 
  , m_mmax    ( 2 * ps.lowEdge() < mmax ? mmax : 2 * ps.highEdge()   )
  , m_psx_aux ( ps.lowEdge() , ps.highEdge() , ps.L () , ps.N () - 1 ) 
  , m_psy_aux ( ps.lowEdge() , ps.highEdge() , ps.L () , ps.N () - 1 ) 
{}
// ===========================================================================
// get the value
// ===========================================================================
double Ostap::Math::PS2DPol3Sym::operator () 
  ( const double x , 
    const double y ) const 
{
  //
  if      ( x < xmin () || x > xmax () ) { return 0 ; }
  else if ( y < ymin () || y > ymax () ) { return 0 ; }
  //
  if ( x + y > m_mmax ) { return 0 ; }
  //
  m_psx_aux.setThresholds ( m_ps->lowEdge() , m_mmax - y ) ;
  m_psy_aux.setThresholds ( m_ps->lowEdge() , m_mmax - x ) ;
  //
  return 0.5 * 
    ( m_ps ( y ) * m_psx_aux ( x ) + m_ps ( x ) * m_psy_aux ( y ) ) ;
}
// ============================================================================
double Ostap::Math::PS2DPol3Sym::integral 
( const double xlow , const double xhigh , 
  const double ylow , const double yhigh ) const 
{
  //
  if      ( s_equal ( xlow , xhigh ) || s_equal ( ylow  , yhigh ) ) { return 0 ; }
  //
  else if ( xlow  > xhigh ) { return -1*integral ( xhigh , xlow  , ylow  , yhigh ) ; }
  else if ( ylow  > yhigh ) { return -1*integral ( xlow  , xhigh , yhigh , ylow  ) ; }
  //
  if      ( xhigh <  xmin () || xlow > xmax () ) { return 0 ; }
  else if ( yhigh <  ymin () || ylow > ymax () ) { return 0 ; }
  //
  const double  x_low  = std::max ( xmin () , xlow  ) ;
  const double  x_high = std::min ( xmax () , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  const double  y_low  = std::max ( ymin () , ylow  ) ;
  const double  y_high = std::min ( ymax () , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  if ( x_low + y_low >= m_mmax ) { return 0 ; }
  //
  // use cubature   
  static const Ostap::Math::GSL::Integrator2D<PS2DPol3Sym> s_cubature{} ;
  static const char s_message[] = "Integral(PS2DPol3Sym)" ;
  const auto F = s_cubature.make_function ( this , x_low , x_high , y_low , y_high ) ;
  //
  int    ierror  =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  std::tie ( ierror , result , error ) = s_cubature.cubature_with_cache 
    ( tag () , &F , 20000 , s_PRECISION , s_PRECISION , s_message , __FILE__ , __LINE__ ) ;
  // ==========================================================================
  return  result ;
}
// ============================================================================
double Ostap::Math::PS2DPol3Sym::integrateY
( const double x    , 
  const double ylow , const double yhigh ) const 
{
  //
  if      ( s_equal ( ylow  , yhigh ) ) { return 0 ; }
  else if ( ylow  > yhigh ) { return -1*integrateY ( x , yhigh , ylow  ) ; }
  //
  if      ( x     <  xmin () || x     > xmax() ) { return 0 ; }
  else if ( ylow  >  ymax () || yhigh < ymin() ) { return 0 ; }
  //
  const double  y_low  = std::max ( ymin () , ylow  ) ;
  const double  y_high = std::min ( ymax (), yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  typedef Ostap::Math::IntegrateY<PS2DPol3Sym> IY ;
  static const Ostap::Math::GSL::Integrator1D<IY> s_integrator ;
  static const char message[] = "IntegrateY(PS2DPol3Sym)" ;
  //
  const IY fy { this , x } ;
  const auto F = s_integrator.make_function ( &fy ) ;
  //  
  int    ierror    =  0   ;
  double result    =  1.0 ;
  double error     = -1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate_with_cache
    ( std::hash_combine ( tag() , 'X' , x )                  , 
      &F                        ,   // the function
      y_low   , y_high          ,   // low & high edges
      workspace ( m_workspace ) ,   // workspace
      s_PRECISION               ,   // absolute precision
      s_PRECISION               ,   // relative precision
      m_workspace.size()                    ,   // maximum number of subintervals
      message                   ,   // message 
      __FILE__  , __LINE__      ) ; // filename & line number 
  // ==========================================================================
  return result ;
}
// ============================================================================
double Ostap::Math::PS2DPol3Sym::integrateX
( const double y                         , 
  const double xlow , const double xhigh ) const 
{ return integrateY ( y , xlow , xhigh ) ; }
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::PS2DPol3Sym::tag () const 
{ return std::hash_combine ( m_ps.tag () , m_mmax ) ; }
// ============================================================================






// ===========================================================================
// constructor from the order
// ===========================================================================
Ostap::Math::ExpoPS2DPol::ExpoPS2DPol
( const Ostap::Math::PhaseSpaceNL&   psy  , 
  const double                       xmin , 
  const double                       xmax , 
  const unsigned short               Nx   ,
  const unsigned short               Ny   ,
  const double                       tau  )      
  : m_positive ( Nx , 
                 Ny , 
                 std::min ( xmin , xmax ) , std::max ( xmin , xmax ) ,
                 psy.lowEdge()            , psy.highEdge()           )
  , m_psy ( psy )
  , m_tau ( tau ) 
{}
// ===========================================================================
// constructor from the order
// ===========================================================================
Ostap::Math::ExpoPS2DPol::ExpoPS2DPol
( const Ostap::Math::PhaseSpaceNL&   psy  , 
  const double                       xmin , 
  const double                       xmax , 
  const unsigned short               Nx   ,
  const unsigned short               Ny   ,
  const double                       ymin , 
  const double                       ymax ,
  const double                       tau  )      
  : m_positive ( Nx , 
                 Ny , 
                 std::min ( xmin , xmax )   , std::max ( xmin , xmax ) ,
                 std::max ( psy. lowEdge () , std::min ( ymin , ymax ) ) , 
                 std::min ( psy.highEdge () , std::max ( ymin , ymax ) ) )
  , m_psy ( psy )
  , m_tau ( tau ) 
{}
// ============================================================================
// set tau-parameter
// ============================================================================
bool Ostap::Math::ExpoPS2DPol::setTau ( const double value )
{
  //
  if ( s_equal ( m_tau , value ) ) { return false ; }
  //
  m_tau = value ;
  //
  return true ;
}
// ===========================================================================
// get the value
// ===========================================================================
double Ostap::Math::ExpoPS2DPol::operator () 
  ( const double x , 
    const double y ) const 
{
  //
  if      ( x < m_positive.xmin () || x > m_positive.xmax () ) { return 0 ; }
  else if ( y < m_psy. lowEdge  () || y < m_positive.ymin () ) { return 0 ; }
  else if ( y > m_psy.highEdge  () || y > m_positive.ymax () ) { return 0 ; }
  //
  return m_positive ( x , y ) * my_exp ( m_tau * x ) * m_psy ( y ) ;
}
// ============================================================================
// helper function to make calculations
// ============================================================================
double Ostap::Math::ExpoPS2DPol::calculate
( const std::vector<double>& fx , 
  const std::vector<double>& fy ) const 
{
  double       result = 0 ;
  const Ostap::Math::Bernstein2D& b2d = m_positive.bernstein() ;
  for  ( unsigned short ix = 0 ; ix <= nX()  ; ++ix )
  {
    for  ( unsigned short iy = 0 ; iy <= nY() ; ++iy )
    { result += b2d.par ( ix , iy ) * fx[ix] * fy[iy] ; } 
  }
  //
  const double scalex = ( nX () + 1 ) / ( xmax () - xmin () ) ;
  const double scaley = ( nY () + 1 ) / ( ymax () - ymin () ) ;
  //
  return result * scalex * scaley ;
}
// ============================================================================
double Ostap::Math::ExpoPS2DPol::integral 
( const double xlow , const double xhigh , 
  const double ylow , const double yhigh ) const 
{
  //
  if      ( s_equal ( xlow , xhigh ) || s_equal ( ylow  , yhigh ) ) { return 0 ; }
  //
  else if ( xlow  > xhigh ) { return -1*integral ( xhigh , xlow  , ylow  , yhigh ) ; }
  else if ( ylow  > yhigh ) { return -1*integral ( xlow  , xhigh , yhigh , ylow  ) ; }
  //
  if      ( xhigh <  m_positive.xmin () ) { return 0 ; }
  else if ( xlow  >  m_positive.xmax () ) { return 0 ; }
  else if ( yhigh <  m_positive.ymin () || yhigh <  m_psy. lowEdge () ) { return 0 ; }
  else if ( ylow  >  m_positive.ymax () || ylow  >  m_psy.highEdge () ) { return 0 ; }
  //
  const double  x_low  = std::max ( m_positive.xmin() , xlow  ) ; 
  const double  x_high = std::min ( m_positive.xmax() , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  const double  y_low  = std::max ( std::max ( m_psy. lowEdge() , m_positive.ymin() ) , ylow  ) ;
  const double  y_high = std::min ( std::min ( m_psy.highEdge() , m_positive.ymax() ) , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  const unsigned short nx  = m_positive.nX() ;
  const unsigned short ny  = m_positive.nY() ;
  //
  const Bernstein2D&   b2d = m_positive.bernstein() ;
  //
  std::vector<double> fy ( ny + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= ny ; ++i ) 
  { fy[i] = _integral_ ( m_psy , b2d.basicY ( i ) , y_low , y_high , m_workspace ) ; }
  //
  std::vector<double> fx ( nx + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nx ; ++i )
  { fx[i] = Ostap::Math::integrate ( b2d.basicX ( i ) , m_tau , x_low , x_high ) ; }
  //
  return calculate ( fx  , fy ) ;
}
// ============================================================================
double Ostap::Math::ExpoPS2DPol::integrateY
( const double x    , 
  const double ylow , const double yhigh ) const 
{
  //
  if      ( s_equal ( ylow  , yhigh ) ) { return 0 ; }
  else if ( ylow  > yhigh ) { return -1*integrateY ( x , yhigh , ylow  ) ; }
  //
  if      ( x     <  m_positive.xmin () ) { return 0 ; }
  else if ( x     >  m_positive.xmax () ) { return 0 ; }
  else if ( yhigh <  m_positive.ymin () || yhigh <  m_psy. lowEdge () ) { return 0 ; }
  else if ( ylow  >  m_positive.ymax () || ylow  >  m_psy.highEdge () ) { return 0 ; }
  //
  const double  y_low  = std::max ( std::max ( m_psy. lowEdge() , m_positive.ymin() ) , ylow  ) ;
  const double  y_high = std::min ( std::min ( m_psy.highEdge() , m_positive.ymax() ) , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  const unsigned short nx  = m_positive.nX() ;
  const unsigned short ny  = m_positive.nY() ;
  //
  const Bernstein2D&   b2d = m_positive.bernstein() ;
  //
  std::vector<double> fy ( ny + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= ny ; ++i ) 
  { fy[i] = _integral_ ( m_psy , b2d.basicY ( i ) , y_low , y_high , m_workspace ) ; }
  //
  std::vector<double> fx ( nx + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nx ; ++i )
  { fx[i] = b2d.basicX ( i ) ( x ) * my_exp (  m_tau * x ) ; }
  //
  return calculate ( fx  , fy ) ;
}
// ============================================================================
double Ostap::Math::ExpoPS2DPol::integrateX 
( const double y    ,
  const double xlow , const double xhigh ) const 
{
  //
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; }
  //
  else if ( xlow  > xhigh ) { return -1*integrateX ( y , xhigh , xlow ) ; }
  //
  if      ( xhigh <  m_positive.xmin () ) { return 0 ; }
  else if ( xlow  >  m_positive.xmax () ) { return 0 ; }
  else if ( y     <  m_positive.ymin () || y <  m_psy. lowEdge () ) { return 0 ; }
  else if ( y     >  m_positive.ymax () || y >  m_psy.highEdge () ) { return 0 ; }
  //
  const double  x_low  = std::max ( m_positive.xmin() , xlow  ) ; 
  const double  x_high = std::min ( m_positive.xmax() , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  const unsigned short nx  = m_positive.nX() ;
  const unsigned short ny  = m_positive.nY() ;
  //
  const Bernstein2D&   b2d = m_positive.bernstein() ;
  //
  std::vector<double> fy ( ny + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= ny ; ++i ) 
  { fy[i] = m_psy ( y ) *  b2d.basicY ( i ) ( y ) ; }
  //
  std::vector<double> fx ( nx + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nx ; ++i )
  { fx[i] = Ostap::Math::integrate ( b2d.basicX ( i ) , m_tau , x_low , x_high ) ; }
  //
  return calculate ( fx  , fy ) ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::ExpoPS2DPol::tag () const 
{ return std::hash_combine ( m_positive.tag() , m_psy.tag () , m_tau ) ; }
// ============================================================================


// ===========================================================================
// constructor from the order
// ===========================================================================
Ostap::Math::Expo2DPol::Expo2DPol
( const double                       xmin , 
  const double                       xmax , 
  const double                       ymin , 
  const double                       ymax , 
  const unsigned short               Nx   ,
  const unsigned short               Ny   ,
  const double                       taux ,
  const double                       tauy )
  : m_positive ( Nx , Ny , 
                 std::min ( xmin , xmax ) , std::max ( xmin , xmax ) ,
                 std::min ( ymin , ymax ) , std::max ( ymin , ymax ) ) 
  , m_tauX ( taux ) 
  , m_tauY ( tauy )
{}
// ===========================================================================
// set tau-parameter
// ============================================================================
bool Ostap::Math::Expo2DPol::setTauX ( const double value )
{
  //
  if ( s_equal ( m_tauX , value ) ) { return false ; }
  //
  m_tauX = value ;
  //
  return true ;
}
// ===========================================================================
// set tau-parameter
// ============================================================================
bool Ostap::Math::Expo2DPol::setTauY ( const double value )
{
  //
  if ( s_equal ( m_tauY , value ) ) { return false ; }
  //
  m_tauY = value ;
  //
  return true ;
}
// ===========================================================================
// get the value
// ===========================================================================
double Ostap::Math::Expo2DPol::operator () 
  ( const double x , const double y ) const 
{
  //
  if      ( x < m_positive.xmin () || x > m_positive.xmax () ) { return 0 ; }
  else if ( y < m_positive.ymin () || y > m_positive.ymax () ) { return 0 ; }
  //
  return m_positive ( x , y ) * my_exp ( m_tauX * x ) * my_exp ( m_tauY * y ) ;
}
// ============================================================================
// helper function to make calculations
// ============================================================================
double Ostap::Math::Expo2DPol::calculate
( const std::vector<double>& fx , 
  const std::vector<double>& fy ) const 
{
  double       result = 0 ;
  const Ostap::Math::Bernstein2D& b2d = m_positive.bernstein() ;
  for  ( unsigned short ix = 0 ; ix <= nX()  ; ++ix )
  {
    for  ( unsigned short iy = 0 ; iy <= nY() ; ++iy )
    { result += b2d.par ( ix , iy ) * fx[ix] * fy[iy] ; } 
  }
  //
  const double scalex = ( nX () + 1 ) / ( xmax () - xmin () ) ;
  const double scaley = ( nY () + 1 ) / ( ymax () - ymin () ) ;
  //
  return result * scalex * scaley ;
}
// ============================================================================
double Ostap::Math::Expo2DPol::integral 
( const double xlow , const double xhigh , 
  const double ylow , const double yhigh ) const 
{
  //
  if      ( s_equal ( xlow , xhigh ) || s_equal ( ylow  , yhigh ) ) { return 0 ; }
  //
  else if ( xlow  > xhigh ) { return -1*integral ( xhigh , xlow  , ylow  , yhigh ) ; }
  else if ( ylow  > yhigh ) { return -1*integral ( xlow  , xhigh , yhigh , ylow  ) ; }
  //
  if      ( xhigh <  m_positive.xmin () ) { return 0 ; }
  else if ( xlow  >  m_positive.xmax () ) { return 0 ; }
  else if ( yhigh <  m_positive.ymin () ) { return 0 ; }
  else if ( ylow  >  m_positive.ymax () ) { return 0 ; }
  //
  const double  x_low  = std::max ( m_positive.xmin() , xlow  ) ; 
  const double  x_high = std::min ( m_positive.xmax() , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  const double  y_low  = std::max ( m_positive.ymin() , ylow  ) ;
  const double  y_high = std::min ( m_positive.ymax() , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  const unsigned short nx  = m_positive.nX() ;
  const unsigned short ny  = m_positive.nY() ;
  //
  const Bernstein2D&   b2d = m_positive.bernstein() ;
  //
  std::vector<double> fy ( ny + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= ny ; ++i ) 
  { fy[i] =  Ostap::Math::integrate ( b2d.basicY ( i ) , m_tauY , y_low , y_high ) ; }
  //
  std::vector<double> fx ( nx + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nx ; ++i ) 
  { fx[i] =  Ostap::Math::integrate ( b2d.basicX ( i ) , m_tauX , x_low , x_high  ) ; }
  //
  return calculate  ( fx  , fy ) ;
}
// ============================================================================
double Ostap::Math::Expo2DPol::integrateY  
( const double x    , 
  const double ylow , const double yhigh ) const 
{
  //
  if      ( s_equal ( ylow  , yhigh ) ) { return 0 ; }
  else if ( ylow  > yhigh ) { return -1*integrateY ( x , yhigh , ylow  ) ; }
  //
  if      ( x     <  m_positive.xmin () ) { return 0 ; }
  else if ( x     >  m_positive.xmax () ) { return 0 ; }
  else if ( yhigh <  m_positive.ymin () ) { return 0 ; }
  else if ( ylow  >  m_positive.ymax () ) { return 0 ; }
  //
  const double  y_low  = std::max ( m_positive.ymin() , ylow  ) ;
  const double  y_high = std::min ( m_positive.ymax() , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  const unsigned short nx  = m_positive.nX() ;
  const unsigned short ny  = m_positive.nY() ;
  //
  const Bernstein2D&   b2d = m_positive.bernstein() ;
  //
  std::vector<double> fy ( ny + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= ny ; ++i ) 
  { fy[i] =  Ostap::Math::integrate ( b2d.basicY ( i ) , m_tauY , y_low , y_high ) ; }
  //
  std::vector<double> fx ( nx + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nx ; ++i ) 
  { fx[i] =  my_exp ( m_tauX * x ) * b2d.basicX ( i ) ( x ) ; }
  //
  return calculate  ( fx  , fy ) ;
}
// ============================================================================
double Ostap::Math::Expo2DPol::integrateX
( const double y    , 
  const double xlow , const double xhigh ) const 
{
  //
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; }
  else if ( xlow  > xhigh ) { return -1*integrateX ( y , xhigh , xlow ) ; }
  //
  if      ( xhigh <  m_positive.xmin () ) { return 0 ; }
  else if ( xlow  >  m_positive.xmax () ) { return 0 ; }
  else if ( y     <  m_positive.ymin () ) { return 0 ; }
  else if ( y     >  m_positive.ymax () ) { return 0 ; }
  //
  const double  x_low  = std::max ( m_positive.xmin() , xlow  ) ; 
  const double  x_high = std::min ( m_positive.xmax() , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  const unsigned short nx  = m_positive.nX() ;
  const unsigned short ny  = m_positive.nY() ;
  //
  const Bernstein2D&   b2d = m_positive.bernstein() ;
  //
  std::vector<double> fy ( ny + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= ny ; ++i ) 
  { fy[i] =  my_exp ( m_tauY * y ) * b2d.basicY ( i ) ( y )  ; }
  //
  std::vector<double> fx ( nx + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nx ; ++i ) 
  { fx[i] =  Ostap::Math::integrate ( b2d.basicX ( i ) , m_tauX , x_low , x_high  ) ; }
  //
  return calculate  ( fx  , fy ) ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Expo2DPol::tag () const 
{ return std::hash_combine ( m_positive.tag () , m_tauX , m_tauY ) ; }
// ============================================================================

// ===========================================================================
// constructor from the order
// ===========================================================================
Ostap::Math::Expo2DPolSym::Expo2DPolSym
( const double                       xmin , 
  const double                       xmax , 
  const unsigned short               N    ,
  const double                       tau  )      
  : m_positive ( N , std::min ( xmin , xmax ) , std::max ( xmin , xmax ) )
  , m_tau ( tau ) 
{}
// ===========================================================================
// set tau-parameter
// ============================================================================
bool Ostap::Math::Expo2DPolSym::setTau ( const double value )
{
  //
  if ( s_equal ( m_tau , value ) ) { return false ; }
  //
  m_tau = value ;
  //
  return true ;
}
// ===========================================================================
// get the value
// ===========================================================================
double Ostap::Math::Expo2DPolSym::operator () 
  ( const double x , const double y ) const 
{
  //
  if      ( x < m_positive.xmin () || x > m_positive.xmax () ) { return 0 ; }
  else if ( y < m_positive.ymin () || y > m_positive.ymax () ) { return 0 ; }
  //
  return m_positive ( x , y ) * my_exp ( m_tau * ( x + y ) ) ;
}
// ============================================================================
// helper function to make calculations
// ============================================================================
double Ostap::Math::Expo2DPolSym::calculate
( const std::vector<double>& fx , 
  const std::vector<double>& fy ) const 
{
  double       result = 0 ;
  const Ostap::Math::Bernstein2DSym& b2d = m_positive.bernstein() ;
  for  ( unsigned short ix = 0 ; ix <= nX()  ; ++ix )
  {
    result   += b2d.par ( ix , ix ) * fx[ix] * fy[ix] ;
    for  ( unsigned short iy = 0 ; iy < ix ; ++iy )
    { result += b2d.par ( ix , iy ) * ( fx[ix] * fy[iy] + fx[iy] * fy[ix] ) ; } 
  }
  //
  const double scalex = ( nX () + 1 ) / ( xmax () - xmin () ) ;
  const double scaley = scalex ;
  //
  return result * scalex * scaley ;
}
// ============================================================================
double Ostap::Math::Expo2DPolSym::integral 
( const double xlow , const double xhigh , 
  const double ylow , const double yhigh ) const 
{
  //
  if      ( s_equal ( xlow , xhigh ) || s_equal ( ylow  , yhigh ) ) { return 0 ; }
  //
  else if ( xlow  > xhigh ) { return -1*integral ( xhigh , xlow  , ylow  , yhigh ) ; }
  else if ( ylow  > yhigh ) { return -1*integral ( xlow  , xhigh , yhigh , ylow  ) ; }
  //
  if      ( xhigh <  m_positive.xmin () ) { return 0 ; }
  else if ( xlow  >  m_positive.xmax () ) { return 0 ; }
  else if ( yhigh <  m_positive.ymin () ) { return 0 ; }
  else if ( ylow  >  m_positive.ymax () ) { return 0 ; }
  //
  const double  x_low  = std::max ( m_positive.xmin() , xlow  ) ; 
  const double  x_high = std::min ( m_positive.xmax() , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  const double  y_low  = std::max ( m_positive.ymin() , ylow  ) ;
  const double  y_high = std::min ( m_positive.ymax() , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  const unsigned short nx  = m_positive.nX() ;
  const unsigned short ny  = m_positive.nY() ;
  //
  const Bernstein2DSym&   b2d = m_positive.bernstein() ;
  //
  std::vector<double> fy ( ny + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= ny ; ++i ) 
  { fy[i] = Ostap::Math::integrate ( b2d.basic ( i ) , m_tau , y_low , y_high ) ; }
  //
  std::vector<double> fx ( nx + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nx ; ++i )
  { fx[i] = Ostap::Math::integrate ( b2d.basic ( i ) , m_tau , x_low , x_high ) ; }
  //
  return calculate ( fx  , fy ) ;
}
// ============================================================================
double Ostap::Math::Expo2DPolSym::integrateY 
( const double x    , 
  const double ylow , const double yhigh ) const 
{
  //
  if      ( s_equal ( ylow  , yhigh ) ) { return 0 ; }
  else if ( ylow  > yhigh ) { return -1*integrateY ( x , yhigh , ylow  ) ; }
  //
  if      ( x     <  m_positive.xmin () ) { return 0 ; }
  else if ( x     >  m_positive.xmax () ) { return 0 ; }
  else if ( yhigh <  m_positive.ymin () ) { return 0 ; }
  else if ( ylow  >  m_positive.ymax () ) { return 0 ; }
  //
  const double  y_low  = std::max ( m_positive.ymin() , ylow  ) ;
  const double  y_high = std::min ( m_positive.ymax() , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  const unsigned short nx  = m_positive.nX() ;
  const unsigned short ny  = m_positive.nY() ;
  //
  const Bernstein2DSym&   b2d = m_positive.bernstein() ;
  //
  std::vector<double> fy ( ny + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= ny ; ++i ) 
  { fy[i] = Ostap::Math::integrate ( b2d.basic ( i ) , m_tau , y_low , y_high ) ; }
  //
  std::vector<double> fx ( nx + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nx ; ++i )
  { fx[i] = my_exp ( m_tau * x ) * b2d.basic ( i ) ( x ) ; }
  //
  return calculate ( fx  , fy ) ;
}
// ============================================================================
double Ostap::Math::Expo2DPolSym::integrateX
( const double y    , 
  const double xlow , const double xhigh ) const 
{ return integrateY ( y , xlow , xhigh ) ; }
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Expo2DPolSym::tag () const 
{ return std::hash_combine ( m_positive.tag () , m_tau ) ; }
// ============================================================================



// ============================================================================
// The END 
// ============================================================================
