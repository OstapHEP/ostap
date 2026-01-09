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
#include "Ostap/Hash.h"
#include "Ostap/Math.h"
#include "Ostap/MoreMath.h"
#include "Ostap/qMath.h"
#include "Ostap/Power.h"
#include "Ostap/Clenshaw.h"
#include "Ostap/StatusCode.h"
#include "Ostap/Models2D.h"
// ============================================================================
//  Local 
// ============================================================================
#include "local_math.h"
#include "local_hash.h"
#include "Integrator1D.h"
#include "Integrator2D.h"
#include "status_codes.h"
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
    PSBERN
    ( const Ostap::Math::PhaseSpaceNL* ps        , 
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
    const std::size_t tag = Ostap::Utils::hash_combiner ( bp.tag () , ps.tag () ) ;
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
    std::tie ( ierror , result , error ) = s_integrator.qag_integrate
      ( tag                  , 
        &F                   ,   // the function
        xlow   , xhigh       ,   // low & high edges
        workspace ( work )   ,   // workspace
        s_APRECISION         ,   // absolute precision
        s_RPRECISION         ,   // relative precision
        work.size ()         ,   // maximum number of subintervals
        message              ,   // message 
        __FILE__  , __LINE__ ) ; // filename & line number 
    //
    return result ;
  }
  // ==========================================================================
  typedef  std::array<short,9> SPLITS ;
  /// split point for 2D Gaussia
  static const SPLITS s_SPLITS = { -10 , -5 , -3 , -1 , 0 , 1 , 3 , 5 , 10 } ;
  // =========================================================================
}
// ===========================================================================
// constructor from the order
// ===========================================================================
Ostap::Math::PS2DPol::PS2DPol
( const Ostap::Math::PhaseSpaceNL&   psx  , 
  const Ostap::Math::PhaseSpaceNL&   psy  , 
  const unsigned short               Nx   ,
  const unsigned short               Ny   )
  : Ostap::Math::PolyFactor2D ( Nx, Ny , 
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
  : Ostap::Math::PolyFactor2D ( Nx, Ny , 
				std::max ( psx. lowEdge () , std::min ( xmin , xmax ) ) , 
				std::min ( psx.highEdge () , std::max ( xmin , xmax ) ) , 
				std::max ( psy. lowEdge () , std::min ( ymin , ymax ) ) , 
				std::min ( psy.highEdge () , std::max ( ymin , ymax ) ) ) 
  , m_psx  ( psx   ) 
  , m_psy  ( psy   )
{}
// ===========================================================================
Ostap::Math::PS2DPol::PS2DPol
( const Ostap::Math::Positive2D& pol ,
  const PhaseSpaceNL&            psx ,
  const PhaseSpaceNL&            psy ) 
  : Ostap::Math::PolyFactor2D ( pol ) 
  , m_psx  ( psx   ) 
  , m_psy  ( psy   )
{
  Ostap::Assert ( m_psx.lowEdge   () < m_positive.xmax () , 
                  "Invalid setting of lowEdge/highEdge/xmin/xmax"   ,
                  "Ostap::Math::PS2DPol" ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;  
  Ostap::Assert ( m_positive.xmin () < m_psx.highEdge() , 
                  "Invalid setting of lowEdge/highEdge/xmin/xmax"   ,
                  "Ostap::Math::PS2DPol" , 
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  Ostap::Assert ( m_psy.lowEdge   () < m_positive.ymax () , 
                  "Invalid setting of lowEdge/highEdge/ymin/ymax"   ,
                  "Ostap::Math::PS2DPol" , 
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  Ostap::Assert ( m_positive.ymin () < m_psy.highEdge() , 
                  "Invalid setting of lowEdge/highEdge/ymin/ymax"   ,
                  "Ostap::Math::PS2DPol" , 
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
}
// ===========================================================================
// get the value
// ===========================================================================
double Ostap::Math::PS2DPol::evaluate 
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
{ return Ostap::Utils::hash_combiner ( m_positive.tag () , m_psx.tag () , m_psy.tag () ) ; }
// ============================================================================


// ===========================================================================
// constructor from the order
// ===========================================================================
Ostap::Math::PS2DPolSym::PS2DPolSym
( const Ostap::Math::PhaseSpaceNL&   ps   , 
  const unsigned short               N    )
  : Ostap::Math::PolyFactor2DSym ( N, ps.lowEdge() , ps.highEdge() )
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
  : Ostap::Math::PolyFactor2DSym
    ( N  , 
      std::max ( ps. lowEdge () , std::min ( xmin , xmax ) ) , 
      std::min ( ps.highEdge () , std::max ( xmin , xmax ) ) ) 
  , m_ps   ( ps   ) 
{}
// ===========================================================================
Ostap::Math::PS2DPolSym::PS2DPolSym
( const Ostap::Math::Positive2DSym& pol ,
  const PhaseSpaceNL&               ps  ) 
  : Ostap::Math::PolyFactor2DSym ( pol ) 
  , m_ps       ( ps   ) 
{
  Ostap::Assert ( m_ps.lowEdge    () < m_positive.xmax () , 
                  "Invalid setting of lowEdge/highEdge/xmin/xmax"   ,
                  "Ostap::Math::PS2DPolSym"    ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  Ostap::Assert ( m_positive.xmin () < m_ps.highEdge   () , 
                  "Invalid setting of lowEdge/highEdge/xmin/xmax"   ,
                  "Ostap::Math::PS2DPolSym"    , 
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
}
// ===========================================================================
// get the value
// ===========================================================================
double Ostap::Math::PS2DPolSym::evaluate 
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
{ return Ostap::Utils::hash_combiner ( m_positive.tag () , m_ps.tag () ) ; }
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
  : Ostap::Math::PolyFactor2D ( Nx, Ny , 
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
  : Ostap::Math::PolyFactor2D ( Nx , 
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
Ostap::Math::PS2DPol2::PS2DPol2
( const Ostap::Math::Positive2D& pol  ,
  const PhaseSpaceNL&            psx  ,
  const PhaseSpaceNL&            psy  , 
  const double                   mmax )
  : Ostap::Math::PolyFactor2D ( pol ) 
  , m_psx  ( psx   ) 
  , m_psy  ( psy   )
  , m_mmax ( psx.lowEdge() + psy.lowEdge() < mmax ? mmax : psx.highEdge() + psy.highEdge() )
  , m_psx_aux ( psx.lowEdge() , psx.highEdge() , psx.L () , psx.N () - 1 ) 
  , m_psy_aux ( psy.lowEdge() , psy.highEdge() , psy.L () , psy.N () - 1 ) 
{
  Ostap::Assert ( m_psx.lowEdge   () < m_positive.xmax () , 
                  "Invalid setting of lowEdge/highEdge/xmin/xmax"   ,
                  "Ostap::Math::PS2DPol2" , 
  		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;  
  Ostap::Assert ( m_positive.xmin () < m_psx.highEdge() , 
                  "Invalid setting of lowEdge/highEdge/xmin/xmax"   ,
                  "Ostap::Math::PS2DPol2" , 
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;  
  Ostap::Assert ( m_psy.lowEdge   () < m_positive.ymax () , 
                  "Invalid setting of lowEdge/highEdge/ymin/ymax"   ,
                  "Ostap::Math::PS2DPol2" , 
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;  
  Ostap::Assert ( m_positive.ymin () < m_psy.highEdge() , 
                  "Invalid setting of lowEdge/highEdge/ymin/ymax"   ,
                  "Ostap::Math::PS2DPol2" , 
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;  
}
// ===========================================================================
// get the value
// ===========================================================================
double Ostap::Math::PS2DPol2::evaluate 
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
  std::tie ( ierror , result , error ) = s_cubature.cubature
    ( tag () , &F , 20000 , s_APRECISION , s_RPRECISION , s_message , __FILE__ , __LINE__ ) ;
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
  typedef Ostap::Math::IntegrateY2<PS2DPol2> IY ;
  static const Ostap::Math::GSL::Integrator1D<IY> s_integrator ;
  static const char message[] = "IntegrateY2(PS2DPol2)" ;
  //
  const IY fy { this , x } ;
  const auto F = s_integrator.make_function ( &fy ) ;
  //  
  int    ierror    =  0   ;
  double result    =  1.0 ;
  double error     = -1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    (Ostap::Utils::hash_combiner ( tag() , 'X' , x )                   , 
      &F                        ,   // the function
      y_low   , y_high          ,   // low & high edges
      workspace ( m_workspace ) ,   // workspace
      s_APRECISION               ,   // absolute precision
      s_RPRECISION               ,   // relative precision
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
  typedef Ostap::Math::IntegrateX2<PS2DPol2> IX ;
  static const Ostap::Math::GSL::Integrator1D<IX> s_integrator ;
  static const char message[] = "IntegrateX2(PS2DPol2)" ;
  //
  const IX fx { this , y } ;
  const auto F = s_integrator.make_function ( &fx ) ;
  //  
  int    ierror    =  0   ;
  double result    =  1.0 ;
  double error     = -1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( Ostap::Utils::hash_combiner ( tag() , 'Y' , y )                   , 
      &F                        ,   // the function
      x_low   , x_high          ,   // low & high edges
      workspace ( m_workspace ) ,   // workspace
      s_APRECISION               ,   // absolute precision
      s_RPRECISION               ,   // relative precision
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
{ return Ostap::Utils::hash_combiner ( m_positive. tag () , 
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
  : Ostap::Math::PolyFactor2DSym ( N, ps.lowEdge() , ps.highEdge() )
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
  : Ostap::Math::PolyFactor2DSym ( N ,
				   std::max ( ps. lowEdge () , std::min ( xmin , xmax ) ) , 
				   std::min ( ps.highEdge () , std::max ( xmin , xmax ) ) ) 
  , m_ps      ( ps   ) 
  , m_mmax    ( 2 * ps.lowEdge() < mmax ? mmax : 2 * ps.highEdge()   )
  , m_psx_aux ( ps.lowEdge() , ps.highEdge() , ps.L () , ps.N () - 1 ) 
  , m_psy_aux ( ps.lowEdge() , ps.highEdge() , ps.L () , ps.N () - 1 ) 
{}
// ===========================================================================
Ostap::Math::PS2DPol2Sym::PS2DPol2Sym
( const Ostap::Math::Positive2DSym& pol ,
  const PhaseSpaceNL&               ps  , 
  const double                      mmax )
  : Ostap::Math::PolyFactor2DSym ( pol ) 
  , m_ps       ( ps   ) 
  , m_mmax    ( 2 * ps.lowEdge() < mmax ? mmax : 2 * ps.highEdge()   )
  , m_psx_aux ( ps.lowEdge() , ps.highEdge() , ps.L () , ps.N () - 1 ) 
  , m_psy_aux ( ps.lowEdge() , ps.highEdge() , ps.L () , ps.N () - 1 ) 
{
  Ostap::Assert ( m_ps.lowEdge    () < m_positive.xmax () , 
                  "Invalid setting of lowEdge/highEdge/xmin/xmax"   ,
                  "Ostap::Math::PS2DPol2Sym"   , 
  		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;  
  Ostap::Assert ( m_positive.xmin () < m_ps.highEdge   () , 
                  "Invalid setting of lowEdge/highEdge/xmin/xmax"   ,
                  "Ostap::Math::PS2DPol2Sym"   , 
  		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;  
}
// ===========================================================================
// get the value
// ===========================================================================
double Ostap::Math::PS2DPol2Sym::evaluate 
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
  std::tie ( ierror , result , error ) = s_cubature.cubature
    ( tag () , &F , 20000 , s_APRECISION , s_RPRECISION , s_message , __FILE__ , __LINE__ ) ;
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
  typedef Ostap::Math::IntegrateY2<PS2DPol2Sym> IY ;
  static const Ostap::Math::GSL::Integrator1D<IY> s_integrator ;
  static const char message[] = "IntegrateY2(PS2DPol2Sym)" ;
  //
  const IY fy { this , x } ;
  const auto F = s_integrator.make_function ( &fy ) ;
  //  
  int    ierror    =  0   ;
  double result    =  1.0 ;
  double error     = -1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( Ostap::Utils::hash_combiner ( tag() , 'X' , x )                   , 
      &F                        ,   // the function
      y_low   , y_high          ,   // low & high edges
      workspace ( m_workspace ) ,   // workspace
      s_APRECISION               ,   // absolute precision
      s_RPRECISION               ,   // relative precision
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
{ return Ostap::Utils::hash_combiner ( m_positive.tag()  , m_ps.tag () , m_mmax ) ; }
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
double Ostap::Math::PS2DPol3::evaluate 
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
  std::tie ( ierror , result , error ) = s_cubature.cubature
    ( tag () , &F , 20000 , s_APRECISION , s_RPRECISION , s_message , __FILE__ , __LINE__ ) ;
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
  typedef Ostap::Math::IntegrateY2<PS2DPol3> IY ;
  static const Ostap::Math::GSL::Integrator1D<IY> s_integrator ;
  static const char message[] = "IntegrateY2(PS2DPol3)" ;
  //
  const IY fy { this , x } ;
  const auto F = s_integrator.make_function ( &fy ) ;
  //  
  int    ierror    =  0   ;
  double result    =  1.0 ;
  double error     = -1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag  ()                   , 
      &F                        ,   // the function
      y_low   , y_high          ,   // low & high edges
      workspace ( m_workspace ) ,   // workspace
      s_APRECISION               ,   // absolute precision
      s_RPRECISION               ,   // relative precision
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
  typedef Ostap::Math::IntegrateX2<PS2DPol3> IX ;
  static const Ostap::Math::GSL::Integrator1D<IX> s_integrator ;
  static const char message[] = "IntegrateX2(PS2DPol3)" ;
  //
  const IX fx { this , y } ;
  const auto F = s_integrator.make_function ( &fx ) ;
  //  
  int    ierror    =  0   ;
  double result    =  1.0 ;
  double error     = -1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    (Ostap::Utils::hash_combiner ( tag() , 'Y' , y )                   , 
      &F                        ,   // the function
      x_low   , x_high          ,   // low & high edges
      workspace ( m_workspace ) ,   // workspace
      s_APRECISION               ,   // absolute precision
      s_RPRECISION               ,   // relative precision
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
{ return Ostap::Utils::hash_combiner ( m_psx.tag () , m_psy.tag () , m_mmax ) ; }
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
double Ostap::Math::PS2DPol3Sym::evaluate 
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
  std::tie ( ierror , result , error ) = s_cubature.cubature
    ( tag () , &F , 20000 , s_APRECISION , s_RPRECISION , s_message , __FILE__ , __LINE__ ) ;
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
  typedef Ostap::Math::IntegrateY2<PS2DPol3Sym> IY ;
  static const Ostap::Math::GSL::Integrator1D<IY> s_integrator ;
  static const char message[] = "IntegrateY2(PS2DPol3Sym)" ;
  //
  const IY fy { this , x } ;
  const auto F = s_integrator.make_function ( &fy ) ;
  //  
  int    ierror    =  0   ;
  double result    =  1.0 ;
  double error     = -1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( Ostap::Utils::hash_combiner ( tag() , 'X' , x )                  , 
      &F                        ,   // the function
      y_low   , y_high          ,   // low & high edges
      workspace ( m_workspace ) ,   // workspace
      s_APRECISION               ,   // absolute precision
      s_RPRECISION               ,   // relative precision
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
{ return Ostap::Utils::hash_combiner ( m_ps.tag () , m_mmax ) ; }
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
  : Ostap::Math::PolyFactor2D
    ( Nx , 
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
  : Ostap::Math::PolyFactor2D 
    ( Nx , 
      Ny , 
      std::min ( xmin , xmax )   , std::max ( xmin , xmax ) ,
      std::max ( psy. lowEdge () , std::min ( ymin , ymax ) ) , 
      std::min ( psy.highEdge () , std::max ( ymin , ymax ) ) )
  , m_psy ( psy )
  , m_tau ( tau ) 
{}
// ============================================================================
// constructor from components
// ============================================================================
Ostap::Math::ExpoPS2DPol::ExpoPS2DPol
( const Ostap::Math::Positive2D&   pol , 
  const Ostap::Math::PhaseSpaceNL& psy ,        
  const double                     tau ) 
  : Ostap::Math::PolyFactor2D ( pol ) 
  , m_psy ( psy )
  , m_tau ( tau ) 
{
  Ostap::Assert ( m_psy.lowEdge   () < m_positive.ymax () , 
                  "Invalid setting of lowEdge/highEdge/ymin/ymax"   ,
                  "Ostap::Math::ExpoPS2DPol"   , 
    		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;  
  Ostap::Assert ( m_positive.ymin () < m_psy.highEdge() , 
                  "Invalid setting of lowEdge/highEdge/ymin/ymax"   ,
                  "Ostap::Math::ExpoPS2DPol"   , 
  		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;  
}
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
double Ostap::Math::ExpoPS2DPol::evaluate 
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
{ return Ostap::Utils::hash_combiner ( m_positive.tag() , m_psy.tag () , m_tau ) ; }
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
  : Ostap::Math::PolyFactor2D
    ( Nx , Ny , 
      std::min ( xmin , xmax ) , std::max ( xmin , xmax ) ,
      std::min ( ymin , ymax ) , std::max ( ymin , ymax ) ) 
  , m_tauX ( taux ) 
  , m_tauY ( tauy )
{}
// ===========================================================================
/// constructor from polynomial 
// ===========================================================================
Ostap::Math::Expo2DPol::Expo2DPol
( const Ostap::Math::Positive2D& pol  , 
  const double                   taux ,
  const double                   tauy ) 
  : Ostap::Math::PolyFactor2D ( pol ) 
  , m_tauX     ( taux ) 
  , m_tauY     ( tauy )
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
double Ostap::Math::Expo2DPol::evaluate  
( const double x ,
  const double y ) const 
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
{ return Ostap::Utils::hash_combiner ( m_positive.tag () , m_tauX , m_tauY ) ; }
// ============================================================================

// ===========================================================================
// constructor from the order
// ===========================================================================
Ostap::Math::Expo2DPolSym::Expo2DPolSym
( const double                       xmin , 
  const double                       xmax , 
  const unsigned short               N    ,
  const double                       tau  )      
  : Ostap::Math::PolyFactor2DSym
    ( N , std::min ( xmin , xmax ) , std::max ( xmin , xmax ) )
  , m_tau ( tau ) 
{}
// ===========================================================================
/// constructor from polynomial 
// ===========================================================================
Ostap::Math::Expo2DPolSym::Expo2DPolSym
( const Ostap::Math::Positive2DSym& pol , 
  const double                      tau ) 
  : Ostap::Math::PolyFactor2DSym ( pol ) 
  , m_tau      ( tau )
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
double Ostap::Math::Expo2DPolSym::evaluate 
( const double x ,
  const double y ) const 
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
{ return Ostap::Utils::hash_combiner ( m_positive.tag () , m_tau ) ; }
// ============================================================================




// ============================================================================
/// constructor 
// ============================================================================
Ostap::Math::Gauss2D::Gauss2D 
( const double  muX    , 
  const double  muY    , 
  const double  sigmaX , 
  const double  sigmaY ,
  const double  theta  ) 
  : m_muX       ( muX )
  , m_muY       ( muY )
  , m_sigmaX    ( std::abs ( sigmaX ) )
  , m_sigmaY    ( std::abs ( sigmaY ) )
  , m_theta     ( theta ) 
  , m_sintheta  ( std::sin ( theta ) ) 
  , m_costheta  ( std::cos ( theta ) )
{}
// ============================================================================
// set mux-parameter
// ============================================================================
bool Ostap::Math::Gauss2D::setMuX ( const double value )
{
  if ( s_equal ( m_muX , value ) ) { return false ; }
  m_muX = value ;
  return true ;
}
// ============================================================================
// set muy-parameter
// ============================================================================
bool Ostap::Math::Gauss2D::setMuY ( const double value )
{
  if ( s_equal ( m_muY , value ) ) { return false ; }
  m_muY = value ;
  return true ;
}
// ============================================================================
// set theta-parameter
// ============================================================================
bool Ostap::Math::Gauss2D::setTheta ( const double value )
{
  if ( s_equal ( m_theta , value ) ) { return false ; }
  m_theta    = value ;
  m_sintheta = std::sin ( m_theta ) ;
  m_costheta = std::cos ( m_theta ) ;
  return true ;
}
// ============================================================================
// set sigmax-parameter
// ============================================================================
bool Ostap::Math::Gauss2D::setSigmaX ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_sigmaX , avalue ) ) { return false ; }
  m_sigmaX = avalue ;
  return true ;
}
// ============================================================================
// set sigmay-parameter
// ============================================================================
bool Ostap::Math::Gauss2D::setSigmaY ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_sigmaY , avalue ) ) { return false ; }
  m_sigmaY = avalue ;
  return true ;
}
// ============================================================================
// get the value
// ============================================================================
double Ostap::Math::Gauss2D::operator () 
  ( const double x , const double y ) const 
{
  const double dx = x - m_muX ;
  const double dy = y - m_muY ;
  //
  const double ct = std::cos ( m_theta ) ;
  const double st = std::sin ( m_theta ) ;
  //
  const double dxp = ( ct * dx + st * dy ) / m_sigmaX ;
  const double dyp = ( ct * dy - st * dx ) / m_sigmaY ;
  //
  return 
    std::exp ( -0.5 * ( dxp * dxp + dyp * dyp ) ) / ( 2 * M_PI * m_sigmaX * m_sigmaY ) ;
}
// ============================================================================
// get the integral over the whole 2D-region
// ============================================================================
double Ostap::Math::Gauss2D::integral () const { return 1 ; }
// ============================================================================
/* get the integral over 2D-region
 *  \f[ \int_{x_low}^{x_high}\int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f]
 *  @param xlow  low  edge in x
 *  @param xhigh high edge in x
 *  @param ylow  low  edge in y
 *  @param yhigh high edge in y
 */
// ============================================================================
double Ostap::Math::Gauss2D::integral 
( const double xlow  ,
  const double xhigh ,
  const double ylow  , 
  const double yhigh ) const 
{ 
  if      ( s_equal ( xlow , xhigh ) || s_equal ( ylow  , yhigh ) ) { return 0 ; }
  //
  else if ( xlow  > xhigh ) { return -1*integral ( xhigh , xlow  , ylow  , yhigh ) ; }
  else if ( ylow  > yhigh ) { return -1*integral ( xlow  , xhigh , yhigh , ylow  ) ; }
  //
  if ( s_zero  ( m_sintheta ) || ( s_equal ( m_sigmaX , m_sigmaY ) ) ) 
  {
    return 
      ( Ostap::Math::gauss_cdf ( xhigh , m_muX , m_sigmaX ) - 
        Ostap::Math::gauss_cdf ( xlow  , m_muX , m_sigmaX ) ) * 
      ( Ostap::Math::gauss_cdf ( yhigh , m_muY , m_sigmaY ) - 
        Ostap::Math::gauss_cdf ( ylow  , m_muY , m_sigmaY ) ) ;
  }
  //
  // very far from the peak?
  // 
  const double sx = std::max ( std::abs ( m_costheta ) *  m_sigmaX , 
                               std::abs ( m_sintheta ) *  m_sigmaY ) ;
  const double sy = std::max ( std::abs ( m_costheta ) *  m_sigmaY , 
                               std::abs ( m_sintheta ) *  m_sigmaX ) ;
  
  if      ( xhigh <= m_muX - 50 * sx ) { return 0 ; }
  else if ( xlow  >= m_muX + 50 * sx ) { return 0 ; }
  else if ( yhigh <= m_muY - 50 * sy ) { return 0 ; }
  else if ( ylow  >= m_muY + 50 * sy ) { return 0 ; }
  //
  // split into smaller regions 
  //
  for ( SPLITS::const_iterator ix = s_SPLITS.begin() ; s_SPLITS.end() != ix  ; ++ix ) 
  {
    const double px = m_muX + (*ix) * sx ;
    if ( xlow < px && px < xhigh ) 
    {
      return 
        integral ( xlow , px    , ylow , yhigh ) +
        integral ( px   , xhigh , ylow , yhigh ) ;
    }
  }
  //
  for ( SPLITS::const_iterator iy = s_SPLITS.begin() ; s_SPLITS.end() != iy  ; ++iy ) 
  {
    const double py = m_muY + (*iy) * sy ;
    if ( ylow < py && py < yhigh ) 
    {
      return 
        integral ( xlow , xhigh , ylow , py    ) +
        integral ( xlow , xhigh , py   , yhigh ) ;
    }
  }
  //
  const bool in_tail = 
    ( xhigh <= m_muX + sx * s_SPLITS.front () ) || 
    ( xlow  >= m_muX + sx * s_SPLITS.back  () ) || 
    ( yhigh <= m_muY + sy * s_SPLITS.front () ) || 
    ( ylow  >= m_muY + sy * s_SPLITS.back  () ) ;
  // 
  // use 2D-cubature   
  //
  static const Ostap::Math::GSL::Integrator2D<Gauss2D> s_cubature{} ;
  static const char s_message[] = "Integral(Gauss2D)" ;
  const auto F = s_cubature.make_function ( this , xlow , xhigh , ylow , yhigh ) ;
  //
  int    ierror  =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  std::tie ( ierror , result , error ) = s_cubature.cubature
    ( tag () , &F , 20000 , 
      in_tail ? s_APRECISION_TAIL : s_APRECISION ,
      in_tail ? s_APRECISION_TAIL : s_RPRECISION , 
      s_message , __FILE__ , __LINE__ ) ;
  return  result ;
}
// ======================================================================
/*  integral over x-dimension
 *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}y\f]
 *  @param y     variable
 *  @param xlow  low  edge in y
 *  @param xhigh high edge in y
 */
// ============================================================================
double Ostap::Math::Gauss2D::integrateX 
( const double y     ,
  const double xlow  , 
  const double xhigh ) const 
{
  //
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; }
  else if ( xlow  > xhigh ) { return -1 * integrateX ( y , xhigh , xlow ) ; }
  //
  if ( s_zero  ( m_sintheta ) || ( s_equal ( m_sigmaX , m_sigmaY ) ) ) 
  {
    return 
      ( Ostap::Math::gauss_cdf ( xhigh , m_muX , m_sigmaX ) - 
        Ostap::Math::gauss_cdf ( xlow  , m_muX , m_sigmaX ) ) 
      * Ostap::Math::gauss_pdf ( y     , m_muY , m_sigmaY ) ;                             
  }
  // very far from the peak?
  //
  const double sx = std::max ( std::abs ( m_costheta ) * m_sigmaX , std::abs ( m_sintheta ) * m_sigmaY ) ;
  const double sy = std::max ( std::abs ( m_costheta ) * m_sigmaY , std::abs ( m_sintheta ) * m_sigmaX ) ;  
  //
  if      ( y     <= m_muY - 50 * sy ) { return 0 ; }
  else if ( y     >= m_muY + 50 * sy ) { return 0 ; }
  else if ( xhigh <= m_muX - 50 * sx ) { return 0 ; }
  else if ( xlow  >= m_muX + 50 * sx ) { return 0 ; }
  //
  //
  // split into smaller regions 
  //
  if ( xlow < m_muX + sx * s_SPLITS.back() || xhigh > m_muX + sx * s_SPLITS.front() ) 
  {
    for ( SPLITS::const_iterator ix = s_SPLITS.begin() ; s_SPLITS.end() != ix  ; ++ix ) 
    {
      const double px = m_muX + (*ix) * sx ;
      if ( xlow < px && px < xhigh ) 
      { return integrateX ( y , xlow , px    ) + integrateX ( y , px   , xhigh ) ; }
    }
  }
  //
  const bool in_tail = 
    ( xhigh <= m_muX + sx * s_SPLITS.front () ) || 
    ( xlow  >= m_muX + sx * s_SPLITS.back  () ) || 
    ( y     <= m_muY + sy * s_SPLITS.front () ) || 
    ( y     >= m_muY + sy * s_SPLITS.back  () ) ;
  // 
  typedef Ostap::Math::IntegrateX2<Gauss2D> IX ;
  static const Ostap::Math::GSL::Integrator1D<IX> s_integrator ;
  static const char message[] = "IntegrateX2(Gauss2D)" ;
  //
  const IX fx { this , y } ;
  const auto F = s_integrator.make_function ( &fx ) ;
  //  
  int    ierror    =  0   ;
  double result    =  1.0 ;
  double error     = -1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( Ostap::Utils::hash_combiner ( tag() , 'Y' , y ) , 
      &F                        ,   // the function
      xlow    , xhigh           ,   // low & high edges
      workspace ( m_workspace ) ,   // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size()        ,   // maximum number of subintervals
      message                   ,   // message 
      __FILE__  , __LINE__      ) ; // filename & line number 
  //
  return result ;
}
// ======================================================================
/* integral over x-dimension
 *  \f[ \int_{x_low}^{x_high} \mathcal{B}(x,y) \mathrm{d}x\f]
 *  @param x     variable
 *  @param ylow  low  edge in x
 *  @param yhigh high edge in x
 */
// ============================================================================
double Ostap::Math::Gauss2D::integrateY 
( const double x     ,
  const double ylow  ,
  const double yhigh ) const 
{
  //
  if      ( s_equal ( ylow , yhigh ) ) { return 0 ; }
  else if ( ylow  > yhigh ) { return -1 * integrateY ( x , yhigh , ylow ) ; }
  //
  if ( s_zero  ( m_sintheta ) || ( s_equal ( m_sigmaX , m_sigmaY ) ) ) 
  {
    return 
      ( Ostap::Math::gauss_cdf ( yhigh , m_muY , m_sigmaY ) - 
        Ostap::Math::gauss_cdf ( ylow  , m_muY , m_sigmaY ) )
      * Ostap::Math::gauss_pdf ( x     , m_muX , m_sigmaX ) ;                             
  }
  // very far from the peak?
  // 
  const double sx = std::max ( std::abs ( m_costheta ) * m_sigmaX , 
                               std::abs ( m_sintheta ) * m_sigmaY ) ;
  const double sy = std::max ( std::abs ( m_costheta ) * m_sigmaY , 
                               std::abs ( m_sintheta ) * m_sigmaX ) ;  
  //
  if      ( x     <= m_muX - 50 * sx ) { return 0 ; }
  else if ( x     >= m_muX + 50 * sx  ) { return 0 ; }
  else if ( yhigh <= m_muY - 50 * sy  ) { return 0 ; }
  else if ( ylow  >= m_muY + 50 * sy  ) { return 0 ; }
  //
  // split into smaller regions 
  //
  if ( ylow < m_muY + sy * s_SPLITS.back() || yhigh > m_muY + sy * s_SPLITS.front() ) 
  {
    for ( SPLITS::const_iterator iy = s_SPLITS.begin() ; s_SPLITS.end() != iy  ; ++iy ) 
    {
      const double py = m_muY + (*iy) * sy  ;
      if ( ylow < py && py < yhigh ) 
      { return integrateY ( x , ylow , py ) + integrateY ( x , py , yhigh ) ; }
    }
  }
  //
  const bool in_tail = 
    ( yhigh <= m_muY + sy * s_SPLITS.front () ) || 
    ( ylow  >= m_muY + sy * s_SPLITS.back  () ) || 
    ( x     <= m_muX + sx * s_SPLITS.front () ) || 
    ( x     >= m_muX + sx * s_SPLITS.back  () ) ;
  // 
  typedef Ostap::Math::IntegrateY2<Gauss2D> IY ;
  static const Ostap::Math::GSL::Integrator1D<IY> s_integrator ;
  static const char message[] = "IntegrateY2(Gauss2D)" ;
  //
  const IY fy { this , x } ;
  const auto F = s_integrator.make_function ( &fy ) ;
  //  
  int    ierror    =  0   ;
  double result    =  1.0 ;
  double error     = -1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( Ostap::Utils::hash_combiner ( tag() , 'X' , x ) , 
      &F                        ,   // the function
      ylow    , yhigh           ,   // low & high edges
      workspace ( m_workspace ) ,   // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size()        ,   // maximum number of subintervals
      message                   ,   // message 
      __FILE__  , __LINE__      ) ; // filename & line number 
  //
  return result ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Gauss2D::tag () const 
{ return Ostap::Utils::hash_combiner ( m_muX    , 
                             m_muY    , 
                             m_sigmaX , 
                             m_sigmaY , 
                             m_theta  ) ; }
// ============================================================================


// ============================================================================
/*  constructor from all parameters 
 *  @param mass particle mass (needed to calculate transverse mass)
 *  @param q    q-parameter of Tsallis  (q=1 corresponds to Boltzman statistics)
 *  @param T    the temperature
 *  @param mu   chemical potential 
 */
// ============================================================================
Ostap::Math::Tsallis2::Tsallis2 
( const double mass ,   // mass 
  const double T    ,   // temperature 
  const double q    ,   // q=1 -> Boltzman statistics 
  const double mu   )   // chemical potential 
  : m_mass ( std::abs ( mass ) ) 
  , m_T    ( std::abs ( T ) ) 
  , m_q    ( std::abs ( q ) ) 
  , m_mu   ( mu ) 
  , m_workspace () 
{}
// ============================================================================
/*  evaluate Tsallis function
 *  @param pt transverse momentum of the particle 
 *  @param y  rapidity of the particle 
 */
// ============================================================================
double Ostap::Math::Tsallis2::evaluate 
( const double pt ,
  const double y  ) const 
{
  if ( pt <= 0 ) { return 0 ; }
  //
  const double mtcy = mT ( pt ) * std::cosh ( y ) ;
  //
  const double f    = pt * mtcy ;
  const double arg  = ( mtcy - m_mu ) / m_T ;
  ///
  const double texp = Ostap::Math::tsallis_qexp ( -arg , m_q ) ;
  ///
  return texp <= 0 ? 0 : f * std::pow ( texp , m_q ) ;
}
// ============================================================================
// set new value for mass
// ============================================================================
bool Ostap::Math::Tsallis2::setMass ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_mass , avalue ) ) { return false ; }
  m_mass = avalue ;
  return true ;
}
// ============================================================================
// set new value for n-parameter
// ============================================================================
bool Ostap::Math::Tsallis2::setQ ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_q , avalue ) ) { return false ; }
  m_q = avalue ;
  return true ;
}
// ============================================================================
// set new value for T-parameter
// ============================================================================
bool Ostap::Math::Tsallis2::setT ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_T , avalue ) ) { return false ; }
  m_T = avalue ;
  return true ;
}
// ============================================================================
// set new value for mu-parameter
// ============================================================================
bool Ostap::Math::Tsallis2::setMu ( const double value )
{
  if ( s_equal ( m_mu , value ) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ============================================================================
//  get Tsallis integrals  
// ============================================================================
double Ostap::Math::Tsallis2::integral 
( const double ptlow  ,
  const double pthigh , 
  const double ylow   , 
  const double yhigh  ) const 
{
  if      ( s_equal ( ptlow  , pthigh ) ) { return 0 ; }
  else if ( s_equal ( ylow   , yhigh  ) ) { return 0 ; }
  else if (           pthigh < ptlow    ) { return - integral ( pthigh , ptlow  , ylow  , yhigh ) ; }
  else if (           yhigh  < ylow     ) { return - integral ( ptlow  , pthigh , yhigh , ylow  ) ; }
  else if ( pthigh <= 0                 ) { return 0 ; }
  //
  const double pt_min = std::max ( ptlow , 0.0 ) ;
  const double pt_max = pthigh ;
  const double y_min  = ylow   ;
  const double y_max  = yhigh  ;
  //
  // use cubature   
  static const Ostap::Math::GSL::Integrator2D<Tsallis2> s_cubature{} ;
  static const char s_message[] = "Integral(Tsallis2)" ;
  const auto F = s_cubature.make_function ( this , pt_min , pt_max , y_min , y_max ) ;
  //
  int    ierror  =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  std::tie ( ierror , result , error ) = s_cubature.cubature
    ( tag () , &F , 20000 , s_APRECISION , s_RPRECISION , s_message , __FILE__ , __LINE__ ) ;
  // ==========================================================================
  return result ;
}
// ============================================================================
//  get the integral between ylow-yhigh for given pt 
// ============================================================================
double Ostap::Math::Tsallis2::integrate_y  
( const double pt     ,
  const double ylow   , 
  const double yhigh  ) const 
{
  //
  if      ( pt <= 0                  ) { return 0 ; }
  else if ( s_equal ( ylow , yhigh ) ) { return 0 ; }
  else if (           yhigh < ylow   ) { return -integrate_y ( pt , yhigh , ylow ) ; }
  //
  typedef Ostap::Math::IntegrateY2<Tsallis2> IY ;
  static const Ostap::Math::GSL::Integrator1D<IY> s_integrator ;
  static const char message[] = "IntegrateY(Tsallis2)" ;
  //
  const IY fy { this , pt } ;
  const auto F = s_integrator.make_function ( &fy ) ;
  //  
  int    ierror    =  0   ;
  double result    =  1.0 ;
  double error     = -1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( Ostap::Utils::hash_combiner ( tag() , 'P' , pt )                   , 
      &F                        ,   // the function
      ylow    , yhigh           ,   // low & high edges
      workspace ( m_workspace ) ,   // workspace
      s_APRECISION              ,   // absolute precision
      s_RPRECISION              ,   // relative precision
      m_workspace.size()        ,   // maximum number of subintervals
      message                   ,   // message 
      __FILE__  , __LINE__      ) ; // filename & line number 
  //
  return result ;
}
// ============================================================================
// get the integral between ptlow-pthigh for given rapidity 
// ============================================================================
double Ostap::Math::Tsallis2::integrate_pt
( const double y      , 
  const double ptlow  ,
  const double pthigh ) const 
{
  //
  if      ( s_equal ( ptlow  , pthigh ) ) { return 0 ; }
  else if (           pthigh < ptlow    ) { return -integrate_pt ( y , pthigh , ptlow ) ; }
  else if ( pthigh <= 0                 ) { return 0 ; }
  //
  const double pt_min = std::max ( ptlow , 0.0 ) ;
  const double pt_max = pthigh ;
  //
  typedef Ostap::Math::IntegrateX2<Tsallis2> IPT ;
  static const Ostap::Math::GSL::Integrator1D<IPT> s_integrator ;
  static const char message[] = "IntegratePT(Tsallis2)" ;
  //
  const IPT fpt { this , y } ;
  const auto F = s_integrator.make_function ( &fpt ) ;
  //  
  int    ierror    =  0   ;
  double result    =  1.0 ;
  double error     = -1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( Ostap::Utils::hash_combiner ( tag() , 'Y' , y )                   , 
      &F                        ,   // the function
      pt_min , pt_max           ,   // low & high edges
      workspace ( m_workspace ) ,   // workspace
      s_APRECISION              ,   // absolute precision
      s_RPRECISION              ,   // relative precision
      m_workspace.size()        ,   // maximum number of subintervals
      message                   ,   // message 
      __FILE__  , __LINE__      ) ; // filename & line number 
  //
  return result ; 
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Tsallis2::tag () const 
{ 
  static const std::string s_name = "Tsallis2" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mass , m_q , m_T , m_mu ) ; 
}
// ============================================================================


// ============================================================================
//                                                                      The END 
// ============================================================================
