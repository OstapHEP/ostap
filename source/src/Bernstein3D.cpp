// ============================================================================
// Include files
// ============================================================================
// STD& STL
// ============================================================================
#include <cassert>
#include <numeric>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Hash.h"
#include "Ostap/Math.h"
#include "Ostap/Bernstein2D.h"
#include "Ostap/Bernstein3D.h"
#include "Ostap/StatusCode.h"
// ============================================================================
// Local
// ============================================================================
#include "local_math.h"
#include "local_hash.h"
// ============================================================================
/** @file
 *  Implementation file for functions, related to Bernstein's polynomnials
 *
 *  @see http://en.wikipedia.org/wiki/Bernstein_polynomial
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2010-04-19
 */
// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::Bernstein3D::Bernstein3D 
( const unsigned short       nX    ,
  const unsigned short       nY    ,
  const unsigned short       nZ    ,
  const double               xmin  ,
  const double               xmax  ,
  const double               ymin  ,
  const double               ymax  ,
  const double               zmin  ,
  const double               zmax  ) 
  : Ostap::Math::Parameters ( ( nX + 1 ) * ( nY + 1 ) * ( nZ + 1 ) )
  , m_nx   ( nX )
  , m_ny   ( nY )
  , m_nz   ( nZ )
    //
  , m_xmin ( std::min ( xmin , xmax ) )
  , m_xmax ( std::max ( xmin , xmax ) )
  , m_ymin ( std::min ( ymin , ymax ) )
  , m_ymax ( std::max ( ymin , ymax ) )
  , m_zmin ( std::min ( zmin , zmax ) )
  , m_zmax ( std::max ( zmin , zmax ) )
    //
  , m_bx   ()
  , m_by   ()
  , m_bz   ()
    //
  , m_fx   ( nX + 1 , 0.0 )
  , m_fy   ( nY + 1 , 0.0 )
  , m_fz   ( nZ + 1 , 0.0 )
{
  //
  m_bx.reserve ( m_nx + 1 ) ;
  m_by.reserve ( m_ny + 1 ) ;
  m_bz.reserve ( m_nz + 1 ) ;
  //
  typedef  Ostap::Math::Bernstein::Basic BB ;
  for ( unsigned short ix = 0 ; ix <= nX ; ++ix )
  { m_bx.push_back ( Bernstein ( BB ( ix , nX ) , xmin , xmax ) ) ; }
  //
  for ( unsigned short iy = 0 ; iy <= nY ; ++iy )
  { m_by.push_back ( Bernstein ( BB ( iy , nY ) , ymin , ymax ) ) ; }
  //
  for ( unsigned short iz = 0 ; iz <= nZ ; ++iz )
  { m_bz.push_back ( Bernstein ( BB ( iz , nZ ) , zmin , zmax ) ) ; }
  //
}
// ============================================================================
// constructor from parameters 
// ============================================================================
Ostap::Math::Bernstein3D::Bernstein3D 
( const std::vector<double>& pars  ,
  const unsigned short       nX    ,
  const unsigned short       nY    ,
  const unsigned short       nZ    ,
  const double               xmin  ,
  const double               xmax  ,
  const double               ymin  ,
  const double               ymax  ,
  const double               zmin  ,
  const double               zmax  )
  : Ostap::Math::Bernstein3D ( nX   , nY   , nZ ,
                               xmin , xmax ,
                               ymin , ymax ,
                               zmin , zmax )
{
  setPars ( pars ) ;
}
// ======================================================================
/*  As a product of three 1D-polynomials:
 *  \f[  B_{n^x,n^y,n^z}(x,y,z) \equiv 
 *      B^{n^x}(x)B^{n^y}(y)B^{n^z}(z) = 
 *  \left(\sum_{i=0}{n^{x}} \alpha_{i} B_{n^{x}}^i(x)]\right)
 *  \left(\sum_{j=0}{n^{y}} \beta_{j} B_{n^{y}}^j(y)]\right) = 
 *  \left(\sum_{k=0}{n^{z}} \gamma_{k} B_{n^{z}}^k(z)]\right) = 
 *    \sum_{i=0}{n^{x}}
 *    \sum_{j=0}{n^{y}} 
 *    \sum_{k=0}{n^{z}} 
 *   \alpha_{i}\beta_{j}\gamma_{k} 
 *    B_{n^{x}}^i(x) B_{n^{y}}^j(y) B_{n^{z}}^k(z) \f]
 */          
// ======================================================================
Ostap::Math::Bernstein3D::Bernstein3D
( const Ostap::Math::Bernstein& bx , 
  const Ostap::Math::Bernstein& by ,
  const Ostap::Math::Bernstein& bz ) 
  : Ostap::Math::Parameters ( ( bx.n  () + 1 ) * ( by.n () + 1 ) * ( bz.n () + 1 ) )
  , m_nx   ( bx. n   () )
  , m_ny   ( by. n   () )
  , m_nz   ( bz. n   () )
  , m_xmin ( bx.xmin () )
  , m_xmax ( bx.xmax () )
  , m_ymin ( by.xmin () )
  , m_ymax ( by.xmax () )
  , m_zmin ( bz.xmin () )
  , m_zmax ( bz.xmax () )
    //
  , m_bx   () 
  , m_by   ()
  , m_bz   ()
    // 
  , m_fx   ( bx.n() + 1 , 0.0 )
  , m_fy   ( by.n() + 1 , 0.0 )
  , m_fz   ( bz.n() + 1 , 0.0 )
{
  //
  m_bx.reserve ( m_nx + 1 ) ;
  m_by.reserve ( m_ny + 1 ) ;
  m_bz.reserve ( m_nz + 1 ) ;
  //
  typedef  Ostap::Math::Bernstein::Basic BB ;
  for ( unsigned short ix = 0 ; ix <= m_nx ; ++ix )
  { m_bx.push_back ( Bernstein ( BB ( ix , m_nx ) , m_xmin , m_xmax ) ) ; }
  //
  for ( unsigned short iy = 0 ; iy <= m_ny ; ++iy )
  { m_by.push_back ( Bernstein ( BB ( iy , m_ny ) , m_ymin , m_ymax ) ) ; }
  //
  for ( unsigned short iz = 0 ; iz <= m_nz ; ++iz )
  { m_bz.push_back ( Bernstein ( BB ( iz , m_nz ) , m_zmin , m_zmax ) ) ; }
  //
  for ( unsigned short ix = 0 ; ix <= m_nx ; ++ix ) 
  { for ( unsigned short iy = 0 ; iy <= m_ny ; ++iy )
    { for ( unsigned short iz = 0 ; iz <= m_nz ; ++iz ) 
      { setPar ( ix , iy , iz , bx.par ( ix ) * by.par ( iy ) * bz.par ( iz ) ) ; } } }
  //
}
// ============================================================================
// constructor from symmetric polynomial
// ============================================================================
Ostap::Math::Bernstein3D::Bernstein3D 
( const Ostap::Math::Bernstein3DSym& right ) 
  : Ostap::Math::Parameters ( ( right.nX() + 1 ) * ( right.nY() + 1 ) * ( right.nZ () + 1 ) )
  , m_nx   ( right.nX() )
  , m_ny   ( right.nY() )
  , m_nz   ( right.nZ() )
    //
  , m_xmin ( right.xmin() )
  , m_xmax ( right.xmax() )
  , m_ymin ( right.ymin() )
  , m_ymax ( right.ymax() )
  , m_zmin ( right.zmin() )
  , m_zmax ( right.zmax() )
    //
  , m_bx   ()
  , m_by   ()
  , m_bz   ()
    //
  , m_fx   ( right.nX() + 1 , 0.0 )
  , m_fy   ( right.nY() + 1 , 0.0 )
  , m_fz   ( right.nZ() + 1 , 0.0 )
{
  //
  m_bx.reserve ( m_nx + 1 ) ;
  m_by.reserve ( m_ny + 1 ) ;
  m_bz.reserve ( m_nz + 1 ) ;
  //
  typedef  Ostap::Math::Bernstein::Basic BB ;
  for ( unsigned short ix = 0 ; ix <= m_nx ; ++ix )
  { m_bx.push_back ( Bernstein ( BB ( ix , m_nx ) , m_xmin , m_xmax ) ) ; }
  //
  for ( unsigned short iy = 0 ; iy <= m_ny ; ++iy )
  { m_by.push_back ( Bernstein ( BB ( iy , m_ny ) , m_ymin , m_ymax ) ) ; }
  //
  for ( unsigned short iz = 0 ; iz <= m_nz ; ++iz )
  { m_bz.push_back ( Bernstein ( BB ( iz , m_nz ) , m_zmin , m_zmax ) ) ; }
  //
}
// ============================================================================
// constructor from mixed symmetry polynomial
// ============================================================================
Ostap::Math::Bernstein3D::Bernstein3D 
( const Ostap::Math::Bernstein3DMix& right ) 
  : Ostap::Math::Parameters ( ( right.nX() + 1 ) * ( right.nY() + 1 ) * ( right.nZ () + 1 ) )
  , m_nx   ( right.nX() )
  , m_ny   ( right.nY() )
  , m_nz   ( right.nZ() )
    //
  , m_xmin ( right.xmin() )
  , m_xmax ( right.xmax() )
  , m_ymin ( right.ymin() )
  , m_ymax ( right.ymax() )
  , m_zmin ( right.zmin() )
  , m_zmax ( right.zmax() )
    //
  , m_bx   ()
  , m_by   ()
  , m_bz   ()
    //
  , m_fx   ( right.nX() + 1 , 0.0 )
  , m_fy   ( right.nY() + 1 , 0.0 )
  , m_fz   ( right.nZ() + 1 , 0.0 )
{
  //
  m_bx.reserve ( m_nx + 1 ) ;
  m_by.reserve ( m_ny + 1 ) ;
  m_bz.reserve ( m_nz + 1 ) ;
  //
  typedef  Ostap::Math::Bernstein::Basic BB ;
  for ( unsigned short ix = 0 ; ix <= m_nx ; ++ix )
  { m_bx.push_back ( Bernstein ( BB ( ix , m_nx ) , m_xmin , m_xmax ) ) ; }
  //
  for ( unsigned short iy = 0 ; iy <= m_ny ; ++iy )
  { m_by.push_back ( Bernstein ( BB ( iy , m_ny ) , m_ymin , m_ymax ) ) ; }
  //
  for ( unsigned short iz = 0 ; iz <= m_nz ; ++iz )
  { m_bz.push_back ( Bernstein ( BB ( iz , m_nz ) , m_zmin , m_zmax ) ) ; }
  //
}
// ============================================================================
// swap  two 3D-polynomials 
// ============================================================================
void Ostap::Math::Bernstein3D::swap ( Ostap::Math::Bernstein3D&  right ) 
{
  Ostap::Math::Parameters::swap ( right ) ;
  //
  std::swap ( m_nx   , right.m_nx    ) ;
  std::swap ( m_ny   , right.m_ny    ) ;
  std::swap ( m_nz   , right.m_nz    ) ;
  std::swap ( m_xmin , right.m_xmin  ) ;
  std::swap ( m_xmax , right.m_xmax  ) ;
  std::swap ( m_ymin , right.m_ymin  ) ;
  std::swap ( m_ymax , right.m_ymax  ) ;
  std::swap ( m_zmin , right.m_zmin  ) ;
  std::swap ( m_zmax , right.m_zmax  ) ;
  std::swap ( m_bx   , right.m_bx    ) ;
  std::swap ( m_by   , right.m_by    ) ;
  std::swap ( m_bz   , right.m_bz    ) ;
  std::swap ( m_fz   , right.m_fx    ) ;
  std::swap ( m_fz   , right.m_fx    ) ;
  std::swap ( m_fz   , right.m_fx    ) ;
}
// ============================================================================
// helper function to make calculations
// ============================================================================
double Ostap::Math::Bernstein3D::calculate
( const std::vector<double>& fx , 
  const std::vector<double>& fy , 
  const std::vector<double>& fz ) const 
{
  double       result = 0 ;
  for  ( unsigned short ix = 0 ; ix <= nX () ; ++ix )
  {
    for  ( unsigned short iy = 0 ; iy <= nY () ; ++iy )
    { 
      for  ( unsigned short iz = 0 ; iz <= nZ () ; ++iz )
      { 
        result += par ( ix , iy , iz ) * fx[ix] * fy[iy] * fz[iz]; 
      }
    }
  }
  //
  const double scalex = ( nX () + 1 ) / ( xmax() - xmin() ) ;
  const double scaley = ( nY () + 1 ) / ( ymax() - ymin() ) ;
  const double scalez = ( nZ () + 1 ) / ( zmax() - zmin() ) ;
  //
  return result * scalex * scaley * scalez ;
}
// ============================================================================
// get the value
// ============================================================================
double Ostap::Math::Bernstein3D::evaluate
( const double x ,
  const double y , 
  const double z ) const
{
  /// the trivial cases
  if ( x < xmin () || x > xmax () ) { return 0.0        ; }
  if ( y < ymin () || y > ymax () ) { return 0.0        ; }
  if ( z < zmin () || z > zmax () ) { return 0.0        ; }
  //
  if      ( 0 == npars ()       ) { return 0.0        ; }
  else if ( 1 == npars ()       ) 
  { 
    const double scalex = ( nX () + 1 ) / ( xmax() - xmin() ) ;
    const double scaley = ( nY () + 1 ) / ( ymax() - ymin() ) ;
    const double scalez = ( nZ () + 1 ) / ( zmax() - zmin() ) ;
    //
    return m_pars [0] * scalex * scaley * scalez ;
  }
  ///
  // std::vector<double> fx ( nX () + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nX ()  ; ++i )
  { m_fx[i] = m_bx[i] ( x )  ; }
  //
  // std::vector<double> fy ( nY ()  + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nY ()  ; ++i )
  { m_fy[i] = m_by[i] ( y )  ; }
  //
  // std::vector<double> fz ( nZ () + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nZ () ; ++i )
  { m_fz[i] = m_bz[i] ( z )  ; }
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
/** get the integral over 3D-region
 *  \f[  x_{min} < x < x_{max}, y_{min}< y< y_{max} , z_{min} < z < z_{max}\f]
 */
// ============================================================================
double Ostap::Math::Bernstein3D::integral() const
{ return std::accumulate ( m_pars.begin() , m_pars.end() , 0.0 ) ; }
// ============================================================================
/* get the integral over 3D-region
 *  \f[ \int_{x_{low}}^{x_{high}}
 *      \int_{y_{low}}^{y_{high}} 
 *      \int_{z_{low}}^{z_{high}} 
 *\mathcal{B}(x,y,z) \mathrm{d}x\mathrm{d}y\mathrm{d}z\f]
 *  @param xlow  low  edge in x
 *  @param xhigh high edge in x
 *  @param ylow  low  edge in y
 *  @param yhigh high edge in y
 *  @param zlow  low  edge in z
 *  @param zhigh high edge in z
 */
// ============================================================================
double Ostap::Math::Bernstein3D::integral
( const double xlow , const double xhigh ,
  const double ylow , const double yhigh ,
  const double zlow , const double zhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ||
            s_equal ( ylow , yhigh ) ||
            s_equal ( zlow , zhigh ) ) { return 0 ; }
  //
  else if ( s_equal ( xlow  , xmin () ) &&
            s_equal ( xhigh , xmax () ) &&
            s_equal ( ylow  , ymin () ) &&
            s_equal ( yhigh , ymax () ) &&  
            s_equal ( zlow  , zmin () ) &&
            s_equal ( zhigh , zmax () ) )  { return integral () ; }
  //
  else if ( xlow  > xhigh ) 
  { return -1*integral ( xhigh , xlow  , ylow  , yhigh , zlow  , zhigh ) ; }
  else if ( ylow  > yhigh ) 
  { return -1*integral ( xlow  , xhigh , yhigh , ylow  , zlow  , zhigh ) ; }
  else if ( zlow  > zhigh ) 
  { return -1*integral ( xlow  , xhigh , ylow  , yhigh , zhigh , zlow  ) ; }
  //
  else if ( xhigh <  xmin () || xlow >  xmax() ) { return 0 ; }
  else if ( yhigh <  ymin () || ylow >  ymax() ) { return 0 ; }
  else if ( zhigh <  zmin () || zlow >  zmax() ) { return 0 ; }
  //
  const double  x_low  = std::max ( xmin() , xlow  ) ;
  const double  x_high = std::min ( xmax() , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  const double  y_low  = std::max ( ymin() , ylow  ) ;
  const double  y_high = std::min ( ymax() , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  const double  z_low  = std::max ( zmin() , zlow  ) ;
  const double  z_high = std::min ( zmax() , zhigh ) ;
  if ( z_low >= z_high ) { return 0 ; }
  //
  // std::vector<double> fx ( nX () + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nX ()  ; ++i )
  { m_fx[i] = m_bx[i].integral ( x_low , x_high ) ; }
  //
  // std::vector<double> fy ( nY () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nY ()  ; ++i )
  { m_fy[i] = m_by[i].integral ( y_low , y_high ) ; }
  //
  // std::vector<double> fz ( nZ ()  + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nZ ()  ; ++i )
  { m_fz[i] = m_bz[i].integral ( z_low , z_high ) ; }
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
/*  integral over x-dimension
 *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y,z) \mathrm{d}y\f]
 *  @param y     variable
 *  @param z     variable
 *  @param ylow  low  edge in y
 *  @param yhigh high edge in y
 */
// ============================================================================
double Ostap::Math::Bernstein3D::integrateX
( const double y    ,
  const double z    ,
  const double xlow , const double xhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; }
  else if ( xlow  > xhigh  ) { return -1*integrateX ( y , z , xhigh , xlow ) ; }
  else if ( xhigh <= xmin () || xlow >= xmax() ) { return 0 ; }
  else if ( y     <  ymin () || y    >  ymax() ) { return 0 ; }
  else if ( z     <  zmin () || z    >  zmax() ) { return 0 ; }
  else if ( s_equal ( xlow  , xmin () ) &&
            s_equal ( xhigh , xmax () )         ) { return integrateX ( y ,  z ) ; }
  //
  const double  x_low  = std::max ( xmin() , xlow  ) ;
  const double  x_high = std::min ( xmax() , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  // std::vector<double> fx ( nX () + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nX () ; ++i )
  { m_fx[i] = m_bx[i].integral ( x_low , x_high ) ; }
  //
  // std::vector<double> fy ( nY () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nY () ; ++i )
  { m_fy[i] = m_by[i] ( y ) ; }
  //
  // std::vector<double> fz ( nZ ()  + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nZ ()  ; ++i )
  { m_fz[i] = m_bz[i] ( z ) ; }
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
/** integral over y-dimension
 *  \f[ \int_{x_low}^{x_high} \mathcal{B}(x,y,z) \mathrm{d}y\f]
 *  @param x     variable
 *  @param z     variable
 *  @param ylow  low  edge in y
 *  @param yhigh high edge in y
 */
// ============================================================================
double Ostap::Math::Bernstein3D::integrateY
( const double x    ,
  const double z    ,
  const double ylow , const double yhigh ) const
{
  if      ( s_equal ( ylow  , yhigh ) ) { return 0 ; }
  else if ( ylow  >  yhigh ) { return -1*integrateY ( x , z , yhigh , ylow  ) ; }
  else if ( x     <  xmin () || x    >  xmax() ) { return 0 ; }
  else if ( z     <  zmin () || z    >  zmax() ) { return 0 ; }
  else if ( yhigh <= ymin () || ylow >= ymax() ) { return 0 ; }
  else if ( s_equal ( ylow  , ymin () ) &&
            s_equal ( yhigh , ymax () )         ) { return integrateY ( x , z ) ; }
  //
  const double  y_low  = std::max ( ymin() , ylow  ) ;
  const double  y_high = std::min ( ymax() , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  // std::vector<double> fx ( nX ()  + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nX ()  ; ++i )
  { m_fx[i] = m_bx[i] ( x ) ; }
  //
  // std::vector<double> fy ( nY () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nY () ; ++i )
  { m_fy[i] = m_by[i].integral ( y_low , y_high ) ; }
  //
  // std::vector<double> fz ( nZ ()  + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nZ ()  ; ++i )
  { m_fz[i] = m_bz[i] ( z ) ; }
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
/** integral over z-dimension
 *  \f[ \int_{z_low}^{z_high} \mathcal{B}(x,y,z) \mathrm{d}z\f]
 *  @param x     variable
 *  @param y     variable
 *  @param zlow  low  edge in z
 *  @param zhigh high edge in z
 */
// ============================================================================
double
Ostap::Math::Bernstein3D::integrateZ
( const double x    ,
  const double y    ,
  const double zlow , const double zhigh ) const
{
  if      ( s_equal ( zlow  , zhigh ) ) { return 0 ; }
  else if ( zlow  >  zhigh ) { return -1*integrateZ ( x , y , zhigh , zlow  ) ; }
  else if ( x     <  xmin () || x    >  xmax() ) { return 0 ; }
  else if ( y     <  ymin () || y    >  ymax() ) { return 0 ; }
  else if ( zhigh <= zmin () || zlow >= zmax() ) { return 0 ; }
  else if ( s_equal ( zlow  , zmin () ) &&
            s_equal ( zhigh , zmax () )         ) { return integrateZ ( x , y ) ; }
  //
  const double  z_low  = std::max ( zmin() , zlow  ) ;
  const double  z_high = std::min ( zmax() , zhigh ) ;
  if ( z_low >= z_high ) { return 0 ; }
  //
  // std::vector<double> fx ( nX () + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nX ()  ; ++i )
  { m_fx[i] = m_bx[i] ( x ) ; }
  //
  // std::vector<double> fy ( nY () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nY ()  ; ++i )
  { m_fy[i] = m_by[i] ( y ) ; }
  //
  // std::vector<double> fz ( nZ () + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nZ () ; ++i ) 
  { m_fz[i] = m_bz[i].integral ( z_low , z_high ) ; }
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
double
Ostap::Math::Bernstein3D::integrateX
( const double y , 
  const double z ) const
{
  if      ( y < ymin () || y > ymax() ) { return 0 ; }
  else if ( z < zmin () || z > zmax() ) { return 0 ; }
  //
  // const std::vector<double> fx ( nX () + 1 , ( xmax() - xmin () ) / ( nX () + 1 ) ) ;
  std::fill ( m_fx.begin() , m_fx.end() , ( xmax() - xmin () ) / ( nX () + 1 ) ) ;
  //
  // std::vector<double> fy ( nY () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nY () ; ++i )
  { m_fy[i] = m_by[i] ( y ) ; }
  //
  // std::vector<double> fz ( nZ () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nZ ()  ; ++i )
  { m_fz[i] = m_bz[i] ( z ) ; }
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
double
Ostap::Math::Bernstein3D::integrateY
( const double x , 
  const double z ) const
{
  if      ( x < xmin () || x > xmax() ) { return 0 ; }
  else if ( z < zmin () || z > zmax() ) { return 0 ; }
  //
  // std::vector<double> fx ( nX () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nX () ; ++i )
  { m_fx[i] = m_bx[i] ( x ) ; }
  //
  // const std::vector<double> fy ( nY () + 1 , ( ymax() - ymin () ) / ( nY () + 1 ) ) ;
  std::fill  ( m_fy.begin () , m_fy.end() , ( ymax() - ymin () ) / ( nY () + 1 ) ) ;
  //
  // std::vector<double> fz ( nZ () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nZ ()  ; ++i )
  { m_fz[i] = m_bz[i] ( z ) ; }
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
double
Ostap::Math::Bernstein3D::integrateZ
( const double x , 
  const double y ) const
{
  if      ( x < xmin () || x > xmax() ) { return 0 ; }
  else if ( y < ymin () || y > ymax() ) { return 0 ; }
  //
  // std::vector<double> fx ( nX () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nX () ; ++i )
  { m_fx[i] = m_bx[i] ( x ) ; }
  //
  // std::vector<double> fy ( nY () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nY () ; ++i )
  { m_fy[i] = m_by[i] ( y ) ; }
  //
  // const std::vector<double> fz ( nZ() + 1 , ( zmax() - zmin () ) / ( nZ () + 1 ) ) ;
  std::fill  ( m_fz.begin () , m_fz.end() , ( zmax() - zmin () ) / ( nZ () + 1 ) ) ;
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
/*  integral over x&y-dimensions
 *  \f[ \int_{x_{low}}^{x_{high}}
 *      \int_{y_{low}}^{y_{high}} \mathcal{B}(x,y,z) \mathrm{d}x\mathrm{d}y\f]
 *  @param z     variable
 *  @param xlow  low  edge in x
 *  @param xhigh high edge in x
 *  @param ylow  low  edge in y
 *  @param yhigh high edge in y
 */
// ============================================================================
double Ostap::Math::Bernstein3D::integrateXY
( const double z    ,                          
  const double xlow , const double xhigh ,
  const double ylow , const double yhigh ) const 
{
  if      ( s_equal ( xlow  , xhigh ) ) { return 0 ; }
  else if ( s_equal ( ylow  , yhigh ) ) { return 0 ; }
  else if ( xlow  >  xhigh ) { return -1*integrateXY ( z , xhigh , xlow  , ylow  , yhigh ) ; }
  else if ( ylow  >  yhigh ) { return -1*integrateXY ( z , xlow  , xhigh , yhigh , ylow  ) ; }
  else if ( z     <  zmin () || z    >  zmax() ) { return 0 ; }
  else if ( s_equal ( xlow  , xmin () ) &&
            s_equal ( xhigh , xmax () ) &&
            s_equal ( ylow  , ymin () ) &&
            s_equal ( yhigh , ymax () ) ) { return integrateXY ( z ) ; }
  //  //
  const double  x_low  = std::max ( xmin() , xlow  ) ;
  const double  x_high = std::min ( xmax() , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  const double  y_low  = std::max ( ymin() , ylow  ) ;
  const double  y_high = std::min ( ymax() , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  // std::vector<double> fx ( nX () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nX ()  ; ++i )
  { m_fx[i] = m_bx[i].integral ( x_low , x_high ) ; }
  //
  // std::vector<double> fy ( nY () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nY () ; ++i )
  { m_fy[i] = m_by[i].integral ( y_low , y_high ) ; }
  //
  // std::vector<double> fz ( nZ () + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nZ () ; ++i )
  { m_fz[i] = m_bz[i] ( z ) ; }
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
/** integral over x&z-dimensions
 *  \f[ \int_{x_{low}}^{x_{high}}
 *      \int_{z_{low}}^{z_{high}} \mathcal{B}(x,y,z) \mathrm{d}x\mathrm{d}z\f]
 *  @param y     variable
 *  @param xlow  low  edge in x
 *  @param xhigh high edge in x
 *  @param zlow  low  edge in y
 *  @param zhigh high edge in y
 */
// ============================================================================
double Ostap::Math::Bernstein3D::integrateXZ
( const double y    ,                          
  const double xlow , const double xhigh ,
  const double zlow , const double zhigh ) const 
{
  if      ( s_equal ( xlow  , xhigh ) ) { return 0 ; }
  else if ( s_equal ( zlow  , zhigh ) ) { return 0 ; }
  else if ( xlow  >  xhigh ) { return -1*integrateXZ ( y , xhigh , xlow  , zlow  , zhigh ) ; }
  else if ( zlow  >  zhigh ) { return -1*integrateXZ ( y , xlow  , xhigh , zhigh , zlow  ) ; }
  else if ( y     <  ymin () || y    >  ymax() ) { return 0 ; }
  else if ( s_equal ( xlow  , xmin () ) &&
            s_equal ( xhigh , xmax () ) &&
            s_equal ( zlow  , zmin () ) &&
            s_equal ( zhigh , zmax () ) ) { return integrateXZ ( y ) ; }
  //
  const double  x_low  = std::max ( xmin() , xlow  ) ;
  const double  x_high = std::min ( xmax() , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  const double  z_low  = std::max ( zmin() , zlow  ) ;
  const double  z_high = std::min ( zmax() , zhigh ) ;
  if ( z_low >= z_high ) { return 0 ; }
  //
  // std::vector<double> fx ( nX () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nX () ; ++i )
  { m_fx[i] = m_bx[i].integral ( x_low , x_high ) ; }
  //
  // std::vector<double> fy ( nY () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nY () ; ++i )
  { m_fy[i] = m_by[i] ( y ) ; }
  //
  // std::vector<double> fz ( nZ () + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nZ () ; ++i )
  { m_fz[i] = m_bz[i].integral ( z_low , z_high ) ; }
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
/** integral over y&z-dimensions
 *  \f[ \int_{y_{low}}^{y_{high}}
 *      \int_{z_{low}}^{z_{high}} \mathcal{B}(x,y,z) \mathrm{d}y\mathrm{d}z\f]
 *  @param x     variable
 *  @param ylow  low  edge in x
 *  @param yhigh high edge in x
 *  @param zlow  low  edge in y
 *  @param zhigh high edge in y
 */
// ============================================================================ 
double Ostap::Math::Bernstein3D::integrateYZ
( const double x    ,                          
  const double ylow , const double yhigh ,
  const double zlow , const double zhigh ) const 
{
  if      ( s_equal ( ylow  , yhigh ) ) { return 0 ; }
  else if ( s_equal ( zlow  , zhigh ) ) { return 0 ; }
  else if ( ylow  >  yhigh ) { return -1*integrateYZ ( x , yhigh , ylow  , zlow  , zhigh ) ; }
  else if ( zlow  >  zhigh ) { return -1*integrateYZ ( x , ylow  , yhigh , zhigh , zlow  ) ; }
  else if ( x     <  xmin () || x    >  xmax() ) { return 0 ; }
  else if ( s_equal ( ylow  , ymin () ) &&
            s_equal ( yhigh , ymax () ) &&
            s_equal ( zlow  , zmin () ) &&
            s_equal ( zhigh , zmax () ) ) { return integrateYZ ( x ) ; }
  //
  const double  y_low  = std::max ( ymin() , ylow  ) ;
  const double  y_high = std::min ( ymax() , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  const double  z_low  = std::max ( zmin() , zlow  ) ;
  const double  z_high = std::min ( zmax() , zhigh ) ;
  if ( z_low >= z_high ) { return 0 ; }
  //
  // std::vector<double> fx ( nX () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nX () ; ++i )
  { m_fx[i] = m_bx[i] ( x ) ; }
  //
  // std::vector<double> fy ( nY () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nY () ; ++i )
  { m_fy[i] = m_by[i].integral ( y_low , y_high ) ; }
  //
  // std::vector<double> fz ( nZ () + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <=  nZ () ; ++i )
  { m_fz[i] = m_bz[i].integral ( z_low , z_high ) ; }
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
/*  integral over x&y-dimensions
 *  \f[ \int_{x_{min}}^{x_{max}}
 *      \int_{y_{min}}^{y_{max}} \mathcal{B}(x,y,z) \mathrm{d}x\mathrm{d}y\f]
 *  @param z     variable
 */
// ============================================================================
double
Ostap::Math::Bernstein3D::integrateXY
( const double z    ) const 
{
  if ( z < zmin () || z > zmax() ) { return 0 ; }
  //
  // const std::vector<double> fx ( nX () + 1 , ( xmax() - xmin () ) / ( nX () + 1 ) ) ;
  // const std::vector<double> fy ( nY () + 1 , ( ymax() - ymin () ) / ( nY () + 1 ) ) ;
  std::fill ( m_fx.begin() , m_fx.end() , ( xmax() - xmin () ) / ( nX () + 1 ) ) ;
  std::fill ( m_fy.begin() , m_fy.end() , ( ymax() - ymin () ) / ( nY () + 1 ) ) ;
  //
  // std::vector<double> fz ( nZ () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nZ () ; ++i )
  { m_fz[i] = m_bz[i] ( z ) ; }
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
/*  integral over x&z-dimensions
 *  \f[ \int_{x_{min}}^{x_{min}}
 *      \int_{z_{max}}^{z_{max}} \mathcal{B}(x,y,z) \mathrm{d}x\mathrm{d}z\f]
 *  @param y     variable
 */
// ============================================================================
double Ostap::Math::Bernstein3D::integrateXZ ( const double y    ) const 
{
  if ( y < ymin () || y > ymax() ) { return 0 ; }
  //
  // const std::vector<double> fx ( nX() + 1 , ( xmax() - xmin () ) / ( nX () + 1 ) ) ;
  std::fill ( m_fx.begin() , m_fx.end() , ( xmax() - xmin () ) / ( nX () + 1 ) ) ;
  //
  // std::vector<double> fy ( nY ()  + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nY () ; ++i )
  { m_fy[i] = m_by[i] ( y ) ; }
  //
  // const std::vector<double> fz ( nZ () + 1 , ( zmax() - zmin () ) / ( nZ () + 1 ) ) ;
  std::fill ( m_fz.begin() , m_fz.end() , ( zmax() - zmin () ) / ( nZ () + 1 ) ) ;
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
/* integral over y&z-dimensions
 *  \f[ \int_{y_{min}}^{y_{max}}
 *      \int_{z_{min}}^{z_{max}} \mathcal{B}(x,y,z) \mathrm{d}y\mathrm{d}z\f]
 *  @param x     variable
 */
// ============================================================================
double Ostap::Math::Bernstein3D::integrateYZ ( const double x    ) const 
{
  if ( x < xmin () || x > xmax() ) { return 0 ; }
  //
  // std::vector<double> fx ( nX () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nX () ; ++i )
  { m_fx[i] = m_bx[i] ( x ) ; }
  //
  // const std::vector<double> fy ( nY () + 1 , ( ymax() - ymin () ) / ( nY () + 1 ) ) ;
  std::fill ( m_fy.begin() , m_fy.end() , ( ymax() - ymin () ) / ( nY () + 1 ) ) ;
  //
  // const std::vector<double> fz ( nZ () + 1 , ( zmax() - zmin () ) / ( nZ () + 1 ) ) ;
  std::fill ( m_fz.begin() , m_fz.end() , ( zmax() - zmin () ) / ( nZ () + 1 ) ) ;
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
/* get the integral 
 *  \f$ \mathcal{B}(z) = 
 *   \int_{x_{min}}^{x_{max}}
 *   \int_{y_{min}}^{y_{max}} \mathcal{B}(x,y,z)dxdy \f$ 
 *  @see Ostap::Math::Bernstein3D::integrateXY  
 */
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Bernstein3D::integralXY () const
{
  for  ( unsigned int iz = 0 ; iz <= m_nz ; ++iz )
  {
    m_fz [ iz ] = 0 ;
    for ( unsigned short ix = 0 ; ix <= m_nx ; ++ix )
    {
      for ( unsigned short iy = 0 ; iy <= m_ny ; ++iy )
      { m_fz [ iz ] += par ( ix , iy , iz ) ; }
    }
  }
  const double scale = 1.0 * ( m_nz + 1 )/ ( m_zmax - m_zmin ) ;
  Ostap::Math::scale ( m_fz.begin(), m_fz.end() , scale ) ;
  return Bernstein ( m_fz , zmin() , zmax() ) ; 
}
// ============================================================================
/* get the integral 
 *  \f$ \mathcal{B}(y) = 
 *   \int_{x_{min}}^{x_{max}}
 *   \int_{z_{min}}^{z_{max}} \mathcal{B}(x,y,z)dxdz \f$ 
 *  @see Ostap::Math::Bernstein3D::integrateXZ  
 */
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Bernstein3D::integralXZ () const
{
  for  ( unsigned int iy = 0 ; iy <= m_ny ; ++iy )
  {
    m_fy [ iy ] = 0 ;
    for ( unsigned short ix = 0 ; ix <= m_nx ; ++ix )
    {
      for ( unsigned short iz = 0 ; iz <= m_nz ; ++iz )
      { m_fy [ iy ] += par ( ix , iy , iz ) ; }
    }
  }
  const double scale = 1.0 * ( m_ny + 1 )/ ( m_ymax - m_ymin ) ;
  Ostap::Math::scale ( m_fy.begin(), m_fy.end() , scale ) ;
  return Bernstein ( m_fy , ymin() , ymax() ) ; 
}
// ============================================================================
/* get the integral 
 *  \f$ \mathcal{B}(x) = 
 *   \int_{y_{min}}^{y_{max}}
 *   \int_{z_{min}}^{z_{max}} \mathcal{B}(x,y,z)dydz \f$ 
 *  @see Ostap::Math::Bernstein3D::integrateYZ  
 */
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Bernstein3D::integralYZ () const
{
  for  ( unsigned int ix = 0 ; ix <= m_nx ; ++ix )
  {
    m_fx [ ix ] = 0 ;
    for ( unsigned short iy = 0 ; iy <= m_ny ; ++iy )
    {
      for ( unsigned short iz = 0 ; iz <= m_nz ; ++iz )
      { m_fx [ ix ] += par ( ix , iy , iz ) ; }
    }
  }
  const double scale = 1.0 * ( m_nx + 1 )/ ( m_xmax - m_xmin ) ;
  Ostap::Math::scale ( m_fx.begin(), m_fx.end() , scale ) ;
  return Bernstein ( m_fx , xmin() , xmax() ) ; 
}
// ============================================================================
/*  get the integral 
 *  \f$ \mathcal{B}(z) = 
 *   \int_{x_{low}}^{x_{high}}
 *   \int_{y_{low}}^{y_{max}} \mathcal{B}(x,y,z)dxdy \f$ 
 *  @see Ostap::Math::Bernstein3D::integrateXY  
 */
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Bernstein3D::integralXY
( const double xlow , const double xhigh ,
  const double ylow , const double yhigh ) const
{
  // get integrals over X 
  for  ( unsigned short ix = 0 ; ix <= m_nx ; ++ix )
  { m_fx [ ix ] = m_bx [ ix ].integral ( xlow , xhigh ) ; }
  // get integrals over Y 
  for  ( unsigned short iy = 0 ; iy <= m_ny ; ++iy )
  { m_fy [ iy ] = m_by [ iy ].integral ( ylow , yhigh ) ; }
  //
  for  ( unsigned int iz = 0 ; iz <= m_nz ; ++iz )
  {
    m_fz [ iz ] = 0 ;
    for ( unsigned short ix = 0 ; ix <= m_nx ; ++ix )
    {
      for ( unsigned short iy = 0 ; iy <= m_ny ; ++iy )
      { m_fz [ iz ] += par ( ix , iy , iz ) * m_fx [ ix ] * m_fy [ iy ] ; }
    }
  }
  const double scale =
    1.0 * ( m_nx  + 1 ) * ( m_ny + 1 ) * ( m_nz + 1 ) /
    ( (  m_xmax - m_xmin ) * ( m_ymax  - m_ymin ) * ( m_zmax - m_zmin ) ) ;
  Ostap::Math::scale ( m_fz.begin(), m_fz.end() , scale ) ;
  return Bernstein ( m_fz , zmin() , zmax() ) ; 
}
// ============================================================================
/*  get the integral 
 *  \f$ \mathcal{B}(y) = 
 *   \int_{x_{low}}^{x_{high}}
 *   \int_{z_{low}}^{z_{high}} \mathcal{B}(x,y,z)dxdz \f$ 
 *  @see Ostap::Math::Bernstein3D::integrateXZ  
 */
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Bernstein3D::integralXZ
( const double xlow , const double xhigh ,
  const double zlow , const double zhigh ) const
{
  // get integrals over X 
  for  ( unsigned short ix = 0 ; ix <= m_nx ; ++ix )
  { m_fx [ ix ] = m_bx [ ix ].integral ( xlow , xhigh ) ; }
  // get integrals over Z 
  for  ( unsigned short iz = 0 ; iz <= m_nz ; ++iz )
  { m_fz [ iz ] = m_bz [ iz ].integral ( zlow , zhigh ) ; }
  //
  for  ( unsigned int iy = 0 ; iy <= m_ny ; ++iy )
  {
    m_fy [ iy ] = 0 ;
    for ( unsigned short ix = 0 ; ix <= m_nx ; ++ix )
    {
      for ( unsigned short iz = 0 ; iz <= m_nz ; ++iz )
      { m_fy [ iy ] += par ( ix , iy , iz ) * m_fx [ ix ] * m_fz [ iz ] ; }
    }
  }
  const double scale =
    1.0 * ( m_nx  + 1 ) * ( m_ny + 1 ) * ( m_nz + 1 ) /
    ( (  m_xmax - m_xmin ) * ( m_ymax  - m_ymin ) * ( m_zmax - m_zmin ) ) ;
  Ostap::Math::scale ( m_fy.begin(), m_fy.end() , scale ) ;
  return Bernstein ( m_fy , ymin() , ymax() ) ; 
}
// ============================================================================
/* get the integral 
 *  \f$ \mathcal{B}(x) = 
 *   \int_{y_{low}}^{y_{high}}
 *   \int_{z_{low}}^{z_{high}} \mathcal{B}(x,y,z)dydz \f$ 
 *  @see Ostap::Math::Bernstein3D::integrateYZ  
 */
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Bernstein3D::integralYZ
( const double ylow , const double yhigh ,
  const double zlow , const double zhigh ) const
{
  // get integrals over Y 
  for  ( unsigned short iy = 0 ; iy <= m_ny ; ++iy )
  { m_fy [ iy ] = m_by [ iy ].integral ( ylow , yhigh ) ; }
  // get integrals over Z 
  for  ( unsigned short iz = 0 ; iz <= m_nz ; ++iz )
  { m_fz [ iz ] = m_bz [ iz ].integral ( zlow , zhigh ) ; }
  //
  for  ( unsigned int ix = 0 ; ix <= m_nx ; ++ix )
  {
    m_fx [ ix ] = 0 ;
    for ( unsigned short iy = 0 ; iy <= m_ny ; ++iy )
    {
      for ( unsigned short iz = 0 ; iz <= m_nz ; ++iz )
      { m_fx [ ix ] += par ( ix , iy , iz ) * m_fy [ iy ] * m_fz [ iz ] ; }
    }
  }
  const double scale =
    1.0 * ( m_nx  + 1 ) * ( m_ny + 1 ) * ( m_nz + 1 ) /
    ( (  m_xmax - m_xmin ) * ( m_ymax  - m_ymin ) * ( m_zmax - m_zmin ) ) ;
  Ostap::Math::scale ( m_fx.begin(), m_fx.end() , scale ) ;
  return Bernstein ( m_fx , xmin() , xmax() ) ; 
}
// ============================================================================
/*  get the integral 
 *  \f$ \mathcal{B}(y,z) = 
 *   \int_{x_{min}}^{x_{max}}
 *   \mathcal{B}(x,y,z)dx \f$ 
 *  @see Ostap::Math::Bernstein3D::integrateX  
 */
// ============================================================================
Ostap::Math::Bernstein2D
Ostap::Math::Bernstein3D::integralX() const
{
  // crete the output object 
  Ostap::Math::Bernstein2D b2 { m_ny   , m_nz   ,
                                m_ymin , m_ymax ,
                                m_zmin , m_zmax };
  // fill it! 
  for ( unsigned short iy = 0 ; iy <= m_ny ; ++iy )
  {
    for ( unsigned short iz = 0 ; iz <= m_nz ; ++iz )
    {
      double pp = 0.0 ;
      for  ( unsigned int ix = 0 ; ix <= m_nx ; ++ix )
      { pp += par ( ix , iy , iz ) ; }
      //
      b2.setPar ( iy , iz , pp ) ;
    }
  }
  return b2 ;
}
// ============================================================================
/*  get the integral 
 *  \f$ \mathcal{B}(x,z) = 
 *   \int_{y_{min}}^{y_{max}}
 *   \mathcal{B}(x,y,z)dy \f$ 
 *  @see Ostap::Math::Bernstein3D::integrateY  
 */
// ============================================================================
Ostap::Math::Bernstein2D
Ostap::Math::Bernstein3D::integralY() const
{
  // crete the output object 
  Ostap::Math::Bernstein2D b2 { m_nx   , m_nz   ,
                                m_xmin , m_xmax ,
                                m_zmin , m_zmax };
  // fill it! 
  for ( unsigned short ix = 0 ; ix <= m_nx ; ++ix )
  {
    for ( unsigned short iz = 0 ; iz <= m_nz ; ++iz )
    {
      double pp = 0.0 ;
      for  ( unsigned int iy = 0 ; iy <= m_ny ; ++iy )
      { pp += par ( ix , iy , iz ) ; }
      //
      b2.setPar ( ix , iz , pp ) ;
    }
  }
  return b2 ;
}
// ============================================================================
/*  get the integral 
 *  \f$ \mathcal{B}(x,y) = 
 *   \int_{z_{min}}^{z_{max}}
 *   \mathcal{B}(x,y,z)dz \f$ 
 *  @see Ostap::Math::Bernstein3D::integrateZ  
 */
// ============================================================================
Ostap::Math::Bernstein2D
Ostap::Math::Bernstein3D::integralZ() const
{
  // crete the output object 
  Ostap::Math::Bernstein2D b2 { m_nx   , m_ny   ,
                                m_xmin , m_xmax ,
                                m_ymin , m_ymax };
  // fill it! 
  for ( unsigned short ix = 0 ; ix <= m_nx ; ++ix )
  {
    for ( unsigned short iy = 0 ; iy <= m_ny ; ++iy )
    {
      double pp = 0.0 ;
      for  ( unsigned int iz = 0 ; iz <= m_nz ; ++iz )
      { pp += par ( ix , iy , iz ) ; }
      //
      b2.setPar ( ix , iy , pp ) ;
    }
  }
  return b2 ;
}
// ============================================================================
/*  get the integral 
 *  \f$ \mathcal{B}(y,z) = 
 *   \int_{x_{low}}^{x_{high}}
 *   \mathcal{B}(x,y,z)dx \f$ 
 *  @see Ostap::Math::Bernstein3D::integrateX  
 */
// ============================================================================
Ostap::Math::Bernstein2D
Ostap::Math::Bernstein3D::integralX
( const double xlow  ,
  const double xhigh ) const 
{
  // perform integration 
  for ( unsigned short ix = 0 ; ix <= m_nx  ; ++ix )
  { m_fx [ ix ] = m_bx [ ix ] .integral ( xlow , xhigh ) ; }
  
  // crete the output object 
  Ostap::Math::Bernstein2D b2 { m_ny   , m_nz   ,
                                m_ymin , m_ymax ,
                                m_zmin , m_zmax };
  // fill it! 
  for ( unsigned short iy = 0 ; iy <= m_ny ; ++iy )
  {
    for ( unsigned short iz = 0 ; iz <= m_nz ; ++iz )
    {
      double pp = 0.0 ;
      for  ( unsigned int ix = 0 ; ix <= m_nx ; ++ix )
      { pp += par ( ix , iy , iz ) * m_fx [ ix ] ; }
      //
      b2.setPar ( iy , iz , pp ) ;
    }
  }
  //
  b2 *= ( m_nx + 1 ) / ( m_xmax - m_xmin ) ;
  return b2 ;
}
// ============================================================================
/*  get the integral 
 *  \f$ \mathcal{B}(x,z) = 
 *   \int_{y_{low}}^{y_{high}}
 *   \mathcal{B}(x,y,z)dy \f$ 
 *  @see Ostap::Math::Bernstein3D::integrateY  
 */
// ============================================================================
Ostap::Math::Bernstein2D
Ostap::Math::Bernstein3D::integralY
( const double ylow  ,
  const double yhigh ) const 
{
  // perform integration 
  for ( unsigned short iy = 0 ; iy <= m_ny  ; ++iy )
  { m_fy [ iy ] = m_by [ iy ] .integral ( ylow , yhigh ) ; }
  
  // crete the output object 
  Ostap::Math::Bernstein2D b2 { m_nx   , m_nz   ,
                                m_xmin , m_xmax ,
                                m_zmin , m_zmax };
  // fill it! 
  for ( unsigned short ix = 0 ; ix <= m_nx ; ++ix )
  {
    for ( unsigned short iz = 0 ; iz <= m_nz ; ++iz )
    {
      double pp = 0.0 ;
      for  ( unsigned int iy = 0 ; iy <= m_ny ; ++iy )
      { pp += par ( ix , iy , iz ) * m_fy [ iy ] ; }
      //
      b2.setPar ( ix , iz , pp ) ;
    }
  }
  //
  b2 *= ( m_ny + 1 ) / ( m_ymax - m_ymin ) ;
  return b2 ;
}
// ============================================================================
/** get the integral 
 *  \f$ \mathcal{B}(x,y) = 
 *   \int_{z_{low}}^{z_{high}}
 *   \mathcal{B}(x,y,z)dz \f$ 
 *  @see Ostap::Math::Bernstein3D::integrateZ  
 */
// ============================================================================
Ostap::Math::Bernstein2D
Ostap::Math::Bernstein3D::integralZ
( const double zlow  ,
  const double zhigh ) const 
{
  // perform integration 
  for ( unsigned short iz = 0 ; iz <= m_nz  ; ++iz )
  { m_fz [ iz ] = m_bz [ iz ] .integral ( zlow , zhigh ) ; }
  
  // crete the output object 
  Ostap::Math::Bernstein2D b2 { m_nx   , m_ny   ,
                                m_xmin , m_xmax ,
                                m_ymin , m_ymax };
  // fill it! 
  for ( unsigned short ix = 0 ; ix <= m_nx ; ++ix )
  {
    for ( unsigned short iy = 0 ; iy <= m_ny ; ++iy )
    {
      double pp = 0.0 ;
      for  ( unsigned int iz = 0 ; iz <= m_nz ; ++iz )
      { pp += par ( ix , iy , iz ) * m_fz [ iz ] ; }
      //
      b2.setPar ( ix , iy , pp ) ;
    }
  }
  //
  b2 *= ( m_nz + 1 ) / ( m_zmax - m_zmin ) ;
  return b2 ;
}



// ============================================================================
/** update Bernstein expansion by addition of one "event" with 
 *  the given weight
 *  @code
 *  Bernstein3D sum = ... ;
 *  for ( auto x : .... ) { sum.fill ( x , y , z ) ; }
 *  @endcode
 *  This is a useful function to make an unbinned parameterization 
 *  of certain distribution and/or efficiency 
 */
// ============================================================================
bool Ostap::Math::Bernstein3D::fill 
( const double x      , 
  const double y      , 
  const double z      , 
  const double weight )
{
  // no update 
  if      ( x < m_xmin || x > m_xmax ) { return false ; }
  else if ( y < m_ymin || y > m_ymax ) { return false ; }
  else if ( z < m_zmin || z > m_zmax ) { return false ; }
  else if ( s_zero ( weight )        ) { return true  ; }
  //
  const long double w = weight * 1.0L / 
    // ( 1.0L * ( m_ymax - m_ymin ) * ( m_xmax - m_xmin ) * ( m_nx + 1 ) * ( m_ny + 1 ) ) ;
    ( 1.0L * ( m_nx + 1 ) * ( m_ny + 1 )* ( m_nz  + 1 )  ) ;
  //
  const double xx  =  tx ( x ) ;
  const double yy  =  ty ( y ) ;
  const double zz  =  tz ( z ) ;
  //
  typedef Ostap::Math::BernsteinDualBasis::Basis Basis ;
  const Basis* basisx = Ostap::Math::BernsteinDualBasis::basis ( m_nx ) ; 
  const Basis* basisy = Ostap::Math::BernsteinDualBasis::basis ( m_ny ) ;
  const Basis* basisz = Ostap::Math::BernsteinDualBasis::basis ( m_nz ) ;
  //
  static const std::string s_MSGX { "Cannot aquire valid Bernstein dual basis/x" } ;
  static const std::string s_MSGY { "Cannot aquire valid Bernstein dual basis/y" } ;
  static const std::string s_MSGZ { "Cannot aquire valid Bernstein dual basis/z" } ;
  static const std::string s_TAG  { "Ostap::Math::Bernstein3D::fill"          } ;
  //
  Ostap::Assert ( basisx && basisx->size() == m_nx + 1 , s_MSGX , s_TAG , 950 ) ;
  Ostap::Assert ( basisy && basisy->size() == m_ny + 1 , s_MSGY , s_TAG , 950 ) ;
  Ostap::Assert ( basisz && basisz->size() == m_nz + 1 , s_MSGZ , s_TAG , 950 ) ;
  //
  std::transform ( basisx->begin () , basisx->end () , m_fx.begin() ,
                   [xx]( const auto& b ) -> double { return b ( xx  ) ; } ) ;
  std::transform ( basisy->begin () , basisy->end () , m_fy.begin() ,
                   [yy]( const auto& b ) -> double { return b ( yy  ) ; } ) ;
  std::transform ( basisz->begin () , basisz->end () , m_fz.begin() ,
                   [zz]( const auto& b ) -> double { return b ( zz  ) ; } ) ;
  //
  for ( unsigned short ix = 0 ; ix <= m_nx  ; ++ix  )
  {
    for ( unsigned short iy = 0 ; iy <= m_ny  ; ++iy  )
    {
      for ( unsigned short iz = 0 ; iz <= m_nz  ; ++iz  )
      {
        const std::size_t k = index  ( ix  , iy , iz ) ;
        m_pars [ k ] += w * m_fx [ ix ] * m_fy [ iy ] * m_fz [ iz  ] ;
      }
    }
  }
  //
  return true ;
}
// ============================================================================
Ostap::Math::Bernstein3D&
Ostap::Math::Bernstein3D::operator+=( const double a )
{
  if   ( s_zero ( a ) ) { return *this ; }
  Ostap::Math::shift ( m_pars , a ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Bernstein3D&
Ostap::Math::Bernstein3D::operator*=( const double a )
{
  if      ( s_equal ( a , 1 ) ) { return *this ; }
  else if ( s_zero  ( a     ) ) { std::fill ( m_pars.begin() , m_pars.end() , 0 ) ; }
  Ostap::Math::scale ( m_pars , a ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Bernstein3D&
Ostap::Math::Bernstein3D::operator-=( const double a )
{
  if ( s_zero ( a ) ) { return *this ; }
  Ostap::Math::shift ( m_pars , -a ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Bernstein3D&
Ostap::Math::Bernstein3D::operator/=( const double a )
{
  if   ( s_equal ( a , 1 ) ) { return *this ; }
  Ostap::Math::scale ( m_pars , 1/a ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Bernstein3D
Ostap::Math::Bernstein3D::operator-() const
{
  Bernstein3D b ( *this ) ;
  Ostap::Math::negate ( b.m_pars ) ;
  return b ;
}
// ============================================================================
// Sum of Bernstein polynomial and a constant
// ============================================================================
Ostap::Math::Bernstein3D
Ostap::Math::Bernstein3D::__add__   ( const double value ) const
{ return (*this) + value ; }
// ============================================================================
// Sum of Bernstein polynomial and a constant
// ============================================================================
Ostap::Math::Bernstein3D
Ostap::Math::Bernstein3D::__radd__  ( const double value ) const
{ return value + (*this) ; }
// ============================================================================
// Product of Bernstein polynomial and a constant
// ============================================================================
Ostap::Math::Bernstein3D
Ostap::Math::Bernstein3D::__mul__   ( const double value ) const
{ return (*this) * value ; }
// ============================================================================
// Product of Bernstein polynomial and a constant
// ============================================================================
Ostap::Math::Bernstein3D
Ostap::Math::Bernstein3D::__rmul__  ( const double value ) const
{ return value * (*this) ; }
// ============================================================================
// Subtract a constant from Benrstein polynomial
// ============================================================================
Ostap::Math::Bernstein3D
Ostap::Math::Bernstein3D::__sub__  ( const double value ) const
{ return (*this) - value ; }
// ============================================================================
// Subtract Bernstein polynomial from a constant
// ============================================================================
Ostap::Math::Bernstein3D
Ostap::Math::Bernstein3D::__rsub__ ( const double value ) const
{ return value - (*this) ; }
// ============================================================================
// Divide Benrstein polynomial by a constant
// ============================================================================
Ostap::Math::Bernstein3D
Ostap::Math::Bernstein3D:: __truediv__   ( const double value ) const
{ return (*this) / value ; }
// ============================================================================
// Negate Bernstein polynomial
// ============================================================================
Ostap::Math::Bernstein3D
Ostap::Math::Bernstein3D::__neg__ ()  const
{ return -(*this); }
// ============================================================================
// get the tag value 
// ============================================================================
std::size_t Ostap::Math::Bernstein3D::tag () const  // get the tag value 
{
  return Ostap::Utils::hash_combiner
    ( Ostap::Utils::hash_range ( m_pars )  ,  
      m_nx   , m_ny   , m_nz ,
      m_xmin , m_xmax ,
      m_ymin , m_ymax , 
      m_zmin , m_zmax ) ; }
// ============================================================================

// ============================================================================
// Add       polynomials (with the same domain!)
// ============================================================================
Ostap::Math::Bernstein3D&
Ostap::Math::Bernstein3D::isum
( const Ostap::Math::Bernstein3D& other ) 
{
  if ( this == &other ) { *this *= 2 ; return *this; }
  //
  Ostap::Assert ( m_nx == other.m_nx && 
                  m_ny == other.m_ny && 
                  m_nz == other.m_nz && 
                  s_equal ( xmin() , other.xmin() ) &&
                  s_equal ( xmax() , other.xmax() ) && 
                  s_equal ( ymin() , other.ymin() ) &&
                  s_equal ( ymax() , other.ymax() ) && 
                  s_equal ( zmin() , other.zmin() ) &&
                  s_equal ( zmax() , other.zmax() )               ,
                  "Cannot isum Bernstein3D with different domains"  , 
                  "Ostap::Math::Bernstei3D"                        , 
                  Ostap::StatusCode ( 527 )                       )  ;
  //
  for ( unsigned short i = 0 ; i < npars() ; ++i ) 
  { m_pars[i] += other.m_pars [ i ] ; }
  //
  return *this ;
}
// ============================================================================
// Subtract polynomials (with the same domain!)
// ============================================================================
Ostap::Math::Bernstein3D&
Ostap::Math::Bernstein3D::isub
( const Ostap::Math::Bernstein3D& other ) 
{
  if ( this == &other ) { *this *= 0 ; return *this; }
  //
  Ostap::Assert ( m_nx == other.m_nx && 
                  m_ny == other.m_ny && 
                  m_nz == other.m_nz && 
                  s_equal ( xmin() , other.xmin() ) &&
                  s_equal ( xmax() , other.xmax() ) && 
                  s_equal ( ymin() , other.ymin() ) &&
                  s_equal ( ymax() , other.ymax() ) && 
                  s_equal ( zmin() , other.zmin() ) &&
                  s_equal ( zmax() , other.zmax() )               ,
                  "Cannot isub Bernstein3D with different domains"  , 
                  "Ostap::Math::Bernstei3D"                        , 
                  Ostap::StatusCode ( 528 )                       )  ;
  //
  for ( unsigned short i = 0 ; i < npars() ; ++i ) 
  { m_pars[i] += other.m_pars [ i ] ; }
  //
  return *this ;
}





// ============================================================================
// 3S symmetric polynomial
// ============================================================================

// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::Bernstein3DSym::Bernstein3DSym 
( const unsigned short       N     ,
  const double               xmin  ,
  const double               xmax  )
  : Ostap::Math::Parameters ( ( N + 1 ) * ( N + 2 ) * ( N + 3 ) / 6  )
  , m_n    ( N )
    //
    //
  , m_xmin ( std::min ( xmin , xmax ) )
  , m_xmax ( std::max ( xmin , xmax ) )
    //
  , m_b    ()
    //
  , m_fx ( N + 1  , 0.0 )
  , m_fy ( N + 1  , 0.0 )
  , m_fz ( N + 1  , 0.0 )    
{
  //
  typedef  Ostap::Math::Bernstein::Basic BB ;
  for ( unsigned short ix = 0 ; ix <= N ; ++ix )
  { m_b.push_back ( Bernstein ( BB ( ix , N ) , xmin , xmax ) ) ; }
  //
}
// ============================================================================
// constructor from parameters 
// ============================================================================
Ostap::Math::Bernstein3DSym::Bernstein3DSym 
( const std::vector<double>& pars  ,
  const unsigned short       N     ,
  const double               xmin  ,
  const double               xmax  )
  : Ostap::Math::Bernstein3DSym ( N , xmin , xmax )
{
  setPars ( pars ) ;
}
// ============================================================================
// swap  two 3D-polynomials 
// ============================================================================
void Ostap::Math::Bernstein3DSym::swap
( Ostap::Math::Bernstein3DSym&  right ) 
{
  Ostap::Math::Parameters::swap ( right ) ;
  //
  std::swap ( m_pars , right.m_pars  ) ;
  std::swap ( m_xmin , right.m_xmin  ) ;
  std::swap ( m_xmax , right.m_xmax  ) ;
  std::swap ( m_b    , right.m_b     ) ;
  std::swap ( m_fx   , right.m_fx    ) ;
  std::swap ( m_fy   , right.m_fy    ) ;
  std::swap ( m_fz   , right.m_fx    ) ;
}
// ============================================================================
// helper function to make calculations
// ============================================================================
double Ostap::Math::Bernstein3DSym::calculate
( const std::vector<double>& fx , 
  const std::vector<double>& fy , 
  const std::vector<double>& fz ) const 
{
  double       result = 0 ;
  for  ( unsigned short ix = 0 ; ix <= nX ()  ; ++ix )
  {
    for  ( unsigned short iy = 0 ; iy <= ix ; ++iy )
    { 
      for  ( unsigned short iz = 0 ; iz <= iy ; ++iz )
      { 
        double r = 0 ;
        if      ( ix == iy && iy == iz ) 
        { 
          r += fx[ix] * fy[iy] * fz[iz]; 
        }
        else if ( ix == iy ) 
        { 
          r += fx[ix] * fy[iy] * fz[iz]
            +  fx[iz] * fy[ix] * fz[iy]
            +  fx[ix] * fy[iz] * fz[iy] ; 
        }
        else if ( iy == iz ) 
        { 
          r += fx[ix] * fy[iy] * fz[iz] 
            +  fx[iy] * fy[ix] * fz[iz]  
            +  fx[iy] * fy[iz] * fz[ix] ; 
        }
        else 
        {
          r += fx[ix] * fy[iy] * fz[iz]
            +  fx[ix] * fy[iz] * fz[iy]
            +  fx[iy] * fy[ix] * fz[iz] 
            +  fx[iy] * fy[iz] * fz[ix] 
            +  fx[iz] * fy[ix] * fz[iy]
            +  fx[iz] * fy[iy] * fz[ix] ;          
        }
        // 
        result += r * par( ix , iy , iz ) ;
      }
    }  
  }
  //
  const double scalex = ( nX () + 1 ) / ( xmax() - xmin() ) ;
  const double scaley = scalex ;
  const double scalez = scalex ;
  //
  return result * scalex * scaley * scalez ;
}
// ============================================================================
// get the value
// ============================================================================
double Ostap::Math::Bernstein3DSym::evaluate
( const double x ,
  const double y , 
  const double z ) const
{
  /// the trivial cases
  if ( x < xmin () || x > xmax () ) { return 0.0        ; }
  if ( y < ymin () || y > ymax () ) { return 0.0        ; }
  if ( z < zmin () || z > zmax () ) { return 0.0        ; }
  //
  if      ( 0 == npars ()       ) { return 0.0        ; }
  else if ( 1 == npars ()       ) 
  { 
    const double scale = ( nX () + 1 ) / ( xmax() - xmin() ) ;
    //
    return m_pars [0] * scale * scale * scale ;
  }
  ///
  // std::vector<double> fx ( nX () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nX () ; ++i )
  { m_fx[i] = m_b [i] ( x )  ; }
  //
  // std::vector<double> fy ( nY () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nY ()  ; ++i )
  { m_fy[i] = m_b [i] ( y )  ; }
  //
  // std::vector<double> fz ( nZ () + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nZ () ; ++i )
  { m_fz[i] = m_b [i] ( z )  ; }
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
/** get the integral over 3D-region
 *  \f[  x_{min} < x < x_{max}, y_{min}< y< y_{max} , z_{min} < z < z_{max}\f]
 */
// ============================================================================
double Ostap::Math::Bernstein3DSym::integral() const
{ 
  double       result = 0 ;
  for  ( unsigned short ix = 0 ; ix <= nX () ; ++ix )
  {
    for  ( unsigned short iy = 0 ; iy <= ix ; ++iy )
    { 
      for  ( unsigned short iz = 0 ; iz <= iy ; ++iz )
      { 
        unsigned short r =
          ( ix == iy && iy == iz ) ? 1 :
          ( ix == iy || iy == iz ) ? 3 : 6 ;
        // 
        result += r  *  par( ix , iy , iz ) ;
      }
    }  
  }
  //
  return result ;
}
// ============================================================================
/* get the integral over 3D-region
 *  \f[ \int_{x_{low}}^{x_{high}}
 *      \int_{y_{low}}^{y_{high}} 
 *      \int_{z_{low}}^{z_{high}} 
 *\mathcal{B}(x,y,z) \mathrm{d}x\mathrm{d}y\mathrm{d}z\f]
 *  @param xlow  low  edge in x
 *  @param xhigh high edge in x
 *  @param ylow  low  edge in y
 *  @param yhigh high edge in y
 *  @param zlow  low  edge in z
 *  @param zhigh high edge in z
 */
// ============================================================================
double Ostap::Math::Bernstein3DSym::integral
( const double xlow , const double xhigh ,
  const double ylow , const double yhigh ,
  const double zlow , const double zhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ||
            s_equal ( ylow , yhigh ) ||
            s_equal ( zlow , zhigh ) ) { return 0 ; }
  //
  else if ( s_equal ( xlow  , xmin () ) &&
            s_equal ( xhigh , xmax () ) &&
            s_equal ( ylow  , ymin () ) &&
            s_equal ( yhigh , ymax () ) &&  
            s_equal ( zlow  , zmin () ) &&
            s_equal ( zhigh , zmax () ) )  { return integral () ; }
  //
  else if ( xlow  > xhigh ) 
  { return -1*integral ( xhigh , xlow  , ylow  , yhigh , zlow  , zhigh ) ; }
  else if ( ylow  > yhigh ) 
  { return -1*integral ( xlow  , xhigh , yhigh , ylow  , zlow  , zhigh ) ; }
  else if ( zlow  > zhigh ) 
  { return -1*integral ( xlow  , xhigh , ylow  , yhigh , zhigh , zlow  ) ; }
  //
  else if ( xhigh <  xmin () || xlow >  xmax() ) { return 0 ; }
  else if ( yhigh <  ymin () || ylow >  ymax() ) { return 0 ; }
  else if ( zhigh <  zmin () || zlow >  zmax() ) { return 0 ; }
  //
  const double  x_low  = std::max ( xmin() , xlow  ) ;
  const double  x_high = std::min ( xmax() , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  const double  y_low  = std::max ( ymin() , ylow  ) ;
  const double  y_high = std::min ( ymax() , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  const double  z_low  = std::max ( zmin() , zlow  ) ;
  const double  z_high = std::min ( zmax() , zhigh ) ;
  if ( z_low >= z_high ) { return 0 ; }
  //
  // std::vector<double> fx ( nX () + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nX () ; ++i )
  { m_fx[i] = m_b [i].integral ( x_low , x_high ) ; }
  //
  // std::vector<double> fy ( nY () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nY () ; ++i )
  { m_fy[i] = m_b [i].integral ( y_low , y_high ) ; }
  //
  // std::vector<double> fz ( nZ () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nZ ()  ; ++i )
  { m_fz[i] = m_b [i].integral ( z_low , z_high ) ; }
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
/*  integral over x-dimension
 *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y,z) \mathrm{d}y\f]
 *  @param y     variable
 *  @param z     variable
 *  @param ylow  low  edge in y
 *  @param yhigh high edge in y
 */
// ============================================================================
double Ostap::Math::Bernstein3DSym::integrateX
( const double y    ,
  const double z    ,
  const double xlow , const double xhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; }
  else if ( xlow  > xhigh  ) { return -1*integrateX ( y , z , xhigh , xlow ) ; }
  else if ( xhigh <= xmin () || xlow >= xmax() ) { return 0 ; }
  else if ( y     <  ymin () || y    >  ymax() ) { return 0 ; }
  else if ( z     <  zmin () || z    >  zmax() ) { return 0 ; }
  else if ( s_equal ( xlow  , xmin () ) &&
            s_equal ( xhigh , xmax () )         ) { return integrateX ( y ,  z ) ; }
  //
  const double  x_low  = std::max ( xmin() , xlow  ) ;
  const double  x_high = std::min ( xmax() , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  // std::vector<double> fx ( nX () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nX () ; ++i )
  { m_fx[i] = m_b [i].integral ( x_low , x_high ) ; }
  //
  // std::vector<double> fy ( nY () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nY ()   ; ++i )
  { m_fy[i] = m_b [i] ( y ) ; }
  //
  // std::vector<double> fz ( nZ () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nZ ()  ; ++i )
  { m_fz[i] = m_b [i] ( z ) ; }
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
double
Ostap::Math::Bernstein3DSym::integrateX
( const double y , 
  const double z ) const
{
  if      ( y < ymin () || y > ymax() ) { return 0 ; }
  else if ( z < zmin () || z > zmax() ) { return 0 ; }
  //
  // const std::vector<double> fx ( nX () + 1 , ( xmax() - xmin () ) / ( nX () + 1 ) ) ;
  std::fill ( m_fx.begin() , m_fx.end() , ( xmax() - xmin () ) / ( nX () + 1 ) ) ;
  //
  // std::vector<double> fy ( nY () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nY () ; ++i )
  { m_fy[i] = m_b [i] ( y ) ; }
  //
  // std::vector<double> fz ( nZ () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nZ () ; ++i )
  { m_fz[i] = m_b [i] ( z ) ; }
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
/*  integral over x&y-dimensions
 *  \f[ \int_{x_{low}}^{x_{high}}
 *      \int_{y_{low}}^{y_{high}} \mathcal{B}(x,y,z) \mathrm{d}x\mathrm{d}y\f]
 *  @param z     variable
 *  @param xlow  low  edge in x
 *  @param xhigh high edge in x
 *  @param ylow  low  edge in y
 *  @param yhigh high edge in y
 */
// ============================================================================
double Ostap::Math::Bernstein3DSym::integrateXY
( const double z    ,                          
  const double xlow , const double xhigh ,
  const double ylow , const double yhigh ) const 
{
  if      ( s_equal ( xlow  , xhigh ) ) { return 0 ; }
  else if ( s_equal ( ylow  , yhigh ) ) { return 0 ; }
  else if ( xlow  >  xhigh ) { return -1*integrateXY ( z , xhigh , xlow  , ylow  , yhigh ) ; }
  else if ( ylow  >  yhigh ) { return -1*integrateXY ( z , xlow  , xhigh , yhigh , ylow  ) ; }
  else if ( z     <  zmin () || z    >  zmax() ) { return 0 ; }
  else if ( s_equal ( xlow  , xmin () ) &&
            s_equal ( xhigh , xmax () ) &&
            s_equal ( ylow  , ymin () ) &&
            s_equal ( yhigh , ymax () ) ) { return integrateXY ( z ) ; }
  //  //
  const double  x_low  = std::max ( xmin() , xlow  ) ;
  const double  x_high = std::min ( xmax() , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  const double  y_low  = std::max ( ymin() , ylow  ) ;
  const double  y_high = std::min ( ymax() , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  // std::vector<double> fx ( nX () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nX () ; ++i )
  { m_fx[i] = m_b [i].integral ( x_low , x_high ) ; }
  //
  // std::vector<double> fy ( nY () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nY () ; ++i )
  { m_fy[i] = m_b [i].integral ( y_low , y_high ) ; }
  //
  // std::vector<double> fz ( nZ () + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nZ () ; ++i )
  { m_fz[i] = m_b [i] ( z ) ; }
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
/*  integral over x&y-dimensions
 *  \f[ \int_{x_{min}}^{x_{max}}
 *      \int_{y_{min}}^{y_{max}} \mathcal{B}(x,y,z) \mathrm{d}x\mathrm{d}y\f]
 *  @param z     variable
 */
// ============================================================================
double
Ostap::Math::Bernstein3DSym::integrateXY
( const double z    ) const 
{
  if ( z < zmin () || z > zmax() ) { return 0 ; }
  //
  // const std::vector<double> fx ( nX() + 1 , ( xmax() - xmin () ) / ( nX () + 1 ) ) ;
  // const std::vector<double> fy ( nY() + 1 , ( ymax() - ymin () ) / ( nY () + 1 ) ) ;
  std::fill ( m_fx.begin() , m_fx.end() , ( xmax() - xmin () ) / ( nX () + 1 ) ) ;
  std::fill ( m_fy.begin() , m_fy.end() , ( ymax() - ymin () ) / ( nY () + 1 ) ) ;
  //
  // std::vector<double> fz ( nZ () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nZ () ; ++i )
  { m_fz[i] = m_b [i] ( z ) ; }
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
Ostap::Math::Bernstein3DSym&
Ostap::Math::Bernstein3DSym::operator+=( const double a )
{
  if   ( s_zero ( a ) ) { return *this ; }
  Ostap::Math::shift ( m_pars , a ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Bernstein3DSym&
Ostap::Math::Bernstein3DSym::operator*=( const double a )
{
  if      ( s_equal ( a , 1 ) ) { return *this ; }
  else if ( s_zero  ( a     ) ) { std::fill ( m_pars.begin() , m_pars.end() , 0 ) ; }
  Ostap::Math::scale ( m_pars , a ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Bernstein3DSym&
Ostap::Math::Bernstein3DSym::operator-=( const double a )
{
  if ( s_zero ( a ) ) { return *this ; }
  Ostap::Math::shift ( m_pars , -a ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Bernstein3DSym&
Ostap::Math::Bernstein3DSym::operator/=( const double a )
{
  if   ( s_equal ( a , 1 ) ) { return *this ; }
  Ostap::Math::scale ( m_pars , 1/a ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Bernstein3DSym
Ostap::Math::Bernstein3DSym::operator-() const
{
  Bernstein3DSym b ( *this ) ;
  Ostap::Math::negate ( b.m_pars ) ;
  return b ;
}
// ============================================================================
// Sum of Bernstein polynomial and a constant
// ============================================================================
Ostap::Math::Bernstein3DSym
Ostap::Math::Bernstein3DSym::__add__   ( const double value ) const
{ return (*this) + value ; }
// ============================================================================
// Sum of Bernstein polynomial and a constant
// ============================================================================
Ostap::Math::Bernstein3DSym
Ostap::Math::Bernstein3DSym::__radd__  ( const double value ) const
{ return value + (*this) ; }
// ============================================================================
// Product of Bernstein polynomial and a constant
// ============================================================================
Ostap::Math::Bernstein3DSym
Ostap::Math::Bernstein3DSym::__mul__   ( const double value ) const
{ return (*this) * value ; }
// ============================================================================
// Product of Bernstein polynomial and a constant
// ============================================================================
Ostap::Math::Bernstein3DSym
Ostap::Math::Bernstein3DSym::__rmul__  ( const double value ) const
{ return value * (*this) ; }
// ============================================================================
// Subtract a constant from Benrstein polynomial
// ============================================================================
Ostap::Math::Bernstein3DSym
Ostap::Math::Bernstein3DSym::__sub__  ( const double value ) const
{ return (*this) - value ; }
// ============================================================================
// Subtract Bernstein polynomial from a constant
// ============================================================================
Ostap::Math::Bernstein3DSym
Ostap::Math::Bernstein3DSym::__rsub__ ( const double value ) const
{ return value - (*this) ; }
// ============================================================================
// Divide Benrstein polynomial by a constant
// ============================================================================
Ostap::Math::Bernstein3DSym
Ostap::Math::Bernstein3DSym:: __truediv__   ( const double value ) const
{ return (*this) / value ; }
// ============================================================================
// Negate Bernstein polynomial
// ============================================================================
Ostap::Math::Bernstein3DSym
Ostap::Math::Bernstein3DSym::__neg__ ()  const
{ return -(*this); }
// ============================================================================
// get the tag value 
// ============================================================================
std::size_t Ostap::Math::Bernstein3DSym::tag () const  // get the hash value 
{
  return Ostap::Utils::hash_combiner 
    ( Ostap::Utils::hash_range ( m_pars ) , 
      m_n    , 
      m_xmin , m_xmax ) ; 
}
// ============================================================================



// ============================================================================
// 3D  with X<->Y symmetry 
// ============================================================================

// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::Bernstein3DMix::Bernstein3DMix
( const unsigned short       N     ,
  const unsigned short       Nz    ,
  const double               xmin  ,
  const double               xmax  ,
  const double               zmin  ,
  const double               zmax  )
  : Ostap::Math::Parameters ( ( N + 1 ) * ( N + 2 ) * ( Nz + 1 ) / 2 )
  , m_n    ( N  )
  , m_nz   ( Nz )
    //
  , m_xmin ( std::min ( xmin , xmax ) )
  , m_xmax ( std::max ( xmin , xmax ) )
  , m_zmin ( std::min ( zmin , zmax ) )
  , m_zmax ( std::max ( zmin , zmax ) )
    //
  , m_b    ()
  , m_bz   ()
    //
  , m_fx ( N  + 1 , 0.0 ) 
  , m_fy ( N  + 1 , 0.0 ) 
  , m_fz ( Nz + 1 , 0.0 ) 
{
  //
  typedef  Ostap::Math::Bernstein::Basic BB ;
  for ( unsigned short ix = 0 ; ix <= N  ; ++ix )
  { m_b .push_back ( Bernstein ( BB ( ix , N  ) , xmin , xmax ) ) ; }
  for ( unsigned short iz = 0 ; iz <= Nz ; ++iz )
  { m_bz.push_back ( Bernstein ( BB ( iz , Nz ) , zmin , zmax ) ) ; }
  //
}
// ============================================================================
// constructor from parameters  
// ============================================================================
Ostap::Math::Bernstein3DMix::Bernstein3DMix
( const std::vector<double>& pars ,
  const unsigned short       N     ,
  const unsigned short       Nz    ,
  const double               xmin  ,
  const double               xmax  ,
  const double               zmin  ,
  const double               zmax  )
  : Ostap::Math::Bernstein3DMix ( N    , Nz   ,
                                  xmin , xmax ,
                                  zmin , zmax )
{
  setPars ( pars ) ;
}
// ============================================================================
// constructor from symmetric polynomial
// ============================================================================
Ostap::Math::Bernstein3DMix::Bernstein3DMix
( const Ostap::Math::Bernstein3DSym& right ) 
  : Ostap::Math::Parameters ( ( right.nZ() + 1 ) * ( right.nX() + 2 ) * ( right.nZ() + 1 ) / 2  )
    //
  , m_n    ( right.nX() )
  , m_nz   ( right.nZ() )
    //
    //
  , m_xmin ( right.xmin () )
  , m_xmax ( right.xmin () )
  , m_zmin ( right.zmin () )
  , m_zmax ( right.zmin () )
    //
  , m_b    ()
  , m_bz   ()
    //
  , m_fx ( right.nX() + 1 , 0.0 ) 
  , m_fy ( right.nY() + 1 , 0.0 ) 
  , m_fz ( right.nZ() + 1 , 0.0 ) 
{
  //
  typedef  Ostap::Math::Bernstein::Basic BB ;
  for ( unsigned short ix = 0 ; ix <= m_n  ; ++ix )
  { m_b .push_back ( Bernstein ( BB ( ix , m_n  ) , m_xmin , m_xmax ) ) ; }
  for ( unsigned short iz = 0 ; iz <= m_nz ; ++iz )
  { m_bz.push_back ( Bernstein ( BB ( iz , m_nz ) , m_zmin , m_zmax ) ) ; }
  //
}
// ============================================================================
// swap  two 3D-polynomials 
// ============================================================================
void Ostap::Math::Bernstein3DMix::swap
( Ostap::Math::Bernstein3DMix&  right ) 
{
  Ostap::Math::Parameters::swap ( right ) ;
  //
  std::swap ( m_n    , right.m_n     ) ;
  std::swap ( m_nz   , right.m_nz    ) ;
  std::swap ( m_xmin , right.m_xmin  ) ;
  std::swap ( m_xmax , right.m_xmax  ) ;
  std::swap ( m_zmin , right.m_zmin  ) ;
  std::swap ( m_zmax , right.m_zmax  ) ;
  std::swap ( m_b    , right.m_b     ) ;
  std::swap ( m_bz   , right.m_bz    ) ;
  std::swap ( m_fx   , right.m_fx    ) ;
  std::swap ( m_fy   , right.m_fy    ) ;
  std::swap ( m_fz   , right.m_fz    ) ;
}
// ============================================================================
// helper function to make calculations
// ============================================================================
double Ostap::Math::Bernstein3DMix::calculate
( const std::vector<double>& fx , 
  const std::vector<double>& fy , 
  const std::vector<double>& fz ) const 
{
  double       result = 0 ;
  for  ( unsigned short ix = 0 ; ix <= nX () ; ++ix )
  {
    for  ( unsigned short iy = 0 ; iy <= ix ; ++iy )
    { 
      for  ( unsigned short iz = 0 ; iz <= nZ () ; ++iz )
      { 
        double r = 0 ;
        if      ( ix == iy ) 
        { 
          r += fx[ix] * fy[iy] * fz[iz] ;
        }
        else 
        { 
          r += fx[ix] * fy[iy] * fz[iz] 
            +  fx[iy] * fy[ix] * fz[iz] ;  
        }
        //
        result += r * par( ix , iy , iz ) ;
      }
    }  
  }
  //
  const double scalex = ( nX () + 1 ) / ( xmax() - xmin() ) ;
  const double scaley = scalex ;
  const double scalez = ( nZ () + 1 ) / ( zmax() - zmin() ) ;
  //
  return result * scalex * scaley * scalez ;
}
// ============================================================================
// get the value
// ============================================================================
double Ostap::Math::Bernstein3DMix::evaluate
( const double x ,
  const double y , 
  const double z ) const
{
  /// the trivial cases
  if ( x < xmin () || x > xmax () ) { return 0.0        ; }
  if ( y < ymin () || y > ymax () ) { return 0.0        ; }
  if ( z < zmin () || z > zmax () ) { return 0.0        ; }
  //
  if      ( 0 == npars ()       ) { return 0.0        ; }
  else if ( 1 == npars ()       ) 
  { 
    const double scalex = ( nX () + 1 ) / ( xmax() - xmin() ) ;
    const double scaley = scalex ;
    const double scalez = ( nZ () + 1 ) / ( zmax() - zmin() ) ;
    //
    return m_pars [0] * scalex * scaley * scalez ;
  }
  ///
  // std::vector<double> fx ( nX () + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nX () ; ++i )
  { m_fx[i] = m_b [i] ( x )  ; }
  //
  // std::vector<double> fy ( nY () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nY () ; ++i )
  { m_fy[i] = m_b [i] ( y )  ; }
  //
  // std::vector<double> fz ( nZ () + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nZ () ; ++i )
  { m_fz[i] = m_bz[i] ( z )  ; }
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
/** get the integral over 3D-region
 *  \f[  x_{min} < x < x_{max}, y_{min}< y< y_{max} , z_{min} < z < z_{max}\f]
 */
// ============================================================================
double Ostap::Math::Bernstein3DMix::integral() const
{ 
  double       result = 0 ;
  for  ( unsigned short ix = 0 ; ix <= nX () ; ++ix )
  {
    for  ( unsigned short iy = 0 ; iy <= ix ; ++iy )
    { 
      for  ( unsigned short iz = 0 ; iz <= nZ () ; ++iz )
      { 
        const unsigned short r =  ix == iy ? 1 : 2 ;
        result += r  *  par( ix , iy , iz ) ;
      }
    }  
  }
  //
  return result ;
}
// ============================================================================
/* get the integral over 3D-region
 *  \f[ \int_{x_{low}}^{x_{high}}
 *      \int_{y_{low}}^{y_{high}} 
 *      \int_{z_{low}}^{z_{high}} 
 *\mathcal{B}(x,y,z) \mathrm{d}x\mathrm{d}y\mathrm{d}z\f]
 *  @param xlow  low  edge in x
 *  @param xhigh high edge in x
 *  @param ylow  low  edge in y
 *  @param yhigh high edge in y
 *  @param zlow  low  edge in z
 *  @param zhigh high edge in z
 */
// ============================================================================
double Ostap::Math::Bernstein3DMix::integral
( const double xlow , const double xhigh ,
  const double ylow , const double yhigh ,
  const double zlow , const double zhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ||
            s_equal ( ylow , yhigh ) ||
            s_equal ( zlow , zhigh ) ) { return 0 ; }
  //
  else if ( s_equal ( xlow  , xmin () ) &&
            s_equal ( xhigh , xmax () ) &&
            s_equal ( ylow  , ymin () ) &&
            s_equal ( yhigh , ymax () ) &&  
            s_equal ( zlow  , zmin () ) &&
            s_equal ( zhigh , zmax () ) )  { return integral () ; }
  //
  else if ( xlow  > xhigh ) 
  { return -1*integral ( xhigh , xlow  , ylow  , yhigh , zlow  , zhigh ) ; }
  else if ( ylow  > yhigh ) 
  { return -1*integral ( xlow  , xhigh , yhigh , ylow  , zlow  , zhigh ) ; }
  else if ( zlow  > zhigh ) 
  { return -1*integral ( xlow  , xhigh , ylow  , yhigh , zhigh , zlow  ) ; }
  //
  else if ( xhigh <  xmin () || xlow >  xmax() ) { return 0 ; }
  else if ( yhigh <  ymin () || ylow >  ymax() ) { return 0 ; }
  else if ( zhigh <  zmin () || zlow >  zmax() ) { return 0 ; }
  //
  const double  x_low  = std::max ( xmin() , xlow  ) ;
  const double  x_high = std::min ( xmax() , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  const double  y_low  = std::max ( ymin() , ylow  ) ;
  const double  y_high = std::min ( ymax() , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  const double  z_low  = std::max ( zmin() , zlow  ) ;
  const double  z_high = std::min ( zmax() , zhigh ) ;
  if ( z_low >= z_high ) { return 0 ; }
  //
  // std::vector<double> fx ( nX () + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nX () ; ++i )
  { m_fx[i] = m_b [i].integral ( x_low , x_high ) ; }
  //
  // std::vector<double> fy ( nY () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nY () ; ++i )
  { m_fy[i] = m_b [i].integral ( y_low , y_high ) ; }
  //
  // std::vector<double> fz ( nZ () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nZ () ; ++i )
  { m_fz[i] = m_bz[i].integral ( z_low , z_high ) ; }
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
/*  integral over x-dimension
 *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y,z) \mathrm{d}y\f]
 *  @param y     variable
 *  @param z     variable
 *  @param ylow  low  edge in y
 *  @param yhigh high edge in y
 */
// ============================================================================
double Ostap::Math::Bernstein3DMix::integrateX
( const double y    ,
  const double z    ,
  const double xlow , const double xhigh ) const
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; }
  else if ( xlow  > xhigh  ) { return -1*integrateX ( y , z , xhigh , xlow ) ; }
  else if ( xhigh <= xmin () || xlow >= xmax() ) { return 0 ; }
  else if ( y     <  ymin () || y    >  ymax() ) { return 0 ; }
  else if ( z     <  zmin () || z    >  zmax() ) { return 0 ; }
  else if ( s_equal ( xlow  , xmin () ) &&
            s_equal ( xhigh , xmax () )         ) { return integrateX ( y ,  z ) ; }
  //
  const double  x_low  = std::max ( xmin() , xlow  ) ;
  const double  x_high = std::min ( xmax() , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  //  std::vector<double> fx ( nX () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nX ()  ; ++i )
  { m_fx[i] = m_b [i].integral ( x_low , x_high ) ; }
  //
  // std::vector<double> fy ( nY () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nY ()  ; ++i )
  { m_fy[i] = m_b [i] ( y ) ; }
  //
  // std::vector<double> fz ( nZ () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nZ () ; ++i )
  { m_fz[i] = m_bz[i] ( z ) ; }
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
/*  integral over z-dimension
 *  \f[ \int_{z_{low}}^{z_{high}} \mathcal{B}(x,y,z) \mathrm{d}z\f]
 *  @param x     variable
 *  @param y     variable
 *  @param zlow  low  edge in z
 *  @param zhigh high edge in z
 */
// ============================================================================
double Ostap::Math::Bernstein3DMix::integrateZ
( const double x    ,
  const double y    ,
  const double zlow , const double zhigh ) const
{
  if      ( s_equal ( zlow  , zhigh ) ) { return 0 ; }
  else if ( zlow  >  zhigh ) { return -1*integrateZ ( x , y , zhigh , zlow  ) ; }
  else if ( x     <  xmin () || x    >  xmax() ) { return 0 ; }
  else if ( y     <  ymin () || y    >  ymax() ) { return 0 ; }
  else if ( zhigh <= zmin () || zlow >= zmax() ) { return 0 ; }
  else if ( s_equal ( zlow  , zmin() ) &&
            s_equal ( zhigh , zmax() )         ) { return integrateZ ( x , y ) ; }
  //
  const double  z_low  = std::max ( zmin() , zlow  ) ;
  const double  z_high = std::min ( zmax() , zhigh ) ;
  if ( z_low >= z_high ) { return 0 ; }
  //
  // std::vector<double> fx ( nX () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nX () ; ++i )
  { m_fx[i] = m_b [i] ( x ) ; }
  //
  // std::vector<double> fy ( nY () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nY () ; ++i )
  { m_fy[i] = m_b [i] ( y ) ; }
  //
  // std::vector<double> fz ( nZ () + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nZ () ; ++i ) 
  { m_fz[i] = m_bz[i].integral ( z_low , z_high ) ; }
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
double
Ostap::Math::Bernstein3DMix::integrateX
( const double y , 
  const double z ) const
{
  if      ( y < ymin () || y > ymax() ) { return 0 ; }
  else if ( z < zmin () || z > zmax() ) { return 0 ; }
  //
  // const std::vector<double> fx ( nX () + 1 , ( xmax() - xmin () ) / ( m_n + 1 ) ) ;
  std::fill ( m_fx.begin() , m_fx.end() , ( xmax() - xmin () ) / ( m_n + 1 ) ) ;
  //
  // std::vector<double> fy ( nY () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nY () ; ++i )
  { m_fy[i] = m_b [i] ( y ) ; }
  //
  // std::vector<double> fz ( nY () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nY ()  ; ++i )
  { m_fz[i] = m_bz[i] ( z ) ; }
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
double Ostap::Math::Bernstein3DMix::integrateZ ( const double x , 
                                                 const double y ) const
{
  if      ( x < xmin () || x > xmax() ) { return 0 ; }
  else if ( y < ymin () || y > ymax() ) { return 0 ; }
  //
  // std::vector<double> fx ( nX () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nX () ; ++i )
  { m_fx[i] = m_b [i] ( x ) ; }
  //
  // std::vector<double> fy ( nY() + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nY () ; ++i )
  { m_fy[i] = m_b [i] ( y ) ; }
  //
  // const std::vector<double> fz ( nZ () + 1  , ( zmax() - zmin () ) / ( nZ () + 1  ) ) ;
  std::fill ( m_fz.begin() , m_fz.end() , ( zmax() - zmin () ) / ( nZ () + 1 ) ) ;
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
/*  integral over x&y-dimensions
 *  \f[ \int_{x_{low}}^{x_{high}}
 *      \int_{y_{low}}^{y_{high}} \mathcal{B}(x,y,z) \mathrm{d}x\mathrm{d}y\f]
 *  @param z     variable
 *  @param xlow  low  edge in x
 *  @param xhigh high edge in x
 *  @param ylow  low  edge in y
 *  @param yhigh high edge in y
 */
// ============================================================================
double Ostap::Math::Bernstein3DMix::integrateXY
( const double z    ,                          
  const double xlow , const double xhigh ,
  const double ylow , const double yhigh ) const 
{
  if      ( s_equal ( xlow  , xhigh ) ) { return 0 ; }
  else if ( s_equal ( ylow  , yhigh ) ) { return 0 ; }
  else if ( xlow  >  xhigh ) { return -1*integrateXY ( z , xhigh , xlow  , ylow  , yhigh ) ; }
  else if ( ylow  >  yhigh ) { return -1*integrateXY ( z , xlow  , xhigh , yhigh , ylow  ) ; }
  else if ( z     <  zmin () || z    >  zmax() ) { return 0 ; }
  else if ( s_equal ( xlow  , xmin () ) &&
            s_equal ( xhigh , xmax () ) &&
            s_equal ( ylow  , ymin () ) &&
            s_equal ( yhigh , ymax () ) ) { return integrateXY ( z ) ; }
  //  //
  const double  x_low  = std::max ( xmin() , xlow  ) ;
  const double  x_high = std::min ( xmax() , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  const double  y_low  = std::max ( ymin() , ylow  ) ;
  const double  y_high = std::min ( ymax() , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  // std::vector<double> fx ( nX () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nX () ; ++i )
  { m_fx[i] = m_b [i].integral ( x_low , x_high ) ; }
  //
  // std::vector<double> fy ( nY () + 1  , 0 ) ;
  for ( unsigned short i = 0 ; i <= nY ()  ; ++i )
  { m_fy[i] = m_b [i].integral ( y_low , y_high ) ; }
  //
  // std::vector<double> fz ( nZ () + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nZ () ; ++i )
  { m_fz[i] = m_bz [i] ( z ) ; }
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
/*  integral over x&z-dimensions
 *  \f[ \int_{x_{low}}^{x_{high}}
 *      \int_{z_{low}}^{z_{high}} \mathcal{B}(x,y,z) \mathrm{d}x\mathrm{d}z\f]
 *  @param x     variable
 *  @param xlow  low  edge in x
 *  @param xhigh high edge in x
 *  @param zlow  low  edge in z
 *  @param zhigh high edge in z
 */
// ============================================================================
double Ostap::Math::Bernstein3DMix::integrateXZ
( const double y    ,                          
  const double xlow , const double xhigh ,
  const double zlow , const double zhigh ) const 
{
  if      ( s_equal ( xlow  , xhigh ) ) { return 0 ; }
  else if ( s_equal ( zlow  , zhigh ) ) { return 0 ; }
  else if ( xlow  >  xhigh ) { return -1*integrateXZ ( y , xhigh , xlow  , zlow  , zhigh ) ; }
  else if ( zlow  >  zhigh ) { return -1*integrateXZ ( y , xlow  , xhigh , zlow  , zhigh ) ; }
  else if ( y     <  ymin () || y    >  ymax() ) { return 0 ; }
  else if ( s_equal ( xlow  , xmin () ) &&
            s_equal ( xhigh , xmax () ) &&
            s_equal ( zlow  , zmin () ) &&
            s_equal ( zhigh , zmax () ) ) { return integrateXZ ( y ) ; }
  //  //
  const double  x_low  = std::max ( xmin() , xlow  ) ;
  const double  x_high = std::min ( xmax() , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  const double  z_low  = std::max ( zmin() , zlow  ) ;
  const double  z_high = std::min ( zmax() , zhigh ) ;
  if ( z_low >= z_high ) { return 0 ; }
  //
  // std::vector<double> fx ( nX () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nX () ; ++i )
  { m_fx[i] = m_b [i].integral ( x_low , x_high ) ; }
  //
  // std::vector<double> fy ( nY () + 1  , 0 ) ;
  for ( unsigned short i = 0 ; i <= nY ()  ; ++i )
  { m_fy[i] = m_b [i] ( y ) ; }
  //
  // std::vector<double> fz ( nZ () + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= nZ () ; ++i )
  { m_fz[i] = m_b [i].integral ( z_low , z_high ) ; }
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
/*  integral over x&y-dimensions
 *  \f[ \int_{x_{min}}^{x_{max}}
 *      \int_{y_{min}}^{y_{max}} \mathcal{B}(x,y,z) \mathrm{d}x\mathrm{d}y\f]
 *  @param z     variable
 */
// ============================================================================
double
Ostap::Math::Bernstein3DMix::integrateXY
( const double z    ) const 
{
  if ( z < zmin () || z > zmax() ) { return 0 ; }
  //
  // const std::vector<double> fx ( nX () + 1 , ( xmax() - xmin () ) / ( nX () + 1 ) ) ;
  // const std::vector<double> fy ( nY () + 1 , ( ymax() - ymin () ) / ( nY () + 1 ) ) ;
  std::fill ( m_fx.begin() , m_fx.end() , ( xmax() - xmin () ) / ( nX() + 1 ) ) ;
  std::fill ( m_fy.begin() , m_fy.end() , ( ymax() - ymin () ) / ( nY() + 1 ) ) ;
  //
  // std::vector<double> fz ( nZ () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nZ () ; ++i )
  { m_fz[i] = m_bz [i] ( z ) ; }
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
/*  integral over x&z-dimensions
 *  \f[ \int_{x_{min}}^{x_{max}}
 *      \int_{z_{min}}^{z_{max}} \mathcal{B}(x,y,z) \mathrm{d}x\mathrm{d}z\f]
 *  @param y     variable
 */
// ============================================================================
double Ostap::Math::Bernstein3DMix::integrateXZ ( const double y    ) const 
{
  if ( y < ymin () || y > ymax() ) { return 0 ; }
  //
  // const std::vector<double> fx ( nX () + 1 , ( xmax() - xmin () ) / ( nX () + 1 ) ) ;
  std::fill ( m_fx.begin() , m_fx.end() , ( xmax() - xmin () ) / ( nX() + 1 ) ) ;
  //
  // std::vector<double> fy ( nY () + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= nY () ; ++i )
  { m_fy[i] = m_b  [i] ( y ) ; }
  //
  // const std::vector<double> fz ( nZ () + 1 , ( zmax() - zmin () ) / ( nZ () + 1 ) ) ;
  std::fill ( m_fz.begin() , m_fz.end() , ( zmax() - zmin () ) / ( nZ() + 1 ) ) ;
  //
  return calculate ( m_fx , m_fy , m_fz ) ;
}
// ============================================================================
Ostap::Math::Bernstein3DMix&
Ostap::Math::Bernstein3DMix::operator+=( const double a )
{
  if   ( s_zero ( a ) ) { return *this ; }
  Ostap::Math::shift ( m_pars , a ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Bernstein3DMix&
Ostap::Math::Bernstein3DMix::operator*=( const double a )
{
  if      ( s_equal ( a , 1 ) ) { return *this ; }
  else if ( s_zero  ( a     ) ) { std::fill ( m_pars.begin() , m_pars.end() , 0 ) ; }
  Ostap::Math::scale ( m_pars , a ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Bernstein3DMix&
Ostap::Math::Bernstein3DMix::operator-=( const double a )
{
  if ( s_zero ( a ) ) { return *this ; }
  Ostap::Math::shift ( m_pars , -a ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Bernstein3DMix&
Ostap::Math::Bernstein3DMix::operator/=( const double a )
{
  if   ( s_equal ( a , 1 ) ) { return *this ; }
  Ostap::Math::scale ( m_pars , 1/a ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Bernstein3DMix
Ostap::Math::Bernstein3DMix::operator-() const
{
  Bernstein3DMix b ( *this ) ;
  Ostap::Math::negate ( b.m_pars ) ;
  return b ;
}
// ============================================================================
// Sum of Bernstein polynomial and a constant
// ============================================================================
Ostap::Math::Bernstein3DMix
Ostap::Math::Bernstein3DMix::__add__   ( const double value ) const
{ return (*this) + value ; }
// ============================================================================
// Sum of Bernstein polynomial and a constant
// ============================================================================
Ostap::Math::Bernstein3DMix
Ostap::Math::Bernstein3DMix::__radd__  ( const double value ) const
{ return value + (*this) ; }
// ============================================================================
// Product of Bernstein polynomial and a constant
// ============================================================================
Ostap::Math::Bernstein3DMix
Ostap::Math::Bernstein3DMix::__mul__   ( const double value ) const
{ return (*this) * value ; }
// ============================================================================
// Product of Bernstein polynomial and a constant
// ============================================================================
Ostap::Math::Bernstein3DMix
Ostap::Math::Bernstein3DMix::__rmul__  ( const double value ) const
{ return value * (*this) ; }
// ============================================================================
// Subtract a constant from Benrstein polynomial
// ============================================================================
Ostap::Math::Bernstein3DMix
Ostap::Math::Bernstein3DMix::__sub__  ( const double value ) const
{ return (*this) - value ; }
// ============================================================================
// Subtract Bernstein polynomial from a constant
// ============================================================================
Ostap::Math::Bernstein3DMix
Ostap::Math::Bernstein3DMix::__rsub__ ( const double value ) const
{ return value - (*this) ; }
// ============================================================================
// Divide Benrstein polynomial by a constant
// ============================================================================
Ostap::Math::Bernstein3DMix
Ostap::Math::Bernstein3DMix:: __truediv__   ( const double value ) const
{ return (*this) / value ; }
// ============================================================================
// Negate Bernstein polynomial
// ============================================================================
Ostap::Math::Bernstein3DMix
Ostap::Math::Bernstein3DMix::__neg__ ()  const
{ return -(*this); }
// ============================================================================
// get the tag value 
// ============================================================================
std::size_t Ostap::Math::Bernstein3DMix::tag () const  // get the tag value 
{
  return Ostap::Utils::hash_combiner 
    ( Ostap::Utils::hash_range ( m_pars ) , 
      m_n    , m_nz   ,
      m_xmin , m_xmax ,
      m_zmin , m_zmax ) ; }
// ============================================================================


// ============================================================================
// 3D-POSITIVE 
// ============================================================================

// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::Positive3D::Positive3D 
( const unsigned short       Nx    ,
  const unsigned short       Ny    ,
  const unsigned short       Nz    ,
  const double               xmin  ,
  const double               xmax  ,
  const double               ymin  ,
  const double               ymax  ,
  const double               zmin  ,
  const double               zmax  )
  : m_bernstein (   Nx , Ny , Nz,  xmin , xmax , ymin , ymax , zmin , zmax ) 
  , m_sphere    ( ( Nx + 1 ) * ( Ny + 1 ) * ( Nz + 1 ) - 1 )
{
  updateBernstein () ;
}
// ============================================================================
// constructor from parameters 
// ============================================================================
Ostap::Math::Positive3D::Positive3D 
( const std::vector<double>& pars  , 
  const unsigned short       Nx    ,
  const unsigned short       Ny    ,
  const unsigned short       Nz    ,
  const double               xmin  ,
  const double               xmax  ,
  const double               ymin  ,
  const double               ymax  ,
  const double               zmin  ,
  const double               zmax  )
  : m_bernstein (   Nx , Ny , Nz,  xmin , xmax , ymin , ymax , zmin , zmax ) 
  , m_sphere    ( ( Nx + 1 ) * ( Ny + 1 ) * ( Nz + 1 ) - 1 )
{
  m_sphere.setPars ( pars ) ;
  updateBernstein () ;
}
// ============================================================================
// swap  two 2D-polynomials 
// ============================================================================
void Ostap::Math::Positive3D::swap ( Ostap::Math::Positive3D&  right ) 
{
  Ostap::Math::swap ( m_bernstein , right.m_bernstein ) ;
  Ostap::Math::swap ( m_sphere    , right.m_sphere    ) ;  
}
// ============================================================================
// set k-parameter
// ============================================================================
bool Ostap::Math::Positive3D::setPar 
( const unsigned int k     , 
  const double       value ,
  const bool         force )
{
  //
  const bool update = m_sphere.setPhase ( k , value , force ) ;
  if ( !update ) { return false ; }   // no actual change 
  //
  return updateBernstein () ;
}
// =============================================================================
// update bernstein coefficients
// =============================================================================
bool Ostap::Math::Positive3D::updateBernstein ()
{
  //
  bool update = false ;
  for ( unsigned int ix = 0 ; ix < m_sphere.nX() ; ++ix ) 
  { 
    const bool updated = m_bernstein.setPar ( ix , m_sphere.x2 ( ix ) ) ;
    update = updated || update ;  
  }
  //
  if ( update ) { m_bernstein /= m_bernstein.integral() ; }  
  //
  return update ;
}
// ============================================================================
/*  get the integral over 3D-region           
 *  \f[ \int_{x_{min}}^{x_{max}}
 *      \int_{y_{min}}^{y_{max}} 
 *      \int_{z_{min}}^{z_{max}} 
 *        \mathcal{B}(x,y,z) \mathrm{d}x\mathrm{d}y\mathrm{d}z\f] 
 */
// ============================================================================
double  Ostap::Math::Positive3D::integral   () const { return 1 ; }
// ============================================================================
/* get the integral over 3D-region 
 *  \f[ \int_{x_{low}}^{x_{high}}
 *      \int_{y_{low}}^{y_{high}} 
 *      \int_{z_{low}}^{z_{high}} 
 *      \mathcal{B}(x,y,z) \mathrm{d}x\mathrm{d}y\mathrm{d}z\f] 
 *  @param xlow  low  edge in x 
 *  @param xhigh high edge in x 
 *  @param ylow  low  edge in y 
 *  @param yhigh high edge in y 
 *  @param zlow  low  edge in z 
 *  @param zhigh high edge in z 
 */
// ============================================================================
double Ostap::Math::Positive3D::integral   
( const double xlow , const double xhigh , 
  const double ylow , const double yhigh ,
  const double zlow , const double zhigh ) const 
{ 
  return
    s_equal ( xlow  , xmin() ) && 
    s_equal ( xhigh , xmax() ) && 
    s_equal ( ylow  , ymin() ) && 
    s_equal ( yhigh , ymax() ) && 
    s_equal ( zlow  , zmin() ) && 
    s_equal ( zhigh , zmax() )  ? 1.0 : 
    m_bernstein.integral ( xlow , xhigh , ylow , yhigh , zlow , zhigh ) ; 
}
// ============================================================================

// ============================================================================
// 3D-SYMMETRIC POSITIVE 
// ============================================================================

// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::Positive3DSym::Positive3DSym 
( const unsigned short       N     ,
  const double               xmin  ,
  const double               xmax  )
  : m_bernstein (   N ,  xmin , xmax ) 
  , m_sphere    ( ( N + 1 ) * ( N + 2 ) * ( N + 3 ) / 6 - 1 ) ///  ????? 
{
  updateBernstein () ;
}
// ============================================================================
// constructor from parameters 
// ============================================================================
Ostap::Math::Positive3DSym::Positive3DSym 
( const std::vector<double>& pars  ,
  const unsigned short       N     ,
  const double               xmin  ,
  const double               xmax  )
  : m_bernstein (   N ,  xmin , xmax ) 
  , m_sphere    ( ( N + 1 ) * ( N + 2 ) * ( N + 3 ) / 6 - 1 ) ///  ????? 
{
  m_sphere.setPars ( pars ) ;
  updateBernstein () ;
}
// ============================================================================
// swap  two 2D-polynomials 
// ============================================================================
void Ostap::Math::Positive3DSym::swap ( Ostap::Math::Positive3DSym&  right ) 
{
  Ostap::Math::swap ( m_bernstein , right.m_bernstein ) ;
  Ostap::Math::swap ( m_sphere    , right.m_sphere    ) ;  
}
// ============================================================================
// set k-parameter
// ============================================================================
bool Ostap::Math::Positive3DSym::setPar 
( const unsigned int k     , 
  const double       value ,
  const bool         force )
{
  //
  const bool update = m_sphere.setPhase ( k , value , force ) ;
  if ( !update ) { return false ; }   // no actual change 
  //
  return updateBernstein () ;
}
// =============================================================================
// update bernstein coefficients
// =============================================================================
bool Ostap::Math::Positive3DSym::updateBernstein ()
{
  //
  bool update = false ;
  for ( unsigned int ix = 0 ; ix < m_sphere.nX() ; ++ix ) 
  { 
    const bool updated = m_bernstein.setPar ( ix , m_sphere.x2 ( ix ) ) ;
    update = updated || update ;  
  }
  //
  if ( update ) { m_bernstein /= m_bernstein.integral() ; }
  //
  return update ;
}
// ============================================================================
/*  get the integral over 3D-region           
 *  \f[ \int_{x_{min}}^{x_{max}}
 *      \int_{y_{min}}^{y_{max}} 
 *      \int_{z_{min}}^{z_{max}} 
 *        \mathcal{B}(x,y,z) \mathrm{d}x\mathrm{d}y\mathrm{d}z\f] 
 */
// ============================================================================
double  Ostap::Math::Positive3DSym::integral   () const { return 1 ; }
// ============================================================================
/* get the integral over 3D-region 
 *  \f[ \int_{x_{low}}^{x_{high}}
 *      \int_{y_{low}}^{y_{high}} 
 *      \int_{z_{low}}^{z_{high}} 
 *      \mathcal{B}(x,y,z) \mathrm{d}x\mathrm{d}y\mathrm{d}z\f] 
 *  @param xlow  low  edge in x 
 *  @param xhigh high edge in x 
 *  @param ylow  low  edge in y 
 *  @param yhigh high edge in y 
 *  @param zlow  low  edge in z 
 *  @param zhigh high edge in z 
 */
// ============================================================================
double Ostap::Math::Positive3DSym::integral   
( const double xlow , const double xhigh , 
  const double ylow , const double yhigh ,
  const double zlow , const double zhigh ) const 
{ 
  return
    s_equal ( xlow  , xmin() ) && 
    s_equal ( xhigh , xmax() ) && 
    s_equal ( ylow  , ymin() ) && 
    s_equal ( yhigh , ymax() ) && 
    s_equal ( zlow  , zmin() ) && 
    s_equal ( zhigh , zmax() )  ? 1.0 : 
    m_bernstein.integral ( xlow , xhigh , ylow , yhigh , zlow , zhigh ) ; 
}
// ============================================================================


// ============================================================================
// 3D-POSITIVE with X<-->Y  SYMMETRY
// ============================================================================

// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::Positive3DMix::Positive3DMix
( const unsigned short       N     , 
  const unsigned short       Nz    ,
  const double               xmin  ,
  const double               xmax  ,
  const double               zmin  ,
  const double               zmax  )
  : m_bernstein (   N , Nz ,  xmin , xmax , zmin , zmax ) 
  , m_sphere    ( ( N + 1 ) * ( N + 2 ) * ( Nz + 1 ) / 2  - 1 ) 
{
  updateBernstein () ;
}
// ============================================================================
// constructor from parameters  
// ============================================================================
Ostap::Math::Positive3DMix::Positive3DMix
( const std::vector<double>& pars  , 
  const unsigned short       N     , 
  const unsigned short       Nz    ,
  const double               xmin  ,
  const double               xmax  ,
  const double               zmin  ,
  const double               zmax  )
  : m_bernstein (   N , Nz ,  xmin , xmax , zmin , zmax ) 
  , m_sphere    ( ( N + 1 ) * ( N + 2 ) * ( Nz + 1 ) / 2  - 1 ) 
{
  m_sphere.setPars  ( pars ) ;
  updateBernstein () ;
}
// ============================================================================
// swap  two 2D-polynomials 
// ============================================================================
void Ostap::Math::Positive3DMix::swap ( Ostap::Math::Positive3DMix&  right ) 
{
  Ostap::Math::swap ( m_bernstein , right.m_bernstein ) ;
  Ostap::Math::swap ( m_sphere    , right.m_sphere    ) ;  
}
// ============================================================================
// set k-parameter
// ============================================================================
bool Ostap::Math::Positive3DMix::setPar 
( const unsigned int k     , 
  const double       value ,
  const bool         force )
{
  //
  const bool update = m_sphere.setPhase ( k , value , force ) ;
  if ( !update ) { return false ; }   // no actual change 
  //
  return updateBernstein () ;
}
// =============================================================================
// update bernstein coefficients
// =============================================================================
bool Ostap::Math::Positive3DMix::updateBernstein ()
{
  //
  bool update = false ;
  for ( unsigned int ix = 0 ; ix < m_sphere.nX() ; ++ix ) 
  { 
    const bool updated = m_bernstein.setPar ( ix , m_sphere.x2 ( ix ) ) ;
    update = updated || update ;  
  }
  //
  if ( update ) { m_bernstein /= m_bernstein.integral() ; }
  //
  return update ;
}
// ============================================================================
/*  get the integral over 3D-region           
 *  \f[ \int_{x_{min}}^{x_{max}}
 *      \int_{y_{min}}^{y_{max}} 
 *      \int_{z_{min}}^{z_{max}} 
 *        \mathcal{B}(x,y,z) \mathrm{d}x\mathrm{d}y\mathrm{d}z\f] 
 */
// ============================================================================
double  Ostap::Math::Positive3DMix::integral   () const { return 1 ; }
// ============================================================================
/* get the integral over 3D-region 
 *  \f[ \int_{x_{low}}^{x_{high}}
 *      \int_{y_{low}}^{y_{high}} 
 *      \int_{z_{low}}^{z_{high}} 
 *      \mathcal{B}(x,y,z) \mathrm{d}x\mathrm{d}y\mathrm{d}z\f] 
 *  @param xlow  low  edge in x 
 *  @param xhigh high edge in x 
 *  @param ylow  low  edge in y 
 *  @param yhigh high edge in y 
 *  @param zlow  low  edge in z 
 *  @param zhigh high edge in z 
 */
// ============================================================================
double Ostap::Math::Positive3DMix::integral   
( const double xlow , const double xhigh , 
  const double ylow , const double yhigh ,
  const double zlow , const double zhigh ) const 
{ 
  return
    s_equal ( xlow  , xmin() ) && 
    s_equal ( xhigh , xmax() ) && 
    s_equal ( ylow  , ymin() ) && 
    s_equal ( yhigh , ymax() ) && 
    s_equal ( zlow  , zmin() ) && 
    s_equal ( zhigh , zmax() )  ? 1.0 : 
    m_bernstein.integral ( xlow , xhigh , ylow , yhigh , zlow , zhigh ) ; 
}
// ============================================================================

// ============================================================================
//                                                                      The END
// ============================================================================
