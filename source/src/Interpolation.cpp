// $Id$
// ===========================================================================
// Include files
// ===========================================================================
// Ostap
// ===========================================================================
#include "Ostap/Interpolation.h"
// ===========================================================================
/** @file 
 *  Implementation file for interpolation functions 
 *  @see Ostap::Math::Interpolation
 *  @see LHCbMath/Interpolation.h
 *  
 *  @date 2016-07-23 
 *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
 *
 *                    $Revision$
 *  Last modification $Date$
 *                 by $author$
 */
// ============================================================================
/*  very simp;e lagrange interpolation 
 *  @param  xs INPUT sequence of abscissas 
 *  @param  ys INPUT sequence of values 
 *  @param  x      INPUT evaluate the polynomial in this point
 *  @return the value of Largange interpolation polynomial at point x
 *
 *  - If ys-vector is shorter than x-vector, it is assumed to be zero padded
 *  - If ys-vector is longer  than x-vector, extra avalues are ignored       
 *  - If xs-vector is empty, polynomial is zero
 *
 *  @warning it could be CPU inefficient
 *  @warning it should *NOT* be applied for very long sequence of points 
 *           (e.g. >20) due to bad numerical  stability
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2016-07-23
 */
// ============================================================================
double Ostap::Math::Interpolation::lagrange 
( const std::vector<double>& xs , 
  const std::vector<double>& ys , 
  const double               x  ) 
{
  return lagrange ( xs.begin () , 
                    xs.end   () , 
                    ys.begin () , 
                    ys.end   () , 
                    x           , 
                    0.0         , 
                    [] ( double x ) { return x ; } , 
                    [] ( double y ) { return y ; } ) ;
}
// ============================================================================
/*  very simple lagrange interpolation 
 *  @param  data INPUT sequence of (x,y)
 *  @param  x      INPUT evaluate the polynomial in this point
 *  @return the value of Largange interpolation polynomial at point x
 *
 *  - If data-vector is empty, polynomial is zero
 *
 *  @warning it could be CPU inefficient
 *  @warning it should *NOT* be applied for very long sequence of points 
 *           (e.g. >20) due to bad numerical  stability
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2016-07-23
 */
// ============================================================================
double Ostap::Math::Interpolation::lagrange 
( const Ostap::Math::Interpolation::DATA& data , 
  const double                            x  ) 
{
  typedef std::pair<double,double>  PP ;
  return lagrange ( data.begin () , 
                    data.end   () , 
                    data.begin () , 
                    data.end   () , 
                    x             , 
                    0.0           , 
                    [] ( const PP& i ) { return i.first  ; } , 
                    [] ( const PP& i ) { return i.second ; } ) ; 
}
// ============================================================================
/*  Simple lagrange interpolation 
 *  - it also evaluate the derivative wity respect to y_i 
 *  @param  xs INPUT sequence of abscissas 
 *  @param  ys INPUT sequence of values 
 *  @param  x  INPUT evaluate the polynomial in this point
 *  @param  it INPUT index of y_i, the derivative shodul be calculated.
 *  @return the value of Largange interpolation polynomial at point x
 *
 *  - If ys-vector is shorter than x-vector, it is assumed to be zero padded
 *  - If ys-vector is longer  than x-vector, extra avalues are ignored       
 *  - If xs-vector is empty, polynomial is zero
 *
 *  @warning it could be CPU inefficient
 *  @warning it should *NOT* be applied for very long sequence of points 
 *           (e.g. >20) due to bad numerical  stability
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2016-07-23
 */
// ============================================================================
std::pair<double,double>
Ostap::Math::Interpolation::lagrange2 
( const std::vector<double>& xs , 
  const std::vector<double>& ys , 
  const double               x  , 
  const unsigned int         iy ) 
{
  const double r = lagrange ( xs.begin () , 
                              xs.end   () , 
                              ys.begin () , 
                              ys.end   () , 
                              x           , 
                              0.0         , 
                              [] ( double x )  { return x ; } , 
                              [] ( double y )  { return y ; } ) ;
  //
  if  ( ys.size() <= iy || xs.size() <= iy ) { return std::make_pair ( r , 0. ) ; } // RETURN
  //
  double  d = 1 ;
  const unsigned int nx = xs.size() ;
  const double       xi = xs[iy]    ;
  for ( unsigned int jx = 0  ; jx < nx ; ++jx ) 
  {
    if ( jx == iy ) { continue ; }
    const double xj = xs[jx] ;
    d *= ( x - xj ) / ( xi - xj ) ;
  }
  return std::make_pair (  r, d ) ;
}
// ============================================================================
/*  Simple lagrange interpolation 
 *  - it also evaluate the derivative wity respect to y_i 
 *  @param  data INPUT sequence of (x,y)
 *  @param  x      INPUT evaluate the polynomial in this point
 *  @param  it INPUT index of y_i, the derivative shodul be calculated.
 *  @return the value of Largange interpolation polynomial at point x
 *
 *  - If data-vector is empty, polynomial is zero
 *
 *  @warning it could be CPU inefficient
 *  @warning it should *NOT* be applied for very long sequence of points 
 *           (e.g. >20) due to bad numerical  stability
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2016-07-23
 */
// ============================================================================
std::pair<double,double>
Ostap::Math::Interpolation::lagrange2 
( const Ostap::Math::Interpolation::DATA& data ,
  const double                            x    , 
  const unsigned int                      iy   ) 
{
  typedef std::pair<double,double>  PP ;
  const double r = lagrange ( data.begin () , 
                              data.end   () , 
                              data.begin () , 
                              data.end   () , 
                              x             , 
                              0.0           , 
                              [] ( const PP&i ) { return i.first  ; } , 
                              [] ( const PP&i ) { return i.second ; } ) ; 
  //
  if  ( data.size() <= iy ) { return std::make_pair ( r , 0. ) ; } // RETURN
  //
  double  d = 1 ;
  const unsigned int nx = data.size()    ;
  const double       xi = data[iy].first ;
  for ( unsigned int jx = 0  ; jx < nx ; ++jx ) 
  {
    if ( jx == iy ) { continue ; }
    const double xj = data[jx].first ;
    d *= ( x - xj ) / ( xi - xj ) ;
  }
  return std::make_pair (  r, d ) ;
}
// ============================================================================
/** very simple Neville's interpolation 
 *  @param  xs INPUT sequence of abscissas 
 *  @param  ys INPUT sequence of values 
 *  @param  x      INPUT evaluate the polynomial in this point
 *  @return the value of Largange interpolation polynomial at point x
 *
 *  - If ys-vector is shorter than x-vector, it is assumed to be zero padded
 *  - If ys-vector is longer  than x-vector, extra avalues are ignored       
 *  - If xs-vector is empty, polynomial is zero
 *
 *  @warning it could be CPU inefficient
 *  @warning it should *NOT* be applied for very long sequence of points 
 *           (e.g. >20) due to bad numerical  stability
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2016-07-23
 */
// ============================================================================
double Ostap::Math::Interpolation::neville 
( const std::vector<double>& xs , 
  const std::vector<double>& ys , 
  const double               x  ) 
{
  return neville ( xs.begin() , xs.end  () ,  
                   ys.begin() , ys.end  () , 
                   x          , 
                   [] ( double x ) { return x ; } , 
                   [] ( double y ) { return y ; } ) ;
}                   
// ============================================================================
/*  very simple Neville interpolation 
 *  @param  data INPUT sequence of (x,y)
 *  @param  x      INPUT evaluate the polynomial in this point
 *  @return the value of Largange interpolation polynomial at point x
 *
 *  - If data-vector is empty, polynomial is zero
 *
 *  @warning it could be CPU inefficient
 *  @warning it should *NOT* be applied for very long sequence of points 
 *           (e.g. >20) due to bad numerical  stability
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2016-07-23
 */
// ============================================================================
double Ostap::Math::Interpolation::neville 
( const Ostap::Math::Interpolation::DATA& data , const double x  ) 
{
  return neville ( data.begin () , 
                   data.end   () , 
                   data.begin () , 
                   data.end   () , 
                   x             , 
                   [] ( const std::pair<double,double>& i ) { return i.first  ; } , 
                   [] ( const std::pair<double,double>& i ) { return i.second ; } );
}
// ============================================================================
/** very simple Neville's interpolation 
 *  -  it evaluates the polynomial and the derivative    
 *  @param  xs INPUT sequence of abscissas 
 *  @param  ys INPUT sequence of values 
 *  @param  x      INPUT evaluate the polynomial in this point
 *  @return the value of Largange interpolation polynomial at point x
 *
 *  - If ys-vector is shorter than x-vector, it is assumed to be zero padded
 *  - If ys-vector is longer  than x-vector, extra avalues are ignored       
 *  - If xs-vector is empty, polynomial is zero
 *
 *  @warning it could be CPU inefficient
 *  @warning it should *NOT* be applied for very long sequence of points 
 *           (e.g. >20) due to bad numerical  stability
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2016-07-23
 */
// ============================================================================
std::pair<double,double> 
Ostap::Math::Interpolation::neville2 
( const std::vector<double>& xs , 
  const std::vector<double>& ys , 
  const double               x  ) 
{
  return neville2 ( xs.begin() , xs.end  () ,  
                    ys.begin() , ys.end  () , 
                    x          , 
                    [] ( double x ) { return x ; } , 
                    [] ( double y ) { return y ; } ) ;
}                   
// ============================================================================
/*  very simple Neville interpolation 
 *  -  it evaluates the polynomial and the derivative    
 *  @param  data INPUT sequence of (x,y)
 *  @param  x      INPUT evaluate the polynomial in this point
 *  @return the value of Largange interpolation polynomial at point x
 *
 *  - If data-vector is empty, polynomial is zero
 *
 *  @warning it could be CPU inefficient
 *  @warning it should *NOT* be applied for very long sequence of points 
 *           (e.g. >20) due to bad numerical  stability
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2016-07-23
 */
// ============================================================================
std::pair<double,double>
Ostap::Math::Interpolation::neville2 
( const Ostap::Math::Interpolation::DATA& data , const double x  ) 
{
  return neville2 ( data.begin () , 
                    data.end   () , 
                    data.begin () , 
                    data.end   () , 
                    x             , 
                    [] ( const std::pair<double,double>& i ) { return i.first  ; } , 
                    [] ( const std::pair<double,double>& i ) { return i.second ; } );
}
// ============================================================================



// ============================================================================
//  The END 
// ============================================================================

