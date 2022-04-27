// ===========================================================================
// Include files
// ===========================================================================
// STD/STL
// ===========================================================================
#include <algorithm>
// ===========================================================================
// Ostap
// ===========================================================================
#include "Ostap/Math.h"
#include "Ostap/Interpolation.h"
#include "Ostap/Interpolants.h"
#include "Ostap/Differences.h"
#include "Ostap/Choose.h"
// ===========================================================================
// Local 
// ===========================================================================
#include "local_math.h"
// ===========================================================================
/** @file 
 *  Implementation file for interpolation functions 
 *  @see Ostap::Math::Interpolation
 *  @see Ostap/Interpolation.h
 *  
 *  @date 2016-07-23 
 *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
 */
// ============================================================================
namespace
{
  // ==========================================================================
  /// "numerically less" 
  const Ostap::Math::NumLess<double>  s_num_less  {} ;
  /// "numerically equal" 
  const Ostap::Math::Equal_To<double> s_num_equal {} ;
  // ==========================================================================
}
// ============================================================================
/*  create the abscissas from vector of abscissas 
 *  @param x       input vector of abscissas (to be sorted internally)
 *  @param sorted  indicate if input data is already sorted 
 */
// ============================================================================
Ostap::Math::Interpolation::Abscissas::Abscissas 
( const Ostap::Math::Interpolation::Abscissas::Data& x      , 
  const bool                                         sorted ) 
  : Ostap::Math::Interpolation::Abscissas  ( x.begin() , x.end() , sorted )
{}
// ============================================================================
/* sort abscissas and eliminate the duplicates  
 * @return number of removed duplicates  
 */
// ============================================================================
unsigned int Ostap::Math::Interpolation::Abscissas::sort ( ) 
{
  const unsigned int l = m_x.size() ;
  /// sort it 
  std::sort ( m_x.begin() , m_x.end() , s_num_less ) ;
  /// remove duplicates if any 
  m_x.erase ( std::unique ( m_x.begin () , m_x.end () , s_num_equal ) , m_x.end () ) ;
  /// get the result 
  return l - m_x.size() ;  
}
// ============================================================================
// make a slice for the given abscissas 
// ============================================================================
Ostap::Math::Interpolation::Abscissas 
Ostap::Math::Interpolation::Abscissas::slice 
( const int i , 
  const int j ) const 
{
  // treat the negative indices 
  const int ii = 
    i < 0     ? 0        : i ; // ignore  negative  i 
  const int jj = 
        j <  0    ? j + n () : 
    j >= n () ?     n () : j ; // adjust a bit not very negative j  
  //
  return 
        0 >= ii && jj >= n() ? (*this) : 
    0 <= ii && ii < jj   ? Abscissas ( m_x.begin() + ii , 
                                       m_x.begin() + jj , true ) :
    Abscissas () ;
}
// ============================================================================
/* special constructor for the given interpoaltion type 
 * @param n    number of interpolation points 
     * @param low  low edge of the interval 
 * @param high high edge of the interval 
 * @param t    interpolation type 
 * @see Ostap::Math::Interpolation::Abscissas::Type 
 */
// ============================================================================
Ostap::Math::Interpolation::Abscissas::Abscissas 
( const unsigned short                               n    ,
  const double                                       low  , 
  const double                                       high  , 
  const Ostap::Math::Interpolation::Abscissas::AType t    ) 
  : m_x     ( n , 0.0 )
  , m_atype ( Generic ) 
  , m_xmin  ( std::min ( low , high ) )
  , m_xmax  ( std::max ( low , high ) )
{  
  const long double mn = m_xmin ;
  const long double mx = m_xmax ;
  //
  if ( 0 == n ) { return ; }                               // QUIT 
  if ( 1 == n ) { m_x [0] = 0.5 * ( mn + mx ) ; return ; } // QUIT 
  //
  switch ( t ) 
  {
    // ========================================================================
  case Chebyshev1 :  // roots of Tn(x)
    // ========================================================================
    for  ( unsigned short i = 0 ; i < n ; ++i ) 
    {
      const long double a = 1.0L * ( 2 * ( n - i ) - 1 ) * M_PIl / ( 2 * n ) ;
      const long double x = std::cos ( a ) ;
      m_x [ i ] = 0.5 * ( ( 1 - x ) * mn + ( 1 + x ) * mx ) ;
    }
    m_atype = Chebyshev1 ;
    break ;
    // ========================================================================    
  case Chebyshev2 :  // extrema of Tn(x)
    // ========================================================================
    for  ( unsigned short i = 0 ; i < n ; ++i ) 
    {
      const long double a = 1.0L * ( ( n - i - 1 ) * M_PIl / ( n - 1 ) ) ;
      const long double x = std::cos ( a ) ;
      m_x [ i ] = 0.5 * ( ( 1 - x ) * mn + ( 1 + x ) * mx ) ;
    }
    m_atype = Chebyshev2 ;
    break ;
    // ========================================================================
  default :  // uniform 
    // ========================================================================
    { 
      const long double in = 1.0L / ( n - 1 ) ;
      for ( unsigned short i = 0 ; i < n ; ++i ) 
      { 
        const long double beta  = i * in       ;
        m_x [ i ] = ( 1 - beta ) * mn + beta * mx ;
      }
      m_atype = Uniform ;
      break ;
    }
  }
  //
}  
// ============================================================================
/* add one more abscissa in the list 
 * @param xnew value to be added 
 * @return -1 if point is NOT added or new index for added points  
 */
// ============================================================================
int Ostap::Math::Interpolation::Abscissas::add ( const double xnew ) 
{
  // where to insert new element?
  Data::iterator ifound = std::lower_bound 
    ( m_x.begin () , 
      m_x.end   () , 
      xnew         ,
      s_num_less   ) ;
  // all elements are less than new element
  if  ( m_x.end() == ifound ) 
  { 
    m_x.push_back ( xnew ) ;
    m_xmin  = m_x.front () ;
    m_xmax  = m_x.back  () ;
    m_atype = Generic ;                           // ATTENTION 
    return m_x.size() - 1 ;                       // RETURN
  } 
  // there exists element with the same value 
  else if ( !s_num_less ( xnew , *ifound ) ) { return - 1 ; } // RETURN
  // insert new element at the given position
  const int index = ifound - m_x.begin() ;
  m_x.insert ( ifound , xnew ) ;
  //
  m_atype = Generic ;                           // ATTENTION! 
  m_xmin = m_x.front () ;
  m_xmax = m_x.back  () ;
  return index ;
}
// ============================================================================
/* remove the point with the given index 
 * @param index poitn with the index 
 * @return true if point is really removed 
 */
// ============================================================================
bool  Ostap::Math::Interpolation::Abscissas::remove ( const unsigned short index ) 
{
  if ( index >= m_x.size() ) { return false ; }
  //
  m_atype = Generic ;                             // ATTENTION!
  m_x.erase ( m_x.begin() + index ) ;
  if ( !empty() ) 
  {
    m_xmin = m_x.front () ;
    m_xmax = m_x.back  () ;
  }
  else 
  {
    m_xmin = 0 ;
    m_xmax = 0 ;
  }
  //
  return true ;
}
// ============================================================================

// ============================================================================
// Interpolation table 
// ============================================================================
/* simple constructor from x&y-lists 
 *  @param x input vector of x 
 *  @param y input vector of y 
 *  @param sorted  indicate that x-list is already sorted 
 *  - if vector of y is longer  than vector x, extra values are ignored 
 *  - if vector of y is shorter than vector x, missing entries are assumed to be zero  
 *  @attention duplicated abscissas will be removed 
 */
// ============================================================================
Ostap::Math::Interpolation::Table::Table
( const Ostap::Math::Interpolation::Abscissas::Data& x      , 
  const Ostap::Math::Interpolation::Abscissas::Data& y      , 
  const bool                                         sorted ) 
  : Ostap::Math::Interpolation::Table ( x.begin () , 
                                        x.end   () , 
                                        y.begin () , 
                                        y.end   () , 
                                        sorted     )
{}
// ======================================================================
/*  the simplest constructor 
 *  @param table  input data 
 *  @param sorted indicate if data already  sorted and duplicates removed  
 */
// ============================================================================
Ostap::Math::Interpolation::Table::Table
( Ostap::Math::Interpolation::TABLE table  , 
  const bool                        sorted ) 
  : m_abscissas () 
  , m_values    () 
{
  /// sort the table and eliminate duplicates (if needed) 
  get_sorted ( sorted , table ) ;
  /// 4) fill abscissas 
  m_abscissas = Abscissas 
    ( table.begin () , 
      table.end   () ,  
      []( const auto& p )-> double { return p.first ; } , true ) ;
  /// 5) fill values 
  m_values.resize( table.size () ) ;
  std::transform 
    ( table.begin    () ,
      table.end      () , 
      m_values.begin () ,
      []( const auto& p )-> double { return p.second ; }       ) ; 
}
// ============================================================================
// sort table and remove duplicated abscissas 
// ============================================================================
void Ostap::Math::Interpolation::Table::get_sorted 
( const bool                         sorted , 
  Ostap::Math::Interpolation::TABLE& table  )
{ 
  // =========================================================================
  if ( !sorted ) 
  {
    std::sort ( table.begin () , 
                table.end   () , 
                [] ( auto a , auto b ) 
                { return s_num_less ( a.first , b.first ) ; } ) ;
  }
  /// remove duplicates if any 
  table.erase ( std::unique ( table.begin () , 
                              table.end   () , 
                              [] ( auto a , auto b ) 
                              { return s_num_equal ( a.first , b.first ) ; } ) , 
                table.end ()  ) ;
  // =========================================================================
}
// ============================================================================
/*  add the point (x,y) into interpolation table 
 *  @param x abscissas of the point to be added 
 *  @param y the value of function at x 
 *  @return the index of new point, or -1  if point if not added 
 */
// ============================================================================
int Ostap::Math::Interpolation::Table::add  
( const double x , 
  const double y ) 
{
  const int index = m_abscissas.add ( x ) ;
  if ( index < 0 ) { return -1 ; }
  //
  m_values.insert ( m_values.begin() + index , y ) ;                 
  return index ;
}
// ============================================================================
//  get the value    for the given index 
// ============================================================================
double Ostap::Math::Interpolation::Table::y    ( const unsigned short index ) const 
{ 
  static const double s_nan = std::numeric_limits<double>::quiet_NaN () ;
  return 
    m_values.empty()        ? s_nan :
    index < m_values.size() ? m_values[ index ] : m_values.back()  ;
}
// ============================================================================
/* remove the point with the given index 
 * @param index poitn with the index 
 * @return true if point is really removed 
 */
// ============================================================================
bool  Ostap::Math::Interpolation::Table::remove ( const unsigned short index ) 
{
  if ( index < size () ) 
  {
    m_abscissas.remove ( index ) ;
    m_values.erase ( m_values.begin() + index ) ;
    return true ;  
  }
  return false ;
}
// ============================================================================
/*  Interpolation using the straightforward Lagrange interpolant 
 *  https://en.wikipedia.org/wiki/Lagrange_polynomial
 *  - it is rather slow O(n^2)
 */
// ============================================================================
double Ostap::Math::Interpolation::Table::lagrange 
( const double x ) const 
{
  return Ostap::Math::Interpolation::lagrange 
    ( m_abscissas.begin () , 
      m_abscissas.end   () , 
      m_values   .begin () , 
      m_values   .end   () , 
      x                    , 
      0.0L                 , 
      [] ( const double x ) -> double { return x ; } , 
      [] ( const double y ) -> double { return y ; } ) ;
}
// ============================================================================
/* interpolation using Neville's algorithm
 *  @see https://en.wikipedia.org/wiki/Neville%27s_algorithm
 *  - it is rather slow O(n^2)
 */
// ============================================================================
double Ostap::Math::Interpolation::Table::neville  
( const double x ) const 
{
  return Ostap::Math::Interpolation::neville 
    ( m_abscissas.begin () , 
      m_abscissas.end   () , 
      m_values   .begin () , 
      m_values   .end   () , 
      x                , 
      [] ( const double x ) -> double { return x ; } , 
      [] ( const double y ) -> double { return y ; } ) ;
}
// ============================================================================
/*  Interpolation using Neville's algorithm
 *  @see https://en.wikipedia.org/wiki/Neville%27s_algorithm
 *  - it is rather slow O(n^2)
 *  @return interpolated value and the derivative 
 */
// ============================================================================
std::pair<double,double> 
Ostap::Math::Interpolation::Table::neville2  
( const double       x ) const 
{
  return Ostap::Math::Interpolation::neville2
    ( m_abscissas.begin () , 
      m_abscissas.end   () , 
      m_values   .begin () , 
      m_values   .end   () , 
      x            , 
      [] ( const double x ) -> double { return x ; } , 
      [] ( const double y ) -> double { return y ; } ) ;
}
// ============================================================================
/*  Very simple lagrange interpolation 
 *  - it also evaluates the derivative with respect to y_i 
 *  - it is rather slow O(n^2) and numerically not very stable 
 *  @param x interpolation point
 *  @param iy index of y_i
 *  @return ( y(x) , dy/d(y_i))
 */
// ============================================================================
std::pair<double,double> 
Ostap::Math::Interpolation::Table::lagrange2 
( const double       x  , 
  const unsigned int iy ) const 
{ return Ostap::Math::Interpolation::lagrange2 ( m_abscissas.x () , 
                                                 m_values         ,
                                                 x , iy ) ; }
// ============================================================================
// interpolation using the 1st Berrut interpolant 
// ============================================================================
double Ostap::Math::Interpolation::Table::berrut1st ( const double x ) const 
{
  //
  if ( empty() ) { return 0 ; }
  //
  long double s1 = 0 ;
  long double s2 = 0 ;
  int         bi = 1 ;
  //
  const unsigned int N = size() ;
  for ( unsigned int i = 0 ; i < N ; ++ i ) 
  {
    const long double xv = m_abscissas [ i ] ;
    const long double yv = m_values    [ i ] ;
    //
    if ( s_equal ( x , xv ) ) { return yv ; } // RETURN 
    //
    const long double w = bi / ( x - xv ) ;
    //
    s1 += w * yv ;
    s2 += w      ;
    //
    bi *= -1     ;
  }
  return s1 / s2 ;
}
// ============================================================================
// interpolation using the 2nd Berrut interpolant 
// ============================================================================
double Ostap::Math::Interpolation::Table::berrut2nd ( const double x ) const 
{
  //
  if ( empty() ) { return 0 ; }
  //
  long double        s1 = 0 ;
  long double        s2 = 0 ;
  int                bi = 1 ;
  //
  const unsigned int N  = size() ;
  for ( unsigned int i = 0 ; i < N ; ++i ) 
  {
    const long double xv = m_abscissas [ i ] ;
    const long double yv = m_values    [ i ] ;
    //
    if ( s_equal ( x , xv ) ) { return yv ; } // RETURN 
    //
    const long double w  = bi / ( x - xv ) * ( ( ( 0 == i ) || ( ( i + 1 ) == N ) ) ? 1 : 2 ) ;
    //
    s1 += w * yv ;
    s2 += w      ;
    //
    bi *= -1 ;
  }
  return s1 / s2 ;
}
// ============================================================================
//  make a slice for the given abscissas 
// ============================================================================
Ostap::Math::Interpolation::Table
Ostap::Math::Interpolation::Table::slice 
( const int i , 
  const int j ) const 
{
  //
  /// allow small  negative numbers a'la Python
  //
  const int ii = i < 0 ? i + size () : i ;
  const int jj = j < 0 ? j + size () : j ;
  //
  if ( ii <  0 || jj < 0 || jj <= ii ) { return Table () ; } // RETURN 
  if ( ii == 0 && jj + 1 == size ()  ) { return *this    ; } // RETURN
  //
  return Table ( m_abscissas .begin () + ii , 
                 m_abscissas .begin () + jj , 
                 m_values    .begin () + ii , 
                 m_values    .begin () + jj , true ) ;
}
// ============================================================================
/*  simple constructor from the interpolation points  
 *  @param p input vector of abscissas  
 */
// ============================================================================
Ostap::Math::Neville::Neville 
( const Ostap::Math::Interpolation::Table& p ) 
  : Ostap::Math::Interpolation::Table ( p )
{}
// ============================================================================
/*  simple constructor from the interpolation points  
 *  @param p input vector of abscissas  
 */
// ============================================================================
Ostap::Math::Neville::Neville 
( Ostap::Math::Interpolation::Table&& p ) 
  : Ostap::Math::Interpolation::Table ( std::move ( p ) ) 
{}
// ============================================================================
/*  simple constructor from the interpolation points  
 *  @param p input vector of abscissas  
 */
// ============================================================================
Ostap::Math::Lagrange::Lagrange
( const Ostap::Math::Interpolation::Table& p ) 
  : Ostap::Math::Interpolation::Table ( p )
{}
// ============================================================================
/*  simple constructor from the interpolation points  
 *  @param p input vector of abscissas  
 */
// ============================================================================
Ostap::Math::Lagrange::Lagrange
( Ostap::Math::Interpolation::Table&& p ) 
  : Ostap::Math::Interpolation::Table ( std::move ( p ) ) 
{}
// ============================================================================
/*  simple constructor from the interpolation points  
 *  @param p input vector of abscissas  
 */
// ============================================================================
Ostap::Math::Berrut1st::Berrut1st
( const Ostap::Math::Interpolation::Table& p ) 
  : Ostap::Math::Interpolation::Table ( p )
{}
// ============================================================================
/*  simple constructor from the interpolation points  
 *  @param p input vector of abscissas  
 */
// ============================================================================
Ostap::Math::Berrut1st::Berrut1st
( Ostap::Math::Interpolation::Table&& p ) 
  : Ostap::Math::Interpolation::Table ( std::move ( p ) ) 
{}
// ============================================================================
/*  simple constructor from the interpolation points  
 *  @param p input vector of abscissas  
 */
// ============================================================================
Ostap::Math::Berrut2nd::Berrut2nd
( const Ostap::Math::Interpolation::Table& p ) 
  : Ostap::Math::Interpolation::Table ( p )
{}
// ============================================================================
/*  simple constructor from the interpolation points  
 *  @param p input vector of abscissas  
 */
// ============================================================================
Ostap::Math::Berrut2nd::Berrut2nd
( Ostap::Math::Interpolation::Table&& p ) 
  : Ostap::Math::Interpolation::Table ( std::move ( p ) ) 
{}




// ============================================================================
// Floater-Hormann rational interpolant 
// ============================================================================
/*  Simple constructor from the interpolation points  
 *  @param p input data 
 *  @param d Floater-Hormann degree parameter 
 */
// ============================================================================
Ostap::Math::FloaterHormann::FloaterHormann
( const Ostap::Math::Interpolation::Table& p ,
  const unsigned short                     d )
  : Ostap::Math::Interpolation::Table ( p )
  , m_d       ( d        ) 
  , m_weights ( p.size() ) 
{ get_weights () ; }
// ============================================================================
/*  simple constructor from the interpolation points  
 *  @param p input data 
 *  @param d Floater-Hormann degree parameter 
 */
// ============================================================================
Ostap::Math::FloaterHormann::FloaterHormann
( Ostap::Math::Interpolation::Table&&     p ,
  const unsigned short                    d )
  : Ostap::Math::Interpolation::Table ( std::move ( p ) ) 
  , m_d       ( d        ) 
  , m_weights ( p.size() ) 
{ get_weights () ; }
// ============================================================================
// calculate weigthts for Floater-Hormann interpolant 
// ============================================================================
void Ostap::Math::FloaterHormann::get_weights() 
{
  /// (1) resize container of weigths 
  m_weights.resize ( size() ) ;
  //
  const unsigned int nn = size ()  ;
  const unsigned int n  = 0 == nn ? 0 : nn - 1 ;
  //
  if ( n < m_d ) { m_d = n ; }
  //
  for ( unsigned int i = 0 ; i <= n ; ++i )
  {    
    const long double xi = x ( i ) ;
    //
    long double ib = 0 ;
    const unsigned int jmin = i       <= m_d ? 0 : i - m_d ;
    const unsigned int jmax = i + m_d <= n   ? i : n - m_d ;
    //
    for ( unsigned int j = jmin ; j <= jmax ; ++j ) 
    {
      const unsigned int kmin = j       ;
      const unsigned int kmax = j + m_d ;
      //
      long double bb = 1  ;
      for ( unsigned int k = kmin ; k <= kmax ; ++k ) 
      {
        const long double xk = x ( k ) ;
        if ( k != i ) { bb *= 1 / std::abs ( xi - xk ) ; }
      }
      ib += bb ;
    }
    /// fill table of weigths 
    m_weights [ i ] = ( i % 2 == 0 ? 1 : -1 ) * ib ;
  }
}
// ============================================================================
// the main method: get the value of Floater-Hormann interpolant 
// ============================================================================
double Ostap::Math::FloaterHormann::evaluate 
( const double x ) const 
{
  if ( empty() ) { return 0 ; }
  //
  long double s1 = 0 ;
  long double s2 = 0 ;
  //
  unsigned int N = size() ;
  for ( unsigned int i = 0  ; i < N  ;  ++i )
  {
    const long double xi = this -> x ( i ) ;
    const long double yi = this -> y ( i ) ;
    //
    if ( s_equal ( x , xi ) ) { return yi ; }  // RETURN
    //
    const double wi = m_weights[i] / ( x - xi ) ;
    //
    s1 += wi * yi ;
    s2 += wi ;
  }
  return s1 / s2 ;
}


  






// ============================================================================
// Newton interpolation polynomial
// ============================================================================
/*  simple contructor from interpolation points  
 *  @param p input data 
 */
// ============================================================================
Ostap::Math::Newton::Newton
( const Ostap::Math::Interpolation::Table& p )
  : Ostap::Math::Interpolation::Table ( p ) 
  , m_diffs ( p.size() ) 
{ get_differences () ; }
// ============================================================================
/*  simple contructor from interpolation points  
 *  @param p input data 
 */
// ============================================================================
Ostap::Math::Newton::Newton
( Ostap::Math::Interpolation::Table&& p )
  : Ostap::Math::Interpolation::Table ( std::move ( p ) ) 
  , m_diffs ( p.size() ) 
{  get_differences () ; }
// ============================================================================
// get the divided differences 
// ============================================================================
void Ostap::Math::Newton::get_differences () // get the divided differences 
{
  m_diffs.resize ( size() )  ;
  const unsigned int N   = m_diffs.size () ;
  for ( unsigned short i = 0 ; i < N ; ++i ) 
  {
    m_diffs [i] = Ostap::Math::Differences::divided 
      ( x_begin ()             , 
        x_begin () + ( i + 1 ) , 
        y_begin ()             , 
        [] ( const double x ) -> double { return x ; } , 
        [] ( const double y ) -> double { return y ; } ) ;  
  } 
}
// ============================================================================
// the main method: get the value of interpolated polynomial 
// ============================================================================
double Ostap::Math::Newton::evaluate ( const double x ) const 
{
  //
  if ( empty() ) { return 0 ; }
  //
  long double        result  = 0 ;
  long double        product = 1 ;
  const unsigned int N       = m_diffs.size () ;
  for ( unsigned short i = 0 ; i < N ; ++i )
  {
    result  += m_diffs [ i ] * product ;
    const long double xi =  this->x ( i ) ;
    product *= ( x - xi ) ;
  }
  return result ;  
}
// ============================================================================


// ============================================================================
// Barycentric Lagrange interpolation
// ============================================================================
// simple constructor from the interpolation table 
// ============================================================================
Ostap::Math::Barycentric::Barycentric
( const Ostap::Math::Interpolation::Table& p ) 
  : Interpolation::Table ( p  ) 
  , m_weights ( p.size() ) 
{ get_weights () ; }
// ============================================================================
// simple constructor from the interpolation table 
// ============================================================================
Ostap::Math::Barycentric::Barycentric
( Ostap::Math::Interpolation::Table&&      p ) 
  : Interpolation::Table ( std::move ( p ) ) 
  , m_weights ( p.size() ) 
{ get_weights () ; }
// ============================================================================
// calculate weigthts for the Barycentric interpolant 
// ============================================================================
void Ostap::Math::Barycentric::get_weights() 
{
  /// (1) resize container of weigths 
  m_weights.resize ( size() ) ;
  //
  const unsigned int N = size() ;
  //
  switch ( atype() ) 
  {
    // ========================================================================
  case Ostap::Math::Interpolation::Abscissas::Uniform : 
    // ========================================================================
    {
      for ( unsigned int i = 0 ; i < N ; ++ i ) 
      { m_weights [i] = ( i % 2 ? 1 : -1 ) * Ostap::Math::choose_double ( N - 1 , i ) ; }
      //
      return ;  // RETURN 
    }
    // ========================================================================
  case Ostap::Math::Interpolation::Abscissas::Chebyshev1 : 
    // ========================================================================
    {
      for ( unsigned int i = 0 ; i < N ; ++ i ) 
      {
        const long double a = 1.0L * ( 2 * ( N - i ) - 1 ) * M_PIl / ( 2 * N ) ;
        const long double x = std::sin ( a ) ;
        m_weights [i] = ( i % 2 ? 1 : -1 ) * x ;
      }
      //
      return ;  // RETURN 
    }
    // ========================================================================
  case Ostap::Math::Interpolation::Abscissas::Chebyshev2 : 
    // ========================================================================    
    {
      for ( unsigned int i = 0 ; i < N ; ++ i ) 
      { m_weights [i] = ( i % 2 ? 1 : -1 ) ; }
      if ( 2 <= m_weights.size () ) 
      {
        m_weights.front () *= 0.5 ;
        m_weights.back  () *= 0.5 ;
      }
      //
      return ;  // RETURN 
    }  
    // ======================================================================
  default :
    break ;
  }
  // ========================================================================
  // Generic case 
  // ========================================================================
  for   ( unsigned short i = 0 ;  i < N ; ++i ) 
  {
    //
    const long double xi  = x ( i ) ;
    //
    long double ww = 1 ;
    for  ( unsigned short j = 0 ; j < N ; ++j  )
    { 
      if ( i != j ) 
      { 
        const long double xj  = x( j ) ; 
        ww *= ( xi - xj ) ; 
      } 
    }
    //   
    m_weights [ i ] = 1 / ww ; 
  }
}
// ============================================================================
// the main method: get the value of Floater-Hormann interpolant 
// ============================================================================
double Ostap::Math::Barycentric::evaluate 
( const double x ) const 
{
  if ( empty() ) { return 0 ; }
  //
  long double s1 = 0 ;
  long double s2 = 0 ;
  //
  const unsigned int N = size () ;
  for ( unsigned int i = 0  ; i < N ; ++i )
  {
    const long double xi =  this-> x ( i ) ;
    const long double yi =  this-> y ( i ) ;
    //
    if ( s_equal ( x , xi ) ) { return yi ; }  // RETURN
    //
    const double wi = m_weights[i] / ( x - xi ) ;
    //
    s1 += wi * yi ;
    s2 += wi ;
  }
  return s1 / s2 ;
}


// ============================================================================
// Tiele rational interpolator
// ============================================================================
// contructor from the interpolation table 
// ============================================================================
Ostap::Math::Thiele::Thiele 
( const Ostap::Math::Interpolation::Table& p ) 
  : Thiele ( p.x_begin() , 
             p.x_end  () , 
             p.y_begin() )
{}
// ============================================================================
// contructor from the interpolation table 
// ============================================================================
Ostap::Math::Thiele::Thiele 
( const Ostap::Math::Interpolation::TABLE& p ) 
  : m_x    ( p.size() ) 
  , m_y    ( p.size() ) 
  , m_rho0 ( ) 
{
  std::transform 
    ( p.begin    () ,
      p.end      () , 
      m_y.begin  () ,
      []( const auto& e )-> double { return e.second ; } ) ;  
  std::transform 
    ( p.begin    () ,
      p.end      () , 
      m_x.begin  () ,
      []( const auto& e )-> double { return e.first ; } ) ; 
  //
  this->calc_rho0 () ;
}
// ============================================================================
// calculate inverted/reciprocal differences 
// ============================================================================
void Ostap::Math::Thiele::calc_rho0 () 
{
  //
  typedef Interpolation::Abscissas::Data R  ;
  typedef std::vector<R>                 RT ;
  //
  const unsigned int N = m_x.size() ;
  //
  RT rho ( N ) ;
  //
  // 1st step 
  for ( unsigned int i = 0 ; i < N ; ++i ) 
  { rho [ i ].insert ( rho [ i ].begin() ,  N - i , m_y [ i ] ) ; }
  //
  // 2nd step
  for ( unsigned int i = 0 ; i + 1 < N ; ++i ) 
  { rho [i] [1] = ( m_x [ i ] - m_x [ i + 1 ] ) /
      ( rho [ i ] [ 0 ] - rho [ i + 1 ] [ 0 ] ) ; }
  //
  // 3rd step 
  const unsigned int NR = rho.size() ;
  for ( unsigned int i = 2 ; i < NR ; ++i ) 
  { for ( unsigned int j = 0 ; i + j < NR ; ++j ) 
    { rho  [ j ] [ i ] = 
        ( m_x [ j ] - m_x [ j + i ] ) / ( rho [ j ] [ i - 1 ] - rho [ j + 1 ] [ i - 1 ] )
        + rho [ j + 1 ] [ i - 2 ] ; } }
  //
  m_rho0 = rho [0] ;
}
// ============================================================================
// the main method: get the value of interpolant 
// ============================================================================
double Ostap::Math::Thiele::evaluate    ( const  double x ) const 
{
  const unsigned int NR = m_rho0.size() ;
  //
  double r   = 0.0 ;
  for ( unsigned int i = NR - 1 ;  1 < i  ; --i ) 
  { r = ( x - m_x [ i - 1 ] ) / ( m_rho0 [ i ] - m_rho0 [ i - 2 ] + r ) ; }
  //
  return m_y [ 0 ] + ( x - m_x [ 0 ] ) / ( m_rho0 [ 1 ] + r ) ;
}
// ============================================================================
//  get the abscissas for the given index 
// ============================================================================
double Ostap::Math::Thiele::x    ( const unsigned short index ) const 
{ 
  static const double s_nan = std::numeric_limits<double>::quiet_NaN () ;
  return 
    m_x.empty() ? s_nan : 
    index < m_x.size() ? m_x [ index ] : m_x.back() ; 
}
// ============================================================================
//  get the value    for the given index 
// ============================================================================
double Ostap::Math::Thiele::y    ( const unsigned short index ) const 
{ 
  static const double s_nan = std::numeric_limits<double>::quiet_NaN () ;
  return 
    m_y.empty()        ? s_nan :
    index < m_y.size() ? m_y [ index ] : m_y.back()  ;
}
// ============================================================================


// ============================================================================
/*  very simple lagrange interpolation 
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
                    [] ( double z ) { return z ; } , 
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
                   [] ( double z ) { return z ; } , 
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
                    [] ( double z ) { return z ; } , 
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

