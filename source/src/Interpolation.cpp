// ===========================================================================
// Include files
// ===========================================================================
// Ostap
// ===========================================================================
#include "Ostap/Math.h"
#include "Ostap/Interpolation.h"
#include "Ostap/Differences.h"
#include "Ostap/Choose.h"
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
  : m_x ( x ) 
{ if  ( !sorted ) { sort() ; } }
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
( const unsigned short                              n    ,
  const double                                      low  , 
  const double                                      high  , 
  const Ostap::Math::Interpolation::Abscissas::Type t    ) 
  : m_x ( n , 0.0 )
{  
  const long double mn = std::min ( low , high ) ;
  const long double mx = std::max ( low , high ) ;
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
  { m_x.push_back ( xnew ) ;                   return m_x.size() - 1 ; } // RETURN
  // there exists element with the same value 
  else if ( !s_num_less ( xnew , *ifound ) ) { return            - 1 ; } // RETURN
  // insert new element at the given position
  const int index = ifound - m_x.begin() ;
  m_x.insert ( ifound , xnew ) ;
  //
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
  m_x.erase ( m_x.begin() + index ) ;
  return true ;
}
// ============================================================================

// ============================================================================
// Interpolatrion table 
// ============================================================================
// sort internal data and remove duplicated abscissas 
// ============================================================================
void Ostap::Math::Interpolation::Table::get_sorted ( const bool sorted ) 
{ 
  // =========================================================================
  if ( !sorted ) 
  {
    std::sort ( m_table.begin () , 
                m_table.end   () , 
                [] ( auto a , auto b ) 
                { return s_num_less ( a.first , b.first ) ; } ) ;
  }
  /// remove duplicates if any 
  m_table.erase ( std::unique ( m_table.begin () , 
                                m_table.end   () , 
                                [] ( auto a , auto b ) 
                                { return s_num_equal ( a.first , b.first ) ; } ) , 
                  m_table.end ()  ) ;
}
// ====================================================================
/* the simplest constructor 
 * @param data   input data 
 * @param sorted indicate if data already sorted 
 */
// ============================================================================
Ostap::Math::Interpolation::Table::Table
( const Ostap::Math::Interpolation::TABLE& data   , 
  const bool                               sorted ) 
  :  m_table ( data )
{ if  ( !sorted ) { get_sorted() ; } }
// ============================================================================
/*  simple constructor from abscissas and y-list 
 *  @param x input vector of abscissas  
 *  @param y input vector of y 
 *  - if vector of y is longer  than vector x, extra values are ignored 
 *  - if vector of y is shorter than vector x, missing entries are assumed to be zero  
 */
// ============================================================================
Ostap::Math::Interpolation::Table::Table
( const Ostap::Math::Interpolation::Abscissas&       x , 
  const Ostap::Math::Interpolation::Abscissas::Data& y ) 
  : Table ( x.begin() , x.end() , y.begin() , y.end() , true ) 
{}
// ============================================================================
/*  simple contructor from x&y-lists 
 *  @param x input vector of x 
 *  @param y input vector of y 
 *  - if vector of y is longer  than vector x, extra values are ignored 
 *  - if vector of y is shorter than vector x, missing entries are assumed to be zero  
 *  @attention duplicated abscissas will be removed 
 */
// ============================================================================
Ostap::Math::Interpolation::Table::Table
( const Ostap::Math::Interpolation::Abscissas::Data& x , 
  const Ostap::Math::Interpolation::Abscissas::Data& y ) 
  : Table ( x.begin() , x.end() , y.begin() , y.end() ) 
{}
// ============================================================================
// get the abscissas 
// ============================================================================
Ostap::Math::Interpolation::Abscissas 
Ostap::Math::Interpolation::Table::abscissas () const 
{
  return Abscissas 
    ( m_table.begin () , 
      m_table.end   () , 
      [] ( const auto& p ) { return p.first ; } , 
      true ) ; 
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
  /// comparison criteria: compare only abscissas 
  auto cmp_x = [] ( const auto& p1 , const auto&p2 ) 
    { return s_num_less ( p1.first , p2.first ) ; } ;
  
  // where to insert new element?
  const TABLE::value_type pnew { std::make_pair ( x , y ) } ;
  TABLE::iterator ifound = std::lower_bound 
    ( m_table.begin () , 
      m_table.end   () , 
      pnew             , cmp_x ) ;
  
  // 1) all elements are less than new element
  if  ( m_table.end() == ifound ) 
  { m_table.push_back ( pnew ) ;          return m_table.size() - 1 ; } // RETURN  
  // 2) there exists element with the same value 
  else if ( !cmp_x ( pnew , *ifound ) ) { return                - 1 ; } // RETURN
  // 3) insert new element at the given position
  const int index = ifound - m_table.begin() ;
  m_table.insert ( ifound , pnew ) ;
  //
  return index ;
}
// ============================================================================
/* remove the point with the given index 
 * @param index poitn with the index 
 * @return true if point is really removed 
 */
// ============================================================================
bool  Ostap::Math::Interpolation::Table::remove ( const unsigned short index ) 
{
  if ( index >= m_table.size() ) { return false ; }
  m_table.erase ( m_table.begin() + index ) ;
  return true ;
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
    ( m_table.begin () , 
      m_table.end   () , 
      m_table.begin () , 
      m_table.end   () , 
      x                , 
      0.0L             , 
      [] ( const auto& p ) { return p.first  ; } , 
      [] ( const auto& p ) { return p.second ; } ) ;
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
    ( m_table.begin () , 
      m_table.end   () , 
      m_table.begin () , 
      m_table.end   () , 
      x                , 
      [] ( const auto& p ) { return p.first  ; } , 
      [] ( const auto& p ) { return p.second ; } ) ;
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
    ( m_table.begin () , 
      m_table.end   () , 
      m_table.begin () , 
      m_table.end   () , 
      x            , 
      [] ( const auto& p ) { return p.first  ; } , 
      [] ( const auto& p ) { return p.second ; } ) ;
}
// ============================================================================
/*  Very simple lagrange interpolation 
 *  - it also evaluates the derivative wity respect to y_i 
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
{ return Ostap::Math::Interpolation::lagrange2 ( m_table , x , iy ) ; }
// ============================================================================
//  make a slice for the given abscissas 
// ============================================================================
Ostap::Math::Interpolation::Table
Ostap::Math::Interpolation::Table::slice ( const int i , const int j ) const 
{
  //
  // treat the negative indices 
  const int ii = 
    i < 0     ? 0        : i ; // ignore  negative  i 
  const int jj = 
    j <  0    ? j + n () : 
    j >= n () ?     n () : j ; // adjust a bit not very negative j  
  //
  if      ( 0 >= ii && n() <= jj ) { return *this ; }
  //  regular slice 
  else if ( 0 <= ii && i   <  jj )
  {
    Table newt {} ;
    newt.m_table.resize ( jj - ii ) ;
    std::copy ( m_table      .begin () + ii , 
                m_table      .begin () + jj ,
                newt.m_table .begin ()      ) ;
    return newt  ; 
  }
  // empty slice 
  return Table() ;
}
// ============================================================================

// ============================================================================
// Barycentric Weights  
// ============================================================================
// constructor from abscissas 
// ============================================================================
Ostap::Math::Interpolation::Weights::Weights 
( const Ostap::Math::Interpolation::Abscissas& a       )
  :  m_a ( a ) 
  ,  m_w ( a.size() , 0.0 ) 
{ get_weights() ; }
// ============================================================================
// constructor from abscissas 
// ============================================================================
Ostap::Math::Interpolation::Weights::Weights 
( const Ostap::Math::Interpolation::Abscissas::Data& a )
  :  m_a ( a ) 
  ,  m_w () 
{
  m_w .resize ( m_a.size() ) ;
  get_weights() ; 
}
// ============================================================================
// constructor from abscissas 
// ============================================================================
Ostap::Math::Interpolation::Weights::Weights 
( const unsigned short                               n    ,
  const double                                       xmin , 
  const double                                       xmax , 
  const  Ostap::Math::Interpolation::Abscissas::Type t    ) 
  : m_a ( n , xmin , xmax , t ) 
  , m_w () 
{
  m_w .resize ( m_a.size() ) ;
  //
  switch ( t ) 
  {
    // ========================================================================
  case Ostap::Math::Interpolation::Abscissas::Uniform : 
    // ========================================================================
    {
      for ( unsigned short i = 0 ; i < n ; ++ i ) 
      { m_w [i] = ( i % 2 ? 1 : -1 ) * Ostap::Math::choose_double ( n - 1 , i ) ; }
      //
      const auto imax = std::max_element ( m_w.begin() , 
                                           m_w.end  () , 
                                           [] ( const double a , const double b ) 
                                           { return std::abs ( a ) < std::abs ( b ) ; } ) ;
      //
      if ( m_w.end() != imax ) 
      {
        std::pair<double,int> a = Ostap::Math::frexp2 ( *imax ) ;
        if ( 2 < std::abs ( a.second ) ) { Ostap::Math::ldexp ( m_w , -a.second ) ; }
      }
      break ;
    }
    // ========================================================================
  case Ostap::Math::Interpolation::Abscissas::Chebyshev1 : 
    // ========================================================================
    {
      for ( unsigned short i = 0 ; i < n ; ++ i ) 
      {
        const long double a = 1.0L * ( 2 * ( n - i ) - 1 ) * M_PIl / ( 2 * n ) ;
        const long double x = std::sin ( a ) ;
        m_w [i] = ( i % 2 ? 1 : -1 ) * x ;
      }
      break ; 
    }
    // ========================================================================
  case Ostap::Math::Interpolation::Abscissas::Chebyshev2 : 
    // ========================================================================    
    {
      for ( unsigned short i = 0 ; i < n ; ++ i ) { m_w [i] = ( i % 2 ? 1 : -1 ) ; }
      m_w.front () *= 0.5 ;
      m_w.back  () *= 0.5 ;
      //
      break ;
    }
    // ========================================================================
  default :  
    // ========================================================================
    get_weights () ;
    // ========================================================================
  }
  // ==========================================================================
}
// ============================================================================
// calculate the weights 
// ============================================================================
void Ostap::Math::Interpolation::Weights::get_weights() 
{
  const unsigned short n = m_a.size() ;
  //
  for   ( unsigned short i = 0 ;  i < n ; ++i ) 
  {
    const long double xi  = m_a.x ( i ) ;
    //
    long double ww = 1 ;
    for  ( unsigned short j = 0 ; j < n ; ++j  )
    { 
      if ( i != j ) {  ww *= ( xi - m_a.x ( j ) ) ; }
    }
    //   
    m_w[ i ] = 1 / ww ;
  }
  //
  const auto imax = 
    std::max_element ( m_w.begin() , 
                       m_w.end  () , 
                       [] ( const double a , const double b ) 
                       { return std::abs ( a ) < std::abs ( b ) ; } ) ;
  //
  std::pair<double,int> a = Ostap::Math::frexp2 ( *imax ) ;
  if ( 2 < std::abs ( a.second ) ) { Ostap::Math::ldexp ( m_w , -a.second ) ; }
}
// ============================================================================
/*  add the point x into  collection 
 *  @param x abscissas of the point to be added 
 *  @return the index of new point, or -1  if point if not added 
 */
// ============================================================================
int Ostap::Math::Interpolation::Weights::add ( const  double x ) 
{
  //
  const int index = m_a.add ( x ) ;
  if ( index < 0 ) { return index ; }
  //
  m_w.insert ( m_w.begin() + index , 0.0 ) ; //  temporatily to keep indices valid 
  const int    N   = m_a.size() ;
  long  double sum = 0 ;
  for ( int i = 0 ; i < N ; ++i ) 
  {
    if ( i == index ) { continue ; }
    const long double xi = m_a.x( i ) ;
    const long double wi = m_w[i] / (  xi  - x ) ;
    m_w[i]  = wi ;
    sum    += wi ;
  }
  /// finally get the new weight 
  m_w [ index ] = -1 * sum ;
  //
  // rescale weights if needed 
  const auto imax = 
    std::max_element ( m_w.begin() , 
                       m_w.end  () , 
                       [] ( const double a , const double b ) 
                       { return std::abs ( a ) < std::abs ( b ) ; } ) ;
  //
  std::pair<double,int> a = Ostap::Math::frexp2 ( *imax ) ;
  if ( 2 < std::abs ( a.second ) ) { Ostap::Math::ldexp ( m_w , -a.second ) ; }
  //
  return index ;
}
// =============================================================================
/*  remove the point with the  given index 
 *  @param index the point to be removed 
 *  @return   true if point is removed 
 */
// =============================================================================
bool Ostap::Math::Interpolation::Weights::remove ( unsigned short index ) 
{
  if ( index >= m_a.size() ) { return false ; }
  const long double xj = m_a.x ( index ) ;
  const unsigned short N = n() ;
  for  ( unsigned short i = 0 ; i < N  ; ++i ) 
  {
    if ( index == i ) { continue ; }
    const long double xi = m_a.x ( i ) ;
    m_w [ i ] *= (  xi - xj ) ;
  }
  m_a.remove  ( index ) ;
  m_w.erase   ( m_w.begin() + index ) ;
  //
  return true ;
}
// ============================================================================
//  make a slice for the given abscissas 
// ============================================================================
Ostap::Math::Interpolation::Weights
Ostap::Math::Interpolation::Weights::slice 
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
    0 <= ii && ii <  jj  ? 
    Weights ( m_a.slice ( ii , jj ) ) : Weights ( ) ;                             
}
// ============================================================================
/*  simple contructor from interpolation points  
 *  @param x input vector of abscissas  
 */
// ============================================================================
Ostap::Math::Neville::Neville
  ( const Ostap::Math::Interpolation::Table& p ) 
  : Ostap::Math::Interpolation::Table ( p )
{}
// ============================================================================
/*  simple contructor from interpolation points  
 *  @param x input vector of abscissas  
 */
// ============================================================================
Ostap::Math::Lagrange::Lagrange
( const Ostap::Math::Interpolation::Table& p ) 
  : Ostap::Math::Interpolation::Table ( p )
{}
// ============================================================================

// ============================================================================
// Newton innterpolation polynomial
// ============================================================================
/** simple contructor from interpolation points  
 *  @param x input vector of abscissas  
 */
// ============================================================================
Ostap::Math::Newton::Newton ( const Ostap::Math::Interpolation::Table& p )
  : m_table ( p        ) 
  , m_diffs ( p.size() ) 
{
  get_differences () ;
} 
// ============================================================================
// get the divided differences 
// ============================================================================
void Ostap::Math::Newton::get_differences () // get the divided differences 
{
  m_diffs.resize ( m_table.size() )  ;
  const unsigned int N   = m_diffs.size () ;
  for ( unsigned short i = 0 ; i < N ; ++i ) 
  {
    m_diffs [i] = Ostap::Math::Differences::divided 
      ( m_table.begin()             , 
        m_table.begin() + ( i + 1 ) , 
        m_table.begin()             , 
        [] ( const auto& p )->double { return p.first  ; } , 
        [] ( const auto& p )->double { return p.second ; } ) ;  
  } 
}
// ============================================================================
// the main method: get the value of interpolated polynomial 
// ============================================================================
double Ostap::Math::Newton::evaluate ( const double x ) const 
{
  long double        result  = 0 ;
  long double        product = 1 ;
  const unsigned int N       = m_diffs.size () ;
  for ( unsigned short i = 0 ; i < N ; ++i )
  {
    result  += m_diffs[i] * product ;
    product *= ( x * 1.0L - m_table.x ( i ) ) ;
  }
  return result ;  
}
// ============================================================================
/** add the point (x,y) into interpolation table 
 *  @param x abscissas of the point to be added 
 *  @param y the value of function at x 
 *  @return the index of new point, or -1  if point if not added 
 */
// ============================================================================
int Ostap::Math::Newton::add ( const  double x ,   const  double y ) 
{
  const int index = m_table.add ( x , y ) ;
  if  ( 0 <= index ) { get_differences () ; } //   recalculate differences 
  return index ;
}
// ============================================================================
/* remove the point with the  given index 
 *  @param index the point to be removed 
 *  @return   true if point is removed 
 */
// ============================================================================
bool Ostap::Math::Newton::remove ( unsigned short index ) 
{
  const bool removed = m_table.remove ( index ) ;
  if ( removed ) { get_differences() ; } // recalculate the differences 
  return removed ;
}
// ============================================================================
// make a slice for the given range of points
// ============================================================================
Ostap::Math::Newton
Ostap::Math::Newton::slice ( const int i , const int j ) const 
{
  //
  Newton tmp{} ;
  tmp.m_table = m_table.slice ( i , j ) ;
  tmp.get_differences() ;
  //
  return tmp ;
}
// ============================================================================
 


// ============================================================================
// Barycentric Lagrange interpolation
// ============================================================================
/* simple constructor from the interpolation table 
 */
// ============================================================================
Ostap::Math::Barycentric::Barycentric
  ( const Ostap::Math::Interpolation::Table& data ) 
  : m_w ( data.begin () , data.end () , [] ( const auto& p ) { return p.first ; } , true )
  , m_y ( data.size() )
{
  std::transform ( data.begin () , 
                   data.end   () ,
                   m_y.begin  () , 
                   [] ( const auto&p ) { return  p.second ; } ) ; 
}
// ============================================================================
/** simple constructor from the interpolation data 
 */
// ============================================================================
Ostap::Math::Barycentric::Barycentric
( const Ostap::Math::Interpolation::TABLE& data   , 
  const bool                               sorted ) 
  : Barycentric ( Ostap::Math::Interpolation::Table ( data , sorted ) ) 
{}
// ============================================================================
/** simple constructor from x&y-lists 
 *  @param x input vector of x 
 *  @param y input vector of y 
 *  - if vector of y is longer  than vector x, extra values are ignored 
 *  - if vector of y is shorter than vector x, missing entries are assumed to be zero  
 */
// ============================================================================O
Ostap::Math::Barycentric::Barycentric
( const Ostap::Math::Interpolation::Abscissas& x , 
  const Ostap::Math::Barycentric::Data&        y ) 
  : m_w ( x ) 
  , m_y ( y )
{
  if      ( m_w.size() > m_y.size() ) 
  { m_y.insert ( m_y.end() , m_w.size() - m_y.size() , 0.0 ) ; }
  else if ( m_w.size() < m_y.size() ) 
  { m_y.erase  ( m_y.begin () + m_w.size() , m_y.end ()   ) ; }
}
// ============================================================================
/** simple constructor from x&y-lists 
 *  @param x input vector of x 
 *  @param y input vector of y 
 *  - if vector of y is longer  than vector x, extra values are ignored 
 *  - if vector of y is shorter than vector x, missing entries are assumed to be zero  
 */
// ============================================================================O
Ostap::Math::Barycentric::Barycentric
( const Ostap::Math::Interpolation::Weights&   x , 
  const Ostap::Math::Barycentric::Data&        y ) 
  : m_w ( x ) 
  , m_y ( y )
{
  if      ( m_w.size() > m_y.size() ) 
  { m_y.insert ( m_y.end() , m_w.size() - m_y.size() , 0.0 ) ; }
  else if ( m_w.size() < m_y.size() ) 
  { m_y.erase  ( m_y.begin () + m_w.size() , m_y.end ()   ) ; }
}
// ============================================================================
/*  Simple constructor from x&y-lists 
 *  @param x input vector of x 
 *  @param y input vector of y 
 *  - if vector of y is longer  than vector x, extra values are ignored 
 *  - if vector of y is shorter than vector x, missing entries are assumed to be zero  
 *  @attention duplicated abscissas will be removed 
 */
// ============================================================================
Ostap::Math::Barycentric::Barycentric
( const Ostap::Math::Barycentric::Data& x , 
  const Ostap::Math::Barycentric::Data& y ) 
  : Barycentric ( Ostap::Math::Interpolation::Table ( x , y ) ) 
{}
// ============================================================================
/** add the point (x,y)  into  collection 
 *  @param x abscissas of the point to be added 
 *  @param y the value of function at x 
 *  @return the index of new point, or -1  if point if not added 
 */
// ============================================================================
int Ostap::Math::Barycentric::add ( const  double x ,   const  double y ) 
{
  const int index = m_w.add ( x ) ;
  if ( index < 0 ) { return index ; }
  m_y.insert ( m_y.begin() + index , y ) ;
  return index ;
}
// ============================================================================
/*  remove the point with the  given index 
 *  @param index the point to be removed 
 *  @return   true if point is removed 
 */
// ============================================================================
bool Ostap::Math::Barycentric::remove ( const unsigned short index ) 
{
  if ( !m_w.remove ( index ) ) { return false ; }
  m_y.erase ( m_y.begin() +  index ) ;
  return false ;
}
// ============================================================================
//  make a slice for the given abscissas 
// ============================================================================
Ostap::Math::Barycentric
Ostap::Math::Barycentric::slice 
( const unsigned short i , 
  const unsigned short j ) const 
{
  return 
    i == 0  && j == n() ? (*this) :
    i <  j  && j <= n() ? 
    Barycentric ( m_w.slice ( i , j ) , Data ( m_y.begin() + i , m_y.begin() + j ) ) :
    Barycentric () ;
}
// ============================================================================
// evaluate the interpolation polynomial 
// ============================================================================
double Ostap::Math::Barycentric::evaluate   ( const double x ) const 
{
  const unsigned short N = n() ;
  long double s1 = 0 ;
  long double s2 = 0 ;
  for ( unsigned short i = 0 ; i < N ; ++i ) 
  {
    //
    const long double xi = m_w . x ( i ) ;
    if ( s_num_equal ( x , xi ) ) { return y ( i ) ; }   // RETURN!!! 
    //
    const long double c  = m_w . w ( i ) / ( x - xi ) ;
    const long double yi = m_y [ i ] ;
    //
    s1  = std::fma ( c , yi , s1 )  ;
    s2 += c ;
    //
  }
  return s1 / s2 ;
}
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

