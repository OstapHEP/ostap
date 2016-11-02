#ifndef OSTAP_INTERPOLATION_H 
#define OSTAP_INTERPOLATION_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <vector>
#include <array>
#include <utility>
#include <algorithm>
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math  
  {
    // ========================================================================
    /** @namespace Gaudi::Math::Interpolation
     *  Collection of simple utilities for various types of interpolation
     *  - lagrange interpolation 
     *  - Neville  interpolation 
     *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
     *  @see https://en.wikipedia.org/wiki/Neville%27s_algorithm
     *
     *  None of this method should be applies for "long" sequence of 
     *  interpolation points ( e.g. >20), especially for uniform grid 
     *  (see https://en.wikipedia.org/wiki/Runge%27s_phenomenon)
     *  
     *  Lagrange interpoaltion is numerically not stable, and 
     *  Neville's algorithm has (a little bit) better numerical stability.
     *
     *  Using Lagrange algorithm it is easy to get derivative with 
     *  respect to the data points, while using Neville's algorithm 
     *  one can easily calcualate the derivative with respect to 
     *  the argument.  
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2016-07-23
     */
    namespace Interpolation
    {
      // ======================================================================
      /// the actual type of "simple" data 
      typedef std::vector<std::pair<double,double> >           DATA    ;
      /// the actual type of "simple" data 
      typedef std::vector<double>                              DATAVCT ;
      // ======================================================================
      /** Very simple lagrange interpolation 
       *
       *  @param  xbegin INPUT start iterator for the sequence of abscissas 
       *  @param  xend   INPUT end   iterator for the sequence of abscissas 
       *  @param  ybegin INPUT start iterator for the sequence of values 
       *  @param  yend   INPUT end   iterator for the sequence of values 
       *  @param  x      INPUT evaluate the polynomial in this point
       *  @param  xvalue INPUT adapter for x-values  
       *  @param  yvalue INPUT adapter for y-values  
       *  @return the value of Largange interpolation polynomial at point x
       *
       *  - If y-vector is shorter than x-vector, it is assumed to be zero padded
       *  - If y-vector is longer  than x-vector, extra avalues are ignored       
       *  - If x-vector is empty, polynomial is zero
       *
       *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
       *  @warning it could be CPU inefficient
       *  @warning it should *NOT* be applied for very long sequence of points 
       *           (e.g. >20) due to bad numerical  stability and Runge phenomenon
       *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
       *  @date 2016-07-23
       */
      template <class XITERATOR, 
                class YITERATOR,
                class RESULT   , 
                class XADAPTER ,
                class YADAPTER >
      inline RESULT lagrange
      ( XITERATOR    xbegin , 
        XITERATOR    xend   , 
        YITERATOR    ybegin , 
        YITERATOR    yend   , 
        const double x      , 
        RESULT       result , 
        XADAPTER     xvalue ,   // adaptor to get y-value
        YADAPTER     yvalue ) ; // adaptor to get x-value
      // ======================================================================
      /** simple interpolation using Neville's algorithm 
       *
       *  In general it should be faster than largange algorithm, 
       *  but it includes a copy on input data, that could affect CPU performance
       *  Numerically it is more stable that Lagrange interpolation, 
       *  but anyhow it also should not be used for very high polynomial 
       *  degrees (e.g.>20), especially for uniform grid (Runge phenomenon)
       *
       *  - If y-vector is shorter than x-vector, it is assumed to be zero padded
       *  - If y-vector is longer  than x-vector, extra avalues are ignored       
       *  - If x-vector is empty, polynomial is zero
       *
       *  @param  xbegin INPUT start iterator for the sequence of abscissas 
       *  @param  xend   INPUT end   iterator for the sequence of abscissas 
       *  @param  ybegin INPUT start iterator for the sequence of values 
       *  @param  yend   INPUT end   iterator for the sequence of values 
       *  @param  x      INPUT evaluate the polynomial in this point
       *  @param  xvalue INPUT adapter for x-values  
       *  @param  yvalue INPUT adapter for y-values  
       *  @return the value of interpolation polynomial at point x
       *  @see https://en.wikipedia.org/wiki/Neville%27s_algorithm
       */
      template <class XITERATOR, 
                class YITERATOR, 
                class XADAPTER ,
                class YADAPTER >
      inline double neville
      ( XITERATOR    xbegin ,
        XITERATOR    xend   ,
        YITERATOR    ybegin ,
        YITERATOR    yend   ,
        const double x      , 
        XADAPTER     xvalue ,   // adaptor to get y-value
        YADAPTER     yvalue ) ; // adaptor to get x-value  
      // ======================================================================
      /** simple interpolation using Neville's algorithm: 
       *    evaluate the interpolation polynomial and also the derivative 
       *
       *  In general it should be faster than largange algorithm, 
       *  but it includes a copy on input data, that could affect CPU performance
       *  Numerically it is more stable that Lagrange interpolation, 
       *  but anyhow it also should not be used for very high polynomial 
       *  degrees (e.g.>20), especially for uniform grid (Runge phenomenon)
       *
       *  - If y-vector is shorter than x-vector, it is assumed to be zero padded
       *  - If y-vector is longer  than x-vector, extra avalues are ignored       
       *  - If x-vector is empty, polynomial is zero
       *
       *  @param  xbegin INPUT start iterator for the sequence of abscissas 
       *  @param  xend   INPUT end   iterator for the sequence of abscissas 
       *  @param  ybegin INPUT start iterator for the sequence of values 
       *  @param  yend   INPUT end   iterator for the sequence of values 
       *  @param  x      INPUT evaluate the polynomial in this point
       *  @param  xvalue INPUT adapter for x-values  
       *  @param  yvalue INPUT adapter for y-values  
       *  @return the value of interpolation polynomial at point x
       *  @see https://en.wikipedia.org/wiki/Neville%27s_algorithm
       */
      template <class XITERATOR, 
                class YITERATOR, 
                class XADAPTER ,
                class YADAPTER >
      inline 
      std::pair<double,double> neville2
      ( XITERATOR    xbegin ,
        XITERATOR    xend   ,
        YITERATOR    ybegin ,
        YITERATOR    yend   ,
        const double x      , 
        XADAPTER     xvalue ,   // adaptor to get y-value
        YADAPTER     yvalue ) ; // adaptor to get x-value  
      // ======================================================================
      /** simple interpolation using Neville's algorithm 
       *
       *  In general it should be faster than largange algorithm, 
       *  but it modified input data!
       *
       *  Numerically it is more stable than Lagrange interpolation, 
       *  but anyhow it also should not be used for very high polynomial 
       *  degrees (e.g.>20), especially for uniform grid (Runge phenomenon)
       *
       *  @attention y-sequence must be "simple", convertible to doubles  
       *  @attention y-sequence is *MODIFIED*
       *
       *  @param  xbegin INPUT  start iterator for the sequence of abscissas 
       *  @param  xend   INPUT  end   iterator for the sequence of abscissas 
       *  @param  ybegin UPDATE start iterator for the sequence of values 
       *  @param  x      INPUT  evaluate the polynomial in this point
       *  @param  xvalue INPUT  adapter for x-values  
       *  @param  yvalue INPUT  adapter for y-values  
       *  @return the value of  interpolation polynomial at point x
       *  @see https://en.wikipedia.org/wiki/Neville%27s_algorithm
       */
      template <class XITERATOR, 
                class YITERATOR, 
                class XADAPTOR >
      inline double 
      neville
      ( XITERATOR    xbegin , 
        XITERATOR    xend   , 
        YITERATOR    ybegin , // NON-const!
        const double x      , 
        XADAPTOR     xvalue ) ;
      // ======================================================================      
      /** simple interpolation using Neville's algorithm with simultaneous 
       *  estimation of the derivative  
       *
       *  Numerically it is more stable than Lagrange interpolation, 
       *  but anyhow it also should not be used for very high polynomial 
       *  degrees (e.g.>20), especially for uniform grid (Runge phenomenon)
       *
       *  @attention y-sequence must be "simple", convertible to doubles  
       *  @attention y-sequence is *MODIFIED*
       *  @attention d-sequence is *MODIFIED*
       *
       *  @param  xbegin INPUT  start iterator for the sequence of abscissas 
       *  @param  xend   INPUT  end   iterator for the sequence of abscissas 
       *  @param  ybegin UPDATE start iterator for the sequence of values 
       *  @param  x      INPUT  evaluate the polynomial in this point
       *  @param  xvalue INPUT  adapter for x-values  
       *  @param  yvalue INPUT  adapter for y-values  
       *  @return the pair (function,derivative) at point x 
       *  @see https://en.wikipedia.org/wiki/Neville%27s_algorithm
       */
      template <class XITERATOR, 
                class YITERATOR, 
                class DITERATOR, 
                class XADAPTOR >
      inline std::pair<double,double>
      neville
      ( XITERATOR    xbegin , 
        XITERATOR    xend   , 
        YITERATOR    ybegin , // NON-const!
        DITERATOR    dbegin , // NON-const!
        const double x      , 
        XADAPTOR     xvalue ) ;
      // ======================================================================      
      /** very simple lagrange interpolation 
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
      double lagrange ( const std::vector<double>& xs , 
                        const std::vector<double>& ys , 
                        const double               x  ) ;
      // ======================================================================      
      /** very simple lagrange interpolation 
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
      double lagrange ( const DATA& data , const double x  ) ;
      // ======================================================================      
      /** Simple lagrange interpolation 
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
      std::pair<double,double>
      lagrange2 ( const std::vector<double>& xs , 
                  const std::vector<double>& ys , 
                  const double               x  , 
                  const unsigned int         iy ) ;
      // ======================================================================      
      /** Simple lagrange interpolation 
       *  - it also evaluate the derivative with respect to y_i 
       *  @param  data INPUT sequence of (x,y)
       *  @param  x    INPUT evaluate the polynomial in this point
       *  @param  it   INPUT index of y_i, the derivative shodul be calculated.
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
      std::pair<double,double>
      lagrange2 ( const DATA&        data , 
                  const double       x    , 
                  const unsigned int iy   ) ;
      // ======================================================================      
      /** very simple Neville's interpolation 
       *  @param  xs INPUT sequence of abscissas 
       *  @param  ys INPUT sequence of values 
       *  @param  x  INPUT evaluate the polynomial in this point
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
      double neville  ( const std::vector<double>& xs ,
                        const std::vector<double>& ys , 
                        const double               x  ) ;
      // ======================================================================      
      /** very simple Neville interpolation 
       *  @param  data INPUT sequence of (x,y)
       *  @param  x    INPUT evaluate the polynomial in this point
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
      double 
      neville  ( const DATA&  data , 
                 const double x    ) ;
      // ======================================================================      
      /** very simple Neville's interpolation 
       *  -  it evalutes the polynomial and the derivative
       *  @param  xs INPUT sequence of abscissas 
       *  @param  ys INPUT sequence of values 
       *  @param  x  INPUT evaluate the polynomial in this point
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
      std::pair<double,double> 
      neville2  ( const std::vector<double>& xs ,
                  const std::vector<double>& ys , 
                  const double               x  ) ;
      // ======================================================================      
      /** very simple Neville interpolation 
       *  -  it evaluates the polynomial and the derivative
       *  @param  data INPUT sequence of (x,y)
       *  @param  x    INPUT evaluate the polynomial in this point
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
      std::pair<double,double> 
      neville2  ( const DATA&  data , 
                  const double x  ) ;
      // ======================================================================      
    } //                            end of namespace Gaudi::Math::Interpolation 
    // ========================================================================
  } //                                             end of namespace Gaudi::Math
  // ==========================================================================
} //                                                     end of namespace Gaudi
// ============================================================================
/*  Very simple lagrange interpolation 
 *
 *  @param  xbegin INPUT start iterator for the sequence of abscissas 
 *  @param  xend   INPUT end   iterator for the sequence of abscissas 
 *  @param  ybegin INPUT start iterator for the sequence of values 
 *  @param  yend   INPUT end   iterator for the sequence of values 
 *  @param  x      INPUT evaluate the polynomial in this point
 *  @param  xvalue INPUT adapter for x-values  
 *  @param  yvalue INPUT adapter for y-values  
 *  @return the value of Largange interpolation polynomial at point x
 *
 *  - If y-vector is shorter than x-vector, it is assumed to be zero padded
 *  - If y-vector is longer  than x-vector, extra avalues are ignored       
 *  - If x-vector is empty, polynomial is zero
 *
 *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
 *  @warning it could be CPU inefficient
 *  @warning it should *NOT* be applied for very long sequence of points 
 *           (e.g. >20) due to bad numerical  stability and Runge phenomenon
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2016-07-23
 */
// ============================================================================
template <class XITERATOR, 
          class YITERATOR, 
          class RESULT   ,
          class XADAPTER ,
          class YADAPTER >
inline RESULT 
Ostap::Math::Interpolation::lagrange
( XITERATOR    xbegin , 
  XITERATOR    xend   , 
  YITERATOR    ybegin , 
  YITERATOR    yend   , 
  const double x      , 
  RESULT       result ,
  XADAPTER     xvalue ,   // adaptor to get y-value
  YADAPTER     yvalue )   // adaptor to get x-value
{
  /// 1)special case no data 
  if ( xbegin == xend || ybegin == yend ) { return result ; }    // RETURN
  const typename std::iterator_traits<XITERATOR>::difference_type nx = std::distance ( xbegin , xend ) ;
  /// 2) special case: constant function 
  if ( 1 == nx ) { return  result + yvalue ( *ybegin )    ; } // RETURN
  const typename std::iterator_traits<XITERATOR>::difference_type ny = std::distance ( ybegin , yend ) ;
  //
  /// redefine yend : skip extra y-values 
  if ( nx < ny ) { std::advance ( yend , nx - ny ) ; }
  //
  unsigned int i  =      0 ;
  XITERATOR    ix = xbegin ;
  YITERATOR    iy = ybegin ;
  //
  for ( ; iy != yend ; ++iy , ++ix , ++i ) 
  {
    const double xi = xvalue ( *ix ) ; // get values 
    //
    double       r = 1 ;
    unsigned int j = 0 ;
    for ( XITERATOR jx = xbegin ; jx != xend ; ++jx, ++j  )
    {
      const double xj = xvalue ( *jx ) ;  // get the value 
      if ( i == j ) { continue ; }
      r *= ( x - xj ) / ( xi - xj ) ;
    }
    //
    result += r * yvalue ( *iy ) ;
  }
  return result ;                            // RETURN 
}
// ============================================================================
/** simple interpolation using Neville's algorithm 
 *
 *  In general it should be faster than largange algorithm, 
 *  but it includes a copy on input data, that could affect CPU performance
 *  Numerically it is more stable that Lagrange interpolation, 
 *  but anyhow it also should not be used for very high polynomial 
 *  degrees (e.g.>20), especially for uniform grid (Runge phenomenon)
 *
 *  - If y-vector is shorter than x-vector, it is assumed to be zero padded
 *  - If y-vector is longer  than x-vector, extra avalues are ignored       
 *  - If x-vector is empty, polynomial is zero
 *
 *  @param  xbegin INPUT start iterator for the sequence of abscissas 
 *  @param  xend   INPUT end   iterator for the sequence of abscissas 
 *  @param  ybegin INPUT start iterator for the sequence of values 
 *  @param  yend   INPUT end   iterator for the sequence of values 
 *  @param  x      INPUT evaluate the polynomial in this point
 *  @param  xvalue INPUT adapter for x-values  
 *  @param  yvalue INPUT adapter for y-values  
 *  @return the value of interpolation polynomial at point x
 *  @see https://en.wikipedia.org/wiki/Neville%27s_algorithm
 */
// ============================================================================
template <class XITERATOR, 
          class YITERATOR, 
          class XADAPTER ,
          class YADAPTER >
inline double 
Ostap::Math::Interpolation::neville
( XITERATOR    xbegin ,
  XITERATOR    xend   ,
  YITERATOR    ybegin ,
  YITERATOR    yend   ,
  const double x      , 
  XADAPTER     xvalue , // adaptor to get y-value
  YADAPTER     yvalue ) // adaptor to get x-value  
{
  /// 1)special case no data 
  if ( xbegin == xend || ybegin == yend ) { return                   0 ; } // RETURN
  const typename std::iterator_traits<XITERATOR>::difference_type NX = std::distance ( xbegin , xend ) ;  
  /// 2) special case: constant function 
  if ( 1 == NX                          ) { return  yvalue ( *ybegin ) ; } // RETURN
  const typename std::iterator_traits<YITERATOR>::difference_type NY = std::distance ( ybegin , yend ) ;  
  // temporary storage 
  std::vector<double> _y ( NX , 0 ) ;
  YITERATOR ylast = yend ;
  if ( NX < NY ) { std::advance ( ylast , NX - NY ) ; }
  std::transform ( ybegin , ylast , _y.begin() , yvalue ) ;
  // simple version of neville 
  return neville ( xbegin , xend , _y.begin() , x , xvalue ) ;
}
// ============================================================================
/** simple interpolation using Neville's algorithm:
 *  - evaluate the interpolation polynomial and the derivative 
 *
 *  In general it should be faster than largange algorithm, 
 *  but it includes a copy on input data, that could affect CPU performance
 *  Numerically it is more stable that Lagrange interpolation, 
 *  but anyhow it also should not be used for very high polynomial 
 *  degrees (e.g.>20), especially for uniform grid (Runge phenomenon)
 *
 *  - If y-vector is shorter than x-vector, it is assumed to be zero padded
 *  - If y-vector is longer  than x-vector, extra avalues are ignored       
 *  - If x-vector is empty, polynomial is zero
 *
 *  @param  xbegin INPUT start iterator for the sequence of abscissas 
 *  @param  xend   INPUT end   iterator for the sequence of abscissas 
 *  @param  ybegin INPUT start iterator for the sequence of values 
 *  @param  yend   INPUT end   iterator for the sequence of values 
 *  @param  x      INPUT evaluate the polynomial in this point
 *  @param  xvalue INPUT adapter for x-values  
 *  @param  yvalue INPUT adapter for y-values  
 *  @return the value of interpolation polynomial and derivative at point x 
 *  @see https://en.wikipedia.org/wiki/Neville%27s_algorithm
 */
// ============================================================================
template <class XITERATOR, 
          class YITERATOR, 
          class XADAPTER ,
          class YADAPTER >
inline std::pair<double,double>
Ostap::Math::Interpolation::neville2
( XITERATOR    xbegin ,
  XITERATOR    xend   ,
  YITERATOR    ybegin ,
  YITERATOR    yend   ,
  const double x      , 
  XADAPTER     xvalue , // adaptor to get y-value
  YADAPTER     yvalue ) // adaptor to get x-value  
{
  /// 1)special case no data 
  if ( xbegin == xend || ybegin == yend ) { return std::make_pair ( 0                  , 0 ) ; } // RETURN
  const typename std::iterator_traits<XITERATOR>::difference_type NX = std::distance ( xbegin , xend ) ;  
  /// 2) special case: constant function 
  if ( 1 == NX                          ) { return std::make_pair ( yvalue ( *ybegin ) , 0 ) ; } // RETURN
  const typename std::iterator_traits<YITERATOR>::difference_type NY = std::distance ( ybegin , yend ) ;  
  // temporary storage 
  std::vector<double> _y ( NX , 0 ) ;
  std::vector<double> _d ( NX , 0 ) ;
  YITERATOR ylast = yend ;
  if ( NX < NY ) { std::advance ( ylast , NX - NY ) ; }
  std::transform ( ybegin , ylast , _y.begin() , yvalue ) ;
  // simple version of neville 
  return neville ( xbegin , xend , _y.begin() , _d.begin() , x , xvalue ) ;
}
// ============================================================================
/** simple interpolation using Neville's algorithm: 
 *
 *
 *  In general it should be faster than largange algorithm, 
 *  but it modifies input data!
 *
 *  Numerically it is more stable than Lagrange interpolation, 
 *  but anyhow it also should not be used for very high polynomial 
 *  degrees (e.g.>20), especially for uniform grid (Runge phenomenon)
 *
 *  - If y-vector is shorter than x-vector, it is assumed to be zero padded
 *  - If y-vector is longer  than x-vector, extra avalues are ignored       
 *  - If x-vector is empty, polynomial is zero
 *
 *  @attention y-sequence must be "simple", convertible to doubles  
 *  @attention y-sequence is *MODIFIED*
 *
 *  @param  xbegin INPUT  start iterator for the sequence of abscissas 
 *  @param  xend   INPUT  end   iterator for the sequence of abscissas 
 *  @param  ybegin UPDATE start iterator for the sequence of values 
 *  @param  yend   UPDATE end   iterator for the sequence of values 
 *  @param  x      INPUT  evaluate the polynomial in this point
 *  @param  xvalue INPUT  adapter for x-values  
 *  @param  yvalue INPUT  adapter for y-values  
 *  @return the value of  interpolation polynomial at point x
 *  @see https://en.wikipedia.org/wiki/Neville%27s_algorithm
 */
// ============================================================================
template <class XITERATOR, 
          class YITERATOR, 
          class XADAPTOR >
inline double 
Ostap::Math::Interpolation::neville
( XITERATOR    xbegin , 
  XITERATOR    xend   , 
  YITERATOR    ybegin , // NON-CONST 
  const double x      , 
  XADAPTOR     xvalue )
{
  if ( xbegin == xend ) { return 0        ; } // RETURN 
  const typename std::iterator_traits<XITERATOR>::difference_type NX = std::distance ( xbegin , xend ) ;
  if ( 1 == NX        ) { return *ybegin  ; } // RETURN
  //
  for ( unsigned int k = 1 ; k < NX ; ++k ) 
  {
    for ( unsigned int i = 0 ; i < NX - k ; ++i ) 
    {
      const double xj  = xvalue ( *(xbegin +(i + k) ) ) ;
      const double xi  = xvalue ( *(xbegin + i      ) ) ;
      const double yi  = *(ybegin + i      ) ;
      const double yi1 = *(ybegin +(i + 1) ) ;
      *(ybegin+i) = ( ( x - xj ) * yi + ( xi - x ) * yi1 ) /( xi - xj ) ;
    }
  }
  return *ybegin ;
}
// ============================================================================
/** simple interpolation using Neville's algorithm with simultaneous 
 *  estimation of the derivative  
 *
 *  Numerically it is more stable than Lagrange interpolation, 
 *  but anyhow it also should not be used for very high polynomial 
 *  degrees (e.g.>20), especially for uniform grid (Runge phenomenon)
 *
 *  @attention y-sequence must be "simple", convertible to doubles  
 *  @attention y-sequence is *MODIFIED*
 *  @attention d-sequence is *MODIFIED*
 *
 *  @param  xbegin INPUT  start iterator for the sequence of abscissas 
 *  @param  xend   INPUT  end   iterator for the sequence of abscissas 
 *  @param  ybegin UPDATE start iterator for the sequence of values 
 *  @param  x      INPUT  evaluate the polynomial in this point
 *  @param  xvalue INPUT  adapter for x-values  
 *  @param  yvalue INPUT  adapter for y-values  
 *  @return the pair (function,derivative) at point x 
 *  @see https://en.wikipedia.org/wiki/Neville%27s_algorithm
 */
// ============================================================================
template <class XITERATOR, 
          class YITERATOR, 
          class DITERATOR, 
          class XADAPTOR >
inline std::pair<double,double>
Ostap::Math::Interpolation::neville
( XITERATOR    xbegin , 
  XITERATOR    xend   , 
  YITERATOR    ybegin , // NON-const!
  DITERATOR    dbegin , // NON-const!
  const double x      , 
  XADAPTOR     xvalue ) 
{
  if ( xbegin == xend ) { return std::make_pair ( 0       , 0 ) ; } // RETURN 
  const typename std::iterator_traits<XITERATOR>::difference_type NX = std::distance ( xbegin , xend ) ;
  if ( 1 == NX        ) { return std::make_pair ( *ybegin , 0 ) ; } // RETURN
  //
  for ( unsigned int k = 1 ; k < NX ; ++k ) 
  {
    for ( unsigned int i = 0 ; i < NX - k ; ++i ) 
    {
      const double xj  = xvalue ( *(xbegin +(i + k) ) ) ;
      const double xi  = xvalue ( *(xbegin + i      ) ) ;
      const double yi  = *(ybegin + i      ) ;
      const double yi1 = *(ybegin +(i + 1) ) ;
      const double di  = *(dbegin + i      ) ;
      const double di1 = *(dbegin +(i + 1) ) ;
      *(ybegin+i) = ( ( x - xj ) * yi      + ( xi - x ) * yi1        ) / ( xi - xj ) ;
      *(dbegin+i) = ( ( x - xj ) * di + yi + ( xi - x ) * di1 - yi1  ) / ( xi - xj ) ;
    }
  }
  return std::make_pair ( *ybegin , *dbegin ) ;
}
// ============================================================================
// The END 
// ============================================================================
#endif // OSTAP_INTERPOLATION_H
// ============================================================================
