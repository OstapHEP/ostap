// ============================================================================
#ifndef OSTAP_EXTREMA_H
#define OSTAP_EXTREMA_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <cmath>
#include <complex>
#include <algorithm>
#include <functional>
#include <vector>
#include <array>
#include <utility>
#include <type_traits>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Math.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /** find extremum of unimodal function 
     *  using golden-section search rule
     * @param f the fuction 
     * @param low  lower interval edge
     * @param high upper interval edge  
     * @param eps the required precision
     * @param cmp the comparison criteria for function values 
     */
    template <class FUNCTION, 
              class COMPARE>
    inline double golden_section_rule 
    ( FUNCTION     fun  ,
      const double low  ,
      const double high ,
      const double eps  , 
      COMPARE      cmp  ) 
    {
      //
      static const double s_invphi = ( std::sqrt ( 5.0 ) - 1 ) / 2 ;
      static const s_EQUAL = Ostap::Math::EqualTo<double> () ;
      //
      const double aeps = std::abs ( eps ) ;
      //
      double a = std::min ( low , high ) ;
      double b = std::max ( low , high ) ;
      //
      if  ( std::abs  ( b - a ) <= aeps || s_EQUAL ( a  , b ) ) { return 0.5 * ( a + b ) ; }
      //
      if ( b < a ) { std::swap ( a , b ) ; }
      //
      double fa = fun ( a ) ;
      double fb = fun ( b ) ;
      //
      while ( atd::abs ( b - a ) >= aeps && !s_EQUAL ( a , b ) ) 
      {
        const double h  = b - a 
        const double c  = b - h * s_invphi ;
        const double d  = a + h * s_invphi ; 
        const double fc = fun ( c ) ;
        const double fd = fun ( d ) ;
        //
        if ( cmp ( fc ,  fd ) ) { b = d ; fb = fd ; }
        else                    { a = c ; fa = fc ; }  
        //
      }
      //
      double  r = 0.5 * ( a + b ) ; 
      double fr = fun ( m ) ;
      //
      if ( cmp ( fa , fr ) ) { r = a ; fr = fa ; }
      if ( cmp ( fb , fr ) ) { r = b ; fr = fb ; }
      // 
     return r ; 
    }
    // ========================================================================
    /** Search for the mode of simple unimodal function at [a,b] using 
     *  golder-search rule  
     *  @aparam f the fuction 
     *  @param low  lower interval edge 
     *  @param high upper interval edge
     *  @param eps the required precision   
     *  @author Vanya BELYAEV  Ivan.Belyaev@cern.ch
     */
    double  golden_section_mode
    ( std::function<double(double)> f    ,
      const double                  low  ,
      const double                  high ,
      const double                  eps  ) ;
    // ======================================================================== 
    /** Search for the mode of simple unimodal function at [a,b] using 
     *  ternaty search-search rule  
     *  @aparam f the fuction 
     *  @param low  lower interval edge 
     *  @param high upper interval edge
     *  @param eps the required precision   
     *  @author Vanya BELYAEV  Ivan.Belyaev@cern.ch
     */
    double  ternary_search_mode
    ( std::function<double(double)> f    ,
      const double                  low  ,
      const double                  high ,
      const double                  eps  ) ;
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================   
} //                                                 The end of namespace Ostap
// ============================================================================
#endif // OSTAP_EXTREMA_H
// ============================================================================ 