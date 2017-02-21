// $Id$
// ============================================================================
#ifndef OSTAP_HISTOINTERPOLATE_H 
#define OSTAP_HISTOINTERPOLATE_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/ValueWithError.h"
// ============================================================================
// forward declarations 
// ============================================================================
class TH1 ; // from ROOT 
class TH2 ; // from ROOT 
class TH3 ; // from ROOT 
// ============================================================================
/** @file HistoInterpolation.h LHCbMath/HistoInterpolation.h
 *  Collection of primitive utilities for hiostorgam interpoaltion 
 *
 *  Originally developed for Ostap in python
 *  Translated to C++ to get some gain in CPU 
 *
 *  @author Vanya Belyaev
 *  @date   2015-10-12
 *
 *  Version           $Revision$
 *  Last mofidication $Date$
 *                 by $Author$
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    /** @class HistoInterpolation
     *  Collection of primitive utilities for hiostorgam interpoaltion 
     *
     *  Originally developed for Ostap in python
     *  Translated to C++ to get some gain in CPU 
     *
     *  @author Vanya Belyaev
     *  @date   2015-10-12
     */
    class HistoInterpolation 
    {
    public:
      // ======================================================================
      enum Type { Nearest   , 
                  Linear    , 
                  Quadratic , 
                  Cubic     } ;
      // ======================================================================
    public: // 1D interpolation 
      // ======================================================================
      /** linear interpolation between two points 
       *  @param x  (INPUT) the x-value
       *  @param x0 (INPUT) x-coordinate of the first  point 
       *  @param x1 (INPUT) x-coordinate of the second point 
       *  @param y0 (INPUT) y-coordinate of the first  point  \f$ y(x_0) \f$ 
       *  @param y1 (INPUT) y-coordinate of the second point  \f$ y(x_1) \f$  
       *  @return result of linear interpolation \f$ y(x) \f$
       */
      static Ostap::Math::ValueWithError interpolate
        ( const double x  , 
          const double x0 , 
          const double x1 ,
          const Ostap::Math::ValueWithError& y0 , 
          const Ostap::Math::ValueWithError& y1 ) ;
      // ======================================================================
      /** quadratic (parabolic)  interpolation between three points 
       *  @param x  (INPUT) the x-value
       *  @param x0 (INPUT) x-coordinate of the first  point 
       *  @param x1 (INPUT) x-coordinate of the second point 
       *  @param x2 (INPUT) x-coordinate of the third  point 
       *  @param y0 (INPUT) y-coordinate of the first  point \f$ y(x_0) \f$  
       *  @param y1 (INPUT) y-coordinate of the second point \f$ y(x_1) \f$ 
       *  @param y2 (INPUT) x-coordinate of the third  point \f$ y(x_2) \f$ 
       *  @return result of  quadratic (parabolic) interpolation  \f$ y(x) \f$
       */
      static Ostap::Math::ValueWithError interpolate
        ( const double                       x  , 
          const double                       x0 , 
          const double                       x1 ,
          const double                       x2 ,
          const Ostap::Math::ValueWithError& y0 , 
          const Ostap::Math::ValueWithError& y1 , 
          const Ostap::Math::ValueWithError& y2 ) ;
      // ======================================================================
      /** qubic interpolation between four points 
       *  @param x  (INPUT) the x-value
       *  @param x0 (INPUT) x-coordinate of the first  point 
       *  @param x1 (INPUT) x-coordinate of the second point 
       *  @param x2 (INPUT) x-coordinate of the third  point 
       *  @param x3 (INPUT) x-coordinate of the fourth point 
       *  @param y0 (INPUT) y-coordinate of the first  point \f$ y(x_0) \f$  
       *  @param y1 (INPUT) y-coordinate of the second point \f$ y(x_1) \f$ 
       *  @param y2 (INPUT) x-coordinate of the third  point \f$ y(x_2) \f$ 
       *  @param y3 (INPUT) x-coordinate of the third  point \f$ y(x_3) \f$ 
       *  @return result of  quadratic (parabolic) interpolation  \f$ y(x) \f$
       */
      static Ostap::Math::ValueWithError interpolate
        ( const double                       x  , 
          const double                       x0 , 
          const double                       x1 ,
          const double                       x2 ,
          const double                       x3 ,
          const Ostap::Math::ValueWithError& y0 , 
          const Ostap::Math::ValueWithError& y1 , 
          const Ostap::Math::ValueWithError& y2 ,
          const Ostap::Math::ValueWithError& y3 ) ;
      // ======================================================================
    public: // 2D interpolation 
      // ======================================================================
      /** bi-linear interpolation on grid 
       *  @param x   (INPUT) the x-value
       *  @param y   (INPUT) the x-value
       *  @param x0  (INPUT) the first  x-coordinate on the grid 
       *  @param x1  (INPUT) the second x-coordinate on the grid 
       *  @param y0  (INPUT) the first  y-coordinate on the grid 
       *  @param y1  (INPUT) the second y-coordinate on the grid 
       *  @param f00 (INPUT) function value for (x0,y0)
       *  @param f10 (INPUT) function value for (x1,y0)
       *  @param f01 (INPUT) function value for (x0,y1)
       *  @param f11 (INPUT) function value for (x1,y1)
       *  @return result of bi-linear interpolation
       */
      static Ostap::Math::ValueWithError interpolate
        ( const double                       x   , 
          const double                       y   , 
          const double                       x0  , 
          const double                       x1  ,
          const double                       y0  , 
          const double                       y1  ,
          const Ostap::Math::ValueWithError& f00 , 
          const Ostap::Math::ValueWithError& f10 ,
          const Ostap::Math::ValueWithError& f01 , 
          const Ostap::Math::ValueWithError& f11 ) ;
      // ======================================================================
      /** bi-quadratic interpolation on grid 
       *  @param x    (INPUT) the x-value
       *  @param y    (INPUT) the x-value
       *  @param x0   (INPUT) the first  x-coordinate on the grid 
       *  @param x1   (INPUT) the second x-coordinate on the grid 
       *  @param x2   (INPUT) the third  x-coordinate on the grid 
       *  @param y0   (INPUT) the first  y-coordinate on the grid 
       *  @param y1   (INPUT) the second y-coordinate on the grid 
       *  @param y2   (INPUT) the third  y-coordinate on the grid 
       *  @param f00  (INPUT) function value for (x0,y0)
       *  @param f10  (INPUT) function value for (x1,y0)
       *  @param f20  (INPUT) function value for (x2,y0)
       *  @param f01  (INPUT) function value for (x0,y1)
       *  @param f11  (INPUT) function value for (x1,y1)
       *  @param f21  (INPUT) function value for (x2,y1)
       *  @param f02  (INPUT) function value for (x0,y2)
       *  @param f12  (INPUT) function value for (x1,y2)
       *  @param f22  (INPUT) function value for (x2,y2)
       *  @return result of bi-quadrate interpolation
       */
      static Ostap::Math::ValueWithError interpolate
        ( const double                       x   , 
          const double                       y   , 
          const double                       x0  , 
          const double                       x1  ,
          const double                       x2  ,
          const double                       y0  , 
          const double                       y1  ,
          const double                       y2  ,
          const Ostap::Math::ValueWithError& f00 , 
          const Ostap::Math::ValueWithError& f10 ,
          const Ostap::Math::ValueWithError& f20 ,
          const Ostap::Math::ValueWithError& f01 , 
          const Ostap::Math::ValueWithError& f11 ,
          const Ostap::Math::ValueWithError& f21 ,
          const Ostap::Math::ValueWithError& f02 , 
          const Ostap::Math::ValueWithError& f12 ,
          const Ostap::Math::ValueWithError& f22 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /** interpolate 1D histogram 
       *  @param h1          (INPUT) input histogram 
       *  @param x           (INPUT) the x-value 
       *  @param t           (INPUT) interpolation type 
       *  @param edges       (INPUT) use the special treatment of edges ? 
       *  @param extrapolate (INPUT) use extrapolation ? 
       *  @param density     (INPUT) rescale to density? 
       *  If "density" flag is activated, actually   the value of 
       *  density function, defined as a ratio of bin content over bin volume 
       *  is interpolated 
       *  @return value of interpolated function/density
       */
      static Ostap::Math::ValueWithError interpolate_1D
        ( const TH1&   h1                   , 
          const double x                    ,
          const Type   t           = Linear , 
          const bool   edges       = true   , 
          const bool   extrapolate = false  , 
          const bool   density     = false  ) ;
      // ======================================================================
      /** interpolate 2D histogram 
       *  @param h2          (INPUT) input histogram 
       *  @param x           (INPUT) the x-value 
       *  @param y           (INPUT) the y-value 
       *  @param tx          (INPUT) interpolation type in x-direction
       *  @param ty          (INPUT) interpolation type in y-direction
       *  @param edges       (INPUT) use the special treatment of edges ? 
       *  @param extrapolate (INPUT) use extrapolation ? 
       *  @param density     (INPUT) rescale to density? 
       *  If "density" flag is activated, actually   the value of 
       *  density function, defined as a ratio of bin content over bin volume 
       *  is interpolated 
       *  @return value of interpolated function/density
       */
      static Ostap::Math::ValueWithError interpolate_2D 
        ( const TH2&   h1                   , 
          const double x                    ,
          const double y                    ,
          const Type   tx          = Linear , 
          const Type   ty          = Linear ,
          const bool   edges       = true   , 
          const bool   extrapolate = false  , 
          const bool   density     = false  ) ;
      // ======================================================================
      /** interpolate 3D histogram 
       *  @param h3          (INPUT) input histogram 
       *  @param x           (INPUT) the x-value 
       *  @param y           (INPUT) the y-value 
       *  @param z           (INPUT) the z-value 
       *  @param tx          (INPUT) interpolation type in x-direction
       *  @param ty          (INPUT) interpolation type in y-direction
       *  @param tz          (INPUT) interpolation type in z-direction
       *  @param edges       (INPUT) use the special treatment of edges ? 
       *  @param extrapolate (INPUT) use extrapolation ? 
       *  @param density     (INPUT) rescale to density? 
       *  If "density" flag is activated, actually   the value of 
       *  density function, defined as a ratio of bin content over bin volume 
       *  is interpolated 
       *  @return value of interpolated function/density
       */
      static Ostap::Math::ValueWithError interpolate_3D 
        ( const TH3&   h3                   , 
          const double x                    ,
          const double y                    ,
          const double z                    ,
          Type         tx          = Linear , 
          Type         ty          = Linear ,
          Type         tz          = Linear ,
          const bool   edges       = true   , 
          const bool   extrapolate = false  , 
          const bool   density     = false  ) ;
      // ======================================================================
    } ;  
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_HISTOINTERPOLATE_H
// ============================================================================
