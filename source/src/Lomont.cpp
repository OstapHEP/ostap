// ============================================================================
// Include files
// ============================================================================
// STD & STL 
// ============================================================================
#include <cstdint>
#include <cmath>
#include <limits>
#include <cassert>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Lomont.h"
// ============================================================================
namespace 
{
  // ==========================================================================
  static_assert( std::numeric_limits<float>          ::is_specialized &&
                 std::numeric_limits<std::int32_t>   ::is_specialized && 
                 std::numeric_limits<std::uint32_t>  ::is_specialized &&
		 sizeof ( float ) == sizeof ( std::int32_t  )         &&
		 sizeof ( float ) == sizeof ( std::uint32_t )         &&  
                 31 == std::numeric_limits<std::int32_t>::digits      && 
                 32 == std::numeric_limits<std::uint32_t>::digits     , 
		 "FAILED FLOAT/INT32 ASSUMPTIONS" ) ;
  // ===========================================================================
  static_assert( std::numeric_limits<double>         ::is_specialized &&
                 std::numeric_limits<std::int64_t>   ::is_specialized && 
                 std::numeric_limits<std::uint64_t>  ::is_specialized &&
		 sizeof ( double ) == sizeof ( std::int64_t  )        &&
		 sizeof ( double ) == sizeof ( std::uint64_t )        &&  
                 63 == std::numeric_limits<std::int64_t>::digits      && 
                 64 == std::numeric_limits<std::uint64_t>::digits     , 
		 "FAILED DOUBLE/INT64 ASSUMPTIONS" ) ;
  // ===========================================================================
  /** @struct Cast_F 
   *  Helper structure to perfrom "cast" between int and float 
   *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
   *  @date 2009-05-23
   */
  struct Cast_F
  {
    // ===========================================================================
  public:
    // ========================================================================
    /// int -> float 
    inline float        i2f ( const std::int32_t i ) { m_f.i = i ; return m_f.f ; } // int   -> float
    /// float -> in
    inline std::int32_t f2i ( const float        f ) { m_f.f = f ; return m_f.i ; } // float -> int 
    // ========================================================================
  private:
    // ========================================================================
    /// Helper union to avoid the reinterpret cast for floats 
    union Float_U                     // Helper union to avoid reinterpret cast 
    {
      float        f ;  // float value 
      std::int32_t i ;  // int   value 
    } ;
    // ========================================================================
  private:
    // ========================================================================
    /// the helper union
    Float_U m_f ;                                           // the helper union
    // ========================================================================
  } ;
  // ==========================================================================
  /** @struct Cast_D 
   *  Helper structure to perfrom "cast" between long and double  
   *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
   *  @date 2009-05-23
   */
  struct Cast_D
  {
    // ========================================================================
  public:
    // ========================================================================
    /// long   -> double 
    inline double       l2d ( const std::int64_t l ) { m_d.l = l ; return m_d.d ; } // long   -> double
    /// double -> long 
    inline std::int64_t d2l ( const double       d ) { m_d.d = d ; return m_d.l ; } // double -> long
    // ========================================================================
  private:
    // ========================================================================
    /// Helper union to avoid the reinterpret cast for floats 
    union Double_U                     // Helper union to avoid reinterpret cast 
    {
      double        d ; // double value 
      std::int64_t  l ; // long   value 
    } ;
    // ========================================================================
  private:
    // ========================================================================
    /// the helper union
    Double_U m_d ;                                          // the helper union
    // ========================================================================
  } ;
  // ==========================================================================
  // kind of "distance" between two floats
  inline std::intmax_t _distance_float_
  ( const float a  , 
    const float b  ) 
  {
    //
    if      ( a == b               ) { return 0 ; }
    else if ( !std::isfinite ( a ) ) { return std::numeric_limits<std::intmax_t>::max() ; } 
    else if ( !std::isfinite ( b ) ) { return std::numeric_limits<std::intmax_t>::max() ; } 
    else if ( a >  b               ) { return -_distance_float_ (  b ,  a ) ; }
    //
    if      ( !b                   ) { return  _distance_float_ (  0 , -a ) ; }
    // both numbers are negative:
    else if ( b < 0                ) { return  _distance_float_ ( -b , -a ) ; }
    // both numbers have different  signs: 
    else if ( a < 0 && 0 < b       ) { return -_distance_float_ (  0 ,  a ) + _distance_float_ ( 0 , b ) ; }
    //
    Cast_F caster{} ;
    //
    // const int ai   = caster.f2i ( af ) ;
    // const int bi   = caster.f2i ( bf ) ;
    // const int test = (((unsigned int)(ai^bi))>>31)-1;
    // return ((( const_min_int - ai ) & (~test)) | ( ai& test )) - bi ;
    //
    const std::intmax_t ai = caster.f2i ( a ) ;
    const std::intmax_t bi = caster.f2i ( b ) ;
    return bi - ai ; 
    // ========================================================================
  }
  // ==========================================================================
  // kind of "distance" between two doubles
  inline std::intmax_t
  _distance_double_
  ( const double a , 
    const double b ) 
  {
    //
    if      (  a == b               ) { return 0 ; }
    else if (  !std::isfinite ( a ) ) { return std::numeric_limits<std::intmax_t>::max() ; }    
    else if (  !std::isfinite ( b ) ) { return std::numeric_limits<std::intmax_t>::max() ; } 
    else if (  a >  b               ) { return - _distance_double_ (  b ,  a ) ; }
    //
    if      ( !b                    ) { return   _distance_double_ (  0 , -a ) ; }
    // both numbers are negative:
    else if (  b <  0               ) { return   _distance_double_ ( -b , -a ) ; }
    // both numbers have different  sign: 
    else if (  a <  0 && 0 <  b     )  { return -_distance_double_ (  0 ,  a ) + _distance_double_ ( 0 , b ) ; }
    //
    Cast_D caster{} ;
    //
    // const Long ai   = caster.d2l ( af ) ;
    // const Long bi   = caster.d2l ( bf ) ;
    // const Long test = (((ULong)(ai^bi))>>63)-1;
    // return ((( const_min_long - ai ) & (~test)) | ( ai& test )) - bi ;
    //
    const std::intmax_t ai = caster.d2l ( a ) ;
    const std::intmax_t bi = caster.d2l ( b ) ;
    return bi - ai ; 
  }
  // ==========================================================================
  inline bool _compare_float_ 
  ( const float          a       , 
    const float          b       , 
    const unsigned short maxULPs ) 
  {
    const std::intmax_t diff = _distance_float_ ( a , b ) ;
    return std::abs ( diff ) <= maxULPs ; 
  }
  // ==========================================================================
  bool _compare_double_
  ( const double       a       , 
    const double       b       , 
    const unsigned int maxULPs ) 
  {
    const std::intmax_t diff = _distance_double_ ( a , b ) ;
    return std::abs ( diff ) <= maxULPs ; 
  }
  // ==========================================================================
  // next  float
  inline float _next_float_
  ( const float a           ,
    const std::int32_t ulps ) 
  {
    if      (  0 == ulps             ) { return a  ; }
    else if (  !std::isfinite ( a  ) ) { return a  ; }
    else if (  0 > a                 ) { return -_next_float_ ( -a , -ulps ) ; }
    //
    if ( 0 > ulps ) 
    {
      const std::int32_t d = _distance_float_ ( a , 0 ) + ulps ;
      if ( d < 0 ) { return _next_float_ ( 0  , d ) ; }
    }
    //    
    Cast_F caster{} ;
    std::int32_t ai  = caster.f2i ( a ) ;
    ai              += ulps ;
    return caster.i2f ( ai );
  }
  // ============================================================================
  // next  double
  inline double _next_double_
  ( const double       a    ,
    const std::int64_t ulps ) 
  {
    //
    if      ( 0 == ulps            ) { return a ; }
    else if ( !std::isfinite ( a ) ) { return a ; }
    else if ( 0 > a                ) { return - _next_double_ ( -a , -ulps ) ; }
    //
    if ( 0 > ulps ) 
    {
      const std::int64_t d = a ?  (_distance_double_ ( 0 , a  ) + ulps ) : ulps ;
      if ( d < 0 ) { return - _next_double_ ( 0 , -d ) ; }
    }    
    //
    Cast_D caster {} ;    
    std::int64_t al  = caster.d2l ( a ) ;
    al              += ulps ;
    return caster.l2d ( al );
  }
  // ============================================================================
} // end of anonymous namespace 
// ============================================================================
/*  equality comparison of float numbers using as the metric the maximal 
 *  number of Units in the Last Place (ULP).
 *  It is a slightly modified version of very efficient implementation 
 *  of the initial Bruce Dawson's algorithm by Chris Lomont.
 *
 *  @see www.lomont.org 
 *  @see http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
 *
 *  Lomont claims the algorithm is factor 2-10 more efficient 
 *  with respect to  Knuth's algorithm from comparisons of floating number 
 *  using the relative precision.
 *
 *  The effective relative difference depends on the choice of 
 *   <c>maxULPS</c>:
 *  - For the case of maxULPs=1, (of cource it is totally unphysical case!!!)
 *  the effective relative precision r = |a-b|/(|a|+|b|)is 
 *  between 3.5e-8 and 5.5e-8 for |a|,|b|>1.e-37, and 
 *  then it quickly goes to ~1 
 *  - For the case of maxULPS=10 
 *  the effective relative precision is 
 *  between 3e-8 and 6e-7 for |a|,|b|>1.e-37, and 
 *  then it quickly goes to ~1 
 *  - For the case of maxULPS=100 
 *  the effective relative precision is 
 *  around ~6e-6 for |a|,|b|>1.e-37, and 
 *  then it quickly goes to ~1 
 *  - For the case of maxULPS=1000 
 *  the effective relative precision is 
 *  around ~6e-5 for |a|,|b|>1.e-37, and 
 *  then it quickly goes to ~1 
 *  
 *  @param  af the first number 
 *  @param  bf the second number 
 *  @param  maxULPS the maximal metric deciation in the terms of 
 *                 maximal number of units in the last place
 *  @author Vanya BELYAEV  Ivan.Belyaev@nikhef.nl
 */
// ============================================================================
bool Ostap::Math::Lomont::compare_float 
( const float          a       , 
  const float          b       , 
  const unsigned short maxULPs ) 
{
  // ==========================================================================
  static_assert ( std::numeric_limits<float>          ::is_specialized &&
		  std::numeric_limits<std::int32_t>   ::is_specialized && 
		  std::numeric_limits<std::uint32_t>  ::is_specialized &&
		  sizeof ( float ) == sizeof ( std::int32_t  )         &&
		  sizeof ( float ) == sizeof ( std::uint32_t )         &&  
		  31 == std::numeric_limits<std::int32_t>::digits      && 
		  32 == std::numeric_limits<std::uint32_t>::digits     , 
		  "FAILED FLOAT/INT32 ASSUMPTIONS" ) ;
  // ===========================================================================
  return a == b || _compare_float_ ( a , b , maxULPs ) ;
  // ==========================================================================
}
// ============================================================================
/*  get the floating number that representation 
 *  is different with respect  to the argument for 
 *  the certain number of "Units in the Last Position".
 *  For ulps=1, it is just next float number, for ulps=-1 is is the 
 *  previous one.
 *
 *  This routine is very convinient to test the parameter maxULPS for
 *  the routine Gaudi::Math::lomont_compare 
 *
 *  @see Gaudi::Math::lomont_compare
 *  @param af the reference number 
 *  @param ulps the bias 
 *  @return the biased float number (on distance "ulps")
 *  @author Vanya BELYAEV  Ivan.Belyaev@nikhef.nl
 *  @date 2008-11-08
 */  
// ============================================================================
float Ostap::Math::Lomont::next_float
( const float a    ,
  const short ulps ) 
{
  // ==========================================================================
  return _next_float_ ( a , ulps ) ;
  // ==========================================================================
}
// ============================================================================
float Ostap::Math::Lomont::prev_float
( const float a    ,
  const short ulps ) { return next_float ( a , -ulps ) ; }
// ============================================================================
/*  equality comparison of float numbers using as the metric the maximal 
 *  number of Units in the Last Place (ULP).
 *  It is a slightly modified version of very efficient implementation 
 *  of the initial Bruce Dawson's algorithm by Chris Lomont.
 *
 *  @see www.lomont.org 
 *  @see http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
 *
 *  C.Lomont claims the algorithm is factor 2-10 more efficient 
 *  with respect to  Knuth's algorithm from comparisons of floating number 
 *  using the relative precision.
 *
 *  The effective relative difference depends on the choice of 
 *   <c>maxULPS</c>:
 *  - For the case of maxULPs=1, (of cource it is totally unphysical case!!!)
 *  the effective relative precision r = |a-b|/(|a|+|b|)is 
 *  ~6e-16 for |a|,|b|>1.e-304, and 
 *  then it quickly goes to ~1 
 *  
 *  @param  af the first number 
 *  @param  bf the second number 
 *  @param  maxULPS the maximal metric deciation in the terms of 
 *                 maximal number of units in the last place
 *  @author Vanya BELYAEV  Ivan.Belyaev@nikhef.nl
 *  @date 2008-11-08
 */
// ============================================================================
bool Ostap::Math::Lomont::compare_double 
( const double       a       , 
  const double       b       , 
  const unsigned int maxULPs ) 
{
  // ==========================================================================
  static_assert( std::numeric_limits<double>         ::is_specialized &&
                 std::numeric_limits<std::int64_t>   ::is_specialized && 
                 std::numeric_limits<std::uint64_t>  ::is_specialized &&
		 sizeof ( double ) == sizeof ( std::int64_t  )        &&
		 sizeof ( double ) == sizeof ( std::uint64_t )        &&  
                 63 == std::numeric_limits<std::int64_t>::digits      && 
                 64 == std::numeric_limits<std::uint64_t>::digits     , 
		 "FAILED DOUBLE/INT64 ASSUMPTIONS" ) ;
  // ===========================================================================
  return a == b || _compare_double_ ( a , b , maxULPs ) ;
  // ===========================================================================  
}
// ============================================================================
/*  Get the floating number that representation 
 *  is different with respect  to the argument for 
 *  the certain number of "Units in the Last Position".
 *  For ulps=1, it is just next float number, for ulps=-1 is is the 
 *  previous one.
 *
 *  This routine is very convinient to test the parameter maxULPS for
 *  the routine LHCb::Math::lomont_compare_float 
 *
 *  @see Gaudi::Math::lomont_compare
 *  @param ad the reference number 
 *  @param ulps the bias 
 *  @return the biased float number (on distance "ulps")
 *  @author Vanya BELYAEV  Ivan.Belyaev@nikhef.nl
 *  @date 2008-11-08
 */  
// ============================================================================
double Ostap::Math::Lomont::next_double
( const double a    ,
  const short  ulps ) 
{
  // ==========================================================================
  return _next_double_ (  a  , ulps ) ;
  // ==========================================================================
}
// ============================================================================
double Ostap::Math::Lomont::prev_double
( const double a    ,
  const short  ulps ) { return next_double ( a , - ulps ) ; } 
// ============================================================================
/*  "distance" in ULPS between two float values 
 *   @param a (INPUT) the first  number 
 *   @param b (INPUT) the second number 
 *   @param "distance" in ULPs
 */
// ============================================================================
std::intmax_t Ostap::Math::Lomont::ulps_distance_float
( const float  a ,
  const float  b ) 
{ return _distance_float_  ( a , b ) ; }
// ============================================================================
/*  "distance" in ULPS between two double values 
 *   @param a (INPUT) the first  number 
 *   @param b (INPUT) the second number 
 *   @param "distance" in ULPs
 */
// ============================================================================
std::intmax_t Ostap::Math::Lomont::ulps_distance_double
( const double a ,
  const double b )
{ return _distance_double_ ( a , b ) ; }
// ============================================================================
// explicit "cast" of float to int32 
// ============================================================================
std::int32_t Ostap::Math::Lomont::float2int
( const float        v )
{
  Cast_F caster {} ;
  return caster.f2i ( v ) ;
}
// ============================================================================
// explicit "cast" of int   to float 
// ============================================================================
float Ostap::Math::Lomont::int2float
( const std::int32_t i ) 
{
  Cast_F caster {} ;
  return caster.i2f ( i ) ;
}
// ============================================================================
// explicit "cast" of double to int64
// ============================================================================
std::int64_t Ostap::Math::Lomont::double2int
( const double       v )
{
  Cast_D caster {} ;
  return caster.d2l ( v ) ;
}
// ============================================================================
// explicit "cast" of int   to float 
// ============================================================================
double Ostap::Math::Lomont::int2double
( const std::int64_t i )
{
  Cast_D caster {} ;
  return caster.l2d ( i ) ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================
  
