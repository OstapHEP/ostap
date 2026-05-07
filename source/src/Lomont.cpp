// ============================================================================
// Include files
// ============================================================================
// STD & STL 
// ============================================================================
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
  // prerequisites for "correct" Float 
  static_assert( std::numeric_limits<float>          ::is_specialized &&
                 std::numeric_limits<int>            ::is_specialized && 
                 std::numeric_limits<unsigned int>   ::is_specialized &&
                 sizeof(float)==sizeof(int)                           && 
                 sizeof(float)==sizeof(unsigned int)                  &&
                 32 == std::numeric_limits<unsigned int>::digits , "FAILED ASSUMPTIONS") ;
  // ==========================================================================
  // define proper double 
  // ==========================================================================
  template <bool I> 
  struct __Longs ;
  template <>
  struct __Longs<true>
  {
    typedef long               Long   ;
    typedef unsigned long      ULong  ;
  } ;
  template <>
  struct __Longs<false>
  {
    typedef long long               Long   ;
    typedef unsigned long long      ULong  ;
  } ;
  struct _Longs : public __Longs<sizeof(double)==sizeof(long)>{};
  /// get the final types:
  typedef _Longs::Long  Long  ;
  typedef _Longs::ULong ULong ;
  // ==========================================================================
  ///  minimal int 
  const int  const_min_int  = std::numeric_limits<int>::min  () ;
  ///  minimal Long
  const Long const_min_long = std::numeric_limits<Long>::min () ;
  // ==========================================================================
  /// the final check
  static_assert( std::numeric_limits<double>  ::is_specialized &&
                 std::numeric_limits<Long>    ::is_specialized && 
                 std::numeric_limits<ULong>   ::is_specialized &&
                 std::numeric_limits<ULong>   ::is_specialized &&
                 std::numeric_limits<Long>    ::is_specialized &&
                 sizeof(double)==sizeof(Long)                  && 
                 sizeof(double)==sizeof(ULong)                 &&
                 64 == std::numeric_limits<ULong>::digits      , "FAILED ASSUMPTIONS") ;
  // ==========================================================================

  // ==========================================================================
  /** @struct Cast_F 
   *  Helper structure to perfrom "cast" between int and float 
   *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
   *  @date 2009-05-23
   */
  struct Cast_F
  {
  public :
    // ========================================================================
    // prerequisites:
    static_assert( std::numeric_limits<float>          ::is_specialized &&
                   std::numeric_limits<int>            ::is_specialized && 
                   std::numeric_limits<unsigned int>   ::is_specialized &&
                   sizeof(float)==sizeof(int)                           && 
                   32 == std::numeric_limits<unsigned int>::digits , "FAILED ASSUMPTIONS") ;
    // ========================================================================
  public:
    // ========================================================================
    /// int -> float 
    float i2f ( const int   i ) { m_f.i = i ; return m_f.f ; } // int   -> float
    /// float -> in
    int   f2i ( const float f ) { m_f.f = f ; return m_f.i ; } // float -> int 
    // ========================================================================
  private:
    // ========================================================================
    /// Helper union to avoid reinterpret cast for floats 
    union Float_U                     // Helper union to avoid reinterpret cast 
    {
      float f ;  // float value 
      int   i ;  // int   value 
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
    /// the final check
    static_assert( std::numeric_limits<double>  ::is_specialized &&
                   std::numeric_limits<Long>    ::is_specialized && 
                   std::numeric_limits<ULong>   ::is_specialized &&
                   std::numeric_limits<ULong>   ::is_specialized &&
                   std::numeric_limits<Long>    ::is_specialized &&
                   sizeof(double)==sizeof(Long)                  && 
                   sizeof(double)==sizeof(ULong)                 &&
                   64 == std::numeric_limits<ULong>::digits      , "FAILED ASSUMPTIONS") ;
    // ========================================================================
  public:
    // ========================================================================
    /// long   -> double 
    double l2d ( const Long   l ) { m_d.l = l ; return m_d.d ; } // long   -> double
    /// double -> long 
    Long   d2l ( const double d ) { m_d.d = d ; return m_d.l ; } // double -> long
    // ========================================================================
  private:
    // ========================================================================
    /// Helper union to avoid reinterpret cast for floats 
    union Double_U                     // Helper union to avoid reinterpret cast 
    {
      double d ; // double value 
      Long   l ; // long   value 
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
  inline int _distance_float_ ( const float af , 
                                const float bf ) 
  {
    //
    if      ( af == bf    ) { return 0 ; }
    else if ( af >  bf    ) { return -_distance_float_ ( bf , af ) ; }
    //  
    // both numbers are negative:
    if      ( bf < 0      ) { return _distance_float_ ( -bf , -af ) ; }
    // both numbers have differrent  signs: 
    else if ( af < 0 && 0 < bf ) 
    {
      return 
        _distance_float_ ( af   , 0.0f ) + 
        _distance_float_ ( 0.0f ,   bf ) ;
    }
    //
    Cast_F caster{} ;
    //
    const int ai   = caster.f2i ( af ) ;
    const int bi   = caster.f2i ( bf ) ;
    const int test = (((unsigned int)(ai^bi))>>31)-1;
    //
    return ((( const_min_int - ai ) & (~test)) | ( ai& test )) - bi ;
  }
  // ==========================================================================
  // kind of "distance" between two doubles
  inline Long  _distance_double_ ( const double af , 
                                   const double bf ) 
  {
    //
    if      ( af == bf    ) { return 0 ; }
    else if ( af >  bf    ) { return - _distance_double_ (  bf ,  af ) ; }
    //  
    // both numbers are negative:
    if      ( bf < 0      ) { return   _distance_double_ ( -bf , -af ) ; }
    // both numbers have differrent  signs: 
    else if ( af < 0 && 0 < bf ) 
    {
      return
        _distance_double_ ( af   , 0.0l ) + 
        _distance_double_ ( 0.0l ,   bf ) ;
    }
    //
    Cast_D caster{} ;
    //
    const Long ai   = caster.d2l ( af ) ;
    const Long bi   = caster.d2l ( bf ) ;
    const Long test = (((ULong)(ai^bi))>>63)-1;
    //
    return ((( const_min_long - ai ) & (~test)) | ( ai& test )) - bi ;
  }
  // ==========================================================================
  inline bool _compare_float_ 
  ( const float          af      , 
    const float          bf      , 
    const unsigned short maxULPs ) 
  {
    //
    const int diff     = _distance_float_ ( af  , bf ) ;
    //
    const int maxDiff_ = maxULPs ;
    //
    const int v1       = maxDiff_ + diff ;
    const int v2       = maxDiff_ - diff ;
    //
    return 0 <= ( v1 | v2 ) ;
  }
  // ==========================================================================
  bool _compare_double_
  ( const double       af      , 
    const double       bf      , 
    const unsigned int maxULPs ) 
  {
    // ==========================================================================
    //
    const Long diff     = _distance_double_ ( af  , bf ) ;
    //
    const Long maxDiff_ = maxULPs ;
    //
    const Long v1       = maxDiff_ + diff ;
    const Long v2       = maxDiff_ - diff ;
    //
    return 0 <= ( v1 | v2 ) ;
  }
  // ==========================================================================
  // next  float
  inline float _next_float_ ( const float af , const short ulps ) 
  {
    if      ( 0 == ulps ) { return af ; }
    else if ( 0 > af    ) { return -_next_float_ ( -af , -ulps ) ; }
    //
    if ( 0 > ulps ) 
    {
      const int d = _distance_float_ ( af , 0.0f ) + ulps ;
      if (  d < 0 ) { return -_next_float_ ( 0.0f , -d ) ; }
    }
    //
    Cast_F caster{} ;
    int ai  = caster.f2i ( af ) ;
    ai     += ulps ;
    return caster.i2f ( ai );
  }
  // ============================================================================
  // next  double
  inline double _next_double_ ( const double ad , const short ulps ) 
  {
    if      ( 0 == ulps ) { return ad ; }
    else if ( 0 > ad    ) { return -_next_double_ ( -ad , -ulps ) ; }
    //
    if ( 0 > ulps ) 
    {
      const Long d = _distance_double_ ( ad , 0 ) + ulps ;
      if (  d < 0 ) { return -_next_double_ ( 0 , -d ) ; }
    }
    //
    Cast_D caster{} ;
    Long al  = caster.d2l ( ad ) ;
    al      += ulps ;
    return caster.l2d ( al );
  }
  // ============================================================================
} // end of anonymous namespace 
// ============================================================================
/*  equality comparion of float numbers using as the metric the maximal 
 *  number of Units in the Last Place (ULP).
 *  It is a slightly modified version of very efficient implementation 
 *  of the initial Bruce Dawson's algorithm by Chris Lomont.
 *
 *  @see www.lomont.org 
 *  @see http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
 *
 *  Lomont claims the algorithm is factor 2-10 more efficient 
 *  with respect to  Knuth's algorithm fomr comparions of floating number 
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
bool Ostap::Math::lomont_compare_float 
( const float          af      , 
  const float          bf      , 
  const unsigned short maxULPs ) 
{
  // ==========================================================================
  // prerequisites:
  static_assert( std::numeric_limits<float>        ::is_specialized &&
                 std::numeric_limits<int>          ::is_specialized && 
                 std::numeric_limits<unsigned int> ::is_specialized &&
                 sizeof(float)==sizeof(int)                         && 
                 sizeof(float)==sizeof(unsigned int)                &&
                 32 == std::numeric_limits<unsigned int>::digits    , "FAILED ASSUMPTIONS") ;
  // ==========================================================================
  return _compare_float_ ( af ,   bf , maxULPs ) ;
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
float Ostap::Math::next_float ( const float af , const short ulps ) 
{
  /// the final check
  static_assert( std::numeric_limits<float> ::is_specialized &&
                 std::numeric_limits<int>   ::is_specialized && 
                 sizeof(float)==sizeof(int) , "FAILED ASSUMPTIONS") ;
  // ==========================================================================
  return _next_float_ ( af , ulps ) ;
  // ==========================================================================
}
// ============================================================================
float Ostap::Math::prev_float ( const float af , const short ulps ) 
{ return next_float ( af , -ulps ) ; }
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
 *  with respect to  Knuth's algorithm fomr comparions of floating number 
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
bool Ostap::Math::lomont_compare_double 
( const double       af      , 
  const double       bf      , 
  const unsigned int maxULPs ) 
{
  // ==========================================================================
  /// the final check
  static_assert( std::numeric_limits<double>  ::is_specialized &&
                 std::numeric_limits<Long>    ::is_specialized && 
                 std::numeric_limits<ULong>   ::is_specialized &&
                 std::numeric_limits<ULong>   ::is_specialized &&
                 std::numeric_limits<Long>    ::is_specialized &&
                 sizeof(double)==sizeof(Long)                  && 
                 sizeof(double)==sizeof(ULong)                 &&
                 64 == std::numeric_limits<ULong>::digits      , "FAILED ASSUMPTIONS") ;
  // ==========================================================================
  return _compare_double_ ( af ,   bf , maxULPs ) ;
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
double Ostap::Math::next_double ( const double ad , const short ulps ) 
{
  // ==========================================================================
  /// the final check
  static_assert( std::numeric_limits<double> ::is_specialized &&
                 std::numeric_limits<Long>   ::is_specialized && 
                 sizeof(double)==sizeof(Long) , "FAILED ASSUMPTIONS") ;
  // ==========================================================================
  return _next_double_ (  ad , ulps ) ;
  // ==========================================================================
}
// ============================================================================
double Ostap::Math::prev_double ( const double ad , const short ulps )
{  return next_double ( ad, - ulps ) ; } 
// ============================================================================
/*  "distance" in ULPS between two float values 
 *   @param a (INPUT) the first  number 
 *   @param b (INPUT) the second number 
 *   @param "distance" in ULPs
 */
// ============================================================================
long Ostap::Math::ulps_distance_float  ( const float  a , const float  b ) 
{ return _distance_float_  ( a , b ) ; }
// ============================================================================
/*  "distance" in ULPS between two double values 
 *   @param a (INPUT) the first  number 
 *   @param b (INPUT) the second number 
 *   @param "distance" in ULPs
 */
// ============================================================================
long Ostap::Math::ulps_distance_double ( const double a , const double b )
{ return _distance_double_ ( a , b ) ; }
// ============================================================================
//                                                                      The END 
// ============================================================================
  
