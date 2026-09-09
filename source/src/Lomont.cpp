// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <type_traits>

// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Lomont.h"

// ============================================================================
namespace 
{
  // ==========================================================================
  // Compile-time checks for float/int32 and IEEE 754 compliance
  // ==========================================================================
  static_assert( std::numeric_limits<float>::is_specialized &&
                 std::numeric_limits<std::int32_t>::is_specialized &&
                 std::numeric_limits<std::uint32_t>::is_specialized &&
                 sizeof( float ) == sizeof( std::int32_t ) &&
                 sizeof( float ) == sizeof( std::uint32_t ) &&
                 31 == std::numeric_limits<std::int32_t>::digits &&
                 32 == std::numeric_limits<std::uint32_t>::digits &&
                 std::numeric_limits<float>::is_iec559,
                 "FAILED FLOAT/INT32 OR IEEE 754 ASSUMPTIONS" );

  // ==========================================================================
  // Compile-time checks for double/int64 and IEEE 754 compliance
  // ==========================================================================
  static_assert( std::numeric_limits<double>::is_specialized &&
                 std::numeric_limits<std::int64_t>::is_specialized &&
                 std::numeric_limits<std::uint64_t>::is_specialized &&
                 sizeof( double ) == sizeof( std::int64_t ) &&
                 sizeof( double ) == sizeof( std::uint64_t ) &&
                 63 == std::numeric_limits<std::int64_t>::digits &&
                 64 == std::numeric_limits<std::uint64_t>::digits &&
                 std::numeric_limits<double>::is_iec559,
                 "FAILED DOUBLE/INT64 OR IEEE 754 ASSUMPTIONS" );

  // ==========================================================================
  // Safe C++17 bit_cast replacement (avoids undefined behavior from unions)
  // ==========================================================================
  template <typename To, typename From>
  inline To bit_cast( const From& src ) noexcept
  {
    static_assert( sizeof( To ) == sizeof( From ), "Sizes must match" );
    static_assert( std::is_trivially_copyable_v<To> && std::is_trivially_copyable_v<From>,
                   "Types must be trivially copyable" );
    To dst;
    std::memcpy( &dst, &src, sizeof( To ) );
    return dst;
  }

  // Metafunction to map float/double to their corresponding integer types
  template <typename T>
  using int_type_t = std::conditional_t<sizeof( T ) == 4, std::int32_t, std::int64_t>;

  // ==========================================================================
  // Unified template for ULP distance calculation (C++17)
  // ==========================================================================
  template <typename T>
  inline std::intmax_t _distance_impl_( const T a, const T b ) noexcept 
  {
    static_assert( std::numeric_limits<T>::is_iec559, "Type must conform to IEEE 754 (IEC 559)" );

    if ( a == b ) { return 0; }
    if ( !std::isfinite( a ) || !std::isfinite( b ) ) {
      return std::numeric_limits<std::intmax_t>::max();
    }
    if ( a > b ) { return -_distance_impl_( b, a ); }
    if ( !b ) { return _distance_impl_( static_cast<T>( 0 ), -a ); }
    if ( b < 0 ) { return _distance_impl_( -b, -a ); }
    if ( a < 0 && 0 < b ) {
      return -_distance_impl_( static_cast<T>( 0 ), a ) +
             _distance_impl_( static_cast<T>( 0 ), b );
    }

    using IntT = int_type_t<T>;
    const auto ai = bit_cast<IntT>( a );
    const auto bi = bit_cast<IntT>( b );
    return static_cast<std::intmax_t>( bi ) - static_cast<std::intmax_t>( ai );
  }

  inline std::intmax_t _distance_float_( const float a, const float b ) noexcept 
  { return _distance_impl_( a, b ); }

  inline std::intmax_t _distance_double_( const double a, const double b ) noexcept 
  { return _distance_impl_( a, b ); }

  // ==========================================================================
  // Unified template for ULP stepping (C++17)
  // ==========================================================================
  template <typename T>
  inline T _next_impl_( const T a, const int_type_t<T> ulps ) noexcept 
  {
    static_assert( std::numeric_limits<T>::is_iec559, "Type must conform to IEEE 754 (IEC 559)" );

    if ( 0 == ulps || !std::isfinite( a ) ) { return a; }
    if ( 0 > a ) { return -_next_impl_( -a, -ulps ); }

    if ( 0 > ulps ) {
      const auto d = a ? ( _distance_impl_( static_cast<T>( 0 ), a ) + ulps ) : ulps;
      if ( d < 0 ) { return -_next_impl_( static_cast<T>( 0 ), -d ); }
    }

    using IntT = int_type_t<T>;
    auto ai = bit_cast<IntT>( a );
    ai += ulps;
    return bit_cast<T>( ai );
  }

  inline float _next_float_( const float a, const std::int32_t ulps ) noexcept 
  { return _next_impl_( a, ulps ); }

  inline double _next_double_( const double a, const std::int64_t ulps ) noexcept 
  { return _next_impl_( a, ulps ); }

  // ==========================================================================
  inline bool _compare_float_ 
  ( const float a, 
    const float b, 
    const unsigned short maxULPs ) noexcept 
  {
    const std::intmax_t diff = _distance_float_( a, b );
    return std::abs( diff ) <= maxULPs;
  }

  // ==========================================================================
  inline bool _compare_double_
  ( const double a, 
    const double b, 
    const unsigned int maxULPs ) noexcept 
  {
    const std::intmax_t diff = _distance_double_( a, b );
    return std::abs( diff ) <= maxULPs;
  }

} // namespace

// ============================================================================
// Implementation of public functions declared in Ostap/Lomont.h
// ============================================================================

bool Ostap::Math::Lomont::compare_float
( const float a, 
  const float b, 
  const unsigned short maxULPs ) 
{ return a == b || _compare_float_( a, b, maxULPs ); }

float Ostap::Math::Lomont::next_float
( const float a, 
  const short ulps ) 
{ return _next_float_( a, ulps ); }

float Ostap::Math::Lomont::prev_float
( const float a, 
  const short ulps )
{ return next_float( a, -ulps ); }

bool Ostap::Math::Lomont::compare_double
( const double a, 
  const double b, 
  const unsigned int maxULPs ) 
{ return a == b || _compare_double_( a, b, maxULPs ); }

double Ostap::Math::Lomont::next_double
( const double a, 
  const short ulps ) 
{ return _next_double_( a, ulps ); }

double Ostap::Math::Lomont::prev_double
( const double a, 
  const short ulps ) 
{ return next_double( a, -ulps ); }

std::intmax_t Ostap::Math::Lomont::ulps_distance_float
( const float a, 
  const float b ) 
{ return _distance_float_( a, b ); }

std::intmax_t Ostap::Math::Lomont::ulps_distance_double
( const double a, 
  const double b ) 
{ return _distance_double_( a, b ); }

std::int32_t Ostap::Math::Lomont::float2int
( const float v ) 
{ return bit_cast<std::int32_t>( v ); }

float Ostap::Math::Lomont::int2float
( const std::int32_t i ) 
{ return bit_cast<float>( i ); }

std::int64_t Ostap::Math::Lomont::double2int
( const double v ) 
{ return bit_cast<std::int64_t>( v ); }

double Ostap::Math::Lomont::int2double
( const std::int64_t i ) 
{ return bit_cast<double>( i ); }

// ============================================================================
// The END
// ============================================================================