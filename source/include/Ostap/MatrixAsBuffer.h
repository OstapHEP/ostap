// ============================================================================
#ifndef OSTAP_MATRIXASBUFFER_H 
#define OSTAP_MATRIXASBUFFER_H 1
// ============================================================================
// Incldue files
// ============================================================================
// ROOT
// ============================================================================
#include "Math/SMatrix.h"
#include "Math/SVector.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Span.h"
#include "Ostap/Buffer.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Utils 
  {
    // ========================================================================
    /// SMatrix as buffer 
    template <class T, unsigned int D1, unsigned int D2, typename R>
    inline Ostap::Utils::Buffer<T> buffer
    ( const ROOT::Math::SMatrix<T,D1,D2,R>& m )
    { return Ostap::Utils::Buffer<T>( m.Array() , R::kSize ) ; }
    /// SVector as buffer 
    template <class T, unsigned int D>
    inline Ostap::Utils::Buffer<T> buffer
    ( const ROOT::Math::SVector<T,D>& v )
    { return Ostap::Utils::Buffer<T>( v.Array() , D ) ; }
    
    // ========================================================================
#if defined(OSTAP_HAS_STD_SPAN) && OSTAP_HAS_STD_SPAN
    // ========================================================================
    /// SMatrix as std::span
    template <class T, unsigned int D1, unsigned int D2, typename R>
    inline std::span<T> span 
    ( const ROOT::Math::SMatrix<T,D1,D2,R>& m )
    { return std::span<T>( m.Array() , R::kSize ) ; }
    // ========================================================================
    /// SVector as ast::span 
    template <class T, unsigned int D>
    inline std::span<T> span 
    ( const ROOT::Math::SVector<T,D>& v )
    { return std::span<T>( v.Array() , D ) ; }
    // ========================================================================    
#endif // OSTAP_HAS_STD_SPAN
    // ========================================================================
    
    // ========================================================================
  } //                                        The end of namespace Ostap::Utils 
// ============================================================================  
} //                                                 The end of namespace Ostap
// ============================================================================
#endif // OSTAP_MATRIX_AS_BUFFER 
// ============================================================================
//                                                                      The END 
// ============================================================================

