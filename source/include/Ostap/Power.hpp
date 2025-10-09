// ============================================================================
#ifndef OSTAP_POWER_HPP 
#define OSTAP_POWER_HPP 1
// ============================================================================
// Include files
// ============================================================================
/** @file Ostap/Power.hpp
 *  @author Juan PALACIOS juan.palacios@cern.ch
 *  @date 2006-10-26 
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // =======================================================================
    /** @struct Power 
     *  Template metafunction to calculate integer powers of integer
     *  and floating point numbers.
     *
     *  @code 
     *   // calculate 50.**10
     *   const double result = Ostap::Math::pow<double, 10> ( 50. ) ;
     *
     *  @code 
     *
     *  @author Juan PALACIOS juan.palacios@cern.ch
     *  @date 2006-10-26 
     */
    template<typename TYPE, int N>
    struct Power 
    {
      static inline TYPE pow ( TYPE __x ) 
      { 
        return Power<TYPE, N-1>::pow(__x)*__x;
      }
    };
    // =======================================================================
    template<typename TYPE>
    struct Power<TYPE, 1> 
    {      
      static inline TYPE pow ( TYPE __x ) { return __x; }
    };
    // =======================================================================    
    template<typename TYPE>
    struct Power<TYPE, 0> 
    {
      static inline TYPE pow ( TYPE /* __x */ ) { return 1; }
    };
    // =======================================================================    
    template<typename TYPE, int N, bool C>
    struct ReturnPolicy 
    {
      typedef TYPE RETURN_TYPE;
    };
    // =======================================================================        
    template<typename TYPE, int N>
    struct ReturnPolicy<TYPE,N,false> 
    {
      typedef double RETURN_TYPE;
    };
    // =======================================================================        
    template<typename TYPE, int N>
    struct InvPower
    {
      static inline typename ReturnPolicy<TYPE, N, (N>=0) >::RETURN_TYPE pow( TYPE __x )
      {
        return 1./Power<TYPE, -1*N>::pow(__x);
      }  
    };
    // =======================================================================        
    template <typename TYPE, int N, bool C>
    struct ImplementationSwitch 
    {
      typedef Ostap::Math::Power<TYPE, N>  powImpl;
    };
    // =======================================================================        
    template <typename TYPE, int N>
    struct ImplementationSwitch<TYPE, N, false> 
    {
      typedef Ostap::Math::InvPower<TYPE, N>  powImpl;
    };
    // =======================================================================        
    template<typename TYPE, int N>
    static inline typename ReturnPolicy<TYPE, N, (N>=0) >::RETURN_TYPE pow ( TYPE __x ) 
    {
      return ImplementationSwitch< TYPE, N, (N>=0) >::powImpl::pow(__x);
    };
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_POWER_HPP
// ============================================================================
