// ============================================================================
#ifndef OSTAP_CHEBYSHEVAPPROXIMATION_H 
#define OSTAP_CHEBYSHEVAPPROXIMATION_H 1
// ============================================================================
//  include files
// ============================================================================
// STD&STL
// ============================================================================
#include <functional>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/ValueWithError.h"
// ============================================================================
// forward declarations 
// ============================================================================
namespace Ostap { namespace Functions { class  PyCallable ; } } // from Ostap 
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    /** @class ChebyshevApproximation Ostap/ChebyshevApproximation.h
     *  Helper class  for Chebyshev Approximation. 
     *  Unlike the similar ROOT's class ROOT::Math::ChebyshedApprox
     *  it could be copied and moved.
     *  @see https://en.wikipedia.org/wiki/Approximation_theory#Chebyshev_approximation
     *  @see https://www.gnu.org/software/gsl/doc/html/cheb.html
     *  @see ROOT::Math::ChebyshedApprox
     *  @author Vanya Belyaev
     *  @date   2019-09-25
     */
    class ChebyshevApproximation 
    {
    public:
      // ======================================================================
      /** constructor from the function, low/high-limits and the approximation order 
       *  @param func the function 
       *  @param a the low-limit 
       *  @param b the high-limit 
       *  @param N the approximation order 
       */
      ChebyshevApproximation 
      ( const std::function<double(double)>& func , 
        const double                         a    , 
        const double                         b    , 
        const unsigned short                 N    ) ;
      // ======================================================================
      /** constructor from the function, low/high-limits and the approximation order 
       *  @param func the function 
       *  @param a the low-limit 
       *  @param b the high-limit 
       *  @param N the approximation order 
       */
      ChebyshevApproximation 
      ( const Ostap::Functions::PyCallable& func , 
        const double             a    , 
        const double             b    , 
        const unsigned short     N    ) ;
      // ======================================================================
      /// copy constructor 
      ChebyshevApproximation ( const ChebyshevApproximation&  right ) ;
      /// move constructor 
      ChebyshevApproximation (       ChebyshevApproximation&& right ) ;
      // ======================================================================
      /// destructor: deallocate GSL structures 
      ~ChebyshevApproximation() ;
      // ======================================================================
    protected:
      // ======================================================================
      // default (protected) constructor 
      ChebyshevApproximation () ;
      // ======================================================================
    public:
      // ======================================================================
      ///  the main method: evaluate the approximation sum 
      double operator () ( const double         x ) const { return evaluate ( x ) ; }
      /** the main method: evaluate the approximation sum 
       *  using at most <code>n</code> terms 
       */
      double operator () ( const double         x , 
                           const unsigned short n ) const { return evaluate ( x , n ) ; }
      // ======================================================================
    public:
      // ======================================================================
      ///  the main method: evaluate the approximation sum 
      double evaluate ( const double         x ) const ;
      /** the main method: evaluate the approximation sum, 
       *  using at most <code>n</code> terms 
       */
      double evaluate ( const double         x , 
                        const unsigned short n ) const ;
      // ======================================================================
    public: // trivial accessors 
      // ======================================================================
      /// get the low edge 
      double         a () const { return m_a ; }
      /// get the high edge 
      double         b () const { return m_b ; }
      /// get the approximation order 
      unsigned short N () const { return m_N ; }
      // ======================================================================
    public:
      // ======================================================================
      double xmin () const { return m_a ; }
      double xmax () const { return m_b ; }
      // ======================================================================
    public: // derivatives and integrals 
      // ======================================================================
      /** the main method: evaluate the approximation sum 
       *  @return the approximation with the error estimate 
       */
      Ostap::Math::ValueWithError
      eval_err ( const double         x ) const ;
      /** the main method: evaluate the approximation sum
       *  using at most <code>n</code> terms 
       *  @return the approximation with the error estimate 
       */
      Ostap::Math::ValueWithError 
      eval_err ( const double         x , 
                 const unsigned short n ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get a derivative  
      ChebyshevApproximation derivative () const ;
      /// get an integral: \f$ F(x) \equiv \int_a^{z} f(t) \deriv t  + C \f$ 
      ChebyshevApproximation integral   ( const double C = 0 ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// swap two objects 
      void swap ( ChebyshevApproximation&  right ) ;
      // ======================================================================
    private :
      // ======================================================================
      /// low edge 
      double m_a ;                                       // low edge 
      /// high edge 
      double m_b ;                                       // high edge 
      /// approximation order
      unsigned short m_N ;                               // approximation order
      // ======================================================================
    private:
      // ======================================================================
      /** the actual GSL structure <code>gsl_cheb_series</code>
       *  <code>char*</code> here instead of <code>void*</code>
       *  to please the dictionary generator 
       *  @see gsl_cheb_series  
       */
      char* m_chebyshev { nullptr } ; // the actual GSL structure
      // ======================================================================
    };
    // ========================================================================
    /// swap two objects 
    inline void swap ( ChebyshevApproximation& a , 
                       ChebyshevApproximation& b ) { a.swap ( b ) ;}
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_CHEBYSHEVAPPROXIMATION_H
// ============================================================================
