// ============================================================================
#ifndef OSTAP_CHI2FIT_H 
#define OSTAP_CHI2FIT_H 1
// ============================================================================
// Include files
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/ValueWithError.h"
// ============================================================================
/** @file Ostap/Chi2Fit.h
 *   Trivial chi2-fit 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2012-05-26
 */  
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /** @class Chi2Fit LHCbMath/Chi2Fit.h
    *   Trivial chi2-fit 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2012-05-26
     */
    class Chi2Fit 
    {
    public: 
      // =======================================================================
      /// the atomic type 
      typedef Ostap::Math::ValueWithError  VE   ;
      /// the data-vector ("histogram")
      typedef std::vector<VE>              DATA ;
      /// the components 
      typedef std::vector<DATA>            CMPS ;
      // ======================================================================
    public: 
      // ======================================================================
      /// fit with one component 
      Chi2Fit ( const DATA& data ,    // the data 
                const DATA& cmp  ) ;  // the component
      /// fit with many component 
      Chi2Fit ( const DATA& data ,    // the data 
                const CMPS& cmps ) ;  // the components
      /// destructor 
      ~Chi2Fit () ; // destructor 
      // ======================================================================
    public:
      // ======================================================================
      std::size_t size () const { return m_cmps.size() ; }
      // ======================================================================
      /// data used in the fit 
      const DATA& data () const { return m_data ; }
      /// components used in the fit 
      const CMPS& cmps () const { return m_cmps ; }      
      /// init values and steps 
      const DATA& init () const { return m_init ; }      
      // ======================================================================
    public : // Fit results 
      // ======================================================================
      StatusCode  status () const { return m_code ; }
      /// get the parametr 
      VE          param  ( const unsigned int index ) const ;
      /// get the covariance matrix elements
      double      cov2   ( const unsigned int i1 , 
                           const unsigned int i2 ) const ;
      /// the function at minimum
      double      chi2   () const { return m_chi2   ; } // the function at minimum
      /// number of function calls 
      std::size_t ncalls () const { return m_calls  ; } // function calls 
      /// number of iterations 
      std::size_t niters () const { return m_iters  ; }
      /// number of points 
      std::size_t points () const { return m_points ; }
      // ======================================================================
      /// get all parameters at once 
      DATA   params () const ; // get all parameters at once 
      // ======================================================================
    public:
      // ======================================================================
      /// output 
      std::ostream& fillStream ( std::ostream& s ) const ;
      /// conversion to string 
      std::string   toString   ( ) const ;
      // ======================================================================
    private:
      // ======================================================================
      DATA m_data ;
      CMPS m_cmps ;
      // the init values and steps 
      DATA m_init ;
      // ======================================================================
    private: // fit results 
      // ======================================================================
      StatusCode  m_code   ;
      void*       m_solu   ;
      void*       m_cov2   ;
      // ======================================================================
      double      m_chi2   ;
      std::size_t m_calls  ;
      std::size_t m_iters  ;
      std::size_t m_points ;
      // ======================================================================
    };
    // ========================================================================
    /// output operator 
    inline std::ostream& operator<< ( std::ostream& s , const Chi2Fit& f ) 
    { return f.fillStream ( s ) ; }
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_CHI2FIT_H
// ============================================================================
