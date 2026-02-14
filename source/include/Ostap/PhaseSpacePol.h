// ============================================================================
#ifndef OSTAP_PHASESPACEPOL_H 
#define OSTAP_PHASESPACEPOL_H 1
// ============================================================================
// Include files
// ============================================================================
// STD &  STL
// ============================================================================
#include <memory>
#include <complex>
#include <vector>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/PhaseSpace.h"
#include "Ostap/Bernstein1D.h"
// ============================================================================
/** @file Ostap/PhaseSpace.h
 *  collection of functions related to the phase space calculations 
 *  @author Vanya Belyaev
 */  
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /** @class PhaseSpacePol
     *  Function to represent the product of l/n-body phase space
     *  and positive polynomial
     *  \f[ \Phi_{l,n}^{(N)(x)} \equiv 
     *      \Phi_{l,n}(x;x_{low},x_{high}) P_{N}(x) \f]
     *  where :
     *  -  \f$  \Phi_{l,n}(x;x_{low},x_{high}) \f$  is a phase space of 
     *     l-particles from n-body decay
     *  -  \f$ P_{N}(x) \f$ is a positive polynomial of degree N
     *  
     *  @see Ostap::Math::PhaseSpaceNL
     *  @see Ostap::Math::Positive
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2011-11-30
     */
    class PhaseSpacePol final : public Ostap::Math::PolyFactor1D
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from thresholds and number of particles
       *  @param threshold_L the low-mass  threshold
       *  @param threshold_H the high-mass threshold
       *  @param l           how many particles we consider
       *  @param n           total number of particles ( n>l!)
       *  @param N           degree of polynomial
       */
      PhaseSpacePol
      ( const double         threshold_L =  0 ,
	const double         threshold_H = 10 ,
	const unsigned short l           =  2 ,
	const unsigned short n           =  3 ,
	const unsigned short N           =  1 ) ; // degree of polynomial
      // =====================================================================
      /** constructor from phase space and polynomial degree
       *  @param ps          phase space factor
       *  @param N           degree of polynomial
       */
      PhaseSpacePol 
      ( const PhaseSpaceNL&  ps      ,
	const unsigned short N  =  1 ) ; // degree of polynomial
      // ======================================================================
      /** constructor from phase space and polynomial degree
       *  @param ps          phase space factor
       *  @param N           degree of polynomial
       */
      PhaseSpacePol 
      ( const PhaseSpaceNL&  ps      ,
	const unsigned short N       ,
	const double         xlow    ,
	const double         xhigh   ) ;
      // ======================================================================
      /// constructor from phase space and polynomial
      PhaseSpacePol 
      ( const PhaseSpaceNL&          ps  ,
	const Ostap::Math::Positive& pol ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate N/L-body modulated phase space
      inline double operator () ( const double x ) const { return evaluate ( x ) ; }
      /// evaluate N/L-body modulated phase space
      double evaluate    ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the phase space 
      const Ostap::Math::PhaseSpaceNL& phasespace () const { return m_phasespace ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral
      double integral () const ;
      /// get the integral between low and high limits
      double integral 
      ( const double low  ,
	const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::PhaseSpaceNL* operator->() const 
      { return &m_phasespace ; }
      // ======================================================================
    private:
      // ======================================================================
      /// the phase space
      Ostap::Math::PhaseSpaceNL   m_phasespace ; // the phase space
      // ======================================================================
    private:
      // ======================================================================
      /// integration workspace
      Ostap::Math::WorkSpace      m_workspace  ;    // integration workspace
      // ======================================================================
    } ;
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_PHASESPACEPOL_H
// ============================================================================
