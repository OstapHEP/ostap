// ============================================================================
#ifndef OSTAP_NSPHERE_H 
#define OSTAP_NSPHERE_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <vector>
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ==========================================================================
    /** @class NSphere Ostap/NSphere.h
     *  "N-sphere" of parameters: useful class to get 
     *   certain combinations of normalized/constrained variables     
     *   parameterized with several angular-like phase parameters       
     *
     *  Useful cases:
     *  - Get n-parameters \f$ \x_i f$, such \f$ \sum x_i^2 = 1 \f$
     *  - Get n-parameters \f$ \x_i f$, such \f$0 \le x_i < 1\f$ and  \f$ \sum x_i = 1 \f$
     *  - Get n-parameters \f$ \x_i f$, such \f$ 0 \le ... \le x_i \le x_{i+1} \le ..  \le 1 \f$ 
     *
     *  @author Vanya BELYAEV   Ivan.Belyaev@itep.ru 
     *  @date   2014-01-21
     */
    class NSphere 
    {
    public: 
      // ======================================================================
      /** Standard constructor for rotated sphere 
       *  @param nPhases  dimensionality of N-sphere 
       */
      NSphere  ( const unsigned short       nPhases = 1 ) ;
      // =======================================================================
      /** Standard constructor
       *  @param nPhases  dimensionality of N-sphere 
       *  @param rotated  rotate it?
       */
      NSphere  ( const unsigned short       nPhases     ,
                 const unsigned short       rotated     ) ;
      // =======================================================================
      /** Standard constructor
       *  @param phases  vector of phases 
       *  @param rotated  rotate it?
       */
      NSphere  ( const std::vector<double>& phases      ,
                 const unsigned short       rotated     ) ;
      // =======================================================================
      /** Standard constructor with deltas 
       *  @param phases  vector of phases 
       *  @param deltas  rotation deltas 
       */
      NSphere  ( const std::vector<double>& phases ,
                 const std::vector<double>& deltas ) ;
      // =======================================================================
      /** constructor from deltas 
       *  @param deltas  rotation deltas 
       */
      NSphere  ( const std::string&       /* fake */ ,
                 const std::vector<double>& deltas ) ;
      // =======================================================================
      /** Standard constructor
       *  @param phases  vector of phases 
       */
      NSphere  ( const std::vector<double>& phases      ) ;
      // =======================================================================
      /// copy
      NSphere  ( const NSphere&  right ) ;
      /// move
      NSphere  (       NSphere&& right ) ;
      /// destructor 
      ~NSphere () ; 
      // ======================================================================
    public:
      // ======================================================================
      unsigned int   nX      () const { return nPhi () + 1      ; } 
      unsigned int   nPhi    () const { return m_sin_phi.size() ; } 
      // ======================================================================
    public:
      // ======================================================================
      /// get x_i coefficient:               0 <= i < nX  
      inline double x        ( const unsigned short index ) const ;
      // get x_i coefficient squared        0 <= i < nX 
      inline double x2       ( const unsigned short index ) const 
      { const double xi = x ( index ) ; return xi * xi ; }
      // get x_i coefficient squared        0 <= i < nX 
      inline double xsquared ( const unsigned short index ) const 
      { const double xi = x ( index ) ; return xi * xi ; }
      // ======================================================================
    public:
      // ======================================================================
      double sin_phi   ( const unsigned short index ) const 
      { return index < nPhi() ? m_sin_phi [index] : 0.0 ; }
      // ======================================================================
      double cos_phi   ( const unsigned short index ) const 
      { return index < nPhi() ? m_cos_phi [index] : 0.0 ; }
      // ======================================================================
    public:
      // ======================================================================
      double phase     ( const unsigned short index ) const 
      { return index < nPhi() ? m_phases  [index] : 0.0 ; }
      // ======================================================================
      /// get all phases 
      const std::vector<double>& phases  () const { return m_phases  ; }
      // get all   sines 
      const std::vector<double>& sines   () const { return m_sin_phi ; }
      // get all cosines 
      const std::vector<double>& cosines () const { return m_cos_phi ; }
      // get all deltas 
      const std::vector<double>& delta   () const { return m_delta   ; }
      // ======================================================================
    public:
      // ======================================================================
      /** set new value for phi(i)      
       *  @param index (input) the index (0 <= index < nPhi)
       *  @param value new value to be set 
       *  @return true is new value is really set 
       */
      bool setPhase   ( const unsigned short index , 
                        const double         value ) ;
      // ======================================================================
    public: // "par"-like interface
      // ======================================================================
      /// number of parameters/phases 
      unsigned short             npars   () const { return nPhi()    ; }
      /// get all parameters/phases 
      const std::vector<double>& pars    () const { return phases()  ; }
      /// get the parameter/phase 
      double par       ( const unsigned short index ) const 
      { return phase ( index ) ; }
      // ======================================================================
      /** set new value for phi(i)      
       *  @param index (input) the index (0 <= index < nPhi)
       *  @param value new value to be set 
       *  @return true is new value is really set 
       */
      bool setPar     ( const unsigned short index , 
                        const double         value ) 
      { return setPhase ( index , value ) ; }
      // ======================================================================
      /** set several/all parameters at once 
       *  @param begin  start itertaor for sequence of coefficients 
       *  @param end    end  titertaor for sequence of coefficients 
       *  @return true if at least one parameter is actually changed 
       */
      template <class ITERATOR>
      bool setPars ( ITERATOR begin  , 
                     ITERATOR end    ) ;
      // ======================================================================
      /** set several/all parameters at once 
       *  @param pars (NIPUT) vector of parameters 
       *  @return true if at least one parameter is actually changed 
       */
      inline bool setPars ( const std::vector<double>& pars ) 
      { return setPars ( pars.begin() , pars.end() ) ; }
      // ======================================================================
    public:
      // ======================================================================
      ///  sphere "size" 
      unsigned int size () const { return nPhi () ; } // sphere size 
      ///  sphere dimension  
      unsigned int dim  () const { return nPhi () ; } // sphere dimension 
      // ======================================================================      
    public:
      // ======================================================================
      /// get biases to equalize the data 
      double delta    ( const unsigned short index ) const  
      { return index < nPhi() ? m_delta [ index ] : 0.0 ; }
      // ======================================================================
      /// get the phase biases 
      const std::vector<double>& deltas () const { return m_delta ; }
      // ======================================================================      
    public:
      // ======================================================================
      /// copy assignement 
      NSphere& operator=( const NSphere&  right ) ;
      /// move assignement 
      NSphere& operator=(       NSphere&& right ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// swap two spheres 
      void swap ( NSphere& right ) ; // swap two spheres 
      // ======================================================================
    public:
      // ======================================================================
      /** convert n-coordinates \f$ x_i \f$ into (n-1) phases \f$ \phi_i\f$  
       *  @param  x vector  in n-dimensional space 
       *  @return spherical phases
       */
      static std::vector<double> phis ( const std::vector<double>& x ) ;
      // ======================================================================
    private:
      // ======================================================================
      /// the phase biases for rotated sphere 
      std::vector<double> m_delta   ; // the phase biases for rotated sphere 
      /// the phases  
      std::vector<double> m_phases  ; // the phases 
      /// vector of sin(phi)
      std::vector<double> m_sin_phi ; // vector of sin(phi)
      /// vector of cos(phi)
      std::vector<double> m_cos_phi ; // vector of cos(phi)
      // ======================================================================
    };
    // ========================================================================
    /// swap two spheres 
    inline void swap ( NSphere& a , NSphere& b ) { a.swap ( b ) ; }
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
// get x_i coefficient:  0 <= i < nX  
// ============================================================================
inline double Ostap::Math::NSphere::x 
( const unsigned short index ) const 
{
  //
  const unsigned int nx = nX () ;
  if      ( nx <= index ) { return 0               ; } // invalid 
  else if ( nx == 1     ) { return 1               ; } // trivial
  else if ( 0  == index ) { return m_cos_phi [ 0 ] ; } // trivial
  //
  const bool last = ( index + 1 == nx ) ;
  //
  long double xi = 1.0 ;
  for ( unsigned short j = 0 ; j < index ; ++j ) { xi *= m_sin_phi[j] ; }
  //
  return last ? xi : xi * m_cos_phi [ index ] ;
}
// ============================================================================
/** set several/all parameters at once 
 *  @param pars (NIPUT) vector of parameters 
 *  @return true if at least one parameter is actually changed 
 */
template <class ITERATOR>
inline bool Ostap::Math::NSphere::setPars ( ITERATOR begin  , 
                                            ITERATOR end    ) 
{
  bool update = false ;
  const unsigned int   N = nPhi ()  ;
  for ( unsigned short k ; k < N && begin != end ;  ++k, ++begin ) 
  { update = setPar ( k  , *begin ) | update ; }
  return update ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_NSPHERE_H
// ============================================================================
