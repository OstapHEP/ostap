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
    /** @class NSphere 
     *  "N-sphere" of parameters 
     *  @author Vanya Belyaev
     *  @date   2014-01-21
     */
    class NSphere 
    {
    public: 
      // ======================================================================
      /** Standard constructor for rotated sphere 
       *  @param nPhases  dimensionality of N-sphere 
       */
      NSphere ( const unsigned short       nPhases = 1 ) ;
      // =======================================================================
      /** Standard constructor
       *  @param nPhases  dimensionality of N-sphere 
       *  @param rotated  rotate it?
       */
      NSphere ( const unsigned short       nPhases     ,
                const unsigned short       rotated     ) ;
      // =======================================================================
      /** Standard constructor
       *  @param phases  vector of phases 
       *  @param rotated  rotate it?
       */
      NSphere ( const std::vector<double>& phases      ,
                const unsigned short       rotated     ) ;
      // =======================================================================
      /** Standard constructor
       *  @param phases  vector of phases 
       */
      NSphere ( const std::vector<double>& phases      ) ;
      // =======================================================================
      /// copy
      NSphere ( const NSphere&  right ) ;
      /// move
      NSphere (       NSphere&& right ) ;
      /// destructor 
      ~NSphere() ; 
      // ======================================================================
    public:
      // ======================================================================
      unsigned int   nX      () const { return nPhi() + 1       ; } 
      unsigned int   nPhi    () const { return m_sin_phi.size() ; } 
      unsigned short rotated () const { return m_rotated        ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get x_i coefficient:               0 <= i < nX  
      inline double x        ( const unsigned short index ) const ;      
      /// get x_i coefficient squared        0 <= i < nX 
      inline double x2       ( const unsigned short index ) const ;
      inline double xsquared ( const unsigned short index ) const 
      { return x2 ( index ) ; }
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
      double par       ( const unsigned short index ) const 
      { return phase ( index ) ; }
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
      /// number of phases 
      unsigned short             npars   () const { return nPhi()    ; }
      /// get all phases 
      const std::vector<double>& pars    () const { return phases()  ; }
      /** set new value for phi(i)      
       *  @param index (input) the index (0 <= index < nPhi)
       *  @param value new value to be set 
       *  @return true is new value is really set 
       */
      bool setPar     ( const unsigned short index , 
                        const double         value ) 
      { return setPhase ( index , value ) ; }
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
      { return m_rotated && ( index < nPhi() ) ? m_delta [index] : 0.0 ; }
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
    private:
      // ======================================================================
      /// bias to equalize the x_i 
      unsigned short      m_rotated ; // rotated sphere ?
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
    /// swap two speheres 
    inline void swap ( NSphere& a , NSphere& b ) { a.swap ( b ) ; }
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
// get x_i coefficient:               0 <= i < nX  
// ============================================================================
inline double Ostap::Math::NSphere::x 
( const unsigned short indx ) const 
{
  //
  const unsigned int nx = nX() ;
  if      ( nx <= indx       ) { return 0                ; } // invalid 
  else if ( nx  == 1         ) { return 1                ; } // trivial
  else if ( nx  == 1u + indx ) { return m_cos_phi[0]     ; } // trivial
  //
  const unsigned int index = nX() - indx - 1 ;
  //
  // get index as phi 
  //
  double xi = 1.0 ;
  for  ( unsigned short j = 0 ; j < index ; ++j ) { xi *= m_sin_phi[j] ; }
  //
  return index == nPhi() ? xi : xi * m_cos_phi[index] ;
}
// ============================================================================
// get x_i coefficient squared        0 <= i < nX 
// ============================================================================
inline double Ostap::Math::NSphere::x2 
( const unsigned short index ) const  
{
  const double xi = x ( index ) ;
  return xi * xi ;
}
// ============================================================================
// The END 
// ============================================================================
#endif // OSTAP_NSPHERE_H
// ============================================================================
