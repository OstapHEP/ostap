#ifndef OSTAP_FOURIER_H 
#define OSTAP_FOURIER_H 1
// ============================================================================
// include files
// ============================================================================
// STD & STL
// ============================================================================
#include <vector>
#include <complex>
// ============================================================================
// Ostap
// ============================================================================
#include "Parameters.h"
// ============================================================================
/** @file Ostap/Fourier.h
 *  set of useful models
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2010-04-19
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    // forward declarations
    class FourierSum ;
    class CosineSum  ;
    // ========================================================================
    /** @class FourierSum
     *  Fourier sum
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-07-26
     */
    class  FourierSum : public Ostap::Math::Parameters 
    {
      // ======================================================================
    public:
      // ======================================================================
      /** @param degree  degree
       *  @param xmin    low  edge
       *  @param xmax    high edge
       *  @param fejer   use fejer summation
       */
      FourierSum 
      ( const unsigned short degree = 0     ,   // degree
        const double         xmin   = 0     ,   // low edge
        const double         xmax   = 1     ,   // high edge
        const bool           fejer  = false );  // use Fejer summation
      /// constructor from cosine serie
      FourierSum ( const CosineSum&  sum ) ;
      /// constructor from Fourier series and fejer flag
      FourierSum ( const FourierSum& sum  , const bool fejer ) ;
      // ======================================================================
      /// copy 
      FourierSum ( const FourierSum&  sum ) ;
      ///  move 
      FourierSum (       FourierSum&& sum ) ;
      /// constructor from odd numbers of parameters 
      FourierSum
      ( const std::vector<double>& pars  ,
        const double               xmin  ,
        const double               xmax  ,
        const double               fejer );
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      inline double operator () ( const double x ) const
      { return m_fejer ? fejer_sum ( x ) : fourier_sum ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// calculate Fourier sum
      double fourier_sum ( const double x ) const ;
      /// calculate Fejer sum
      double fejer_sum   ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get lower edge
      double xmin  () const { return m_xmin  ; }
      /// get upper edge
      double xmax  () const { return m_xmax  ; }
      /// use Fejer summation?
      bool   fejer () const { return m_fejer ; }
      // ======================================================================
    public:
      // ======================================================================
      double x ( const double t ) const { return    t / m_scale   + m_delta ; }
      double t ( const double x ) const { return  ( x - m_delta ) * m_scale ; }
      // ======================================================================
    public:
      // ======================================================================
      /// maximal trigonometric index 
      unsigned short degree () const { return ( m_pars.size() - 1 ) / 2 ; }
      /// get k-th cos-parameter
      double a ( const unsigned short k ) const { return par ( 2 * k     ) ; }
      /// get k-th sin-parameter
      double b ( const unsigned short k ) const
      { return 1 <= k ? par   ( 2 * k - 1 ) : 0 ; }
      // set cosine terms
      bool setA ( const unsigned short k , const double value )
      { return setPar ( 2 * k , value ) ; }
      // set cosine terms
      bool setB ( const unsigned short k , const double value )
      { return 1<= k ? setPar ( 2 * k - 1 , value ) : false ; }
      /** get the magnitude of nth-harmonic
       *  \f$m_k = \sqrt( a^2_k + b^2_k) \f$
       */
      double mag    ( const unsigned short k ) const ;
      /// get the phase for n-th harmonic
      double phase  ( const unsigned short k ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get Fejer sum
      FourierSum fejer_sum   () const ;                       // get Fejer sum
      // ======================================================================
    public:
      // ======================================================================
      /// get the derivative at point x
      double     derivative   ( const double x ) const ;
      /// get the derivative as function
      FourierSum derivative   ( ) const ;
      /// get nth derivative as function
      FourierSum derivative_n ( const unsigned short n ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get integral between low and high
      double     integral   ( const double low , const double high ) const ;
      /** get integral as function
       *  @param c0  integration constant
       */
      FourierSum integral   ( const double c0 = 0 ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /** convolute with gaussian
       *  @param sigma resoltuion parameter for gaussian
       *  @return convolution witgh gaussian
       */
      FourierSum   convolve 
      ( const double sigma     ) const ;
      /** deconvolute with optional regularization
       *  @param sigma sigma of gaussian
       *  @param delta parameter of Tikhonov's regularization
       *  for delta<=0, no regularization
       *  @return regularised deconvolution
       */
      FourierSum deconvolve  
      ( const double sigma     ,
        const double delta = 0 ) const ;
      /**  get the effective cut-off (==number of effective harmonics)
       *   of Tikhonov's regularization
       *   \f$ n \equiv  \sqrt{2 \ln \delta} \frac{\pi\sigma}{L} \f$
       *   @param sigma  gaussian resoltuion
       *   @param delta  regularization parameter
       *   @return number of effective harmonic
       */
      double     regularization
      ( const double sigma     ,
        const double delta     ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// simple  manipulations with polynoms: scale it!
      FourierSum& operator *= ( const double a ) ;     // scale it!
      /// simple  manipulations with polynoms: scale it!
      FourierSum& operator /= ( const double a ) ;     // scale it!
      /// simple  manipulations with polynoms: add constant
      FourierSum& operator += ( const double a ) ;     // add constant
      /// simple  manipulations with polynoms: subtract constant
      FourierSum& operator -= ( const double a ) ;     // subtract constant
      // ======================================================================
    public:
      // ======================================================================
      /** sum of two Fourier series (with the same interval!)
       *  @param other the first fourier sum
       *  @return the sum of two Fourier series
       */
      FourierSum sum ( const FourierSum& other ) const ;
      // ======================================================================
      /** get "shifted" fourier sum
       *  \f$ g(x) \equiv f ( x - a ) \f$
       *  @param a the bias aprameter
       *  @return the shifted fourier sum
       */
      FourierSum shift ( const double a ) const ;
      // ======================================================================
      /// negate it! 
      FourierSum operator-() const ;
      // ======================================================================
    public: // for python
      // ======================================================================
      FourierSum __add__     ( const double value ) const ;
      FourierSum __mul__     ( const double value ) const ;
      FourierSum __sub__     ( const double value ) const ;
      FourierSum __truediv__ ( const double value ) const ;
      FourierSum __div__     ( const double value ) const { return __truediv__ ( value ) ; }
      FourierSum __radd__    ( const double value ) const ;
      FourierSum __rmul__    ( const double value ) const ;
      FourierSum __rsub__    ( const double value ) const ;
      // ======================================================================
      FourierSum __add__   ( const FourierSum& right ) const ;
      FourierSum __sub__   ( const FourierSum& right ) const ;
      // ======================================================================
      FourierSum __neg__   () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// swap  it! 
      void swap ( FourierSum& right ) ;
      // ======================================================================
    private:
      // ======================================================================
      /// low edge
      double m_xmin  ;             // the low edge
      /// high edge
      double m_xmax  ;             // the high edge
      /// scale factor
      double m_scale ;             // scale factor
      /// delta
      double m_delta ;             // delta
      /// summation algorithm
      bool   m_fejer ;             // summation algorithm
       // ======================================================================
    } ;
    // ========================================================================
    ///  swap  two objects 
    inline void swap ( FourierSum& a , FourierSum& b ) { a.swap ( b ) ; }
    /// some operators 
    inline FourierSum operator+( const FourierSum& a , const FourierSum& b ) 
    { return a.sum ( b ) ; }
    inline FourierSum operator-( const FourierSum& a , const FourierSum& b ) 
    { return a + (-b); }
    inline FourierSum operator+( const FourierSum& a , const double      b ) 
    { return FourierSum (  a )+= b ; }
    inline FourierSum operator*( const FourierSum& a , const double      b ) 
    { return FourierSum (  a )*= b ; }
    inline FourierSum operator/( const FourierSum& a , const double      b ) 
    { return FourierSum (  a )/= b ; }
    inline FourierSum operator+( const double      b , const FourierSum& a )
    { return FourierSum (  a )+= b ; }
    inline FourierSum operator-( const double      b , const FourierSum& a )
    { return  -a + b ; }
    inline FourierSum operator*( const double      b , const FourierSum& a )
    { return FourierSum ( a )*= b ; }
    // ========================================================================
    /** @class CosineSum
     *  Fourier sum over cosines
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-07-26
     */
    class  CosineSum : public Ostap::Math::Parameters 
    {
    public:
      // ======================================================================
      /** @param degree  degree
       *  @param xmin    low  edge
       *  @param xmax    high edge
       *  @param fejer   use fejer summation
       */
      CosineSum 
      ( const unsigned short degree = 0     ,    // degree
        const double         xmin   = 0     ,    // low edge
        const double         xmax   = 1     ,    // high edge
        const bool           fejer  = false ) ;  // use Fejer summation
      /// constructor from Fourier sum
      CosineSum ( const FourierSum&    sum            ) ;
      /// constructor from Fourier series and fejer flag
      CosineSum ( const CosineSum&     sum  , const bool fejer ) ;
      // ======================================================================
      /// copy 
      CosineSum ( const CosineSum&  sum ) ;
      ///  move 
      CosineSum (       CosineSum&& sum ) ;
      /// constructor from on-empty list of parameters
      CosineSum 
      ( const std::vector<double>& pars  ,
        const double               xmin  ,
        const double               xmax  ,
        const double               fejer );
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      double operator () ( const double x ) const
      { return m_fejer ? fejer_sum ( x ) : fourier_sum ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// calculate Fourier sum
      double fourier_sum ( const double x ) const ;
      /// calculate Fejer sum
      double fejer_sum   ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get lower edge
      double xmin  () const { return m_xmin  ; }
      /// get upper edge
      double xmax  () const { return m_xmax  ; }
      /// use Fejer summation?
      bool   fejer () const { return m_fejer ; }
      // ======================================================================
    public:
      // ======================================================================
      double x ( const double t ) const { return    t / m_scale  + m_xmin  ; }
      double t ( const double x ) const { return  ( x - m_xmin ) * m_scale ; }
      // ======================================================================
    public:
      // ======================================================================
      /// degree  of polynomial
      unsigned short degree () const { return m_pars.size() - 1 ; }
      /// get k-th cos-parameter
      double a    ( const unsigned short k ) const { return par ( k     ) ; }
      // set cosine terms
      bool   setA ( const unsigned short k , const double value )
      { return setPar ( k , value ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get Fejer sum
      CosineSum fejer_sum   () const ;                         // get Fejer sum
      // ======================================================================
    public:
      // ======================================================================
      /// get the derivative at point x
      double     derivative   ( const double x ) const ;
      /// get the derivative as function
      FourierSum derivative   ( ) const ;
      /// get nth derivative as function
      FourierSum derivative_n ( const unsigned short n ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get integral between low and high
      double     integral   ( const double low , const double high ) const ;
      /** get integral as function
       *  @param c0  integration constant
       */
      FourierSum integral   ( const double c0 = 0 ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /** convolute with gaussian
       *  @param sigma resoltuion parameter for gaussian
       *  @return convolution witgh gaussian
       */
      CosineSum   convolve
      ( const double sigma     ) const ;
      /** deconvolute with optional regularization
       *  @param sigma sigma of gaussian
       *  @param delta parameter of Tikhonov's regularization
       *  for delta<=0, no regularization
       *  @return regularised deconvolution
       */
      CosineSum deconvolve
      ( const double sigma     ,
        const double delta = 0 ) const ;
      /** get the effective cut-off (==number of terms/harmonics)
       *  of Tikhonov's regularization
       *  \f$ n \equiv  \sqrt{2 \ln \delta} \frac{\pi\sigma}{L} \f$
       *  @param sigma  gaussian resoltuion
       *  @param delta  regularization parameter
       *  @return number of effective harmonic
       */
      double    regularization
      ( const double sigma     ,
        const double delta     ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// simple  manipulations with polynoms: scale it!
      CosineSum& operator *= ( const double a ) ;     // scale it!
      /// simple  manipulations with polynoms: scale it!
      CosineSum& operator /= ( const double a ) ;     // scale it!
      /// simple  manipulations with polynoms: add constant
      CosineSum& operator += ( const double a ) ;     // add constant
      /// simple  manipulations with polynoms: subtract constant
      CosineSum& operator -= ( const double a ) ;     // subtract constant
      // ======================================================================
    public:
      // ======================================================================
      /** sum of two Fourier series (with the same interval!)
       *  @param other the first fourier sum
       *  @return the sum of two Fourier series
       */
      CosineSum sum ( const CosineSum& other ) const ;
      // ======================================================================
      /// negate it! 
      CosineSum operator-() const ;
      // ======================================================================
    public:
      // ======================================================================
      /// swap  it! 
      void swap ( CosineSum& right ) ;
      // ======================================================================
    public: // for python
      // ======================================================================
      CosineSum __add__     ( const double value ) const ;
      CosineSum __mul__     ( const double value ) const ;
      CosineSum __sub__     ( const double value ) const ;
      CosineSum __truediv__ ( const double value ) const ;
      CosineSum __div__     ( const double value ) const { return __truediv__ ( value ) ; }
      // ======================================================================
      CosineSum __radd__    ( const double value ) const ;
      CosineSum __rmul__    ( const double value ) const ;
      CosineSum __rsub__    ( const double value ) const ;
      // ======================================================================
      CosineSum __add__     ( const CosineSum& right ) const ;
      CosineSum __sub__     ( const CosineSum& right ) const ;
      // ======================================================================
      CosineSum __neg__     () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// low edge
      double m_xmin  ;             // the low edge
      /// high edge
      double m_xmax  ;             // the high edge
      /// scale factor
      double m_scale ;             // scale factor
      /// summation algorithm
      bool m_fejer   ;             // summation algorithm
      // ======================================================================
    } ;
    // ========================================================================
    ///  swap  two objects 
    inline void swap ( CosineSum& a , CosineSum& b ) { a.swap ( b ) ; }
    // ========================================================================
    /// some operators 
    inline CosineSum operator+( const CosineSum& a , const CosineSum& b ) 
    { return a.sum ( b ) ; }
    inline CosineSum operator-( const CosineSum& a , const CosineSum& b ) 
    { return a + (-b); }
    inline CosineSum operator+( const CosineSum& a , const double     b ) 
    { return CosineSum ( a )+=  b ; }
    inline CosineSum operator*( const CosineSum& a , const double     b ) 
    { return CosineSum ( a )*=  b ; }
    inline CosineSum operator/( const CosineSum& a , const double     b ) 
    { return CosineSum ( a )/=  b ; }
    inline CosineSum operator+( const double     b , const CosineSum& a )
    { return CosineSum ( a )+=  b ; }
    inline CosineSum operator-( const double     b , const CosineSum& a )
    { return -a + b ; }
    inline CosineSum operator*( const double     b , const CosineSum& a )
    { return CosineSum ( a )*=  b ; }
    // ========================================================================    
    /** make a sum of two fourier series (with the same interval!)
     *  @param s1 the first fourier sum
     *  @param s2 the first fourier sum
     *  @return s1+s2
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2016-06-26
     */
    FourierSum sum ( const FourierSum& s1 , const FourierSum& s2 ) ;
    // ========================================================================
    /** make a sum of two fourier cosine series (with the same interval!)
     *  @param s1 the first fourier cosine sum
     *  @param s2 the first fourier cosine sum
     *  @return s1+s2
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2016-06-26
     */
    CosineSum sum ( const CosineSum& s1 , const CosineSum& s2 ) ;
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_FOURIER_H
// ============================================================================
