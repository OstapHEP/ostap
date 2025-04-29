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
    class SineSum    ;
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
       */
      FourierSum 
      ( const unsigned short N    =  0 ,   // degree
        const double         xmin = -1 ,   // low edge
        const double         xmax =  1 ) ; // high edge
      /// copy 
      FourierSum ( const FourierSum&  sum ) = default ;
      ///  move 
      FourierSum (       FourierSum&& sum ) = default ;
      /// constructor from odd numbers of parameters 
      FourierSum
      ( const std::vector<double>& pars  ,
        const double               xmin  ,
        const double               xmax  ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      inline double operator () ( const double x ) const
      { return evaluate ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// calculate Fourier sum
      double evaluate ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get lower edge
      inline double xmin  () const { return m_xmin  ; }
      /// get upper edge
      inline double xmax  () const { return m_xmax  ; }
      /// get the x0 point 
      inline double x0    () const { return m_delta ; }
      /// get the frequence of the first/base harmonic
      inline double omega () const { return m_scale ; }
      // ======================================================================
    public:
      // ======================================================================
      /// t -> x transformartion 
      double x ( const double t ) const { return    t / m_scale   + m_delta ; }
      /// x -> t transformartion 
      double t ( const double x ) const { return  ( x - m_delta ) * m_scale ; }
      // ======================================================================
    public:
      // ======================================================================
      /// maximal trigonometric index 
      unsigned short N () const { return ( m_pars.size() - 1 ) / 2 ; }
      /// get k-th cos-parameter
      double a ( const unsigned short k ) const { return par ( 2 * k ) ; }
      /// get k-th sin-parameter
      double b ( const unsigned short k ) const
      { return 1 <= k ? par   ( 2 * k - 1 ) : 0 ; }
      // set cosine terms
      bool setA ( const unsigned short k , const double value )
      { return setPar ( 2 * k , value ) ; }
      // set   sine terms
      bool setB ( const unsigned short k , const double value )
      { return 1 <= k ? setPar ( 2 * k - 1 , value ) : false ; }
      /** get the magnitude of nth-harmonic
       *  \f$m_k = \sqrt( a^2_k + b^2_k) \f$
       */
      double mag    ( const unsigned short k ) const ;
      /// get the phase for n-th harmonic
      double phase  ( const unsigned short k ) const ;
      // ======================================================================
    public: // integrals & derivatives 
      // ======================================================================
      /// get the derivative at point x
      double     derivative     ( const double x ) const ;
      /// get derivative ad function object
      FourierSum the_derivative () const ;
      // ======================================================================        
    public: // integrals & derivatives 
      // ======================================================================
      // get integral between low and high
      double     integral
      ( const double low ,
        const double high ) const ;
      // ======================================================================
      // get the integral between x0 and x 
      inline double integral ( const double x ) const
      { return integral ( x0 () , x ) ; }
      // ======================================================================
      /**  Get integral as function object           
       *   @attention The linear term p0/2*x is not included! 
       *   It needs to be added explicitely: 
       *   @code
       *   FourierSum fs = ... ; 
       *   FourierSum fi = fs.indefinite_integral();
       *   x  = 0.12 ;
       *   x0 = fs.x0() ; // NOTE:   fi ( x0 ) == 0 
       *   // integral between x0 an X; 
       *   double result = fi ( x ) + 0.5 * cs[0] *  ( x - x0 )  ; 
       *   @endcode
       */
      FourierSum the_integral ( const double C = 0 ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// make convolution with Gaussian function with sigmna 
      FourierSum convolute ( const double sigma ) const ; 
      // ======================================================================
    public:
      // ======================================================================
      /** Calcuate the series using Cesaro methdod rder k 
       *  @param   k  (INPUT)  th esummation order
       *  @return the sum that will be calcualetd usnig Cesaro method        
       */
      FourierSum cesaro ( const unsigned short k ) const ; 
      // ======================================================================
    public:
      // ======================================================================
      /// simple  manipulations with sum : scale it!
      FourierSum& operator *= ( const double a ) ;     // scale it!
      /// simple  manipulations with sum : scale it!
      FourierSum& operator /= ( const double a ) ;     // scale it!
      /// simple  manipulations with sum : add constant
      FourierSum& operator += ( const double a ) ;     // add constant
      /// simple  manipulations with sum : subtract constant
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
      FourierSum __add__     ( const double      value ) const ; 
      FourierSum __mul__     ( const double      value ) const ; 
      FourierSum __sub__     ( const double      value ) const ;  
      FourierSum __truediv__ ( const double      value ) const ;
      FourierSum __div__     ( const double      value ) const { return __truediv__ ( value ) ; }
      FourierSum __radd__    ( const double      value ) const ;
      FourierSum __rmul__    ( const double      value ) const ;
      FourierSum __rsub__    ( const double      value ) const ;
      // ======================================================================
      FourierSum __add__     ( const FourierSum& right ) const ;
      FourierSum __sub__     ( const FourierSum& right ) const ;
      // ======================================================================
      FourierSum __neg__     () const ;
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
      // ======================================================================
    private:
      // ======================================================================
      mutable std::vector<double> m_aux {}  ; // array for derivatives&intgrals
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
       */
      CosineSum 
      ( const unsigned short degree = 0 ,    // degree
        const double         xmin   = 0 ,   // low edge
        const double         xmax   = 1 ) ; // high edge
      /// constructor from Fourier sum
      CosineSum ( const FourierSum& sum ) ;
      /// copy 
      CosineSum ( const CosineSum&  sum ) = default ;
      ///  move 
      CosineSum (       CosineSum&& sum ) = default ;
      /// constructor from non-empty list of parameters
      CosineSum 
      ( const std::vector<double>& pars  ,
        const double               xmin  ,
        const double               xmax  ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      inline double operator () ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// calculate Fourier sum
      double evaluate ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get lower edge
      inline double xmin  () const { return m_xmin  ; }
      /// get upper edge
      inline double xmax  () const { return m_xmax  ; }
      /// get the x0 point 
      inline double x0    () const { return m_xmin  ; }
      /// get the frequence of the first/base harmonic
      inline double omega () const { return m_scale ; }
      // ======================================================================
    public:
      // ======================================================================
      double x ( const double t ) const { return    t / m_scale  + m_xmin  ; }
      double t ( const double x ) const { return  ( x - m_xmin ) * m_scale ; }
      // ======================================================================
    public:
      // ======================================================================
      /// maximal trigonometric index 
      unsigned short N () const { return m_pars.size() - 1 ; }
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
      /// get the derivative at point x
      double  derivative     ( const double x ) const ;
      /// get derivative ad function object
      SineSum the_derivative () const ;
      // ======================================================================
    public:
      // ======================================================================
      // get integral between low and high
      double integral
      ( const double low  ,
        const double high ) const ;
      // ======================================================================
      // get the integral between x0 and x 
      inline double integral ( const double x ) const
      { return integral ( x0 () , x ) ; }
      // ======================================================================
      /**  Get integral as function object           
       *   @attention The term p0/2*x is not included! 
       *   It needs to be added explicitely: 
       *   @code
       *   CosineSum cs = ... ; 
       *   SineSum   ci = cs.the_integral();
       *   x = 0.12     ;
       *   x0 = cs.x0() ; // NOTE  ci ( x0 ) == 0 
       *   double result = ci ( x ) + 0.5 * cs[0] * ( x - x0 ) ; 
       *   @endcode
       */
      SineSum the_integral () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// make convolution with Gaussian function with sigmna 
      CosineSum convolute ( const double sigma ) const ;
      // ======================================================================
    public :
      // ======================================================================
      /** Calcuate the series using Cesaro methdod rder k 
       *  @param   k  (INPUT)  th esummation order
       *  @return the sum that will be calcualetd usnig Cesaro method        
       */
      CosineSum cesaro ( const unsigned short k ) const ; 
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
      CosineSum __add__     ( const double     value ) const ;
      CosineSum __mul__     ( const double     value ) const ;
      CosineSum __sub__     ( const double     value ) const ;
      CosineSum __truediv__ ( const double     value ) const ;
      CosineSum __div__     ( const double     value ) const { return __truediv__ ( value ) ; }
      CosineSum __radd__    ( const double     value ) const ;
      CosineSum __rmul__    ( const double     value ) const ;
      CosineSum __rsub__    ( const double     value ) const ;
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
      // ======================================================================
    private:
      // ======================================================================
      mutable std::vector<double> m_aux {}  ; // array for derivatives&intgrals
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
    /** @class SineSum
     *  Fourier sum over sines
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2025-04-26[
     */
    class  SineSum : public Ostap::Math::Parameters 
    {
    public:
      // ======================================================================
      /** @param degree  degree
       *  @param xmin    low  edge
       *  @param xmax    high edge
       */
      SineSum 
      ( const unsigned short degree = 0 ,    // degree
        const double         xmin   = 0 ,    // low edge
        const double         xmax   = 1 ) ;     // high edge
      /// copy 
      SineSum ( const SineSum&     sum ) = default ;
      ///  move 
      SineSum (       SineSum&&    sum ) = default ;
      /// constructor from on-empty list of parameters
      SineSum 
      ( const std::vector<double>& pars  ,
        const double               xmin  ,
        const double               xmax  ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      inline double operator () ( const double x ) const
      { return evaluate ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// calculate Fourier sum
      double evaluate ( const double x ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get lower edge
      inline double xmin  () const { return m_xmin  ; }
      /// get upper edge
      inline double xmax  () const { return m_xmax  ; }
      /// get the x0 point 
      inline double x0    () const { return m_xmin  ; }
      /// get the frequence of the first/base harmonic
      inline double omega () const { return m_scale ; }
      // ======================================================================
    public:
      // ======================================================================
      double x ( const double t ) const { return    t / m_scale  + m_xmin  ; }
      double t ( const double x ) const { return  ( x - m_xmin ) * m_scale ; }
      // ======================================================================
    public:
      // ======================================================================
      /// maximal trigonometric index 
      unsigned short N () const { return m_pars.size() ; }
      /// get k-th sin-parameter
      double a    ( const unsigned short k ) const { return par ( k     ) ; }
      // set cosine terms
      bool   setA ( const unsigned short k , const double value )
      { return setPar ( k , value ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the derivative at point x
      double    derivative     ( const double x ) const ;
      /// get derivative as function object
      CosineSum the_derivative () const ;
      // ======================================================================
    public:
      // ======================================================================
      // get integral between low and high
      double    integral ( const double low , const double high ) const ;
      // ======================================================================
      // get the integral between x0 and x 
      inline double integral ( const double x ) const
      { return integral ( x0 () , x ) ; }
      // ======================================================================
      /// get integral as function object + ingegrtaion constant 
      CosineSum the_integral ( const double C = 0 ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// make convolution with Gaussian function with sigmna 
      SineSum convolute ( const double sigma ) const ;
      // ======================================================================
    public :
      // ======================================================================
      /** Calcuate the series using Cesaro methdod rder k 
       *  @param   k  (INPUT)  th esummation order
       *  @return the sum that will be calcualetd usnig Cesaro method        
       */
      SineSum cesaro ( const unsigned short k ) const ; 
      // ======================================================================
    public:
      // ======================================================================
      /// simple  manipulations with the sum: scale it!
      SineSum& operator *= ( const double a ) ;     // scale it!
      /// simple  manipulations with the sum : scale it!
      SineSum& operator /= ( const double a ) ;     // scale it!
      // ======================================================================
    public:
      // ======================================================================
      /** sum of two Fourier series (with the same interval!)
       *  @param other the first fourier sum
       *  @return the sum of two Fourier series
       */
      SineSum sum ( const SineSum& other ) const ;
      // ======================================================================
      /// negate it! 
      SineSum operator-() const ;
      // ======================================================================
    public:
      // ======================================================================
      /// swap  it! 
      void swap ( SineSum& right ) ;
      // ======================================================================
    public: // for python
      // ======================================================================
      SineSum __mul__     ( const double   value ) const ;
      SineSum __truediv__ ( const double   value ) const ;
      SineSum __div__     ( const double   value ) const { return __truediv__ ( value ) ; }
      // ======================================================================
      SineSum __rmul__    ( const double   value ) const ;
      // ======================================================================
      SineSum __add__     ( const SineSum& right ) const ;
      SineSum __sub__     ( const SineSum& right ) const ;
      // ======================================================================
      SineSum __neg__     () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// low edge
      double m_xmin  ;             // the low edge
      /// high edge
      double m_xmax  ;             // the high edge
      /// scale factor
      double m_scale ;             // scale factor
      // ======================================================================
    private:
      // ======================================================================
      mutable std::vector<double> m_aux {}  ; // array for derivatives&intgrals
      // ======================================================================
    } ;
    // ========================================================================
    ///  swap  two objects 
    inline void swap ( SineSum& a , SineSum& b ) { a.swap ( b ) ; }
    // ========================================================================
    /// some operators 
    inline SineSum operator+( const SineSum& a , const SineSum& b ) 
    { return a.sum ( b ) ; }
    inline SineSum operator-( const SineSum& a , const SineSum& b ) 
    { return a + (-b); }
    inline SineSum operator*( const SineSum& a , const double   b ) 
    { return SineSum ( a ) *=  b ; }
    inline SineSum operator/( const SineSum& a , const double   b ) 
    { return SineSum ( a ) /=  b ; }
    inline SineSum operator*( const double   b , const SineSum& a )
    { return SineSum ( a ) *=  b ; }
    // ========================================================================


    
    // // ========================================================================
    // namespace Fourier
    // {
    //   // ======================================================================
    //   /** get the derivative of Fourier sum 
    //    *  @param f input the function
    //    *  @param x point x 
    //    *  @reutrn the value of derivative at point x 
    //    */
    //   double derivative
    //   ( const FourierSum& f ,
    //     const double      x ) ;
    //   // ======================================================================
    //   /** get the derivative of Cosine sum 
    //    *  @param f input the function
    //    *  @param x point x 
    //    *  @reutrn the value of derivative at point x 
    //    */
    //   double derivative
    //   ( const CosineSum& f ,
    //     const double     x ) ;
    //   // =======================================================================
    //   /** get the derivative of Sine sum 
    //    *  @param f input the function
    //    *  @param x point x 
    //    *  @reutrn the value of derivative at point x 
    //    */
    //   double derivative
    //   ( const SineSum& f ,
    //     const double   x ) ;
    //   // ======================================================================
    //   /** get the integral for Fouriersum 
    //    *  @param f input the function
    //    *  @param low       lower limitof integration 
    //    *  @param high      lower lumit of integration 
    //    *  @reutrn the integral 
    //    */
    //   double integral 
    //   ( const Ostap::Math::FourierSum& f    ,
    //     const double                   low  , 
    //     const double                   high ) ;
    //   // ======================================================================
    //   /** get the integral for CosineSum 
    //    *  @param f input the function
    //    *  @param low       lower limitof integration 
    //    *  @param high      lower lumit of integration 
    //    *  @reutrn the integral 
    //    */
    //   double integral 
    //   ( const Ostap::Math::CosineSum& f    ,
    //     const double                  low  , 
    //     const double                  high ) ;
    //   // ======================================================================
    //   /** get the integral for SineSum 
    //    *  @param f input the function
    //    *  @param low       lower limitof integration 
    //    *  @param high      lower lumit of integration 
    //    *  @reutrn the integral 
    //    */
    //   double integral 
    //   ( const Ostap::Math::SineSum& f ,
    //     const double                low  , 
    //     const double                high ) ;
    //   // ======================================================================                 
    // }
    
      
    // // ========================================================================
    // /// get the derivative as function
    // FourierSum derivative ( const FourunerSum& ) ;
    // /// get the derivative as function
    // CoSineSum  derivative ( const SineSum&     ) ;
    // /// get the derivative as function
    // SineSum    derivative ( const CosineSum&      ) ;
    // // ========================================================================
    // ///  get a integral  as function 
    // FourierSum integral   ( const FourunerSum& ) ;
    // ///  get a integral  as function 
    // CoSineSum  integral   ( const SineSum&     ) ;
    // ///  get a integral  as function 
    // SineSum    integral   ( const CosineSum&   ) ;
    // // ========================================================================
    // /**  get Gaussian convolution as function
    //  *  @param f (INOUT) the function 
    //  *  @param sigma Gaussian convoluieon 
    //  */
    // FourierSum convolve   ( const FourierSum& , const double sigma ) ;       
    // /**  get Gaussian convolution as function
    //  *  @param f (INOUT) the function 
    //  *  @param sigma Gaussian convoluieon 
    //  */
    // FourierSum convolve   ( const FourierSum& , const double sigma ) ;
    // /**  get Gaussian convolution as function
    //  *  @param f (INOUT) the function 
    //  *  @param sigma Gaussian convoluieon 
    //  */
    // FourierSum convolve   ( const FourierSum& , const double sigma ) ;
    // // ========================================================================
    // /** get Gaussian de-convolution as function
    //  *  @param f (INOUT) the function 
    //  *  @param sigma Gaussian convoluieon 
    //  *  @param delta de-convolution parameter 
    //  */
    // FourierSum deconvolve ( const FourierSum& f , const double sigma , const double delta ) ;       
    // /** get Gaussian de-convolution as function
    //  *  @param f (INOUT) the function 
    //  *  @param sigma Gaussian convoluieon 
    //  *  @param delta de-convolution parameter 
    //  */
    // FourierSum deconvolve ( const FourierSum& f , const double sigma , const double delta ) ;
    // // ========================================================================
    // /** get Gaussian de-convolution as function
    //  *  @param f (INOUT) the function 
    //  *  @param sigma Gaussian convoluieon 
    //  *  @param delta de-convolution parameter 
    //  */
    // FourierSum deconvolve ( const FourierSum& f , const double sigma , const double delta ) ;
    // // ========================================================================


    
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_FOURIER_H
// ============================================================================
