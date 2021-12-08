// ============================================================================
#ifndef OSTAP_INTERPOLANTS_H 
#define OSTAP_INTERPOLANTS_H 1
// ============================================================================
// Include files 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/Interpolation.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math  
  {
    // ========================================================================
    /** @class Neville 
     *  Simple interpolation polynomial using Neville's algorithm 
     *  @see https://en.wikipedia.org/wiki/Neville%27s_algorithm
     *  @attention  it is not CPU efficient!
     */
    class Neville : public Interpolation::Table 
    {
    public:
      // ====================================================================== 
      /** simple contsructor from the interpolation points  
       *  @param p input vector of abscissas  
       */
      Neville 
      ( const Interpolation::Table& p  ) ;
      // ====================================================================
      /** simple contsructor from the interpolation points  
       *  @param p input vector of abscissas  
       */
      Neville 
      (       Interpolation::Table&& p ) ;
      // ====================================================================
      // other constructors from the base class 
      // ====================================================================
      using Interpolation::Table::Table ;
      /// default constructor
      Neville () = default ;
      // =====================================================================
    public:
      // ======================================================================
      /// the main method: get the value of interpolated polynomial 
      double evaluate    ( const  double x ) const { return neville  ( x ) ; }
      /// the main method: get the value of interpolated polynomial 
      double operator () ( const  double x ) const { return evaluate ( x ) ; }
      // ======================================================================
      /// get the derivative   (dy/dx) at point x 
      double derivative  ( const  double x ) const { return neville2( x ).second  ; }
      // ======================================================================
    public : 
      // ======================================================================
      /// swap two interpolators 
      void exchange ( Neville& right ) 
      { Interpolation::Table::exchange ( right ) ; }
      // ======================================================================      
    } ;    
    // ========================================================================
    /// swap two interpolators 
    inline void swap ( Neville& a , Neville& b ) { a.exchange ( b ) ; }
    // ========================================================================
    /** @class Lagrange  
     *  https://en.wikipedia.org/wiki/Lagrange_polynomial
     *  @attention  it is not CPU efficient and numerically not stable 
     */
    class Lagrange : public Interpolation::Table 
    {
    public:
      // ====================================================================== 
      /** simple constructor from the interpolation points  
       *  @param p input vector of abscissas  
       */
      Lagrange 
      ( const  Interpolation::Table& p ) ;
      // ====================================================================
      /** simple constructor from the interpolation points  
       *  @param p input vector of abscissas  
       */
      Lagrange 
      (        Interpolation::Table&& p ) ;
      // ====================================================================
      // other constructors from the base class 
      // ====================================================================
      using Interpolation::Table::Table ;
      // ======================================================================
      /// default constructor
      Lagrange () = default ;
      // =====================================================================
    public:
      // ======================================================================
      /// the main method: get the value of interpolated polynomial 
      double evaluate    ( const  double x ) const { return lagrange ( x ) ; }
      /// the main method: get the value of interpolated polynomial 
      double operator () ( const  double x ) const { return evaluate ( x ) ; }
      // ======================================================================
      /// get the drivative with respect to i-th parameter (dy/d(y_i)) at point x 
      double derivative  ( const double       x , 
                           const unsigned int iy ) const 
      { return lagrange2 ( x , iy ).second  ; }
      // ======================================================================
    public : 
      // ======================================================================
      /// swap two interpolators 
      void exchange ( Lagrange& right ) 
      { Interpolation::Table::exchange ( right ) ; }
      // ======================================================================      
    } ;    
    // ========================================================================
    /// swap two Lagrange interpolators 
    inline void swap ( Lagrange& a , Lagrange& b ) { a.exchange ( b ) ; }
    // ========================================================================
    /** @class Berrut1st
     *  Very efficient rational 1st Berrut's interpolant
     *  \f[ F_n(x) = \frac{ \sum_i \frac{\beta_i}{x-x_i} f_i}
     *                    { \sum_i \frac{\beta_i}{x-x_i} } \f]
     *  where \f$ \beta_i = (-1)^i\f$.
     *  It is barycentric-like rational interpolation.
     *  For odd number of points it is barycentric, 
     * - performance is about O(n)
     *  @date 2021-11-25
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     */
    class Berrut1st : public Interpolation::Table 
    {
    public:
      // ======================================================================
      /** simple constructor from the interpolation points  
       *  @param p input vector of abscissas  
       */
      Berrut1st 
      ( const  Interpolation::Table& p ) ;
      // ====================================================================
      /** simple constructor from the interpolation points  
       *  @param p input vector of abscissas  
       */
      Berrut1st 
      (        Interpolation::Table&& p ) ;
      // ====================================================================
      // other constructors from the base class 
      // ====================================================================
      using Interpolation::Table::Table ;
      // ======================================================================
    public:
      // ======================================================================
      /// the main method: get the value of interpolant 
      double evaluate    ( const  double x ) const { return berrut1st ( x ) ; }
      /// the main method: get the value of interpolant 
      double operator () ( const  double x ) const { return evaluate  ( x ) ; }
      // ======================================================================
      /// get the weight 
      double weight ( const unsigned short index ) const
      { return size () <= index ? 0 : ( 0 == index % 2 ) ? 1 : -1 ; }
      // ======================================================================
    public : 
      // ======================================================================
      /// swap two interpolators 
      void exchange ( Berrut1st& right ) 
      { Interpolation::Table::exchange ( right ) ; }
      // ======================================================================      
    } ;    
    // ========================================================================
    /// swap two Berrut1st interpolators 
    inline void swap ( Berrut1st& a , Berrut1st& b ) { a.exchange ( b ) ; }
    // ========================================================================
    /** @class Berrut2nd
     *  Very efficient 2nd Berrut's rational interpolant
     *  \f[ F_n(x) = \frac{ \sum_i \frac{\beta_i}{x-x_i} f_i}
     *                    { \sum_i \frac{\beta_i}{x-x_i} } \f]
     *  where \f$ \beta_i = \alpha (-1)^i\f$ (2nd Berrut inetrpolant)
     *  - \f$ \alpha = 1\f$ fo the first and last points
     *  - \f$ \alpha = 2\f$ otherwise 
     *  It is barycentric rational interpolation.
     * - performance is about O(n)
     *  @date 2021-11-25
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     */
    class Berrut2nd : public Interpolation::Table 
    {
    public:
      // ======================================================================
      /** simple constructor from the interpolation points  
       *  @param p input vector of abscissas  
       */
      Berrut2nd 
      ( const  Interpolation::Table& p ) ;
      // ====================================================================
      /** simple constructor from the interpolation points  
       *  @param p input vector of abscissas  
       */
      Berrut2nd 
      (        Interpolation::Table&& p ) ;
      // ====================================================================
      // other constructors from the base class 
      // ======================================================================
      using Interpolation::Table::Table ;
      // ======================================================================
    public:
      // ======================================================================
      /// the main method: get the value of interpolant 
      double evaluate    ( const  double x ) const { return berrut2nd ( x ) ; }      
      /// the main method: get the value of interpolant 
      double operator () ( const  double x ) const { return evaluate  ( x ) ; }
      // ======================================================================
      /// get the weight 
      double weight ( const unsigned short index ) const
      {
        const unsigned int s = size() ;
        return s <= index ? 0 :  
          ( ( 0 == index || index + 1 == s ) ? 1  : 2 ) * ( 0 == index % 2 ? 1 : -1 ) ;
      }
      // ======================================================================
    public : 
      // ======================================================================
      /// swap two interpolators 
      void exchange ( Berrut2nd& right )
      { Interpolation::Table::exchange ( right ) ; }
      // ======================================================================      
    } ;    
    // ========================================================================
    /// swap two Berrut2nd interpolators 
    inline void swap ( Berrut2nd& a , Berrut2nd& b ) { a.exchange ( b ) ; }
    // ========================================================================
    /** @class FloaterHormann
     *  Efficient Floater-Hormann interpolant for barycentric rational interpolation
     *  \f[ F^d_n(x) = \frac{ \sum_i \frac{\beta_i}{x-x_i} f_i}
     *                      { \sum_i \frac{\beta_i}{x-x_i} } \f]
     *  - for \f$d=0\f$, it corresponds to 1st Berrut interpolant 
     *  - for \f$d\ge n\f$ , it corresponds to true barycentric polymnomial onterpolation
     *  - for small values of d it behaves reasonably even for weird interpolaiton mesh 
     *  - performance is about O(n), initialization about \f$ O(nd^2) \f$ 
     *  @date 2021-11-16
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     */
    class FloaterHormann : public Interpolation::Table 
    {
      // ======================================================================
    public:
      // ======================================================================
      /** simple constructor from the interpolation points  
       *  @param p input data 
       *  @param d Floater-Hormann degree parameter 
       */
      FloaterHormann
      ( const Interpolation::Table& p     , 
        const unsigned short        d = 3 ) ;
      // ======================================================================
      /** simple constructor from the interpolation points  
       *  @param p input data 
       *  @param d Floater-Hormann degree parameter 
       */
      FloaterHormann
      ( Interpolation::Table&&      p     , 
        const unsigned short        d = 3 ) ;
      // ======================================================================
      /// default constructor
      FloaterHormann () = default ;
      // =====================================================================
    public:
      // ======================================================================
      /// the main method: get the value of Floater-Hormann interpolant 
      double evaluate    ( const  double x ) const ;
      /// the main method: get the value of interpolant 
      double operator () ( const  double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the degree for the  Floater-Hormann interpolant 
      unsigned short d () const { return m_d ; }
      // ======================================================================
      //// get the weights 
      double weight ( const unsigned short index ) const
      { return index < m_weights.size() ? m_weights[index] : 0 ; }
      // ======================================================================
      /// get the weights 
      const Interpolation::Abscissas::Data& weights () const 
      { return m_weights ; }
      // ======================================================================
    public : 
      // ======================================================================
      /// swap two interpolators 
      void exchange ( FloaterHormann& right ) 
      { 
        Interpolation::Table::exchange ( right ) ;
        swap      ( m_weights , right.m_weights ) ;        
        std::swap ( m_d       , right.m_d       ) ;
      }
      // ======================================================================      
    private :
      // ======================================================================
      /// calculate weigthts for Floater-Hormann interpolant 
      void get_weights() ;
      // ======================================================================      
    private :
      // ======================================================================
      /// degree for the Floater-Hormann interpolant  
      unsigned short                  m_d       ; // Floater Hormann parameter 
      /// vector of weights 
      Interpolation::Abscissas::Data  m_weights ; // vector of weights 
      // ======================================================================
    } ;
    // ========================================================================
    /// swap two Floater-Hormann interpolators 
    inline void swap ( FloaterHormann& a , FloaterHormann& b ) { a.exchange ( b ) ; }
    // ========================================================================
    /** @class Barycentric 
     *  Very efficient (true) Barycentric Lagrange interpolation
     *  For intialization the barycentric weigths are calculated as:
     *  - O(n) for <code>Chebyshev</code> and <code>Chebvyshev2</code>  abscissas 
     *  - O(n) for <code>Uniform</code> abscissas, but much slower 
     *  - O(n^2) for general case 
     *  @see Ostap::Math::Interpolation::Abscissas
     *  @see Ostap::Math::Interpolation::Abscissas::AType
     *  For evaluation it takes O(n) - that is very fast!
     *  @see Jean-Paul Berrut and Lloyd N. Trefethen, 
     *       Barycentric Lagrange Interpolation, SIAM Rev., 46(3), 501–517.
     *       ISSN (print): 0036-1445
     *       ISSN (online): 1095-7200
     *  @see https://doi.org/10.1137/S0036144502417715
     *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
     *  @see https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
     */
    class Barycentric : public Interpolation::Table 
    {
    public:
      // ======================================================================
      /** simple constructor from the interpolation points  
       *  @param p input vector of abscissas  
       */
      Barycentric 
      ( const Interpolation::Table& p ) ;
      /** simple constructor from the interpolation points  
       *  @param p input data 
       */
      Barycentric 
      ( Interpolation::Table&&      p ) ;
      // ======================================================================
      /// default constructor 
      Barycentric () = default ;
      // ======================================================================
    public:
      // ======================================================================
      /// the main method: get the value of interpolant 
      double evaluate    ( const  double x ) const;
      /// the main method: get the value of interpolant 
      double operator () ( const  double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the weights 
      double weight ( const unsigned short index ) const
      { return index < m_weights.size() ? m_weights[index] : 0.0 ; }
      // ======================================================================
      /// get the weights 
      const Interpolation::Abscissas::Data& weights () const 
      { return m_weights ; }
      // ======================================================================
    public : 
      // ======================================================================
      /// swap two interpolators 
      void exchange ( Barycentric& right ) 
      { 
        Interpolation::Table::exchange ( right ) ;
        std::swap ( m_weights , right.m_weights ) ;        
      }
      // ======================================================================      
    private :
      // ======================================================================
      /// calculate weigthts for Barycentric interpolant 
      void get_weights() ;
      // ======================================================================      
    private :
      // ======================================================================
      /// vector of weights
      Interpolation::Abscissas::Data m_weights {} ; // vector of weights
      // ======================================================================
    } ;
    // ========================================================================
    /// swap two Floater-Hormann interpolators 
    inline void swap ( Barycentric& a , Barycentric& b ) 
    { a.exchange ( b ) ; }
    // ========================================================================
    /** @class Newton 
     *  Interpolation polynomial
     *  https://en.wikipedia.org/wiki/Newton_polynomial
     *  @attention it is efficient  and relatively stable numerically
     */
    class Newton : public Interpolation::Table 
    {
    public:
      // ======================================================================
      /** simple constructor from the interpolation points  
       *  @param x input vector of abscissas  
       */
      Newton
      ( const Interpolation::Table& p ) ;
      /** simple constructor from the interpolation points  
       *  @param x input vector of abscissas  
       */
      Newton
      ( Interpolation::Table&&      p ) ;
      // ======================================================================
      /// default constructor
      Newton () = default ;
      // ======================================================================
    public:
      // ======================================================================
      /// the main method: get the value of interpolant 
      double evaluate    ( const  double x ) const ;
      /// the main method: get the value of interpolant 
      double operator () ( const  double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public : 
      // ======================================================================
      /// swap two interpolators 
      void exchange ( Newton& right ) 
      { 
        Interpolation::Table::exchange ( right ) ;
        swap ( m_diffs , right.m_diffs ) ;        
      }
      // ======================================================================      
    private :
      // ======================================================================
      /// calculate differences for Newton interpolation
      void get_differences () ;
      // ======================================================================      
    private :
      // ======================================================================
      /// vector of weights/differences 
      Interpolation::Abscissas::Data m_diffs {} ; // vector of weights 
      // ======================================================================
    } ;
    // ========================================================================
    /// swap two Newton interpolators 
    inline void swap ( Newton& a , Newton& b ) { a.exchange ( b ) ; }
    // ========================================================================
    namespace Interpolation 
    {
      // ======================================================================
      /** Efficient (true) Barycentric Lagrange Interpolation: O(n)
       *  @code 
       *  auto fun             = [] ( double x ) { return std::sin(x) ; } ;
       *  Abscissas abscissas  = ... ;
       *  auto interpolant     = lagrange ( fun , abscissas ) ;
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param abscissas interpolation abscissas  
       *  @see Jean-Paul Berrut and Lloyd N. Trefethen, 
       *       Barycentric Lagrange Interpolation, SIAM Rev., 46(3), 501–517.
       *       ISSN (print): 0036-1445
       *       ISSN (online): 1095-7200
       *  @see https://doi.org/10.1137/S0036144502417715
       *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
       *  @see https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
       *  @see Ostap::Math::Barycentric
       */
      template <class FUNCTION>
      inline  
      Ostap::Math::Barycentric
      lagrange
      ( FUNCTION         func      ,  
        const Abscissas& abscissas )
      { return Ostap::Math::Barycentric ( Table ( func , abscissas ) ) ; }
      // ======================================================================
      /** Efficient (true) Barycentric Lagrange Interpolation: 
       *  @code 
       *  auto fun             = [] ( double x ) { return std::sin(x) ; } ;
       *  auto interpolant     = lagrange ( fun , 12 , 0. , 1. , Abscissas::Chebyshev ) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param N         number of interpolation absicssas
       *  @param low       low  edge of interpolating ranage 
       *  @param high      high edge of interpolating ranage 
       *  @param i         interpolation type       
       *  @see Jean-Paul Berrut and Lloyd N. Trefethen, 
       *       Barycentric Lagrange Interpolation, SIAM Rev., 46(3), 501–517.
       *       ISSN (print): 0036-1445
       *       ISSN (online): 1095-7200
       *  @see https://doi.org/10.1137/S0036144502417715
       *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
       *  @see https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
       *  @see Ostap::Math::Barycentric
       */
      template <class FUNCTION>
      inline 
      Ostap::Math::Barycentric
      lagrange
      ( FUNCTION               func ,  
        const unsigned int     N    , 
        const double           low  ,   
        const double           high , 
        const Abscissas::AType t    )        
      { return Ostap::Math::Barycentric ( Table ( func , N , low , high , t ) )  ; }      
      // ======================================================================
      /** Efficient (true) Barycentric Lagrange Interpolation: 
       *  @code 
       *  auto fun             = [] ( double x ) { return std::sin(x) ; } ;
       *  auto interpolant     = lagrange ( fun , {{ 0.0 , 0.1 , 0.2,  0.3, 0.7, 1.0 }} ) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param x         interpolation abscissas 
       *  @see Jean-Paul Berrut and Lloyd N. Trefethen, 
       *       Barycentric Lagrange Interpolation, SIAM Rev., 46(3), 501–517.
       *       ISSN (print): 0036-1445
       *       ISSN (online): 1095-7200
       *  @see https://doi.org/10.1137/S0036144502417715
       *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
       *  @see https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
       *  @see Ostap::Math::Barycentric
       */      
      template <class FUNCTION>
      inline 
      Ostap::Math::Barycentric
      lagrange 
      ( FUNCTION                                            func ,  
        const  Ostap::Math::Interpolation::Abscissas::Data& x    )
      { return Ostap::Math::Barycentric ( Table ( func , x ) ) ; }
      // ======================================================================
      /** Efficient (true) Barycentric Lagrange Interpolation: 
       *  @code 
       *  Table          table = ... ; // interpolation table 
       *  auto interpolant     = lagrange ( table ) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param table    interpolation table 
       *  @see Jean-Paul Berrut and Lloyd N. Trefethen, 
       *       Barycentric Lagrange Interpolation, SIAM Rev., 46(3), 501–517.
       *       ISSN (print): 0036-1445
       *       ISSN (online): 1095-7200
       *  @see https://doi.org/10.1137/S0036144502417715
       *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
       *  @see https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
       *  @see Ostap::Math::Barycentric
       */      
      inline 
      Ostap::Math::Barycentric
      lagrange
      ( const Table& data )   
      { return Ostap::Math::Barycentric ( data ) ; }
      // ======================================================================      
      /** Efficient (true) Barycentric Lagrange Interpolation: 
       *  @code 
       *  TABLE          table = ... ; // interpolation table 
       *  auto interpolant     = lagrange ( table ) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param table    interpolation table 
       *  @param sorted indicate if data is sorted and duplicates are removed 
       *  @see Jean-Paul Berrut and Lloyd N. Trefethen, 
       *       Barycentric Lagrange Interpolation, SIAM Rev., 46(3), 501–517.
       *       ISSN (print): 0036-1445
       *       ISSN (online): 1095-7200
       *  @see https://doi.org/10.1137/S0036144502417715
       *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
       *  @see https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
       *  @see Ostap::Math::Barycentric
       */
      inline 
      Ostap::Math::Barycentric
      lagrange
      ( const TABLE& data           , 
        const bool   sorted = false )   
      { return Ostap::Math::Barycentric ( Table ( data , sorted ) ) ; }
      // ======================================================================
      /** Efficient (true) Barycentric Lagrange Interpolation: 
       *  @code
       *  typdef std::map<double,double> MAP ;
       *  MAP data ;
       *  data[1.0] = std::sin ( 1.0 ) ;
       *  data[1.5] = std::sin ( 1.5 ) ;
       *  data[2.0] = std::sin ( 2.0 ) ;
       *  data[2.5] = std::sin ( 2.5 ) ;
       *  auto interpolant     = lagrange ( data) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param data      interpolation table 
       *  @see Jean-Paul Berrut and Lloyd N. Trefethen, 
       *       Barycentric Lagrange Interpolation, SIAM Rev., 46(3), 501–517.
       *       ISSN (print): 0036-1445
       *       ISSN (online): 1095-7200
       *  @see https://doi.org/10.1137/S0036144502417715
       *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
       *  @see https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
       *  @see Ostap::Math::Barycentric
       */
      template <class KEY,class VALUE>
      inline 
      Ostap::Math::Barycentric
      lagrange 
      ( const  std::map<KEY,VALUE>& data )   
      { return Ostap::Math::Barycentric ( Table ( data ) ) ; }
      // ======================================================================
      /** Efficient (true) Barycentric Lagrange Interpolation: 
       *  @code
       *  typedef std::vector<double> DATA ;
       *  DATA data {{ std::sin( 0.1 ) , std::sin( 0.1 ) , std::sin( 0.3) }} ;
       *  Abscissas a          = ... ;
       *  auto interpolant     = lagrange ( a , data ) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param x     interpolation abscissas 
       *  @param y     interpolation data 
       *  @see Jean-Paul Berrut and Lloyd N. Trefethen, 
       *       Barycentric Lagrange Interpolation, SIAM Rev., 46(3), 501–517.
       *       ISSN (print): 0036-1445
       *       ISSN (online): 1095-7200
       *  @see https://doi.org/10.1137/S0036144502417715
       *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
       *  @see https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
       *  @see Ostap::Math::Barycentric
       */
      inline   
      Ostap::Math::Barycentric
      lagrange
      ( const  Abscissas&                                   x , 
        const  Ostap::Math::Interpolation::Abscissas::Data& y )
      { return Ostap::Math::Barycentric ( Table ( x , y ) ) ; }
      // ======================================================================      
      /** Efficient (true) Barycentric Lagrange Interpolation: 
       *  @code
       *  typedef std::vector<double> DATA ;
       *  DATA xx {{           0.1   ,           0.2   ,            0.3   }} ;
       *  DATA yy {{ std::sin( 0.1 ) , std::sin( 0.2 ) , std::sin ( 0.3 ) }} ;
       *  auto interpolant     = lagrange ( xx , yy) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param x     interpolation abscissas 
       *  @param y     interpolation data 
       *  @see Jean-Paul Berrut and Lloyd N. Trefethen, 
       *       Barycentric Lagrange Interpolation, SIAM Rev., 46(3), 501–517.
       *       ISSN (print): 0036-1445
       *       ISSN (online): 1095-7200
       *  @see https://doi.org/10.1137/S0036144502417715
       *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
       *  @see https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
       *  @see Ostap::Math::Barycentric
       */
      inline 
      Ostap::Math::Barycentric
      lagrange 
      ( const  Ostap::Math::Interpolation::Abscissas::Data& x ,
        const  Ostap::Math::Interpolation::Abscissas::Data& y )
      { return Ostap::Math::Barycentric ( Table ( x , y ) ) ; }
      // ======================================================================
      
      // ======================================================================      
      //  Newton interpolation 
      // ======================================================================      
      /** Newton interpolation
       *  @code 
       *  auto fun             = [] ( double x ) { return std::sin(x) ; } ;
       *  Abscissas abscissas  = ... ;
       *  auto interpolant     = newton ( fun , abscissas ) ;
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param abscissas interpolation abscissas  
       *  @see https://en.wikipedia.org/wiki/Newton_polynomial
       *  @see Ostap::Math::Newton 
       */
      template <class FUNCTION>
      inline 
      Ostap::Math::Newton
      newton 
      ( FUNCTION         func      ,  
        const Abscissas& abscissas )
      { return Ostap::Math::Newton ( Table ( func , abscissas ) ) ; }
      // =======================================================================
      /** Newton interpolation
       *  @code 
       *  auto fun             = [] ( double x ) { return std::sin(x) ; } ;
       *  auto interpolant     = newton ( fun , 12 , 0. , 1. , Abscissas::Chebyshev ) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param N         number of interpolation absicssas
       *  @param low       low  edge of interpolating ranage 
       *  @param high      high edge of interpolating ranage 
       *  @param i         interpolation type       
       *  @see https://en.wikipedia.org/wiki/Newton_polynomial
       *  @see Ostap::Math::Newton 
       */
      template <class FUNCTION>
      inline 
      Ostap::Math::Newton
      newton 
      ( FUNCTION               func ,  
        const unsigned int     N    , 
        const double           low  ,   
        const double           high , 
        const Abscissas::AType t    )        
      { return Ostap::Math::Newton ( Table ( func , N , low , high , t ) )  ; }
      // ======================================================================      
      /** Newton interpolation
       *  @code 
       *  auto fun             = [] ( double x ) { return std::sin(x) ; } ;
       *  auto interpolant     = newton ( fun , {{ 0.0 , 0.1 , 0.2,  0.3, 0.7, 1.0 }} ) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param x         interpolation abscissas 
       *  @see https://en.wikipedia.org/wiki/Newton_polynomial
       *  @see Ostap::Math::Newton 
       */
      template <class FUNCTION>
      inline 
      Ostap::Math::Newton
      newton 
      ( FUNCTION               func ,  
        const Abscissas::Data& x    )
      { return Ostap::Math::Newton ( Table ( func , x ) ) ; }
      // ======================================================================      
      /** Newton interpolation
       *  @code 
       *  Table          table = ... ; // interpolation table 
       *  auto interpolant     = newton ( table ) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param table    interpolation table  
       *  @see https://en.wikipedia.org/wiki/Newton_polynomial
       *  @see Ostap::Math::Newton 
       */
      inline 
      Ostap::Math::Newton
      Newton 
      ( const Table& data )   
      { return Ostap::Math::Newton ( data ) ; }
      // ======================================================================      
      /** Newton interpolation
       *  @code 
       *  TABLE          table = ... ; // interpolation table 
       *  auto interpolant     = newton ( table ) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param table    interpolation table 
       *  @param sorted indicate if data is sorted and duplicates are removed 
       *  @see https://en.wikipedia.org/wiki/Newton_polynomial
       *  @see Ostap::Math::Newton 
       */
      inline 
      Ostap::Math::Newton
      newton
      ( const TABLE& data , 
        const bool sorted = false )   
      { return Ostap::Math::Newton ( Table ( data , sorted ) ) ; }
      // ======================================================================
      /** Newton interpolation
       *  @code
       *  typdef std::map<double,double> MAP ;
       *  MAP data ;
       *  data[1.0] = std::sin ( 1.0 ) ;
       *  data[1.5] = std::sin ( 1.5 ) ;
       *  data[2.0] = std::sin ( 2.0 ) ;
       *  data[2.5] = std::sin ( 2.5 ) ;
       *  auto interpolant     = newton ( data) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param data      interpolation table 
       *  @see https://en.wikipedia.org/wiki/Newton_polynomial
       *  @see Ostap::Math::Newton 
       */
      template <class KEY,class VALUE>
      inline 
      Ostap::Math::Newton
      newton 
      ( const std::map<KEY,VALUE>& data )   
      { return Ostap::Math::Newton ( Table ( data ) ) ; }
      // ======================================================================
      /** Newton interpolation
       *  @code
       *  typedef std::vector<double> DATA ;
       *  DATA data {{ std::sin( 0.1 ) , std::sin( 0.1 ) , std::sin( 0.3) }} ;
       *  Abscissas a          = ... ;
       *  auto interpolant     = lagrange ( a , data ) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param x     interpolation abscissas 
       *  @param y     interpolation data 
       *  @see https://en.wikipedia.org/wiki/Newton_polynomial
       *  @see Ostap::Math::Newton 
       */
      inline   
      Ostap::Math::Newton
      newton
      ( const Abscissas&       x , 
        const Abscissas::Data& y )
      { return Ostap::Math::Newton ( Table ( x , y ) ) ; }
      // ======================================================================
      /** Newton interpolation
       *  @code
       *  typedef std::vector<double> DATA ;
       *  DATA xx {{           0.1   ,           0.2   ,            0.3   }} ;
       *  DATA yy {{ std::sin( 0.1 ) , std::sin( 0.2 ) , std::sin ( 0.3 ) }} ;
       *  Weights w          = ... ;
       *  auto interpolant     = lagrange ( xx , yy) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param x     interpolation abscissas 
       *  @param y     interpolation data 
       *  @see https://en.wikipedia.org/wiki/Newton_polynomial
       *  @see Ostap::Math::Newton 
       */
      inline 
      Ostap::Math::Newton
      newton
      ( const Abscissas::Data& x , 
        const Abscissas::Data& y )
      { return Ostap::Math::Newton ( Table ( x , y ) ) ; }
      // ======================================================================
    }
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_INTERPOLANTS_H
// ============================================================================
