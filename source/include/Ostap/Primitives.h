// ============================================================================
#ifndef OSTAP_PRIMITIVES_H 
#define OSTAP_PRIMITIVES_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <functional>
// ============================================================================
namespace  Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /** @class Linear
     *  Linear combination of two functions 
     *   \f[ f(x) =  c_1 f_1(x) + c_2 * f_2 ( x  ) \f] 
     *  With the approproate choice of \f$ c_1 \f$  and \f$c_2 \f$  one gest
     *   - sum
     *   - difference 
     *   - scaling 
     *   - bias 
     */
    class Linear
    {
    public :
      // ======================================================================
      /** constructor from two functions and two scale factors 
       *   \f[ f(x) =  c_1 f_1(x) + c_2 * f_2 ( x  ) \f] 
       *  @param f1  the first function 
       *  @param f1  the secons function 
       *  @param c1 the first  scale parameter
       *  @param c2 the second  scale parameter
       */
      template <class FUNCTION1, class FUNCTION2>
      Linear ( FUNCTION1    f1     ,
               FUNCTION2    f2     ,
               const double c1 = 1 ,
               const double c2 = 1 )
        : m_fun1 ( f1 )
        , m_fun2 ( f2 )
        , m_c1   ( c1 ) 
        , m_c2   ( c2 ) 
      {}
      // ======================================================================
      template <class FUNCTION1>
      Linear ( FUNCTION1    f1     ,
               const double f2     ,
               const double c1 = 1 ,
               const double c2 = 1 )
        : m_fun1 ( f1 )
        , m_fun2 ( [f2]( const double /* x */ ) -> double { return f2 ; } )
        , m_c1   ( c1 ) 
        , m_c2   ( c2 )          
      {}
      // ======================================================================
      template <class FUNCTION2>
      Linear ( const double f1     ,
               FUNCTION2    f2     , 
               const double c1 = 1 ,
               const double c2 = 1 )
        : m_fun1 ( [f1]( const double /* x */ ) -> double { return f1 ; } )
        , m_fun2 ( f2 )
        , m_c1   ( c1 ) 
        , m_c2   ( c2 )          
      {}
      // ======================================================================
      /// the main method
      inline  double operator() ( const double x ) const
      { return m_c1 * m_fun1 ( x ) + m_c2 * m_fun2 ( x ) ; }
      // ======================================================================
    public:
      //  =====================================================================
      template <typename ...  ARGS>
      inline static Linear
      create ( ARGS... args ) { return Linear  ( args ... ) ; }
      // ======================================================================
    private :
      //  =====================================================================
      /// the first function 
      std::function<double(double)> m_fun1 {   } ; // the first fuction
      /// the second function 
      std::function<double(double)> m_fun2 {   } ; // the second  fuction
      /// c1-parameter
      double                        m_c1   { 1 } ;  // c1-parameter
      /// c2-parameter
      double                        m_c2   { 1 } ;  // c2-parameter
      //  =====================================================================      
    } ;
    // ========================================================================
    /** @class Compose 
     *  Composition of two functions 
     *   \f[ f(x) =  c_1 f_1(x) ( c_2 * f_2 ( x ) ) \f] 
     */
    class Compose 
    {
    public :
      // ======================================================================
      template <class FUNCTION1, class FUNCTION2>
      Compose ( FUNCTION1    f1     ,
                FUNCTION2    f2     ,
                const double c1 = 1 ,
                const double c2 = 1 )
        : m_fun1 ( f1 )
        , m_fun2 ( f2 )
        , m_c1   ( c1 ) 
        , m_c2   ( c2 ) 
      {}
      // ======================================================================
      /// the main method
      inline  double operator() ( const double x ) const
      { return m_c1 * m_fun1 ( m_c2 * m_fun2 ( x ) ) ; }
      // ======================================================================
    public:
      // ======================================================================
      template <typename ...  ARGS>
      inline static Compose 
      create ( ARGS... args ) { return Compose  ( args ... ) ; }
      // ======================================================================
    private :
      //  =====================================================================
      /// the first function 
      std::function<double(double)> m_fun1 {   } ; // the first fuction
      /// the second function 
      std::function<double(double)> m_fun2 {   } ; // the second  fuction
      /// c1-parameter
      double                        m_c1   { 1 } ;  // c1-parameter
      /// c2-parameter
      double                        m_c2   { 1 } ;  // c2-parameter
      //  =====================================================================      
    } ;
    // ========================================================================
    /** @class Multiply 
     *  Multiplication of two functions 
     *   \f[ f(x) =  c f_1(x) f_2 ( x  ) \f] 
     */
    class Multiply 
    {
    public :
      // ======================================================================
      template <class FUNCTION1, class FUNCTION2>
      Multiply ( FUNCTION1    f1     ,
                 FUNCTION2    f2     ,
                 const double c  = 1 ) 
        : m_fun1 ( f1 )
        , m_fun2 ( f2 )
        , m_c    ( c  ) 
      {}
      // ======================================================================
      /// the main method
      inline  double operator() ( const double x ) const
      { return m_c * m_fun1 ( x ) * m_fun2 ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      template <typename ...  ARGS>
      inline static Multiply 
      create ( ARGS... args ) { return Multiply  ( args ... ) ; }
      // ======================================================================
    private :
      //  =====================================================================
      /// the first function 
      std::function<double(double)> m_fun1 ; // the first fuction
      /// the second function 
      std::function<double(double)> m_fun2 ; // the second  fuction
      /// c-parameter
      double                        m_c   { 1 } ;  // c-parameter
      //  =====================================================================      
    } ;
    // ========================================================================
    /** @class Divide 
     *  Division of two functions 
     *   \f[ f(x) =  c \frac{ f_1(x)}{ f_2 ( x ) }  \f] 
     */
    class Divide 
    {
    public :
      // ======================================================================
      template <class FUNCTION1, class FUNCTION2>
      Divide   ( FUNCTION1    f1     ,
                 FUNCTION2    f2     ,
                 const double c  = 1 ) 
        : m_fun1 ( f1 )
        , m_fun2 ( f2 )
        , m_c    ( c  ) 
      {}
      // ======================================================================
      /// the main method
      inline  double operator() ( const double x ) const
      { return m_c * m_fun1 ( x ) / m_fun2 ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      template <typename ...  ARGS>
      inline static Divide 
      create ( ARGS... args ) { return Divide  ( args ... ) ; }
      // ======================================================================
    private :
      //  =====================================================================
      /// the first function 
      std::function<double(double)> m_fun1 ; // the first fuction
      /// the second function 
      std::function<double(double)> m_fun2 ; // the second  fuction
      /// c-parameter
      double                        m_c   { 1 } ;  // c-parameter
      //  =====================================================================      
    } ;
    // ========================================================================
    /** @class Moebius
     *  Moebius tranformation  \f$ f(x) = \frac{ax+b}{cx+d}\f$,with
     *  \f$ ad-bx \neq 0\f$  
     *  Depending on the parameters one gets:
     *   - linear function, \f$ c=0\f$
     *   - scaling
     *   - constant function     
     *   - hyperbola 
     *  @see https://en.wikipedia.org/wiki/M%C3%B6bius_transformation
     */
    class Moebius
    {
      // ======================================================================
    public:
      // ======================================================================
      Moebius ( const double a = 1  ,
                const double b = 0  ,
                const double c = 0  ,
                const double d = 1  ) ;
      // ======================================================================
      /// the only important method 
      inline double operator()  ( const double x ) const
      { return ( m_a *  x + m_b ) / ( m_c  * x  + m_d ) ;}        
      // ======================================================================
    private :
      // ======================================================================
      /// a-parameter
      double m_a { 1.0 } ;  // a-parameter
      /// b-parameter
      double m_b { 0.0 } ;  // b-parameter
      /// c-parameter
      double m_c { 0.0 } ;  // c-parameter
      /// d-parameter
      double m_d { 1.0 } ;  // d-parameter
      // ======================================================================
    };
    // ========================================================================
    /** @class Step
     *  Heaviside step function
     *   \f[ h(x) =  \left \begin{array}{lcl} 1 & \mathrm{for} & ax+b\ge 0 \\
     *                                        0 & \mathrm{for} & ax+b<0   
     *                       \end{array} \right\} \f] 
     *  @see https://en.wikipedia.org/wiki/Heaviside_step_function
     */
    class Step
    {
    public:
      // ======================================================================
      Step ( const double a = 1 ,
             const double b = 0 )
        : m_a ( a ) 
        , m_b ( 0 )
      {}
      // ======================================================================
      /// the only one method
      inline double operator() ( const double x ) const
      { return  0 <= m_a * x + m_b ? 1 : 0 ; }
      // ======================================================================
    private:
      // ======================================================================
      /// a-parameter
      double m_a { 1 } ; // a-parameter 
      /// b-parameter
      double m_b { 0 } ; // b-parameter 
      // ======================================================================
    } ;
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_PRIMITIVES_H
// ============================================================================
