// ============================================================================
#ifndef OSTAP_PRIMITIVES_H 
#define OSTAP_PRIMITIVES_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <algorithm>
#include <functional>
#include <cmath>
// ============================================================================
namespace  Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /** @class Const
     *  Constant "function": \f$ f(x) \equiv  c \f$ 
     */
    class Const
    {
    public :
      // ======================================================================
      Const
      ( const double c = 0 )
        : m_c ( c )
      {}
      // ======================================================================
      /// the main method
      inline  double operator() ( const double /* x */ ) const { return m_c ; }
      // ======================================================================
    private :
      //  =====================================================================
      /// c-parameter
      double m_c   { 0 } ;  // c-parameter
      //  =====================================================================      
    } ;
    // ========================================================================
    /** @class Id
     *  \f$ f(x) \equiv  x \f$ 
     */
    class Id 
    {
    public :
      // ======================================================================
      /// the main method
      inline  double operator() ( const double  x ) const { return x ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class F1 
     *  holder for one simple function
     */
    class F1
    {
    public :
      // ======================================================================
      /// constructor from the function
      F1  ( std::function<double(double)> f )
        : m_fun ( f )
      {}
      // ======================================================================
      /// constructor from the constant
      F1  ( const double v )
        : m_fun ( Const ( v )  )
      {}
      // ======================================================================
      /// template constructor 
      template <class FUNCTION>
      F1  ( FUNCTION f )
        : m_fun ( f )
      {}
      // ======================================================================
    protected:
      // ======================================================================      
      /// the function itself
      std::function<double(double)> m_fun ; // the function itself
      // ======================================================================      
    } ;
    // ========================================================================
    /** @class F2 
    *   holder for two simple functions
    */
    class F2
    {
    public :
      // ======================================================================
      /// constructor form the function 
      F2
      ( std::function<double(double)> fun1 ,
        std::function<double(double)> fun2 )
        : m_fun1 ( fun1 )
        , m_fun2 ( fun2 )
      {}      
      /// constructor from the function and value  
      F2
      ( std::function<double(double)> fun1 ,
        const double                  fun2 )
        : m_fun1 (         fun1   )
        , m_fun2 ( Const ( fun2 ) )  
      {}
      /// constructor from the function and value  
      F2
      ( const double                  fun1 ,
        std::function<double(double)> fun2 )
        : m_fun1 ( Const ( fun1 ) )  
        , m_fun2 (         fun2   )
      {}
      /// constructor from the values  
      F2
      ( const double                  fun1 ,
        const double                  fun2 )
        : m_fun1 ( Const ( fun1 ) )  
        , m_fun2 ( Const ( fun2 ) )
      {}      
      // ======================================================================
      /// template constructor 
      template <class FUNCTION1>
      F2
      ( FUNCTION1    fun1 ,
        const double fun2  )
        : m_fun1 ( fun1 )
        , m_fun2 ( Const ( fun2 ) ) 
      {}
      // ======================================================================
      /// template constructor 
      template <class FUNCTION2>
      F2
      ( const double fun1 , 
        FUNCTION2    fun2 ) 
        : m_fun1 ( Const ( fun1 ) )  
        , m_fun2 (         fun2   ) 
      {}
      // ======================================================================
      /// template constructor 
      template <class FUNCTION1,
                class FUNCTION2>
      F2
      ( FUNCTION1 fun1 ,
        FUNCTION2 fun2 ) 
        : m_fun1 ( fun1 )
        , m_fun2 ( fun2 )          
      {}
      // ======================================================================
    protected:
      // ======================================================================      
      /// the first function 
      std::function<double(double)> m_fun1 ; // the ferst  function
      /// the second function 
      std::function<double(double)> m_fun2 ; // the second function
      // ======================================================================      
    } ;
    // ========================================================================
    /** @class Linear
     *  Linear combination of two functions 
     *   \f[ f(x) =  c_1 f_1(x) + c_2 * f_2 ( x  ) \f] 
     *  With the approproate choice of \f$ c_1 \f$  and \f$c_2 \f$  one gets
     *   - sum
     *   - difference 
     *   - scaling 
     *   - bias 
     */
    class Linear : public F2 
    {
    public :
      // ======================================================================
      /** constructor from two functions and two scale factors of 1  
       *   \f[ f(x) =  c_1 f_1(x) + c_2 * f_2 ( x  ) \f] 
       *  @param f1  the first function 
       *  @param f2  the secons function 
       */
      Linear
      ( std::function<double(double)> f1 , 
        std::function<double(double)> f2 )
        : F2   ( f1 , f2 ) 
        , m_c1 ( 1  ) 
        , m_c2 ( 1  ) 
      {}
      // ======================================================================
      /** constructor from two functions and two scale factors 
       *  \f[ f(x) =  c_1 f_1(x) + c_2 * f_2 ( x  ) \f] 
       *  @param f1 the first function 
       *  @param c1 the first   scale parameter
       *  @param f2 the second function 
       *  @param c2 the second  scale parameter
       */
      Linear
      ( std::function<double(double)> f1 , 
        const double                  c1 ,
        std::function<double(double)> f2 , 
        const double                  c2 )
        : F2   ( f1 , f2 ) 
        , m_c1 ( c1 ) 
        , m_c2 ( c2 ) 
      {}
      // ======================================================================
      /** constructor for linear function
       *  \f[ f(x) =  ax + b  +\f] 
       *  @param a    \f$ a\f$  coefficient
       *  @param b    \f$ b\f$ coefficient 
       */
      Linear
      ( const double  a , 
        const double  b ) 
        : Linear ( Id() , a , Const(1) , b )
      {}
      // ======================================================================
      /** constructor for linear function passing points 
       *  \f$ (x_1,y_1) \f$ and \f$ (x_2,y_2) \f$
       *  @param x1    \f$ x_1\f$ 
       *  @param y1    \f$ y_1\f$ 
       *  @param x2    \f$ x_2\f$ 
       *  @param y2    \f$ y_2\f$ 
       */
      Linear
      ( const double x1 ,
        const double y1 ,               
        const double x2 ,  
        const double y2 )
        : Linear ( ( y2 - y1 ) / ( x2 - x1 ) , ( x2 * y1 - y2 * x1 ) / ( x2 - x1 ) )
      {}
      // ======================================================================
      /** constructor from two functions and two scale factors of 1  
       *   \f[ f(x) =  c_1 f_1(x) + c_2 * f_2 ( x  ) \f] 
       *  @param f1  the first function 
       *  @param f2  the secons function 
       */
      template <class FUNCTION1, class FUNCTION2>
      Linear
      ( FUNCTION1 f1 , 
        FUNCTION2 f2 )
        : F2   ( f1 , f2 ) 
        , m_c1 ( 1  ) 
        , m_c2 ( 1  ) 
      {}
      // ======================================================================
      /** constructor from two functions and two scale factors 
       *  \f[ f(x) =  c_1 f_1(x) + c_2 * f_2 ( x  ) \f] 
       *  @param f1 the first function 
       *  @param c1 the first   scale parameter
       *  @param f2 the second function 
       *  @param c2 the second  scale parameter
       */
      template <class FUNCTION1, class FUNCTION2>
      Linear
      ( FUNCTION1    f1 ,
        const double c1 ,
        FUNCTION2    f2 ,
        const double c2 )
        : F2   ( f1 , f2 ) 
        , m_c1 ( c1 ) 
        , m_c2 ( c2 ) 
      {}
      // ======================================================================
      /** constructor from functions and scale factors 
       *  \f[ f(x) =  \sum_i c_i f_i(x) +\f] 
       *  @param f1   the first function 
       *  @param c1   the first  scale parameter
       *  @param f2   the second function 
       *  @param c2   the second  scale parameter
       *  @param args other functions/parameters 
       */
      template <class FUNCTION1, 
                class FUNCTION2, 
                typename ... ARGS>
      Linear
      ( FUNCTION1    f1   ,
        const double c1   ,
        FUNCTION2    f2   ,
        const double c2   , 
        ARGS...      args )
        : Linear ( Linear ( f1 , c1 , f2 , c2 ) , 1.0 , args... )
      {}
      // ======================================================================
    public:
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
    class Compose : public F2 
    {
    public :
      // ======================================================================
      Compose
      ( std::function<double(double)> f1     , 
        std::function<double(double)> f2     ,
        const double                  c1 = 1 ,
        const double                  c2 = 1 )
        : F2   ( f1 , f2 ) 
        , m_c1 ( c1 ) 
        , m_c2 ( c2 ) 
      {}
      // ======================================================================
      template <class FUNCTION1, class FUNCTION2>
      Compose
      ( FUNCTION1    f1     ,
        FUNCTION2    f2     ,
        const double c1 = 1 ,
        const double c2 = 1 )
        : F2   ( f1 , f2 ) 
        , m_c1 ( c1 ) 
        , m_c2 ( c2 ) 
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
      /// c1-parameter
      double                        m_c1   { 1 } ;  // c1-parameter
      /// c2-parameter
      double                        m_c2   { 1 } ;  // c2-parameter
      //  =====================================================================      
    } ;
    // ========================================================================
    /** @class Multiply 
     *  Multiplication of two functions 
     *   \f[ f(x) =  f_1(x) \times f_2 ( x  ) \f] 
     */
    class Multiply : public F2 
    {
    public :
      // ======================================================================
      Multiply
      ( std::function<double(double)> f1 , 
        std::function<double(double)> f2 )
        : F2  ( f1 , f2 ) 
      {}
      // ======================================================================
      template <class FUNCTION1, 
                class FUNCTION2>
      Multiply
      ( FUNCTION1    f1     ,
        FUNCTION2    f2     ) 
        : F2  ( f1 , f2 ) 
      {}
      // ====================================================================== 
      /// Constructor with several functions
      template <class FUNCTION1, class FUNCTION2, typename ... ARGS >
      Multiply
      ( FUNCTION1 f1   ,
        FUNCTION2 f2   , 
        ARGS ...  args ) 
        : Multiply ( Multiply ( f1 , f2 ) , args ... )
      {}
      // ======================================================================
      /// the main method
      inline  double operator() ( const double x ) const
      { return m_fun1 ( x ) * m_fun2 ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      template <typename ...  ARGS>
      inline static Multiply 
      create ( ARGS... args ) { return Multiply  ( args ... ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Divide 
     *  Division of two functions 
     *   \f[ f(x) =  \frac{ f_1(x)}{ f_2 ( x ) }  \f] 
     */
    class Divide : public F2  
    {
    public :
      // ======================================================================
      Divide
      ( std::function<double(double)> f1 , 
        std::function<double(double)> f2 )
        : F2 ( f1 , f2 ) 
      {}
      // ======================================================================
      template <class FUNCTION1, class FUNCTION2>
      Divide
      ( FUNCTION1 f1 ,
        FUNCTION2 f2 ) 
        : F2 ( f1 , f2 ) 
      {}
      /// Constructor with several functions
      template <class FUNCTION1, class FUNCTION2, typename ... ARGS >
      Divide
      ( FUNCTION1 f1   ,
        FUNCTION2 f2   , 
        ARGS ...  args ) 
        : Divide ( f1  , Multiply ( f2  , args ... ) )
      {}
      // ======================================================================
      /// the main method
      inline  double operator() ( const double x ) const
      { return m_fun1 ( x ) / m_fun2 ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      template <typename ...  ARGS>
      inline static Divide 
      create ( ARGS... args ) { return Divide  ( args ... ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Sum
     *  simple sum of several functions 
     *  \f[ f(x) = \sum_i f_i(x)  \f] 
     */
    class Sum : public F2 
    {
    public :
      // ======================================================================
      /** constructor from two functions
       *   \f[ f(x) =  f_1(x) + f_2 ( x  ) \f] 
       *  @param f1 the first  function 
       *  @param f1 the second function 
       */
      Sum
      ( std::function<double(double)> f1 , 
        std::function<double(double)> f2 )
        : F2 ( f1 , f2 ) 
      {}
      // ======================================================================
      /** constructor from two functions
       *   \f[ f(x) =  f_1(x) + f_2 ( x  ) \f] 
       *  @param f1 the first  function 
       *  @param f1 the second function 
       */
      template <class FUNCTION1, class FUNCTION2>
      Sum
      ( FUNCTION1    f1 ,
        FUNCTION2    f2 ) 
        : F2 ( f1 , f2 ) 
      {}
      // ====================================================================== 
      /// Constructor with several functions
      template <class FUNCTION1, class FUNCTION2, typename ... ARGS >
      Sum
      ( FUNCTION1 f1   ,
        FUNCTION2 f2   , 
        ARGS ...  args ) 
        : Sum ( Sum ( f1 , f2 ) , args ... )
      {}
      // ======================================================================
    public: 
      // ======================================================================
      /// the main method
      inline  double operator() ( const double x ) const
      { return m_fun1 ( x ) + m_fun2 ( x ) ; }
      // ======================================================================
    public:
      //  =====================================================================
      template <typename ...  ARGS>
      inline static Sum
      create ( ARGS... args ) { return Sum ( args ... ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Subtract
     *  Subtract two or more functions 
     *  \f[ f(x) = f_1 (x) - f_2(x) \f] 
     *  \f[ f(x) = f_1 (x) - \sum_{i=2}} f_i(x) \f] 
     */
    class Subtract: public F2 
    {
    public :
      // ======================================================================
      /** constructor from two functions
       *   \f[ f(x) =  f_1(x) - f_2 ( x  ) \f] 
       *  @param f1 the first  function 
       *  @param f1 the second function 
       */
      Subtract
      ( std::function<double(double)> f1 , 
        std::function<double(double)> f2 )
        : F2 ( f1 , f2 ) 
      {}
      // ======================================================================
      /** constructor from two functions
       *   \f[ f(x) =  f_1(x) - f_2 ( x  ) \f] 
       *  @param f1 the first  function 
       *  @param f1 the second function 
       */
      template <class FUNCTION1, class FUNCTION2>
      Subtract
      ( FUNCTION1    f1 ,
        FUNCTION2    f2 ) 
        : F2 ( f1 , f2 ) 
      {}
      // ====================================================================== 
      /// Constructor with several functions
      template <class FUNCTION1, class FUNCTION2, typename ... ARGS >
      Subtract
      ( FUNCTION1 f1   ,
        FUNCTION2 f2   , 
        ARGS ...  args ) 
        : Subtract ( f1 , Sum ( f2 , args ... ) ) 
      {}
      // ======================================================================
    public: 
      // ======================================================================
      /// the main method
      inline  double operator() ( const double x ) const
      { return m_fun1 ( x ) - m_fun2 ( x ) ; }
      // ======================================================================
    public:
      //  =====================================================================
      template <typename ...  ARGS>
      inline static Subtract
      create ( ARGS... args ) { return Subtract ( args ... ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Moebius
     *  Moebius tranformation  \f$ f(x) = \frac{ax+b}{cx+d}\f$,
     *  with \f$ ad - bc \neq 0\f$  
     *  Depending on the parameters one gets:
     *   - linear function, \f$ c=0\f$
     *   - scaling
     *   - hyperbola , \f$ a =0 \f$ 
     *  @see https://en.wikipedia.org/wiki/M%C3%B6bius_transformation
     */
    class Moebius
    {
      // ======================================================================
    public:
      // ======================================================================
      Moebius
      ( const double a = 1  ,
        const double b = 0  ,
        const double c = 0  ,
        const double d = 1  ) ;
      // ======================================================================
      /// the only important method 
      inline double operator()  ( const double x ) const
      { return ( m_a * x + m_b ) / ( m_c  * x  + m_d ) ; }        
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
     *   \f[ h(x) =  \left \begin{array}{lcl} 1 & \mathrm{for} & ax+b\ge 0 \ \
     *                                        0 & \mathrm{for} & ax+b<0   
     *                       \end{array} \right\} \f] 
     *  @see https://en.wikipedia.org/wiki/Heaviside_step_function
     */
    class Step
    {
    public:
      // ======================================================================
      Step
      ( const double a = 1 ,
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
    /** @class Max
     *  Maximal of two functions 
     *   \f[ f(x) =  max ( f_1(x) , f_2(x) ) \f] 
     */
    class Max : public F2 
    {
    public :
      // ======================================================================
      /** constructor from two functions 
       *   \f[ f(x) =  max  ( f_1(x) , f_2(x) )\f] 
       *  @param f1  the first  function 
       *  @param f2  the second function 
       */
      Max
      ( std::function<double(double)> f1 , 
        std::function<double(double)> f2 )
        : F2 ( f1 , f2 ) 
      {}
      // ======================================================================
      /** constructor from two functions 
       *   \f[ f(x) =  max  ( f_1(x) , f_2(x) )\f] 
       *  @param f1  the first  function 
       *  @param f2  the second function 
       */
      template <class FUNCTION1, class FUNCTION2>
      Max
      ( FUNCTION1    f1 , 
        FUNCTION2    f2 )
        : F2 ( f1 , f2 ) 
      {}
      // ======================================================================
      /** constructor from two functions 
       *   \f[ f(x) =  max  ( f_1(x) , f_2(x) )\f] 
       *  @param f1  the first  function 
       *  @param f2  the second function 
       */
      template <class FUNCTION1>
      Max
      ( FUNCTION1    f1 , 
        const double f2 )
        : F2 ( f1 , Const ( f2 ) ) 
      {}
      // ======================================================================
      /** constructor from two functions 
       *   \f[ f(x) =  max  ( f_1(x) , f_2(x) )\f] 
       *  @param f1  the first  function 
       *  @param f2  the second function 
       */
      template <class FUNCTION2>
      Max
      ( const double f1 , 
        FUNCTION2    f2 ) 
        : F2 ( Const ( f1 ) , f2 ) 
      {}
      // ======================================================================
      /// Constructor with several functions
      template <class FUNCTION1, class FUNCTION2, typename ... ARGS >
      Max
      ( FUNCTION1 f1   ,
        FUNCTION2 f2   , 
        ARGS ...  args ) 
        : Max ( Max ( f1 , f2 ) , args ... )
      {}
      // ======================================================================
    public:
      // ======================================================================
      /// the main method
      inline  double operator() ( const double x ) const
      { return std::max  ( m_fun1 ( x ) , m_fun2 ( x ) ) ; }
      // ======================================================================
    public:
      //  =====================================================================
      template <typename ...  ARGS>
      inline static Max
      create ( ARGS... args ) { return Max  ( args ... ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Min
     *  Minimal value  of two functions 
     *   \f[ f(x) =  min ( f_1(x) , f_2(x) ) \f] 
     */
    class Min : public F2 
    {
    public :
      // ======================================================================
      /** constructor from two functions 
       *   \f[ f(x) =  min  ( f_1(x) , f_2(x) )\f] 
       *  @param f1  the first function 
       *  @param f1  the secons function 
       */
      Min 
      ( std::function<double(double)> f1 , 
        std::function<double(double)> f2 )
        : F2 ( f1  , f2 ) 
      {}
      // ======================================================================
      /** constructor from two functions 
       *   \f[ f(x) =  min  ( f_1(x) , f_2(x) )\f] 
       *  @param f1  the first function 
       *  @param f1  the secons function 
       */
      template <class FUNCTION1, class FUNCTION2>
      Min
      ( FUNCTION1    f1 , 
        FUNCTION2    f2 )
        : F2 ( f1  , f2 ) 
      {}
      // ======================================================================
      /** constructor from two functions 
       *   \f[ f(x) =  min  ( f_1(x) , f_2(x) )\f] 
       *  @param f1  the first function 
       *  @param f1  the secons function 
       */
      template <class FUNCTION1>
      Min
      ( FUNCTION1    f1 ,
        const double f2 )
        : F2 ( f1  , Const ( f2 ) ) 
      {}
      // ======================================================================
      /** constructor from two functions 
       *   \f[ f(x) =  min  ( f_1(x) , f_2(x) )\f] 
       *  @param f1  the first function 
       *  @param f1  the secons function 
       */
      template <class FUNCTION2>
      Min
      ( const double f1 , 
        FUNCTION2    f2 )
        : F2 ( Const ( f1 ) , f2 ) 
      {}
      // ======================================================================
      /// Constructor with several functions
      template <class FUNCTION1, class FUNCTION2, typename ... ARGS >
      Min
      ( FUNCTION1 f1   ,
        FUNCTION2 f2   , 
        ARGS ...  args ) 
        : Min ( Min ( f1 , f2 ) , args ... )
      {}
      // ======================================================================
    public:
      // ======================================================================      
      /// the main method
      inline  double operator() ( const double x ) const
      { return std::min  ( m_fun1 ( x ) , m_fun2 ( x ) ) ; }
      // ======================================================================
    public:
      //  =====================================================================
      template <typename ...  ARGS>
      inline static Min
      create ( ARGS... args ) { return Min  ( args ... ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Apply
     *  keep and apply arbitrary function
     */
    class Apply : public F1 
    {
    public:
      // ======================================================================
      /// constructor  from the function 
      Apply
      ( std::function<double(double)> f )
        : F1 ( f  )
      {}
      // ======================================================================
      /// constrtuctor  from the function 
      template <class FUNCTION>
      Apply ( FUNCTION f )
        : F1 ( f  )
      {}
      /// constant function 
      Apply ( const double a ) 
        : F1 ( a )
      {}
      /// copy constructor 
      Apply ( const Apply&  ) = default ;
      /// move constructor
      Apply (       Apply&& ) = default ;
      // ======================================================================
    public:
      // ======================================================================
      template <class FUNCTION>
      static inline Apply 
      create ( FUNCTION f ) { return Apply ( f ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// the only one important method 
      inline double operator() ( const double x ) const 
      { return m_fun ( x ) ; }
      // ======================================================================        
    } ;
    // ========================================================================
    /** @class Abs
     *  \f$  F(x) = \left| f(x) \right| \f$ 
     */
    class Abs : public F1 
    {
    public:
      // ======================================================================
      Abs 
      ( std::function<double(double)> f )
        : F1 ( f  )
      {}
      // ======================================================================      
      template <class FUNCTION>
      Abs ( FUNCTION f )
        : F1 ( f  )
      {}
      // ======================================================================
      /// default constructor 
      Abs () : Abs ( Id () ) {}
      // ======================================================================
    public:
      // ======================================================================      
      template <class FUNCTION>
      static inline Abs
      create ( FUNCTION f ) { return Abs ( f ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// the only one important method 
      inline double operator() ( const double x ) const
      { return std::abs ( m_fun ( x ) ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Sqrt 
     *  \f$  F(x) = \sqrt { f(x) } \f$ 
     */
    class Sqrt : public F1 
    {
    public:
      // ======================================================================
      Sqrt 
      ( std::function<double(double)> f )
        : F1 ( f  )
      {}
      // ======================================================================
      template <class FUNCTION>
      Sqrt ( FUNCTION f )
        : F1 ( f  )
      {}
      // ======================================================================
      /// default contsructor 
      Sqrt () : Sqrt ( Id () ) {}
      // ======================================================================
    public:
      // ======================================================================      
      template <class FUNCTION>
      static inline Sqrt
      create ( FUNCTION f ) { return Sqrt ( f ) ; }
      // ======================================================================
    public:
      /// the only one important method 
      inline double operator() ( const double x ) const
      { return std::sqrt ( m_fun ( x ) ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Cbrt
     *  \f$  F(x) = \sqrt[3] { f(x) } \f$ 
     */
    class Cbrt : public F1 
    {
    public:
      // ======================================================================
      Cbrt 
      ( std::function<double(double)> f )
        : F1 ( f  )
      {}
      // ======================================================================
      template <class FUNCTION>
      Cbrt ( FUNCTION f )
        : F1 ( f  )
      {}
      // ======================================================================
      /// default contsructor 
      Cbrt () : Cbrt ( Id () ) {}
      // ======================================================================
    public:
      // ======================================================================      
      template <class FUNCTION>
      static inline Cbrt
      create ( FUNCTION f ) { return Cbrt ( f ) ; }
      // ======================================================================
    public:
      /// the only one important method 
      inline double operator() ( const double x ) const
      { return std::cbrt ( m_fun ( x ) ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Exp
     *  \f$  F(x) = e^{f(x)} \f$ 
     */
    class Exp : public F1 
    {
    public:
      // ======================================================================
      Exp
      ( std::function<double(double)> f )
        : F1 ( f  )
      {}
      // ======================================================================
      template <class FUNCTION>
      Exp ( FUNCTION f )
        : F1 ( f  )
      {}
      // ======================================================================
      /// default contsructor 
      Exp () : Exp ( Id () ) {}
      // ======================================================================
    public:
      // ======================================================================      
      template <class FUNCTION>
      static inline Exp
      create ( FUNCTION f ) { return Exp ( f ) ; }
      // ======================================================================
    public:
      /// the only one important method 
      inline double operator() ( const double x ) const
      { return std::exp ( m_fun ( x ) ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Log
     *  \f$  F(x) = \log f(x) \f$ 
     */
    class Log : public F1 
    {
    public:
      // ======================================================================
      Log
      ( std::function<double(double)> f )
        : F1 ( f  )
      {}
      // ======================================================================
      template <class FUNCTION>
      Log ( FUNCTION f )
        : F1 ( f  )
      {}
      // ======================================================================
      /// default contsructor 
      Log () : Log ( Id () ) {}
      // ======================================================================
    public:
      // ======================================================================      
      template <class FUNCTION>
      static inline Log
      create ( FUNCTION f ) { return Log ( f ) ; }
      // ======================================================================
    public:
      /// the only one important method 
      inline double operator() ( const double x ) const
      { return std::log ( m_fun ( x ) ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Log10
     *  \f$  F(x) = \log_{10} f(x) \f$ 
     */
    class Log10 : public F1 
    {
    public:
      // ======================================================================
      Log10
      ( std::function<double(double)> f )
        : F1 ( f  )
      {}
      // ======================================================================
      template <class FUNCTION>
      Log10 ( FUNCTION f )
        : F1 ( f  )
      {}
      // ======================================================================
      /// default contsructor 
      Log10 () : Log10 ( Id () ) {}
      // ======================================================================
    public:
      // ======================================================================      
      template <class FUNCTION>
      static inline Log10
      create ( FUNCTION f ) { return Log10 ( f ) ; }
      // ======================================================================
    public:
      /// the only one important method 
      inline double operator() ( const double x ) const
      { return std::log10 ( m_fun ( x ) ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Erf
     *  \f$  F(x) = erf f(x) \f$ 
     */
    class Erf : public F1 
    {
    public:
      // ======================================================================
      Erf
      ( std::function<double(double)> f )
        : F1 ( f  )
      {}
      // ======================================================================
      template <class FUNCTION>
      Erf ( FUNCTION f )
        : F1 ( f  )
      {}
      // ======================================================================
      /// default contsructor 
      Erf () : Erf ( Id () ) {}
      // ======================================================================
    public:
      // ======================================================================      
      template <class FUNCTION>
      static inline Erf
      create ( FUNCTION f ) { return Erf ( f ) ; }
      // ======================================================================
    public:
      /// the only one important method 
      inline double operator() ( const double x ) const
      { return std::erf ( m_fun ( x ) ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Erfc
     *  \f$  F(x) = erfc f(x) \f$ 
     */
    class Erfc : public F1 
    {
    public:
      // ======================================================================
      Erfc
      ( std::function<double(double)> f )
        : F1 ( f  )
      {}
      // ======================================================================
      template <class FUNCTION>
      Erfc ( FUNCTION f )
        : F1 ( f  )
      {}
      // ======================================================================
      /// default contsructor 
      Erfc () : Erfc ( Id () ) {}
      // ======================================================================
    public:
      // ======================================================================      
      template <class FUNCTION>
      static inline Erfc
      create ( FUNCTION f ) { return Erfc ( f ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// the only one important method 
      inline double operator() ( const double x ) const
      { return std::erfc ( m_fun ( x ) ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class TGamma
     *  \f$  F(x) = \Gamma f(x) \f$ 
     */
    class TGamma : public F1 
    {
    public:
      // ======================================================================
      TGamma
      ( std::function<double(double)> f )
        : F1 ( f  )
      {}
      // ======================================================================
      template <class FUNCTION>
      TGamma ( FUNCTION f )
        : F1 ( f  )
      {}
      // ======================================================================
      /// default contsructor 
      TGamma () : TGamma ( Id () ) {}
      // ======================================================================
    public:
      // ======================================================================      
      template <class FUNCTION>
      static inline TGamma
      create ( FUNCTION f ) { return TGamma ( f ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// the only one important method 
      inline double operator() ( const double x ) const
      { return std::tgamma ( m_fun ( x ) ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class LGamma
     *  \f$  F(x) = \log \Gamma f(x) \f$ 
     */
    class LGamma : public F1 
    {
    public:
      // ======================================================================
      LGamma
      ( std::function<double(double)> f )
        : F1 ( f  )
      {}
      // ======================================================================
      template <class FUNCTION>
      LGamma ( FUNCTION f )
        : F1 ( f  )
      {}
      // ======================================================================
      /// default contsructor 
      LGamma () : LGamma ( Id () ) {}
      // ======================================================================
    public:
      // ======================================================================      
      template <class FUNCTION>
      static inline LGamma
      create ( FUNCTION f ) { return LGamma ( f ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// the only one important method 
      inline double operator() ( const double x ) const
      { return std::lgamma ( m_fun ( x ) ) ; }
      // ======================================================================
    } ;
    
    // ========================================================================
    /** @class Sin
     *  \f$  F(x) = \sin f(x) \f$ 
     */
    class Sin : public F1 
    {
    public:
      // ======================================================================
      Sin
      ( std::function<double(double)> f )
        : F1 ( f  )
      {}
      // ======================================================================
      template <class FUNCTION>
      Sin ( FUNCTION f )
        : F1 ( f  )
      {}
      // ======================================================================
      /// default contsructor 
      Sin () : Sin ( Id () ) {}
      // ======================================================================
    public:
      // ======================================================================      
      template <class FUNCTION>
      static inline Sin
      create ( FUNCTION f ) { return Sin ( f ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// the only one important method 
      inline double operator() ( const double x ) const
      { return std::sin ( m_fun ( x ) ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Cos
     *  \f$  F(x) = \cos f(x) \f$ 
     */
    class Cos : public F1 
    {
    public:
      // ======================================================================
      Cos
      ( std::function<double(double)> f )
        : F1 ( f  )
      {}
      // ======================================================================
      template <class FUNCTION>
      Cos ( FUNCTION f )
        : F1 ( f  )
      {}
      // ======================================================================
      /// default contsructor 
      Cos () : Cos ( Id () ) {}
      // ======================================================================
    public:
      // ======================================================================      
      template <class FUNCTION>
      static inline Cos
      create ( FUNCTION f ) { return Cos ( f ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// the only one important method 
      inline double operator() ( const double x ) const
      { return std::cos ( m_fun ( x ) ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Tan
     *  \f$  F(x) = \tan f(x) \f$ 
     */
    class Tan : public F1 
    {
    public:
      // ======================================================================
      Tan
      ( std::function<double(double)> f )
        : F1 ( f  )
      {}
      // ======================================================================
      template <class FUNCTION>
      Tan ( FUNCTION f )
        : F1 ( f  )
      {}
      // ======================================================================
      /// default contsructor 
      Tan () : Tan ( Id () ) {}
      // ======================================================================
    public:
      // ======================================================================      
      template <class FUNCTION>
      static inline Tan
      create ( FUNCTION f ) { return Tan ( f ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// the only one important method 
      inline double operator() ( const double x ) const
      { return std::tan ( m_fun ( x ) ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class ASin
     *  \f$  F(x) = \asin f(x) \f$ 
     */
    class ASin : public F1 
    {
    public:
      // ======================================================================
      ASin
      ( std::function<double(double)> f )
        : F1 ( f  )
      {}
      // ======================================================================
      template <class FUNCTION>
      ASin ( FUNCTION f )
        : F1 ( f  )
      {}
      // ======================================================================
      /// default contsructor 
      ASin () : ASin ( Id () ) {}
      // ======================================================================
    public:
      // ======================================================================      
      template <class FUNCTION>
      static inline ASin
      create ( FUNCTION f ) { return ASin ( f ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// the only one important method 
      inline double operator() ( const double x ) const
      { return std::asin ( m_fun ( x ) ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class ACos
     *  \f$  F(x) = \acos f(x) \f$ 
     */
    class ACos : public F1 
    {
    public:
      // ======================================================================
      ACos
      ( std::function<double(double)> f )
        : F1 ( f  )
      {}
      // ======================================================================
      template <class FUNCTION>
      ACos ( FUNCTION f )
        : F1 ( f  )
      {}
      // ======================================================================
      /// default contsructor 
      ACos () : ACos ( Id () ) {}
      // ======================================================================
    public:
      // ======================================================================      
      template <class FUNCTION>
      static inline ACos
      create ( FUNCTION f ) { return ACos ( f ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// the only one important method 
      inline double operator() ( const double x ) const
      { return std::acos ( m_fun ( x ) ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class ATan
     *  \f$  F(x) = \atan f(x) \f$ 
     */
    class ATan : public F1 
    {
    public:
      // ======================================================================
      ATan
      ( std::function<double(double)> f )
        : F1 ( f  )
      {}
      // ======================================================================
      template <class FUNCTION>
      ATan ( FUNCTION f )
        : F1 ( f  )
      {}
      // ======================================================================
      /// default contsructor 
      ATan () : ATan ( Id () ) {}
      // ======================================================================
    public:
      // ======================================================================      
      template <class FUNCTION>
      static inline ATan
      create ( FUNCTION f ) { return ATan ( f ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// the only one important method 
      inline double operator() ( const double x ) const
      { return std::atan ( m_fun ( x ) ) ; }
      // ======================================================================
    } ;
    
    // ========================================================================
    /** @class Sinh
     *  \f$  F(x) = \sinh f(x) \f$ 
     */
    class Sinh : public F1 
    {
    public:
      // ======================================================================
      Sinh
      ( std::function<double(double)> f )
        : F1 ( f  )
      {}
      // ======================================================================
      template <class FUNCTION>
      Sinh ( FUNCTION f )
        : F1 ( f  )
      {}
      // ======================================================================
      /// default contsructor 
      Sinh () : Sinh ( Id () ) {}
      // ======================================================================
    public:
      // ======================================================================      
      template <class FUNCTION>
      static inline Sinh
      create ( FUNCTION f ) { return Sinh ( f ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// the only one important method 
      inline double operator() ( const double x ) const
      { return std::sinh ( m_fun ( x ) ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Cosh
     *  \f$  F(x) = \cos f(x) \f$ 
     */
    class Cosh : public F1 
    {
    public:
      // ======================================================================
      Cosh
      ( std::function<double(double)> f )
        : F1 ( f  )
      {}
      // ======================================================================
      template <class FUNCTION>
      Cosh ( FUNCTION f )
        : F1 ( f  )
      {}
      // ======================================================================
      /// default contsructor 
      Cosh () : Cosh ( Id () ) {}
      // ======================================================================
    public:
      // ======================================================================      
      template <class FUNCTION>
      static inline Cosh
      create ( FUNCTION f ) { return Cosh ( f ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// the only one important method 
      inline double operator() ( const double x ) const
      { return std::cosh ( m_fun ( x ) ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Tanh
     *  \f$  F(x) = \tanh f(x) \f$ 
     */
    class Tanh : public F1 
    {
    public:
      // ======================================================================
      Tanh
      ( std::function<double(double)> f )
        : F1 ( f  )
      {}
      // ======================================================================
      template <class FUNCTION>
      Tanh ( FUNCTION f )
        : F1 ( f  )
      {}
      // ======================================================================
      /// default contsructor 
      Tanh () : Tanh ( Id () ) {}
      // ======================================================================
    public:
      // ======================================================================      
      template <class FUNCTION>
      static inline Tanh
      create ( FUNCTION f ) { return Tanh ( f ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// the only one important method 
      inline double operator() ( const double x ) const
      { return std::tanh ( m_fun ( x ) ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class ASinh
     *  \f$  F(x) = \asinh f(x) \f$ 
     */
    class ASinh : public F1 
    {
    public:
      // ======================================================================
      ASinh
      ( std::function<double(double)> f )
        : F1 ( f  )
      {}
      // ======================================================================
      template <class FUNCTION>
      ASinh ( FUNCTION f )
        : F1 ( f  )
      {}
      // ======================================================================
      /// default constructor 
      ASinh () : ASinh ( Id () ) {}
      // ======================================================================
    public:
      // ======================================================================      
      template <class FUNCTION>
      static inline ASinh
      create ( FUNCTION f ) { return ASinh ( f ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// the only one important method 
      inline double operator() ( const double x ) const
      { return std::asinh ( m_fun ( x ) ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class ACosh
     *  \f$  F(x) = \acosh f(x) \f$ 
     */
    class ACosh : public F1 
    {
    public:
      // ======================================================================
      ACosh
      ( std::function<double(double)> f )
        : F1 ( f  )
      {}
      // ======================================================================
      template <class FUNCTION>
      ACosh ( FUNCTION f )
        : F1 ( f  )
      {}
      // ======================================================================
      /// default contsructor 
      ACosh () : ACosh ( Id () ) {}
      // ======================================================================
    public:
      // ======================================================================      
      template <class FUNCTION>
      static inline ACosh
      create ( FUNCTION f ) { return ACosh ( f ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// the only one important method 
      inline double operator() ( const double x ) const
      { return std::acosh ( m_fun ( x ) ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class ATanh
     *  \f$  F(x) = \atanh f(x) \f$ 
     */
    class ATanh : public F1 
    {
    public:
      // ======================================================================
      ATanh
      ( std::function<double(double)> f )
        : F1 ( f  )
      {}
      // ======================================================================
      template <class FUNCTION>
      ATanh ( FUNCTION f )
        : F1 ( f  )
      {}
      // ======================================================================
      /// default contsructor 
      ATanh () : ATanh ( Id () ) {}
      // ======================================================================
    public:
      // ======================================================================      
      template <class FUNCTION>
      static inline ATanh
      create ( FUNCTION f ) { return ATanh ( f ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// the only one important method 
      inline double operator() ( const double x ) const
      { return std::atanh ( m_fun ( x ) ) ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Pow 
     *  pow for the function 
     */
    class Pow : public F2 
    {
    public:
      // ======================================================================
      /// the type 
      enum Type { Integer = 0 ,
                  Double      ,
                  Function    } ;
      // ======================================================================        
    public:
      // ======================================================================
      /// constructor  from the function and degree 
      Pow
      ( std::function<double(double)> f     , 
        const double                  order ) 
        : F2 ( f , Const ( order ) ) 
        , m_iorder ( 0      )
        , m_dorder ( order  )
        , m_type   ( Double ) 
      {}      
      /// constructor  from the function and degree 
      Pow
      ( std::function<double(double)> f     , 
        const int                     order ) 
        : F2 ( f , Const ( 1.0 * order ) ) 
        , m_iorder ( order   )
        , m_dorder ( order   )
        , m_type   ( Integer ) 
      {}
      /// constructor  from two functions
      Pow 
      ( std::function<double(double)> f1 , 
        std::function<double(double)> f2 ) 
        : F2 ( f1 , f2   ) 
        , m_iorder ( 0        )
        , m_dorder ( 0        )
        , m_type   ( Function ) 
      {}
      // =======================================================================
      /// constructor  from the function and degree 
      template <class FUNCTION>
      Pow 
      ( FUNCTION     f     ,
        const double order ) 
        : F2 ( f , Const ( order ) ) 
        , m_iorder ( 0      )
        , m_dorder ( order  )
        , m_type   ( Double ) 
      {}
      /// constructor  from the function and degree 
      template <class FUNCTION>
      Pow
      ( FUNCTION     f     ,
        const int    order ) 
        : F2 ( f , Const ( 1.0 * order ) ) 
        , m_iorder ( order   )
        , m_dorder ( order   )
        , m_type   ( Integer ) 
      {}
      
      /// constructor  from two functions
      template <class FUNCTION1,
                class FUNCTION2>
      Pow 
      ( FUNCTION1    f1    ,
        FUNCTION2    f2    ) 
        : F2 ( f1 , f2   ) 
        , m_iorder ( 0        )
        , m_dorder ( 0        )
        , m_type   ( Function ) 
      {}

      /// constructor  from degree 
      Pow ( const double order ) 
        : Pow ( Id () ,  order )
      {}
      /// constructor  from degree 
      Pow ( const int    order ) 
        : Pow ( Id () ,  order )
      {}
      /// copy constructor 
      Pow ( const Pow&  ) = default ;
      /// move constructor
      Pow (       Pow&& ) = default ;
      // ======================================================================
    public:
      // ======================================================================
      template <class FUNCTION>
      static inline Pow 
      create
      ( FUNCTION     f ,
        const double o ) { return Pow ( f , o ) ; }
      // ======================================================================
      template <class FUNCTION>
      static inline Pow 
      create
      ( FUNCTION     f ,
        const int    o ) { return Pow ( f , o ) ; }
      // ======================================================================
      template <class FUNCTION1,
                class FUNCTION2>
      static inline Pow 
      create
      ( FUNCTION1    f1 ,
        FUNCTION2    f2 ) { return Pow ( f1 , f2 ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// the only one important method 
      inline double operator() ( const double x ) const
      {
        return
          m_type == Integer ? std::pow ( m_fun1 ( x ) , m_iorder ) :
          m_type == Double  ? std::pow ( m_fun1 ( x ) , m_dorder ) :
          std::pow  ( m_fun1 ( x  ) , m_fun2  ( x ) ) ;
      }
      // ======================================================================        
    private :
      // ======================================================================     
      /// integer order
      int                            m_iorder { 0       } ; // integer order
      /// non-integer order
      double                         m_dorder { 0       } ; // non-integer order
      /// type 
      Type                           m_type   { Integer } ; // type 
      // ======================================================================
    } ;
    // ========================================================================
    /** @class SmoothTransition
     *  smooth transition between two functions 
     *  @see Ostap::Math::smoooth_transition 
     */
    class SmoothTransition : public F2 
    {
    public:
      // ======================================================================
      /// constructor from two functions
      SmoothTransition
      ( std::function<double(double)> f1    , 
        std::function<double(double)> f2    ,
        const double                  a = 0 , 
        const double                  b = 1 ) 
        : F2 ( f1 , f2 ) 
        , m_a   ( std::min ( a , b ) ) 
        , m_b   ( std::max ( a , b ) ) 
      {}
      /// constructor  from two functions 
      template <class FUNCTION1, 
                class FUNCTION2>
      SmoothTransition
      ( FUNCTION1 f1       ,
        FUNCTION2 f2       ,
        const double a = 0 , 
        const double b = 1 ) 
        : F2 ( f1 , f2 ) 
        , m_a   ( std::min ( a , b ) ) 
        , m_b   ( std::max ( a , b ) ) 
      {}
      /// copy constructor 
      SmoothTransition ( const SmoothTransition&  ) = default ;
      /// move constructor
      SmoothTransition (       SmoothTransition&& ) = default ;
      // ======================================================================
    public:
      // ======================================================================
      template <class FUNCTION1,
                class FUNCTION2>
      static inline SmoothTransition
      create
      ( FUNCTION1    f1    ,
        FUNCTION2    f2    ,
        const double a = 0 , 
        const double b = 1 )
      { return SmoothTransition ( f1 , f2 , a , b ) ; }
      // ======================================================================      
    public: // transition function 
      // ======================================================================
      inline double transition 
      ( const double x ) const  
      {
        if      ( x <= m_a ) { return 0 ; }
        else if ( x >= m_b ) { return 1 ; }
        return phi ( ( x - m_a ) / ( m_b - m_a ) ) ;
      }          
      // ======================================================================
    private:
      // ======================================================================
      inline double psi 
      ( const double t )  const  
      { return t <= 0 ? 0.0 : std::exp ( -1/t ) ;}
      // ======================================================================
      inline double phi 
      ( const double t )  const  
      {
        const double p1 = psi (     t ) ;
        const double p2 = psi ( 1 - t ) ;
        return p1 / ( p1 + p2 ) ;
      }
      // ======================================================================      
    public:
      // ======================================================================
      /// the only one important  method 
      inline double operator() 
      ( const double x ) const
      {
        if       ( x <= m_a ) { return m_fun1 ( x ) ; }
        else if  ( x >= m_b ) { return m_fun2 ( x ) ; }
        //
        const double p = transition ( x ) ;
        return p * m_fun2 ( x ) + ( 1 - p ) * m_fun1 ( x ) ;
      }
      // ======================================================================        
    private :
      // ======================================================================     
      /// low x  
      double m_a { 0 } ; // low x 
      /// high x  
      double m_b { 1 } ; // high x 
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Apply2
     *  keep and apply arbitrary 2D-function
     */
    class Apply2
    {
    public:
      // ======================================================================
      /// constructor  from the function 
      Apply2
      ( std::function<double(double,double)> f )
        : m_fun ( f ) 
      {}
      // ======================================================================
      /// constrtuctor  from the function 
      template <class FUNCTION>
      Apply2 ( FUNCTION f )
        : m_fun ( f  )
      {}
      /// copy constructor 
      Apply2 ( const Apply2&  ) = default ;
      /// move constructor
      Apply2 (       Apply2&& ) = default ;
      // ======================================================================
    public:
      // ======================================================================
      template <class FUNCTION>
      static inline Apply2 
      create ( FUNCTION f ) { return Apply2 ( f ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// the only one important method 
      inline double operator() 
      ( const double x , 
        const double y ) const 
      { return m_fun ( x , y ) ; }
      // ======================================================================        
    protected:
      // ======================================================================      
      /// the function itself
      std::function<double(double,double)> m_fun ; // the function itself
      // ======================================================================      
    } ;
    // ========================================================================
    /** @class Apply3
     *  keep and apply arbitrary 3D-function
     */
    class Apply3
    {
    public:
      // ======================================================================
      /// constructor  from the function 
      Apply3
      ( std::function<double(double,double,double)> f )
        : m_fun ( f ) 
      {}
      // ======================================================================
      /// constrtuctor  from the function 
      template <class FUNCTION>
      Apply3 ( FUNCTION f )
        : m_fun ( f  )
      {}
      /// copy constructor 
      Apply3 ( const Apply3&  ) = default ;
      /// move constructor
      Apply3 (       Apply3&& ) = default ;
      // ======================================================================
    public:
      // ======================================================================
      template <class FUNCTION>
      static inline Apply3 
      create ( FUNCTION f ) { return Apply3 ( f ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// the only one important method 
      inline double operator() 
      ( const double x , 
        const double y , 
        const double z ) const 
      { return m_fun ( x , y , z ) ; }
      // ======================================================================        
    protected:
      // ======================================================================      
      /// the function itself
      std::function<double(double,double,double)> m_fun ; // the function itself
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
