// ============================================================================
#ifndef OSTAP_DALITZINTEGRAL_H 
#define OSTAP_DALITZINTEGRAL_H 1
// ============================================================================
// Include files 
// ============================================================================
//  STD&STL
// ============================================================================
#include  <functional>
// ============================================================================
// forward declaratinos 
// ============================================================================
namespace  Ostap
{
  // ==========================================================================
  namespace  Kinematics
  {
    // ========================================================================
    class Dalitz0 ;
    class Dalitz  ;
    // ========================================================================
  }
  // ==========================================================================
  namespace Math
  {
    class WorkSpace ;
  }
  // ==========================================================================
}
// ============================================================================
namespace Ostap
{
  // =========================================================================
  namespace Math
  {
    // =======================================================================
    /** @class DalitzIntegrator DalitzIntegrtor.h Ostap/DalitzIntegrator.h
     *  Helper class to make the integration over Dalitz plane easy 
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2019-11-02
     */
    class DalitzIntegrator 
    {
    public: 
      // ======================================================================
      /// the actual type of 1D-function to integrate 
      typedef std::function<double(double)>               function1 ;
      /// the actual type of 2D-function to integrate 
      typedef std::function<double(double,double)>        function2 ;
      /// the actual type of 3D-function to integrate 
      typedef std::function<double(double,double,double)> function3 ;
      // ======================================================================
    public : // 1D integrations  (no workspace) 
      // ======================================================================
      /** evaluate integral over \f$s\f$ for \f$ f(s,s_1,s_2) \f$
       *  \f[ F(s_1,s_2)  = \int_{s_{min}}^{s_{max}} ds f(s,s_1,s_2) \f]
       *  @param f3  the fun§ãtion \f$  f(s,s_1,s_2)\f$
       *  @param s1    value of \f$ s_1\f$
       *  @param s2    value of \f$ s_2\f$
       *  @param smax  upper inntegration limit for  \f$s\f$
       *  @param d     helper Dalitz-object 
       */
      static double integrate_s
      ( function3                         f3   ,
        const double                      s1   ,
        const double                      s2   ,
        const double                      smax , 
        const Ostap::Kinematics::Dalitz0& d    ) ;
      // ======================================================================
      /** evaluate integral over \f$s_1\f$ for \f$ f(s,s_1,s_2) \f$
       *  \f[ F(s,s_2)  = \int  ds_1 f(s,s_1,s_2) \f]
       *  @param f3  the function \f$  f(s,s_1,s_2)\f$
       *  @param s     value of \f$ s\f$
       *  @param s2    value of \f$ s_2\f$
       *  @param d     helper Dalitz-object 
       */
      static double integrate_s1
      ( function3                         f3 ,
        const double                      s  ,
        const double                      s2 ,
        const Ostap::Kinematics::Dalitz0& d  ) ;
      // ======================================================================
      /** evaluate integral over \f$s_1\f$ for \f$ f(s_1,s_2) \f$
       *  \f[ F(s_2)  = \int  ds_1 f(s_1,s_2) \f]
       *  @param f2  the function \f$  f(s_1,s_2)\f$
       *  @param s2    value of \f$ s_2\f$
       *  @param d     helper Dalitz-object 
       */
      static double integrate_s1
      ( function2                         f2 ,
        const double                      s2 ,
        const Ostap::Kinematics::Dalitz&  d  ) ;
      // ======================================================================
      /** evaluate integral over \f$s_1\f$ for \f$ f(s_1,s_2) \f$
       *  \f[ F(s_2)  = \int  ds_1 f(s_1,s_2) \f]
       *  @param f2    the function \f$  f(s_1,s_2)\f$
       *  @param s     value of \f$ s \f$
       *  @param s2    value of \f$ s_2\f$
       *  @param d     helper Dalitz-object 
       */
      static double integrate_s1
      ( function2                         f2 ,
        const double                      s  , 
        const double                      s2 ,
        const Ostap::Kinematics::Dalitz0& d  ) ;
      // ==========================a============================================
    public : // 1D integrations  (with workspace) 
      // ======================================================================
      /** evaluate integral over \f$s\f$ for \f$ f(s,s_1,s_2) \f$
       *  \f[ F(s_1,s_2)  = \int_{s_{min}}^{s_{max}} ds f(s,s_1,s_2) \f]
       *  @param f3  the fun§ãtion \f$  f(s,s_1,s_2)\f$
       *  @param s1    value of \f$ s_1\f$
       *  @param s2    value of \f$ s_2\f$
       *  @param smax  upper inntegration limit for  \f$s\f$
       *  @param d     helper Dalitz-object 
       *  @param ws    integration workspace  
       */
      static double integrate_s
      ( function3                         f3   ,
        const double                      s1   ,
        const double                      s2   ,
        const double                      smax , 
        const Ostap::Kinematics::Dalitz0& d    ,
        const Ostap::Math::WorkSpace&     ws   ) ;
      // ======================================================================
      /** evaluate integral over \f$s_1\f$ for \f$ f(s,s_1,s_2) \f$
       *  \f[ F(s,s_2)  = \int  ds_1 f(s,s_1,s_2) \f]
       *  @param f3  the function \f$  f(s,s_1,s_2)\f$
       *  @param s     value of \f$ s\f$
       *  @param s2    value of \f$ s_2\f$
       *  @param d     helper Dalitz-object 
       *  @param ws    integration workspace  
       */
      static double integrate_s1
      ( function3                         f3 ,
        const double                      s  ,
        const double                      s2 ,
        const Ostap::Kinematics::Dalitz0& d  ,
        const Ostap::Math::WorkSpace&     ws ) ;
      // ======================================================================
      /** evaluate integral over \f$s_1\f$ for \f$ f(s_1,s_2) \f$
       *  \f[ F(s_2)  = \int  ds_1 f(s_1,s_2) \f]
       *  @param f2  the function \f$  f(s_1,s_2)\f$
       *  @param s2    value of \f$ s_2\f$
       *  @param d     helper Dalitz-object 
       *  @param ws    integration workspace  
       */
      static double integrate_s1
      ( function2                         f2 ,
        const double                      s2 ,
        const Ostap::Kinematics::Dalitz&  d  ,
        const Ostap::Math::WorkSpace&     ws ) ;
      // ======================================================================
      /** evaluate integral over \f$s_1\f$ for \f$ f(s_1,s_2) \f$
       *  \f[ F(s_2)  = \int  ds_1 f(s_1,s_2) \f]
       *  @param f2    the function \f$  f(s_1,s_2)\f$
       *  @param s     value of \f$ s \f$
       *  @param s2    value of \f$ s_2\f$
       *  @param d     helper Dalitz-object 
       *  @param ws    integration workspace  
       */
      static double integrate_s1
      ( function2                         f2 ,
        const double                      s  , 
        const double                      s2 ,
        const Ostap::Kinematics::Dalitz0& d  ,
        const Ostap::Math::WorkSpace&     ws ) ;
      // ==========================a============================================
    public : // 1D integrations with cache
      // ======================================================================
      /** evaluate integral over \f$s\f$ for \f$ f(s,s_1,s_2) \f$
       *  \f[ F(s_1,s_2)  = \int_{s_{min}}^{s_{max}} ds f(s,s_1,s_2) \f]
       *  @param tag tag that indicate theuniquness of function
       *  @param f3  the fun§ãtion \f$  f(s,s_1,s_2)\f$
       *  @param s1    value of \f$ s_1\f$
       *  @param s2    value of \f$ s_2\f$
       *  @param smax  upper inntegration limit for  \f$s\f$
       *  @param d     helper Dalitz-object 
       *  @param ws    integration workspace  
       */
      static double integrate_s
      ( const std::size_t                 tag ,
        function3                         f3   ,
        const double                      s1   ,
        const double                      s2   ,
        const double                      smax , 
        const Ostap::Kinematics::Dalitz0& d    ,
        const Ostap::Math::WorkSpace&     ws   ) ;
      // ======================================================================
      /** evaluate integral over \f$s_§í\f$ for \f$ f(s,s_1,s_2) \f$
       *  \f[ F(s,s_2)  = \int  ds_1 f(s,s_1,s_2) \f]
       *  @param tag tag that indicate theuniquness of function
       *  @param f3  the fun§ãtion \f$  f(s,s_1,s_2)\f$
       *  @param s     value of \f$ s\f$
       *  @param s2    value of \f$ s_2\f$
       *  @param d     helper Dalitz-object 
       *  @param ws    integration workspace  
       */
      static double integrate_s1
      ( const std::size_t                 tag ,
        function3                         f3 ,
        const double                      s  ,
        const double                      s2 ,
        const Ostap::Kinematics::Dalitz0& d  ,
        const Ostap::Math::WorkSpace&     ws ) ;
      // ======================================================================
      /** evaluate integral over \f$s_§í\f$ for \f$ f(s_1,s_2) \f$
       *  \f[ F(s_2)  = \int  ds_1 f(s_1,s_2) \f]
       *  @param tag tag that indicate theuniquness of function
       *  @param f2  the fun§ãtion \f$  f(s_1,s_2)\f$
       *  @param s2    value of \f$ s_2\f$
       *  @param d     helper Dalitz-object 
       *  @param ws    integration workspace  
       */
      static double integrate_s1
      ( const std::size_t                 tag ,
        function2                         f2 ,
        const double                      s2 ,
        const Ostap::Kinematics::Dalitz&  d  ,
        const Ostap::Math::WorkSpace&     ws ) ;
      // ======================================================================
      /** evaluate integral over \f$s_§í\f$ for \f$ f(s_1,s_2) \f$
       *  \f[ F(s_2)  = \int  ds_1 f(s_1,s_2) \f]
       *  @param tag tag that indicate the uniquness of function
       *  @param f2    the fun§ãtion \f$  f(s_1,s_2)\f$
       *  @param s     value of \f$ s \f$
       *  @param s2    value of \f$ s_2\f$
       *  @param d     helper Dalitz-object 
       *  @param ws    integration workspace  
       */
      static double integrate_s1
      ( const std::size_t                 tag ,
        function2                         f2 ,
        const double                      s  , 
        const double                      s2 ,
        const Ostap::Kinematics::Dalitz0& d  ,
        const Ostap::Math::WorkSpace&     ws ) ;
      // ======================================================================
    public : // 2D integrations
      // ======================================================================
      /** evaluate the integral over \f$s_1\f$ , \f$s_2\f$ variables 
       *  \f[ \int\int ds_1 ds_2 f(s, s_1,s_2) = 
       *  \int_{s_1^{min}}^{s_1^{max}} ds_1 
       *  \int_{s_2^{min}(s_1)}^{s_2^{max}(s_1)} ds_2 f(s,s_1,s_2)  \f] 
       *  @param f3 the function \f$ f(s, s_1,s_2)\f$ 
       *  @param d  helper Dalitz-object 
       *  @return integral over Dalitz plot
       */
      static double integrate_s1s2
      ( function3                        f3 ,
        const Ostap::Kinematics::Dalitz& d  ) ;
      // ======================================================================
      /** evaluate the integral over \f$s_1\f$ , \f$s_2\f$ variables 
       *  \f[ \int\int ds_1 ds_2 f(s, s_1,s_2) = 
       *  \int_{s_1^{min}}^{s_1^{max}} ds_1 
       *  \int_{s_2^{min}(s_1)}^{s_2^{max}(s_1)} ds_2 f(s,s_1,s_2)  \f] 
       *  @param f3 the function \f$ f(s, s_1,s_2)\f$ 
       *  @param d  helper Dalitz-object 
       *  @return integral over Dalitz plot
       */
      static double integrate_s1s2
      ( function3                        f3 ,
        const double                      s ,
        const Ostap::Kinematics::Dalitz0& d ) ;
      // ======================================================================
      /** evaluate the integral over \f$s_1\f$ , \f$s_2\f$ variables 
       *  \f[ \int\int ds_1 ds_2 f( s_1,s_2) = 
       *  \int_{s_1^{min}}^{s_1^{max}} ds_1 
       *  \int_{s_2^{min}(s_1)}^{s_2^{max}(s_1)} ds_2 f(s,s_1,s_2)  \f] 
       *  @param f2 the function \f$ f(s_1,s_2)\f$ 
       *  @param s  the \f$s=M^2\f$
       *  @param d  helper Dalitz-object 
       *  @return integral over Dalitz plot
       */
      static double integrate_s1s2
      ( function2                         f2 ,
        const double                      s  ,
        const Ostap::Kinematics::Dalitz0& d  ) ;
      // =======================================================================      
      /** evaluate the integral over \f$s\f$ , \f$s_1\f$ variables 
       *  \f[ \int\int f( s, s_1,s_2) ds ds_1 \f]  
       *  @param f3 the function \f$  f(s,s_1,s_2) \f$
       *  @param s2 fixed value of s2 
       *  @param smax upper-edge fo inntegrato ove \f$ s \f$
       *  @param d  helper Dalitz-object 
       *  @return integral over \f$ s, s_1\f$
       */
      static double integrate_ss1
      ( function3                         f3   ,
        const double                      s2   ,
        const double                      smax ,
        const Ostap::Kinematics::Dalitz0& d    ) ;
      // =======================================================================
    public: // integrations with cache
      // =======================================================================
      /** evaluate the integral over \f$s_1\f$ , \f$s_2\f$ variables 
       *  \f[ \int\int ds_1 ds_2 f(s, s_1,s_2) = 
       *  \int_{s_1^{min}}^{s_1^{max}} ds_1 
       *  \int_{s_2^{min}(s_1)}^{s_2^{max}(s_1)} ds_2 f(s,s_1,s_2)  \f] 
       *  @param tag tag that indicate the uniquness of function
       *  @param f3 the function \f$ f(s, s_1,s_2)\f$ 
       *  @param d  helper Dalitz-object 
       *  @return integral over Dalitz plot
       */
      static double integrate_s1s2
      ( const std::size_t                tag ,
        function3                        f3  ,
        const Ostap::Kinematics::Dalitz& d   ) ;
      // ======================================================================
      /** evaluate the integral over \f$s_1\f$ , \f$s_2\f$ variables 
       *  \f[ \int\int ds_1 ds_2 f(s, s_1,s_2) = 
       *  \int_{s_1^{min}}^{s_1^{max}} ds_1 
       *  \int_{s_2^{min}(s_1)}^{s_2^{max}(s_1)} ds_2 f(s,s_1,s_2)  \f] 
       *  @param tag tag that indicate the uniquness of function
       *  @param f3 the function \f$ f(s, s_1,s_2)\f$ 
       *  @param d  helper Dalitz-object 
       *  @return integral over Dalitz plot
       */
      static double integrate_s1s2
      ( const std::size_t                tag ,
        function3                        f3 ,
        const double                      s ,
        const Ostap::Kinematics::Dalitz0& d ) ;
      // ======================================================================
      /** evaluate the integral over \f$s_1\f$ , \f$s_2\f$ variables 
       *  \f[ \int\int ds_1 ds_2 f( s_1,s_2) = 
       *  \int_{s_1^{min}}^{s_1^{max}} ds_1 
       *  \int_{s_2^{min}(s_1)}^{s_2^{max}(s_1)} ds_2 f(s,s_1,s_2)  \f] 
       *  @param tag tag that indicate the uniquness of function
       *  @param f2 the function \f$ f(s_1,s_2)\f$ 
       *  @param s  the \f$s=M^2\f$
       *  @param d  helper Dalitz-object 
       *  @return integral over Dalitz plot
       */
      static double integrate_s1s2
      ( const std::size_t                tag ,
        function2                         f2 ,
        const double                      s  ,
        const Ostap::Kinematics::Dalitz0& d  ) ;
      // =======================================================================
      /** evaluate the integral over \f$s\f$ , \f$s_1\f$ variables 
       *  \f[ \int\int f( s, s_1,s_2) ds ds_1 \f]  
       *  @param tag tag that indicate the uniquness of function
       *  @param f3 the function \f$  f(s,s_1,s_2) \f$
       *  @param s2 fixed value of s2 
       *  @param smax upper-edge fo inntegrato ove \f$ s \f$
       *  @param d  helper Dalitz-object 
       *  @return integral over \f$ s, s_1\f$
       */
      static double integrate_ss1
      ( const std::size_t                tag ,
        function3                         f3   ,
        const double                      s2   ,
        const double                      smax ,
        const Ostap::Kinematics::Dalitz0& d    ) ;
      // =======================================================================
    public:
      // =======================================================================
      /** integrate the function \f$ f(s,s_1,s_2)\f$ over \f$ s \f$
       *  \f[ F(s_1,s_2)  = \int_{s_{min}}^{s_{max}} ds f(s,s_1,s_2) \f]
       *  @see Ostap::Math::DalitzIntegrator::integrate_s
       */
      template <class FUNCTION, typename... ARGS>
      static inline double 
      integral_s ( FUNCTION        f    , 
                   const double    s1   , 
                   const double    s2   , 
                   const double    smax , 
                   const ARGS& ... args )
      { return integrate_s ( std::cref ( f ) , s1 , s2  , smax , args... ) ; }
      // =====================================================================      
      /** integrate the function \f$ f(s,s_1,s_2)\f$ or \f$ f(s_1,s_2) 
       *  over \f$ \f$ s_1 \f$
       *  @see Ostap::Math::DalitzIntegrator::integrate_s1
       */
      template <class FUNCTION, typename... ARGS>
      static inline double 
      integral_s1 ( FUNCTION        f    , 
                    const double    x    , 
                    const ARGS& ... args )
      { return integrate_s1 ( std::cref ( f ) , x , args... ) ; }
      // =====================================================================
      /** integrate the function \f$ f(s,s_1,s_2)\f$ or \f$ f(s_1,s_2) 
       *  over \f$ s\f$  and \f$s_1\f$
       *  @see Ostap::Math::DalitzIntegrator::integrate_ss1
       */
      template <class FUNCTION, typename... ARGS>
      static inline double 
      integral_ss1 ( FUNCTION        f    , 
                     const ARGS& ... args )
      { return integrate_ss1 ( std::cref ( f ) , args... ) ; }
      // =====================================================================
    public:
      // =====================================================================
      /** integrate the function \f$ f(s,s_1,s_2)\f$ over \f$ s \f$
       *  \f[ F(s_1,s_2)  = \int_{s_{min}}^{s_{max}} ds f(s,s_1,s_2) \f]
       *  @see Ostap::Math::DalitzIntegrator::integrate_s
       */
      template <class FUNCTION, typename... ARGS>
      static inline double 
      integral_s ( const std::size_t tag  , 
                   FUNCTION          f    , 
                   const double      s1   , 
                   const double      s2   , 
                   const double      smax , 
                   const ARGS& ...   args )
      { return integrate_s ( tag , std::cref ( f ) , s1 , s2  , smax , args... ) ; }
      // =======================================================================
      /** integrate the function \f$ f(s,s_1,s_2)\f$ or \f$ f(s_1,s_2) 
       *  over \f$ \f$ s_1 \f$
       *  @see Ostap::Math::DalitzIntegrator::integrate_s1
       */
      template <class FUNCTION, typename... ARGS>
      static inline double 
      integral_s1 ( const std::size_t tag  ,
                    FUNCTION          f    , 
                    const double      x    , 
                    const ARGS& ...   args )
      { return integrate_s1 ( tag , std::cref ( f ) , x , args... ) ; }
      // ======================================================================
      /** integrate the function \f$ f(s,s_1,s_2)\f$ or \f$ f(s_1,s_2) 
       *  over \f$ s\f$  and \f$s_1\f$
       *  @see Ostap::Math::DalitzIntegrator::integrate_ss1
       */
      template <class FUNCTION, typename... ARGS>
      static inline double 
      integral_ss1 ( const std::size_t tag  , 
                     FUNCTION          f    , 
                     const ARGS& ...   args )
      { return integrate_ss1 ( tag , std::cref ( f ) , args... ) ; }
      // ======================================================================
    public: // as functions of energy 
      // ======================================================================
      /** evaluate the integral over \f$e_2\f$ , \f$e_3\f$ variables 
       *  \f[ \int\int de_2 de_3 f ( M ,e_2,e_3) = 
       *  \int_{e_2^{min}}^{e_2^{max}} de_2 
       *  \int_{e_3^{min}(e_2)}^{e_3^{max}(e_2)} de_3 f(M,e_2,e_3)  \f] 
       *  @param  f the function \f$ f(M,e_2,e_3)\f$ 
       *  @param d  helper Dalitz-object 
       *  @return the integral over Dalitz plot
       */
      static double integrate_e2e3
      ( function3                        f3 ,
        const Ostap::Kinematics::Dalitz& d  ) ;
      // =======================================================================
      /** evaluate the integral over \f$e_2\f$ , \f$e_3\f$ variables 
       *  \f[ \int\int de_2 de_3 f(e_2,e_3) = 
       *  \int_{e_2^{min}}^{e_2^{max}} de_2 
       *  \int_{e_3^{min}(e_2)}^{e_3^{max}(e_2)} de_3 f(e_2,e_3)  \f] 
       *  @param  f the function \f$ f(e_2,e_3)\f$ 
       *  @param d  helper Dalitz-object 
       *  @return the integral over dalitz plot
       */
      static double integrate_e2e3
      ( function2                        f2 ,
        const Ostap::Kinematics::Dalitz& d  ) ;
      // ======================================================================
    };
    // ========================================================================
  } //                                         The end of namespace Ostap::Math 
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_DALITZINTEGRAL_H
// ============================================================================
