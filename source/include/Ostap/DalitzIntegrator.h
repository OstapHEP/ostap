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
// Ostap
// ============================================================================
#include "Ostap/Dalitz.h"
#include "Ostap/Workspace.h"
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
    public:
      // ======================================================================
      /** constructor from Dalitz configurtaion and the integration workspace
       *  @param dalitz  Dalitz configuration 
       *  @param size    size of the integrtaion woekjspace 
       *  @see Ostap::Math::WorkSpace
       *  @see Ostap::Kinematics::Dalitz0 
       */
      DalitzIntegrator
      ( const Ostap::Kinematics::Dalitz0& dalitz   , 
        const std::size_t                 size = 0 ) ;
      // ======================================================================
      /// get the Dalitz configuration 
      const Ostap::Kinematics::Dalitz0& dalitz () const { return m_dalitz ; }
      // ======================================================================
    public: // 1D integrals 
      // ======================================================================
      /** evaluate integral over \f$s\f$ for \f$ f(s,s_1,s_2) \f$
       *  \f[ F(s_1,s_2)  = \int_{s_{min}}^{s_{max}} ds f(s,s_1,s_2) \f]
       *  @param f3  the fun§ãtion \f$  f(s,s_1,s_2)\f$
       *  @param s1    value of \f$ s_1\f$
       *  @param s2    value of \f$ s_2\f$
       *  @param smax  upper inntegration limit for  \f$s\f$
       */      
      template <class FUNCTION3, 
                typename = std::enable_if<std::is_convertible<FUNCTION3,function3>::value> >      
      double integrate_s 
      ( FUNCTION3         f3      , 
        const double      s1      ,  
        const double      s2      ,
        const double      smax    , 
        const std::size_t tag = 0 ) const 
      { return integrate_s ( std::cref ( f3 ) , s2 , s2 , smax , m_dalitz , m_workspace , tag ) ; }
      // ======================================================================
      /** evaluate integral over \f$s_1\f$ for \f$ f(s,s_1,s_2) \f$ of \f$ f(s_1,s_2) \f$
       *  \f[ F(s,s_2)  = \int  ds_1 f(s,s_1,s_2) \f] or 
       *  \f[ F(s,s_2)  = \int  ds_1 f(s_1,s_2) \f] 
       *  @param f3  the function \f$  f(s,s_1,s_2)\f$ or \f$ f(s_1,s_2)\f$
       *  @param s     value of \f$ s\f$
       *  @param s2    value of \f$ s_2\f$
       */
      template <class FUNCTION23,
                typename = std::enable_if<std::is_convertible<FUNCTION23,function2>::value ||
                                          std::is_convertible<FUNCTION23,function3>::value > > 
        double integrate_s1
        ( FUNCTION23        f23     ,
          const double      s       ,
          const double      s2      , 
          const std::size_t tag = 0 ) const 
      { return integrate_s1 ( std::cref ( f23 ) , s , s2 ,  m_dalitz , m_workspace , tag ) ; }
      // ======================================================================
      // integrate over s2  
      // ======================================================================
      template <class FUNCTION2,
                typename = std::enable_if<std::is_convertible<FUNCTION2,function2>::value> >
      double integrate_s2 
      ( FUNCTION2         f2      , 
        const double      s       ,
        const double      s1      , 
        const std::size_t tag = 0 ) const 
      {
        // swap arguments   
        auto fc = std::cref ( f2 ) ;
        auto ff = [fc]( const double s_1 , const double s_2 )-> double
          { return fc ( s_2 , s_1 ) ; } ;
        return integrate_s1 ( std::cref ( ff) , s , s1 ,  m_dalitz321 , m_workspace , tag ) ; 
      }
      // ======================================================================
      // integrate over s2  
      // ======================================================================
      template <class FUNCTION3, 
                typename = std::enable_if<std::is_convertible<FUNCTION3,function3>::value> >      
      double integrate3_s2 
      ( FUNCTION3         f3      ,
        const double      s       ,
        const double      s1      , 
        const std::size_t tag = 0 ) const 
      {
        // swap arguments   
        auto fc = std::cref ( f3 ) ;
        auto ff = [fc]( const double s , const double s_1 , const double s_2 )-> double
          { return fc ( s , s_2 , s_1 ) ; } ;
        return integrate_s1 ( std::cref ( ff ) , s , s1 ,  m_dalitz321 , m_workspace , tag ) ; 
      }
      // ======================================================================
    public:  // 2D integrals 
      // ======================================================================
      /** evaluate the integral over \f$s_1\f$ , \f$s_2\f$ variables 
       *  \f[ \int\int ds_1 ds_2 f(s, s_1,s_2) = 
       *  \int_{s_1^{min}}^{s_1^{max}} ds_1 
       *  \int_{s_2^{min}(s_1)}^{s_2^{max}(s_1)} ds_2 f(s,s_1,s_2)  \f] 
       *  @param f3 the function \f$ f(s, s_1,s_2)\f$ 
       *  @return integral over Dalitz plot
       */
      template <class FUNCTION23,
                typename = std::enable_if<std::is_convertible<FUNCTION23,function2>::value ||
                                          std::is_convertible<FUNCTION23,function3>::value > >
        double integrate_s1s2
        ( FUNCTION23           f23     ,
          const double         s       , 
          const std::size_t    tag = 0 , 
          const unsigned short n1  = 0 , 
          const unsigned short n2  = 0 ) const 
      { return integrate_s1s2 ( std::cref ( f23 ) , s , m_dalitz , tag , n1 , n2 ) ; }
      // ======================================================================
      /** evaluate the integral over \f$s\f$ , \f$s_1\f$ variables 
       *  \f[ \int\int f( s, s_1,s_2) ds ds_1 \f]  
       *  @param f3 the function \f$  f(s,s_1,s_2) \f$
       *  @param s2 fixed value of s2 
       *  @param smax upper-edge for integration over \f$ s \f$
       *  @param d  helper Dalitz-object 
       *  @return integral over \f$ s, s_1\f$
       */
      template <class FUNCTION3 ,
                typename = std::enable_if<std::is_convertible<FUNCTION3,function3>::value> >      
      double integrate_ss1
      ( FUNCTION3            f3      ,
        const double         s2      ,
        const double         smax    , 
        const std::size_t    tag = 0 ,
        const unsigned short n1  = 0 , 
        const unsigned short n2  = 0 ) const 
      { return integrate_ss1 ( std::cref ( f3 ) , s2 , smax , m_dalitz , tag , n1 , n2 ) ; }
      // ===================================================================== 
      /** evaluate the integral over \f$s\f$ , \f$s_1\f$ variables 
       *  \f[ \int\int f( s, s_1,s_2) ds ds_1 \f]  
       *  @param f3 the function \f$  f(s,s_1,s_2) \f$
       *  @param s2 fixed value of s2 
       *  @param smin lower-edge for integration over \f$ s \f$
       *  @param smax upper-edge for integration over \f$ s \f$
       *  @param d  helper Dalitz-object 
       *  @return integral over \f$ s, s_1\f$
       */
      template <class FUNCTION3 ,
                typename = std::enable_if<std::is_convertible<FUNCTION3,function3>::value> >
      double integrate_ss1
      ( FUNCTION3            f3      ,
        const double         s2      ,
        const double         smin    ,
        const double         smax    ,
        const std::size_t    tag = 0 , 
        const unsigned short n1  = 0 , 
        const unsigned short n2  = 0 ) const 
      { return integrate_ss1 ( std::cref ( f3 ) , s2 , smin , smax , m_dalitz , tag , n1 , n2 ) ; }
      // ====================================================================== 
      /** evaluate the integral over \f$s\f$ , \f$s_2\f$ variables 
       *  \f[ \int\int f( s, s_1,s_2) ds ds_s \f]  
       *  @param f3 the function \f$  f(s,s_1,s_2) \f$
       *  @param s1 fixed value of s1 
       *  @param smax upper-edge for integration over \f$ s \f$
       *  @param d  helper Dalitz-object 
       *  @return integral over \f$ s, s_1\f$
       */
      template <class FUNCTION3 ,
                typename = std::enable_if<std::is_convertible<FUNCTION3,function3>::value> >
      double integrate_ss2
      ( FUNCTION3            f3      ,
        const double         s1      ,
        const double         smax    ,
        const std::size_t    tag = 0 ,
        const unsigned short n1  = 0 , 
        const unsigned short n2  = 0 ) const 
      {
        // swap arguments   
        auto fc = std::cref ( f3 ) ;
        auto ff = [fc]( const double s , const double s_1 , const double s_2 )-> double
                  { return fc ( s , s_2 , s_1 ) ; } ;
        // invoke integration over s,s1 with swapped arguments 
        return integrate_ss1 ( std::cref ( ff ) , s1 , smax , m_dalitz321 , tag , n1 , n2 ) ;
      }
      // ===================================================================== 
      /** evaluate the integral over \f$s\f$ , \f$s_s\f$ variables 
       *  \f[ \int\int f( s, s_1,s_2) ds ds_2 \f]  
       *  @param f3 the function \f$  f(s,s_1,s_2) \f$
       *  @param s1 fixed value of s1 
       *  @param smin lowe-edge for integration over \f$ s \f$
       *  @param smax upper-edge for integration over \f$ s \f$
       *  @param d  helper Dalitz-object 
       *  @return integral over \f$ s, s_1\f$
       */
      template <class FUNCTION3 ,
                typename = std::enable_if<std::is_convertible<FUNCTION3,function3>::value> >
      double integrate_ss2
      ( FUNCTION3            f3      ,
        const double         s1      ,
        const double         smin    ,
        const double         smax    ,
        const std::size_t    tag = 0 , 
        const unsigned short n1  = 0 , 
        const unsigned short n2  = 0 ) const 
      {
        // swap arguments   
        auto fc = std::cref ( f3 ) ;
        auto ff = [fc]( const double s , const double s_1 , const double s_2 )-> double
          { return fc ( s , s_2 , s_1 ) ; } ;
        // invoke integration over s,s1 with swapped arguments 
        return integrate_ss1 ( std::cref ( ff ) , s1 , smin , smax , m_dalitz321 , tag , n1 , n2 ) ;
      }
      // ======================================================================
    public : // 1D integrations  (with workspace) 
      // ======================================================================
      /** evaluate integral over \f$s\f$ for \f$ f(s,s_1,s_2) \f$
       *  \f[ F(s_1,s_2)  = \int_{s_{min}}^{s_{max}} ds f(s,s_1,s_2) \f]
       *  @param f3  the fun§ãtion \f$  f(s,s_1,s_2)\f$
       *  @param s1    value of \f$ s_1\f$
       *  @param s2    value of \f$ s_2\f$
       *  @param smax  upper integration limit for  \f$s\f$
       *  @param d     helper Dalitz-object 
       *  @param ws    integration workspace  
       */
      static double integrate_s
      ( function3                         f3        ,
        const double                      s1        ,
        const double                      s2        ,
        const double                      smax      , 
        const Ostap::Kinematics::Dalitz0& d         ,
        const Ostap::Math::WorkSpace&     ws        , 
        const std::size_t                 tag = 0 ) ;
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
      ( function3                         f3      ,
        const double                      s       ,
        const double                      s2      ,
        const Ostap::Kinematics::Dalitz0& d       ,
        const Ostap::Math::WorkSpace&     ws      ,
        const std::size_t                 tag = 0 ) ;
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
      ( function2                         f2      ,
        const double                      s       , 
        const double                      s2      ,
        const Ostap::Kinematics::Dalitz0& d       ,
        const Ostap::Math::WorkSpace&     ws      ,
        const std::size_t                 tag = 0 ) ;
      // ==========================a============================================
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
      ( function3                         f3      ,
        const double                      s       ,
        const Ostap::Kinematics::Dalitz0& d       ,
        const std::size_t                 tag = 0 ,
        const unsigned short              nx  = 0 , 
        const unsigned short              ny  = 0 ) ;
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
      ( function2                         f2      ,
        const double                      s       ,
        const Ostap::Kinematics::Dalitz0& d       ,
        const std::size_t                 tag = 0 ,
        const unsigned short              nx  = 0 , 
        const unsigned short              ny  = 0 ) ;
      // =======================================================================      
      /** evaluate the integral over \f$s\f$ , \f$s_1\f$ variables 
       *  \f[ \int\int f( s, s_1,s_2) ds ds_1 \f]  
       *  @param f3 the function \f$  f(s,s_1,s_2) \f$
       *  @param s2 fixed value of s2 
       *  @param smax upper-edge fo inntegrato ove \f$ s \f$
       *  @param d  helper Dalitz-object 
       *  @param tag unique label/tag  (for caching)
       *  @return integral over \f$ s, s_1\f$
       */
      static double integrate_ss1
      ( function3                         f3      ,
        const double                      s2      ,
        const double                      smax    ,
        const Ostap::Kinematics::Dalitz0& d       ,
        const std::size_t                 tag = 0 , 
        const unsigned short              nx  = 0 , 
        const unsigned short              ny  = 0 ) ;
      // =======================================================================
      /** evaluate the integral over \f$s\f$ , \f$s_1\f$ variables 
       *  \f[ \int\int f( s, s_1,s_2) ds ds_1 \f]  
       *  @param f3 the function \f$  f(s,s_1,s_2) \f$
       *  @param s2 fixed value of s2 
       *  @param smin lower-edge for integration over \f$ s \f$
       *  @param smax upper-edge for integration over \f$ s \f$
       *  @param d  helper Dalitz-object 
       *  @return integral over \f$ s, s_1\f$
       */
      static double integrate_ss1
      ( function3                         f3      ,
        const double                      s2      ,
        const double                      smin    ,
        const double                      smax    ,
        const Ostap::Kinematics::Dalitz0& d       ,
        const std::size_t                 tag = 0 ,
        const unsigned short              n1  = 0 , 
        const unsigned short              n2  = 0 ) ;
      // =======================================================================
    public:
      // =======================================================================
      /** integrate the function \f$ f(s,s_1,s_2)\f$ over \f$ s \f$
       *  \f[ F(s_1,s_2)  = \int_{s_{min}}^{s_{max}} ds f(s,s_1,s_2) \f]
       *  @see Ostap::Math::DalitzIntegrator::integrate_s
       */
      template <class FUNCTION3  , 
                typename... ARGS ,
                typename = std::enable_if<std::is_convertible<FUNCTION3,function3>::value> >
      static inline double 
      integral_s
      ( FUNCTION3        f    , 
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
      template <class FUNCTION23 ,  
                typename... ARGS ,
                typename = std::enable_if<std::is_convertible<FUNCTION23,function2>::value ||
                                          std::is_convertible<FUNCTION23,function3>::value > >
        static inline double 
        integral_s1 
        ( FUNCTION23      f    , 
          const double    x    , 
          const ARGS& ... args )
      { return integrate_s1 ( std::cref ( f ) , x , args... ) ; }
      // ======================================================================
     /** integrate the function \f$ f(s,s_1,s_2)\f$ or \f$ f(s_1,s_2) 
       *  over \f$ \f$ s_1,s_2 \f$
       *  @see Ostap::Math::DalitzIntegrator::integrate_s1s2
       */
      template <class FUNCTION23 , 
                typename... ARGS ,
                typename = std::enable_if<std::is_convertible<FUNCTION23,function2>::value ||
                                          std::is_convertible<FUNCTION23,function3>::value > >
        static inline double 
        integral_s1s2 
        ( FUNCTION23      f , 
          const double    s , 
          const ARGS& ... args )
      { return integrate_s1s2 ( std::cref ( f ) , s , args... ) ; }
      // =====================================================================
      /** integrate the function \f$ f(s,s_1,s_2)\f$ or \f$ f(s_1,s_2) 
       *  over \f$ s\f$  and \f$s_1\f$
       *  @see Ostap::Math::DalitzIntegrator::integrate_ss1
       */
      template <class FUNCTION3  , 
                typename... ARGS ,
                typename = std::enable_if<std::is_convertible<FUNCTION3,function3>::value> >                
      static inline double 
      integral_ss1 
      ( FUNCTION3       f    , 
        const ARGS& ... args )
      { return integrate_ss1 ( std::cref ( f ) , args... ) ; }
      // =====================================================================
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
      ( function3                        f3      ,
        const Ostap::Kinematics::Dalitz& d       , 
        const std::size_t                tag = 0 ,
        const unsigned short              n1 = 0 , 
        const unsigned short              n2 = 0 ) ;
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
      ( function2                        f2      ,
        const Ostap::Kinematics::Dalitz& d       , 
        const std::size_t                tag = 0 ,
        const unsigned short              n1 = 0 , 
        const unsigned short              n2 = 0 ) ;
      // ======================================================================
    private:
      // ======================================================================
      /// Dalitz configuration
      Ostap::Kinematics::Dalitz0 m_dalitz    {} ; //  Dalitz configuration
      /// Dalitz configuration (rotated) 3-2-1, \f$ s_1 \leftrghtarrow s_2\f$ 
      Ostap::Kinematics::Dalitz0 m_dalitz321 {} ; //  Dalitz configuration (3,2,1)
      /// Dalitz configuration (rotated) 1-3-2, \f$ s_1 \leftrghtarrow s_3\f$ 
      Ostap::Kinematics::Dalitz0 m_dalitz132 {} ; //  Dalitz configuration (1,3,2)
      /// integration  workspace  
      Ostap::Math::WorkSpace     m_workspace {} ; // integration workspace
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
