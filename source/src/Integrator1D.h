// ============================================================================
#ifndef OSTAP_INTEGRATOR1D_H 
#define OSTAP_INTEGRATOR1D_H 1
// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/GSL_utils.h"
// ============================================================================
// GSL
// ============================================================================
#include "gsl/gsl_errno.h"
#include "gsl/gsl_integration.h"
// ============================================================================
// local
// ============================================================================
#include "GSL_sentry.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    namespace GSL 
    { 
      // ======================================================================
      /** @typedef    Result 
       *  the  type for  result of numerica intehgrtaion routines 
       */
      typedef std::tuple<int,double,double>      Result ;
      // ======================================================================
      /** @class Integrator1D  Integrator1D.h 
       *  Helper class to simplify operations 
       *  with GSL numerical integraion functions 
       *
       *  Typical usage 
       *  @code 
       *  Integrator1D<MYOBJECT> integrator {} ;
       *  ...
       *  auto F = s_i.make_function( this ) ;
       *  int    ierror ;
       *  double result ;
       *  double error  ;
       *  std::tie( ierror, result , error ) = 
       *    integrator.gaq_integrate ( &F    , 
       *                               0 , 1 , // low & high edges 
       *                               workspace ( *this ) ) ;
       *  @endcode
       *  @author Vanya Belyaev
       *  @date   2018-09-21
       */
      template <class FUNCTION>
      class Integrator1D
      {
      public :
        //  ===================================================================
        /// make a function for integration 
        gsl_function make_function ( const FUNCTION* f ) const 
       { 
          gsl_function F ;
          F.params   = const_cast<FUNCTION*>( f )  ;
          F.function = &adapter ;
          return F ; 
        }
        // ====================================================================
      public :
        //  ===================================================================
        /// adaptive integrator 
        Result gaq_integrate   
        ( const gsl_function*        func                 ,       // the function
          const double               xlow                 ,       // low integration edge 
          const double               xhigh                ,       // high integration edge 
          gsl_integration_workspace* workspace            ,       // workspace
          const double               aprecision = 1.e-8   ,       // absolute precision
          const double               rprecision = 1.e-8   ,       // relative precision
          int                        limit      = -1      ,       // limit 
          const char*                reason     = nullptr ,       // message 
          const char*                file       = nullptr ,       // file name 
          const unsigned long        line       = 0       ) const // line number 
        {
          // setup GSL 
          Ostap::Math::GSL::GSL_Error_Handler sentry ;
          //
          double    result =  1.0 ;
          double    error  = -1.0 ;
          if      ( limit <= 0                ) { limit = workspace->limit ; }
          else if ( limit >  workspace->limit ) { limit = workspace->limit ; }
          //
          const int ierror = gsl_integration_qag 
            ( func               ,   // the function
              xlow   , xhigh     ,   // low & high edges
              aprecision         ,   // absolute precision
              rprecision         ,   // relative precision
              limit              ,   // maximum number of subintervals
              GSL_INTEG_GAUSS51  ,   // integration rule
              workspace          ,   // workspace
              &result            ,   // the result
              &error             ) ; // the error in result
          if ( ierror ) { gsl_error ( reason , file , line , ierror ) ; }
          //
          return Result { ierror , result , error } ;  
        }
        // ====================================================================
        /// adaptive integrator 
        Result gaqiu_integrate   
        ( const gsl_function*        func                 ,       // the function
          const double               xlow                 ,       // low integration edge 
          gsl_integration_workspace* workspace            ,       // workspace
          const double               aprecision = 1.e-8   ,       // absolute precision
          const double               rprecision = 1.e-8   ,       // relative precision
          int                        limit      = -1      ,       // limit 
          const char*                reason     = nullptr ,       // message 
          const char*                file       = nullptr ,       // file name 
          const unsigned long        line       = 0       ) const // line number 
        {
          // setup GSL 
          Ostap::Math::GSL::GSL_Error_Handler sentry ;
          //
          double    result =  1.0 ;
          double    error  = -1.0 ;
          if      ( limit <= 0                ) { limit =  workspace->limit() ; }
          else if ( limit >  workspace->limit ) { limit =  workspace->limit() ; }
          //
          const int ierror = gsl_integration_qagiu 
            ( func               ,   // the function
              xlow               ,   // low 
              aprecision         ,   // absolute precision
              rprecision         ,   // relative precision
              limit              ,   // maximum number of subintervals
              workspace          ,   // workspace
              &result            ,   // the result
              &error             ) ; // the error in result
          if ( ierror ) { gsl_error ( reason , file , line , ierror ) ; }
          //
          return Result { ierror , result , error } ;  
        }
        // ====================================================================
        /// adaptive integrator 
        Result gaqil_integrate   
        ( const gsl_function*        func                 ,       // the function
          const double               xhigh                ,       // high integration edge 
          gsl_integration_workspace* workspace            ,       // workspace
          const double               aprecision = 1.e-8   ,       // absolute precision
          const double               rprecision = 1.e-8   ,       // relative precision
          int                        limit      = -1      ,       // limit 
          const char*                reason     = nullptr ,       // message 
          const char*                file       = nullptr ,       // file name 
          const unsigned long        line       = 0       ) const // line number 
        {
          // setup GSL 
          Ostap::Math::GSL::GSL_Error_Handler sentry ;
          //
          double    result =  1.0 ;
          double    error  = -1.0 ;
          if      ( limit <= 0                ) { limit =  workspace->limit() ; }
          else if ( limit >  workspace->limit ) { limit =  workspace->limit() ; }
          //
          const int ierror = gsl_integration_qagil 
            ( func               ,   // the function
              xhigh              ,   // high edges
              aprecision         ,   // absolute precision
              rprecision         ,   // relative precision
              limit              ,   // maximum number of subintervals
              workspace          ,   // workspace
              &result            ,   // the result
              &error             ) ; // the error in result
         if ( ierror ) { gsl_error ( reason , file , line , ierror ) ; }
          //
          return Result { ierror , result , error } ;  
        }
        // ====================================================================
      public:
        // ====================================================================
        /// the actual adapter for GSL 
        static double adapter ( double x , void* params ) 
        {
          const FUNCTION* f = (FUNCTION*) params ;
          return (*f) ( x ) ;
        }
        // ====================================================================
      };  
      // ======================================================================
    } //                                  The end of namespace Ostap::Math::GSL
    // ========================================================================
    /** @class IntegrateX 
     *  helper class to perform X-integration of 2D-funtion
     *  @author Vanya Belyaev
     *  @date   2018-09-21
     */
    template <class FUNCTION2D> 
    class IntegrateX  
    {
    public:
      // ======================================================================
      IntegrateX ( const FUNCTION2D* f2d , 
                   const double      y   ) 
        : m_f2d  ( f2d ) 
        , m_y    ( y   ) 
      {}
      IntegrateX() =delete ;
      // ======================================================================
      double operator() ( const double x ) const 
      { return (*m_f2d) ( x , m_y ) ; }
      // ======================================================================
      const FUNCTION2D* m_f2d ;
      double            m_y   ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class IntegrateY 
     *  helper class to perform Y-integration of 2D-funtion
     *  @author Vanya Belyaev
     *  @date   2018-09-21
     */
    template <class FUNCTION2D> 
    class IntegrateY  
    {
    public:
      // ======================================================================
      IntegrateY ( const FUNCTION2D* f2d , 
                   const double      x   ) 
        : m_f2d  ( f2d ) 
        , m_x    ( x   ) 
      {}
      IntegrateY() = delete ;
      // ======================================================================
      double operator() ( const double y ) const 
      { return (*m_f2d) ( m_x , y ) ; }
      // ======================================================================
      const FUNCTION2D* m_f2d ;
      double            m_x   ;
      // ======================================================================
    } ;
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_INTEGRATOR1D_H
// ============================================================================
