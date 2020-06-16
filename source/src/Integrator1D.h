// ============================================================================
#ifndef OSTAP_INTEGRATOR1D_H 
#define OSTAP_INTEGRATOR1D_H 1
// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================
#include <map>
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
#include "local_gsl.h"
#include "local_hash.h"   // hash_combine 
#include "syncedcache.h"  // the cache 
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
      typedef std::tuple<int,double,double> Result ;
      // ======================================================================
      /** @class Integrator1D  Integrator1D.h 
       *  Helper class to simplify operations 
       *  with GSL numerical integraion functions 
       *
       *  Typical usage 
       *  @code 
       *  Integrator1D<MYOBJECT> integrator {} ;
       *  ...
       *  auto F = integrator.make_function( this ) ;
       *  int    ierror ;
       *  double result ;
       *  double error  ;
       *  std::tie( ierror, result , error ) = 
       *    integrator.gaq_integrate ( &F    , 
       *                               0 , 1 , // low & high edges 
       *                               workspace ( *this ) ) ;
       *  @endcode
       *  @see https://www.gnu.org/software/gsl/doc/html/integration.html
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
          const unsigned long        line       = 0       ,       // line number 
          const int                  rule       = GSL_INTEG_GAUSS51 ) const // integration rule 
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
              rule               ,   // integration rule
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
          if      ( limit <= 0                ) { limit =  workspace->limit ; }
          else if ( limit >  workspace->limit ) { limit =  workspace->limit ; }
          //
          const int ierror = gsl_integration_qagiu 
            ( const_cast<gsl_function*> ( func ) ,   // the function
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
          if      ( limit <= 0                ) { limit =  workspace->limit ; }
          else if ( limit >  workspace->limit ) { limit =  workspace->limit ; }
          //
          const int ierror = gsl_integration_qagil 
            ( const_cast<gsl_function*> ( func ) ,   // the function
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
      public: // integration with cache 
        // ====================================================================
        /// adaptive integrator with cache 
        Result gaq_integrate_with_cache   
        ( const std::size_t          tag                  ,  
          const gsl_function*        func                 ,       // the function
          const double               xlow                 ,       // low integration edge 
          const double               xhigh                ,       // high integration edge 
          gsl_integration_workspace* workspace            ,       // workspace
          const double               aprecision = 1.e-8   ,       // absolute precision
          const double               rprecision = 1.e-8   ,       // relative precision
          int                        limit      = -1      ,       // limit 
          const char*                reason     = nullptr ,       // message 
          const char*                file       = nullptr ,       // file name 
          const unsigned long        line       = 0       ,       // line number 
          const int                  rule       = GSL_INTEG_GAUSS51 ) const // integration rule 
        {
          // ==================================================================
          const std::size_t key = std::hash_combine 
            ( tag , func->params , xlow , xhigh ,  
              aprecision , rprecision , 
              limit      , reason , file , line , rule ) ;
          // ==================================================================
          { // look into the cache ============================================
            CACHE::Lock lock { s_cache.mutex() } ;
            auto it = s_cache->find  ( key ) ;
            if ( s_cache->end() != it ) {  return it->second ; }  // AVOID calculation
            // ================================================================
          } // ================================================================
          // ==================================================================
          // perform numerical integration using GSL 
          Result result = gaq_integrate ( func ,
                                          xlow ,  xhigh , 
                                          workspace     , 
                                          aprecision    , rprecision  ,
                                          limit         , 
                                          reason        , file , line , rule ) ;
          // ==================================================================
          { // update the cache ===============================================
            CACHE::Lock lock  { s_cache.mutex() } ;
            // clear the cache is too large
            if ( s_CACHESIZE < s_cache->size() ) { s_cache->clear() ; }
            // update the cache
            s_cache->insert ( std::make_pair ( key , result ) ) ;
          } // ================================================================
          // ==================================================================
          return result ;
          // ==================================================================
        }
        // ====================================================================
        /// adaptive integrator 
        Result gaqiu_integrate_with_cache 
        ( const std::size_t          tag                  ,
          const gsl_function*        func                 ,       // the function
          const double               xlow                 ,       // low integration edge 
          gsl_integration_workspace* workspace            ,       // workspace
          const double               aprecision = 1.e-8   ,       // absolute precision
          const double               rprecision = 1.e-8   ,       // relative precision
          int                        limit      = -1      ,       // limit 
          const char*                reason     = nullptr ,       // message 
          const char*                file       = nullptr ,       // file name 
          const unsigned long        line       = 0       ) const // line number 
        {
          //
          // ==================================================================
          const std::size_t key = std::hash_combine 
            ( tag , func->params , xlow ,  
              aprecision , rprecision , 
              limit      , reason , file , line ) ;
          // ==================================================================
          { // look into the cache ============================================
            CACHE::Lock lock { s_cache.mutex() } ;
            auto it = s_cache->find  ( key ) ;
            if ( s_cache->end() != it ) {  return it->second ; }  // AVOID calculation
            // ================================================================
          } // ================================================================
          // ==================================================================
          // perform numerical inntegration using GSL 
          Result result = gaqiu_integrate ( func       ,
                                            xlow       ,
                                            workspace  , 
                                            aprecision , rprecision  , 
                                            limit      , 
                                            reason     , file , line ) ;
          // ==================================================================
          { // update the cache ===============================================
            CACHE::Lock lock  { s_cache.mutex() } ;
            // clear the cache is too large
            if ( s_CACHESIZE < s_cache->size() ) { s_cache->clear() ; }
            // update the cache
            s_cache->insert ( std::make_pair ( key , result ) ) ;
          } // ================================================================
          // ==================================================================
          return result ;
          // ==================================================================
        }
        // ====================================================================
        /// adaptive integrator 
        Result gaqil_integrate_with_cache 
        ( const std::size_t          tag                  ,
          const gsl_function*        func                 ,       // the function
          const double               xhigh                ,       // upper integration edge 
          gsl_integration_workspace* workspace            ,       // workspace
          const double               aprecision = 1.e-8   ,       // absolute precision
          const double               rprecision = 1.e-8   ,       // relative precision
          int                        limit      = -1      ,       // limit 
          const char*                reason     = nullptr ,       // message 
          const char*                file       = nullptr ,       // file name 
          const unsigned long        line       = 0       ) const // line number 
        {
          //
          // ==================================================================
          const std::size_t key = std::hash_combine 
            ( tag , func->params , xhigh ,  
              aprecision , rprecision , 
              limit      , reason , file , line ) ;
          // ==================================================================
          { // look into the cache ============================================
            CACHE::Lock lock { s_cache.mutex() } ;
            auto it = s_cache->find  ( key ) ;
            if ( s_cache->end() != it ) {  return it->second ; }  // AVOID calculation
            // ================================================================
          } // ================================================================
          // ==================================================================
          // perform numerical inntegration using GSL 
          Result result = gaqil_integrate ( func       ,
                                            xhigh      ,
                                            workspace  , 
                                            aprecision , rprecision  , 
                                            limit      , 
                                            reason     , file , line ) ;
          // ==================================================================
          { // update the cache ===============================================
            CACHE::Lock lock  { s_cache.mutex() } ;
            // clear the cache is too large
            if ( s_CACHESIZE < s_cache->size() ) { s_cache->clear() ; }
            // update the cache
            s_cache->insert ( std::make_pair ( key , result ) ) ;
          } // ================================================================
          // ==================================================================
          return result ;
          // ==================================================================
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
      private:
        // ====================================================================
        typedef std::map<std::size_t,Result>  MAP   ;
        typedef SyncedCache<MAP>              CACHE ;
        /// the actual integrtaion cache 
        static CACHE              s_cache     ; // integration cache 
        static const unsigned int s_CACHESIZE ; // cache size 
        // ====================================================================
      };  
      // ======================================================================
      template <class FUNCTION>
      typename Integrator1D<FUNCTION>::CACHE 
      Integrator1D<FUNCTION>::s_cache = Integrator1D<FUNCTION>::CACHE{} ;
      // ======================================================================
      template <class FUNCTION>
      const unsigned int Integrator1D<FUNCTION>::s_CACHESIZE = 10000 ;
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
