// ============================================================================
#ifndef OSTAP_INTEGRATOR1D_H 
#define OSTAP_INTEGRATOR1D_H 1
// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================
#include <map>
#include <vector>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Hash.h"
#include "Ostap/LinAlg.h"
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
      /** @class Integrator1D  Integrator1D.h 
       *  Helper class to simplify operations 
       *  with GSL numerical integration functions 
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
       *    integrator.qag_integrate ( &F    , 
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
        /// make a function for numerical integration 
        gsl_function 
        make_function 
        ( const FUNCTION* f ) const 
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
        Result qag_integrate   
        ( const gsl_function*        func                 ,           // the function
          const double               xlow                 ,           // low integration edge 
          const double               xhigh                ,           // high integration edge 
          gsl_integration_workspace* workspace            ,           // workspace
          const double               aprecision = s_APRECISION_QAG ,  // absolute precision
          const double               rprecision = s_RPRECISION_QAG ,  // relative precision
          int                        limit      = -1      ,           // limit 
          const char*                reason     = nullptr ,           // message 
          const char*                file       = nullptr ,           // file name 
          const unsigned long        line       = 0                 , // line number          
          const int                  rule       = GSL_INTEG_GAUSS61 , // integration rule 
          const std::size_t          tag        = 0                 ) const  // tag/label 
        {
          // cache? 
          if ( 0 != tag ) { return qag_integrate ( tag        , 
                                                   func       , 
                                                   xlow       , 
                                                   xhigh      , 
                                                   workspace  , 
                                                   aprecision , 
                                                   rprecision ,
                                                   limit      , 
                                                   reason     ,
                                                   file       , 
                                                   line       , 
                                                   rule       ) ; }
          //
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
        Result qagi_integrate   
        ( const gsl_function*        func                 ,           // the function
          gsl_integration_workspace* workspace            ,           // workspace
          const double               aprecision = s_APRECISION_QAGI , // absolute precision
          const double               rprecision = s_RPRECISION_QAGI , // relative precision
          int                        limit      = -1      ,           // limit 
          const char*                reason     = nullptr ,           // message 
          const char*                file       = nullptr ,           // file name 
          const unsigned long        line       = 0       ,           // line number 
          const std::size_t          tag        = 0       ) const     // tag/label 
        {          
          // cache? 
          if ( 0 != tag ) { return qagi_integrate ( tag        , 
                                                    func       , 
                                                    workspace  ,
                                                    aprecision , 
                                                    rprecision ,
                                                    limit      , 
                                                    reason     ,
                                                    file       , 
                                                    line       ) ; }
          //
          // setup GSL 
          Ostap::Math::GSL::GSL_Error_Handler sentry ;
          //
          double    result =  1.0 ;
          double    error  = -1.0 ;
          if      ( limit <= 0                ) { limit =  workspace->limit ; }
          else if ( limit >  workspace->limit ) { limit =  workspace->limit ; }
          //
          const int ierror = gsl_integration_qagi 
            ( const_cast<gsl_function*> ( func ) ,   // the function
              aprecision                         ,   // absolute precision
              rprecision                         ,   // relative precision
              limit                              ,   // maximum number of subintervals
              workspace                          ,   // workspace
              &result                            ,   // the result
              &error                             ) ; // the error in result
          if ( ierror ) { gsl_error ( reason , file , line , ierror ) ; }
          //
          return Result { ierror , result , error } ;  
        }
        // ====================================================================
        /// adaptive integrator 
        Result qagiu_integrate   
        ( const gsl_function*        func                 ,            // the function
          const double               xlow                 ,            // low integration edge 
          gsl_integration_workspace* workspace            ,            // workspace
          const double               aprecision = s_APRECISION_QAGIU , // absolute precision
          const double               rprecision = s_RPRECISION_QAGIU , // relative precision
          int                        limit      = -1      ,            // limit 
          const char*                reason     = nullptr ,            // message 
          const char*                file       = nullptr ,            // file name 
          const unsigned long        line       = 0       ,            // line number 
          const std::size_t          tag        = 0       ) const      // tag/label 
        {          
          // cache? 
          if ( 0 != tag ) { return qagiu_integrate ( tag        , 
                                                     func       , 
                                                     xlow       ,
                                                     workspace  , 
                                                     aprecision , 
                                                     rprecision ,
                                                     limit      , 
                                                     reason     ,
                                                     file       , 
                                                     line       ) ; }
          //
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
        Result qagil_integrate   
        ( const gsl_function*        func                 ,            // the function
          const double               xhigh                ,            // high integration edge 
          gsl_integration_workspace* workspace            ,            // workspace
          const double               aprecision = s_APRECISION_QAGIL , // absolute precision
          const double               rprecision = s_RPRECISION_QAGIL , // relative precision
          int                        limit      = -1      ,            // limit 
          const char*                reason     = nullptr ,            // message 
          const char*                file       = nullptr ,            // file name 
          const unsigned long        line       = 0       ,            // line number 
          const std::size_t          tag        = 0       ) const      // tag/label 
        {
          // cache? 
          if ( 0 != tag ) { return qagil_integrate ( tag        , 
                                                     func       , 
                                                     xhigh      ,
                                                     workspace  ,
                                                     aprecision , 
                                                     rprecision ,
                                                     limit      , 
                                                     reason     ,
                                                     file       , 
                                                     line       ) ; }
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
        /// adaptive integrator 
        Result qagp_integrate   
        ( const gsl_function*        func                 ,           // the  function
          const double               xlow                 ,           // low  integration edge 
          const double               xhigh                ,           // high integration edge 
          const std::vector<double>& pnts                 ,           // known singular points 
          gsl_integration_workspace* workspace            ,           // workspace
          const double               aprecision = s_APRECISION_QAGP , // absolute precision
          const double               rprecision = s_RPRECISION_QAGP , // relative precision
          int                        limit      = -1      ,           // limit 
          const char*                reason     = nullptr ,           // message 
          const char*                file       = nullptr ,           // file name 
          const unsigned long        line       = 0       ,           // line number 
          const std::size_t          tag        = 0       ) const     // tag/label 
        {
          // cache? 
          if ( 0 != tag ) { return qagp_integrate ( tag        , 
                                                    func       , 
                                                    xlow       ,
                                                    xhigh      ,
                                                    pnts       , 
                                                    workspace  ,
                                                    aprecision , 
                                                    rprecision ,
                                                    limit      , 
                                                    reason     ,
                                                    file       , 
                                                    line       ) ; }
          
          // setup GSL 
          Ostap::Math::GSL::GSL_Error_Handler sentry ;
          //
          double    result =  1.0 ;
          double    error  = -1.0 ;
          if      ( limit <= 0                ) { limit =  workspace->limit ; }
          else if ( limit >  workspace->limit ) { limit =  workspace->limit ; }
          //
          std::vector<double> pts {} ; pts.reserve ( pnts.size() + 2 ) ;
          pts.push_back ( xlow  ) ; 
          for ( auto p : pnts ) { if ( xlow < p && p < xhigh ) { pts.push_back ( p ) ; } }
          pts.push_back ( xhigh ) ;
          //
          const int ierror = gsl_integration_qagp
            ( const_cast<gsl_function*> ( func ) ,   // the function
              &pts[0]    , pts.size ()           ,   // known singular points 
              aprecision                         ,   // absolute precision
              rprecision                         ,   // relative precision
              limit                              ,   // maximum number of subintervals
              workspace                          ,   // workspace
              &result                            ,   // the result
              &error                             ) ; // the error in result
          if ( ierror ) { gsl_error ( reason , file , line , ierror ) ; }
          //
          return Result { ierror , result , error } ;  
        }
        // ====================================================================
        /// Cauchy principal value adaptive integrator 
        Result qawc_integrate   
        ( const gsl_function*        func                 ,           // the  function
          const double               xlow                 ,           // low  integration edge 
          const double               xhigh                ,           // high integration edge 
          const double               c                    ,           // Cauchy's point
          gsl_integration_workspace* workspace            ,           // workspace
          const double               aprecision = s_APRECISION_QAWC , // absolute precision
          const double               rprecision = s_RPRECISION_QAWC , // relative precision
          int                        limit      = -1      ,           // limit 
          const char*                reason     = nullptr ,           // message 
          const char*                file       = nullptr ,           // file name 
          const unsigned long        line       = 0       ,           // line number 
          const std::size_t          tag        = 0       ) const     // tag/label 
        {
          // cache? 
          if ( 0 != tag ) { return qawc_integrate ( tag        , 
                                                    func       , 
                                                    xlow       ,
                                                    xhigh      ,
                                                    c          , 
                                                    workspace  ,
                                                    aprecision , 
                                                    rprecision ,
                                                    limit      , 
                                                    reason     ,
                                                    file       , 
                                                    line       ) ; }
          // setup GSL 
          Ostap::Math::GSL::GSL_Error_Handler sentry ;
          //
          double    result =  1.0 ;
          double    error  = -1.0 ;
          if      ( limit <= 0                ) { limit =  workspace->limit ; }
          else if ( limit >  workspace->limit ) { limit =  workspace->limit ; }
          //
          const int ierror = gsl_integration_qawc 
            ( const_cast<gsl_function*> ( func ) ,   // the function
              xlow   , xhigh                     ,   // low & high edges
              c                                  ,   // Cauchy point 
              aprecision                         ,   // absolute precision
              rprecision                         ,   // relative precision
              limit                              ,   // maximum number of subintervals
              workspace                          ,   // workspace
              &result                            ,   // the result
              &error                             ) ; // the error in result
          if ( ierror ) { gsl_error ( reason , file , line , ierror ) ; }
          //
          return Result { ierror , result , error } ;  
        }
        // ====================================================================
        /// double-adaptive integration CQUAD 
        Result cquad_integrate   
        ( const gsl_function*              func                 ,           // the  function
          const double                     xlow                 ,           // low  integration edge 
          const double                     xhigh                ,           // high integration edge 
          gsl_integration_cquad_workspace* workspace            ,           // workspace
          const double                     aprecision = s_APRECISION_CQUAD , // absolute precision
          const double                     rprecision = s_RPRECISION_CQUAD , // relative precision
          const char*                      reason     = nullptr ,           // message 
          const char*                      file       = nullptr ,           // file name 
          const unsigned long              line       = 0       ,           // line number 
          const std::size_t                tag        = 0       ) const     // tag/label 
        {
          // cache? 
          if ( 0 != tag ) { return cquad_integrate ( tag        , 
                                                     func       , 
                                                     xlow       ,
                                                     xhigh      ,
                                                     workspace  ,
                                                     aprecision , 
                                                     rprecision ,
                                                     reason     ,
                                                     file       , 
                                                     line       ) ; }
          // setup GSL 
          Ostap::Math::GSL::GSL_Error_Handler sentry ;
          //
          double      result =  1.0 ;
          double      error  = -1.0 ;
          std::size_t nevals = 0    ;
          //
          const int ierror = gsl_integration_cquad 
            ( const_cast<gsl_function*> ( func ) ,   // the function
              xlow   , xhigh                     ,   // low & high edges
              aprecision                         ,   // absolute precision
              rprecision                         ,   // relative precision
              workspace                          ,   // workspace
              &result                            ,   // the result
              &error                             ,   // the error in result
              &nevals                            ) ; // number of function evaluations 
            if ( ierror ) { gsl_error ( reason , file , line , ierror ) ; }
          //
          return Result { ierror , result , error } ;  
        }
        // ====================================================================
        /// Romberg integration 
        Result romberg_integrate   
        ( const gsl_function*                func                 ,           // the  function
          const double                       xlow                 ,           // low  integration edge 
          const double                       xhigh                ,           // high integration edge 
          gsl_integration_romberg_workspace* workspace            ,           // workspace
          const double                       aprecision = s_APRECISION_ROMBERG , // absolute precision
          const double                       rprecision = s_RPRECISION_ROMBERG , // relative precision
          const char*                        reason     = nullptr ,           // message 
          const char*                        file       = nullptr ,           // file name 
          const unsigned long                line       = 0       ,           // line number 
          const std::size_t                  tag        = 0       ) const     // tag/label 
        {
          // cache? 
          if ( 0 != tag ) { return romberg_integrate ( tag        , 
                                                       func       , 
                                                       xlow       ,
                                                       xhigh      ,
                                                       workspace  ,
                                                       aprecision , 
                                                       rprecision ,
                                                       reason     ,
                                                       file       , 
                                                       line       ) ; }
          // setup GSL 
          Ostap::Math::GSL::GSL_Error_Handler sentry ;
          //
          double      result =  1.0 ;
          std::size_t nevals = 0    ;
          //
          const int ierror = gsl_integration_romberg 
            ( const_cast<gsl_function*> ( func ) ,   // the function
              xlow   , xhigh                     ,   // low & high edges
              aprecision                         ,   // absolute precision
              rprecision                         ,   // relative precision
              &result                            ,   // the result
              &nevals                            ,   // number of function evaluations 
              workspace                          ) ; // workspace
          if ( ierror ) { gsl_error ( reason , file , line , ierror ) ; }
          //
          // imitate the adaptive integration
          if ( GSL_EMAXITER == ierror )
          {
            // split into 3 intervals 
            const double x1 = xlow  + 0.35 * ( xhigh - xlow ) ;
            const double x2 = xhigh - 0.35 * ( xhigh - xlow ) ;
            //
            Result r1 = romberg_integrate ( func , xlow , x1         , workspace , 
                                            0.5 * aprecision  , rprecision , 
                                            reason , file , line     , tag ) ;
            Result r2 = romberg_integrate ( func , x1   , x2         , workspace , 
                                            0.5 * aprecision  , rprecision , 
                                            reason , file , line     , tag ) ;
            Result r3 = romberg_integrate ( func , x2   , xhigh      , workspace , 
                                            0.5 * aprecision  , rprecision , 
                                            reason      , file       , line , tag ) ;
            //
            return Result { std::max ( std::get<0> ( r1 ) , std::max ( std::get<0> ( r2 ) , std::get<0> ( r3 ) ) ) , 
                std::get<1> ( r1 ) + std::get<1> ( r2 ) + std::get<1> ( r3 ) , 
                std::get<2> ( r1 ) + std::get<2> ( r2 ) + std::get<2> ( r3 ) } ;
          }
          //
          /// fake to keep the interface more or less coherent 
          const double error = std::max ( abs ( aprecision          ) , 
                                          abs ( rprecision * result ) ) ;
          //
          return Result { ierror , result , error } ;  
        }
        // ====================================================================
      public: // integration with cache 
        // ====================================================================
        /// adaptive integrator with cache 
        Result qag_integrate   
        ( const std::size_t          tag                  ,  
          const gsl_function*        func                 ,          // the function
          const double               xlow                 ,          // low integration edge 
          const double               xhigh                ,          // high integration edge 
          gsl_integration_workspace* workspace            ,          // workspace
          const double               aprecision = s_APRECISION_QAG , // absolute precision
          const double               rprecision = s_RPRECISION_QAG , // relative precision
          int                        limit      = -1      ,          // limit 
          const char*                reason     = nullptr ,          // message 
          const char*                file       = nullptr ,          // file name 
          const unsigned long        line       = 0       ,          // line number 
          const int                  rule       = GSL_INTEG_GAUSS61 ) const // integration rule 
        {
          // ==================================================================
          static const std::string s_QAG { "QAG" } ;
          const std::size_t key = Ostap::Utils::hash_combiner 
            ( tag , func->params , xlow , xhigh ,  s_QAG , 
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
          Result result = qag_integrate ( func ,
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
        /// adaptive integrator with cache 
        Result qagi_integrate
        ( const std::size_t          tag                  ,
          const gsl_function*        func                 ,           // the function
          gsl_integration_workspace* workspace            ,           // workspace
          const double               aprecision = s_APRECISION_QAGI , // absolute precision
          const double               rprecision = s_RPRECISION_QAGI , // relative precision
          int                        limit      = -1      ,           // limit 
          const char*                reason     = nullptr ,           // message 
          const char*                file       = nullptr ,           // file name 
          const unsigned long        line       = 0       ) const     // line number 
        {
          //
          // ==================================================================
          static const std::string s_QAGI { "QAGI" } ;
          const std::size_t key = Ostap::Utils::hash_combiner 
            ( tag , func->params      ,  s_QAGI     , 
              aprecision , rprecision , 
              limit      , reason     , file , line ) ;
          // ==================================================================
          { // look into the cache ============================================
            CACHE::Lock lock { s_cache.mutex() } ;
            auto it = s_cache->find  ( key ) ;
            if ( s_cache->end() != it ) {  return it->second ; }  // AVOID calculation
            // ================================================================
          } // ================================================================
          // ==================================================================
          // perform numerical inntegration using GSL 
          Result result = qagi_integrate ( func       ,
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
        /// adaptive integrator with cache 
        Result qagiu_integrate
        ( const std::size_t          tag                  ,
          const gsl_function*        func                 ,            // the function
          const double               xlow                 ,            // low integration edge 
          gsl_integration_workspace* workspace            ,            // workspace
          const double               aprecision = s_APRECISION_QAGIU , // absolute precision
          const double               rprecision = s_RPRECISION_QAGIU , // relative precision
          int                        limit      = -1      ,            // limit 
          const char*                reason     = nullptr ,            // message 
          const char*                file       = nullptr ,            // file name 
          const unsigned long        line       = 0       ) const      // line number 
        {
          //
          // ==================================================================
          static const std::string s_QAGIU { "QAGIU" } ;
          const std::size_t key = Ostap::Utils::hash_combiner 
            ( tag , func->params , xlow ,  s_QAGIU , 
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
          Result result = qagiu_integrate ( func       ,
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
        Result qagil_integrate
        ( const std::size_t          tag                  ,
          const gsl_function*        func                 ,            // the function
          const double               xhigh                ,            // upper integration edge 
          gsl_integration_workspace* workspace            ,            // workspace
          const double               aprecision = s_APRECISION_QAGIL , // absolute precision
          const double               rprecision = s_RPRECISION_QAGIL , // relative precision
          int                        limit      = -1      ,            // limit 
          const char*                reason     = nullptr ,            // message 
          const char*                file       = nullptr ,            // file name 
          const unsigned long        line       = 0       ) const      // line number 
        {
          //
          // ==================================================================
          static const std::string s_QAGIL { "QAGIL" } ;
          const std::size_t key = Ostap::Utils::hash_combiner 
            ( tag , func->params , xhigh , s_QAGIL ,  
              aprecision , rprecision    , 
              limit      , reason        ,  file , line ) ;
          // ==================================================================
          { // look into the cache ============================================
            CACHE::Lock lock { s_cache.mutex() } ;
            auto it = s_cache->find  ( key ) ;
            if ( s_cache->end() != it ) {  return it->second ; }  // AVOID calculation
            // ================================================================
          } // ================================================================
          // ==================================================================
          // perform numerical inntegration using GSL 
          Result result = qagil_integrate ( func       ,
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
        /// adaptive integrator with cache 
        Result qagp_integrate   
        ( const std::size_t          tag                  ,           // tag/label
          const gsl_function*        func                 ,           // the  function
          const double               xlow                 ,           // low  integration edge 
          const double               xhigh                ,           // high integration edge 
          const std::vector<double>& pnts                 ,           // knowns singular points 
          gsl_integration_workspace* workspace            ,           // workspace
          const double               aprecision = s_APRECISION_QAGP , // absolute precision
          const double               rprecision = s_RPRECISION_QAGP , // relative precision
          int                        limit      = -1      ,           // limit 
          const char*                reason     = nullptr ,           // message 
          const char*                file       = nullptr ,           // file name 
          const unsigned long        line       = 0       ) const     // line number 
        {
          // ==================================================================
          static const std::string s_QAGP { "QAGP" } ;
          const std::size_t key = Ostap::Utils::hash_combiner 
            ( tag  , func->params       , s_QAGP ,  
              xlow , xhigh , pnts       ,   
              aprecision   , rprecision , 
              limit        , reason     ,  file , line ) ;
          // ==================================================================
          { // look into the cache ============================================
            CACHE::Lock lock { s_cache.mutex() } ;
            auto it = s_cache->find  ( key ) ;
            if ( s_cache->end() != it ) {  return it->second ; }  // AVOID calculation
            // ================================================================
          } // ================================================================
          // ==================================================================
          // perform numerical inntegration using GSL 
          Result result = qagp_integrate ( func       ,
                                           xlow       ,
                                           xhigh      ,
                                           pnts       ,
                                           workspace  , 
                                           aprecision , 
                                           rprecision , 
                                           limit      , 
                                           reason     , 
                                           file       , 
                                           line       ) ;
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
        /// Cauchy principal value adaptive integrator 
        Result qawc_integrate   
        ( const std::size_t          tag                  ,           // tag/label    
          const gsl_function*        func                 ,           // the  function
          const double               xlow                 ,           // low  integration edge 
          const double               xhigh                ,           // high integration edge 
          const double               c                    ,           // Cauchy's point
          gsl_integration_workspace* workspace            ,           // workspace
          const double               aprecision = s_APRECISION_QAWC , // absolute precision
          const double               rprecision = s_RPRECISION_QAWC , // relative precision
          int                        limit      = -1      ,           // limit 
          const char*                reason     = nullptr ,           // message 
          const char*                file       = nullptr ,           // file name 
          const unsigned long        line       = 0       ) const     // line number 
        {
          // ==================================================================
          static const std::string s_GAWC { "GAWC" } ;
          const std::size_t key = Ostap::Utils::hash_combiner 
            ( tag  , func->params       , s_GAWC ,  
              xlow , xhigh , c          ,   
              aprecision   , rprecision , 
              limit        , reason     ,  file , line ) ;
          // ==================================================================
          { // look into the cache ============================================
            CACHE::Lock lock { s_cache.mutex() } ;
            auto it = s_cache->find  ( key ) ;
            if ( s_cache->end() != it ) {  return it->second ; }  // AVOID calculation
            // ================================================================
          } // ================================================================
          // ==================================================================
          // perform numerical integration using GSL 
          Result result = qawc_integrate ( func       ,
                                           xlow       ,
                                           xhigh      ,
                                           c          ,
                                           workspace  , 
                                           aprecision , 
                                           rprecision , 
                                           limit      , 
                                           reason     ,
                                           file       , 
                                           line       ) ;
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
        /// double-adaptive integration CQUAD 
        Result cquad_integrate   
        ( const std::size_t                tag                  ,            // tag/label    
          const gsl_function*              func                 ,            // the  function
          const double                     xlow                 ,            // low  integration edge 
          const double                     xhigh                ,            // high integration edge 
          gsl_integration_cquad_workspace* workspace            ,            // workspace
          const double                     aprecision = s_APRECISION_CQUAD , // absolute precision
          const double                     rprecision = s_RPRECISION_CQUAD , // relative precision
          const char*                      reason     = nullptr ,            // message 
          const char*                      file       = nullptr ,            // file name 
          const unsigned long              line       = 0       ) const      // line number 
        {
          // ==================================================================
          static const std::string s_CQUAD { "CQUAD" } ;
          const std::size_t key = Ostap::Utils::hash_combiner 
            ( tag  , func->params       , s_CQUAD ,  
              xlow , xhigh ,    
              aprecision   , rprecision , 
              reason       , file       , line    ) ;
          // ==================================================================
          { // look into the cache ============================================
            CACHE::Lock lock { s_cache.mutex() } ;
            auto it = s_cache->find  ( key ) ;
            if ( s_cache->end() != it ) {  return it->second ; }  // AVOID calculation
            // ================================================================
          } // ================================================================
          // ==================================================================
          // perform numerical integration using GSL 
          Result result = cquad_integrate ( func       ,
                                            xlow       ,
                                            xhigh      ,
                                            workspace  , 
                                            aprecision , 
                                            rprecision , 
                                            reason     ,
                                            file       , 
                                            line       ) ;
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
        /// Romberg integration
        Result romberg_integrate   
        ( const std::size_t                  tag                  ,            // tag/label    
          const gsl_function*                func                 ,            // the  function
          const double                       xlow                 ,            // low  integration edge 
          const double                       xhigh                ,            // high integration edge 
          gsl_integration_romberg_workspace* workspace            ,            // workspace
          const double                       aprecision = s_APRECISION_ROMBERG , // absolute precision
          const double                       rprecision = s_RPRECISION_ROMBERG , // relative precision
          const char*                        reason     = nullptr ,            // message 
          const char*                        file       = nullptr ,            // file name 
          const unsigned long                line       = 0       ) const      // line number 
        {
          // ==================================================================
          static const std::string s_ROMBERG { "ROMBERG" } ;
          const std::size_t key = Ostap::Utils::hash_combiner 
            ( tag  , func->params       , s_ROMBERG ,  
              xlow , xhigh ,    
              aprecision   , rprecision , 
              reason       , file       , line    ) ;
          // ==================================================================
          { // look into the cache ============================================
            CACHE::Lock lock { s_cache.mutex() } ;
            auto it = s_cache->find  ( key ) ;
            if ( s_cache->end() != it ) {  return it->second ; }  // AVOID calculation
            // ================================================================
          } // ================================================================
          // ==================================================================
          // perform numerical integration using GSL 
          Result result = romberg_integrate ( func       ,
                                              xlow       ,
                                              xhigh      ,
                                              workspace  , 
                                              aprecision , 
                                              rprecision , 
                                              reason     ,
                                              file       , 
                                              line       ) ;
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
      private :
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
        /// the actual integration cache 
        static CACHE              s_cache     ; // integration cache 
        /// integration cache size 
        static const unsigned int s_CACHESIZE ; // cache size 
        // ====================================================================
      };  
      // ======================================================================
      template <class FUNCTION>
      typename Integrator1D<FUNCTION>::CACHE 
      Integrator1D<FUNCTION>::s_cache = Integrator1D<FUNCTION>::CACHE{} ;
      // ======================================================================
      template <class FUNCTION>
      const unsigned int Integrator1D<FUNCTION>::s_CACHESIZE = 50000 ;
      // ======================================================================
    } //                                  The end of namespace Ostap::Math::GSL
    // ========================================================================
    /** @class IntegrateX2 
     *  helper class to perform X-integration of 2D-funtion
     *  @author Vanya Belyaev
     *  @date   2018-09-21
     */
    template <class FUNCTION2D> 
    class IntegrateX2  
    {
    public:
      // ======================================================================
      IntegrateX2 
      ( const FUNCTION2D* f2d , 
        const double      y   ) 
        : m_f2d  ( f2d ) 
        , m_y    ( y   ) 
      {}
      IntegrateX2 () =delete ;
      // ======================================================================
      double operator() ( const double x ) const 
      { return (*m_f2d) ( x , m_y ) ; }
      // ======================================================================
      const FUNCTION2D* m_f2d ;
      double            m_y   ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class IntegrateY2 
     *  helper class to perform Y-integration of 2D-funtion
     *  @author Vanya Belyaev
     *  @date   2018-09-21
     */
    template <class FUNCTION2D> 
    class IntegrateY2  
    {
    public:
      // ======================================================================
      IntegrateY2
      ( const FUNCTION2D* f2d , 
        const double      x   ) 
        : m_f2d  ( f2d ) 
        , m_x    ( x   ) 
      {}
      IntegrateY2 () = delete ;
      // ======================================================================
      double operator() ( const double y ) const 
      { return (*m_f2d) ( m_x , y ) ; }
      // ======================================================================
      const FUNCTION2D* m_f2d ;
      double            m_x   ;
      // ======================================================================
    } ;
    // =======================================================================
    /** @class IntegrateX3 
     *  helper class to perform X-integration of 3D-funtion
     *  @author Vanya Belyaev
     *  @date   2018-09-21
     */
    template <class FUNCTION3D> 
    class IntegrateX3  
    {
    public:
      // ======================================================================
      IntegrateX3
      ( const FUNCTION3D* f3d ,
        const double       y  , 
        const double       z  ) 
        : m_f3d  ( f3d ) 
        , m_y    ( y   ) 
        , m_z    ( z   ) 
      {}
      IntegrateX3 () =delete ;
      // ======================================================================
      double operator() ( const double x ) const 
      { return (*m_f3d) ( x , m_y , m_z ) ; }
      // ======================================================================
      const FUNCTION3D* m_f3d ;
      double            m_y   ;
      double            m_z   ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class IntegrateY3 
     *  helper class to perform Y-integration of 3D-funtion
     *  @author Vanya Belyaev
     *  @date   2018-09-21
     */
    template <class FUNCTION3D> 
    class IntegrateY3  
    {
    public:
      // ======================================================================
      IntegrateY3 
      ( const FUNCTION3D* f3d , 
        const double      x   , 
        const double      z   ) 
        : m_f3d  ( f3d ) 
        , m_x    ( x   ) 
        , m_z    ( z   ) 
      {}
      IntegrateY3 () = delete ;
      // ======================================================================
      double operator() ( const double y ) const 
      { return (*m_f3d) ( m_x , y , m_z ) ; }
      // ======================================================================
      const FUNCTION3D* m_f3d ;
      double            m_x   ;
      double            m_z   ;
      // ======================================================================
    } ;
    // =======================================================================
    /** @class IntegrateZ3 
     *  helper class to perform Z-integration of 3D-funtion
     *  @author Vanya Belyaev
     *  @date   2018-09-21
     */
    template <class FUNCTION3D> 
    class IntegrateZ3  
    {
    public:
      // ======================================================================
      IntegrateZ3
      ( const FUNCTION3D* f3d , 
        const double      x   , 
        const double      y   ) 
        : m_f3d  ( f3d ) 
        , m_x    ( x   ) 
        , m_y    ( y   ) 
      {}
      IntegrateZ3 () = delete ;
      // ======================================================================
      double operator() ( const double z ) const 
      { return (*m_f3d) ( m_x , m_y , z ) ; }
      // ======================================================================
      const FUNCTION3D* m_f3d ;
      double            m_x   ;
      double            m_y   ;
      // ======================================================================
    } ;
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
#endif // OSTAP_INTEGRATOR1D_H
// ============================================================================
//                                                                      The END 
// ============================================================================
