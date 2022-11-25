#ifndef OSTAP_INTEGRATOR3D_H 
#define OSTAP_INTEGRATOR3D_H 1
// ============================================================================
// Include  files
// ============================================================================
// STD&STL
// ============================================================================
#include <map>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Hash.h"
// ============================================================================
// Local 
// ============================================================================
#include "Integrator1D.h"     // GSL-integrator 
#include "cubature.h"         // cubature 
#include "syncedcache.h"      // the cache 
#include "local_hash.h"       // hash_combine 
#include "local_gsl.h"        // hash_combine 
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
      /** @class Integrator3D  Integrator3D.h 
       *  Helper class to simplify operations with integration of 3D-functions 
       *  @see https://www.gnu.org/software/gsl/doc/html/integration.html
       *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
       *  @date   2022-06-23
       */
      template <class FUNCTION>
      class Integrator3D
      {
      public:
        // ====================================================================
        struct Fun 
        {
          integrand  fun     ;
          void*      fdata   ;
          double     min [3] ;
          double     max [3] ;
        } ;  
        // ====================================================================
        Fun make_function ( const FUNCTION* f                     , 
                            const double xmin , const double xmax , 
                            const double ymin , const double ymax , 
                            const double zmin , const double zmax ) const 
        {
          Fun F ;
          F.fdata   = const_cast<FUNCTION*>( f ) ;
          F.fun     = &adapter3d ;
          F.min [0] = xmin ;
          F.min [1] = ymin ;
          F.min [2] = zmin ;
          F.max [0] = xmax ;
          F.max [1] = ymax ;
          F.max [2] = zmax ;
          return F ;
        } ;
        // ====================================================================
      public:
        // ====================================================================
        Result cubature 
        ( const Fun*          fun                  , 
          const unsigned      maxcalls   = 50000   ,
          const double        aprecision = s_APRECISION , 
          const double        rprecision = s_RPRECISION ,
          const char*         reason     = nullptr ,       // message 
          const char*         file       = nullptr ,       // file name 
          const unsigned long line       = 0       ,       // line number 
          const std::size_t   tag        = 0       ) const // tag/label
        {
          if ( 0 != tag ) { return cubature ( tag        , 
                                              fun        , 
                                              maxcalls   , 
                                              aprecision , 
                                              rprecision , 
                                              reason     , 
                                              file       , 
                                              line       ) ; }
          //
          double result =  1 ;        
          double error  = -1 ;
          const int ierror = hcubature 
            ( 1 , fun -> fun   , fun->fdata , // f-dimension, function  & data 
              3 , fun -> min   , fun->max   , // dimension and integration range 
              maxcalls         ,              // maximal number of  function calls 
              aprecision       ,              // absolute precision 
              rprecision       ,              // relative precision
              ERROR_INDIVIDUAL ,              // error norm 
              &result, &error  ) ;            // output: result&error
          //
          if ( ierror )  { gsl_error( reason , file  , line , ierror ) ; }
          return Result { ierror , result , error } ;
        }
        // ====================================================================
        Result cubature
        ( const std::size_t   tag                  ,
          const Fun*          fun                  , 
          const unsigned      maxcalls   = 50000   ,
          const double        aprecision = s_APRECISION , 
          const double        rprecision = s_RPRECISION ,
          const char*         reason     = nullptr ,       // message 
          const char*         file       = nullptr ,       // file name 
          const unsigned long line       = 0       ) const // line number 
        {
          //
          const std::size_t key = Ostap::Utils::hash_combiner 
            ( tag  , fun->fdata  , 
              fun->min[0] , fun->min[1] , fun->min[2] , 
              fun->max[0] , fun->max[1] , fun->max[2] , 
              maxcalls    , aprecision  , rprecision  , 
              reason      , file        , line        ) ;
          // ==================================================================
          { // look into the cache ============================================
            CACHE::Lock lock { s_cache.mutex() } ;
            auto it = s_cache->find ( key ) ;
            if ( s_cache->end() != it ) {  return it->second ; }  // AVOID calculation
            // ================================================================
          } // ================================================================
          // ==================================================================
          // perform the numerical integration via the cubature method 
          Result result = cubature ( fun      , 
                                     maxcalls , aprecision , rprecision  , 
                                     reason   , file       , line        ) ;
          // ==================================================================
          { // update the cache ===============================================
            CACHE::Lock lock  { s_cache.mutex() } ;
            // clear the cache if too large
            if ( s_CACHESIZE < s_cache->size() ) { s_cache->clear() ; }
            // update the cache
            auto r = s_cache->insert ( std::make_pair ( key , result ) ) ;
          } // ================================================================
          // ==================================================================
          return result ;
          // ==================================================================
        }
        // ====================================================================
      public:
        // ====================================================================
        /// the actual adapter for cubature 
        static int adapter3d 
        ( unsigned      ndim  , 
          const double* x     , 
          void*         fdata ,
          unsigned      fdim  , 
          double*       fval  )   
        {
          if ( 1       != fdim  || 
               3       != ndim  || 
               nullptr == x     || 
               nullptr == fdata || 
               nullptr == fval  ) { return 1 ; }
          const FUNCTION* f = (FUNCTION*) fdata  ; 
          fval [ 0 ] = (*f) ( x [ 0 ] , x [ 1 ] , x [ 2 ] ) ;
          return 0 ;
        }
        // ====================================================================
      private:
        // ====================================================================
        typedef std::map<std::size_t,Result>  MAP   ;
        typedef SyncedCache<MAP>              CACHE ;
        /// the actual integration cache 
        static CACHE              s_cache     ; // integration cache 
        static const unsigned int s_CACHESIZE ; // cache size 
        // ====================================================================
      };  
      // ======================================================================
      template <class FUNCTION>
      typename Integrator3D<FUNCTION>::CACHE 
      Integrator3D<FUNCTION>::s_cache = Integrator3D<FUNCTION>::CACHE{} ;
      // ======================================================================
      template <class FUNCTION>
      const unsigned int Integrator3D<FUNCTION>::s_CACHESIZE = 75000 ;
      // ======================================================================
    } //                                  The end of namespace Ostap::Math::GSL 
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_INTEGRATOR3D_H
// ============================================================================
