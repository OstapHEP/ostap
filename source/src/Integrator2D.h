// ============================================================================
#ifndef OSTAP_INTEGRATOR2D_H 
#define OSTAP_INTEGRATOR2D_H 1
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
      /** @class Integrator2D  Integrator2D.h 
       *  Helper class to simplify operations with integration of 2D-functions 
       *  @see https://www.gnu.org/software/gsl/doc/html/integration.html
       *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
       *  @date   2018-09-21
       */
      template <class FUNCTION>
      class Integrator2D
      {
      public:
        // ====================================================================
        struct Fun 
        {
          integrand  fun     ;
          void*      fdata   ;
          double     min [2] ;
          double     max [2] ;
        } ;  
        // ====================================================================
        Fun make_function 
        ( const FUNCTION* f                     , 
          const double xmin , const double xmax , 
          const double ymin , const double ymax ) const 
        {
          Fun F ;
          F.fdata   = const_cast<FUNCTION*>( f ) ;
          F.fun     = &adapter2d ;
          F.min [0] = xmin ;
          F.min [1] = ymin ;
          F.max [0] = xmax ;
          F.max [1] = ymax ;
          return F ;
        } ;
        // ====================================================================
      public:
        // ====================================================================
        Result cubature 
        ( const Fun*          fun                  , 
          const unsigned      maxcalls   = 50000   ,
          const double        aprecision = s_APRECISION_CUBE2D , 
          const double        rprecision = s_RPRECISION_CUBE2D ,
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
              2 , fun -> min   , fun->max   , // dimension and integration range 
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
          const double        aprecision = s_APRECISION_CUBE2D , 
          const double        rprecision = s_RPRECISION_CUBE2D ,
          const char*         reason     = nullptr ,       // message 
          const char*         file       = nullptr ,       // file name 
          const unsigned long line       = 0       ) const // line number 
        {
          //
          const std::size_t key = Ostap::Utils::hash_combiner 
            ( tag  , fun->fdata  , 
              fun->min[0] , fun->min[1] , fun->max[0] , fun->max[1] ,
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
        static int adapter2d 
        ( unsigned      ndim  , 
          const double* x     , 
          void*         fdata ,
          unsigned      fdim  , 
          double*       fval  )   
        {
          if ( 1       != fdim  || 
               2       != ndim  || 
               nullptr == x     || 
               nullptr == fdata || 
               nullptr == fval  ) { return 1 ; }
          const FUNCTION* f = (FUNCTION*) fdata  ; 
          fval [ 0 ] = (*f) ( x [ 0 ] , x [ 1 ] ) ;
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
      typename Integrator2D<FUNCTION>::CACHE 
      Integrator2D<FUNCTION>::s_cache = Integrator2D<FUNCTION>::CACHE{} ;
      // ======================================================================
      template <class FUNCTION>
      const unsigned int Integrator2D<FUNCTION>::s_CACHESIZE = 50000 ;
      // ======================================================================
    } //                                  The end of namespace Ostap::Math::GSL 
    // ========================================================================
    /** @class IntegrateXY 
     *  helper class to perform XY-integration of 3D-funtion
     *  @author Vanya Belyaev
     *  @date   2018-09-21
     */
    template <class FUNCTION3D> 
    class IntegrateXY  
    {
    public:
      // ======================================================================
      IntegrateXY ( const FUNCTION3D* f3d , 
                    const double      z   ) 
        : m_f3d  ( f3d ) 
        , m_z    ( z   ) 
      {}
      IntegrateXY () =delete ;
      // ======================================================================
      double operator() ( const double x , const double y ) const 
      { return (*m_f3d) ( x , y ,  m_z ) ; }
      // ======================================================================
      const FUNCTION3D* m_f3d ;
      double            m_z   ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class IntegrateXZ 
     *  helper class to perform XZ-integration of 3D-funtion
     *  @author Vanya Belyaev
     *  @date   2018-09-21
     */
    template <class FUNCTION3D> 
    class IntegrateXZ  
    {
    public:
      // ======================================================================
      IntegrateXZ ( const FUNCTION3D* f3d , 
                    const double      y   ) 
        : m_f3d  ( f3d ) 
        , m_y    ( y   ) 
      {}
      IntegrateXZ () =delete ;
      // ======================================================================
      double operator() ( const double x , const double z ) const 
      { return (*m_f3d) ( x , m_y ,  z ) ; }
      // ======================================================================
      const FUNCTION3D* m_f3d ;
      double            m_y   ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class IntegrateYZ 
     *  helper class to perform YZ-integration of 3D-funtion
     *  @author Vanya Belyaev
     *  @date   2018-09-21
     */
    template <class FUNCTION3D> 
    class IntegrateYZ  
    {
    public:
      // ======================================================================
      IntegrateYZ ( const FUNCTION3D* f3d , 
                    const double      x   ) 
        : m_f3d  ( f3d ) 
        , m_x    ( x   ) 
      {}
      IntegrateYZ () =delete ;
      // ======================================================================
      double operator() ( const double y , const double z ) const 
      { return (*m_f3d) ( m_x , y ,  z ) ; }
      // ======================================================================
      const FUNCTION3D* m_f3d ;
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
#endif // OSTAP_INTEGRATOR2D_H
// ============================================================================
