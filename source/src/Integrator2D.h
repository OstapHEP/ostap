// ============================================================================
#ifndef OSTAP_INTEGRATOR2D_H 
#define OSTAP_INTEGRATOR2D_H 1
// ============================================================================
// Include  files
// ============================================================================
// Local 
// ============================================================================
#include "Integrator1D.h" //  GSL-integrator 
#include "cubature.h"     //  cubature 
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
       *  @endcode
       *  @author Vanya Belyaev
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
        Fun make_function ( const FUNCTION* f                     , 
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
          const unsigned      maxcalls   = 20000   ,
          const double        aprecision = 1.e-8   , 
          const double        rprecision = 1.e-8   ,
          const char*         reason     = nullptr ,       // message 
          const char*         file       = nullptr ,       // file name 
          const unsigned long line       = 0       ) const // line number 
        {
          double result =  1 ;        
          double error  = -1 ;
          const int ierror = hcubature 
            ( 1 , fun -> fun   , fun->fdata , // f-dimension, function  & data 
              2 , fun -> min   , fun->max   , //   dimension and integration range 
              maxcalls         ,              // maximal number of  function calls 
              aprecision       ,              // absolute precision 
              rprecision       ,              // relative precision
              ERROR_INDIVIDUAL ,              // error norm 
              &result, &error  ) ;            // output: result& error
          //
          if ( ierror )  { gsl_error( reason , file  , line , ierror ) ; }
          return Result { ierror , result , error } ;
        }
        // ====================================================================
      public:
        // ====================================================================
        /// the actual adapter for cubature 
        static int adapter2d ( unsigned      ndim  , 
                               const double* x     , 
                               void*         fdata ,
                               unsigned      fdim  , 
                               double*       fval  )   
        {
          if ( 1  != fdim || 2 != ndim || 
               nullptr == x || nullptr == fdata || nullptr == fval ) { return 1 ; }
          const FUNCTION* f = (FUNCTION*) fdata  ; 
          fval[0] = (*f) ( x[0] , x[1] ) ;
          return 0 ;
        }
        // ====================================================================
      };  
      // ======================================================================
    } //                                  The end of namespace Ostap::Math::GSL 
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_INTEGRATOR2D_H
// ============================================================================
