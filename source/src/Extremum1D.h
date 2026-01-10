// ============================================================================
#ifndef OSTAP_EXTREMUM1D_H 
#define OSTAP_EXTREMUM1D_H 1
// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Math.h"
// ============================================================================
// GSL
// ============================================================================
#include "gsl/gsl_errno.h"
#include "gsl/gsl_min.h"
// ============================================================================
// local
// ============================================================================
#include "local_gsl.h"
// ============================================================================
#include <iostream> 
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    namespace GSL 
    { 
      // ======================================================================
      /** @class Minimizer 
       *  Helper class to allocate/deallocate GSL minimizer 
       */
      class Minimizer 
      {
      public :
		// ====================================================================
		/// constructor: allocate & initialize the minimizer 
		Minimizer 
		( const gsl_function* fun   , 
		  const double        guess ,
		  const double        low   , 
		  const double        high  ,   
		  const gsl_min_fminimizer_type* mtype = gsl_min_fminimizer_brent ) ;
		// ====================================================================
		/// destructor: deallocate the minimier 
		~Minimizer  () ;
		// ====================================================================
      public:
		// ====================================================================
		/// get the minimizer 
		inline gsl_min_fminimizer* minimizer () const { return m_minimizer ; }
		/// conversion to GSL minimizer 
		inline operator gsl_min_fminimizer*  () const { return m_minimizer ; }
		// ===================================================================
      private:
		// ====================================================================
		gsl_min_fminimizer* m_minimizer { nullptr } ;
		// ====================================================================
      } ; //                       The end of class Ostap::Math::GSL::MinSentry 
      // ======================================================================
      /** @class Extremum1D  Extremum1D.h 
       *  Helper class to simplify operations 
       *  with GSL 1D-minimizer/optimizer 
       *
       *  Typical usage 
       *  @code 
       *  Extremum1D<MYOBJECT> extremum {} ;
       *  ...
       *  auto F = integrator.make_function_min ( this ) ;
       *  int    ierror ;
       *  double result ;
       *  double error  ;
       *  std::tie( ierror, result , error ) = 
       *    extremum.optimize_brent ( &F  , 
       *                               0   ,   // low-edge 
       *                               1   ,   // high-edge 
       *                               0.5 ) ; // initial guess 
       *  @endcode
       *  @see https://www.gnu.org/software/gsl/doc/html/integration.html
       *  @author Vanya Belyaev Ivan.Belayev@cern.ch 
       *  @date   2018-09-21
       */
      template <class FUNCTION>
      class Extremum1D
      {
      public :
        //  ===================================================================
        /// make a function for minimization 
        gsl_function 
        make_function_min  
        ( const FUNCTION* f ) const 
        { 
          gsl_function F ;
          F.params   = const_cast<FUNCTION*>( f )  ;
          F.function = &adapter_min ;
          return F ; 
		}
        // ====================================================================
        /// make a function for maximization 
        gsl_function 
        make_function_max
        ( const FUNCTION* f ) const 
        { 
          gsl_function F ;
          F.params   = const_cast<FUNCTION*>( f )  ;
          F.function = &adapter_max ;
          return F ; 
		}
        // ====================================================================
      public : 
        // ====================================================================
		/** perform optimization/minimizaton  of GSL function 
	 	*  using Brent method 
	 	*/
		Result optimize_brent 
		( const gsl_function*      fun        ,
	  	  const double             low        ,
	      const double             high       ,
	      const double             guess      ,
	      const double             aprecision ,
	      const double             rprecision ,
	      const int                limit      = -1      ,        // max iterations 
          const char*              reason     = nullptr ,        // message 
          const char*              file       = nullptr ,        // file name 
          const unsigned long      line       = 0       ) const  // line number                
	   { return optimize
	     ( fun       ,
	       gsl_min_fminimizer_brent , 
	       low        ,
	       high       ,
	       guess      ,
	       aprecision ,
	       rprecision ,
	       limit      ,        // max iterations 
	       reason     ,        // message 
	       file       ,        // file name 
	       line       ) ; } 
    	// ====================================================================
		/** perform optimization/minimizaton  of GSL function 
	 	 *  using golden section search 
	 	 */
		Result optimize_goldensection
		( const gsl_function*      fun        ,
	      const double             low        ,
	      const double             high       ,
	      const double             guess      ,
	      const double             aprecision ,
	      const double             rprecision ,
	      const int                limit      = -1      ,        // max iterations 
          const char*              reason     = nullptr ,        // message 
          const char*              file       = nullptr ,        // file name 
          const unsigned long      line       = 0       ) const  // line number                
		{ return optimize
	      ( fun       ,
	      gsl_min_fminimizer_goldensection , 
	      low        ,
	      high       ,
	      guess      ,
	      aprecision ,
	      rprecision ,
	      limit      ,        // max iterations 
	      reason     ,        // message 
	      file       ,        // file name 
	      line       ) ; } 
    	// ====================================================================
		/** perform optimization/minimizaton  of GSL function 
	 	 *  using variant of Brent’s algorithm which uses 
	 	 *  the safe-guarded step-length algorithm of Gill and Murray.
	     */
	    Result optimize_quad_golden
	    ( const gsl_function*      fun        ,
	      const double             low        ,
	      const double             high       ,
	      const double             guess      ,
	      const double             aprecision ,
	      const double             rprecision ,
	      const int                limit      = -1      ,        // max iterations 
          const char*              reason     = nullptr ,        // message 
          const char*              file       = nullptr ,        // file name 
          const unsigned long      line       = 0       ) const  // line number                
	    { return optimize
	      ( fun       ,
	        gsl_min_fminimizer_quad_golden , 
	        low        ,
	        high       ,
	        guess      ,
	        aprecision ,
	        rprecision ,
	        limit      ,        // max iterations 
	        reason     ,        // message 
	        file       ,        // file name 
	        line       ) ; } 
        // ====================================================================
      private : 
        // ====================================================================
		/** perform optimization/minimizaton  of GSL function 
	 	 *  usinhg the specified GSL algorithms 
	     */
	     Result optimize
		( const gsl_function*            fun        ,
	  	  const gsl_min_fminimizer_type* mtype      , 
	      const double                   low        ,
	      const double                   high       ,
	      const double                   guess      ,
	      const double                   aprecision ,
	      const double                   rprecision ,
	      const int                      limit      = -1      ,        // max iterations 
          const char*                    reason     = nullptr ,        // message 
          const char*                    file       = nullptr ,        // file name 
          const unsigned long            line       = 0       ) const  // line number                
		{
	      /// limits and initial guess 
	      double a = std::min ( low , high ) ;
	  	  double b = std::max ( low , high ) ;
	      double m = a <= guess && guess <= b ? guess : 0.5 * ( a + b ) ;

	  	  /// allocate&initialize  the minimizer
	      const Minimizer minimizer { fun , m , a , b , mtype } ;
	  
		  /// estimaete the maximal number of golden-section steps 
	  	  static const double s_phi  = ( std::sqrt ( 5.0 ) + 1 ) / 2 ;
	  	  static const double s_r    = s_phi - 1 ;
	      static const double s_logr = std::abs ( std::log ( s_r ) )  ; 	  
	  
	  	  const double ap = std::abs ( aprecision ) ;
	      const double rp = std::abs ( rprecision ) ;
	  
	      const double d1 = ap + rp * std::max ( std::abs ( a ) , std::abs ( b ) ) ;
	      const double d2 = std::abs ( std::log ( ( b - a ) / d1 ) / s_logr ) ;
	  
	  	  unsigned short Nmax = Ostap::Math::round ( 2 * d2 + 2 ) ;
	  	  if ( 0 < limit && limit < Nmax ) { Nmax = limit ; }

	  	  /// start iterations
	      int ierror = GSL_SUCCESS ;
	      for ( unsigned int i = 0 ; i < Nmax ; ++i ) 
	    	{
	      	   ierror = gsl_min_fminimizer_iterate ( minimizer  ) ;
	      	   if ( GSL_SUCCESS != ierror ) { break ; }                 // BREAK! 
	      	   //
	      	   m = gsl_min_fminimizer_x_minimum ( minimizer ) ;
	      	   a = gsl_min_fminimizer_x_lower   ( minimizer ) ;
	           b = gsl_min_fminimizer_x_upper   ( minimizer ) ;
	           //
	      	   ierror = gsl_min_test_interval ( a , b , ap , rp ) ;
	      	   if      ( GSL_SUCCESS  == ierror ) { break    ; }	        // BREAK! 
	           else if ( GSL_CONTINUE == ierror ) { continue ; }	        // CONTINUE! 
	    	}
	  	  ///
	      const double result = m ;
	      const double error  = std::min ( std::abs ( a - m ) , std::abs ( b - m ) ) ;
	      //
          if ( ierror ) { gsl_error ( reason , file , line , ierror ) ; }
          //
	      return Result { ierror , result , error } ;  
         }
        // ====================================================================	
      public:
        // ====================================================================
        /// the actual adapter for GSL minimization  
        static double adapter_min ( double x , void* params ) 
        {
          const FUNCTION* f = (FUNCTION*) params ;
          return   (*f) ( x ) ;
        }
        // ====================================================================
        /// the actual adapter for GSL maximization   
        static double adapter_max ( double x , void* params ) 
        {
          const FUNCTION* f = (FUNCTION*) params ;
          return - (*f) ( x ) ;
        }
        // ====================================================================
      } ; //                        The end of class Ostap::Math:GSL:EXtremum1D   
      // ======================================================================
    } //                                  The end of namespace Ostap::Math::GSL
    // ========================================================================
	template <class FUNCTION>
	inline double mode 
	( const FUNCTION*     fun        , 
	  const double        low        , 
	  const double        high       ,
	  const double        guess      ,
	  const double        aprecision = -1      , 
	  const double        rprecision = -1      , 
	  const int           limit      = 1000    , 
	  const char*         message    = nullptr ,
	  const char*         file       = nullptr , 
	  const unsigned long line       = 0 )
	  {
		/// 
		static const Ostap::Math::GSL::Extremum1D<FUNCTION> s_extremum {} ;
  		//
  		const auto F    = s_extremum.make_function_max ( fun  ) ;
  		//
  		int    ierror   =  0 ;
  		double result   =  1 ;
  		double error    = -1 ;
  		//
		/// the width of search window 
		const double w = std::abs ( high - low ) ;

		// adjust absolute and relative precisions 
		const double rp = 0 < rprecision && rprecision < 0.01     ? rprecision : 1.e-6          ;
		const double ap = 0 < aprecision && aprecision < 0.01 * w ? aprecision : w * rprecision ;   

		std::tie ( ierror , result , error ) = s_extremum.optimize_quad_golden 
    	( &F        , // the function 
      	  low       , // low_value 
          high      , // high edge
          guess     , // initial guess 
          ap        , // absolute precision
          rp        , // relative precision
          limit     , // limit on number of iterations
          message   , // 
		  file      , //
		  line      ) ;   
  		//
  		return result ;
	  }  
	// ========================================================================	
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_EXTREMUM1D_H
// ============================================================================
