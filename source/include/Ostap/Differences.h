// ============================================================================
#ifndef OSTAP_DIFFERENCES_H 
#define OSTAP_DIFFERENCES_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <type_traits>
#include <functional>
#include <iterator>
// ============================================================================
/** @file Ostap/Differences.h
 *  Collection of classes and functions to deal with the finite differences 
 */
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /** @namespace Ostap::Math::Differences
     *  Collection of classes and functions to deal with the finite differences 
     *  @see https://en.wikipedia.org/wiki/Divided_differences 
     *  @author Vanya Belyaev
     *  @date   2018-07-30
     */
    namespace Differences 
    {
      // ======================================================================
      // Functions for divided differences :
      // ======================================================================
      /** Divided Forward differences of order-0 from the function 
       *  @code
       *  auto  fun = []  ( double x ) { return std::sin(x) ; }
       *  double d0 = divided ( fun , 0.5 ) ;
       *  @endcode
       *  @see https://en.wikipedia.org/wiki/Divided_differences 
       */ 
      template <class FUNCTION, 
		typename std::enable_if<std::is_invocable<FUNCTION,double>::value,int>::type = 0 >      
      inline double divided 
      ( const FUNCTION& fun , 
        const double    x ) { return fun  ( x ) ; }
      // ======================================================================
      /** Divided Forward differences of high order from the function 
       *  @code
       *  auto  fun = []  ( double x ) { return std::sin(x) ; }
       *  double d0 = divided ( fun , 0.5 ) ;
       *  double d1 = divided ( fun , 0.5 , 0.6 ) 
       *  double d2 = divided ( fun , 0.5 , 0.6 , 0.7 )
       *  double d3 = divided ( fun , 0.5 , 0.6 , 0.7 , 0.8 )       
       *  double d4 = divided ( fun , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 )       
       *  @endcode
       *  @see https://en.wikipedia.org/wiki/Divided_differences 
       */ 
      template <class FUNCTION,
		typename... Args, 
		typename std::enable_if<std::is_invocable<FUNCTION,double>::value,int>::type = 0 >
      inline double divided
      ( const FUNCTION& fun  , 
        const double    x0   , 
        const Args...   args , 
        const double    xn   ) 
      {
        return ( divided ( fun , args... , xn ) - 
                 divided ( fun , x0 , args... ) ) /  ( xn -  x0 ) ;
      }
      // ======================================================================
      /** Divided forward differences from two sequences
       *  @param xbegin the start of sequence of abscissas 
       *  @param xend   the end   of sequence of abscissas 
       *  @param ybegin the start of sequence of function values 
       *  @return divided  difference claculated from these sequences
       */
      template <class XITERATOR, 
                class YITERATOR> 
      inline double divided
      ( XITERATOR xbegin , 
        XITERATOR xend   , 
        XITERATOR ybegin )
      {
        const auto N = std::distance ( xbegin , xend ) ;
        //
        if       ( 0 >= N ) { return 0.0     ; } // RETURN
        else  if ( 1 == N ) { return *ybegin ; } // RETURN
        //
        const long double x0 = * ( xbegin          ) ;
        const long double xn = * ( xbegin +  N - 1 ) ;
        //
        return 
          ( divided ( xbegin + 1 , xend           , ybegin + 1 ) * 1.0L - 
            divided ( xbegin     , xbegin + N - 1 , ybegin     ) ) / ( xn - x0 ) ;
      }
      // ======================================================================
      /** Divided forward differences from two sequnces
       *  @param xbegin the start of sequence of abscissas 
       *  @param xend   the end   of sequence of abscissas 
       *  @param ybegin the start of sequence of function values 
       *  @param xvalue adapter to get x-value from dereferenced x-iterator 
       *  @param yvalue adapter to get y-value from dereferenced y-iterator 
       *  @return divided  difference calculated from these sequences
       */
      template <class XITERATOR, 
                class YITERATOR, 
                class XADAPTER ,
                class YADAPTER >
      inline double divided 
      ( XITERATOR xbegin , 
        XITERATOR xend   , 
        YITERATOR ybegin , 
        XADAPTER  xvalue , 
        YADAPTER  yvalue ) 
      {
        const auto N = std::distance ( xbegin , xend ) ;
        //
        if       ( 0 >= N ) { return 0.0                ; } // RETURN
        else  if ( 1 == N ) { return yvalue ( *ybegin ) ; } // RETURN
        //
        const long double x0 = xvalue ( * ( xbegin          ) ) ;
        const long double xn = xvalue ( * ( xbegin +  N - 1 ) ) ;
        //
        const auto cxv = std::cref ( xvalue ) ;
        const auto cyv = std::cref ( yvalue ) ;
        //
        return 
          ( divided ( xbegin + 1 , xend            , ybegin + 1 , cxv , cyv ) * 1.0L - 
            divided ( xbegin     , xbegin + N - 1  , ybegin     , cxv , cyv ) ) 
          / ( xn - x0 ) ;
      }
      // ======================================================================
      // Finite differences
      // ======================================================================
      /** @class Difference
       *  Finite difference
       *  https://en.wikipedia.org/wiki/Finite_difference
       *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
       */
      class FiniteDifference_
      {
	// ====================================================================
      public : // N = 0 
	// ====================================================================
	/// forward difference
	template <class FUNCTION,
		  const unsigned short N, 
		  typename std::enable_if<std::is_invocable<FUNCTION,double>::value,int>::type = 0 ,	  
		  typename std::enable_if<(N==0),int>::type = 0 > 
	static inline double forward
	( const FUNCTION& fun ,
	  const double    x   ,
	  const double    h   = 1 ) { return fun ( x ) ; }
	// ====================================================================
	/// backward difference
	template <class FUNCTION,
		  const unsigned short N, 
		  typename std::enable_if<std::is_invocable<FUNCTION,double>::value,int>::type = 0 ,	  
		  typename std::enable_if<(N==0),int>::type = 0 > 
	static inline double backward 
	( const FUNCTION& fun     ,
	  const double    x       ,
	  const double    h   = 1 ) { return fun ( x ) ; }
	// ====================================================================
	/// central difference
	template <class FUNCTION,
		  const unsigned short N, 
		  typename std::enable_if<std::is_invocable<FUNCTION,double>::value,int>::type = 0 ,	  
		  typename std::enable_if<(N==0),int>::type = 0 > 
	static inline double central 
	( const FUNCTION& fun     ,
	  const double    x       ,
	  const double    h   = 1 ) { return fun ( x ) ; }
	// ====================================================================
      public : // N == 1 
	// ====================================================================
	/// forward difference 
	template <class FUNCTION,
		  const unsigned short N, 
		  typename std::enable_if<std::is_invocable<FUNCTION,double>::value,int>::type = 0 ,	  
		  typename std::enable_if<(N==1),int>::type = 0 > 
	static inline double forward
	( const FUNCTION& fun     ,
	  const double    x       ,
	  const double    h   = 1 ) { return !h ? 0.0 : fun ( x + 1.0L * h ) - 1.0L * fun ( x ) ; }
	// ====================================================================
	/// backward difference
	template <class FUNCTION,
		  const unsigned short N, 
		  typename std::enable_if<std::is_invocable<FUNCTION,double>::value,int>::type = 0 ,	  
		  typename std::enable_if<(N==1),int>::type = 0 > 
	static inline double backward 
	( const FUNCTION& fun     ,
	  const double    x       ,
	  const double    h   = 1 ) { return !h ? 0.0 : fun ( x ) - 1.0L * fun ( x - 1.0L * h ) ; }
	// ====================================================================
	/// central difference
	template <class FUNCTION,
		  const unsigned short N, 
		  typename std::enable_if<std::is_invocable<FUNCTION,double>::value,int>::type = 0 ,	  
		  typename std::enable_if<(N==1),int>::type = 0 > 
	static inline double central 
	( const FUNCTION& fun     ,
	  const double    x       ,
	  const double    h   = 1 ) { return !h ? 0.0 : fun ( x + 0.5L * h ) - 1.0L * fun ( x - 0.5L * h ) ; }
	// ====================================================================
      public : // N == 2 
	// ====================================================================
	/// forward difference 
	template <class FUNCTION,
		  const unsigned short N, 
		  typename std::enable_if<std::is_invocable<FUNCTION,double>::value,int>::type = 0 ,	  
		  typename std::enable_if<(N==2),int>::type = 0 > 
	static inline double forward
	( const FUNCTION& fun     ,
	  const double    x       ,
	  const double    h   = 1 )
	{ return !h ? 0.0 : fun ( x ) - 2.0L * fun ( x + 1.0L * h ) + fun ( x + 2.0L * h )  ; }
	// ====================================================================
	/// backward difference
	template <class FUNCTION,
		  const unsigned short N, 
		  typename std::enable_if<std::is_invocable<FUNCTION,double>::value,int>::type = 0 ,	  
		  typename std::enable_if<(N==2),int>::type = 0 > 
	static inline double backward 
	( const FUNCTION& fun     ,
	  const double    x       ,
	  const double    h   = 1 )
	{ return !h ? 0.0 : fun ( x ) - 2.0L * fun ( x - 1.0L * h ) + fun ( x - 2.0L * h )  ; }
	// ====================================================================
	/// central difference
	template <class FUNCTION,
		  const unsigned short N, 
		  typename std::enable_if<std::is_invocable<FUNCTION,double>::value,int>::type = 0 ,	  
		  typename std::enable_if<(N==2),int>::type = 0 > 
	static inline double central 
	( const FUNCTION& fun     ,
	  const double    x       ,
	  const double    h   = 1 ) 
	{ return !h ? 0.0 : fun ( x + 1.0L * h ) - 2.0L * fun ( x ) + fun ( x - 1.0L * h )  ; }
	// ====================================================================
      public : // N == 3 
	// ====================================================================
	/// forward difference 
	template <class FUNCTION,
		  const unsigned short N, 
		  typename std::enable_if<std::is_invocable<FUNCTION,double>::value,int>::type = 0 ,	  
		  typename std::enable_if<(N==3),int>::type = 0 > 
	static inline double forward
	( const FUNCTION& fun     ,
	  const double    x       ,
	  const double    h   = 1 )
	{ return !h ? 0.0 : - fun ( x ) + 3.0L * fun ( x + 1.0L * h ) - 3.0L * fun (  x + 2.0L * h ) + fun ( x + 3.0L * h ) ; }
	// ====================================================================
	/// backward difference
	template <class FUNCTION,
		  const unsigned short N, 
		  typename std::enable_if<std::is_invocable<FUNCTION,double>::value,int>::type = 0 ,	  
		  typename std::enable_if<(N==3),int>::type = 0 > 
	static inline double backward 
	( const FUNCTION& fun     ,
	  const double    x       ,
	  const double    h   = 1 )
	{ return !h ? 0.0 : + fun ( x ) - 3.0L * fun ( x - 1.0L * h ) + 3.0L * fun (  x - 2.0L * h ) - fun ( x - 3.0L * h ) ; }
	// ====================================================================
	/// central difference
	template <class FUNCTION,
		  const unsigned short N, 
		  typename std::enable_if<std::is_invocable<FUNCTION,double>::value,int>::type = 0 ,	  
		  typename std::enable_if<(N==3),int>::type = 0 > 
	static inline double central 
	( const FUNCTION& fun     ,
	  const double    x       ,
	  const double    h   = 1 ) 
	{ return !h ? 0.0 : + fun ( x + 1.5L * h ) - 3.0L * fun ( x + 0.5L * h ) + 3.0L * fun ( x - 0.5L * h ) - fun ( x - 1.5L * h ) ; } 
	// ====================================================================
      public : // N > 3 
	// ====================================================================
	/// forward difference 
	template <class FUNCTION,
		  const unsigned short N, 
		  typename std::enable_if<std::is_invocable<FUNCTION,double>::value,int>::type = 0 ,	  
		  typename std::enable_if<(N>3),int>::type = 0 > 
	static inline double forward
	( const FUNCTION& fun     ,
	  const double    x       ,
	  const double    h   = 1 )	
	{
	  if ( !h ) { return 0 ; } 
	  long double r = fun ( x ) * ( 0 == N % 2 ? 1 : -1 ) ;
	  long double t = 1 ;
	  for ( unsigned short j = 1 ; j <= N ; ++j )
	  {
	    t *= - ( N + 1 - j ) ;
	    t /= j ;
	    r += t * fun ( x + 1.0L * j * h ) ; 
	  }
	  return r ;
	}
	// ====================================================================
	/// backward difference
	template <class FUNCTION,
		  const unsigned short N, 
		  typename std::enable_if<std::is_invocable<FUNCTION,double>::value,int>::type = 0 ,	  
		  typename std::enable_if<(N>3),int>::type = 0 > 
	static inline double backward 
	( const FUNCTION& fun     ,
	  const double    x       ,
	  const double    h   = 1 )
	{
	  if ( !h ) { return 0 ; } 
	  long double r = fun ( x ) ;
	  long double t = 1 ;
	  for ( unsigned short j = 1 ; j <= N ; ++j )
	  {
	    t *= - ( N + 1 - j ) ;
	    t /= j ;
	    r += t * fun ( x - 1.0L * j * h ) ; 
	  }
	  return r ;
	}	
	// ====================================================================
	/// central difference
	template <class FUNCTION,
		  const unsigned short N, 
		  typename std::enable_if<std::is_invocable<FUNCTION,double>::value,int>::type = 0 ,	  
		  typename std::enable_if<(N>3),int>::type = 0 > 
	static inline double central 
	( const FUNCTION& fun     ,
	  const double    x       ,
	  const double    h   = 1 )
	{
	  if ( !h ) { return 0 ; } 
	  long double r = fun ( x ) ;
	  long double t = 1 ;
	  for ( unsigned short j = 1 ; j <= N ; ++j )
	  {
	    t *= - ( N + 1 - j ) ;
	    t /= j ;
	    r += t * fun ( x + ( 0.5L * N - j ) * h ) ; 
	  }
	  return r ;
	}	  
	// ====================================================================	
      };
      // ======================================================================
      /** @class Forward
       *  simple evaluator of the Nth forward difference
       *  @see https://en.wikipedia.org/wiki/Finite_difference
       */
      template <unsigned short N> 
      class Forward
      {
      public:
        // ====================================================================
        /** get Nth difference
         *  @param fun the fnuction
         *  @param x   x-value 
         *  @param h   step 
         *  @return Nth difference 
         */
	template <class FUNCTION,
		  typename std::enable_if<std::is_invocable<FUNCTION,double>::value,int>::type = 0>
        static inline double evaluate 
        ( const FUNCTION& fun     , 
          const double    x       , 
          const double    h   = 1 )
	{ return FiniteDifference_::forward<N> ( fun , x , h ) ; }
        // ====================================================================
      } ;
      // ======================================================================
      /** @class Backward
       *  simple evaluator of the Nth backward difference
       *  @see https://en.wikipedia.org/wiki/Finite_difference
       */
      template <unsigned short N> 
      class Backward
      {
      public:
        // ====================================================================
        /** get Nth difference
         *  @param fun the fnuction
         *  @param x   x-value 
         *  @param h   step 
         *  @return Nth difference 
         */
	template <class FUNCTION,
		  typename std::enable_if<std::is_invocable<FUNCTION,double>::value,int>::type = 0>
        static inline double evaluate 
        ( const FUNCTION& fun     , 
          const double    x       , 
          const double    h   = 1 )
	{ return FiniteDifference_::backward<N> ( fun , x , h ) ; }
        // ====================================================================
      } ;
      // ======================================================================
      /** @class Central
       *  simple evaluator of the Nth central difference
       *  @see https://en.wikipedia.org/wiki/Finite_difference
       */
      template <unsigned short N> 
      class Central 
      {
      public:
        // ==================================================================== 
        /** get Nth difference
         *  @param fun the fnuction
         *  @param x   x-value 
         *  @param h   step 
         *  @return Nth difference 
         */
	template <class FUNCTION,
		  typename std::enable_if<std::is_invocable<FUNCTION,double>::value,int>::type = 0>
        static inline double evaluate 
        ( const FUNCTION& fun     , 
          const double    x       , 
          const double    h   = 1 )
	{ return FiniteDifference_::central<N> ( fun , x , h ) ; }
        // ====================================================================
      } ;
      // =======================================================================
      /** Evaluate N-th forward difference of function <code>fun</code>
       *  @param fun the function 
       *  @param N   the order of the difference 
       *  @param x   the point 
       *  @param h   the step 
       *  @return N-th forward dirrerence 
       */
      template <class FUNCTION,
		typename std::enable_if<std::is_invocable<FUNCTION,double>::value,int>::type = 0>
      inline double forward_ 
      ( const FUNCTION&      fun     , 
        const unsigned short N       , 
        const double         x       , 
        const double         h   = 1 ) 
      {
	//
	if      ( 0 == N ) { return fun ( x ) ; }
	else if ( !h     ) { return 0         ; }
	else if ( 1 == N ) { return fun ( x + h ) - 1.0L * fun ( x ) ; }
	else if ( 2 == N ) { return fun ( x )     - 2.0L * fun ( x + h ) + fun ( x + 2 * h ) ; }
	//
	long double r = fun ( x ) * ( 0 == N % 2 ? 1 : -1 ) ;
	long double t = 1 ;
	for ( unsigned short j = 1 ; j <= N ; ++j )
	{
	  t *= - ( N + 1 - j ) ;
	  t /= j ;
	  r += t * fun ( x + j * h ) ; 
	}
	return r ;	
      }
      // ======================================================================
      /** Evaluate N-th backward difference of function <code>fun</code>
       *  @param fun the function 
       *  @param N   the order of the difference 
       *  @param x   the point 
       *  @param h   the step 
       *  @return N-th backward dirrerence 
       */
      template <class FUNCTION,
		typename std::enable_if<std::is_invocable<FUNCTION,double>::value,int>::type = 0>
      inline double backward_ 
      ( const FUNCTION&      fun     , 
        const unsigned short N       , 
        const double         x       , 
        const double         h   = 1 ) 
      {
	//
	if      ( 0 == N ) { return fun ( x ) ; }
	else if ( !h     ) { return 0         ; }
	else if ( 1 == N ) { return fun ( x ) - 1.0L * fun ( x - h ) ; }
	else if ( 2 == N ) { return fun ( x ) - 2.0L * fun ( x - h ) + fun ( x - 2 * h ) ; }
	//
	long double r = fun ( x ) ;
	long double t = 1 ;
	for ( unsigned short j = 1 ; j <= N ; ++j )
	{
	  t *= - ( N + 1 - j ) ;
	  t /= j ;
	  r += t * fun ( x - j * h ) ; 
	}
	return r ;	
      }
      // ======================================================================
      /** Evaluate N-th central difference of function <code>fun</code>
       *  @param fun the function 
       *  @param N   the order of the difference 
       *  @param x   the point 
       *  @param h   the step 
       *  @return N-th central dirrerence 
       */
      template <class FUNCTION,
		typename std::enable_if<std::is_invocable<FUNCTION,double>::value,int>::type = 0>
      inline double central_
      ( const FUNCTION&      fun     , 
        const unsigned short N       , 
        const double         x       , 
        const double         h   = 1 ) 
      {
	if      ( 0 == N ) { return fun ( x ) ; }
	else if ( !h     ) { return 0         ; }
	else if ( 1 == N ) { return fun ( x + 0.5L * h ) - 1.0L * fun ( x - 0.5L * h ) ; }
	else if ( 2 == N ) { return fun ( x + 1.0L * h ) - 2.0L * fun ( x ) + fun ( x - 1.0L * h ) ; }
	//
	long double r = fun ( x ) ;
	long double t = 1 ;
	for ( unsigned short j = 1 ; j <= N ; ++j )
	  {
	  t *= - ( N + 1 - j ) ;
	  t /= j ;
	  r += t * fun ( x - j * h ) ; 
	}
	return r ;
      }
      // ======================================================================
      /** Evaluate N-th forward difference of function <code>fun</code>
       *  @param fun the function 
       *  @param N   the order of the difference 
       *  @param x   the point 
       *  @param h   the step 
       *  @return N-th forward dirrerence 
       */
      double forward 
      ( std::function<double(double)> fun     , 
        const unsigned short          N       , 
        const double                  x       , 
        const double                  h   = 1 ) ;
      // ======================================================================
      /** Evaluate N-th backward difference of function <code>fun</code>
       *  @param fun the function 
       *  @param N   the order of the difference 
       *  @param x   the point 
       *  @param h   the step 
       *  @return N-th backward dirrerence 
       */
      double backward 
      ( std::function<double(double)> fun , 
        const unsigned short          N     , 
        const double                  x     , 
        const double                  h = 1 ) ;
      // ======================================================================
      /** Evaluate N-th central difference of function <code>fun</code>
       *  @param fun the function 
       *  @param N   the order of the difference 
       *  @param x   the point 
       *  @param h   the step 
       *  @return N-th central dirrerence 
       */
      double central
      ( std::function<double(double)> fun   , 
        const unsigned short          N     , 
        const double                  x     , 
        const double                  h = 1 ) ;
      // ======================================================================
      /** @class FiniteDifference 
       */
      class FiniteDifference
      {
      public :
	// ===================================================================
	/** construcructor for the order 
	 *  @param N  the order of the finite differenece 
	 */
	FiniteDifference
	( const unsigned short N = 1 ) ;
	// =================================================================
      public :
	// =================================================================
	/// get forward difference 
	template <class FUNCTION, 
		  typename std::enable_if<std::is_invocable<FUNCTION,double>::value,int>::type = 0 >       	
	inline double forward 
	( const FUNCTION& f     ,
	  const double    x     ,
	  const double    h = 1 ) const
	{ return Ostap::Math::Differences::forward_ ( f , m_N , x , h ) ;  }
	// =================================================================
	/// get backward difference 
	template <class FUNCTION, 
		  typename std::enable_if<std::is_invocable<FUNCTION,double>::value,int>::type = 0 >       	
	inline double backward 
	( const FUNCTION& f     ,
	  const double    x     ,
	  const double    h = 1 ) const
	{ return Ostap::Math::Differences::backward_ ( f , m_N , x , h ) ;  }
	// =================================================================
	/// get central difference 
	template <class FUNCTION, 
		  typename std::enable_if<std::is_invocable<FUNCTION,double>::value,int>::type = 0 >       	
	inline double central 
	( const FUNCTION& f     ,
	  const double    x     ,
	  const double    h = 1 ) const
	{ return Ostap::Math::Differences::central_ ( f , m_N , x , h ) ;  }
	// =================================================================	
      public :
	// =================================================================
	/// get forward difference 
	double forward
	( std::function<double(double)> f     , 
	  const double                  x     ,
	  const double                  h = 1 ) const ;
	// =================================================================
	/// get backward difference 
	double backward 
	( std::function<double(double)> f     , 
	  const double                  x     ,
	  const double                  h = 1 ) const ;
	// =================================================================
	/// get central difference 
	double central 
	( std::function<double(double)> f     , 
	  const double                  x     ,
	  const double                  h = 1 ) const ;
	// =================================================================	
      public :
	// =================================================================
	/// get the order 
	inline unsigned short N () const { return m_N ; }
	// =================================================================
      private :
	// =================================================================
	/// the order :
	unsigned short m_N { 1 } ; // the order 
	// =================================================================	
      } ;		
      // ======================================================================      
    } //                          The end of namespace Ostap::Math::Differences  
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_DIFFERENCES_H
// ============================================================================
