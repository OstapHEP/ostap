// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================
#include <tuple>
#include <functional>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Extrema.h"
// ============================================================================
// Local stuff 
// ============================================================================
#include "Extremum1D.h" 
// ============================================================================
/** @file
 *  Implementation file function from Ostap/Extrema.h
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch 
 *  @date 2019-05-14
 */
// ============================================================================
/* Search for the minimum of function at [low,high] using 
 *  golden-section-search rule  
 *  @aparam f the function 
 *  @param low    lower interval edge 
 *  @param high   upper interval edge
 *  @param guess  initial guess 
 *  @param abseps the required absolute precision   
 *  @param releps the required relative precision   
 *  @author Vanya BELYAEV  Ivan.Belyaev@cern.ch
 */
// ============================================================================
double  Ostap::Math::minimum_golden_section
( std::function<double(double)> f      ,
  const double                  low    ,
  const double                  high   ,
  const double                  guess  ,
  const double                  abseps ,
  const double                  releps ) 
{
  //
  typedef std::function<double(double)> FUNCTION ;
  static const Ostap::Math::GSL::Extremum1D<FUNCTION> extremum {} ;
  //
  auto F = extremum.make_function_min ( &f ) ;
  int    ierror = 0 ;
  double result = 0 ;
  double error  = 0 ;
  std::tie ( ierror , result , error ) =
    extremum.optimize_goldensection ( &F      ,
				      std::min ( low , high ) ,
				      std::max ( low , high ) ,
				      guess    ,  
				      abseps   , 
				      releps   ,
				      -1       ,
				      "Ostap::Math::minimum_golden_section" ,
				      __FILE__ ,
				      __LINE__ ) ;
  return result ;
}
// ============================================================================
/** Search for the maximum of function at [a,b] using 
 *  golden-section-search rule  
 *  @aparam f the function 
 *  @param low    lower interval edge 
 *  @param high   upper interval edge
 *  @param guess  initial guess 
 *  @param abseps the required absolute precision   
 *  @param releps the required relative precision   
 *  @author Vanya BELYAEV  Ivan.Belyaev@cern.ch
 */
// ============================================================================
double Ostap::Math::maximum_golden_section
( std::function<double(double)> f      ,
  const double                  low    ,
  const double                  high   ,
  const double                  guess  ,
  const double                  abseps ,
  const double                  releps )
{
  //
  typedef std::function<double(double)> FUNCTION ;
  static const Ostap::Math::GSL::Extremum1D<FUNCTION> extremum {} ;
  //
  auto F = extremum.make_function_max ( &f ) ;
  int    ierror = 0 ;
  double result = 0 ;
  double error  = 0 ;
  std::tie ( ierror , result , error ) =
    extremum.optimize_goldensection ( &F      ,
				      std::min ( low , high ) ,
				      std::max ( low , high ) ,
				      guess    ,  
				      abseps   , 
				      releps   ,
				      -1       ,
				      "Ostap::Math::maximum_golden_section" ,
				      __FILE__ ,
				      __LINE__ ) ;
  return result ;
}
// ============================================================================
/** Search for the minimum of function at [a,b] using 
 *  Brent' rule  
 *  @aparam f the function 
 *  @param low    lower interval edge 
 *  @param high   upper interval edge
 *  @param guess  initial guess 
 *  @param abseps the required absolute precision   
 *  @param releps the required relative precision   
 *  @author Vanya BELYAEV  Ivan.Belyaev@cern.ch
 */
// ============================================================================
double Ostap::Math::minimum_brent 
( std::function<double(double)> f      ,
  const double                  low    ,
  const double                  high   ,
  const double                  guess  ,
  const double                  abseps ,
  const double                  releps )
{
  //
  typedef std::function<double(double)> FUNCTION ;
  static const Ostap::Math::GSL::Extremum1D<FUNCTION> extremum {} ;
  //
  auto F = extremum.make_function_min ( &f ) ;
  int    ierror = 0 ;
  double result = 0 ;
  double error  = 0 ;
  std::tie ( ierror , result , error ) =
    extremum.optimize_brent ( &F      ,
			      std::min ( low , high ) ,
			      std::max ( low , high ) ,
			      guess    ,  
			      abseps   , 
			      releps   ,
			      -1       ,
			      "Ostap::Math::minimum_brent" ,
			      __FILE__ ,
			      __LINE__ ) ;
  return result ;
}
// ============================================================================
/* Search for the maximum of function at [a,b] using 
 *  Brent' rule  
 *  @aparam f the function 
 *  @param low    lower interval edge 
 *  @param high   upper interval edge
 *  @param guess  initial guess 
 *  @param abseps the required absolute precision   
 *  @param releps the required relative precision   
 *  @author Vanya BELYAEV  Ivan.Belyaev@cern.ch
 */
// ============================================================================
double Ostap::Math::maximum_brent 
( std::function<double(double)> f      ,
  const double                  low    ,
  const double                  high   ,
  const double                  guess  ,
  const double                  abseps ,
  const double                  releps )
{
  //
  typedef std::function<double(double)> FUNCTION ;
  static const Ostap::Math::GSL::Extremum1D<FUNCTION> extremum {} ;
  //
  auto F = extremum.make_function_max ( &f ) ;
  int    ierror = 0 ;
  double result = 0 ;
  double error  = 0 ;
  std::tie ( ierror , result , error ) =
    extremum.optimize_brent ( &F      ,
			      std::min ( low , high ) ,
			      std::max ( low , high ) ,
			      guess    ,  
			      abseps   , 
			      releps   ,
			      -1       ,
			      "Ostap::Math::maximum_brent" ,
			      __FILE__ ,
			      __LINE__ ) ;
  return result ;
}
// ============================================================================
/*  Search for the minimum of function at [a,b] using 
 *  variant of Brent’s algorithm which uses 
 *  the safe-guarded step-length algorithm of Gill and Murray.
 *  @aparam f the function 
 *  @param low    lower interval edge 
 *  @param high   upper interval edge
 *  @param guess  initial guess 
 *  @param abseps the required absolute precision   
 *  @param releps the required relative precision   
 *  @author Vanya BELYAEV  Ivan.Belyaev@cern.ch
 */
// ============================================================================
double Ostap::Math::minimum_quad_golden
( std::function<double(double)> f      ,
  const double                  low    ,
  const double                  high   ,
  const double                  guess  ,
  const double                  abseps ,
  const double                  releps )
{
  //
  typedef std::function<double(double)> FUNCTION ;
  static const Ostap::Math::GSL::Extremum1D<FUNCTION> extremum {} ;
  //
  auto F = extremum.make_function_min ( &f ) ;
  int    ierror = 0 ;
  double result = 0 ;
  double error  = 0 ;
  std::tie ( ierror , result , error ) =
    extremum.optimize_quad_golden ( &F      ,
				    std::min ( low , high ) ,
				    std::max ( low , high ) ,
				    guess    ,  
				    abseps   , 
				    releps   ,
				    -1       ,
				    "Ostap::Math::minimum_quad_golden" ,
				    __FILE__ ,
				    __LINE__ ) ;
  return result ;
}
// ============================================================================
/*  Search for the maximum of function at [a,b] using     
 *  variant of Brent’s algorithm which uses 
 *  the safe-guarded step-length algorithm of Gill and Murray.
 *  @aparam f the function 
 *  @param low    lower interval edge 
 *  @param high   upper interval edge
 *  @param guess  initial guess 
 *  @param abseps the required absolute precision   
 *  @param releps the required relative precision   
 *  @author Vanya BELYAEV  Ivan.Belyaev@cern.ch
 */
// ============================================================================
double Ostap::Math::maximum_quad_golden
( std::function<double(double)> f      ,
  const double                  low    ,
  const double                  high   ,
  const double                  guess  ,
  const double                  abseps ,
  const double                  releps )
{
  //
  typedef std::function<double(double)> FUNCTION ;
  static const Ostap::Math::GSL::Extremum1D<FUNCTION> extremum {} ;
  //
  auto F = extremum.make_function_max ( &f ) ;
  int    ierror = 0 ;
  double result = 0 ;
  double error  = 0 ;
  std::tie ( ierror , result , error ) =
    extremum.optimize_quad_golden ( &F      ,
				    std::min ( low , high ) ,
				    std::max ( low , high ) ,
				    guess    ,  
				    abseps   , 
				    releps   ,
				    -1       ,
				    "Ostap::Math::maximum_quad_golden" ,
				    __FILE__ ,
				    __LINE__ ) ;
  return result ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================
