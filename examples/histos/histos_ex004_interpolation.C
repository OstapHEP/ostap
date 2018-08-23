// ================================================================================
// STD&STL
// ================================================================================
#include <iostream>
#include <iomanip>
#include <random>
#include <functional>
// ================================================================================
// ROOT 
// ================================================================================
#include "TH1.h"
// ================================================================================
// Ostap
// ================================================================================
#include "Ostap/HistoInterpolation.h"
// ================================================================================
/** @file examples/histos/histos_ex004_interpolation.C
 *  Example dealing with various local interpolation functions in C++ 
 *  @see  Ostap::Math::HistoInterpolation
 *  @see  Ostap::Math::HistoInterpolation::interpolate_1D
 *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
 */
// ================================================================================
int histos_ex004_interpolation () 
{
  
  using HI = Ostap::Math::HistoInterpolation ;
  
  TH1F h1 { "h1" , "title" , 10 , 0 , 10 } ;
  
  for ( unsigned short ibin = 1 ; ibin <= h1.GetNbinsX () ;  ++ibin )
  {
    h1.SetBinContent ( ibin , ibin );
    h1.SetBinError   ( ibin , ibin ); 
  }
  
  auto random = std::bind ( std::uniform_real_distribution<double>{ 0.0 , 10.0 } , 
                            std::default_random_engine            {            } ) ;
  
  for ( unsigned short i = 0 ; i < 10 ; ++i ) 
  {
    const double x  = random() ;
    
    const double v  = HI::interpolate_1D ( h1 , x     ) ; // default interpolation
    const double v0 = HI::interpolate_1D ( h1 , x , HI::Nearest   ) ; // no interpolation
    const double v1 = HI::interpolate_1D ( h1 , x , HI::Linear    ) ; // linear 
    const double v2 = HI::interpolate_1D ( h1 , x , HI::Quadratic ) ; // parabolic
    const double v3 = HI::interpolate_1D ( h1 , x , HI::Cubic     ) ; // cubic
    
    std::cout 
      << " x="   << std::setprecision ( 6 ) << std::fixed << x 
      << " \tv=" << std::setprecision ( 6 ) << std::fixed << v  
      << " v0/v1/v2/v3=" 
      << v0 << "/"
      << v1 << "/"
      << v2 << "/"
      << v3 << "/" << std::endl ;
  }
  return 0 ;
}
// ================================================================================
//                                                                          The END 
// ================================================================================
