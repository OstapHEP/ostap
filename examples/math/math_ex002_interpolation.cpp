//usr/bin/env root.exe -l ${0}\(\""${0}"\",\""${*}"\"\); exit $?
// ================================================================================
// STD&STL
// ================================================================================
#include <iostream>
#include <cmath>
#include <map>
// ================================================================================
// Ostap
// ================================================================================
#include "Ostap/Interpolation.h"
#include "Ostap/Interpolants.h"
#include "Ostap/Bernstein.h"
#include "Ostap/StatEntity.h"
// ================================================================================
/** @file examples/math/math_ex002_interpolation.cpp
 *  Example dealing with various interpolation functions in C++ 
 *  @see  Ostap::Math::Interpolation
 *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
 */
// ================================================================================
int math_ex002_interpolation () 
{
  
  using Ostap::Math::Interpolation::Abscissas  ;
  using Ostap::Math::Interpolation::lagrange   ;
  using Ostap::Math::Interpolation::bernstein  ;
  using Ostap::Math::Interpolation::bernstein_ ;
  using Ostap::Math::Interpolation::newton     ;
  using Ostap::Math::Interpolation::TABLE      ; 
  
  // the function to be interpolated 
  auto fun = [] ( double x ) { return std::sin ( x ) ; } ;
  
  // number of interpolation points 
  const unsigned short N    = 10 ;
  // low edge of interpolation region 
  const double         low  =  0 ;
  // high edge of interpolation region 
  const double         high =  4 ;

  
  // ==========================================================================
  // (A) get interpolation data directly the function 
  // ==========================================================================

  // 0) Barycentric innterpolant with uniform abscissas 
  auto l0 = lagrange   ( fun , Abscissas ( N , low , high , Abscissas::Uniform   ) ) ;
  
  // 1) Barycentric interpolant with Chebyshev abscissas 
  auto l1 = lagrange   ( fun , Abscissas ( N , low , high , Abscissas::Chebyshev ) ) ;
  
  // 2) Barycentric interpolant with Lobatto abscissas 
  auto l2 = lagrange   ( fun , Abscissas ( N , low , high , Abscissas::Lobatto   ) ) ;
  
  // 3) Barycentric intepolant with given abscissas 
  auto l3 = lagrange   ( fun , { 0.0 , 0.3 , 0.6 , 0.8 , 1.5 , 2.0 , 2.5 , 3.0 , 3.5 , 4.0 } ) ;
  
  // 4) Bernstein interpolant with given abscissas 
  auto l4 = bernstein_ ( fun , { 0.0 , 0.3 , 0.6 , 0.8 , 1.5 , 2.0 , 2.5 , 3.0 , 3.5 , 4.0 } , low , high ) ;

  // ==========================================================================
  // (B,C) get interpolation data in a form of std::map or Interpolation::Table
  // ==========================================================================

  using MAP = std::map<double,double> ;

  TABLE table ;
  MAP   map   ;

  const std::array<double,10> xarr {{ 0.0 , 0.3 , 0.6 , 0.8 , 1.5 , 
                                      2.0 , 2.5 , 3.0 , 3.5 , 4.0 }} ;
  for ( auto x : xarr) 
  {
    const  double fx = fun ( x )  ;
    table.push_back ( std::make_pair ( x , fx ) )  ;
    map [ x ] = fx ;
  }
  
  // 5) Barycentric interpolant from std::map  
  auto l5 = lagrange ( map ) ;

  // 6) Barycentric interpolant from the interpolation table 
  auto l6 = lagrange ( table ) ;

  // 7) Newton interpolant with given abscissas  
  auto l7 = newton   ( fun , { 0.0 , 0.3 , 0.6 , 0.8 , 1.5 , 2.0 , 2.5 , 3.0 , 3.5 , 4.0 } ) ;


  // ==========================================================================

  // counters 
  Ostap::StatEntity c0 {} ;
  Ostap::StatEntity c1 {} ;
  Ostap::StatEntity c2 {} ;
  Ostap::StatEntity c3 {} ;
  Ostap::StatEntity c4 {} ;
  Ostap::StatEntity c5 {} ;
  Ostap::StatEntity c6 {} ;
  Ostap::StatEntity c7 {} ;
  
  const unsigned int nsteps = 100000 ;
  const double       dx     = ( high - low ) / nsteps ;
  
  for ( unsigned  int i = 0 ; i <= nsteps ; ++i ) 
  {
    const double x = low + i * dx ;
    //
    // true value of function 
    const double f  = fun ( x ) ;
    // l0-interpolant
    const double f0 = l0  ( x ) ;
    // l1-interpolant
    const double f1 = l1  ( x ) ;
    // l2-interpolant
    const double f2 = l2  ( x ) ;
    // fixed abscissas 
    const double f3 = l3  ( x ) ;
    // berinstein interpolant 
    const double f4 = l4  ( x ) ;
    // interpolant5from table 
    const double f5 = l5  ( x ) ;
    // interpolant from map
    const double f6 = l6  ( x ) ;
    // newton interpolant 
    const double f7 = l7  ( x ) ;

    //
    
    c0 += std::abs ( f0 - f ) ;
    c1 += std::abs ( f1 - f ) ;
    c2 += std::abs ( f2 - f ) ;
    c3 += std::abs ( f3 - f ) ;
    c4 += std::abs ( f4 - f ) ;
    c5 += std::abs ( f5 - f ) ;
    c6 += std::abs ( f6 - f ) ;
    c7 += std::abs ( f7 - f ) ;
    //
  }
  const double scale = 1.e-8 ;
  //
  std::cout << "Uniform   : [" << scale 
            << "] mean = " << c0.mean()/scale << "+-" << c0.rms()/scale 
            <<  "  \tmax=" << c0.max ()/scale << std::endl ;
  std::cout << "Chebyshev : [" << scale 
            << "] mean = " << c1.mean()/scale << "+-" << c1.rms()/scale 
            <<  "  \tmax=" << c1.max ()/scale << std::endl ;
  std::cout << "Lobatto   : [" << scale  
            << "] mean = " << c2.mean()/scale << "+-" << c2.rms()/scale 
            <<  "  \tmax=" << c2.max ()/scale << std::endl ;
  std::cout << "Fixed     : [" << scale 
            << "] mean = " << c3.mean()/scale << "+-" << c3.rms()/scale 
            <<  "  \tmax=" << c3.max ()/scale << std::endl ;
  std::cout << "Bernstein : [" << scale 
            << "] mean = " << c4.mean()/scale << "+-" << c4.rms()/scale 
            <<  "  \tmax=" << c4.max ()/scale << std::endl ;
  std::cout << "Map       : [" << scale 
            << "] mean = " << c5.mean()/scale << "+-" << c5.rms()/scale 
            <<  "  \tmax=" << c5.max ()/scale << std::endl ;
  std::cout << "Table     : [" << scale 
            << "] mean = " << c6.mean()/scale << "+-" << c6.rms()/scale 
            <<  "  \tmax=" << c6.max ()/scale << std::endl ;
  std::cout << "Newton    : [" << scale 
            << "] mean = " << c7.mean()/scale << "+-" << c7.rms()/scale 
            <<  "  \tmax=" << c7.max ()/scale << std::endl ;
  //
  return 0 ;
}
// ================================================================================
// The END 
// ================================================================================
