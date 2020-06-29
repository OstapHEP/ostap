// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================
#include <functional>
#include <algorithm>
#include <vector>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Piecewise.h"
// ============================================================================
// local
// ============================================================================
#include "Exception.h"
// ============================================================================
/** @file
 *  Implementation file for class Ostap::Math::Piecewise
 *  @date 2020-06-29 
 *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
 */
// ============================================================================
namespace
{
  // ==========================================================================
  const std::string s_ERROR   { "Cannot insert new value!"    } ;
  const std::string s_METHOD  { "Ostap::Math::Piecewise::add" } ;
  const std::string s_METHOD3 { "Ostap::Math::Piecewise(3)"   } ;
  const std::string s_METHOD4 { "Ostap::Math::Piecewise(4)"   } ;
  // ==========================================================================
}
// ============================================================================
// Standard constructor with a single function 
// ============================================================================
Ostap::Math::Piecewise::Piecewise 
( const  std::function<double(double)>& f1 )
  : m_edges (        ) 
  , m_funcs ( 0 , f1 ) 
{}
// ============================================================================
/*  Standard constructor with two functions 
 *  @param f1 function to bve uses for    \f$ x<   x_1\f$
 *  @param x1 the vale of \f$x_1\f$ 
 *  @param f2 function to bve uses for    \f$ x\ge x_1\f$
 */
// ============================================================================
Ostap::Math::Piecewise::Piecewise 
( const std::function<double(double)>& f1 ,
  const double                         x1 , 
  const std::function<double(double)>& f2 ) 
  : m_edges ( 1 , x1 ) 
  , m_funcs {{ f1, f2 }}
{}
// ============================================================================
/*  Standard constructor with three functions 
 *  @param f1 function to be used for    \f$ x<  x_1\f$
 *  @param x1 the value of \f$x_1\f$ 
 *  @param f2 function to be used for    \f$ x_1 \le  < x_2\f$
 *  @param x2 the value of \f$x_2\f$ 
 *  @param f3 function to be used for    \f$ x_2 \le      \f$
 */
// ============================================================================
Ostap::Math::Piecewise::Piecewise 
( const std::function<double(double)>& f1 ,
  const double                         x1 , 
  const std::function<double(double)>& f2 ,
  const double                         x2 , 
  const std::function<double(double)>& f3 ) 
  : m_edges {{ x1 , x2      }}
  , m_funcs {{ f1 , f2 , f3 }}
{
  Ostap::Assert ( x1 < x2 , s_ERROR , s_METHOD3 ) ;
}
// ============================================================================
/*  Standard constructor with four functions 
 *  @param f1 function to be used for    \f$ x<  x_1\f$
 *  @param x1 the value of \f$x_1\f$ 
 *  @param f2 function to be used for    \f$ x_1 \le  < x_2\f$
 *  @param x2 the value of \f$x_2\f$ 
 *  @param f3 function to be used for    \f$ x_2 \le  < x_3 \f$
 *  @param x3 the value of \f$x_3\f$ 
 *  @param f4 function to be used for    \f$ x_3 \le      \f$
 */
// ============================================================================
Ostap::Math::Piecewise::Piecewise 
( const std::function<double(double)>& f1 ,
  const double                         x1 , 
  const std::function<double(double)>& f2 ,
  const double                         x2 , 
  const std::function<double(double)>& f3 ,
  const double                         x3 , 
  const std::function<double(double)>& f4 ) 
  : m_edges {{ x1 , x2 , x3      }}
  , m_funcs {{ f1 , f2 , f3 , f4 }}
{
  Ostap::Assert ( ( x1 < x2 ) && ( x2 < x3 ) , s_ERROR , s_METHOD4 ) ;
}
// ============================================================================
/*  add new function, defined for  \f$ x\ge x_i\$ 
 *  @attention xi must be larger than any previosly aded ranges!
 *  @param xi value of \f$ x_i \f$
 *  @param fi the function to be used for \f$x\ge x_i\f$ 
 */
// ============================================================================
void Ostap::Math::Piecewise::add 
( const double                         xi , 
  const std::function<double(double)>& fi ) 
{
  //
  Ostap::Assert ( m_edges.empty () || xi > m_edges.back() , s_ERROR , s_METHOD ) ;
  //
  m_edges.push_back ( xi ) ; 
  m_funcs.push_back ( fi ) ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================

