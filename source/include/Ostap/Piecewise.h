// ============================================================================
#ifndef OSTAP_PIECEWISE_H 
#define OSTAP_PIECEWISE_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <functional>
#include <algorithm>
#include <vector>
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    /** @class Piecewise Ostap/Piecewise.h
     *  Simple Piecewise-function
     *  @author Vanya Belyaev
     *  @date   2020-06-29
     */
    class Piecewise 
    {
    public: 
      // ======================================================================
      /// Standard constructor with a single function 
      Piecewise ( const  std::function<double(double)>& f1 ) ;
      // ======================================================================
      /** Standard constructor with two functions 
       *  @param f1 function to be used for    \f$ x<   x_1\f$
       *  @param x1 the value of \f$x_1\f$ 
       *  @param f2 function to be used for    \f$ x\ge x_1\f$
       */
      Piecewise ( const std::function<double(double)>& f1 ,
                  const double                         x1 , 
                  const std::function<double(double)>& f2 ) ;
      // ======================================================================
      /** Standard constructor with three functions 
       *  @param f1 function to be used for    \f$ x<  x_1\f$
       *  @param x1 the value of \f$x_1\f$ 
       *  @param f2 function to be used for    \f$ x_1 \le  < x_2\f$
       *  @param x2 the value of \f$x_2\f$ 
       *  @param f3 function to be used for    \f$ x_2 \le      \f$
       */
      Piecewise ( const std::function<double(double)>& f1 ,
                  const double                         x1 , 
                  const std::function<double(double)>& f2 ,
                  const double                         x2 , 
                  const std::function<double(double)>& f3 ) ;
      // ======================================================================
      /** Standard constructor with four functions 
       *  @param f1 function to be used for    \f$ x<  x_1\f$
       *  @param x1 the value of \f$x_1\f$ 
       *  @param f2 function to be used for    \f$ x_1 \le  < x_2\f$
       *  @param x2 the value of \f$x_2\f$ 
       *  @param f3 function to be used for    \f$ x_2 \le  < x_3 \f$
       *  @param x3 the value of \f$x_3\f$ 
       *  @param f4 function to be used for    \f$ x_3 \le      \f$
       */
      Piecewise ( const std::function<double(double)>& f1 ,
                  const double                         x1 , 
                  const std::function<double(double)>& f2 ,
                  const double                         x2 , 
                  const std::function<double(double)>& f3 ,
                  const double                         x3 , 
                  const std::function<double(double)>& f4 ) ;
      // ======================================================================
    public: // the main method 
      // ======================================================================
      /// the main method 
      double operator() ( const double x ) const 
      { return m_funcs [ index ( x ) ] ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get all edges 
      const std::vector<double>&                        edges    () const 
      { return m_edges ; }
      /// get all functgions 
      const std::vector<std::function<double(double)> > functions() const 
      { return m_funcs ; }
      /// find the proper interval 
      std::size_t index (  const double x ) const 
      {
        auto found = std::upper_bound ( m_edges.begin() , m_edges.end()  , x ) ;
        return m_edges.end() == found ? m_edges.size () : found - m_edges.begin() ;
      }
      // ======================================================================
    public:
      // ======================================================================
      /** add new function, defined for  \f$ x\ge x_i\$ 
       *  @attention xi must be larger than any previosly added ranges!
       *  @param xi value of \f$ x_i \f$
       *  @param fi the function to be used for \f$x\ge x_i\f$ 
       */
      void add ( const double                         xi , 
                 const std::function<double(double)>& fi ) ;
      // ======================================================================
      /** add new function, defined for  \f$ x\ge x_i\$ 
       *  @attention xi must be larger than any previosly added ranges!
       *  @param fi the function to be used for \f$x\ge x_i\f$ 
       *  @param xi value of \f$ x_i \f$
       */
      void add ( const std::function<double(double)>& fi , 
                 const double                         xi ) { add ( xi , fi ) ; }
      // ======================================================================
    private:
      // ======================================================================
      std::vector<double>                         m_edges ;
      std::vector<std::function<double(double)> > m_funcs ;
      // ======================================================================
    };
    // ========================================================================
  } //                                         The end of namespace Ostap::Math  
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_PIECEWISE_H
// ============================================================================
