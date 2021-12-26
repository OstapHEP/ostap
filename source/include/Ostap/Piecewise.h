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
      typedef std::pair< std::function<double(double)>,double> FPAIR ;
      // ======================================================================
    public: 
      // ======================================================================
      /** Standard constructor with a single function 
       *  @param f1 the functon to be used for \f$ -\infty < x < +\infty \f$ 
       *  @param s1 the scale factor
       */
      template <class FUNCTION> 
      Piecewise
      ( FUNCTION     f1     ,
        const double s1 = 1 ) 
        : m_edges () 
        , m_funcs ( 1 , FPAIR ( f1 , s1 ) )
      {}
      // =======================================================================
      /// constructor with several functions 
      template <class FUNCTION, typename... ARGS> 
      Piecewise
      ( FUNCTION     f1   ,
        const double s1   , 
        ARGS...      args ) 
        : m_edges () 
        , m_funcs ( 1 , FPAIR ( f1 , s1 ) )
      {
        this->add ( args... ) ;
      }
      // ======================================================================
      /// default constructor : create a constant function
      Piecewise ( const double value = 0 ) ;
      // ======================================================================
      Piecewise ( const Piecewise&  ) = default ;
      Piecewise (       Piecewise&& ) = default ;
      // ======================================================================
    public: // the main method 
      // ======================================================================
      /// the main method 
      inline double operator() ( const double x ) const 
      { 
        const FPAIR& p = m_funcs [ index( x  ) ] ;
        return p.first ( x ) * p.second ;
      }
      // ======================================================================
    public:
      // ======================================================================
      /// get all edges 
      const std::vector<double>& edges    () const 
      { return m_edges ; }
      /// get all (function,scale) pairs 
      const std::vector<FPAIR> functions() const { return m_funcs ; }
      /// find the proper interval 
      std::size_t index (  const double x ) const 
      {
        auto found = std::upper_bound ( m_edges.begin () , m_edges.end () , x ) ;
        return m_edges.end() == found ? m_edges.size () : found - m_edges.begin() ;
      }
      // ======================================================================
    protected:
      // ======================================================================
      /** add (function,scale) pair into the list into the list 
       *  - check that x is larger than all previous values
       */
      void add_
      ( const double                  x     ,  
        std::function<double(double)> f     , 
        const double                  s = 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// template creator 
      template <class FUNCTION, typename... ARGS>
      static inline Piecewise
      create
      ( FUNCTION     f1   ,
        const double s1   , 
        ARGS...      args ) { return Piecewise ( f1 , s1 , args... ) ; }  
      // ======================================================================
    public:
      // ======================================================================
      /** add a new function, defined for  \f$ x\ge x_i\$ 
       *  @attention xi must be larger than any previosly added ranges!
       *  @param xi value of \f$ x_i \f$
       *  @param fi the function to be used for \f$x\ge x_i\f$ 
       *  @param si scale factor to be used 
       */
      template <class FUNCTION> 
      void add
      ( const double  xi       , 
        FUNCTION      fi       , 
        const double  si = 1.0 ) { this->add_ ( xi , fi , si ) ; }
      // ======================================================================
      /** add new function, defined for  \f$ x\ge x_i\$ 
       *  @attention xi must be larger than any previosly added ranges!
       *  @param fi the function to be used for \f$x\ge x_i\f$ 
       *  @param xi value of \f$ x_i \f$
       *  @param si scale efactor to be used 
       */
      template <class FUNCTION, typename... ARGS>
      void add
      ( const double xi   , 
        FUNCTION     fi   ,
        const double si   , 
        ARGS...      args )
      {
        this->add_ ( xi , fi , si ) ;
        this->add  ( args...      ) ;
      }
      // ======================================================================
    public:
      // ======================================================================
      /// scale the function 
      Piecewise& operator*= ( const double value ) ;
      /// scale the function 
      Piecewise& operator/= ( const double value ) ;
      // ======================================================================
    private:
      // ====================================================================== 
      /// list of edges 
      std::vector<double> m_edges ; // list of edges 
      /// (function,scale) pairs 
      std::vector<FPAIR>  m_funcs ; // (function,scale) pairs 
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
