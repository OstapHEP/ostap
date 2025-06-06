// ============================================================================
#ifndef OSTAP_PROGRESSBAR_H 
#define OSTAP_PROGRESSBAR_H 1
// ============================================================================
// STD&STL
// ============================================================================
#include <string>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/ProgressConf.h"
// ============================================================================
/** @file Ostap/Utils.h
 *  collection of various C++ utilities 
 *  @author Vanya Belyaev
 *  @date   2018-03-23
 */
namespace Ostap
{
  // ==========================================================================
  namespace Utils 
  {
    // ========================================================================
    /** @class ProgressBar
     *  simple pgrepgress bar
     */
    class ProgressBar : public ProgressConf 
    {
      // ======================================================================
    public:
      // ======================================================================
      /** Constructor from configuration and maximal conut 
       *  @param maxcount   maximal count (==0 : no display bar)
       *  @param conf       configuration 
       */
      ProgressBar 
      ( const ProgressConf&      conf          , 
        const unsigned long long maxcount  = 0 ) ;
      // ======================================================================
      /** Constructor from configuration and maximal conut 
       *  @param conf       configuration 
       *  @param maxcount   maximal count (==0 : no display bar)
       */
      ProgressBar 
      ( const unsigned long long maxcount ,  
        const ProgressConf&      conf     ) ;
      // ======================================================================
      /** full constructor for print to std::cout
       *  @param maxcount   maximal count (==0 : no display bar)
       *  @param width      the width     (==0: no display bar)
       *  @param symbol     symbol to show 
       *  @param empty      empty symbol
       *  @param left       left/prefix 
       *  @param right      right/prefix 
       *  @param what       description/prefix 
       *  @param use_timer  use the timer?
       *  @param atty       isatty ? 
       */
      ProgressBar 
      ( const unsigned long long maxcount  = 0      ,
        const unsigned short     width     = 80     ,
        const std::string&       symbol    = "#"    ,
        const std::string&       empty     = " "    ,
        const std::string&       left      = "[ "   ,
        const std::string&       right     = " ]"   ,
        const std::string&       what      = ""     ,
        const bool               timer     = true   ,
        const bool               atty      = true   ) ;
      //
      // ======================================================================
      // destructor
      ~ProgressBar () ;
      // ======================================================================
    public:
      // ======================================================================
      /// current count 
      unsigned long long count     () const { return m_count    ; } 
      /// maximal count 
      unsigned long long max_count () const { return m_maxcount ; }
      // ======================================================================
    public:
      // ======================================================================
      /// Is this progress bar enabled? 
      inline bool enabled   () const { return 0 != m_maxcount && 0 != width () ;  }
      /// Is this progress bar disabled? 
      inline bool disabled  () const { return 0 == m_maxcount || 0 == width () ; }
      // =====================================================================
    public:
      // ======================================================================
      inline ProgressBar& operator++() { return operator+= ( 1 ) ; }
      inline ProgressBar& operator+=( const unsigned int increment ) 
      {
        m_count += increment ;
        return
          !increment                    ? *this : 
          // show nothing 
          disabled ()                   ? *this :
          // last count?
          m_maxcount   == m_count ? show_bar () :
          // right moment to show ? 
          m_next_count <= m_count ? show_bar () :
          // explicitly show the first five counts
          100 <= m_maxcount && m_count <= 5  ? show_bar () : 
          // explicitely show the last five counts 
          100 <= m_maxcount && m_count <= m_maxcount && ( m_maxcount <= m_count + 5 ) ? show_bar () 
          : *this ;
      }
      // ======================================================================
    private:
      // ======================================================================
      /// show the bar
      ProgressBar& show_bar ( const bool show_eta = true ) ;
      // ======================================================================
    private:
      // ======================================================================
      /// maximal count 
      unsigned long long m_maxcount   { 0  }   ;
      /// current count 
      unsigned long long m_count      { 0  }   ;
      /// next counte 
      unsigned long long m_next_count { 0  }   ;
      /// total width of the line
      unsigned int       m_wtot       { 80 }   ;
      // start time 
      unsigned long long m_start      { 0  }   ;
      // ======================================================================
    } ;
    // ========================================================================
  } //                                        The edn of namespace Ostap::Utils 
  // ==========================================================================
} //                                                The  end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_PROGRESSBAR_H
// ============================================================================
