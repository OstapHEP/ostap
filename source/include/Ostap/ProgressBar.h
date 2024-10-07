// ============================================================================
#ifndef OSTAP_PROGRESSBAR_H 
#define OSTAP_PROGRESSBAR_H 1
// ============================================================================
// STD&STL
// ============================================================================
#include <string>
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
    /** @class ProgressConf
     *  configurtation for the progressbar
     *  @see Ostap::Utils::ProgressBar 
     */
    class ProgressConf
    {
      // ======================================================================
    public:
      // ======================================================================
      /** full constructor
       *  @param width      the width  (zero width disables the progress bar!)  
       *  @param symbol     symbol to show as "done" 
       *  @param empty      symbol to show as "not-yet"
       *  @param left       left/prefix 
       *  @param right      right/suffix 
       *  @param what       description/prefix 
       *  @param use_timer  use the timer?
       */
      ProgressConf
      ( const unsigned short     width     = 80     ,        
        const std::string&       symbol    = "#"    ,
        const std::string&       empty     = " "    ,
        const std::string&       left      = "[ "   ,
        const std::string&       right     = " ]"   ,
        const std::string&       what      = ""     , 
        const bool               use_timer = true   ) ;
      // ======================================================================
    public: // ghetters 
      // ======================================================================
      /// "done" symbol 
      const std::string& symbol    () const { return m_symbol    ; }
      ///  "not-yet" symbol 
      const std::string& empty     () const { return m_empty     ; }
      /// left prefix 
      const std::string& left      () const { return m_left      ; }
      /// rigth suffix 
      const std::string& right     () const { return m_right     ; }
      /// description
      const std::string& what      () const { return m_what      ; }
      /// effective width/lngth  of the bar (number of symbols+empty steps )
      unsigned int       width     () const { return m_width     ; }
      /// use the timer (show ETA) ?
      bool               use_timer () const { return m_use_timer ; }
      // ======================================================================
    public : // setters 
      // ======================================================================
      /// symbol for "done"
      void setSymbol   ( const std::string&   value ) { m_symbol = value ; }
      /// symbol for "not-yet"
      void setEmpty    ( const std::string&   value ) { m_empty  = value ; }
      /// left 
      void setLeft     ( const std::string&   value ) { m_left   = value ; }
      /// right 
      void setRight    ( const std::string&   value ) { m_right  = value ; }
      /// width 
      void setWidth    ( const unsigned short value ) ; 
      /// use timer ? 
      void setUseTimer ( const bool           value ) ;
      // ======================================================================
    private:
      // ======================================================================
      /// the width (zero width disables the progress bar!)
      unsigned int       m_width      { 80   } ;
      /// symbol
      std::string        m_symbol     { "#"  } ;
      /// empty symbol
      std::string        m_empty      { " "  } ;
      /// symbols 
      std::string        m_left       { "[ " } ;
      /// symbols 
      std::string        m_right      { " ]" } ;
      /// symbols 
      std::string        m_what       { ""   } ;
      /// use timer ?
      bool               m_use_timer  { true } ;
      // ======================================================================
    } ;
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
       */
      ProgressBar 
      ( const unsigned long long maxcount  = 0      ,
        const unsigned short     width     = 80     ,
        const std::string&       symbol    = "#"    ,
        const std::string&       empty     = " "    ,
        const std::string&       left      = "[ "   ,
        const std::string&       right     = " ]"   ,
        const std::string&       what      = ""     ,
        const bool               timer     = true   ) ;
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
          // rigth moment to show ? 
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
