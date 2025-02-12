// ============================================================================
#ifndef OSTAP_PROGRESSCONF_H 
#define OSTAP_PROGRESSCONF_H 1
// ============================================================================
// STD&STL
// ============================================================================
#include <string>
// ============================================================================
/** @file Ostap/ProgressConf.h
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
       *  @param atty       isatty ? 
       */
      ProgressConf
      ( const unsigned short     width     = 80     ,        
        const std::string&       symbol    = "#"    ,
        const std::string&       empty     = " "    ,
        const std::string&       left      = "[ "   ,
        const std::string&       right     = " ]"   ,
        const std::string&       what      = ""     , 
        const bool               use_timer = true   ,
        const bool               atty      = true   ) ;
      // ======================================================================
      /** onstructor
       *  @param show show default progress bar? 
       */
      // ProgressConf ( const bool show = true ) ;         
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
      /// isatty?
      bool               atty      () const { return m_atty      ; }
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
      /// isatty ?
      bool               m_atty       { true } ;
      // ======================================================================
    } ;
    // ========================================================================
  } //                                        The edn of namespace Ostap::Utils 
  // ==========================================================================
} //                                                The  end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_PROGRESSCONF_H
// ============================================================================
