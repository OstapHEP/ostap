#ifndef OSTAP_HISTODUMP_H 
#define OSTAP_HISTODUMP_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL 
// ============================================================================
#include <string>
// ============================================================================
/// forward declarations:
// ============================================================================
class TH1      ;                                                        // ROOT 
class TProfile ;                                                        // ROOT 
// ============================================================================
/** @file Ostap/HistoDump.h
 *  dump the histogtams 
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Utils 
  {
    // ========================================================================
    namespace Histos 
    {
      // ======================================================================
      /** dump the text representation of the Profile 
       *  @param histo  (INPUT) the histogram 
       *  @param stream (OUTUT) the stream  
       *  @param width  (INPUT) the maximal column width 
       *  @param height (INPUT) the proposed coulmn height 
       *  @return the stream 
       *  @author Vanya BELYAEV  Ivan.BElyaev@itep.ru
       *  @date 2009-09-19
       */ 
      std::ostream& histoDump_
      ( const TProfile*           histo           , 
        std::ostream&             stream          ,
        const std::size_t         width   = 80    , 
        const std::size_t         height  = 50    ) ;
      // ====================================================================
      /** dump the text representation of the histogram 
       *  @param histo  (INPUT) the histogram 
       *  @param width  (INPUT) the maximal column width 
       *  @param height (INPUT) the propsoed coulmn height 
       *  @return string representation of the histogram       
       *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
       *  @date 2009-09-19
       */ 
      std::string   histoDump 
      ( const TProfile*           histo          , 
        const std::size_t         width  = 80    , 
        const std::size_t         height = 50    ) ;
      // ====================================================================
      /** dump the text representation of the histogram 
       *  @param histo  (INPUT) the histogram 
       *  @param stream (OUTUT) the stream  
       *  @param width  (INPUT) the maximal column width 
       *  @param height (INPUT) the proposed coulmn height 
       *  @param errors (INPUT) print/plot errors
       *  @return the stream 
       *  @author Vanya BELYAEV  Ivan.BElyaev@itep.ru
       *  @date 2009-09-19
       */ 
      std::ostream& histoDump_
      ( const TH1*                histo           , 
        std::ostream&             stream          ,
        const std::size_t         width   = 80    , 
        const std::size_t         height  = 50    ,  
        const bool                errors  = false ) ;
      // ====================================================================
      /** dump the text representation of the histogram 
       *  @param histo  (INPUT) the histogram 
       *  @param width  (INPUT) the maximal column width 
       *  @param height (INPUT) the propsoed coulmn height 
       *  @param errors (INPUT) print/plot errors
       *  @return string representation of the histogram       
       *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
       *  @date 2009-09-19
       */ 
      std::string   histoDump 
      ( const TH1*                histo          , 
        const std::size_t         width  = 80    , 
        const std::size_t         height = 50    ,  
        const bool                errors = false ) ;
      // ======================================================================
    } //                                  end of namespace Ostap::Utils::Histos 
    // ========================================================================
  } //                                            end of namespace Ostap::Utils 
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_HISTODUMP_H
// ============================================================================
