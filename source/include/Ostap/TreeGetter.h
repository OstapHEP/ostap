// ============================================================================
#ifndef OSTAP_TREEGETTER_H 
#define OSTAP_TREEGETTER_H 1
// ============================================================================
// Include files 
// ============================================================================
/// STD&STL
// ============================================================================
#include <memory>
#include <vector>
#include <string>
#include <map>
// ============================================================================
// ROOT 
// ============================================================================
#include "TObject.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Notifier.h"
#include "Ostap/Formula.h"
#include "Ostap/StatusCode.h"
// ============================================================================
// Forward declarations
// ============================================================================
class TTree ; // ROOT 
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Trees
  { 
    // ========================================================================
    /** @class Getter 
     *  Helper class to get value for the certains expressions 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2021-09-22
     */
    class Getter : public TObject 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Trees::Getter,1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the tree and expression 
      Getter
      ( const std::string& expression               ,
        TTree*             tree        = nullptr    ) ;
      /// constructor from the tree and expression 
      Getter
      ( const std::string& expression1              ,
        const std::string& expression2              ,
        TTree*             tree        = nullptr    ) ;
      /// constructor from the tree and expression 
      Getter ( const std::string& expression1              ,
               const std::string& expression2              ,
               const std::string& expression3              ,
               TTree*             tree        = nullptr    ) ;
      Getter ( const std::string& expression1              ,
               const std::string& expression2              ,
               const std::string& expression3              ,
               const std::string& expression4              ,
               TTree*             tree        = nullptr    ) ;
      Getter ( const std::string& expression1              ,
               const std::string& expression2              ,
               const std::string& expression3              ,
               const std::string& expression4              ,
               const std::string& expression5              ,
               TTree*             tree        = nullptr    ) ;
      /// constuctor from the list of expressions and TTree 
      Getter
      ( const std::vector<std::string>& expressions ,
        TTree* = nullptr                            ) ; 
      /// constuctor from the list of expressions and TTree 
      Getter
      ( const std::map<std::string,std::string>& expressions ,
        TTree* = nullptr                            ) ; 
      // ======================================================================
      /// copy constructor 
      Getter ( const Getter&  ) = delete ;
      // ======================================================================
    public:
      // ======================================================================
      /// actually we delete all formulae instead of notification.
      Bool_t Notify   () override ; 
      // ======================================================================
    public:
      // ======================================================================
      /// get the tree 
      TTree* tree () const { return m_tree ; }
      // ======================================================================
      /// get the results 
      Ostap::StatusCode  eval
      ( std::vector<double>& result           ,
        TTree*               tree   = nullptr ) const ;
      /// =====================================================================
      /// get the results 
      Ostap::StatusCode  eval
      ( std::map<std::string,double>& result           ,
        TTree*                        tree   = nullptr ) const ;
      // ======================================================================
      /// get the results
      inline 
      Ostap::StatusCode  eval
      ( TTree&                        tree   , 
        std::vector<double>&          result ) const { return eval ( result , &tree ) ; }
      // ======================================================================
      /// get the results
      inline 
      Ostap::StatusCode  eval
      ( TTree&                        tree   , 
        std::map<std::string,double>& result ) const { return eval ( result , &tree ) ; }
      // ======================================================================
    protected : 
      // ======================================================================
      /// make formulae 
      Ostap::StatusCode make_formulae () const ; // make formulae
      // is everything OK ? 
      Ostap::StatusCode ok            ( TTree* tree = nullptr ) const ; 
      // ======================================================================
    private : 
      // ======================================================================
      typedef std::map<std::string,std::string>                      SMAP ;
      typedef std::map<std::string,std::unique_ptr<Ostap::Formula> > FMAP ;
      typedef SMAP::const_iterator                                   SIT  ;
      typedef FMAP::const_iterator                                   FIT  ;
      /// The tree 
      mutable TTree* m_tree     { nullptr } ;
      // list of formulae
      SMAP           m_map      {} ;
      // list of formulae
      mutable FMAP   m_formulae {} ;         
      // ======================================================================
    } ; //                                 The end of class Ostap::Trees:Getter 
    // ========================================================================
  } //                                        The end of namespace Ostap::Trees 
  // ==========================================================================
} //                                                 The end of namesapce Ostap 
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_TREEGETTER_H
// ============================================================================

