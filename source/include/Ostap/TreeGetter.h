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
class TTree            ; // ROOT
class RooAbsCollection ; // ROOT/Roofit  
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
      typedef std::map<std::string,double>                               RMAP ;      
      typedef std::vector<double>                                        RVCT ;      
      typedef std::map<std::string,std::string>                          DCT  ;
      // ======================================================================
    public: 
      // ======================================================================
      ClassDefOverride(Ostap::Trees::Getter,1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constuctor from the list of expressions and TTree 
      Getter
      ( const std::vector<std::string>&          expressions           ,
        const TTree*                             tree        = nullptr ) ; 
      /// constuctor from the map of expressions and TTree 
      Getter
      ( const std::map<std::string,std::string>& expressions           ,
        const TTree*                             tree        = nullptr ) ; 
      // ======================================================================
      /// copy constructor 
      Getter ( const Getter&   ) ;
      /// move  constructor 
      Getter (       Getter&&  ) = default ;
      /// clone method 
      Getter* Clone ( const char* newname = "" ) const override ;
      // ======================================================================
      /// virtual destructor 
      virtual ~Getter() ;
      // ======================================================================      
    public:
      // ======================================================================
      /// get the mapping
      const std::map<std::string,std::string>& mapping () const { return m_map ; }
      /** add entry mapping 
       *  @attention no replacement!!! 
       */
      bool add
      ( const std::string& item            ,
        const std::string& expression = "" ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tree 
      const TTree* tree () const { return m_tree ; }
      // ======================================================================
      /// get the results 
      Ostap::StatusCode  eval
      ( RVCT&             result           ,
        const TTree*      tree   = nullptr ) const ;
      /// =====================================================================
      /// get the results 
      Ostap::StatusCode  eval
      ( RMAP&             result           ,
        const TTree*      tree   = nullptr ) const ;
      // ======================================================================
      /// get the results
      inline 
      Ostap::StatusCode  eval
      ( const TTree&      tree   , 
        RVCT&             result ) const { return eval ( result , &tree ) ; }
      // ======================================================================
      /// get the results
      inline 
      Ostap::StatusCode  eval
      ( const TTree&      tree   , 
        RMAP&             result ) const { return eval ( result , &tree ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// actually we delete all formulae instead of notification.
      Bool_t Notify   () override ; 
      // ======================================================================
    protected : 
      // ======================================================================
      /// make formulae 
      Ostap::StatusCode make_formulae () const ; // make formulae
      // is everything OK ? 
      Ostap::StatusCode ok            ( const TTree* tree = nullptr ) const ; 
      // ======================================================================
    protected : 
      // ======================================================================
      typedef std::map<std::string,std::string>                      SMAP ;
      typedef std::map<std::string,std::unique_ptr<Ostap::Formula> > FMAP ;
      typedef SMAP::const_iterator                                   SIT  ;
      typedef FMAP::const_iterator                                   FIT  ;
      /// The tree 
      mutable const TTree* m_tree     { nullptr } ;
      // list of formulae
      SMAP                 m_map      {} ;
      // list of formulae
      mutable FMAP         m_formulae {} ;         
      // ======================================================================
    } ; //                                 The end of class Ostap::Trees:Getter
    // ========================================================================
    /** @class RooGetter 
     *  Helper class to get value for the certains expressions 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2021-09-22
     */
    class RooGetter : public Getter
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Trees::RooGetter,1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constuctor from the list of expressions and TTree 
      RooGetter
      ( const std::vector<std::string>&          expressions           ,
        const TTree*                             tree        = nullptr ) ; 
      /// constuctor from the map of expressions and TTree 
      RooGetter
      ( const std::map<std::string,std::string>& expressions           ,
        const TTree*                             tree        = nullptr ) ; 
      // ======================================================================
      /// copy constructor 
      RooGetter ( const RooGetter&   ) ;
      /// move  constructor 
      RooGetter (       RooGetter&&  ) = default ;
      /// clone method 
      RooGetter* Clone ( const char* newname = "" ) const override ;
      // ======================================================================
      /// virtual destructor 
      virtual ~RooGetter() ;
      // ======================================================================      
    public:
      // ======================================================================
      /// get the results 
      Ostap::StatusCode  assign
      ( RooAbsCollection& result           ,
        const TTree*      tree   = nullptr ) const ;
      /// get the results
      Ostap::StatusCode  assign
      ( const TTree&      tree             , 
        RooAbsCollection& result           ) const ;
      // ======================================================================
    } ; //                              The end of class Ostap::Trees:RooGetter    
    // ========================================================================
  } //                                        The end of namespace Ostap::Trees 
  // ==========================================================================
} //                                                 The end of namesapce Ostap 
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_TREEGETTER_H
// ============================================================================

