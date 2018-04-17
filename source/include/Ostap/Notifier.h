#ifndef OSTAP_NOTIFIER_H 
#define OSTAP_NOTIFIER_H 1
// ============================================================================
// Include files  
// ============================================================================
// ROOT 
// ============================================================================
#include "TObject.h"
// ============================================================================
// Forward declarationns
// ============================================================================
class TTree ;
// ============================================================================
namespace  Ostap
{
  // ==========================================================================
  namespace  Utils
  {
    // ========================================================================
    /** @class Notifier
     *  Local helper class to keep the proper notifications for TTree
     *  @date 2013-10-13 
     *  @author Vanya BELYAEV Ivan.Brlyaev@itep.ru
     */
    class Notifier : public TObject 
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Utils::Notifier,1) ;
      // ======================================================================
    public:
      // ======================================================================
      Notifier 
      ( TTree*   tree = 0 , 
        TObject* obj0 = 0 , 
        TObject* obj1 = 0 , 
        TObject* obj2 = 0 , 
        TObject* obj3 = 0 , 
        TObject* obj4 = 0 ,
        TObject* obj5 = 0 , 
        TObject* obj6 = 0 , 
        TObject* obj7 = 0 , 
        TObject* obj8 = 0 , 
        TObject* obj9 = 0 ) ;
      // templated constructor 
      template <class ITERATOR>
      Notifier ( ITERATOR  begin ,
                 ITERATOR  end   , 
                 TTree*    tree  ) ;
      /// virtual destructor 
      virtual        ~Notifier () ; // virtual destructor 
      /// the main method 
      virtual Bool_t  Notify   () ;
      // ======================================================================
    public:
      // ======================================================================
      // add object to the notification list 
      inline bool add  ( TObject* o )
      {
        if ( nullptr == o || this == o ) { return false ;}
        m_objects.push_back ( o )  ;
        return true ;
      }
      // exit from  noification 
      bool exit() ;
      // ======================================================================
    private:
      // ======================================================================
      Notifier ( const Notifier & ) ;
      // ======================================================================
    private:
      // ======================================================================
      void _pre_action   () ;
      void _post_action  () ;      
      // ======================================================================
    private:
      // ======================================================================
      TTree*   m_tree ;                  //! the tree 
      TObject* m_old  ;                  //! old notifier  
      // list of fobject to be notified 
      std::vector<TObject*> m_objects ;  //! list of objects 
      // ========================================================================
    } ;
    // ========================================================================
    template <class ITERATOR>
    Notifier::Notifier ( ITERATOR  begin ,
                         ITERATOR  end   , 
                         TTree*    tree  ) 
      : TObject   () 
      , m_tree    ( tree    ) 
      , m_old     ( nullptr )
      , m_objects () 
    {
      this->_pre_action  () ;
      for ( ; begin != end ; ++begin ) { this->add ( *begin ) ; }
      this->_post_action () ;
    }
    // ========================================================================
  } //                                        The end of namespace Ostap::Utils
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                     The END 
// ============================================================================
#endif // OSTAP_NOTIFIER_H
// ============================================================================
