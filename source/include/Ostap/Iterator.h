// ============================================================================
#ifndef ANALYSIS_ITERATOR_H 
#define ANALYSIS_ITERATOR_H 1
// ============================================================================
// Include files
// ============================================================================
//  STD & STL 
// ============================================================================
#include <memory>
// ============================================================================
// Forward deslaration 
// ============================================================================
class TObject          ;  // ROOT 
class TIterator        ;  // ROOT 
class TCollection      ;  // ROOT 
class RooLinkedList    ;  // RooFit
class RooAbsCollection ;  // RooFit
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Utils 
  {
    // ========================================================================
    /** @class Iterator  Ostap/Iterator.h
     *  helper class to simplify manipulations with ROOT-iterator
     *  @see TIterator 
     *  @author Vanay BELYAEV  Ivan.Belyaev@.itep.ru
     *  @date 2015-02-015
     */
    class Iterator 
    {
    public: 
      // ======================================================================
      /// standard constructor: create and keep the ietrator 
      Iterator  ( const RooAbsCollection& collection ) ;    
      /// standard constructor: create and keep the ietrator 
      Iterator  ( const TCollection&      collection ) ;    
      /// standard constructor: create and keep the ietrator 
      Iterator  ( const RooLinkedList&    collection ) ;    
      /// standard constructor: create from existing iterator  
      Iterator  (       TIterator*        iterator   ) ;    
      // ======================================================================
    public:
      // ======================================================================
      /// invoke TIterator::Next
      TObject* next  () const ;                      // invoke TIterator::Next
      /// invoke TIterator::Reset 
      bool     reset () const ;                      // invoke TIterator::Reset 
      // ======================================================================
    public: // aliases
      // ======================================================================
      // alias:
      TObject* Next  () const { return next  () ; }  // invoke TIterator::Next    
      // alias:
      bool     Reset () const { return reset () ; }  // invoke TIterator::Reset
      // ======================================================================
    public: // some pointer alchemistry 
      // ======================================================================
      /// valid iterator 
      bool valid ()             const { return nullptr != m_iterator  ; }
      /// invalid ?
      bool operator!         () const { return !valid () ; }
      /// explicit bool 
      explicit operator bool () const { return  valid () ; }
      /// conversion to underlying type 
      inline TIterator* operator->   () const { return m_iterator.get() ; }
      // ======================================================================
    public:
      // =========================================================================    
      /** iterate and cast to specific type
       */ 
      template <class T>
      T* dynamic_next() const 
      {
        TObject* r = this->Next() ;
        return ( nullptr == r ) ? nullptr : dynamic_cast<T*> ( r ) ; // RETURN
      }
      // ========================================================================
      /** iterate and cast to specific type
       */ 
      template <class T>
      T* static_next() const 
      {
        TObject* r = this->Next() ;
        return ( nullptr == r ) ? nullptr : static_cast<T*> ( r ) ; // RETURN
      }
      // =========================================================================    
    private:
      // ======================================================================
      /// copy constructor is disabled 
      Iterator           ( const Iterator& ) ; // copy constructor is disabled 
      /// assignement is disabled 
      Iterator& operator=( const Iterator& ) ; // assignement is disabled
      // ======================================================================    
    private:
      // ======================================================================
      /// iterator itself 
      std::unique_ptr<TIterator> m_iterator ; // iterator itself 
      // ======================================================================    
    };
    // ========================================================================
  } //                                            end of namespace Ostap::Utils 
  // ==========================================================================
}
// ============================================================================
//                                                                      The END 
// ===========================================================================
#endif // ANALYSIS_ITERATOR_H
// ============================================================================
