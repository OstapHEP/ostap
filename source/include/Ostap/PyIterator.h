// ============================================================================
#ifndef OSTAP_PYITERATOR_H 
#define OSTAP_PYITERATOR_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <limits>
#include <memory>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/ProgressBar.h"
// ============================================================================
// forward declarations 
// ============================================================================
class TTree ; // from ROOT 
class TCut  ;
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  class Formula ;
  // ==========================================================================
  /** @class PyIterator Ostap/PyIterator.h
   *  Helper class for fast iterator over TTree in python.
   *  Iteration in python for large tree can be very time consuming, 
   *  for such cases this "iterator-with-cuts" is much faster
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date   2013-05-06
   */
  class PyIterator 
  {
  public:
    // ========================================================================
    /// constructor 
    PyIterator 
    ( TTree*              tree      , 
      const std::string&  cuts      , 
      const unsigned long first = 0 , 
      const unsigned long last  = std::numeric_limits<unsigned long>::max() ) ;
    //
    PyIterator 
    ( TTree*              tree      , 
      const TCut&         cuts      , 
      const unsigned long first = 0 , 
      const unsigned long last  = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /// constructor 
    PyIterator 
    ( TTree*                            tree      , 
      const Ostap::Utils::ProgressConf& progress  ,
      const std::string&                cuts      , 
      const unsigned long               first = 0 , 
      const unsigned long last  = std::numeric_limits<unsigned long>::max() ) ;
    //
    PyIterator 
    ( TTree*                            tree      , 
      const Ostap::Utils::ProgressConf& progress  ,
      const TCut&                       cuts      , 
      const unsigned long               first = 0 , 
      const unsigned long last  = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
  public:
    // ========================================================================
    /// go to next item 
    TTree* next    () const ;                      // go to next item
    /// get the tree 
    TTree* tree    () const { return m_tree   ; }  // get the tree
    /// get the weight
    double weight  () const { return m_weight ; }  // the weight 
    /// check if formula is ok 
    bool   ok      () const ;
    /// get formula 
    const Ostap::Formula* formula() const { return m_formula.get()  ; }    
    // ========================================================================
  private:
    // ========================================================================
    PyIterator () ;
    PyIterator           ( const PyIterator& ) ;
    PyIterator& operator=( const PyIterator& ) ;    
    /// intialize cuts 
    void init  ( const std::string& cuts ) ;  // intialize cuts 
    // ========================================================================
  private:
    // ========================================================================
    TTree*                            m_tree     { nullptr } ;
    std::unique_ptr<Ostap::Formula>   m_formula  {   }       ;
    unsigned long                     m_last     { 0 }       ;
    mutable double                    m_weight   { 1 }       ; // weight 
    mutable unsigned long             m_current  { 0 }       ;
    mutable Ostap::Utils::ProgressBar m_progress {   }       ;
    // ========================================================================
  };
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_PYITERATOR_H
// ============================================================================
