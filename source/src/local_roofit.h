// ============================================================================
#ifndef LOCAL_ROOFIT_H 
#define LOCAL_ROOFIT_H 1
// ============================================================================
// Include files
// ============================================================================
// ROOT&RooFit
// ============================================================================
#include "RVersion.h"
#include "RooArgList.h"
#include "RooListProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Iterator.h"
// ============================================================================
// Local
// ============================================================================
#include "Exception.h"
// ============================================================================
/** @file local_roofit.h
 *  Collection of simple utilities for RooFit objects 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2020-03-06
 */
// ============================================================================
namespace
{
  // ==========================================================================
  /// size of RooArgList 
  inline std::size_t size ( const RooAbsCollection& lst ) 
  {
#if ROOT_VERSION_CODE < ROOT_VERSION(6,18,0)
    return lst.getSize () ;
#else 
    return lst.size    () ;
#endif
  }
  // ==========================================================================
  /** copy from list to list
   *  @param from objects to be copied from this list 
   *  @param to   objects to be copied to this proxy 
   */
  inline void copy 
  ( const RooAbsCollection&  from ,
    RooAbsCollection&        to   ) 
  {
    //
#if ROOT_VERSION_CODE < ROOT_VERSION(6,18,0)
    //
    Ostap::Utils::Iterator tmp ( from ) ; // only for ROOT < 6.18
    RooAbsArg* c = 0 ;
    while ( c = (RooAbsArg*) tmp.next() ) { to.add ( *c ) ; }
    //
#else
    //
    for ( auto* c : from ) { to.add ( *c ) ; }
    //
#endif 
    //
  }
  // ==========================================================================
  /** copy RooAbsReal from list to list
   *  @param from objects to be copied from this list 
   *  @param to   objects to be copied to this proxy 
   */
  inline unsigned   int copy_real 
  ( const RooAbsCollection&  from                                    ,
    RooAbsCollection&        to                                      , 
    const std::string&       message = "Variable is not RooAbsReal!" ,
    const std::string&       tag     = "Ostap::copy_real"            )
  {
    //
#if ROOT_VERSION_CODE < ROOT_VERSION(6,18,0)
    //
    Ostap::Utils::Iterator tmp ( from ) ; // only for ROOT < 6.18 
    RooAbsArg* c = 0 ;
    while ( c = (RooAbsArg*) tmp.next() )
    {
      Ostap::Assert ( dynamic_cast<RooAbsReal*> ( c ) != nullptr , message , tag , 510 ) ;
      to.add ( *c ) ;
    }
    //
#else
    //
    unsigned ii = 0 ;
    for ( auto* c : from ) 
    {
      Ostap::Assert ( dynamic_cast<RooAbsReal*> ( c ) != nullptr , message , tag , 510 ) ;
      to.add ( *c ) ;   
      ++ii ;
    }
    //
#endif 
    //
    return ::size ( from ) ;
    //
  }
  // ==========================================================================
  /** copy RooAbsReal from list to list
   *  @param from objects to be copied from this list 
   *  @param to   objects to be copied to this proxy 
   */
  inline unsigned   int copy_real 
  ( const RooAbsCollection&  from                                    ,
    RooArgSet&               to                                      , 
    const std::string&       message = "Variable is not RooAbsReal!" ,
    const std::string&       tag     = "Ostap::copy_real"            )
  {
    //
#if ROOT_VERSION_CODE < ROOT_VERSION(6,18,0)
    //
    Ostap::Utils::Iterator tmp ( from ) ; // only for ROOT < 6.18 
    RooAbsArg* c = 0 ;
    while ( c = (RooAbsArg*) tmp.next() )
    {
      Ostap::Assert ( dynamic_cast<RooAbsReal*> ( c ) != nullptr , message , tag , 510 ) ;
      to.add ( *c ) ;
    }
    //
#else
    //
    unsigned ii = 0 ;
    for ( auto* c : from ) 
    {
      Ostap::Assert ( dynamic_cast<RooAbsReal*> ( c ) != nullptr , message , tag , 510 ) ;
      to.add ( *c ) ;   
      ++ii ;
    }
    //
#endif 
    //
    return ::size ( from ) ;
    //
  }
  // ==========================================================================
  /// get parameter from RooListProxy
  inline double get_par 
  ( const unsigned short index , 
    const RooListProxy&  lst   ) 
  {
    const RooAbsArg* v    = lst.at ( index ) ;
    if ( nullptr ==  v ) { return 0 ; }
    const RooArgSet* nset = lst.nset() ;
    //
    const RooAbsReal* r = static_cast<const RooAbsReal*>( v ) ;
    return r->getVal ( nset ) ;
  }
  // ==========================================================================
  template <class OBJECT>
  inline void set_pars 
  ( const RooListProxy&  lst      , 
    OBJECT&              obj      , 
    const unsigned short shift = 0 ) 
  {
    //
    const RooArgSet* nset  = lst.nset() ;
    //
#if ROOT_VERSION_CODE < ROOT_VERSION(6,18,0)
    //
    Ostap::Utils::Iterator it  ( lst  ) ; // only for ROOT < 6.18 
    RooAbsArg*   p = 0 ;
    unsigned int k = 0 ;
    while ( ( p = (RooAbsArg*) it.next() ) )
    {
      const RooAbsReal* r = static_cast<const RooAbsReal*>( p ) ;
      obj.setPar ( k + shift , r->getVal ( nset ) ) ;  
      ++k ;
    }
    //
#else
    //
    const unsigned int N = lst.size() ;
    for ( unsigned int k = 0 ; k < N ; ++k ) 
    {
      const RooAbsReal& r = static_cast<const RooAbsReal&>( lst [ k ] ) ;
      obj.setPar ( k + shift , r.getVal ( nset ) ) ;  
    }
    //
#endif
    //
  }
  // =======================================================================================
  inline void set_pars 
  ( const RooListProxy&  lst   , 
    std::vector<double>& vct   )  
  {
    //
    const RooArgSet* nset  = lst.nset() ;
    vct.resize ( ::size ( lst ) ) ;
    //
#if ROOT_VERSION_CODE < ROOT_VERSION(6,18,0)
    //
    Ostap::Utils::Iterator it  ( lst  ) ; // only for ROOT < 6.18 
    RooAbsArg*   p = 0 ;
    unsigned int k = 0 ;
    while ( ( p = (RooAbsArg*) it.next() ) )
      {
	const RooAbsReal* r = static_cast<const RooAbsReal*>( p ) ;
	vct [ k ] = r->getVal ( nset ) ;  
	++k ;
      }
    //
#else
    //
    const unsigned int N = lst.size() ;
    for ( unsigned int k = 0 ; k < N ; ++k ) 
      {
	const RooAbsReal& r = static_cast<const RooAbsReal&>( lst [ k ] ) ;
	vct [ k ] = r.getVal ( nset ) ;  
      }
    //
#endif
    //
  }
  // ==========================================================================
  void assign 
(       RooAbsCollection& to   ,
  const RooAbsCollection& from ) 
{
  if ( &from == &to ) { return ; } // no self-assignment
  // ==========================================================================
#if ROOT_VERSION(6,26,0)  <= ROOT_VERSION_CODE 
  // ==========================================================================  
  to.assign ( from ) ;
  // ==========================================================================
#else 
  // ==========================================================================
  to.assignValueOnly ( from ) ;
  // ==========================================================================
#endif 
}
  // ==========================================================================
  inline long getValue ( const RooAbsCategory& c )
  {
    // ========================================================================
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE // ==============================
    // ========================================================================
    return  c.getCurrentIndex() ;
    // ========================================================================
#else
    // ========================================================================
    return  c.getIndex() ;    
    // ========================================================================    
#endif
    // ========================================================================
  }
  // ==========================================================================
  inline std::string getLabel ( const RooAbsCategory & c )
  {
    // ========================================================================
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE // ==============================
    // ========================================================================
    return  c.getCurrentLabel () ;
    // ========================================================================
#else
    // ========================================================================
    return  c.getLabel () ;    
    // ========================================================================    
#endif
    // ========================================================================
  }

  

  
  // ==========================================================================
} //                                             The end of anynymous namespace 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // LOCAL_ROOFIT_H
// ============================================================================
