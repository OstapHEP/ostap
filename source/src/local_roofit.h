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
#include "status_codes.h"
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
  inline std::size_t size ( const RooAbsCollection& lst )  { return lst.size () ; }
  // ==========================================================================
  /** copy from list to list
   *  @param from objects to be copied from this list 
   *  @param to   objects to be copied to this proxy 
   */
  inline void copy 
  ( const RooAbsCollection&  from ,
    RooAbsCollection&        to   ) 
  { for ( auto* c : from ) { to.add ( *c ) ; } }
  // ==========================================================================
  /** copy RooAbsReal from list to list
   *  @param from objects to be copied from this list 
   *  @param to   objects to be copied to this proxy 
   */
  inline std::size_t copy_real 
  ( const RooAbsCollection&  from                                    ,
    RooAbsCollection&        to                                      , 
    const std::string&       message = "Variable is not RooAbsReal!" ,
    const std::string&       tag     = "Ostap::copy_real"            ,
    const char*              file    = nullptr                       ,
    const std::size_t        line    = 0                             )
  {
    //
    for ( auto* c : from ) 
      {
        Ostap::Assert ( nullptr != dynamic_cast<RooAbsReal*> ( c ) ,
                        message , tag , INVALID_ABSREAL , file , line ) ;
        to.add ( *c ) ;   
      }
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
    const std::string&       tag     = "Ostap::copy_real"            ,
    const char*              file    = nullptr                       ,
    const std::size_t        line    = 0                             )
  {
    //
    for ( auto* c : from ) 
      {
        Ostap::Assert ( nullptr != dynamic_cast<RooAbsReal*> ( c )    , 
                        message , tag , INVALID_ABSREAL , file , line ) ;
        to.add ( *c ) ;   
      }
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
    const unsigned int N = lst.size() ;
    for ( unsigned int k = 0 ; k < N ; ++k ) 
      {
        const RooAbsReal& r = static_cast<const RooAbsReal&>( lst [ k ] ) ;
        obj.setPar ( k + shift , r.getVal ( nset ) ) ;  
      }
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
    const unsigned int N = lst.size() ;
    for ( unsigned int k = 0 ; k < N ; ++k ) 
      {
	const RooAbsReal& r = static_cast<const RooAbsReal&>( lst [ k ] ) ;
	vct [ k ] = r.getVal ( nset ) ;  
      }
    //
  }
  // ==========================================================================
  void assign 
  ( RooAbsCollection&       to   ,
    const RooAbsCollection& from ) 
  {
    if ( &from == &to ) { return ; } // no self-assignment
    // ==========================================================================
#if ROOT_VERSION(6,26,0)  <= ROOT_VERSION_CODE 
    // ==========================================================================  
    to.assign          ( from ) ;
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
    return  c.getCurrentIndex() ;
    // ========================================================================
  }
  // ==========================================================================
  inline std::string getLabel ( const RooAbsCategory & c )
  {
    // ========================================================================
    return  c.getCurrentLabel () ;
    // ========================================================================
  }
  // ==========================================================================
} //                                             The end of anynymous namespace 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // LOCAL_ROOFIT_H
// ============================================================================
