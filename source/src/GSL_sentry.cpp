// $Id$
// ============================================================================
// Include files 
// ============================================================================
// STD & STL 
// ============================================================================
#include <iostream>
#include <sstream>
#include <map>
#include <tuple>
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/Error2Exception.h"
#include "Ostap/StatusCode.h"
#include "Ostap/ToStream.h"
// ============================================================================
// local
// ============================================================================
#include "GSL_sentry.h"
#include "syncedcache.h"
// ============================================================================
/** @file 
 *  Implementation file for class Ostap::Math::GSL::GSL_Error_Handler
 *  
 *  @see Ostap::Math::GSL::GSL_Error_Handler
 *  @date 2012-05-27 
 *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
 */
// ============================================================================
namespace
{
  // ==========================================================================
  class Cache 
  {
  public: 
    // ========================================================================
    static void add 
    ( const char * reason  ,
      const char * file    ,
      int          line    ,
      int          errcode ) 
    {
      CACHE::Lock lock  { s_CACHE.mutex() } ;
      const ITEM key = std::make_tuple ( reason , file , line , errcode ) ;
      s_CACHE.container() [ key ] += 1 ;
    }
    // =========================================================================    
    static std::size_t size() 
    {
      CACHE::Lock lock  { s_CACHE.mutex() } ;
      return s_CACHE->size () ;
    }
    // =========================================================================    
    /// get summary of GSL errors 
    static Ostap::Utils::GslCount::Table table () 
    {
      CACHE::Lock lock  { s_CACHE.mutex() } ;
      //
      Ostap::Utils::GslCount::Table _table ;
      _table.reserve ( s_CACHE->size () ) ;
      //
      for ( auto it = s_CACHE->begin() ; s_CACHE->end() != it ; ++it ) 
      {
        const std::string reason  = std::get<0>  ( it -> first  ) ;
        const std::string file    = std::get<1>  ( it -> first  ) ;
        const int         line    = std::get<2>  ( it -> first  ) ;
        const int         errcode = std::get<3>  ( it -> first  ) ;
        const std::string message = gsl_strerror ( errcode )      ;
        const int         number  = it -> second                  ;
        //
        Ostap::Utils::GslCount::Row row ; row.reserve ( 6 ) ;
        //
        row.push_back ( Ostap::Utils::toString ( number  ) ) ;
        row.push_back ( Ostap::Utils::toString ( errcode ) ) ;
        row.push_back ( message ) ;        
        row.push_back ( reason  ) ;
        row.push_back ( file    ) ;
        row.push_back ( Ostap::Utils::toString ( line    ) ) ;
        //
        _table.push_back ( row ) ;
      }
      //
      return _table ;
    }
    // =========================================================================    
    static std::size_t clear() 
    {
      CACHE::Lock lock  { s_CACHE.mutex() } ;
      const std::size_t size = s_CACHE->size() ;
      s_CACHE->clear() ;
      return size ;
    }
    // =========================================================================
    ~Cache () 
    {
      CACHE::Lock lock { s_CACHE.mutex() } ;
      if ( 0 < s_CACHE->size() ) 
      {
        std::cerr << "Summary of GSL Errors " << std::endl ;
        for ( auto it = s_CACHE->begin() ; s_CACHE->end() != it ; ++it ) 
        {
          const std::string reason  = std::get<0> ( it->first ) ;
          const std::string file    = std::get<1> ( it->first ) ;
          const int         line    = std::get<2> ( it->first ) ;
          const int         errcode = std::get<3> ( it->first ) ;
          //
          std::cerr 
            << " GSL_ERROR : "   
            << "#"  << it->second << "  "
            << "\t" << errcode << "/'" << gsl_strerror ( errcode ) << "'"
            << "\t reason '"     
            << reason    << "' "
            << "\t file/line '"  
            << file      << "'/" << line 
            << std::endl ;
        }
      }      
      s_CACHE->clear () ;
    }
    // =========================================================================
  private :
    // =========================================================================
    // item in the map 
    typedef std::tuple<std::string,std::string,int,int> ITEM  ;
    // map type itself
    typedef std::map<ITEM,unsigned long>                MAP   ;  
    // synched map 
    typedef SyncedCache<MAP>                            CACHE ;
    // =========================================================================
    // cache itself 
    static CACHE s_CACHE ;
    // =========================================================================
  } ;
  // ===========================================================================
  Cache::CACHE Cache::s_CACHE = CACHE{} ;
  // ===========================================================================
  /** @var s_cache 
   *  The actual cache object 
   */
  Cache s_cache ;
  // ==========================================================================


  // ==========================================================================
  /// print errors  
  void GSL_print_error
  ( const char * reason    ,
    const char * file      ,
    int          line      ,
    int          gsl_errno ) 
  {
    std::cerr 
      << " GSL_ERROR : "   
      << gsl_errno << "/'" << gsl_strerror ( gsl_errno ) << "'"
      << "\t reason '"     
      << reason    << "' "
      << "\t file/line '"  
      << file      << "'/" << line 
      << std::endl ;  
  }
  // ==========================================================================
  /// ignore errors 
  void GSL_ignore_error
  ( const char * /* reason    */ ,
    const char * /* file      */ ,
    int          /* line      */ ,
    int          /* gsl_errno */ ) {}
  // ==========================================================================
  /// convert errors to exceptions 
  void GSL_exception_error
  ( const char * reason    ,
    const char * file      ,
    int          line      ,
    int          gsl_errno ) 
  {
    std::string tag = "GSL/Error" ;
    std::ostringstream ss ;
    ss << gsl_strerror ( gsl_errno )
       << "(" << gsl_errno << ") "
       << reason ; 
    Ostap::Assert ( false                            ,
                    tag + ": " + ss.str()            , 
                    tag                              , 
                    100000 + gsl_errno , file , line ) ; 
  }
  // ==========================================================================
  /// silently count errors 
  void GSL_count_error 
  ( const char * reason    ,
    const char * file      ,
    int          line      ,
    int          gsl_errno ) 
  {
    s_cache.add ( reason , file , line , gsl_errno ) ;
  }
  // ==========================================================================
}
// ============================================================================
// constructor: make use of Gsl Error Handler: print error to stderr 
// ============================================================================
Ostap::Utils::GslError::GslError ( Ostap::Utils::GslError::handler* h     , 
                                   const bool                       force )
  : m_previous ( gsl_set_error_handler ( h ) ) 
  , m_force    ( force ) 
{ 
  //
  if ( m_previous && !m_force ) { gsl_set_error_handler ( m_previous ) ; }
  static_assert( std::is_same<handler,gsl_error_handler_t>::value  ,
                 "``handler'' type is not ``gsl_error_handler_t''" ) ;
}
// ============================================================================
// destructor: stop using the error  handler 
// ============================================================================
Ostap::Utils::GslError::~GslError() { gsl_set_error_handler ( m_previous ) ; }
// ============================================================================
// constructor: make use of Gsl Error Handler: print error to stderr 
// ============================================================================
Ostap::Utils::GslError::GslError 
( const bool force ) 
  : GslError ( &GSL_print_error , force  ) {}
// ============================================================================
// constructor: make use of Gsl Error Handler: print error to stderr 
// ============================================================================
Ostap::Utils::GslIgnore::GslIgnore 
( const bool force ) 
  : GslError( &GSL_ignore_error , force  ) {}
// ============================================================================
// constructor: make use of Gsl Error Handler: print error to stderr 
// ============================================================================
Ostap::Utils::GslCount::GslCount 
( const bool force ) 
  : GslError( &GSL_count_error , force  ) {}
// ============================================================================
// get total number of errors 
// ============================================================================
std::size_t  Ostap::Utils::GslCount::size ()
{ return Cache::size () ; }
// ============================================================================
// clear summary of errors  
// ============================================================================
std::size_t  Ostap::Utils::GslCount::clear ()
{ return Cache::clear () ; }
// ============================================================================
// get all errors in a form of the  table 
// ============================================================================
Ostap::Utils::GslCount::Table 
Ostap::Utils::GslCount::table() 
{ return Cache::table () ; }
// ============================================================================
  
// ============================================================================
// constructor: make use of Gsl Error Handler: print error to stderr 
// ============================================================================
Ostap::Utils::GslException::GslException
( const bool force ) 
  : GslError( &GSL_exception_error , force ) {}
// ============================================================================
//                                                                      The END 
// ============================================================================
