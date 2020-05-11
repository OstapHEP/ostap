// ============================================================================
// Include files 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/MoreRooFit2.h"
// ============================================================================
/** @file
 *  implementaton of various small additions to RooFit 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
 *  @date 2020-05-11
 */
// ============================================================================
namespace
{
  // ==========================================================================
  inline std::string name_  ( const std::string& name  ,
                              const std::string& op    ,
                              const RooAbsReal&  v     )
  { return name.empty  () ? ( op + "_" + v.GetName()       ) : name   ; }
  // ==========================================================================
  inline std::string title_ ( const std::string& title ,
                              const std::string& op    ,
                              const RooAbsReal&  v     )
  { return title.empty () ? ( op + "(" + v.GetName() + ")" ) : title  ; }
  // ==========================================================================
}
// ============================================================================
ClassImp(Ostap::MoreRooFit::M2Q ) ;
ClassImp(Ostap::MoreRooFit::Q2M ) ;
// ============================================================================
// constructor
// ============================================================================
Ostap::MoreRooFit::M2Q::M2Q
( const std::string& name  ,
  const std::string& title , 
  RooAbsReal&        m     ,
  const double       m1    ,
  const double       m2    )
  : RooAbsReal ( name_  ( name  , "m2q" , m ).c_str() ,
                 title_ ( title , "m2q" , m ).c_str() )
  , m_ps ( std::abs ( m1 ) , std::abs ( m2 ) )
  , m_m  ( "!m" , "mass" , this , m  ) 
{}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::MoreRooFit::M2Q::M2Q
( const Ostap::MoreRooFit::M2Q& right ,
  const char*                   name  )
  : RooAbsReal ( right , name )
  , m_ps ( right.m_ps )
  , m_m  ( "!m" , this , right.m_m )
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::MoreRooFit::M2Q::~M2Q(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::MoreRooFit::M2Q*
Ostap::MoreRooFit::M2Q::clone ( const char* newname ) const
{ return new M2Q ( *this , newname ) ; }
// ============================================================================
// the actual evaluation of the result 
// ============================================================================
Double_t Ostap::MoreRooFit::M2Q::evaluate () const 
{
  const double m = m_m ;
  return m_ps.q_ ( m ) ;
}
// ============================================================================


// ============================================================================
// constructor
// ============================================================================
Ostap::MoreRooFit::Q2M::Q2M 
( const std::string& name  , 
  const std::string& title ,
  RooAbsReal&        q     ,
  const double       m1    ,
  const double       m2    )
  : M2Q ( name_  ( name  , "q2m" , q ) ,
          title_ ( title , "q2m" , q ) , q , m1 , m2 )
{}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::MoreRooFit::Q2M::Q2M
( const Ostap::MoreRooFit::Q2M& right ,
  const char*                   name  )
  : M2Q ( right , name )
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::MoreRooFit::Q2M::~Q2M(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::MoreRooFit::Q2M*
Ostap::MoreRooFit::Q2M::clone ( const char* newname ) const
{ return new Q2M( *this , newname ) ; }
// ============================================================================
// the actual evaluation of the result 
// ============================================================================
Double_t Ostap::MoreRooFit::Q2M::evaluate () const 
{
  const double q = m_m ;
  return q <= 0 ? 0.0 : m_ps.m_ ( q ) ;
}
// ============================================================================


// ============================================================================
//                                                                      The END 
// ============================================================================


