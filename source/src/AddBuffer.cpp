// ============================================================================
// Include files 
// ============================================================================
#include <string>
#include <typeinfo> 
// ============================================================================
// STD/STL
// ============================================================================
#include <limits>
// ============================================================================
// ROOT
// ============================================================================
#include "TTree.h"
#include "TBranch.h"
// ============================================================================
// ROOT/RooFit 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatusCode.h"
#include "Ostap/ToStream.h"
#include "Ostap/ProgressBar.h"
#include "Ostap/AddBranch.h"
#include "Ostap/AddBuffer.h"
// ============================================================================
// Local stuff 
// ============================================================================
#include "Exception.h" 
#include "status_codes.h" 
// ============================================================================
/** @file
 *  Implementation file for function Ostap::Trees::add_branch 
 *  @see Ostap::Trees::add_branch 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2019-05-14
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  static_assert ( std::numeric_limits<long double> ::is_specialized ,
                  "std::numeric_limits<double>      is not specialized!" ) ;
  static_assert ( std::numeric_limits<long double> ::is_specialized ,
                  "std::numeric_limits<long double> is not specialized!" ) ;
  // =========================================================================
  const long double s_dmax =  std::numeric_limits<double>::max () ; 
  const long double s_dmin = -std::numeric_limits<double>::max () ;
  // =========================================================================
  template <class DATA>
  inline Ostap::StatusCode 
  _add_buffer_ 
  ( TTree*                            tree     , 
    const Ostap::Utils::ProgressConf& progress , 
    const std::string&                name     ,                
    const std::string&                vtype    ,
    const Ostap::Trees::Buffer<DATA>  buffer   ) 
  {
    //
    if ( !tree   ) { return INVALID_TREE   ; }
    //
    Ostap::Assert ( Ostap::Trees::valid_name_for_branch ( name )             ,
                    "Invalid name for branch:\"" + name + "\"" ,
                    "Ostap::Trees::add_branch_"                ,
                    INVALID_BRANCH_NAME , __FILE__ , __LINE__  ) ;
    //    
    DATA bvalue  { buffer.value() } ;
    TBranch* branch = tree->Branch ( name.c_str() , &bvalue , ( name + vtype ).c_str() );
    Ostap::Assert ( branch ,
                    "Cannot create branch: " + name +
                    " for " + std::string ( typeid ( buffer ).name() ) , 
                    "Ostap::Trees::add_branch"                 ,
                    CANNOT_CREATE_BRANCH , __FILE__ , __LINE__ ) ;
    //
    const Long64_t total    = tree -> GetEntries() ;
    const Long64_t nentries = buffer.size() < total ? buffer.size () : total ;
    //
    for ( Long64_t entry = 0 ; entry < nentries ; ++entry ) 
      {
        if ( tree->GetEntry ( entry ) < 0 ) { break ; };
        bvalue =  buffer [ entry ] ;
        branch -> Fill() ;
      }
    // =========================================================================
    // The rest (if any) is a constant .
    bvalue =  buffer.value () ;
    Ostap::Utils::ProgressBar bar ( nentries , progress ) ;
    for ( Long64_t entry = nentries ; entry < total  ; ++entry , ++bar ) 
      {
        if ( tree->GetEntry ( entry ) < 0 ) { break ; };
        branch->Fill() ;
      }
    //
    return Ostap::StatusCode::SUCCESS ;  
  }
  // ==========================================================================
  inline
  double clamp ( const long double value )
  {
    if ( value <= s_dmin ) { return static_cast<double> ( s_dmin ) ;}
    if ( value >= s_dmax ) { return static_cast<double> ( s_dmax ) ;}
    return  static_cast<double> ( value ) ;
  }
  // ==========================================================================
  inline Ostap::StatusCode 
  _add_buffer_ 
  ( TTree*                                  tree     , 
    const Ostap::Utils::ProgressConf&       progress , 
    const std::string&                      name     ,                
    const Ostap::Trees::Buffer<long double> buffer   ) 
  {
    //
    if ( !tree   ) { return INVALID_TREE   ; }
    //
    Ostap::Assert ( Ostap::Trees::valid_name_for_branch ( name )             ,
                    "Invalid name for branch:\"" + name + "\"" ,
                    "Ostap::Trees::add_branch_"                ,
                    INVALID_BRANCH_NAME , __FILE__ , __LINE__  ) ;
    //
    // ATTENTONO: double here!!! 
    double bvalue  { clamp ( buffer.value() ) } ;
    TBranch* branch = tree->Branch ( name.c_str() , &bvalue , ( name + "/D" ).c_str() );
    Ostap::Assert ( branch ,
                    "Cannot create branch: " + name +
                    " for " + std::string ( typeid ( buffer ).name() ) , 
                    "Ostap::Trees::add_branch"                 ,
                    CANNOT_CREATE_BRANCH , __FILE__ , __LINE__ ) ;
    //
    const Long64_t total    = tree -> GetEntries() ;
    const Long64_t nentries = buffer.size() < total ? buffer.size () : total ;
    //
    for ( Long64_t entry = 0 ; entry < nentries ; ++entry ) 
      {
        if ( tree->GetEntry ( entry ) < 0 ) { break ; };
        bvalue =  clamp ( buffer [ entry ] ) ;
        branch -> Fill() ;
      }
    // =========================================================================
    // The rest (if any) is a constant .
    bvalue = clamp ( buffer.value () ) ;
    Ostap::Utils::ProgressBar bar ( nentries , progress ) ;
    for ( Long64_t entry = nentries ; entry < total  ; ++entry , ++bar ) 
      {
        if ( tree->GetEntry ( entry ) < 0 ) { break ; };
        branch->Fill() ;
      }
    //
    return Ostap::StatusCode::SUCCESS ;  
  }
  // ==========================================================================
} // ==========================================================================
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_buffer
( TTree*                                 tree     ,
  const std::string&                     name     ,
  const Ostap::Trees::Buffer<Double_t>&  buffer   ,
  const Ostap::Utils::ProgressConf&      progress )                  
{ return ::_add_buffer_ ( tree , progress , name , "/D" , buffer ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_buffer
( TTree*                                 tree     ,
  const std::string&                     name     ,
  const Ostap::Trees::Buffer<Float_t>&   buffer   ,
  const Ostap::Utils::ProgressConf&      progress )                  
{ return ::_add_buffer_ ( tree , progress , name , "/F" , buffer ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_buffer
( TTree*                                 tree     ,
  const std::string&                     name     ,
  const Ostap::Trees::Buffer<Char_t>&    buffer   ,
  const Ostap::Utils::ProgressConf&      progress )                  
{ return ::_add_buffer_ ( tree , progress , name , "/B" , buffer ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_buffer
( TTree*                                 tree     ,
  const std::string&                     name     ,
  const Ostap::Trees::Buffer<UChar_t>&   buffer   ,
  const Ostap::Utils::ProgressConf&      progress )                  
{ return ::_add_buffer_ ( tree , progress , name , "/b" , buffer ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_buffer
( TTree*                                 tree     ,
  const std::string&                     name     ,
  const Ostap::Trees::Buffer<Short_t>&   buffer   ,
  const Ostap::Utils::ProgressConf&      progress )                  
{ return ::_add_buffer_ ( tree , progress , name , "/S" , buffer ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_buffer
( TTree*                                 tree     ,
  const std::string&                     name     ,
  const Ostap::Trees::Buffer<UShort_t>&  buffer   ,
  const Ostap::Utils::ProgressConf&      progress )                  
{ return ::_add_buffer_ ( tree , progress , name , "/s" , buffer ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_buffer
( TTree*                                 tree     ,
  const std::string&                     name     ,
  const Ostap::Trees::Buffer<Int_t>&     buffer   ,
  const Ostap::Utils::ProgressConf&      progress )                  
{ return ::_add_buffer_ ( tree , progress , name , "/I" , buffer ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_buffer
( TTree*                                 tree     ,
  const std::string&                     name     ,
  const Ostap::Trees::Buffer<UInt_t>&    buffer   ,
  const Ostap::Utils::ProgressConf&      progress )                  
{ return ::_add_buffer_ ( tree , progress , name , "/i" , buffer ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_buffer
( TTree*                                 tree     ,
  const std::string&                     name     ,
  const Ostap::Trees::Buffer<Long64_t>&  buffer   ,
  const Ostap::Utils::ProgressConf&      progress )                  
{ return ::_add_buffer_ ( tree , progress , name , "/L" , buffer ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_buffer
( TTree*                                 tree     ,
  const std::string&                     name     ,
  const Ostap::Trees::Buffer<ULong64_t>& buffer   ,
  const Ostap::Utils::ProgressConf&      progress )                  
{ return ::_add_buffer_ ( tree , progress , name , "/l" , buffer ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_buffer
( TTree*                                 tree     ,
  const std::string&                     name     ,
  const Ostap::Trees::Buffer<Long_t>&    buffer   ,
  const Ostap::Utils::ProgressConf&      progress )                  
// { return ::_add_buffer_ ( tree , progress , name , "/G" , buffer ) ; }
{ return ::_add_buffer_ ( tree , progress , name , "/L" , buffer ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_buffer
( TTree*                                 tree     ,
  const std::string&                     name     ,
  const Ostap::Trees::Buffer<ULong_t>&   buffer   ,
  const Ostap::Utils::ProgressConf&      progress )                  
// { return ::_add_buffer_ ( tree , progress , name , "/g" , buffer ) ; }
{ return ::_add_buffer_ ( tree , progress , name , "/l" , buffer ) ; }
// ============================================================================
/** add <code>long dobule</code> buffer to TTree
 *  @attentoon actually doubel valeus are stored 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_buffer
( TTree*                                   tree     ,
  const std::string&                       name     ,
  const Ostap::Trees::Buffer<long double>& buffer   ,
  const Ostap::Utils::ProgressConf&      progress )                  
{ return ::_add_buffer_ ( tree , progress , name , buffer ) ; }
// ============================================================================


// ============================================================================
//                                                                      The END 
// ============================================================================
