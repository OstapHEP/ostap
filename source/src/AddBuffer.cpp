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
#include "Ostap/ProgressBar.h"
#include "Ostap/AddBranch.h"
#include "Ostap/AddBuffer.h"
// ============================================================================
// Local stuff 
// ============================================================================
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
    const Ostap::Utils::Buffer<DATA>  buffer   ) 
  {
    //
    if ( !tree   ) { return INVALID_TREE   ; }
    //
    Ostap::Assert ( Ostap::Trees::valid_name_for_branch ( name )             ,
                    "Invalid name for branch:\"" + name + "\"" ,
                    "Ostap::AddBuffer::add_branch_"                ,
                    INVALID_BRANCH_NAME , __FILE__ , __LINE__  ) ;
    //    
    DATA bvalue  { buffer.value() } ;
    TBranch* branch = tree->Branch ( name.c_str() , &bvalue , ( name + vtype ).c_str() );
    Ostap::Assert ( branch ,
                    "Cannot create branch: " + name +
                    " for " + std::string ( typeid ( buffer ).name() ) , 
                    "Ostap::AddBuffer::add_branch"                 ,
                    CANNOT_CREATE_BRANCH , __FILE__ , __LINE__ ) ;
    //
    const Long64_t total    = tree -> GetEntries() ;
    const Long64_t nentries = buffer.size() < total ? buffer.size () : total ;
    //
    Ostap::Utils::ProgressBar bar ( total , progress ) ;
    for ( Long64_t entry = 0 ; entry < nentries ; ++entry , ++bar ) 
      {
        if ( tree->GetEntry ( entry ) < 0 ) { break ; };
        bvalue =  buffer [ entry ] ;
        branch -> Fill() ;
      }
    // =========================================================================
    // The rest (if any) is a constant .
    bvalue =  buffer.value () ;
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
    const Ostap::Utils::Buffer<long double> buffer   ) 
  {
    //
    if ( !tree   ) { return INVALID_TREE   ; }
    //
    Ostap::Assert ( Ostap::Trees::valid_name_for_branch ( name )             ,
                    "Invalid name for branch:\"" + name + "\"" ,
                    "Ostap::AddBuffer::add_branch_"                ,
                    INVALID_BRANCH_NAME , __FILE__ , __LINE__  ) ;
    //
    // ATTENTion: double here!!! 
    double bvalue  { clamp ( buffer.value() ) } ;
    TBranch* branch = tree->Branch ( name.c_str() , &bvalue , ( name + "/D" ).c_str() );
    Ostap::Assert ( branch ,
                    "Cannot create branch: " + name +
                    " for " + std::string ( typeid ( buffer ).name() ) , 
                    "Ostap::AddBuffer::add_branch"                 ,
                    CANNOT_CREATE_BRANCH , __FILE__ , __LINE__ ) ;
    //
    const Long64_t total    = tree -> GetEntries() ;
    const Long64_t nentries = buffer.size() < total ? buffer.size () : total ;
    //
    Ostap::Utils::ProgressBar bar ( total , progress ) ;
    for ( Long64_t entry = 0 ; entry < nentries ; ++entry , ++bar ) 
      {
        if ( tree->GetEntry ( entry ) < 0 ) { break ; };
        bvalue =  clamp ( buffer [ entry ] ) ;
        branch -> Fill() ;
      }
    // =========================================================================
    // The rest (if any) is a constant .
    bvalue = clamp ( buffer.value () ) ;
    for ( Long64_t entry = nentries ; entry < total  ; ++entry , ++bar ) 
      {
        if ( tree->GetEntry ( entry ) < 0 ) { break ; };
        branch->Fill() ;
      }
    //
    return Ostap::StatusCode::SUCCESS ;  
  }
  // ==========================================================================
  /// add several buffers at once 
  template <class DATA>
  inline Ostap::StatusCode 
  _add_buffers_ 
  ( TTree*                             tree     , 
    const Ostap::Utils::ProgressConf&  progress , 
    const std::string&                 vtype    ,
    const Ostap::Utils::Buffers<DATA>  buffers  ) 
  {
    //
    if ( !tree   ) { return INVALID_TREE ; }
    //
    const std::size_t N = buffers.size() ;
    if ( !N      ) { return Ostap::StatusCode::SUCCESS ; } // nothing to add
    //
    /// vector of values 
    std::vector<double>   dvalues  ( N , 0.0     ) ;
    /// vector of Branches 
    std::vector<TBranch*> branches ( N , nullptr ) ;
    ///
    std::size_t index = 0 ;
    for ( const auto& v : buffers )
      {
        const std::string name { v.first } ;
        Ostap::Assert ( Ostap::Trees::valid_name_for_branch ( name )             ,
                        "Invalid name for branch:\"" + name + "\"" ,
                        "Ostap::AddBuffer::add_buffer_"                ,
                        INVALID_BRANCH_NAME , __FILE__ , __LINE__  ) ;
        //
        TBranch* branch = tree->Branch ( name.c_str() , &dvalues [ index ] , ( name + vtype ).c_str() );
        Ostap::Assert ( branch ,
                        "Cannot create branch: " + name +
                        " for " + std::string ( typeid ( v.second ).name() ) , 
                        "Ostap::AddBuffer::add_branch"                       ,
                        CANNOT_CREATE_BRANCH , __FILE__ , __LINE__ ) ;
        //
        dvalues  [ index ] = v.second.value() ; 
        branches [ index ] = branch           ;
        ++index ;
        //
      }
    //
    const Long64_t total = tree -> GetEntries() ;
    Ostap::Utils::ProgressBar bar ( total , progress ) ;
    //
    for ( Long64_t entry = 0 ; entry < total ; ++entry , ++bar ) 
      {
        if ( tree->GetEntry ( entry ) < 0 ) { break ; };
        //
        std::size_t index = 0 ;
        for ( const auto& v : buffers )
          {
            dvalues  [ index ] = v.second [ entry ] ;
            branches [ index ] -> Fill ()           ;
            ++index ;
          }
      }
    //
    return Ostap::StatusCode::SUCCESS ;  
  }
} // ==========================================================================
// ============================================================================
// constructor
// ============================================================================
Ostap::AddBuffer::AddBuffer
( const Ostap::Utils::ProgressConf& progress )
  : m_progress ( progress )
{}
// ============================================================================
Ostap::StatusCode
Ostap::AddBuffer::add_buffer
( TTree*                                 tree     ,
  const std::string&                     name     ,
  const Ostap::Utils::Buffer<Double_t>&  buffer   ) const 
{ return ::_add_buffer_ ( tree , m_progress , name , "/D" , buffer ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::AddBuffer::add_buffer
( TTree*                                 tree     ,
  const std::string&                     name     ,
  const Ostap::Utils::Buffer<Float_t>&   buffer   ) const 
{ return ::_add_buffer_ ( tree , m_progress , name , "/F" , buffer ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::AddBuffer::add_buffer
( TTree*                                 tree     ,
  const std::string&                     name     ,
  const Ostap::Utils::Buffer<Short_t>&   buffer   ) const 
{ return ::_add_buffer_ ( tree , m_progress , name , "/S" , buffer ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::AddBuffer::add_buffer
( TTree*                                 tree     ,
  const std::string&                     name     ,
  const Ostap::Utils::Buffer<UShort_t>&  buffer   ) const 
{ return ::_add_buffer_ ( tree , m_progress , name , "/s" , buffer ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::AddBuffer::add_buffer
( TTree*                                 tree     ,
  const std::string&                     name     ,
  const Ostap::Utils::Buffer<Int_t>&     buffer   ) const 
{ return ::_add_buffer_ ( tree , m_progress , name , "/I" , buffer ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::AddBuffer::add_buffer
( TTree*                                 tree     ,
  const std::string&                     name     ,
  const Ostap::Utils::Buffer<UInt_t>&    buffer   ) const 
{ return ::_add_buffer_ ( tree , m_progress , name , "/i" , buffer ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::AddBuffer::add_buffer
( TTree*                                 tree     ,
  const std::string&                     name     ,
  const Ostap::Utils::Buffer<Long64_t>&  buffer   ) const 
{ return ::_add_buffer_ ( tree , m_progress , name , "/L" , buffer ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::AddBuffer::add_buffer
( TTree*                                 tree     ,
  const std::string&                     name     ,
  const Ostap::Utils::Buffer<ULong64_t>& buffer   ) const 
{ return ::_add_buffer_ ( tree , m_progress , name , "/l" , buffer ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::AddBuffer::add_buffer
( TTree*                                 tree     ,
  const std::string&                     name     ,
  const Ostap::Utils::Buffer<Long_t>&    buffer   ) const 
// { return ::_add_buffer_ ( tree , progress , name , "/G" , buffer ) ; }
{ return ::_add_buffer_ ( tree , m_progress , name , "/L" , buffer ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::AddBuffer::add_buffer
( TTree*                                 tree     ,
  const std::string&                     name     ,
  const Ostap::Utils::Buffer<ULong_t>&   buffer   ) const 
// { return ::_add_buffer_ ( tree , progress , name , "/g" , buffer ) ; }
{ return ::_add_buffer_ ( tree , m_progress , name , "/l" , buffer ) ; }
// ============================================================================
/** add <code>long dobule</code> buffer to TTree
 *  @attentoon actually doubel valeus are stored 
 */
// ============================================================================
Ostap::StatusCode
Ostap::AddBuffer::add_buffer
( TTree*                                   tree     ,
  const std::string&                       name     ,
  const Ostap::Utils::Buffer<long double>& buffer   ) const 
{ return ::_add_buffer_ ( tree , m_progress , name , buffer ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::AddBuffer::add_buffer
( TTree*                                     tree     ,
  const std::string&                         name     ,
  const Ostap::Utils::Buffer<signed char>&   buffer   ) const 
{ return ::_add_buffer_ ( tree , m_progress , name , "/B" , buffer ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::AddBuffer::add_buffer
( TTree*                                     tree     ,
  const std::string&                         name     ,
  const Ostap::Utils::Buffer<unsigned char>& buffer   ) const 
{ return ::_add_buffer_ ( tree , m_progress , name , "/b" , buffer ) ; }
// ============================================================================



// ============================================================================
Ostap::StatusCode
Ostap::AddBuffer::add_buffers
( TTree*                                  tree     ,
  const Ostap::Utils::Buffers<Double_t>&  buffers  ) const 
{ return ::_add_buffers_ ( tree , m_progress , "/D" , buffers ) ; }
// ============================================================================
Ostap::StatusCode
Ostap::AddBuffer::add_buffers
( TTree*                                  tree     ,
  const Ostap::Utils::Buffers<Float_t>&   buffers  ) const 
{ return ::_add_buffers_ ( tree , m_progress , "/F" , buffers ) ; }
// ============================================================================


// ============================================================================
//                                                                      The END 
// ============================================================================
