// ============================================================================
#ifndef OSTAP_ADDBUFFER_H 
#define OSTAP_ADDBUFFER_H 1
// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================ 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/IFuncs.h"
#include "Ostap/StatusCode.h"
#include "Ostap/ProgressBar.h"
#include "Ostap/Buffer.h"
// ============================================================================
// Forward declarations from ROOT 
// ============================================================================
class TTree ; 
// ============================================================================
namespace Ostap
{  
  // ==========================================================================
  namespace Trees
  { 
    // ========================================================================
    /** valid name for branch ?
     *  - not empty 
     *  - no blanks 
     *  - no special symbols 
     */
    bool valid_name_for_branch ( const std::string& name ) ;
    // ========================================================================
    /// add buffer to tree 
    Ostap::StatusCode
    add_buffer
    ( TTree*                                tree             ,
      const std::string&                    name             ,
      const Ostap::Utils::Buffer<Double_t>& buffer           , 
      const Ostap::Utils::ProgressConf&     progress = false ) ;                 
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree* tree                                           ,
      const std::string&                   name             ,
      const Ostap::Utils::Buffer<Float_t>& buffer           ,
      const Ostap::Utils::ProgressConf&    progress = false ) ; 
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                               tree             ,
      const std::string&                   name             ,
      const Ostap::Utils::Buffer<Short_t>& buffer           , 
      const Ostap::Utils::ProgressConf&    progress = false ) ;
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                                tree             ,
      const std::string&                    name             ,
      const Ostap::Utils::Buffer<UShort_t>& buffer           , 
      const Ostap::Utils::ProgressConf&     progress = false ) ; 
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                                tree             ,
      const std::string&                    name             ,
      const Ostap::Utils::Buffer<Int_t>&    buffer           , 
      const Ostap::Utils::ProgressConf&     progress = false ) ; 
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                              tree             ,
      const std::string&                  name             ,
      const Ostap::Utils::Buffer<UInt_t>& buffer           , 
      const Ostap::Utils::ProgressConf&   progress = false ) ; 
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                                tree             ,
      const std::string&                    name             ,
      const Ostap::Utils::Buffer<Long64_t>& buffer           , 
      const Ostap::Utils::ProgressConf&     progress = false ) ; 
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                                 tree             ,
      const std::string&                     name             ,
      const Ostap::Utils::Buffer<ULong64_t>& buffer           , 
      const Ostap::Utils::ProgressConf&      progress = false ) ;
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                              tree             ,
      const std::string&                  name             ,
      const Ostap::Utils::Buffer<Long_t>& buffer           , 
      const Ostap::Utils::ProgressConf&   progress = false ) ; 
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                               tree             ,
      const std::string&                   name             ,
      const Ostap::Utils::Buffer<ULong_t>& buffer           , 
      const Ostap::Utils::ProgressConf&    progress = false ) ;
    // ========================================================================
    /** add <code>long dobule</code> buffer to TTree
     *  @attentoon actually doubel valeus are stored 
     */ 
    Ostap::StatusCode
    add_buffer
    ( TTree*                                   tree             ,
      const std::string&                       name             ,
      const Ostap::Utils::Buffer<long double>& buffer           , 
      const Ostap::Utils::ProgressConf&        progress = false ) ;                 
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                                   tree             ,
      const std::string&                       name             ,
      const Ostap::Utils::Buffer<signed char>& buffer           , 
      const Ostap::Utils::ProgressConf&        progress = false ) ; 
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                                     tree             ,
      const std::string&                         name             ,
      const Ostap::Utils::Buffer<unsigned char>& buffer           , 
      const Ostap::Utils::ProgressConf&          progress = false ) ; 
    // ========================================================================
  } //                                        The end of namespace Ostap::Trees 
  // ==========================================================================
} //                                                 The end of namesapce Ostap 
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_ADDBUFFER_H
// ============================================================================
