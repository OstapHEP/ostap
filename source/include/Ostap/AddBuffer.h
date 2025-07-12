// ============================================================================
#ifndef OSTAP_ADDBUFFER_H 
#define OSTAP_ADDBUFFER_H 1
// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================
#include <string>
// ============================================================================
// ROOT 
// ============================================================================
#include "RtypesCore.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/IFuncs.h"
#include "Ostap/StatusCode.h"
#include "Ostap/ProgressConf.h"
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
  } // The en do fnamespace Ostap::Trees
  // =========================================================================
  /** @class AddBuffer 
   *  Add buffer to TTree
   */
  class AddBuffer
  {    
    // ========================================================================
  public:
    // ========================================================================
    /// constructor with the progress flag
    AddBuffer  ( const Ostap::Utils::ProgressConf& progress = false) ;
    // ========================================================================    
  public: 
    // ========================================================================
    /// add buffer to tree 
    Ostap::StatusCode
    add_buffer
    ( TTree*                                tree             ,
      const std::string&                    name             ,
      const Ostap::Utils::Buffer<Double_t>& buffer           ) const ;
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree* tree                                           ,
      const std::string&                   name             ,
      const Ostap::Utils::Buffer<Float_t>& buffer           ) const ;
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                               tree             ,
      const std::string&                   name             ,
      const Ostap::Utils::Buffer<Short_t>& buffer           ) const ;
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                                tree             ,
      const std::string&                    name             ,
      const Ostap::Utils::Buffer<UShort_t>& buffer           ) const ;
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                                tree             ,
      const std::string&                    name             ,
      const Ostap::Utils::Buffer<Int_t>&    buffer           ) const ; 
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                              tree             ,
      const std::string&                  name             ,
      const Ostap::Utils::Buffer<UInt_t>& buffer           ) const ;
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                                tree             ,
      const std::string&                    name             ,
      const Ostap::Utils::Buffer<Long64_t>& buffer           ) const ;
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                                 tree             ,
      const std::string&                     name             ,
      const Ostap::Utils::Buffer<ULong64_t>& buffer           ) const ;
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                              tree             ,
      const std::string&                  name             ,
      const Ostap::Utils::Buffer<Long_t>& buffer           ) const ;
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                               tree             ,
      const std::string&                   name             ,
      const Ostap::Utils::Buffer<ULong_t>& buffer           ) const ;
    // ========================================================================
    /** add <code>long double</code> buffer to TTree
     *  @attentoon actually doubel values are stored 
     */ 
    Ostap::StatusCode
    add_buffer
    ( TTree*                                   tree             ,
      const std::string&                       name             ,
      const Ostap::Utils::Buffer<long double>& buffer           ) const ; 
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                                   tree             ,
      const std::string&                       name             ,
      const Ostap::Utils::Buffer<signed char>& buffer           ) const ; 
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                                     tree             ,
      const std::string&                         name             ,
      const Ostap::Utils::Buffer<unsigned char>& buffer           ) const ;
    // ========================================================================
  public:
    // ========================================================================
    /// add several (same type) buffers at once 
    Ostap::StatusCode
    add_buffers
    ( TTree*                                     tree             ,
      const Ostap::Utils::Buffers<Double_t>&     buffers          ) const ;
    // ========================================================================
    /// add several (same type) buffers at once 
    Ostap::StatusCode
    add_buffers
    ( TTree*                                     tree             ,
      const Ostap::Utils::Buffers<Float_t>&      buffers          ) const ;
    // ========================================================================
  public: 
    // ========================================================================
    /// congfiguration of the progress bar 
    inline const Ostap::Utils::ProgressConf& progress () const
    { return m_progress ; }
    // ========================================================================
  private :
    // ========================================================================
    /// congfiguration of the progress bar 
    Ostap::Utils::ProgressConf m_progress { false } ; 
    // ========================================================================
  } ; //                                      The end of class Ostap::AddBuffer 
  // ==========================================================================
} //                                                 The end of namesapce Ostap 
// ============================================================================
#endif // OSTAP_ADDBUFFER_H
// ============================================================================
//                                                                      The END
// ============================================================================
