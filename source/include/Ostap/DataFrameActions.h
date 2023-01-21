// ============================================================================
#ifndef OSTAP_DATAFRAMEACTIONS_H 
#define OSTAP_DATAFRAMEACTIONS_H 1
// ============================================================================
// Inclodue files
// ============================================================================
// ROOT 
// ============================================================================
#include "RVersion.h"   // ROOT 
#include "RtypesCore.h" // ROOT 
// ============================================================================
#if ROOT_VERSION(6,14,0) <= ROOT_VERSION_CODE 
#include "ROOT/RDataFrame.hxx"  // ROOT 
#else 
#include "ROOT/TDataFrame.hxx"  // ROOT
#endif
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/DataFrame.h"
#include "Ostap/StatEntity.h"
#include "Ostap/WStatEntity.h"
// ============================================================================
/// ONLY starting from ROOT 6.16
// ============================================================================
#if ROOT_VERSION(6,16,0) <= ROOT_VERSION_CODE
// ============================================================================
namespace ROOT 
{
  // ==========================================================================
  namespace Detail 
  {
    // ========================================================================
    namespace RDF 
    {
      // ======================================================================
      /** @class StatVar
       *  Helper class to get statitsics for the column in DataFrame 
       *  @see Ostap::StatEntity 
       *  @see Ostap::DataFrame 
       */
      class StatVar : public RActionImpl<StatVar> 
      {
      public:
        // ====================================================================
        /// define the resutl type 
        using Result_t = Ostap::StatEntity ;
        // ====================================================================
      public:
        // ====================================================================
        /// default constructor 
        StatVar () ;
        /// Move constructor 
        StatVar (       StatVar&& ) = default ;
        /// Copy constructor is disabled 
        StatVar ( const StatVar&  ) = delete ;
        // ====================================================================
      public:
        // ====================================================================
        /// initialize (empty) 
        void InitTask   ( TTreeReader * , unsigned int ) {} ;
        /// initialize (empty) 
        void Initialize () {} ;
        /// finalize : sum over the slots 
        void Finalize   () ;
        /// who am I ?
        std::string GetActionName() { return "StatVar" ; }
        // ====================================================================
      public:
        // ====================================================================
        /// The basic method: increment the counter 
        void Exec ( unsigned int slot , double value ) 
        { m_slots [ slot % m_N ] += value ; } 
        // ====================================================================
        /// The basic method: increment the counter for the vector-like columns       
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif 
        void Exec ( unsigned int slot , const T &vs )
        { Ostap::StatEntity& e = m_slots [ slot % m_N ] ; for ( const auto & v : vs ) { e += v ; } }
        // ===============================================================================
      public:
        // ===============================================================================
        /// Get the result 
        std::shared_ptr<Result_t> GetResultPtr () const { return m_result ; }
        /// get partial result for the given slot 
        Result_t& PartialUpdate ( unsigned int slot ) { return m_slots [ slot % m_N ] ; }
        // ===============================================================================
      private:
        // ====================================================================
        /// the final result 
        const std::shared_ptr<Ostap::StatEntity> m_result {}    ;
        /// size of m_slots 
        unsigned long                            m_N      { 1 } ;
        /// (current) results per  slot 
        std::vector<Ostap::StatEntity>           m_slots  { 1 } ;
        // ===================================================================
      } ; //                        The end of class ROOT::Detail::RDF::StatVar
      // ======================================================================
      /** @class WStatVar
       *  Helper class to get weighted statitsics for the column in DataFrame 
       *  @see Ostap::WStatEntity 
       *  @see Ostap::DataFrame 
       */
      class WStatVar : public RActionImpl<WStatVar> 
      {
      public:
        // ====================================================================
        /// define the resutl type 
        using Result_t = Ostap::WStatEntity ;
        // ====================================================================
      public:
        // ====================================================================
        /// default constructor 
        WStatVar () ;
        /// Move constructor 
        WStatVar (       WStatVar&& ) = default ;
        /// Copy constructor is disabled 
        WStatVar ( const WStatVar&  ) = delete ;
        // ====================================================================
      public:
        // ====================================================================
        /// initialize (empty) 
        void InitTask   ( TTreeReader * , unsigned int ) {} ;
        /// initialize (empty) 
        void Initialize () {} ;
        /// finalize : sum over the slots 
        void Finalize   () ;
        /// who am I ?
        std::string GetActionName() { return "WStatVar" ; }
        // ====================================================================
      public:
        // ====================================================================
        /// The basic method: increment the counter 
        void Exec ( unsigned int slot , double value , double weight = 1 ) 
        { m_slots [ slot % m_N ].add ( value , weight ) ; } 
        // ====================================================================
        /// The basic method: increment the counter for the vector-like column of values        
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif 
        void Exec ( unsigned int slot , const T &vs , const double weight = 1 )
        {
          Ostap::WStatEntity& e = m_slots [ slot % m_N ] ; for ( const auto & v : vs ) { e.add ( v , weight ) ; } 
        }
        // ====================================================================
        /// The basic method: increment the counter for the vector-like column of weight 
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif
        void Exec ( unsigned int slot , const double value , const T &ws )
        { Ostap::WStatEntity& e = m_slots [ slot % m_N ] ; for ( const auto & w : ws ) { e.add ( value , w   ) ; } }
        // ====================================================================
      public:
        // ====================================================================
        /// Get the result 
        std::shared_ptr<Result_t> GetResultPtr () const { return m_result ; }
        /// get partial result for the given slot 
        Result_t& PartialUpdate ( unsigned int slot ) { return m_slots [ slot % m_N ] ; }
        // ====================================================================
      private:
        // ====================================================================
        /// the final result 
        const std::shared_ptr<Ostap::WStatEntity> m_result {}    ;
        /// size of m_slots 
        unsigned long                             m_N      { 1 } ;
        /// (current) results per  slot 
        std::vector<Ostap::WStatEntity>           m_slots  { 1 } ;
        // ===================================================================
      } ; //                        The end of class ROOT::Detail::RDF::StatVar 
      // ======================================================================
    } //                                 The end of namespace ROOT::Detail::RDF
    // ========================================================================
  } //                                        The end of namespace ROOT::Detail
  // ==========================================================================
} //                                                  The end of namespace ROOT 
// ============================================================================
namespace Ostap 
{
  // ==========================================================================
  namespace Actions 
  {
    // ========================================================================
    using StatVar  = ROOT::Detail::RDF::StatVar  ;
    using WStatVar = ROOT::Detail::RDF::WStatVar ;
    // ========================================================================
  }
  // ==========================================================================
}
// ============================================================================
#endif // #if ROOT_VERSION_CODE >= ROOT_VERSION(6,16,0)
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_DATAFRAMEACTIONS_H
// ============================================================================
