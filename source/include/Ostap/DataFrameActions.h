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
#include "Ostap/Polynomials.h"
#include "Ostap/Bernstein.h"
#include "Ostap/Bernstein2D.h"
#include "Ostap/Bernstein3D.h"
#include "Ostap/Parameterization.h"
#include "Ostap/Moments.h"
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
        /// define the result type 
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
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto & v : vs ) { e += v ; } }
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
        const std::shared_ptr<Result_t> m_result {   } ;
        /// size of m_slots 
        unsigned long                   m_N      { 1 } ;
        /// (current) results per  slot 
        std::vector<Result_t>           m_slots  { 1 } ;
        // ====================================================================
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
          Result_t& e = m_slots [ slot % m_N ] ; for ( const auto & v : vs ) { e.add ( v , weight ) ; } 
        }
        // ====================================================================
        /// The basic method: increment the counter for the vector-like column of weight 
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif
        void Exec ( unsigned int slot , const double value , const T &ws )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto & w : ws ) { e.add ( value , w   ) ; } }
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
        const std::shared_ptr<Result_t> m_result {}    ;
        /// size of m_slots 
        unsigned long                   m_N      { 1 } ;
        /// (current) results per  slot 
        std::vector<Result_t>           m_slots  { 1 } ;
        // ====================================================================
      } ; //                       The end of class ROOT::Detail::RDF::WStatVar
      // ======================================================================
      /** @class Moment_
       *  Helper class to get statitsics for the column in DataFrame 
       *  @see Ostap::Math::Moment_
       */
      template <unsigned short N>
      class Moment_ : public RActionImpl<Moment_<N> > 
      {
      public:
        // ====================================================================
        /// define the result type 
        using Result_t = Ostap::Math::Moment_<N> ;
        // ====================================================================
      public:
        // ====================================================================
        /// default constructor 
        Moment_ () 
          : m_result ( std::make_shared<Result_t>() ) 
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
          , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetThreadPoolSize     () ) : 1u )
#else 
          , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetImplicitMTPoolSize () ) : 1u )
#endif
          , m_slots  ( this->m_N ) 
        {}
        /// Move constructor 
        Moment_ (       Moment_&& ) = default ;
        /// Copy constructor is disabled 
        Moment_ ( const Moment_&  ) = delete ;
        // ====================================================================
      public:
        // ====================================================================
        /// initialize (empty) 
        void InitTask   ( TTreeReader * , unsigned int ) {} ;
        /// initialize (empty) 
        void Initialize () {} ;
        /// finalize : sum over the slots 
        void Finalize   () 
        { 
          Result_t sum { m_slots [ 0 ] } ;
          for ( unsigned int i = 1 ; i < m_N ; ++i ) { sum += m_slots [ i ] ; }
          *m_result = sum ;
        }
        /// who am I ?
        std::string GetActionName() { return "Moment_" ; }
        // ====================================================================
      public:
        // ====================================================================
        /// The basic method: increment the counter 
        void Exec ( unsigned int slot , double value ) 
        { m_slots [ slot % m_N ].add ( value ) ; } 
        // ====================================================================
        /// The basic method: increment the counter for the vector-like columns       
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif 
        void Exec ( unsigned int slot , const T &vs )
        { Result_t& m = m_slots [ slot % m_N ] ; for ( const auto & v : vs ) { m.add ( v ) ; } }
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
        const std::shared_ptr<Result_t> m_result {   } ;
        /// size of m_slots 
        unsigned long                   m_N      { 1 } ;
        /// (current) results per  slot 
        std::vector<Result_t>           m_slots  { 1 } ;
        // ====================================================================
      } ; //                        The end of class ROOT::Detail::RDF::Moment_
      // ======================================================================
      /** @class WMoment_
       *  Helper class to get statitsics for the column in DataFrame 
       *  @see Ostap::Math::WMoment_
       */
      template <unsigned short N>
      class WMoment_ : public RActionImpl<WMoment_<N> > 
      {
      public:
        // ====================================================================
        /// define the result type 
        using Result_t = Ostap::Math::WMoment_<N> ;
        // ====================================================================
      public:
        // ====================================================================
        /// default constructor 
        WMoment_ () 
          : m_result ( std::make_shared<Result_t>() ) 
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
          , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetThreadPoolSize     () ) : 1u )
#else 
          , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetImplicitMTPoolSize () ) : 1u )
#endif
          , m_slots  ( this->m_N ) 
        {}
        /// Move constructor 
        WMoment_ (       WMoment_&& ) = default ;
        /// Copy constructor is disabled 
        WMoment_ ( const WMoment_&  ) = delete ;
        // ====================================================================
      public:
        // ====================================================================
        /// initialize (empty) 
        void InitTask   ( TTreeReader * , unsigned int ) {} ;
        /// initialize (empty) 
        void Initialize () {} ;
        /// finalize : sum over the slots 
        void Finalize   () 
        { 
          Result_t sum { m_slots [ 0 ] } ;
          for ( unsigned int i = 1 ; i < m_N ; ++i ) { sum += m_slots [ i ] ; }
          *m_result = sum ;
        }
        /// who am I ?
        std::string GetActionName() { return "WMoment_" ; }
        // ====================================================================
      public:
        // ====================================================================
        /// The basic method: increment the counter 
        void Exec ( unsigned int slot , double value , double weight ) 
        { m_slots [ slot % m_N ].add ( value , weight ) ; } 
        // ====================================================================
        /// The basic method: increment the counter for the vector-like columns       
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif 
        void Exec ( unsigned int slot , const T &vs , double weight )
        { Result_t& m = m_slots [ slot % m_N ] ; for ( const auto & v : vs ) { m.add ( v ) ; } }
        // ====================================================================
        /// The basic method: increment the counter for the vector-like column of weight 
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif
        void Exec ( unsigned int slot , const double value , const T &ws )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto & w : ws ) { e.add ( value , w   ) ; } }
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
        const std::shared_ptr<Result_t> m_result {   } ;
        /// size of m_slots 
        unsigned long                   m_N      { 1 } ;
        /// (current) results per  slot 
        std::vector<Result_t>           m_slots  { 1 } ;
        // ====================================================================
      } ; //                        The end of class ROOT::Detail::RDF::Moment_
      // ======================================================================
      /** @class LegendrePoly
       *  Helper class to parameterise data as 1D Legendre polynomial
       *  @see Ostap::Math::LegendreSum 
       *  @see Ostap::DataFrame 
       */
      class LegendrePoly : public RActionImpl<LegendrePoly> 
      {
      public:
        // ====================================================================
        /// define the resutl type 
        using Result_t = Ostap::Math::LegendreSum;
        // ====================================================================
      public:
        // ====================================================================
        /// constructor 
        LegendrePoly 
        ( const unsigned short N    , 
          const double         xmin ,
          const double         xmax ) ;
        /// constructor from 
        LegendrePoly 
        ( const Ostap::Math::LegendreSum&  ) ;
        /// Move constructor 
        LegendrePoly (       LegendrePoly&& ) = default ;
        /// Copy constructor is disabled 
        LegendrePoly ( const LegendrePoly&  ) = delete ;
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
        std::string GetActionName() { return "LegendrePoly" ; }
        // ====================================================================
      public:
        // ====================================================================
        /// The basic method: increment the counter 
        void Exec ( unsigned int slot , double value , double weight = 1 ) 
        { m_slots [ slot % m_N ].Fill ( value , weight ) ; } 
        // ====================================================================
        /// The basic method: increment the counter for the vector-like column of values        
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif 
        void Exec ( unsigned int slot , const T &vs , const double weight = 1 )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto & v : vs ) { e.Fill ( v , weight ) ; } }
        // ====================================================================
        /// The basic method: increment the counter for the vector-like column of weight 
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif
        void Exec ( unsigned int slot , const double value , const T &ws )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto & w : ws ) { e.Fill ( value , w   ) ; } }
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
        const std::shared_ptr<Result_t> m_result {}    ;
        /// size of m_slots 
        unsigned long                   m_N      { 1 } ;
        /// (current) results per  slot 
        std::vector<Result_t>           m_slots  { 1 } ;
        // ====================================================================
      } ; //                   The end of class ROOT::Detail::RDF::LegendrePoly
      // ======================================================================
      /** @class ChebyshePoly
       *  Helper class to parameterise data as 1D Chebyshev polynomial
       *  @see Ostap::Math::ChebyshevSum 
       *  @see Ostap::DataFrame 
       */
      class ChebyshevPoly : public RActionImpl<ChebyshevPoly> 
      {
      public:
        // ====================================================================
        /// define the resutl type 
        using Result_t = Ostap::Math::ChebyshevSum;
        // ====================================================================
      public:
        // ====================================================================
        /// constructor
        ChebyshevPoly 
        ( const unsigned short N    , 
          const double         xmin ,
          const double         xmax ) ;
        /// constructor from 
        ChebyshevPoly 
        ( const Ostap::Math::ChebyshevSum&  ) ;
        /// Move constructor 
        ChebyshevPoly (       ChebyshevPoly&& ) = default ;
        /// Copy constructor is disabled 
        ChebyshevPoly ( const ChebyshevPoly&  ) = delete ;
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
        std::string GetActionName() { return "ChebyshevPoly" ; }
        // ====================================================================
      public:
        // ====================================================================
        /// The basic method: increment the counter 
        void Exec ( unsigned int slot , double value , double weight = 1 ) 
        { m_slots [ slot % m_N ].Fill ( value , weight ) ; } 
        // ====================================================================
        /// The basic method: increment the counter for the vector-like column of values        
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif 
        void Exec ( unsigned int slot , const T &vs , const double weight = 1 )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto & v : vs ) { e.Fill ( v , weight ) ; } }
        // ====================================================================
        /// The basic method: increment the counter for the vector-like column of weight 
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif
        void Exec ( unsigned int slot , const double value , const T &ws )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto & w : ws ) { e.Fill ( value , w   ) ; } }
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
        const std::shared_ptr<Result_t> m_result {}    ;
        /// size of m_slots 
        unsigned long                   m_N      { 1 } ;
        /// (current) results per  slot 
        std::vector<Result_t>           m_slots  { 1 } ;
        // ====================================================================
      } ; //                   The end of class ROOT::Detail::RDF::LegendrePoly 
      // ======================================================================
      /** @class BernsteinPoly
       *  Helper class to parameterise data as 1D Bernstein polynomial
       *  @see Ostap::Math::Bernstein 
       *  @see Ostap::DataFrame 
       */
      class BernsteinPoly : public RActionImpl<BernsteinPoly> 
      {
      public:
        // ====================================================================
        /// define the resutl type 
        using Result_t = Ostap::Math::Bernstein;
        // ====================================================================
      public:
        // ====================================================================
        /// constructor
        BernsteinPoly 
        ( const unsigned short N    , 
          const double         xmin ,
          const double         xmax ) ;
        /// constructor from 
        BernsteinPoly 
        ( const Ostap::Math::Bernstein&  ) ;
        /// Move constructor
        BernsteinPoly (       BernsteinPoly&& ) = default ;
        /// Copy constructor is disabled 
        BernsteinPoly ( const BernsteinPoly&  ) = delete ;
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
        std::string GetActionName() { return "BernsteinPoly" ; }
        // ====================================================================
      public:
        // ====================================================================
        /// The basic method: increment the counter 
        void Exec ( unsigned int slot , double value , double weight = 1 ) 
        { m_slots [ slot % m_N ].Fill ( value , weight ) ; } 
        // ====================================================================
        /// The basic method: increment the counter for the vector-like column of values        
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif 
        void Exec ( unsigned int slot , const T &vs , const double weight = 1 )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto & v : vs ) { e.Fill ( v , weight ) ; } }
        // ====================================================================
        /// The basic method: increment the counter for the vector-like column of weight 
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif
        void Exec ( unsigned int slot , const double value , const T &ws )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto & w : ws ) { e.Fill ( value , w   ) ; } }
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
        const std::shared_ptr<Result_t> m_result {}    ;
        /// size of m_slots 
        unsigned long                   m_N      { 1 } ;
        /// (current) results per  slot 
        std::vector<Result_t>           m_slots  { 1 } ;
        // ====================================================================
      } ; //                  The end of class ROOT::Detail::RDF::BernsteinPoly
      // ======================================================================
      /** @class LegendrePoly2
       *  Helper class to parameterise data as 2D Legendre polynomial
       *  @see Ostap::Math::LegendreSum2
       *  @see Ostap::DataFrame 
       */
      class LegendrePoly2 : public RActionImpl<LegendrePoly2> 
      {
      public:
        // ====================================================================
        /// define the resutl type 
        using Result_t = Ostap::Math::LegendreSum2;
        // ====================================================================
      public:
        // ====================================================================
        /// constructor 
        LegendrePoly2 
        ( const unsigned short NX   , 
          const unsigned short NY   , 
          const double         xmin ,
          const double         xmax ,
          const double         ymin ,
          const double         ymax ) ;
        /// constructor from 
        LegendrePoly2 
        ( const Ostap::Math::LegendreSum2&  ) ;
        /// Move constructor 
        LegendrePoly2 (       LegendrePoly2&& ) = default ;
        /// Copy constructor is disabled 
        LegendrePoly2 ( const LegendrePoly2&  ) = delete ;
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
        std::string GetActionName() { return "LegendrePoly2" ; }
        // ====================================================================
      public:
        // ====================================================================
        /// The basic method: increment the counter 
        void Exec ( unsigned int slot , double x , double y  , double weight = 1 ) 
        { m_slots [ slot % m_N ].Fill ( x , y , weight ) ; } 
        // ====================================================================
        /// The basic method: increment the counter for the vector-like column of values
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif 
        void Exec ( unsigned int slot , const T &xs , const double y , const double weight = 1 )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto &x : xs ) { e.Fill ( x , y , weight ) ; } }
        // ====================================================================
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif 
        void Exec ( unsigned int slot , const double x , const T &ys , const double weight = 1 )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto &y : ys ) { e.Fill ( x , y , weight ) ; } }
        // ====================================================================
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif
        void Exec ( unsigned int slot , const double x, const double y , const T &ws )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto & w : ws ) { e.Fill ( x , y , w   ) ; } }
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
        const std::shared_ptr<Result_t> m_result {}    ;
        /// size of m_slots 
        unsigned long                   m_N      { 1 } ;
        /// (current) results per  slot 
        std::vector<Result_t>           m_slots  { 1 } ;
        // ====================================================================
      } ; //                  The end of class ROOT::Detail::RDF::LegendrePoly2
      // ======================================================================
      /** @class BernsteinPoly2
       *  Helper class to parameterise data as 2D Bernstein polynomial
       *  @see Ostap::Math::Bernstein2D
       *  @see Ostap::DataFrame 
       */
      class BernsteinPoly2 : public RActionImpl<BernsteinPoly2> 
      {
      public:
        // ====================================================================
        /// define the resutl type 
        using Result_t = Ostap::Math::Bernstein2D;
        // ====================================================================
      public:
        // ====================================================================
        /// constructor 
        BernsteinPoly2 
        ( const unsigned short NX   , 
          const unsigned short NY   , 
          const double         xmin ,
          const double         xmax ,
          const double         ymin ,
          const double         ymax ) ;
        /// constructor from 
        BernsteinPoly2 
        ( const Ostap::Math::Bernstein2D&  ) ;
        /// Move constructor 
        BernsteinPoly2 (       BernsteinPoly2&& ) = default ;
        /// Copy constructor is disabled 
        BernsteinPoly2 ( const BernsteinPoly2&  ) = delete ;
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
        std::string GetActionName() { return "BernsteinPoly2" ; }
        // ====================================================================
      public:
        // ====================================================================
        /// The basic method: increment the counter 
        void Exec ( unsigned int slot , double x , double y  , double weight = 1 ) 
        { m_slots [ slot % m_N ].Fill ( x , y , weight ) ; } 
        // ====================================================================
        /// The basic method: increment the counter for the vector-like column of values
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif 
        void Exec ( unsigned int slot , const T &xs , const double y , const double weight = 1 )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto &x : xs ) { e.Fill ( x , y , weight ) ; } }
        // ====================================================================
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif 
        void Exec ( unsigned int slot , const double x , const T &ys , const double weight = 1 )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto &y : ys ) { e.Fill ( x , y , weight ) ; } }
        // ====================================================================
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif
        void Exec ( unsigned int slot , const double x, const double y , const T &ws )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto & w : ws ) { e.Fill ( x , y , w   ) ; } }
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
        const std::shared_ptr<Result_t> m_result {}    ;
        /// size of m_slots 
        unsigned long                   m_N      { 1 } ;
        /// (current) results per  slot 
        std::vector<Result_t>           m_slots  { 1 } ;
        // ====================================================================
      } ; //                The end of class ROOT::Detail::RDF::BernnsteinPoly2
      // ======================================================================
      /** @class LegendrePoly3
       *  Helper class to parameterise data as 3D Legendre polynomial
       *  @see Ostap::Math::LegendreSum3 
       *  @see Ostap::DataFrame 
       */
      class LegendrePoly3 : public RActionImpl<LegendrePoly3> 
      {
      public:
        // ====================================================================
        /// define the resutl type 
        using Result_t = Ostap::Math::LegendreSum3;
        // ====================================================================
      public:
        // ====================================================================
        /// constructor 
        LegendrePoly3 
        ( const unsigned short NX   , 
          const unsigned short NY   , 
          const unsigned short NZ   , 
          const double         xmin ,
          const double         xmax ,
          const double         ymin ,
          const double         ymax ,
          const double         zmin ,
          const double         zmax ) ;
        /// constructor from 
        LegendrePoly3 
        ( const Ostap::Math::LegendreSum3&  ) ;
        /// Move constructor 
        LegendrePoly3 (       LegendrePoly3&& ) = default ;
        /// Copy constructor is disabled 
        LegendrePoly3 ( const LegendrePoly3&  ) = delete ;
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
        std::string GetActionName() { return "LegendrePoly3" ; }
        // ====================================================================
      public:
        // ====================================================================
        /// The basic method: increment the counter 
        void Exec ( unsigned int slot , double x , double y , double z , double weight = 1 ) 
        { m_slots [ slot % m_N ].Fill ( x , y , z , weight ) ; } 
        // ====================================================================
        /// The basic method: increment the counter for the vector-like column of values
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif 
        void Exec ( unsigned int slot , const T &xs , const double y , const double z , const double weight = 1 )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto &x : xs ) { e.Fill ( x , y , z , weight ) ; } }
        // ====================================================================
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif 
        void Exec ( unsigned int slot , const double x , const T &ys , const double z , const double weight = 1 )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto &y : ys ) { e.Fill ( x , y , z , weight ) ; } }
        // ====================================================================
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif 
        void Exec ( unsigned int slot , const double x , const double y , const T &zs , const double weight = 1 )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto &z : zs ) { e.Fill ( x , y , z , weight ) ; } }
        // ====================================================================
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif
        void Exec ( unsigned int slot , const double x, const double y , const double z , const T &ws )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto & w : ws ) { e.Fill ( x , y , z , w   ) ; } }
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
        const std::shared_ptr<Result_t> m_result {}    ;
        /// size of m_slots 
        unsigned long                   m_N      { 1 } ;
        /// (current) results per  slot 
        std::vector<Result_t>           m_slots  { 1 } ;
        // ====================================================================
      } ; //                  The end of class ROOT::Detail::RDF::LegendrePoly3
      // ======================================================================
      /** @class BernsteinPoly3
       *  Helper class to parameterise data as 3D Bernstein polynomial
       *  @see Ostap::Math::Bernstein3D 
       *  @see Ostap::DataFrame 
       */
      class BernsteinPoly3 : public RActionImpl<BernsteinPoly3> 
      {
      public:
        // ====================================================================
        /// define the resutl type 
        using Result_t = Ostap::Math::Bernstein3D;
        // ====================================================================
      public:
        // ====================================================================
        /// constructor 
        BernsteinPoly3 
        ( const unsigned short NX   , 
          const unsigned short NY   , 
          const unsigned short NZ   , 
          const double         xmin ,
          const double         xmax ,
          const double         ymin ,
          const double         ymax ,
          const double         zmin ,
          const double         zmax ) ;
        /// constructor from 
        BernsteinPoly3 
        ( const Ostap::Math::Bernstein3D&  ) ;
        /// Move constructor 
        BernsteinPoly3 (       BernsteinPoly3&& ) = default ;
        /// Copy constructor is disabled 
        BernsteinPoly3 ( const BernsteinPoly3&  ) = delete ;
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
        std::string GetActionName() { return "BernsteinPoly3" ; }
        // ====================================================================
      public:
        // ====================================================================
        /// The basic method: increment the counter 
        void Exec ( unsigned int slot , double x , double y , double z , double weight = 1 ) 
        { m_slots [ slot % m_N ].Fill ( x , y , z , weight ) ; } 
        // ====================================================================
        /// The basic method: increment the counter for the vector-like column of values
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif 
        void Exec ( unsigned int slot , const T &xs , const double y , const double z , const double weight = 1 )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto &x : xs ) { e.Fill ( x , y , z , weight ) ; } }
        // ====================================================================
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif 
        void Exec ( unsigned int slot , const double x , const T &ys , const double z , const double weight = 1 )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto &y : ys ) { e.Fill ( x , y , z , weight ) ; } }
        // ====================================================================
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif 
        void Exec ( unsigned int slot , const double x , const double y , const T &zs , const double weight = 1 )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto &z : zs ) { e.Fill ( x , y , z , weight ) ; } }
        // ====================================================================
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif
        void Exec ( unsigned int slot , const double x, const double y , const double z , const T &ws )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto & w : ws ) { e.Fill ( x , y , z , w   ) ; } }
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
        const std::shared_ptr<Result_t> m_result {}    ;
        /// size of m_slots 
        unsigned long                   m_N      { 1 } ;
        /// (current) results per  slot 
        std::vector<Result_t>           m_slots  { 1 } ;
        // ====================================================================
      } ; //                 The end of class ROOT::Detail::RDF::BernsteinPoly3
      // ======================================================================
      /** @class LegendrePoly4
       *  Helper class to parameterise data as 4D Legendre polynomial
       *  @see Ostap::Math::LegendreSum4
       *  @see Ostap::DataFrame 
       */
      class LegendrePoly4 : public RActionImpl<LegendrePoly4> 
      {
      public:
        // ====================================================================
        /// define the resutl type 
        using Result_t = Ostap::Math::LegendreSum4;
        // ====================================================================
      public:
        // ====================================================================
        /// constructor 
        LegendrePoly4 
        ( const unsigned short NX   , 
          const unsigned short NY   , 
          const unsigned short NZ   , 
          const unsigned short NU   , 
          const double         xmin ,
          const double         xmax ,
          const double         ymin ,
          const double         ymax ,
          const double         zmin ,
          const double         zmax ,
          const double         umin ,
          const double         umax ) ;
        /// constructor from 
        LegendrePoly4 
        ( const Ostap::Math::LegendreSum4&  ) ;
        /// Move constructor 
        LegendrePoly4 (       LegendrePoly4&& ) = default ;
        /// Copy constructor is disabled 
        LegendrePoly4 ( const LegendrePoly4&  ) = delete ;
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
        std::string GetActionName() { return "LegendrePoly4" ; }
        // ====================================================================
      public:
        // ====================================================================
        /// The basic method: increment the counter 
        void Exec ( unsigned int slot , double x , double y , double z , double u , double weight = 1 ) 
        { m_slots [ slot % m_N ].Fill ( x , y , z , u , weight ) ; } 
        // ====================================================================
        /// The basic method: increment the counter for the vector-like column of values
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif 
        void Exec ( unsigned int slot , const T &xs , const double y , const double z , const double u , const double weight = 1 )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto &x : xs ) { e.Fill ( x , y , z , u , weight ) ; } }
        // ====================================================================
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif 
        void Exec ( unsigned int slot , const double x , const T &ys , const double z , const double u, const double weight = 1 )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto &y : ys ) { e.Fill ( x , y , z , u , weight ) ; } }

        // ====================================================================
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif 
        void Exec ( unsigned int slot , const double x , const double y , const T &zs , const double u , const double weight = 1 )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto &z : zs ) { e.Fill ( x , y , z , u , weight ) ; } }
        // ====================================================================
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif 
        void Exec ( unsigned int slot , const double x , const double y , const double z , const T &us , const double weight = 1 )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto &u : us ) { e.Fill ( x , y , z , u   , weight ) ; } }
        // ====================================================================
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif
        void Exec ( unsigned int slot , const double x, const double y , const double z , const double u , const T &ws )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto & w : ws ) { e.Fill ( x , y , z , w   ) ; } }
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
        const std::shared_ptr<Result_t> m_result {}    ;
        /// size of m_slots 
        unsigned long                   m_N      { 1 } ;
        /// (current) results per  slot 
        std::vector<Result_t>           m_slots  { 1 } ;
        // ====================================================================
      } ; //                  The end of class ROOT::Detail::RDF::LegendrePoly3
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
    using StatVar        = ROOT::Detail::RDF::StatVar        ;
    using WStatVar       = ROOT::Detail::RDF::WStatVar       ;
    //  
    using LegendrePoly   = ROOT::Detail::RDF::LegendrePoly   ;
    using ChebyshevPoly  = ROOT::Detail::RDF::ChebyshevPoly  ;
    using BernsteinPoly  = ROOT::Detail::RDF::BernsteinPoly  ;
    //
    using LegendrePoly2  = ROOT::Detail::RDF::LegendrePoly2  ;
    using BernsteinPoly2 = ROOT::Detail::RDF::BernsteinPoly2 ;
    //
    using LegendrePoly3  = ROOT::Detail::RDF::LegendrePoly3  ;
    using BernsteinPoly3 = ROOT::Detail::RDF::BernsteinPoly3 ;
    //
    using LegendrePoly4  = ROOT::Detail::RDF::LegendrePoly4  ;
    //
    template <unsigned short N> 
    using Moment_        = ROOT::Detail::RDF::Moment_<N>     ;
    template <unsigned short N> 
    using WMoment_       = ROOT::Detail::RDF::WMoment_<N>    ;
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
