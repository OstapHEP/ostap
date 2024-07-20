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
      /** @class StatAction
       *  Helper class to get statitsics for the column in DataFrame 
       *  using soem COUNTER class
       *  Requirements for COUNTER:
       *  -  counter.add ( value ) 
       *  -  counter += counter 
       *  @see Ostap::StatEntity
       *  @see Ostap::Math::Moment_
       */
      template <class COUNTER> 
      class StatAction : public RActionImpl< StatAction<COUNTER> > 
      {
      public:
        // ====================================================================
        /// define the result type 
        using Result_t = COUNTER ;
        // ====================================================================
      public:
        // ====================================================================
        /// constructor 
	template <typename ...Args>
        StatAction ( const Args& ...args ) 
          : m_result ( std::make_shared<Result_t>( args ... ) ) 
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
          , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetThreadPoolSize     () ) : 1u )
#else 
          , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetImplicitMTPoolSize () ) : 1u )
#endif
	  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
        {}
        /// Move constructor 
        StatAction (       StatAction&& ) = default ;
        /// Copy constructor is disabled 
        StatAction ( const StatAction&  ) = delete ;
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
        std::string GetActionName() { return "StatAction" ; }
        // ====================================================================
      public:
        // ====================================================================
        /// The basic method: increment the counter 
        void Exec ( unsigned int slot , const double value ) 
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
        const std::shared_ptr<Result_t> m_result {} ;
        /// size of m_slots 
        unsigned long                   m_N      {} ;
        /// (current) results per  slot 
        std::vector<Result_t>           m_slots  {} ;
        // ====================================================================
      } ; //                    The end of class ROOT::Detail::RDF::StatCounter
      // ======================================================================
      /** @class WStatAction
       *  Helper class to get statitsics for the column in DataFrame 
       *  using soem COUNTER class
       *  Requirements for COUNTER:
       *  -  counter.add ( value , weight ) 
       *  -  counter += counter 
       *  @see Ostap::WStatEntity
       *  @see Ostap::Math::WMoment_
       */
      template <class Counter>
      class WStatAction : public RActionImpl< WStatAction<Counter> > 
      {
      public:
        // ====================================================================
        /// define the resutl type 
        using Result_t = Counter ;
        // ====================================================================
      public:
        // ====================================================================
        /// default constructor 
	template <typename ...Args>
        WStatAction ( const Args& ...args ) 
	  : m_result ( std::make_shared<Result_t>( args ... ) ) 
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
	  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetThreadPoolSize     () ) : 1u )
#else 
	  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetImplicitMTPoolSize () ) : 1u )
#endif
	  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
        {}
	/// Move constructor 
        WStatAction  (       WStatAction&& ) = default ;
        /// Copy constructor is disabled 
        WStatAction  ( const WStatAction&   ) = delete ;
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
        std::string GetActionName() { return "WStatAction" ; }
        // ====================================================================
      public:
        // ====================================================================
        /// The basic method: increment the counter 
        void Exec ( unsigned int slot , const double value , const double weight ) 
        { m_slots [ slot % m_N ].add ( value , weight ) ; } 
        // ====================================================================
        /// The basic method: increment the counter for the vector-like columns       
#if ROOT_VERSION(6,22,0) <= ROOT_VERSION_CODE
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
#else 
        template <typename T, typename std::enable_if<IsContainer<T>::value, int>::type = 0>
#endif 
        void Exec ( unsigned int slot , const T &vs , const double weight = 1 )
        { Result_t& m = m_slots [ slot % m_N ] ; for ( const auto & v : vs ) { m.add ( v , weight ) ; } }
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
      } ; //                    The end of class ROOT::Detail::RDF::WStatAction
      // ======================================================================
      /** @class Poly1Action
       *  Helper class to parameterise data as 1D polynomial
       *  @see Ostap::Math::LegendreSum 
       *  @see Ostap::Math::ChebyshevSum 
       *  @see Ostap::Math::BernsteinSum 
       *  @see Ostap::DataFrame 
       *  Requirements for class POLYNOMIAL
       *   - copy constructor 
       *   - summation : <code>polynomial += polynomial </code>
       *   - scaling   : <code>polynomial *= scale </code>
       *   - fill:       <code> polynomial.Fill ( value , weight ) </code>
       */
      template <class POLYNOMIAL>
      class Poly1Action : public RActionImpl<Poly1Action<POLYNOMIAL> > 
      {
      public:
        // ====================================================================
        /// define the resutl type 
        using Result_t = POLYNOMIAL ;
        // ====================================================================
      public:
        // ====================================================================
        /// constructor
	template <typename ...Args>
        Poly1Action 
        ( const Args& ...args )
	  : m_result ( std::make_shared<Result_t>( args ... ) ) 
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
	  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetThreadPoolSize     () ) : 1u )
#else 
	  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetImplicitMTPoolSize () ) : 1u )
#endif
	  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
	{
	  (*m_result) *= 0.0 ; // reset the polynomial
	  for ( auto& i : m_slots ) { i *= 0.0 ; }
	}  
        /// Move constructor 
        Poly1Action (       Poly1Action&& ) = default ;
        /// Copy constructor is disabled 
        Poly1Action ( const Poly1Action&  ) = delete ;
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
	  (*m_result) *= 0.0 ;
	  for ( unsigned int i = 0 ; i < m_N ; ++i ) { *m_result += m_slots [ i ] ; }
	}
        /// who am I ?
        std::string GetActionName() { return "Poly1Action" ; }
        // ====================================================================
      public:
        // ====================================================================
        /// The basic method: increment the counter 
        void Exec ( unsigned int slot       ,
		    const double value      ,
		    const double weight = 1 ) 
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
      } ; //                    The end of class ROOT::Detail::RDF::Poly1Action
      // ======================================================================
      /** @class Poly2Action
       *  Helper class to parameterise data as 2D polynomial 
       *  @see Ostap::Math::LegendreSum2
       *  @see Ostap::Math::BErnstein2D 
       *  @see Ostap::DataFrame 
       *  Requirements for class POLYNOMIAL
       *   - copy constructor 
       *   - summation : <code>polynomial += polynomial </code>
       *   - scaling   : <code>polynomial *= scale </code>
       *   - fill:       <code> polynomial.Fill ( x , y , weight ) </code>
       */
      template <class POLYNOMIAL>
      class Poly2Action : public RActionImpl<Poly2Action<POLYNOMIAL>> 
      {
      public:
        // ====================================================================
        /// define the resutl type 
        using Result_t = POLYNOMIAL ;
	  // ====================================================================
      public:
        // ====================================================================
        /// constructor 
	template <typename ...Args>
        Poly2Action 
        ( const Args& ...args )
	  : m_result ( std::make_shared<Result_t>( args ... ) ) 
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
	  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetThreadPoolSize     () ) : 1u )
#else 
	  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetImplicitMTPoolSize () ) : 1u )
#endif
	  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
	{
	  (*m_result) *= 0.0 ; // reset the polynomial
	  for ( auto& i : m_slots ) { i *= 0.0 ; }
	}  
        /// Move constructor 
        Poly2Action (       Poly2Action&& ) = default ;
        /// Copy constructor is disabled 
        Poly2Action ( const Poly2Action&  ) = delete ;
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
	  (*m_result) *= 0.0 ;
	  for ( unsigned int i = 0 ; i < m_N ; ++i ) { *m_result += m_slots [ i ] ; }
	}
        /// who am I ?
        std::string GetActionName() { return "Poly2Action" ; }
        // ====================================================================
      public:
        // ====================================================================
        /// The basic method: increment the counter 
        void Exec ( unsigned int slot , const double x , const double y  , double weight = 1 ) 
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
      } ; //                    The end of class ROOT::Detail::RDF::Poly2Action
      // ======================================================================
      /** @class Poly3Action
       *  Helper class to parameterise data as 3D polynomial 
       *  @see Ostap::Math::LegendreSum3
       *  @see Ostap::Math::Bernstein3D
       *  @see Ostap::DataFrame 
       *  Requirements for class POLYNOMIAL
       *   - copy constructor 
       *   - summation : <code>polynomial += polynomial </code>
       *   - scaling   : <code>polynomial *= scale </code>
       *   - fill:       <code> polynomial.Fill ( x , y , z , weight ) </code>
       */
      template <class POLYNOMIAL>
      class Poly3Action : public RActionImpl<Poly3Action<POLYNOMIAL>> 
      {
      public:
        // ====================================================================
        /// define the resutl type 
        using Result_t = POLYNOMIAL ;
	  // ====================================================================
      public:
        // ====================================================================
	template <typename ...Args>
        Poly3Action 
        ( const Args& ...args )
	  : m_result ( std::make_shared<Result_t>( args... ) ) 
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
	  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetThreadPoolSize     () ) : 1u )
#else 
	  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetImplicitMTPoolSize () ) : 1u )
#endif
	  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
	{
	  (*m_result) *= 0.0 ; // reset the polynomial
	  for ( auto& i : m_slots ) { i *= 0.0 ; }
	}  
        /// Move constructor 
        Poly3Action (       Poly3Action&& ) = default ;
        /// Copy constructor is disabled 
        Poly3Action ( const Poly3Action&  ) = delete ;
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
	  (*m_result) *= 0.0 ;
	  for ( unsigned int i = 0 ; i < m_N ; ++i ) { *m_result += m_slots [ i ] ; }
	}
        /// who am I ?
        std::string GetActionName() { return "Poly3Action" ; }
        // ====================================================================
      public:
        // ====================================================================
        /// The basic method: increment the counter 
        void Exec ( unsigned int slot ,
		    const double x          ,
		    const double y          ,
		    const double z          ,
		    const double weight = 1 ) 
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
        // =====================================================================
      } ; //                     The end of class ROOT::Detail::RDF::Poly3Action
      // =======================================================================
      /** @class Poly4Action
       *  Helper class to parameterise data as 4D polynomial 
       *  @see Ostap::Math::LegendreSum4
       *  @see Ostap::DataFrame 
       *  Requirements for class POLYNOMIAL
       *   - copy constructor 
       *   - summation : <code>polynomial += polynomial </code>
       *   - scaling   : <code>polynomial *= scale </code>
       *   - fill:       <code> polynomial.Fill (  x, y , z ,v  , weight ) </code>
       */
      template <class POLYNOMIAL>
      class Poly4Action : public RActionImpl<Poly4Action<POLYNOMIAL>> 
      {
      public:
        // ====================================================================
        /// define the resutl type 
        using Result_t = POLYNOMIAL ;
	  // ====================================================================
      public:
        // ====================================================================
	template <typename ...Args>
        Poly4Action 
        ( const Args& ...args )
	  : m_result ( std::make_shared<Result_t>( args ... ) ) 
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
	  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetThreadPoolSize     () ) : 1u )
#else 
	  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetImplicitMTPoolSize () ) : 1u )
#endif
	  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
	{
	  (*m_result) *= 0.0 ; // reset the polynomial
	  for ( auto& i : m_slots ) { i *= 0.0 ; }
	}  
        /// Move constructor 
        Poly4Action (       Poly4Action&& ) = default ;
        /// Copy constructor is disabled 
        Poly4Action ( const Poly4Action&  ) = delete ;
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
	  (*m_result) *= 0.0 ;
	  for ( unsigned int i = 0 ; i < m_N ; ++i ) { *m_result += m_slots [ i ] ; }
	}
        /// who am I ?
        std::string GetActionName() { return "Poly4Action" ; }
        // ====================================================================
      public:
        // ====================================================================
        /// The basic method: increment the counter 
        void Exec ( unsigned int slot       ,
		    const double x          ,
		    const double y          ,
		    const double z          ,
		    const double u          ,
		    const double weight = 1 ) 
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
      } ; //                    The end of class ROOT::Detail::RDF::Poly4Action
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
    
    template <class COUNTER>
    using StatAction  = ROOT::Detail::RDF::StatAction<COUNTER>   ;    
    template <class COUNTER>
    using WStatAction = ROOT::Detail::RDF::WStatAction<COUNTER>  ;

    using StatVar  = StatAction <Ostap::StatEntity>  ;    
    using WStatVar = WStatAction<Ostap::WStatEntity> ;

    template <unsigned short N> 
    using Moment_  = StatAction <typename Ostap::Math::Moment_<N>  > ;
    template <unsigned short N> 
    using WMoment_ =  WStatAction<typename Ostap::Math::WMoment_<N> > ;

    using GeometricMean   = StatAction<Ostap::Math::GeometricMean>    ;    
    using ArithmeticMean  = StatAction<Ostap::Math::ArithmeticMean>   ;    
    using HarmonicMean    = StatAction<Ostap::Math::HarmonicMean>     ;    
    using PowerMean       = StatAction<Ostap::Math::PowerMean>        ;    
    using LehmerMean      = StatAction<Ostap::Math::LehmerMean>       ;    
    
    using WGeometricMean  = WStatAction<Ostap::Math::WGeometricMean>  ;    
    using WArithmeticMean = WStatAction<Ostap::Math::WArithmeticMean> ;    
    using WHarmonicMean   = WStatAction<Ostap::Math::WHarmonicMean>   ;    
    using WPowerMean      = WStatAction<Ostap::Math::WPowerMean>      ;    
    using WLehmerMean     = WStatAction<Ostap::Math::WLehmerMean>     ;    
    
    template <class POLYNOMIAL>
    using Poly1Action    = ROOT::Detail::RDF::Poly1Action<POLYNOMIAL> ;    
    template <class POLYNOMIAL>
    using Poly2Action    = ROOT::Detail::RDF::Poly2Action<POLYNOMIAL> ;    
    template <class POLYNOMIAL>
    using Poly3Action    = ROOT::Detail::RDF::Poly4Action<POLYNOMIAL> ;    
    template <class POLYNOMIAL>
    using Poly4Action    = ROOT::Detail::RDF::Poly4Action<POLYNOMIAL> ;    
    
    using LegendrePoly   = Poly1Action<Ostap::Math::LegendreSum>  ;
    using ChebyshevPoly  = Poly1Action<Ostap::Math::ChebyshevSum> ;
    using BernsteinPoly  = Poly1Action<Ostap::Math::Bernstein>    ;

    using LegendrePoly2  = Poly2Action <Ostap::Math::LegendreSum2> ;
    using BernsteinPoly2 = Poly2Action <Ostap::Math::Bernstein2D>  ; 

    using LegendrePoly3  = Poly3Action <Ostap::Math::LegendreSum3> ;
    using BernsteinPoly3 = Poly3Action <Ostap::Math::Bernstein3D>  ; 

    using LegendrePoly4  = Poly4Action <Ostap::Math::LegendreSum4> ;
   
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
