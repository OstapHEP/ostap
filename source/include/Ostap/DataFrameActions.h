// ============================================================================
#ifndef OSTAP_DATAFRAMEACTIONS_H 
#define OSTAP_DATAFRAMEACTIONS_H 1
// ============================================================================
// Inclodue files
// ============================================================================
// ROOT 
// ============================================================================
#include "RtypesCore.h" // ROOT 
// ============================================================================
#include "ROOT/RDataFrame.hxx"  // ROOT 
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
#include "Ostap/ECDF.h"
// ============================================================================
#include "Ostap/DataFrameUtils.h"
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
      /** @class Stat1Action
       *  Helper class to get statitsics for the column in DataFrame 
       *  using soem COUNTER class
       *  Requirements for COUNTER:
       *  -  counter.add ( value ) 
       *  -  counter += counter 
       *  -  default constructor 
       *  @see Ostap::StatEntity
       *  @see Ostap::Math::Moment_
       */
      template <class COUNTER> 
      class Stat1Action : public RActionImpl< Stat1Action<COUNTER> > 
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
        Stat1Action ( const Args& ...args ) 
          : m_result ( std::make_shared<Result_t>( args ... ) ) 
          , m_N      ( std::max ( 1u , Ostap::Utils::mt_pool_size () ) )
	  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
        {}
        /// Move constructor 
        Stat1Action (       Stat1Action&& ) = default ;
        /// Copy constructor is disabled 
        Stat1Action ( const Stat1Action&  ) = delete ;
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
        std::string GetActionName() { return "Stat1Action" ; }
        // ====================================================================
      public:
        // ====================================================================
        /// The basic method: increment the counter 
        void Exec ( unsigned int slot , const double value ) 
        { m_slots [ slot % m_N ].add ( value ) ; } 
        // ====================================================================
        /// The basic method: increment the counter for the vector-like columns       
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
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
      /** @class Stat2Action
       *  Helper class to get statitsics for the column in DataFrame 
       *  using soem COUNTER class
       *  Requirements for COUNTER:
       *  -  counter.add ( value , weight ) 
       *  -  counter += counter 
       *  -  defaut contructor 
       *  @see Ostap::WStatEntity
       *  @see Ostap::Math::Covariance
       *  @see Ostap::Math::WMoment_
       */
      template <class Counter>
      class Stat2Action : public RActionImpl< Stat2Action<Counter> > 
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
        Stat2Action ( const Args& ...args ) 
	  : m_result ( std::make_shared<Result_t>( args ... ) ) 
          , m_N      ( std::max ( 1u , Ostap::Utils::mt_pool_size () ) )
	  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
        {}
	/// Move constructor 
        Stat2Action  (       Stat2Action&& ) = default ;
        /// Copy constructor is disabled 
        Stat2Action  ( const Stat2Action&   ) = delete ;
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
        std::string GetActionName() { return "Stat2Action" ; }
        // ====================================================================
      public:
        // ====================================================================
        /// The basic method: increment the counter 
        void Exec ( unsigned int slot , const double value , const double weight ) 
        { m_slots [ slot % m_N ].add ( value , weight ) ; } 
        // ====================================================================
        /// The basic method: increment the counter for the vector-like columns       
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
        void Exec ( unsigned int slot , const T &vs , const double weight = 1 )
        { Result_t& m = m_slots [ slot % m_N ] ; for ( const auto & v : vs ) { m.add ( v , weight ) ; } }
        // ====================================================================
        /// The basic method: increment the counter for the vector-like column of weight 
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
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
      } ; //                    The end of class ROOT::Detail::RDF::Stat2Action
      // ======================================================================      
      /** @class Stat3Action
       *  Helper class to get statitsics for the column in DataFrame 
       *  using some COUNTER class
       *  Requirements for COUNTER:
       *  -  counter.add ( v1 , v2 , weight ) 
       *  -  counter += counter 
       *  @see Ostap::Math::WCovariance
       */
      template <class Counter>
      class Stat3Action : public RActionImpl< Stat3Action<Counter> > 
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
        Stat3Action ( const Args& ...args ) 
	  : m_result ( std::make_shared<Result_t>( args ... ) ) 
          , m_N      ( std::max ( 1u , Ostap::Utils::mt_pool_size () ) )
	  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
        {}
	/// Move constructor 
        Stat3Action  (       Stat3Action&& ) = default ;
        /// Copy constructor is disabled 
        Stat3Action  ( const Stat3Action&   ) = delete ;
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
        std::string GetActionName() { return "Stat3Action" ; }
        // ====================================================================
      public:
        // ====================================================================
        /// The basic method: increment the counter 
        void Exec
	( unsigned int slot   ,
	  const double v1     ,
	  const double v2     ,
	  const double weight ) 
        { m_slots [ slot % m_N ].add ( v1 , v2 , weight ) ; } 
        // ====================================================================
        /// The basic method: increment the counter for the vector-like columns       
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
        void Exec
	( unsigned int slot       ,
	  const T&     vs         ,
	  const double v2         ,
	  const double weight = 1 )
        { Result_t& m = m_slots [ slot % m_N ] ; for ( const auto& v1 : vs ) { m.add ( v1 , v2  , weight ) ; } }
        // ====================================================================
        /// The basic method: increment the counter for the vector-like columns       
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
        void Exec
	( unsigned int slot       ,
	  const double v1         ,
	  const T&     vs         ,
	  const double weight = 1 )
        { Result_t& m = m_slots [ slot % m_N ] ; for ( const auto& v2 : vs ) { m.add ( v1 , v2  , weight ) ; } }
        // ====================================================================
        /// The basic method: increment the counter for the vector-like column of weight 
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
        void Exec
	( unsigned int slot ,
	  const double v1   ,
	  const double v2   ,		    
	  const T&     ws   )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto& w : ws ) { e.add ( v1 , v2  , w   ) ; } }
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
      } ; //                    The end of class ROOT::Detail::RDF::Stat2Action
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
          , m_N      ( std::max ( 1u , Ostap::Utils::mt_pool_size () ) )
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
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
        void Exec ( unsigned int slot , const T &vs , const double weight = 1 )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto & v : vs ) { e.Fill ( v , weight ) ; } }
        // ====================================================================
        /// The basic method: increment the counter for the vector-like column of weight 
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
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
          , m_N      ( std::max ( 1u , Ostap::Utils::mt_pool_size () ) )
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
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
        void Exec ( unsigned int slot , const T &xs , const double y , const double weight = 1 )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto &x : xs ) { e.Fill ( x , y , weight ) ; } }
        // ====================================================================
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
        void Exec ( unsigned int slot , const double x , const T &ys , const double weight = 1 )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto &y : ys ) { e.Fill ( x , y , weight ) ; } }
        // ====================================================================
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
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
          , m_N      ( std::max ( 1u , Ostap::Utils::mt_pool_size () ) )
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
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
        void Exec ( unsigned int slot , const T &xs , const double y , const double z , const double weight = 1 )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto &x : xs ) { e.Fill ( x , y , z , weight ) ; } }
        // ====================================================================
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
        void Exec ( unsigned int slot , const double x , const T &ys , const double z , const double weight = 1 )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto &y : ys ) { e.Fill ( x , y , z , weight ) ; } }
        // ====================================================================
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
        void Exec ( unsigned int slot , const double x , const double y , const T &zs , const double weight = 1 )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto &z : zs ) { e.Fill ( x , y , z , weight ) ; } }
        // ====================================================================
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
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
          , m_N      ( std::max ( 1u , Ostap::Utils::mt_pool_size () ) )
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
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
        void Exec ( unsigned int slot , const T &xs , const double y , const double z , const double u , const double weight = 1 )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto &x : xs ) { e.Fill ( x , y , z , u , weight ) ; } }
        // ====================================================================
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
        void Exec ( unsigned int slot , const double x , const T &ys , const double z , const double u, const double weight = 1 )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto &y : ys ) { e.Fill ( x , y , z , u , weight ) ; } }
        // ====================================================================
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
        void Exec ( unsigned int slot , const double x , const double y , const T &zs , const double u , const double weight = 1 )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto &z : zs ) { e.Fill ( x , y , z , u , weight ) ; } }
        // ====================================================================
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
        void Exec ( unsigned int slot , const double x , const double y , const double z , const T &us , const double weight = 1 )
        { Result_t& e = m_slots [ slot % m_N ] ; for ( const auto &u : us ) { e.Fill ( x , y , z , u   , weight ) ; } }
        // ====================================================================
        template <typename T, typename std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value, int>::type = 0>
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
    using Stat1Action  = ROOT::Detail::RDF::Stat1Action<COUNTER>   ;    
    template <class COUNTER>
    using Stat2Action = ROOT::Detail::RDF::Stat2Action<COUNTER>  ;

    using StatVar  = Stat1Action <Ostap::StatEntity> ;    
    using WStatVar = Stat2Action<Ostap::WStatEntity> ;

    using ECDF     = Stat1Action <Ostap::Math::ECDF> ;    
    using WECDF    = Stat2Action<Ostap::Math::WECDF> ;
    
    template <unsigned short N> 
    using Moment_  = Stat1Action <typename Ostap::Math::Moment_<N>  > ;
    template <unsigned short N> 
    using WMoment_ =  Stat2Action<typename Ostap::Math::WMoment_<N> > ;

    using GeometricMean   = Stat1Action<Ostap::Math::GeometricMean>    ;    
    using ArithmeticMean  = Stat1Action<Ostap::Math::ArithmeticMean>   ;    
    using HarmonicMean    = Stat1Action<Ostap::Math::HarmonicMean>     ;    
    using PowerMean       = Stat1Action<Ostap::Math::PowerMean>        ;    
    using LehmerMean      = Stat1Action<Ostap::Math::LehmerMean>       ;    
    
    using WGeometricMean  = Stat2Action<Ostap::Math::WGeometricMean>  ;    
    using WArithmeticMean = Stat2Action<Ostap::Math::WArithmeticMean> ;    
    using WHarmonicMean   = Stat2Action<Ostap::Math::WHarmonicMean>   ;    
    using WPowerMean      = Stat2Action<Ostap::Math::WPowerMean>      ;    
    using WLehmerMean     = Stat2Action<Ostap::Math::WLehmerMean>     ;    
    
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
//                                                                      The END 
// ============================================================================
#endif // OSTAP_DATAFRAMEACTIONS_H
// ============================================================================
