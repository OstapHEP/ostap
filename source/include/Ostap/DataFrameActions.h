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
#include "Ostap/Statistic.h"

#include "Ostap/StatEntity.h"
#include "Ostap/WStatEntity.h"


// #include "Ostap/Polynomials.h"
// #include "Ostap/Bernstein.h"
// #include "Ostap/Bernstein2D.h"
// #include "Ostap/Bernstein3D.h"
// #include "Ostap/Parameterization.h"
// #include "Ostap/Moments.h"
// #include "Ostap/ECDF.h"
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
      /** @class StatAction1
       *  Helper class to get statitsics for the column in DataFrame 
       *  using some COUNTER class
       *  Requirements for COUNTER:
       *  -  counter.update ( value ) 
       *  -  couner .reset  () 
       *  -  counter += counter 
       *  -  copy constructor 
       */
      template <class COUNTER> 
      class StatAction1 : public ROOT::Detail::RDF::RActionImpl< StatAction1<COUNTER> > 
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
        StatAction1  ( const Args& ...args ) 
          : m_result ( std::make_shared<Result_t>( args ... ) ) 
          , m_N      ( std::max ( 1u , Ostap::Utils::mt_pool_size () ) )
	  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
        {}
        /// constructor 
        StatAction1  ( const COUNTER& cnt ) 
          : m_result ( std::make_shared<Result_t> ( cnt ) ) 
          , m_N      ( std::max ( 1u , Ostap::Utils::mt_pool_size () ) )
	  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
        {}
        /// Move constructor 
        StatAction1 (       StatAction1&& ) = default ;
        /// Copy constructor is disabled 
        StatAction1 ( const StatAction1&  ) = delete ;
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
          for ( unsigned int i = 1 ; i < m_N ; ++i ) { m_slots [ 0 ] += m_slots [ i ] ; }
	  *m_result = m_slots [ 0 ] ;	  
        }
        /// who am I ?
        std::string GetActionName() { return "StatAction1" ; }
        // ====================================================================
      public:
        // ====================================================================
        /// The basic method: increment the counter 
        inline void Exec
	( unsigned int slot ,
	  const double value ) 
        { m_slots [ slot % m_N ].update ( value ) ; } 
        // ====================================================================
      public: // 1 vector column 
	// ====================================================================
        /// The basic method: increment the counter for the vector-like columns       
        template <typename T,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value>::type>
        inline void Exec
	( unsigned int slot ,
	  const T&     vs   )
        {
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v : vs ) { m.update ( v ) ; }
	}
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
      } ; //                     The end of class ROOT::Detail::RDF::StatAction1
      // ======================================================================      
      /** @class StatAction2
       *  Helper class to get statitsics for the column in DataFrame 
       *  using soem COUNTER class
       *  Requirements for COUNTER:
       *  -  counter.update ( value1 , value2 ) 
       *  -  counter.reset () 
       *  -  counter += counter 
       *  -  copy contructor 
       */
      template <class COUNTER>
      class StatAction2 : public ROOT::Detail::RDF::RActionImpl< StatAction2< COUNTER > > 
      {
      public:
        // ====================================================================
        /// define the resutl type 
        using Result_t = COUNTER ;
        // ====================================================================
      public:
        // ====================================================================
        /// default constructor 
	template <typename ...Args>
        StatAction2 ( const Args& ...args ) 
	  : m_result ( std::make_shared<Result_t> ( args ... ) ) 
          , m_N      ( std::max ( 1u , Ostap::Utils::mt_pool_size () ) )
	  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
        {}
	// ====================================================================
        /// constructor 
        StatAction2 ( const COUNTER& cnt  ) 
	  : m_result ( std::make_shared<Result_t> ( cnt ) ) 
          , m_N      ( std::max ( 1u , Ostap::Utils::mt_pool_size () ) )
	  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
        {}	
	/// Move constructor 
        StatAction2      (       StatAction2&& ) = default ;
        /// Copy constructor is disabled 
        StatAction2      ( const StatAction2&  ) = delete ;
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
          for ( unsigned int i = 1 ; i < m_N ; ++i ) { m_slots [ 0 ] += m_slots [ i ] ; }
          *m_result = m_slots [ 0 ] ;
        }
        /// who am I ?
        std::string GetActionName() { return "StatAction2" ; }
        // ====================================================================
      public:
        // ====================================================================
        /// The basic method: increment the counter 
        inline void Exec
	( unsigned int slot   ,
	  const double value1 ,
	  const double value2 ) 
        { m_slots [ slot % m_N ].update ( value1 , value2 ) ; } 
        // ====================================================================
      public: // 1 vector column 
	// ====================================================================
        /// The basic method: increment the counter for the vector-like columns       
        template <typename T,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value>::type>
	inline void Exec
	( unsigned int slot   ,
	  const T&     vs     ,
	  const double value2 )
        {
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v : vs ) { m.update ( v , value2 ) ; }
	}
        /// The basic method: increment the counter for the vector-like columns       
        template <typename T,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value>::type>
	inline void Exec
	( unsigned int slot   ,
	  const double value1 , 
	  const T&     vs     ) 
        {
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v : vs ) { m.update ( value1 , v ) ; }
	}
        // ====================================================================
      public: // 2 vector columns 
	// ====================================================================
        /// The basic method: increment the counter for the vector-like columns       
        template <typename T1,
		  typename T2,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,  
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type> 
	inline void Exec
	( unsigned int slot   , 
	  const T1&    vs1    ,
	  const T2&    vs2    )	 
        {
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto &v1 : vs1 )
	    { for ( const auto &v2 : vs2 )
		{ m.update ( v1 , v2 ) ; } }
	} ; 
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
      } ; //                    The end of class ROOT::Detail::RDF::StatAction2
      // ======================================================================      
      /** @class StatAction3
       *  Helper class to get statitsics for the column in DataFrame 
       *  using some COUNTER class
       *  Requirements for COUNTER:
       *  - counter.update ( v1 , v2 , v3 )
       *  - counter.reset ()  
       *  - counter += counter 
       *  - copy constructor
       */
      template <class COUNTER>
      class StatAction3 : public ROOT::Detail::RDF::RActionImpl< StatAction3<COUNTER> > 
      {
      public:
	// ====================================================================
	/// define the resutl type 
	using Result_t = COUNTER ;
	// ====================================================================
      public:
	// ====================================================================
	/// default constructor 
	template <typename ...Args>
	StatAction3  ( const Args& ...args ) 
	  : m_result ( std::make_shared<Result_t> ( args ... ) ) 
	  , m_N      ( std::max ( 1u , Ostap::Utils::mt_pool_size () ) )
	  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
	{}
	/// default constructor 
	StatAction3  ( const COUNTER& cnt ) 
	  : m_result ( std::make_shared<Result_t> ( cnt      ) ) 
	  , m_N      ( std::max ( 1u , Ostap::Utils::mt_pool_size () ) )
	  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
	{}
	/// Move constructor 
	StatAction3 (       StatAction3&& ) = default ;
	/// Copy constructor is disabled 
	StatAction3 ( const StatAction3&  ) = delete ;
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
	  for ( unsigned int i = 1 ; i < m_N ; ++i ) { m_slots [ 0 ] += m_slots [ i ] ; }
	  *m_result = m_slots [ 0 ] ;
	}
	/// who am I ?
	std::string GetActionName() { return "StatAction3" ; }
	// ====================================================================
      public:
	// ====================================================================
	/// The basic method: increment the counter 
	inline void Exec
	( unsigned int slot   ,
	  const double v1     ,
	  const double v2     ,
	  const double v3     ) 
	{ m_slots [ slot % m_N ].update ( v1 , v2 , v3 ) ; } 
	// ====================================================================
      public: // 1 vector columns 
	// ====================================================================
	/// The basic method: increment the counter for the vector-like columns       
	template <typename T,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value>::type> 
	inline void Exec
	( unsigned int slot ,
	  const T&     vs   ,
	  const double v2   ,
	  const double v3   )
	{
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v1 : vs ) { m.update ( v1 , v2  , v3 ) ; }
	}
	// ====================================================================
	/// The basic method: increment the counter for the vector-like columns
	template <typename T,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value>::type> 
	inline void Exec
	( unsigned int slot ,
	  const double v1   ,
	  const T&     vs   ,
	  const double v3   )
	{
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v2 : vs ) { m.update ( v1 , v2 , v3 ) ; }
	}
	// ====================================================================
	/// The basic method: increment the counter for the vector-like columns
	template <typename T,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value>::type> 
	inline void Exec
	( unsigned int slot ,
	  const double v1   ,
	  const double v2   , 
	  const T&     vs   )
	{
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v3 : vs ) { m.update ( v1 , v2 , v3 ) ; }
	}
	// ====================================================================
      public: // 2 vector columns 
	// ====================================================================
	/// The basic method: increment the counter for the vector-like column of weight 
	template <typename T1,
		  typename T2,	
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type, 
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type> 
	inline void Exec
	( unsigned int slot ,
	  const T1&    vs1  ,
	  const T2&    vs2  ,
	  const double v3   )		    
	{
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v1 : vs1 )
	    { for ( const auto v2 : vs2 )
		{ m.update ( v1 , v2 , v3 ) ; } }
	} ;
	// ====================================================================
	/// The basic method: increment the counter for the vector-like column of weight 
	template <typename T1,
		  typename T2,		  
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type, 
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type>
	inline void Exec
	( unsigned int slot ,
	  const T1&    vs1  ,
	  const double v2   , 	    
	  const T2&    vs3  )
	{
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v1 : vs1 )
	    { for ( const auto v3 : vs3 )
		{ m.update ( v1 , v2 , v3 ) ; } }
	} ;
	// ====================================================================
	/// The basic method: increment the counter for the vector-like column of weight 
	template <typename T1,
		  typename T2,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type, 
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type>
	inline void Exec
	( unsigned int slot ,
	  const double v1   ,  		    
	  const T1&    vs2  ,
	  const T2&    vs3  )
	{
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v2 : vs2 )
	    { for ( const auto v3 : vs3 )
		{ m.update ( v1 , v2 , v3 ) ; } }
	} ;
	// ====================================================================
      public: // 3 vector columns 
	// ====================================================================
	/// The basic method: increment the counter for the vector-like column of weight 
	template <typename T1,
		  typename T2,	
		  typename T3,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type, 
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type, 
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T3>::value>::type>
	inline void Exec
	( unsigned int slot ,
	  const T1&    vs1  ,
	  const T2&    vs2  ,
	  const T3&    vs3  )
	{
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v1 : vs1 )
	    { for ( const auto v2 : vs2 )
		{ for ( const auto v3 : vs3 )
		    { m.update ( v1 , v2 , v3 ) ; } } }
	} ;
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
        // =====================================================================
      } ; //                     The end of class ROOT::Detail::RDF::StatAction3
      // =======================================================================
      /** @class StatAction4
       *  Helper class to get statitsics for the column in DataFrame 
       *  using some COUNTER class
       *  Requirements for COUNTER:
       *  - counter.update ( v1 , v2 , v3 , v4 )
       *  - counter.reset ()  
       *  - counter += counter 
       *  - copy constructor
       */
      template <class COUNTER>
      class StatAction4 : public ROOT::Detail::RDF::RActionImpl< StatAction4<COUNTER> > 
      {
      public:
	// ====================================================================
	/// define the resutl type 
	using Result_t = COUNTER ;
	// ====================================================================
      public:
	// ====================================================================
	/// default constructor 
	template <typename ...Args>
	StatAction4 ( const Args& ...args ) 
	  : m_result ( std::make_shared<Result_t> ( args ... ) ) 
	  , m_N      ( std::max ( 1u , Ostap::Utils::mt_pool_size () ) )
	  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
	{}
	/// default constructor 
	StatAction4 ( const COUNTER& cnt ) 
	  : m_result ( std::make_shared<Result_t> ( cnt ) ) 
	  , m_N      ( std::max ( 1u , Ostap::Utils::mt_pool_size () ) )
	  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
	{}
	/// Move constructor 
	StatAction4 (       StatAction4&& ) = default ;
	/// Copy constructor is disabled 
	StatAction4 ( const StatAction4&  ) = delete ;
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
	  for ( unsigned int i = 1 ; i < m_N ; ++i ) { m_slots [ 0 ] += m_slots [ i ] ; }
	  *m_result = m_slots [ 0 ] ;
	}
	/// who am I ?
	std::string GetActionName() { return "StatAction4" ; }
	// ====================================================================
      public:
	// ====================================================================
	/// The basic method: increment the counter 
	inline void Exec
	( unsigned int slot   ,
	  const double v1     ,
	  const double v2     ,
	  const double v3     , 
	  const double v4     ) 
	{ m_slots [ slot % m_N ].update ( v1 , v2 , v3 , v4 ) ; } 
	// ====================================================================
      public: // 1 vector column 
	// ====================================================================
	/// The basic method: increment the counter for the vector-like columns       
	template <typename T,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value>::type>
	inline void Exec
	( unsigned int slot ,
	  const T&     vs   ,
	  const double v2   ,
	  const double v3   , 
	  const double v4   )
	{
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v1 : vs ) { m.update ( v1 , v2 , v3 , v4 ) ; }
	}
	// ====================================================================
	/// The basic method: increment the counter for the vector-like columns
	template <typename T,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value>::type>
	inline void Exec
	( unsigned int slot ,
	  const double v1   ,
	  const T&     vs   ,
	  const double v3   , 
	  const double v4   )
	{
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v2 : vs ) { m.update ( v1 , v2 , v3 , v4 ) ; }
	}
	// ====================================================================
	/// The basic method: increment the counter for the vector-like columns
	template <typename T,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value>::type>
	void Exec
	( unsigned int slot ,
	  const double v1   ,
	  const double v2   , 
	  const T&     vs   , 
	  const double v4   )
	{
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v3 : vs ) { m.update ( v1 , v2 , v3 , v4 ) ; }
	}
	// ====================================================================
	/// The basic method: increment the counter for the vector-like columns
	template <typename T,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value>::type>
	void Exec
	( unsigned int slot ,
	  const double v1   ,
	  const double v2   , 
	  const double v3   , 
	  const T&     vs   )
	{
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v4 : vs ) { m.update ( v1 , v2 , v3 , v4 ) ; }
	}
	// ====================================================================
      public: // 2 vector column 
	// ====================================================================
	/// The basic method: increment the counter for the vector-like column of weight 
	template <typename T1,
		  typename T2,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type>
	void Exec
	( unsigned int slot ,
	  const T1&    vs1  ,
	  const T2&    vs2  ,
	  const double v3   , 
	  const double v4   )		    
	{
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v1 : vs1 )
	    { for ( const auto v2 : vs2 )
		{ m.update ( v1 , v2 , v3 , v4 ) ; } }
	} ;
	// ====================================================================
	/// The basic method: increment the counter for the vector-like column of weight 
	template <typename T1,
		  typename T2,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type>
	void Exec
	( unsigned int slot ,
	  const T1&    vs1  ,
	  const double v2   ,		    
	  const T2&    vs3  ,
	  const double v4   )
	{
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v1 : vs1 )
	    { for ( const auto v3 : vs3 )
		{ m.update ( v1 , v2 , v3 , v4 ) ; } }
	} ;
	// ====================================================================
	/// The basic method: increment the counter for the vector-like column of weight 
	template <typename T1,
		  typename T2,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type>
	void Exec
	( unsigned int slot ,
	  const T1&    vs1  ,
	  const double v2   ,		    
	  const double v3   ,		    
	  const T2&    vs4  )
	{
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v1 : vs1 )
	    { for ( const auto v4 : vs4 )
		{ m.update ( v1 , v2 , v3 , v4 ) ; } }
	} ;
	// ====================================================================
	/// The basic method: increment the counter for the vector-like column of weight 
	template <typename T1,
		  typename T2,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type>
	void Exec
	( unsigned int slot ,
	  const double v1   ,		    
	  const T1&    vs2  ,
	  const T2&    vs3  , 
	  const double v4   )
	{
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v2 : vs2 )
	    { for ( const auto v3 : vs3 )
		{ m.update ( v1 , v2 , v3 , v4 ) ; } }
	} ;
	// ====================================================================
	/// The basic method: increment the counter for the vector-like column of weight 
	template <typename T1,
		  typename T2,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type>
	void Exec
	( unsigned int slot ,
	  const double v1   ,		    
	  const T1&    vs2  ,
	  const double v3   , 
	  const T2&    vs4  )
	{
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v2 : vs2 )
	    { for ( const auto v4 : vs4 )
		{ m.update ( v1 , v2 , v3 , v4 ) ; } }
	} ;
	// ====================================================================
	/// The basic method: increment the counter for the vector-like column of weight 
	template <typename T1,
		  typename T2,	
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type>
	void Exec
	( unsigned int slot ,
	  const double v1   ,		    
	  const double v2   , 
	  const T1&    vs3  ,
	  const T2&    vs4  )
	{
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v3 : vs3 )
	    { for ( const auto v4 : vs4 )
		{ m.update ( v1 , v2 , v3 , v4 ) ; } }
	} ;
	// ====================================================================
      public: // 3 vector columns 
	// ====================================================================
	template <typename T1,
		  typename T2,	
		  typename T3,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T3>::value>::type>
	void Exec
	( unsigned int slot ,
	  const double v1   ,		    
	  const T1&    vs2  ,
	  const T2&    vs3  ,
	  const T3&    vs4  )
	{
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v2 : vs2 )
	    { for ( const auto v3 : vs3 )
		{ for ( const auto v4 : vs4 )
		    { m.update ( v1 , v2 , v3 , v4 ) ; } } }
	} ;	
	// ====================================================================
	template <typename T1,
		  typename T2,	
		  typename T3,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T3>::value>::type>
	void Exec
	( unsigned int slot ,
	  const T1&    vs1  ,
	  const double v2   ,		    
	  const T2&    vs3  ,
	  const T3&    vs4  )
	{
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v1 : vs1 )
	    { for ( const auto v3 : vs3 )
		{ for ( const auto v4 : vs4 )
		    { m.update ( v1 , v2 , v3 , v4 ) ; } } }
	} ;
	// ====================================================================
	template <typename T1,
		  typename T2,	
		  typename T3,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T3>::value>::type>
	void Exec
	( unsigned int slot ,
	  const T1&    vs1  ,
	  const T2&    vs2  ,
	  const double v3   ,		    
	  const T3&    vs4  )
	{
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v1 : vs1 )
	    { for ( const auto v2 : vs2 )
		{ for ( const auto v4 : vs4 )
		    { m.update ( v1 , v2 , v3 , v4 ) ; } } }
	} ;
	// ====================================================================
	template <typename T1,
		  typename T2,	
		  typename T3,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T3>::value>::type>
	void Exec
	( unsigned int slot ,
	  const T1&    vs1  ,
	  const T2&    vs2  ,
	  const T3&    vs3  , 
	  const double v4   )		    
	{
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v1 : vs1 )
	    { for ( const auto v2 : vs2 )
		{ for ( const auto v3 : vs3 )
		    { m.update ( v1 , v2 , v3 , v4 ) ; } } }
	} ;
	// ====================================================================
      public: // 4 vector columns 
	// ====================================================================
	template <typename T1,
		  typename T2,	
		  typename T3,	
		  typename T4,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T3>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T4>::value>::type>
	void Exec
	( unsigned int slot ,
	  const T1&    vs1  ,
	  const T2&    vs2  ,
	  const T3&    vs3  , 
	  const T4&    vs4  ) 
	{
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v1 : vs1 )
	    { for ( const auto v2 : vs2 )
		{ for ( const auto v3 : vs3 )
		    { for ( const auto v4 : vs4 )
			{ m.update ( v1 , v2 , v3 , v4 ) ; } } } } 
	} ;
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
        // =====================================================================
      } ; //                     The end of class ROOT::Detail::RDF::StatAction4
      // =======================================================================
      
      // =======================================================================
      // Actions with weight 
      // =======================================================================

      // ======================================================================
      /** @class StatAction1w
       *  Helper class to get statitsics for the column in DataFrame 
       *  using some COUNTER class
       *  Requirements for COUNTER:
       *  -  counter.update ( value , weight ) 
       *  -  couner .reset  () 
       *  -  counter += counter 
       *  -  copy constructor 
       */
      template <class COUNTER> 
      class StatAction1w : public ROOT::Detail::RDF::RActionImpl< StatAction1w<COUNTER> > 
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
        StatAction1w  ( const Args& ...args ) 
          : m_result  ( std::make_shared<Result_t>( args ... ) ) 
          , m_N       ( std::max ( 1u , Ostap::Utils::mt_pool_size () ) )
	  , m_slots   ( this->m_N , *(this->m_result.get() ) ) 
        {}
        /// constructor 
        StatAction1w  ( const COUNTER& cnt ) 
          : m_result  ( std::make_shared<Result_t> ( cnt ) ) 
          , m_N       ( std::max ( 1u , Ostap::Utils::mt_pool_size () ) )
	  , m_slots   ( this->m_N , *(this->m_result.get() ) ) 
        {}
        /// Move constructor 
        StatAction1w (       StatAction1w&& ) = default ;
        /// Copy constructor is disabled 
        StatAction1w ( const StatAction1w&  ) = delete ;
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
          for ( unsigned int i = 1 ; i < m_N ; ++i ) { m_slots [ 0 ] += m_slots [ i ] ; }
	  *m_result = m_slots [ 0 ] ;	  
        }
        /// who am I ?
        std::string GetActionName() { return "StatAction1w" ; }
        // ====================================================================
      public:
        // ====================================================================
        /// The basic method: increment the counter 
        inline void Exec
	( unsigned int slot       ,
	  const double value      ,
	  const double weight = 1 ) 
        {
	  if ( !weight ) { return ; } 
	  m_slots [ slot % m_N ].update ( value , weight ) ;
	} 
        // ====================================================================
      public: // 1 vector column 
	// ====================================================================
        /// The basic method: increment the counter for the vector-like columns       
        template <typename T,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value>::type>
        inline void Exec
	( unsigned int slot       ,
	  const T&     vs         ,
	  const double weight = 1 ) 
        {
	  if ( !weight ) { return ; } 
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v : vs ) { m.update ( v , weight ) ; }
	}
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
      } ; //                   The end of class ROOT::Detail::RDF::StatAction1w
      // ======================================================================      
      /** @class StatAction2w
       *  Helper class to get statitsics for the column in DataFrame 
       *  using soem COUNTER class
       *  Requirements for COUNTER:
       *  -  counter.update ( value1 , value2 , weight = 1 ) 
       *  -  counter.reset () 
       *  -  counter += counter 
       *  -  copy contructor 
       */
      template <class COUNTER>
      class StatAction2w : public ROOT::Detail::RDF::RActionImpl< StatAction2w< COUNTER > > 
      {
      public:
        // ====================================================================
        /// define the resutl type 
        using Result_t = COUNTER ;
        // ====================================================================
      public:
        // ====================================================================
        /// default constructor 
	template <typename ...Args>
        StatAction2w ( const Args& ...args ) 
	  : m_result ( std::make_shared<Result_t> ( args ... ) ) 
          , m_N      ( std::max ( 1u , Ostap::Utils::mt_pool_size () ) )
	  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
        {}
	// ====================================================================
        /// constructor 
        StatAction2w ( const COUNTER& cnt  ) 
	  : m_result ( std::make_shared<Result_t> ( cnt ) ) 
          , m_N      ( std::max ( 1u , Ostap::Utils::mt_pool_size () ) )
	  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
        {}	
	/// Move constructor 
        StatAction2w      (       StatAction2w&& ) = default ;
        /// Copy constructor is disabled 
        StatAction2w      ( const StatAction2w&  ) = delete ;
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
          for ( unsigned int i = 1 ; i < m_N ; ++i ) { m_slots [ 0 ] += m_slots [ i ] ; }
          *m_result = m_slots [ 0 ] ;
        }
        /// who am I ?
        std::string GetActionName() { return "StatAction2w" ; }
        // ====================================================================
      public:
        // ====================================================================
        /// The basic method: increment the counter 
        inline void Exec
	( unsigned int slot       ,
	  const double value1     ,
	  const double value2     , 
	  const double weight = 1 ) 
        {
	  if ( !weight ) { return ; } 
	  m_slots [ slot % m_N ].update ( value1 , value2 , weight ) ;
	} 
        // ====================================================================
      public: // 1 vector column 
	// ====================================================================
        /// The basic method: increment the counter for the vector-like columns       
        template <typename T,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value>::type>
	inline void Exec
	( unsigned int slot       ,
	  const T&     vs         ,
	  const double value2     , 
	  const double weight = 1 ) 
        {
	  if ( !weight ) { return ; } 
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v : vs ) { m.update ( v , value2 , weight ) ; }
	}
        /// The basic method: increment the counter for the vector-like columns       
        template <typename T,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value>::type>
	inline void Exec
	( unsigned int slot       ,
	  const double value1     , 
	  const T&     vs         , 
	  const double weight = 1 ) 
        {
	  if ( !weight ) { return ; } 
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v : vs ) { m.update ( value1 , v , weight ) ; }
	}
        // ====================================================================
      public: // 2 vector columns 
	// ====================================================================
        /// The basic method: increment the counter for the vector-like columns       
        template <typename T1,
		  typename T2,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type>
	inline void Exec
	( unsigned int slot       , 
	  const T1&    vs1        ,
	  const T2&    vs2        ,	 
	  const double weight = 1 ) 
        {
	  if ( !weight ) { return ; } 
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto &v1 : vs1 )
	    { for ( const auto &v2 : vs2 )
		{ m.update ( v1 , v2 , weight ) ; } }
	} ; 
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
      } ; //                    The end of class ROOT::Detail::RDF::StatAction2
      // ======================================================================                 
      /** @class StatAction3w
       *  Helper class to get statitsics for the column in DataFrame 
       *  using some COUNTER class
       *  Requirements for COUNTER:
       *  - counter.update ( v1 , v2 , v3 , weight = 1 )
       *  - counter.reset ()  
       *  - counter += counter 
       *  - copy constructor
       */
      template <class COUNTER>
      class StatAction3w : public ROOT::Detail::RDF::RActionImpl< StatAction3w<COUNTER> > 
      {
      public:
	// ====================================================================
	/// define the resutl type 
	using Result_t = COUNTER ;
	// ====================================================================
      public:
	// ====================================================================
	/// default constructor 
	template <typename ...Args>
	StatAction3w ( const Args& ...args ) 
	  : m_result ( std::make_shared<Result_t> ( args ... ) ) 
	  , m_N      ( std::max ( 1u , Ostap::Utils::mt_pool_size () ) )
	  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
	{}
	/// default constructor 
	StatAction3w ( const COUNTER& cnt ) 
	  : m_result ( std::make_shared<Result_t> ( cnt      ) ) 
	  , m_N      ( std::max ( 1u , Ostap::Utils::mt_pool_size () ) )
	  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
	{}
	/// Move constructor 
	StatAction3w (       StatAction3w&& ) = default ;
	/// Copy constructor is disabled 
	StatAction3w ( const StatAction3w&  ) = delete ;
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
	  for ( unsigned int i = 1 ; i < m_N ; ++i ) { m_slots [ 0 ] += m_slots [ i ] ; }
	  *m_result = m_slots [ 0 ] ;
	}
	/// who am I ?
	std::string GetActionName() { return "StatAction3w" ; }
	// ====================================================================
      public:
	// ====================================================================
	/// The basic method: increment the counter 
	inline void Exec
	( unsigned int slot       ,
	  const double v1         ,
	  const double v2         ,
	  const double v3         , 
	  const double weight = 1 ) 
	{
	  if ( !weight ) { return ; } 
	  m_slots [ slot % m_N ].update ( v1 , v2 , v3 , weight ) ;
	} 
	// ====================================================================
      public: // 1 vector columns 
	// ====================================================================
	/// The basic method: increment the counter for the vector-like columns       
	template <typename T,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value>::type>
	inline void Exec
	( unsigned int slot       ,
	  const T&     vs         ,
	  const double v2         ,
	  const double v3         , 
	  const double weight = 1 ) 
	{
	  if ( !weight ) { return ; } 
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v1 : vs ) { m.update ( v1 , v2  , v3 , weight ) ; }
	}
	// ====================================================================
	/// The basic method: increment the counter for the vector-like columns
	template <typename T,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value>::type>
	inline void Exec
	( unsigned int slot       ,
	  const double v1         ,
	  const T&     vs         ,
	  const double v3         , 
	  const double weight = 1 ) 
	{
	  if ( !weight ) { return ; } 
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v2 : vs ) { m.update ( v1 , v2 , v3 , weight ) ; }
	}
	// ====================================================================
	/// The basic method: increment the counter for the vector-like columns
	template <typename T,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value>::type>
	inline void Exec
	( unsigned int slot       ,
	  const double v1         ,
	  const double v2         , 
	  const T&     vs         , 
	  const double weight = 1 ) 
	{
	  if ( !weight ) { return ; } 
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v3 : vs ) { m.update ( v1 , v2 , v3 , weight ) ; }
	}
	// ====================================================================
      public: // 2 vector columns 
	// ====================================================================
	/// The basic method: increment the counter for the vector-like column of weight 
	template <typename T1,
		  typename T2,	
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type>
	inline void Exec
	( unsigned int slot       ,
	  const T1&    vs1        ,
	  const T2&    vs2        ,
	  const double v3         ,  
	  const double weight = 1 ) 
	{
	  if ( !weight ) { return ; } 
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v1 : vs1 )
	    { for ( const auto v2 : vs2 )
		{ m.update ( v1 , v2 , v3 , weight ) ; } }
	} ;
	// ====================================================================
	/// The basic method: increment the counter for the vector-like column of weight 
	template <typename T1,
		  typename T2,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type>
	inline void Exec
	( unsigned int slot       ,
	  const T1&    vs1        ,
	  const double v2         , 		    
	  const T2&    vs3        ,
	  const double weight = 1 ) 
	{
	  if ( !weight ) { return ; } 
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v1 : vs1 )
	    { for ( const auto v3 : vs3 )
		{ m.update ( v1 , v2 , v3 , weight ) ; } }
	} ;
	// ====================================================================
	/// The basic method: increment the counter for the vector-like column of weight 
	template <typename T1,
		  typename T2,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type>
	inline void Exec
	( unsigned int slot       ,
	  const double v1         ,  		    
	  const T1&    vs2        ,
	  const T2&    vs3        , 
	  const double weight = 1 ) 
	{
	  if ( !weight ) { return ; } 
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v2 : vs2 )
	    { for ( const auto v3 : vs3 )
		{ m.update ( v1 , v2 , v3 , weight ) ; } }
	} ;
	// ====================================================================
      public: // 3 vector columns 
	// ====================================================================
	/// The basic method: increment the counter for the vector-like column of weight 
	template <typename T1,
		  typename T2,	
		  typename T3,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T3>::value>::type>
	inline void Exec
	( unsigned int slot       ,
	  const T1&    vs1        ,
	  const T2&    vs2        ,
	  const T3&    vs3        , 
	  const double weight = 1 ) 
	{
	  if ( !weight ) { return ; } 
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v1 : vs1 )
	    { for ( const auto v2 : vs2 )
		{ for ( const auto v3 : vs3 )
		    { m.update ( v1 , v2 , v3 ) ; } } }
	} ;
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
        // =====================================================================
      } ; //                    The end of class ROOT::Detail::RDF::StatAction3w
      // =======================================================================
      /** @class StatAction4w
       *  Helper class to get statitsics for the column in DataFrame 
       *  using some COUNTER class
       *  Requirements for COUNTER:
       *  - counter.update ( v1 , v2 , v3 , v4 , weight = 1 )
       *  - counter.reset ()  
       *  - counter += counter 
       *  - copy constructor
       */
      template <class COUNTER>
      class StatAction4w : public ROOT::Detail::RDF::RActionImpl< StatAction4w<COUNTER> > 
      {
      public:
	// ====================================================================
	/// define the resutl type 
	using Result_t = COUNTER ;
	// ====================================================================
      public:
	// ====================================================================
	/// default constructor 
	template <typename ...Args>
	StatAction4w ( const Args& ...args ) 
	  : m_result ( std::make_shared<Result_t> ( args ... ) ) 
	  , m_N      ( std::max ( 1u , Ostap::Utils::mt_pool_size () ) )
	  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
	{}
	/// default constructor 
	StatAction4w ( const COUNTER& cnt ) 
	  : m_result ( std::make_shared<Result_t> ( cnt ) ) 
	  , m_N      ( std::max ( 1u , Ostap::Utils::mt_pool_size () ) )
	  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
	{}
	/// Move constructor 
	StatAction4w (       StatAction4w&& ) = default ;
	/// Copy constructor is disabled 
	StatAction4w ( const StatAction4w&  ) = delete ;
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
	  for ( unsigned int i = 1 ; i < m_N ; ++i ) { m_slots [ 0 ] += m_slots [ i ] ; }
	  *m_result = m_slots [ 0 ] ;
	}
	/// who am I ?
	std::string GetActionName() { return "StatAction4w" ; }
	// ====================================================================
      public:
	// ====================================================================
	/// The basic method: increment the counter 
	inline void Exec
	( unsigned int slot       ,
	  const double v1         ,
	  const double v2         ,
	  const double v3         , 
	  const double v4         , 
	  const double weight = 1 )  
	{
	  if ( !weight ) { return ; } 
	  m_slots [ slot % m_N ].update ( v1 , v2 , v3 , v4 , weight ) ;
	} 
	// ====================================================================
      public: // 1 vector column 
	// ====================================================================
	/// The basic method: increment the counter for the vector-like columns       
	template <typename T,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value>::type>
	inline void Exec
	( unsigned int slot       ,
	  const T&     vs         ,
	  const double v2         ,
	  const double v3         , 
	  const double v4         , 
	  const double weight = 1 )  
	{
	  if ( !weight ) { return ; } 
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v1 : vs ) { m.update ( v1 , v2 , v3 , v4 , weight ) ; }
	}
	// ====================================================================
	/// The basic method: increment the counter for the vector-like columns
	template <typename T,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value>::type>
	inline void Exec
	( unsigned int slot       ,
	  const double v1         ,
	  const T&     vs         ,
	  const double v3         , 
	  const double v4         , 
	  const double weight = 1 )  
	{
	  if ( !weight ) { return ; } 
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v2 : vs ) { m.update ( v1 , v2 , v3 , v4 , weight ) ; }
	}
	// ====================================================================
	/// The basic method: increment the counter for the vector-like columns
	template <typename T,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value>::type>
	void Exec
	( unsigned int slot       ,
	  const double v1         ,
	  const double v2         , 
	  const T&     vs         , 
	  const double v4         , 
	  const double weight = 1 )  
	{
	  if ( !weight ) { return ; } 
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v3 : vs ) { m.update ( v1 , v2 , v3 , v4 , weight ) ; }
	}
	// ====================================================================
	/// The basic method: increment the counter for the vector-like columns
	template <typename T,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T>::value>::type>
	void Exec
	( unsigned int slot       ,
	  const double v1         ,
	  const double v2         , 
	  const double v3         , 
	  const T&     vs         , 
	  const double weight = 1 )  
	{
	  if ( !weight ) { return ; } 
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v4 : vs ) { m.update ( v1 , v2 , v3 , v4 , weight ) ; }
	}
	// ====================================================================
      public: // 2 vector column 
	// ====================================================================
	/// The basic method: increment the counter for the vector-like column of weight 
	template <typename T1,
		  typename T2,	
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type>
	void Exec
	( unsigned int slot       ,
	  const T1&    vs1        ,
	  const T2&    vs2        ,
	  const double v3         , 
	  const double v4         , 		    
	  const double weight = 1 )  
	{
	  if ( !weight ) { return ; } 
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v1 : vs1 )
	    { for ( const auto v2 : vs2 )
		{ m.update ( v1 , v2 , v3 , v4 , weight ) ; } }
	} ;
	// ====================================================================
	/// The basic method: increment the counter for the vector-like column of weight 
	template <typename T1,
		  typename T2,	
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type>
	void Exec
	( unsigned int slot       ,
	  const T1&    vs1        ,
	  const double v2         ,		    
	  const T2&    vs3        ,
	  const double v4         ,
	  const double weight = 1 )  
	{
	  if ( !weight ) { return ; } 
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v1 : vs1 )
	    { for ( const auto v3 : vs3 )
		{ m.update ( v1 , v2 , v3 , v4 , weight ) ; } }
	} ;
	// ====================================================================
	/// The basic method: increment the counter for the vector-like column of weight 
	template <typename T1,
		  typename T2,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type>
	void Exec
	( unsigned int slot       ,
	  const T1&    vs1        ,
	  const double v2         ,		    
	  const double v3         ,		    
	  const T2&    vs4        , 
	  const double weight = 1 )  
	{
	  if ( !weight ) { return ; } 
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v1 : vs1 )
	    { for ( const auto v4 : vs4 )
		{ m.update ( v1 , v2 , v3 , v4 , weight ) ; } }
	} ;
	// ====================================================================
	/// The basic method: increment the counter for the vector-like column of weight 
	template <typename T1,
		  typename T2,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type>
	void Exec
	( unsigned int slot       ,
	  const double v1         ,		    
	  const T1&    vs2        ,
	  const T2&    vs3        , 
	  const double v4         ,
	  const double weight = 1 )  
	{
	  if ( !weight ) { return ; } 
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v2 : vs2 )
	    { for ( const auto v3 : vs3 )
		{ m.update ( v1 , v2 , v3 , v4 , weight ) ; } }
	} ;
	// ====================================================================
	/// The basic method: increment the counter for the vector-like column of weight 
	template <typename T1,
		  typename T2,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type>
	void Exec
	( unsigned int slot       ,
	  const double v1         ,		    
	  const T1&    vs2        ,
	  const double v3         ,  
	  const T2&    vs4        , 
	  const double weight = 1 )  
	{
	  if ( !weight ) { return ; } 
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v2 : vs2 )
	    { for ( const auto v4 : vs4 )
		{ m.update ( v1 , v2 , v3 , v4 , weight ) ; } }
	} ;
	// ====================================================================
	/// The basic method: increment the counter for the vector-like column of weight 
	template <typename T1,
		  typename T2,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type>
	void Exec
	( unsigned int slot       ,
	  const double v1         ,		    
	  const double v2         , 
	  const T1&    vs3        ,
	  const T2&    vs4        , 
	  const double weight = 1 )  
	{
	  if ( !weight ) { return ; } 
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v3 : vs3 )
	    { for ( const auto v4 : vs4 )
		{ m.update ( v1 , v2 , v3 , v4 , weight ) ; } }
	} ;
	// ====================================================================
      public: // 3 vector columns 
	// ====================================================================
	template <typename T1,
		  typename T2,	
		  typename T3,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T3>::value>::type>
	void Exec
	( unsigned int slot       ,
	  const double v1         ,		    
	  const T1&    vs2        ,
	  const T2&    vs3        ,
	  const T3&    vs4        ,
	  const double weight = 1 )  
	{
	  if ( !weight ) { return ; } 
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v2 : vs2 )
	    { for ( const auto v3 : vs3 )
		{ for ( const auto v4 : vs4 )
		    { m.update ( v1 , v2 , v3 , v4 , weight ) ; } } }
	} ;	
	// ====================================================================
	template <typename T1,
		  typename T2,	
		  typename T3,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T3>::value>::type>
	void Exec
	( unsigned int slot       ,
	  const T1&    vs1        ,
	  const double v2         ,		    
	  const T2&    vs3        ,
	  const T3&    vs4        , 
	  const double weight = 1 )  
	{
	  if ( !weight ) { return ; } 
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v1 : vs1 )
	    { for ( const auto v3 : vs3 )
		{ for ( const auto v4 : vs4 )
		    { m.update ( v1 , v2 , v3 , v4 , weight ) ; } } }
	} ;
	// ====================================================================
	template <typename T1,
		  typename T2,	
		  typename T3,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T3>::value>::type>
	void Exec
	( unsigned int slot       ,
	  const T1&    vs1        ,
	  const T2&    vs2        ,
	  const double v3         ,		    
	  const T3&    vs4        , 
	  const double weight = 1 )  
	{
	  if ( !weight ) { return ; } 
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v1 : vs1 )
	    { for ( const auto v2 : vs2 )
		{ for ( const auto v4 : vs4 )
		    { m.update ( v1 , v2 , v3 , v4 , weight ) ; } } }
	} ;
	// ====================================================================
	template <typename T1,
		  typename T2,	
		  typename T3,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T3>::value>::type>
	void Exec
	( unsigned int slot       ,
	  const T1&    vs1        ,
	  const T2&    vs2        ,
	  const T3&    vs3        , 
	  const double v4         ,		    
	  const double weight = 1 )  
	{
	  if ( !weight ) { return ; } 
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v1 : vs1 )
	    { for ( const auto v2 : vs2 )
		{ for ( const auto v3 : vs3 )
		    { m.update ( v1 , v2 , v3 , v4 , weight ) ; } } }
	} ;
	// ====================================================================
      public: // 4 vector columns 
	// ====================================================================
	template <typename T1,
		  typename T2,	
		  typename T3,	
		  typename T4,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T1>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T2>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T3>::value>::type,
		  std::enable_if<ROOT::Internal::RDF::IsDataContainer<T4>::value>::type>
	void Exec
	( unsigned int slot       ,
	  const T1&    vs1        ,
	  const T2&    vs2        ,
	  const T3&    vs3        , 
	  const T4&    vs4        , 
	  const double weight = 1 )  
	{
	  if ( !weight ) { return ; } 
	  Result_t& m = m_slots [ slot % m_N ] ;
	  for ( const auto v1 : vs1 )
	    { for ( const auto v2 : vs2 )
		{ for ( const auto v3 : vs3 )
		    { for ( const auto v4 : vs4 )
			{ m.update ( v1 , v2 , v3 , v4 , weight ) ; } } } } 
	} ;
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
      } ; //                    The end of class ROOT::Detail::RDF::StatAction4
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
    /** 
    /// fake helper class to allow namespace loading & lookup. 
    class Action {} ;
    // ========================================================================
    template <class COUNTER>
    using StatAction1  = ROOT::Detail::RDF::StatAction1  <COUNTER> ;
    // ========================================================================
    template <class COUNTER>
    using StatAction1w = ROOT::Detail::RDF::StatAction1w <COUNTER> ;
    // ========================================================================
    template <class COUNTER>
    using StatAction2  = ROOT::Detail::RDF::StatAction2  <COUNTER> ;
    // ========================================================================
    template <class COUNTER>
    using StatAction2w = ROOT::Detail::RDF::StatAction2w <COUNTER> ;
    // ========================================================================
    template <class COUNTER>
    using StatAction3  = ROOT::Detail::RDF::StatAction3  <COUNTER> ;
    // ========================================================================
    template <class COUNTER>
    using StatAction3w = ROOT::Detail::RDF::StatAction3w <COUNTER> ;
    // ========================================================================
    template <class COUNTER>
    using StatAction4  = ROOT::Detail::RDF::StatAction4  <COUNTER> ;
    // ========================================================================
    template <class COUNTER>
    using StatAction4w = ROOT::Detail::RDF::StatAction4w <COUNTER> ;
    */

    /** 
    // ========================================================================
    template <class COUNTER,
	      std::enable_if<std::is_convertible<COUNTER,Ostap::Math::Statistic>::value,bool>::type>      
    inline
    ROOT::Detail::RDF::StatAction1<COUNTER>
    action1
    ( const COUNTER& cnt )
    { return ROOT::Detail::RDF::StatAction1<COUNTER>  ( cnt ) ; }
    
    // ========================================================================
    template <class COUNTER, 
	      std::enable_if<std::is_convertible<COUNTER,Ostap::Math::WStatistic>::value,bool>::type>      
    inline
    ROOT::Detail::RDF::StatAction1w<COUNTER>
    action1w
    ( const COUNTER& cnt )
    { return ROOT::Detail::RDF::StatAction1w<COUNTER> ( cnt ) ; }
    */

    
    /** 
     // ========================================================================
     template <class COUNTER,
     std::enable_if<std::is_convertible<COUNTER,Ostap::Math::Statistic2>::value>::type>      
     inline
     ROOT::Detail::RDF::StatAction2<COUNTER>
     action2
    ( const COUNTER& cnt )
    { return ROOT::Detail::RDF::StatAction2<COUNTER>  ( cnt ) ; }
    // ========================================================================
    template <class COUNTER,
	      std::enable_if<std::is_convertible<COUNTER,Ostap::Math::WStatistic2>::value>::type>      
    inline
    ROOT::Detail::RDF::StatAction2w<COUNTER>
    action2w
    ( const COUNTER& cnt )
    { return ROOT::Detail::RDF::StatAction2w<COUNTER> ( cnt ) ; }
    // ========================================================================
    template <class COUNTER,
	      std::enable_if<std::is_convertible<COUNTER,Ostap::Math::Statistic3>::value>::type>      
    inline
    ROOT::Detail::RDF::StatAction3<COUNTER>
    action3
    ( const COUNTER& cnt )
    { return ROOT::Detail::RDF::StatAction2<COUNTER>  ( cnt ) ; }
    // ========================================================================
    template <class COUNTER,
	      std::enable_if<std::is_convertible<COUNTER,Ostap::Math::WStatistic3>::value>::type>      
    inline
    ROOT::Detail::RDF::StatAction3w<COUNTER>
    action3w
    ( const COUNTER& cnt )
    { return ROOT::Detail::RDF::StatAction3w<COUNTER> ( cnt ) ; }
    // ========================================================================
    template <class COUNTER,
	      std::enable_if<std::is_convertible<COUNTER,Ostap::Math::Statistic4>::value>::type>      
    inline
    ROOT::Detail::RDF::StatAction4<COUNTER>
    action4
    ( const COUNTER& cnt )
    { return ROOT::Detail::RDF::StatAction4<COUNTER>  ( cnt ) ; }
    // ========================================================================
    template <class COUNTER,
    std::enable_if<std::is_convertible<COUNTER,Ostap::Math::WStatistic4>::value>::type>      
    inline
    ROOT::Detail::RDF::StatAction4w<COUNTER>
    action4w
    ( const COUNTER& cnt )
    { return ROOT::Detail::RDF::StatAction4w<COUNTER> ( cnt ) ; }
    // ========================================================================

    */

    
    /** 
    
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

    */
    
    // ========================================================================

  }
  // ==========================================================================
}
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_DATAFRAMEACTIONS_H
// ============================================================================
