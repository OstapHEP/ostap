// ============================================================================
#ifndef OSTAP_PARAMETERS_H 
#define OSTAP_PARAMETERS_H 1
// ============================================================================
// Include files
// ============================================================================
// STD& STL 
// ============================================================================
#include <initializer_list>
#include <vector>
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /** @class Parameters Ostap/Parameneters.h
     *  Holder for parameters 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-02-10
     */
    class Parameters 
    {
    public:
      // =======================================================================
      /// actual type of parameters 
      typedef std::vector<double>        PARAMETERS ;
      /// iterator type 
      typedef PARAMETERS::const_iterator iterator   ;
      // =======================================================================
    public:
      // =======================================================================
      /// constructor from number of parameters 
      Parameters ( const std::size_t           np = 1 ) ; 
      /// constructor from  the list of parameters 
      Parameters ( const PARAMETERS&  pars   ) ;
      /// constructor from  the list of parameters 
      Parameters (       PARAMETERS&& pars   ) ;
      /// templated constructor from the sequence of parameters 
      template <typename ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename std::enable_if<std::is_convertible<value_type,long double>::value,bool>::type = true >           
      Parameters
      ( ITERATOR begin , 
        ITERATOR end   )
        : m_pars ( begin , end )
      {}
      // ======================================================================
    public:
      // ======================================================================
      /// number of parameters 
      inline std::size_t npars  () const { return m_pars.size()     ; }
      // ======================================================================
      /// get the parameter value
      inline double  par          ( const std::size_t k ) const
      { return ( k < m_pars.size() ) ? m_pars [ k ] : 0.0 ; }
      // ======================================================================
      /// get all parameters:
      const PARAMETERS& pars () const { return m_pars ; }
      // ======================================================================
      /** set k-parameter
       *  @param k index
       *  @param value new value 
       *  @return true if parameter is actually changed 
       */
      inline bool setPar          
      ( const std::size_t  k             , 
        const double       value         ,
        const bool         force = false ) 
      { return k < m_pars.size() ?  _setPar ( k , value , force ) : false ; }
      // ======================================================================
      /** set several/all parameters at once 
       *  @param begin  start iterator for the sequence of coefficients 
       *  @param end    end   iterator for the sequence of coefficients 
       *  @return true if at least one parameter is actually changed 
       */
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type ,                
                typename std::enable_if<std::is_convertible<value_type,long double>::value,bool>::type = true >           
      inline bool setPars 
      ( ITERATOR   begin         , 
        ITERATOR   end           ,
        const bool force = false ) 
      {
        bool update = false ;
        const unsigned int   N = m_pars.size()  ;
        for ( unsigned short k ; k < N && begin != end ;  ++k, ++begin ) 
          { update = this->_setPar ( k  , *begin , force ) ? true : update ; }
        return update ;
      }
      // ======================================================================
      /** set several/all parameters at once 
       *  @param pars (NIPUT) vector of parameters 
       *  @return true if at least one parameter is actually changed 
       */
      inline bool setPars
      ( const PARAMETERS& pars          ,
        const bool        force = false ) 
      { return setPars ( pars.begin() , pars.end() , force ) ; }
      // ======================================================================
      /// all parameters are zero ?
      bool           zero   () const ;
      // ======================================================================
    public: // Reset 
      // ======================================================================
      /// reset all parameters to zero 
      void reset() ;
      /// reset all parameters to zero 
      void Reset() { reset () ; }
      // ======================================================================
    public:
      // ======================================================================
      /** Filter out very small terms/parameters 
       *  the term is considered to be very small if
       *   - it is numerically zero: \f$ c_k \approx 0 \f% 
       *   - or if epsilon  > 0: \f$ \left| c_k \right| \le \epsilon \f$ 
       *   - or if scale   != 0: \f$ \left| s \right| + \left| c_k \right| \approx \left| s \right| \f$ 
       *  @param  epsilon  parameter to define "smalness" of terms
       *  @param  scale    parameter to define "smalness" of terms
       *  @return number of nullified terms
       */
      std::size_t remove_noise 
      ( const double epsilon = 0 , 
        const double scale   = 0 ) ;
      // ======================================================================
    public: // simple  manipulations with parameters 
      // ======================================================================
      /// simple  manipulations with parameters: scale it! 
      /// Parameters& operator *= ( const double a ) ;     // scale it! 
      /// simple  manipulations with parameters: scale it! 
      /// Parameters& operator /= ( const double a ) ;     // scale it! 
      // ======================================================================
    public: // expose constant iterators 
      // ======================================================================
      /// begin iterator
      iterator begin () const { return m_pars.begin () ; }
      /// end   iterator
      iterator end   () const { return m_pars.end   () ; }
      // ======================================================================
    protected:
      // ======================================================================
      /// swap two parameter sets 
      void swap ( Parameters& right ) ;
      // ======================================================================
    private:
      // ======================================================================
      /** set k-parameter
       *  @param k index
       *  @param value new value 
       *  @return true if parameter is actually changed 
       */
    bool _setPar 
    ( const std::size_t k             , 
      const double      value         ,
      const bool        force = false ) ;
      // ======================================================================
    protected :
      // ======================================================================
      /// parameters 
      PARAMETERS m_pars ; //  vector of parameters 
      // ======================================================================
    } ;
    // ========================================================================
  } //                                             end of namespace Ostap::Math 
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_PARAMETERS_H
// ============================================================================
