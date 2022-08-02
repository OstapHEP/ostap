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
      /// constructor from number of parameters 
      Parameters ( const std::size_t           np = 1 ) ; 
      /// constructor from  the list of parameters 
      Parameters ( const std::vector<double>&  pars   ) ;
      /// constructor from  the list of parameters 
      Parameters (       std::vector<double>&& pars   ) ;
      /// templated constructor from the sequnce of parameters 
      template <typename ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename = std::enable_if<std::is_convertible<value_type,long double>::value> >
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
      { return ( k < m_pars.size() ) ? m_pars[k] : 0.0 ; }
      // ======================================================================
      /// get all parameters:
      const std::vector<double>& pars () const { return m_pars ; }
      // ======================================================================
      /** set k-parameter
       *  @param k index
       *  @param value new value 
       *  @return true if parameter is actually changed 
       */
      inline bool setPar          
      ( const std::size_t  k     , 
        const double       value ) 
      { return k < m_pars.size() ?  _setPar ( k , value ) : false ; }
      // ======================================================================
      /** set several/all parameters at once 
       *  @param begin  start iterator for the sequence of coefficients 
       *  @param end    end   iterator for the sequence of coefficients 
       *  @return true if at least one parameter is actually changed 
       */
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type ,
                typename = std::enable_if<std::is_convertible<value_type,long double>::value> >
      inline bool setPars 
      ( ITERATOR begin  , 
        ITERATOR end    ) ;
      // ======================================================================
      /** set several/all parameters at once 
       *  @param pars (NIPUT) vector of parameters 
       *  @return true if at least one parameter is actually changed 
       */
      inline bool setPars ( const std::vector<double>& pars ) 
      { return setPars ( pars.begin() , pars.end() ) ; }
      // ======================================================================
      /// all parameters are zero ?
      bool           zero   () const ;
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
      typedef std::vector<double>::const_iterator iterator ;
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
      ( const std::size_t k     , 
        const double      value ) ;
      // ======================================================================
    protected :
      // ======================================================================
      /// parameters 
      std::vector<double> m_pars ; //  vector of parameters 
      // ======================================================================
    } ;
    // ========================================================================
    /** set several/all parameters at once 
     *  @param pars (NIPUT) vector of parameters 
     *  @return true if at least one parameter is actually changed 
     */
    template <class ITERATOR,
              typename value_type = typename std::iterator_traits<ITERATOR>::value_type ,
              typename = std::enable_if<std::is_convertible<value_type,long double>::value> >
    inline bool Parameters::setPars 
    ( ITERATOR begin  , 
      ITERATOR end    ) 
    {
      bool update = false ;
      const unsigned int   N = m_pars.size()  ;
      for ( unsigned short k ; k < N && begin != end ;  ++k, ++begin ) 
      { update = _setPar ( k  , *begin ) ? true : update ; }
      return update ;
    }
    // ========================================================================
  } //                                             end of namespace Ostap::Math 
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_POLYNOMIALS_H
// ============================================================================
