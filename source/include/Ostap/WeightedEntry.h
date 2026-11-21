// ============================================================================
#ifndef OSTAP_WEIGHTEDENTRY_H 
#define OSTAP_WEIGHTEDENTRY_H 1
// ============================================================================
// STD&STL
// ============================================================================
#include <cmath>
#include <utility>
// ============================================================================
namespace Ostap
{
  // ===========================================================================
  namespace  Math
  {
    // =========================================================================
    /** @class WeightedEntry
     */
    class WeightedEntry
    {
    public:
      // =======================================================================
      typedef double                 Data   ;
      typedef double                 Weight ;
      typedef std::pair<Data,Weight> Pair   ; 
      // =======================================================================
    public :
      // =======================================================================
      WeightedEntry
        ( const Data   data   = 0 ,
          const Weight weight = 1 )
        : m_pair ( data , weight )
        {}
      // =======================================================================
      /// constructor from the pair 
      WeightedEntry ( const Pair& data )
        : m_pair ( data )
        {}
      // =======================================================================
    public :
      // =======================================================================
      /// data 
      inline Data          data    () const { return m_pair.first  ; }
      /// data 
      inline Data          value   () const { return m_pair.first  ; }
      /// weight 
      inline Weight        weight  () const { return m_pair.second ; }
      /// weight squared 
      inline Weight        weight2 () const { return m_pair.second * m_pair.second ; }
      /// conversion to the pair
      inline operator const Pair&  () const { return m_pair        ; }
      // =======================================================================      
    public :
      // =======================================================================
      /// both data and weight are finite ?
      inline bool isfinite () const
      { return std::isfinite ( m_pair.first  ) && std::isfinite ( m_pair.second ) ; } 
      /// valid entry: data and weigth are finite and weight is not zero ?
      inline bool valid    () const
      { return this->isfinite () && m_pair.second ; }            
      /// good entry == valid entry 
      inline bool good     () const { return this->valid () ; }            
      // =======================================================================
    public : 
      // =======================================================================
      /// swap two objects 
      inline void swap ( WeightedEntry& entry )
      { std::swap ( this->m_pair , entry.m_pair ) ;  }
      // =======================================================================      
    private :
      // =======================================================================      
      /// actual data is here       
      Pair m_pair { 0 , 1 } ;
      // =======================================================================      
    } ;
    // =========================================================================
    /// swap two objects
    inline void swap ( WeightedEntry& a , WeightedEntry& b ) { a.swap ( b ) ; }    
    // =========================================================================    
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                The  end of namespace Ostap 
// ============================================================================
#endif // OSTAP_WEIGHTEDENTRY_H
// ============================================================================
//                                                                      The END 
// ============================================================================
