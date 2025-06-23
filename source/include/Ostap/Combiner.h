// ============================================================================
#ifndef OSTAP_Combiner_H
#define OSTAP_Combiner_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <algorithm>
#include <vector>
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/Range.h"
// ============================================================================
/** @file
 *
 *  This file is copied from the LoKi project -
 *    "C++ ToolKit  for Smart and Friendly Physics Analysis"
 *
 *  The package has been designed with the kind help from
 *  Galina PAKHLOVA and Sergey BARSUK.  Many bright ideas,
 *  contributions and advices from G.Raven, J.van Tilburg,
 *  A.Golutvin, P.Koppenburg have been used in the design.
 *
 *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
 *  @date 2001-01-23
 */
// ============================================================================
namespace Ostap
 {
  // ==========================================================================
  /** @class Combiner_ Combiner.h LoKi/Combiner.h
   *
   *  Definition of (multy) container with (multy)iterator
   *  Helpful class for implemenattion of variouse looping
   *  techniques
   *
   *  It is in the spirit of former CLHEP Chooser/Combiner_ class,
   *  A. Alexandresku's Loki library and and T. Glebe's (Hera-B)
   *  GCombiner package
   *
   *  The class allows combine "everything" - the content of the
   *  containers is irelevant, it could be objects, pointers, etc.,
   *  even the primitive scalars are allowed.
   *  E.g the examples below shows the loop over all unique
   *  combinations of integers from 2 containers ("double loop").
   *
   *  @code
   *
   *  void  example ()
   *  {
   *
   *  using namespace Ostap;
   *  typedef std::vector<int> Vect;
   *  typedef Combiner_<Vect>  Comp;
   *  typedef Range_<Vect>     Range;
   *
   *  int a1[] = { 1  , 2 , 3 , 4 , 5   };
   *  int a2[] = { 6  , 7 , 8 , 9 , 10  };
   *
   *  Vect v1( a1 , a1 + 5 ) ;
   *  Vect v2( a2 , a2 + 5 ) ;
   *
   *  Comp Combiner;
   *
   *  Range r1( v1.begin() , v1.end() ) ;
   *  Range r2( v2.begin() , v2.end() ) ;
   *
   *  Combiner.add( r1 ) ;
   *  Combiner.add( r2 ) ;
   *
   *  std::cout
   *    << " Combiner parameters "
   *    << " dim  = "       << Combiner.dim()
   *    << " total size = " << Combiner.size() << std::endl ;
   *
   *  int unique = 0 ;
   *  for ( ; Combiner.valid() ; ++Combiner  )
   *  {
   *    if( !Combiner.unique() ) { continue ; }
   *
   *    const Comp::Select& sel = Combiner.current() ;
   *    std::cout
   *      << "  Unique combination #"
   *      << ++unique << " \t selected : \t" ;
   *    for ( Comp::Select::const_iterator it = sel.begin() ;
   *           sel.end() != it ; ++it )
   *    { std::cout << (**it) << " \t " ; }
   *    //
   *    std::cout << std::endl;
   *  }
   *
   *}; // end of the function exmaple
   *
   *  @endcode
   *
   *  @attention The input data are not ownered by
   *  combiner, it operates only with iterators !
   *
   *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
   *  @date   2002-07-11
   */
  template <class CONTAINER>
  class Combiner_ 
  {
  public:
    // ========================================================================
    /// type for the range
    typedef Ostap::Range_<CONTAINER>      Range     ;    
    /// type for actual contained iterator
    typedef typename Range::Container     Container ; 
    /// iterator type
    typedef typename Range::const_iterator iterator ;
    /// definition of (multy)iterator
    typedef std::vector<iterator>          Select   ;
    // ========================================================================
  private:
    // ========================================================================
    /// the Combiner_ itself
    typedef std::vector<Range> Ranges;
    /// iterator over selected combination
    typedef typename Select::iterator it;
    /// iterator over selected combination (const version)
    typedef typename Select::const_iterator cit;
    /// iterator over the Combiner_
    typedef typename Ranges::iterator rit;
    /// iterator over the Combiner_ (cont version)
    typedef typename Ranges::const_iterator rcit;
    // ========================================================================
  public: 
    // ========================================================================
    template <class ITERATOR>
    Combiner_ 
    ( ITERATOR first , 
      ITERATOR last  )
      : m_ranges  () 
      , m_current ()
      , m_size    ( 0 )
      , m_index   ( 0 )  
    { 
      this->add ( first , last ) ;
    }
    // ========================================================================
  public: 
    // ========================================================================
    /// current dimentions (=number of components)  of the Combiner_
    std::size_t N     () const { return m_ranges.size(); }
    // ========================================================================
    /// the size of the Combiner_  (number of ALL combintions)
    std::size_t size () const { return m_size ; }
    // ========================================================================
    template <class ITERATOR>
    Combiner_& add 
    ( ITERATOR first , 
      ITERATOR last  )
    { 
      for  ( ; first != last ; ++first ) { add ( *first ) ; } 
      return *this ;
    }
    // ========================================================================
    /// add one more container 
    Combiner_& add ( const Container& c )
    { return add ( Range ( c.begin () , c.end() ) ) ; } 
    // ========================================================================
    /** add one more "container" to the existing Combiner_.
     *  @param range  the sequence range
     *  return self reference
     */
    Combiner_& add ( const Range& range ) 
    {
      // extend the Combiner_
      m_ranges.push_back  ( range         ) ;
      m_current.push_back ( range.begin() ) ;
      // reset current  (multy) iterator
      return reset() ;
    }
    // ========================================================================
    /** pre-increment operator for the Combiner_ (advance current iterator)
     *  @return self reference
     */
    Combiner_& operator++() 
    {
      /// advance the 'current' iterator
      next();
      return *this;
    }
    // ========================================================================
    /** post-increment operator for the Combiner_ (advance current iterator)
     *  @attention the same as pre-increment
     *  @return self reference
     */
    Combiner_& operator++( int ) { return ++( *this ); }
    // ========================================================================
    /** reset current (multy)iterator
     * go to initial state 
     *  @return self-reference
     */
    Combiner_& reset() 
    {
      // reset current  (multi) iterator
      m_current.resize ( m_ranges.size() ) ;  
      std::transform   ( m_ranges.begin(), m_ranges.end(), m_current.begin(), []( const Range& r ) { return r.begin(); } );
      // reset the current index
      m_index = 0 ;
      m_size  = 1 ;
      for ( const auto&r : m_ranges ) { m_size *= r.size(); }
      // return
      return *this;
    }
    // ========================================================================
    /** access to the current (multy)iterator (const version)
     *  @return 'current' (multy)iterator
     */
    const Select& current() const { return m_current; }
    // ========================================================================
    /** check the validity of current (multi)iterator
     *  @return true if no "actual" sub-iterators are invalid
     */
    bool valid() const { return m_index < m_size ; }
    // ========================================================================
    /** advance 'current' (multi)iterator.
     *  I guess it the most tricky funtion of the whole class.
     *  It is the most primitive one, but I've spent few hours
     *  trying to debug it :-)
     *  @return 'current' (multi)iterator after the advance
     */
    inline const Select& next() 
    {
      // valid index?
      if ( m_index < m_size ) { ++m_index; }
      // the last combination
      if ( m_index == m_size ) { return invalidate(); }
      // evaluate the value of the current iterator from the current index
      size_t prev  = 1;
      auto   range = m_ranges  .begin();
      auto   curr  = m_current .begin();
      for ( ; m_current.end() != curr; ++curr, ++range ) 
      {
        size_t index = ( m_index / prev ) % range->size();
        prev *= range->size();
        *curr = std::next( range->begin(), index );
      }
      // the end
      return m_current;
    }
    // ========================================================================
  public:
    // ========================================================================
    /// copy the content of current iterator (with dereferencing!)
    template <class OUTPUT>
    OUTPUT current ( OUTPUT out ) const 
    {
      // copy only the valid combinations
      if ( !valid() ) { return out; }
      std::transform( m_current.begin(), m_current.end(), out, []( const iterator& i ) { return *i; } );
      return out; // RETURN
    }
    // ========================================================================
  private:
    // ========================================================================
    /** invalidate the current iterator (and index)
     *  @return corrent iterator
     */
    const Select& invalidate() 
    {
      // invalidate the current index
      m_index = m_size;
      // invalidate current  (multy) iterator to global 'end' position
      std::transform( m_ranges.begin(), m_ranges.end(), m_current.begin(), []( const Range& r ) { return r.end(); } );
      // return invalid iterator
      return m_current;
    }
    // ========================================================================
  private:
    // ========================================================================
    /// Combiner_ itself
    Ranges m_ranges   {   } ; // Combiner_ itself
    /// "current" iterator
    Select m_current  {   } ; // "current" iterator
    /// total number of combinations
    size_t m_size     { 0 } ; // total number of combinations
    /// index of current combination
    size_t m_index    { 0 } ; // index of current combination
    // ========================================================================
  }; //                                       The end of class Ostap::Combiner_
  // ==========================================================================
} //                                                 The enf of namespace Ostap
// ============================================================================
// The END
// ============================================================================
#endif // LOKI_Combiner_H
