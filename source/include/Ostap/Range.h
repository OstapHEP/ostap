#ifndef OSTAP_RANGE_H
#define OSTAP_RANGE_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <type_traits>
#include <algorithm>
#include <utility>
#include <vector>
// ============================================================================
/** @file Ostap/Ranhe.h Range.h
 *
 *  This file has been imported from
 *  <a href="http://savannah.cern.ch/projects/loki">LoKi project</a>
 *  <a href="http://cern.ch/lhcb%2Dcomp/Analysis/Loki">
 *  "C++ ToolKit  for Smart and Friendly Physics Analysis"</a>
 *
 *  The package has been designed with the kind help from
 *  Galina PAKHLOVA and Sergey BARSUK.  Many bright ideas,
 *  contributions and advises from G.Raven, J.van Tilburg,
 *  A.Golutvin, P.Koppenburg have been used in the design.
 *
 *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
 *  @date 2001-01-23
 */
// ============================================================================
namespace Ostap  
{
  // ==========================================================================
  /** @class Range_ Range.h Ostap/Range.h
   *
   *  Useful class for representation of "sequence" of the objects
   *  through the range of valid iterators.
   *
   *  The range could be created over *ALL* container types which
   *  supports at least bidirectional iterators.
   *
   *  The minimum requirements from the container type:
   *    - support the concept of "CONTAINER::value_type"
   *    - support the concept of "CONTAINER::const_iterator"
   *    - support the concept of "CONTAINER::const_reference"
   *    - support the concept of "CONTAINER::const_reverse_iterator"
   *    - the iterator should be ok for "std::distance" and "std::advance"
   *    - support for "const_iterator         CONTAINER::begin  () const"
   *    - support for "const_iterator         CONTAINER::end    () const"
   *    - support for "const_reverse_iterator CONTAINER::rbegin () const"
   *    - support for "const_reverse_iterator CONTAINER::rend   () const"
   *
   *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
   *  @date   2002-07-12
   */
  template <class CONTAINER>
  class Range_ 
  {
  public:
    // =========================================================================
    /// type for actual contained iterator
#if defined ( __cplusplus ) && defined ( __cpp_lib_remove_cvref )
    // =========================================================================
    typedef typename std::remove_cvref<CONTAINER>::type  Container      ;
    // =========================================================================
#else // =======================================================================
    // =========================================================================
    typedef typename CONTAINER                           Container      ;
    // =========================================================================
#endif // ======================================================================
    // =========================================================================
    typedef typename Container::const_iterator           iterator       ;
    typedef typename Container::const_iterator           const_iterator ;
    // =========================================================================
  private:
    // =========================================================================
    typedef typename std::iterator_traits<iterator> iter_traits;
    // =========================================================================
  public:
    // =========================================================================
    typedef typename iter_traits::value_type value_type;
    typedef typename iter_traits::reference  reference;
    typedef typename iter_traits::reference  const_reference;
    // =========================================================================
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<iterator> const_reverse_iterator;
    // ========================================================================
  public:
    // ========================================================================
    /// default constructor
    Range_() = default;
    // =========================================================================
    /** Constructor
     *  @param ibegin  iterator to begin of the sequence
     *  @param iend    iterator to end   of the sequence
     */
    template <typename InputIterator>
    Range_( InputIterator first, InputIterator last ) 
    : m_begin ( first ) 
    , m_end   ( last  ) 
    {
      if ( 0 > std::distance ( m_begin , m_end ) ) { std::swap ( m_begin , m_end ) ; }
    }
    // =========================================================================
    /** constructor from the container
     *  @param cnt  reference to the container
     */
    Range_( const Container& cnt ) 
    : m_begin ( cnt.begin () )
    , m_end   ( cnt.end   () ) 
    {}
    // =========================================================================
    /* constructor of empty range/sequence
     * @param ibegin  iterator to begin of empty sequence
     */
    Range_( iterator ibegin ) 
    : m_begin ( ibegin )
    , m_end   ( ibegin ) 
    {}
    // =========================================================================
    public:
    // =========================================================================
    /// empty sequence ?
    inline bool           empty    () const { return m_begin == m_end ; }
    /// size of the sequence (number of elements)
    inline std::size_t    size     () const { return std::distance ( m_begin , m_end ); }
    /// access to begin of the sequence (const version )
    inline iterator       begin    () const { return m_begin ; }
    /// access to end   of the sequence (const version)
    inline iterator       end      () const { return m_end ; }
    /// access to begin of the sequence (const version )
    inline const_iterator cbegin   () const { return m_begin ; }
    /// access to end   of the sequence (const version)
    inline const_iterator cend     () const { return m_end ; }
    /// access to begin of the reversed sequence (const)
    inline reverse_iterator rbegin () const { return reverse_iterator( end() ); }
    /// access to begin of the reversed sequence (const)
    inline reverse_iterator rend   () const { return reverse_iterator( begin() ); }
    /// access for the first element (only for non-empty ranges!)
    inline const_reference front   () const { return *begin(); }
    /// access for the back  element (only for non-empty ranges!)
    inline const_reference back    () const { return *std::prev ( end() ); }
    // ========================================================================
    /// get a "slice" of a range, in Python style
    inline Range_ slice ( long index1, long index2 ) const 
    {
      // trivial cases
      if ( empty() || index1 == index2 ) { return Range_(); } // RETURN
      // adjust indices
      if ( index1 < 0 ) { index1 += size(); }
      if ( index2 < 0 ) { index2 += size(); }
      // check
      if ( index1 < 0 ) { return Range_(); }      // RETURN
      if ( index2 < index1 ) { return Range_(); } // RETURN

      if ( index1 > (long)size() ) { return Range_(); } // RETURN
      if ( index2 > (long)size() ) { index2 = size(); }

      // construct the slice
      return Range_( std::next( begin(), index1 ), std::next( begin(), index2 ) ); // RETURN
    }
    // ========================================================================
    /** non-checked access to the elements by index
     *  (valid only for non-empty sequences)
     *  @param index the index of the lement to be accessed
     */
    inline const_reference operator()( const std::size_t index ) const 
    { return *std::next( begin(), index ); }
    // =====================================================================================
    /** non-checked access to the elements by index
     *  (valid only for non-empty sequences)
     *  @param index the index of the lement to be accessed
     */
    inline const_reference operator[]( const std::size_t index ) const 
    { return ( *this )( index ); }
    // =====================================================================================
    /** Checked access to the elements by index
     *  (valid for all sequences)
     *  @param index the index of the element to be accessed
     */
    inline const_reference at ( const std::size_t index ) const 
    {
      if ( size ()  <= index )  { throw std::out_of_range ( "Ostap::Range-" ) ; }
      return ( *this )( index );
    }
    // ========================================================================
  public:
    // ========================================================================
    /// compare with another range
    template <class C>
    bool operator<( const Range_<C>& right ) const 
    {
      return std::lexicographical_compare( begin(), end(), right.begin(), right.end() );
    }
    /// compare with another container
    template <class ANOTHERCONTAINER>
    bool operator<( const ANOTHERCONTAINER& right ) const 
    {
      return std::lexicographical_compare( begin(), end(), right.begin(), right.end() );
    }
    // ========================================================================
  public:
    // ========================================================================
    /// equality with another range
    bool operator==( const Range_& right ) const 
    {
      if ( &right == this ) { return true; } // RETURN
      return right.size() == size() && std::equal( begin(), end(), right.begin() );
    }
    /// equality with another range type
    template <class CNT>
    bool operator==( const Range_<CNT>& right ) const 
    {
      return right.size() == size() && std::equal( begin(), end(), right.begin() );
    }
    /// compare with another container
    template <class ANOTHERCONTAINER>
    bool operator==( const ANOTHERCONTAINER& right ) const 
    {
      return right.size() == size() && std::equal( begin(), end(), right.begin() );
    }
    // ========================================================================
  public:
    // ========================================================================
    /// empty sequence?
    bool operator!() const { return empty(); }
    /// non-empty sequence?
    explicit operator bool() const { return !empty(); }
    // ========================================================================
  private:
    // ========================================================================
    /// begin-iterator 
    iterator m_begin ;  // begin-iterator 
    /// end-iterator 
    iterator m_end   ; // end-iterator 
    // ========================================================================
  }; // end of class Range_
  // ==========================================================================
} //                                                The end of namespace Ostap 
// ============================================================================
#endif // OSTAP_RANGE_H
// ============================================================================
//                                                                      The END
// ============================================================================= 
