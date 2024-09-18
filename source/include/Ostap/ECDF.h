// ============================================================================
#ifndef OSTAP_ECDF_H 
#define OSTAP_ECDF_H 1
// ============================================================================
// Include files
// ============================================================================
// STD& STD
// ============================================================================
#include <type_traits>
#include <algorithm>
#include <vector>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/ValueWithError.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ==========================================================================
    /** @class ECDF Ostap/ECDF.h
     *  Empirical cumulative distribution function 
     *  @author Vanya Belyaev
     *  @date   2023-09-0167
     */
    class ECDF
    {
    public: 
      // ======================================================================
      /// the actual type pof data 
      typedef std::vector<double>  Data ;
      // ======================================================================
    public: 
      // ======================================================================
      /** Constructor from  data
       *  data must be non-empty!
       */ 
      ECDF
      ( const Data&  data                  ,
        const bool   complementary = false ) ;
      // =======================================================================
      /** Constructor from data
       *  data must be non-empty!
       */ 
      template <class ITERATOR,
                typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                typename            = std::enable_if<std::is_convertible<value_type,double>::value> >
      ECDF  ( ITERATOR   begin                 ,
              ITERATOR   end                   ,
              const bool complementary = false ) 
        : m_data          ( begin , end   )
        , m_complementary ( complementary ) 
      { this->sort_me() ; }
      /// constructoer to create complementary/oridnary ECDF
      ECDF
      ( const ECDF&  right         ,
        const bool   complementary ) ;
      /// copy constructor
      ECDF ( const ECDF&  right ) = default ;
      /// move constructor 
      ECDF (       ECDF&& right ) = default ;
      // ======================================================================
    public: // the main method 
      // ======================================================================
      /// the main method 
      double        evaluate   ( const double x ) const ;
      /// the main method 
      inline double operator() ( const double x ) const
      { return evaluate ( x ) ; }
      // ======================================================================
      /** get the value with binomial estimate for uncertainty
       *  @see Ostap::Math::binomEff 
       */ 
      Ostap::Math::ValueWithError estimate ( const double x ) const ; 
      // ======================================================================
    public:
      // ======================================================================
      /// add a value to data container  
      unsigned long add ( const double value  ) ;
      /// add more values to data constainer 
      unsigned long add ( const Data&  values ) ; 
      /// add more values to data container 
      unsigned long add ( const ECDF&  values ) ;
      /// add a bunch oif values 
      template <class ITERATOR,
        typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
        typename            = std::enable_if<std::is_convertible<value_type,double>::value> >
        unsigned long add
        ( ITERATOR begin ,
          ITERATOR end   )
        {
          /// copy input data 
          Data tmp1 ( begin , end ) ;
          /// sort input data 
          std::sort ( tmp1.begin  () , tmp1.end  ()  ) ;
          /// prepare the output 
          Data tmp2 ( m_data.size () + tmp1.size () )  ;
          /// merge two sorted continers 
          std::merge ( m_data.begin() ,
                       m_data.end  () ,
                       tmp1.begin  () ,
                       tmp1.end    () ,
                       tmp2.begin  () ) ;
          /// swap the output with own data 
          std::swap ( m_data , tmp2 ) ;
          return m_data.size () ;
        }
      // ======================================================================
    public :
      // ======================================================================
      ECDF& __iadd__ ( const double x ) { add ( x ) ; return *this ; }
      ECDF& __iadd__ ( const Data&  x ) { add ( x ) ; return *this ; }
      ECDF& __iadd__ ( const ECDF&  x ) { add ( x ) ; return *this ; }      
      // ======================================================================
    public:
      // ======================================================================
      /// data size 
      inline unsigned long N             () const { return m_data.size () ; } 
      /// data size 
      inline unsigned long size          () const { return m_data.size () ; }
      /// access to data
      inline const Data&   data          () const { return m_data         ; }
      // ======================================================================
      /// complementary?
      inline bool          complementary () const { return m_complementary ; }
      /// minimal x-value
      inline double        xmin          () const { return m_data.front () ; } 
      /// maximal x-value
      inline double        xmax          () const { return m_data.back  () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// swap two objects 
      void swap ( ECDF& right ) ;
      // ======================================================================
    private:
      // ======================================================================
      void sort_me () ;
      // ======================================================================
    private:
      // ======================================================================
      /// container of sorted data 
      Data m_data           {}       ; // container of sorted data
      /// complementary CDF?
      bool m_complementary { false } ; // complementary CDF ? 
      // ======================================================================
    };
    // ========================================================================
    /// swap two objects 
    inline void swap ( ECDF& left , ECDF& right ) { left.swap ( right ) ; }
    // =======================================================================
  } //                                            end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap 
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_ECDF_H
// ============================================================================
