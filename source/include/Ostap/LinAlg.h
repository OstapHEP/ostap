// ============================================================================
#ifndef OSTAP_GSL_LINALG_H 
#define OSTAP_GSL_LINALG_H 1
// ============================================================================
// Include files
// ============================================================================
// STD/STD
// ============================================================================
#include <ostream>
// ============================================================================
// GSL 
// ============================================================================
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_permutation.h"
// =============================================================================
namespace Ostap 
{
  // ==========================================================================
  namespace GSL
  {
    // =========================================================================
    /** @class Ostap::GSL::Matrix
     *  Internal class to  hold GSL-Matrix
     */
    class Matrix
    {
    public : 
      // =====================================================================
      struct Zero     {} ;
      struct Identity {} ;
      // ======================================================================
    public : 
      // ======================================================================
      // General rectangular matrix 
      // ======================================================================
      /// allocate GSL-matrix 
      Matrix
      ( const unsigned int N1        , 
        const unsigned int N2        ) ;
      /// allocate GSL-matrix and initialize all elements to value 
      Matrix
      ( const unsigned int N1        , 
        const unsigned int N2        , 
        const double       value     ) ;
      /// allocate GSL-matrix and initialize all elements to zero 
      Matrix
      ( const unsigned int N1        , 
        const unsigned int N2        , 
        const Zero         /* zero */ ) ;
      /// allocate "identity" GSL-matrix
      Matrix
      ( const unsigned int N1        , 
        const unsigned int N2        , 
        const Identity        /* id  */ ) ;
      // =====================================================================
      // Square matrix 
      // ======================================================================
      /// allocate squared GSL-matrix 
      Matrix
      ( const unsigned int N ) ;      
      /// allocate square GSL-matrix and initialize all elements to zero 
      Matrix
      ( const unsigned int N         , 
        const Zero         zero      ) ;
      /// allocate square identity GSL-matrix 
      Matrix
      ( const unsigned int N    , 
        const Identity     id    ) ;
      /// copy constructor 
      Matrix  ( const Matrix&  right ) ;
      /// move constructor 
      Matrix  (       Matrix&& right ) ;
      ///  destructor: free  GSL-matrix 
      ~Matrix () ;
      // ========================================================================
      /// no default constructor 
      Matrix  () = delete ;
      /// copy assignement! 
      Matrix& operator= ( const Matrix&  ) ;
      /// move assignement! 
      Matrix& operator= (       Matrix&& ) ;
      // ========================================================================
    public:
      // ========================================================================
      // get the matrix
      inline       gsl_matrix* matrix  ()       { return m_matrix ; }
      // get the matrix
      inline const gsl_matrix* matrix  () const { return m_matrix ; }
      // ========================================================================
    public :
      // ========================================================================
      /// get the matrix element 
      double      get
      ( const unsigned int n1 , 
        const unsigned int n2 ) const 
      { return gsl_matrix_get ( m_matrix , n1 , n2 ) ; }
      /// set the matrix element 
      void        set
      ( const unsigned int n1    , 
        const unsigned int n2    , 
        const double       value )
      {        gsl_matrix_set ( m_matrix , n1 , n2 , value ) ; }
      // ========================================================================
    public:
      // ========================================================================
      inline double operator ()
      ( const unsigned int i ,
        const unsigned int j ) const { return get ( i , j ) ; }
      // ========================================================================
    public:
      // ========================================================================
      // #of rows 
      unsigned int nRows () const { return m_matrix->size1 ; }
      // #of columns 
      unsigned int nCols () const { return m_matrix->size2 ; }      
      // ========================================================================
    public:
      // ========================================================================
      /// resize/reset matriz
      Matrix& resize
      ( const unsigned int n1     ,
        const unsigned int n2     ) ;
      /// resize/reset matriz
      Matrix& resize
      ( const unsigned int n1     ,
        const unsigned int n2     ,
        const double       value  ) ;
      /// resize/reset matriz
      Matrix& resize
      ( const unsigned int n1     ,
        const unsigned int n2     ,
        const Zero     /* zero */ ) ;
      /// resize/reset matriz
      Matrix& resize
      ( const unsigned int n1     ,
        const unsigned int n2     ,
        const Identity /* id */   ) ;         
      // ========================================================================
    public:
      // ========================================================================
      /// swap two matrices 
      void swap ( Matrix& right ) ; // swap two matrices 
      // ========================================================================
    private:
      // ========================================================================
      /// the  actual pointer to GSL-matrix 
      gsl_matrix* m_matrix { nullptr } ; // the  actual pointer to GSL-matrix
      // ========================================================================
    };
    // ==========================================================================
    /** @class Ostap::GSL::Vector
     *  Internal class to  hold GSL-Vector
     */
    class Vector
    {
    public: 
      // ========================================================================
      typedef Ostap::GSL::Matrix::Zero Zero ;
      // ========================================================================
    public: 
      // ========================================================================
      /// allocate vector 
      Vector
      ( const unsigned int N     ) ;
      /// allocate vector 
      Vector
      ( const unsigned int N     , 
        const double       value ) ;
      /// allocate vector 
      Vector
      ( const unsigned int N     , 
        const Zero     /* zero */  ) ;
      /// copy constructor 
      Vector
      ( const Vector&  right ) ;
      /// move constructor 
      Vector
      (       Vector&& right ) ;
      /// destructor : free GSL-Vector
      ~Vector () ;
      // =======================================================================
      /// no default constructor 
      Vector() = delete ;
      /// copy assignement! 
      Vector& operator= ( const Vector&  ) ;
      /// move assignement! 
      Vector& operator= (       Vector&& ) ;
      // ========================================================================
    public:
      // ========================================================================
      // get the vector
      inline       gsl_vector* vector ()       { return m_vector ; }
      // get the vector
      inline const gsl_vector* vector () const { return m_vector ; }
      // ========================================================================
      /// get the vector element 
      double      get
      ( const unsigned int n ) const 
      { return gsl_vector_get ( m_vector , n ) ; }
      /// set the vector element 
      void        set
      ( const unsigned int n , 
        const double       value    )
      {        gsl_vector_set ( m_vector , n , value ) ; }
      // ========================================================================
    public:
      // ========================================================================
      // get the element 
      inline double operator() ( const unsigned int i ) const { return get ( i ) ; }
      // get the element 
      inline double operator[] ( const unsigned int i ) const { return get ( i ) ; }
      /// the size of the vector 
      std::size_t   size    () const { return m_vector->size ; } 
      // ========================================================================
    public:
      // ========================================================================
      // resie the vector
      Vector& resize
      ( const unsigned int n     ) ; 
      // resie the vector
      Vector& resize
      ( const unsigned int n     ,
        const double       value ) ; 
      // resie the vector
      Vector& resize
      ( const unsigned int n     ,
        const Zero    /* zero */ ) ;
      // ========================================================================
    public:
      // ========================================================================
      /// swap two vectors 
      void swap ( Vector& right ) ; // swap two vetcors 
      // ========================================================================
    private:
      // ========================================================================
      /// the  actual pointer to GSL-vector 
      gsl_vector* m_vector { nullptr } ; // the  actual pointer to GSL-vector
      // ========================================================================
    };
    // ==========================================================================
    /** @class Ostap::GSL::Permutation
     *  Internal class to keep GSL-permuation
     */
    class  Permutation
    {
    public:
      // ========================================================================
      /// constructor: allocate the permutation 
      Permutation
      ( const unsigned int N ) ;
      /// destructor: free permutation 
      ~Permutation() ;
      // ========================================================================
      Permutation () = delete ;
      Permutation ( const Permutation&  ) ;
      Permutation (       Permutation&& ) ;
      /// copy assignement! 
      Permutation& operator= ( const Permutation&  ) ;
      /// move assignement! 
      Permutation& operator= (       Permutation&& ) ;
      // =======================================================================
    public:
      // =======================================================================
      // get the permutation 
      inline       gsl_permutation* permutation ()       { return m_permutation ; }
      // get the permutation 
      inline const gsl_permutation* permutation () const { return m_permutation ; }
      // ========================================================================
    public:
      // ========================================================================
      /// get the vector element 
      double      get
      ( const unsigned int n ) const 
      { return gsl_permutation_get ( m_permutation , n ) ; }
      // ========================================================================
    public:
      // ========================================================================
      // get the element 
      inline double operator() ( const unsigned int i ) const { return get ( i ) ; }
      // get the element 
      inline double operator[] ( const unsigned int i ) const { return get ( i ) ; }
      /// the size of the vector 
      std::size_t   size    () const { return m_permutation->size ; } 
      // ========================================================================
    public:
      // ========================================================================
      /// swap two permutations 
      void swap ( Permutation& right ) ; // swap two permutations 
      // ========================================================================
    private :
      // =====================================================================
      /// the  actual pointer to GSL-permutation
      gsl_permutation* m_permutation { nullptr } ; // the  actual pointer to GSL-vector
      // ======================================================================
    };
    // ========================================================================
    /// swap two matrices 
    inline void swap ( Matrix& a      , Matrix&      b ) { a.swap ( b ) ; } 
    /// swap two vectors 
    inline void swap ( Vector& a      , Vector&      b ) { a.swap ( b ) ; } 
    /// swap two permuattions 
    inline void swap ( Permutation& a , Permutation& b ) { a.swap ( b ) ; } 
    // ========================================================================
    
    // ========================================================================
    // Linear Algebra 
    // ========================================================================
    
    // ========================================================================
    /** "in-place" LU decomposition  
     *  \f$ PA = LU \f$, where 
     *   - A is \f$ M \times N \f$ matrix 
     *   - P is \f$ M \times M \f$ permutation matrix  
     *   - L is \f$ M \times \min (M,N)\f$  lower trianhular matris 
     *   - U is \f$ \min (M,N) \times N \f$ upper trianhular matris 
     * 
     * For square matrices:
     *   - L is a lower unit triangular matrix
     *   - U is upper triangular
     * 
     * For \f$ M>N \f$: 
     *   - L is a unit lower trapezoidal matrix of size \f$ M\timex N \f$ 
     * 
     * For \f$ M < N \f$: 
     *  - U is upper trapezoidal of size \f$ M \times M \f$  
     *
     *  For square matrices this decomposition can be used to convert the linear 
     *  system \f$ Ax=b\f$  into a pair of triangular systems, 
     *  \f$ Lu=Pb\f$ and  \f$ Ux=y\f$, which can be solved by forward and 
     *   back-substitution. 
     *  Note that the LU decomposition is also valid for singular matrices.
     *
     *  The matrix at the end contans two matrices: 
     * 
     *  On output the diagonal and upper triangular (or trapezoidal) part of the 
     *  input matrix A contain the matrix U. 
     *  The lower triangular (or trapezoidal) part of the input matrix (excluding 
     *  the diagonal) contains L. The diagonal elements of U are unity, and are not stored.
     *
     *  The permutation matrix P is encoded in the permutation p on output.
     *  The j-th column of the matrix P is given by the j-th column of the 
     *   identity matrix, where \f$ k = p_j \f%  
     *  the j-th element of the permutation vector. 
     * 
     *  @param  A  (update) input/update MxN  marix 
     *  @return M-permutation 
     *  @see gsl_linalg_LU_decomp 
     */
    Permutation PLU ( Matrix& A ) ;
    // ========================================================================
    /** perfom LU decomposition  
     *  @param  A   (INOUT)         input matrix 
     *  @param  LU  (UPDATE/OUTPUT) output LU matrix 
     *  @return M-permutation 
     *  @see gsl_linalg_LU_decomp 
     */
    Permutation PLU ( const Matrix& A , Matrix& LU ) ;
    // ========================================================================
    /** perfom LU decomposition  
     *  @param  A   (INOUT)         input matrix 
     *  @param  L   (UPDATE/OUTPUT) lower triangular matrix 
     *  @param  U   (UPDATE/OUTPUT) upper triangular matrix 
     *  @return M-permutation 
     *  @see gsl_linalg_LU_decomp 
     */
    Permutation PLU ( const Matrix& A , Matrix& L, Matrix& U  ) ;
    
    
    // ========================================================================
  } //                                          The end of namespace Ostap::GSL
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /// get the element with maxina absolute value 
    double maxabs_element ( const Ostap::GSL::Matrix&      m ) ;
    /// get the element with maxina absolute value 
    double maxabs_element ( const Ostap::GSL::Vector&      v ) ;
    /// get the element with maxina absolute value 
    double maxabs_element ( const Ostap::GSL::Permutation& p ) ;
    // ========================================================================
  } //                                        The end of namespaxce Ostap::Math 
  // ==========================================================================
  namespace Utils 
  {
    // ========================================================================
    /** print GSL-vector to the stream 
     *  @param v the vector 
     *  @param s the stream 
     *  @return the stream 
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2012-05-28
     */
    std::ostream& toStream 
    ( const gsl_vector&  v , 
      std::ostream&      s ) ;
    // ========================================================================
    /** print GSL-matrix to the stream 
     *  @param m the matrix 
     *  @param s the stream 
     *  @return the stream 
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2012-05-28
     */
    std::ostream& toStream 
    ( const gsl_matrix&  m , 
      std::ostream&      s ) ;
    // ========================================================================
    /** print GSL-permutation to the stream 
     *  @param p the permutation
     *  @param s the stream 
     *  @return the stream 
     */    
    std::ostream& toStream 
    ( const gsl_permutation& p , 
      std::ostream&          s ) ;
    // ========================================================================
    /** print GSL-matrix to the stream 
     *  @param m the matrix 
     *  @param s the stream 
     *  @return the stream 
     */
    std::ostream& toStream
    ( const Ostap::GSL::Matrix& m ,
      std::ostream&             s ) ;
    // ========================================================================
    /** print GSL-vector to the stream 
     *  @param v the vector 
     *  @param s the stream 
     *  @return the stream 
     */
    std::ostream& toStream
    ( const Ostap::GSL::Vector& v ,
      std::ostream&             s ) ;
    // ========================================================================
    /** print GSL-permutation to the stream 
     *  @param p the permutation
     *  @param s the stream 
     *  @return the stream 
     */
    std::ostream& toStream
    ( const Ostap::GSL::Permutation& p ,
      std::ostream&                  s ) ;
    // ========================================================================
  } //                                       The end of name space Ostap::Utils 
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
/// print operator 
inline std::ostream& operator<<
( std::ostream&     s ,
  const gsl_vector& v ) 
{ return Ostap::Utils::toStream ( v , s ) ; }
// ============================================================================
/// print operator 
inline std::ostream& operator<<
( std::ostream&     s ,
  const gsl_matrix& m ) 
{ return Ostap::Utils::toStream ( m , s ) ; }
// ============================================================================
/// print operator 
inline std::ostream& operator<<
( std::ostream&             s ,
  const Ostap::GSL::Matrix& m ) 
{ return Ostap::Utils::toStream ( m , s ) ; }
// ============================================================================
/// print operator 
inline std::ostream& operator<<
( std::ostream&             s ,
  const Ostap::GSL::Vector& v ) 
{ return Ostap::Utils::toStream ( v , s ) ; }
// ============================================================================
#endif // OSTAP_GSL_LINALG_H
// ============================================================================
//                                                                      The END 
// ============================================================================
