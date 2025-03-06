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
    class Matrix      ;
    class Vector      ;
    class Permutation ;
    // =========================================================================
    // GSL version major 
    std::size_t GSL_version_major () ;
    /// GSL version minor
    std::size_t GSL_version_minor () ;
    /// GSL versionmajor  x 1000 + GAL version minor  
    std::size_t GSL_version_int   () ;
    /// GSL version as string
    std::string GSL_version       () ;
    // =========================================================================
    /** @class Ostap::GSL::Matrix
     *  Internal class to  hold GSL-Matrix
     */
    class Matrix
    {
    public : 
      // =====================================================================
      struct Zero {} ;
      struct Id   {} ;
      // ======================================================================
    public : 
      // ======================================================================
      // General rectangular matrix 
      // ======================================================================
      /// allocate GSL-matrix 
      Matrix
      ( const std::size_t N1       , 
        const std::size_t N2       ) ;
      /// allocate GSL-matrix and initialize all elements to value 
      Matrix
      ( const std::size_t N1       , 
        const std::size_t N2       , 
        const double      value    ) ;
      /// allocate GSL-matrix and initialize all elements to zero 
      Matrix
      ( const std::size_t  N1      , 
        const std::size_t  N2      , 
        const Zero      /* zero */ ) ;
      /// allocate "identity" GSL-matrix
      Matrix
      ( const std::size_t  N1      , 
        const std::size_t  N2      , 
        const Id        /* id  */  ) ;
      // =====================================================================
      // Square matrix 
      // ======================================================================
      /// allocate squared GSL-matrix 
      explicit Matrix
      ( const std::size_t  N ) ;      
      /// allocate square GSL-matrix and initialize all elements to zero 
      Matrix
      ( const std::size_t  N        , 
        const Zero      /* zero  */ ) ;
      /// allocate square identity GSL-matrix 
      Matrix
      ( const std::size_t  N        ,  
        const Id        /* id  */    ) ;
      /// make a permutation matrix
      explicit Matrix ( const Permutation& ) ;
      // make diagonal matrix
      explicit Matrix ( const Vector&      ) ; 
      // =======================================================================      
    public : 
      // =======================================================================    
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
      Matrix& operator=( const Matrix&  ) ;
      /// move assignement! 
      Matrix& operator=(       Matrix&& ) ;
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
      inline double       get
      ( const std::size_t n1 , 
        const std::size_t n2 ) const 
      { return gsl_matrix_get ( m_matrix , n1 , n2 ) ; }
      /// set the matrix element 
      inline void         set
      ( const std::size_t n1    , 
        const std::size_t n2    , 
        const double      value )
      {        gsl_matrix_set ( m_matrix , n1 , n2 , value ) ; }
      // ========================================================================
    public:
      // ========================================================================
      inline double operator ()
      ( const std::size_t i ,
        const std::size_t j ) const { return get ( i , j ) ; }
      // ========================================================================
    public:
      // ========================================================================
      // #of rows 
      inline std::size_t nRows () const { return m_matrix->size1 ; }
      // #of columns 
      inline std::size_t nCols () const { return m_matrix->size2 ; }      
      // ========================================================================
    public:
      // ========================================================================
      /// resize/reset matriz
      Matrix& resize
      ( const std::size_t n1     ,
        const std::size_t n2     ) ;
      /// resize/reset matriz
      Matrix& resize
      ( const std::size_t n1     ,
        const std::size_t n2     ,
        const double      value  ) ;
      /// resize/reset matriz
      Matrix& resize
      ( const std::size_t n1     ,
        const std::size_t n2     ,
        const Zero     /* zero */ ) ;
      /// resize/reset matriz
      Matrix& resize
      ( const std::size_t n1     ,
        const std::size_t n2     ,
        const Id      /* id */   ) ;         
      // ========================================================================
    public: // simplest mathops  
      // ========================================================================
      /// scale matrix
      Matrix&        imul ( const double  value ) ;
      /// multiply matrices using CBLAS dgemm function 
      Matrix&        imul ( const Matrix& value ) ;
      /// add matrix 
      Matrix&        iadd ( const Matrix& value ) ;
      /// add I*vlue matrix 
      Matrix&        iadd ( const double  value ) ;   
      /// subtract matrix  
      Matrix&        isub ( const Matrix& value ) ;      
      /// subtract I*value  matrix 
      inline Matrix& isub ( const double  value ) { return iadd ( -value  ) ; }
      /// scale matrix
      inline Matrix& idiv ( const double  value ) { return imul ( 1/value ) ; } 
      // ========================================================================
    public:
      // ========================================================================
      /// multiply  matrices using CBLAS dgemm function 
      Matrix multiply ( const Matrix&      right ) const ;
      /// multiply matrix & vector using CBLAS dgemv function  
      Vector multiply ( const Vector&      right ) const ;
      /// multiply matrix & permutation 
      Matrix multiply ( const Permutation& right ) const ;
      // ========================================================================
    public: // simplest mathops  
      // ========================================================================
      /// add matrix 
      inline Matrix& operator+=( const Matrix& right ) { return iadd (   right ) ; }
      /// subtract matrix 
      inline Matrix& operator-=( const Matrix& right ) { return isub (   right ) ; } 
      /// multiply matrix 
      inline Matrix& operator*=( const Matrix& right ) { return imul (   right ) ; }
      /// scale matrix 
      inline Matrix& operator*=( const double  right ) { return imul (   right ) ; } 
      /// scale matrix 
      inline Matrix& operator/=( const double  right ) { return imul ( 1/right ) ; } 
      /// add      "right*identity" matrix
      inline Matrix& operator+=( const double  right ) { return iadd (   right ) ; } 
      /// subtract "right*identity" matrix
      inline Matrix& operator-=( const double  right ) { return isub (   right ) ; } 
      // ========================================================================
    public:
      // ========================================================================
      /// transpose the matrix
      Matrix T () const ; 
      /// transose the matrix
      inline Matrix transpose () const { return T() ; }
      // ========================================================================
    public:
      // ========================================================================
      /// swap two rows in the matrix 
      Matrix& swap_rows
      ( const std::size_t i1 ,
        const std::size_t i2 ) ;      
      /// swap two columns in the matrix 
      Matrix& swap_cols 
      ( const std::size_t i1 ,
        const std::size_t i2 ) ;      
      // ========================================================================
    public: 
      // ========================================================================
      /// permute the rows     of the ematrix according to permutation 
      Matrix& permute_rows ( const Permutation& p ) ;
      /// permute the coluimns of the ematrix according to permutation 
      Matrix& permute_cols ( const Permutation& p ) ;      
      // ========================================================================
    public:
      // ========================================================================
      /// Are all elements finite ? 
      bool isfinite () const ; // Are all elements finite ? 
      /// Are all elements numerically equal to zero?      
      bool iszero   () const ; // Are all elements numerically equal to zero?
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
    } ; 
    //                                    The end of class Ostap::GSL:Matrix 
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
      ( const std::size_t N     ) ;
      /// allocate vector 
      Vector
      ( const std::size_t N     , 
        const double       value ) ;
      /// allocate vector 
      Vector
      ( const std::size_t N     , 
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
      inline double      get
      ( const std::size_t n ) const 
      { return gsl_vector_get ( m_vector , n ) ; }
      /// set the vector element 
      inline void        set
      ( const std::size_t n , 
        const double       value    )
      {        gsl_vector_set ( m_vector , n , value ) ; }
      // ========================================================================
    public:
      // ========================================================================
      // get the element 
      inline double      operator() ( const std::size_t i ) const { return get ( i ) ; }
      // get the element 
      inline double      operator[] ( const std::size_t i ) const { return get ( i ) ; }
      /// the size of the vector 
      inline std::size_t size    () const { return m_vector -> size ; } 
      // ========================================================================
    public:
      // ========================================================================
      // resie the vector
      Vector& resize
      ( const std::size_t n     ) ; 
      // resie the vector
      Vector& resize
      ( const std::size_t n     ,
        const double       value ) ; 
      // resie the vector
      Vector& resize
      ( const std::size_t n     ,
        const Zero    /* zero */ ) ;
      // ========================================================================
    public: // some simple math oprations 
      // ========================================================================
      /// add a vector of the same size 
      Vector& iadd      ( const Vector&  value ) ;
      /// add a constant 
      Vector& iadd      ( const double   value ) ;      
      /// add a vector of the same size 
      Vector& isub      ( const Vector&  value ) ;
      /// scale vector 
      Vector& imul      ( const double   value ) ;
      /// multiply by the matrix 
      Vector& imul      ( const Matrix&  value ) ;      
      /// scale vector a
      inline Vector& idiv ( const double value ) { return imul ( 1/value ) ; }
      /// subtact a constant 
      inline Vector& isub ( const double value ) { return iadd ( -value  ) ; }      
      // ========================================================================
    public : 
      // ========================================================================
      // multiplu by matrix 
      Vector multiply ( const Matrix& value ) const ;
      // ========================================================================
    public : 
      // ========================================================================      
      /// dot-product   of two vectors
      double dot      ( const Vector& right ) const ;
      /// cross-product of two vectors 
      Matrix cross    ( const Vector& right ) const ;      
      // ========================================================================      
    public: // some simple math oprations 
      // ========================================================================
      /// add vector 
      inline Vector& operator+=( const Vector& value ) { return iadd (  value ) ; }
      /// add constant 
      inline Vector& operator+=( const double  value ) { return iadd (  value ) ; }
      /// subtract vector 
      inline Vector& operator-=( const Vector& value ) { return isub (  value ) ; }
      /// subtract constant 
      inline Vector& operator-=( const double  value ) { return isub (  value ) ; }
      /// multiply with marix 
      inline Vector& operator*=( const Matrix& value ) { return imul (  value ) ; }
      /// scale it 
      inline Vector& operator*=( const double  value ) { return imul (  value ) ; }
      /// scale it 
      inline Vector& operator/=( const double  value ) { return idiv (  value ) ; }
      // ========================================================================
    public:
      // ========================================================================
      /// Are all elements finite ? 
      bool isfinite () const ; // Are all elements finite ? 
      /// Are all elements numerically equal to zero?      
      bool iszero   () const ; // Are all elements numerically equal to zero?
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
      ( const std::size_t N ) ;
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
      inline std::size_t get
      ( const std::size_t n ) const 
      { return gsl_permutation_get ( m_permutation , n ) ; }
      // ========================================================================
    public:
      // ========================================================================
      /// valid permutation ? 
      bool valid () const ; 
      // ========================================================================
    public:
      // ========================================================================
      /// apply permutation to the matrix 
      Matrix apply ( const Matrix& value ) const ;
      // ========================================================================
    public:
      // ========================================================================
      // get the element 
      inline std::size_t operator() ( const std::size_t i ) const { return get ( i ) ; }
      // get the element 
      inline std::size_t operator[] ( const std::size_t i ) const { return get ( i ) ; }
      /// the size of the vector 
      inline std::size_t size () const { return m_permutation -> size ; } 
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
    /// add two matrices 
    inline Matrix operator+( const Matrix& a , const Matrix& b )
    { Matrix c { a } ; c += b ; return c ; }
    
    /// subtract two matrices 
    inline Matrix operator-( const Matrix& a , const Matrix& b )
    { Matrix c { a } ; c -= b ; return c ; }
    
    /// add     b*identity 
    inline Matrix operator+( const Matrix& a , const double& b  )
    { Matrix c { a } ; c += b ; return c ; }
    
    /// subtract  b*identity 
    inline Matrix operator-( const Matrix& a , const double& b  )
    { Matrix c { a } ;  c-= b ; return c ; }
    
    /// multiply two matrices 
    inline Matrix operator*( const Matrix& a , const Matrix& b )
    { return a.multiply ( b ) ; }
    
    /// multiply matrix and vector 
    inline Vector operator*( const Matrix& a , const Vector& b )
    { return a.multiply ( b ) ; }

    /// scale the matrix
    inline Matrix operator*( const Matrix& a , const double b  )
    { Matrix c { a } ; c*= b ; return c ;; }
    
    /// scale the matrix
    inline Matrix operator/( const Matrix& a , const double b  )
    { Matrix c { a } ; c/= b ; return c ;; }
    
    // "right"  forms 

    /// scale the matrix from the right 
    inline Matrix operator*( const double  b , const Matrix& a )
    { Matrix c { a } ; c*= b ; return c ;; }    
    /// add     b*identity from rrght 
    inline Matrix operator+( const double& b , const Matrix& a )
    { Matrix c { a } ; c += b ; return c ; }


    /// add two vectors 
    inline Vector operator+( const Vector& a , const Vector& b )
    { Vector c { a } ; c += b ; return c ; }

    /// subtracy two vectors 
    inline Vector operator-( const Vector& a , const Vector& b )
    { Vector c { a } ; c -= b ; return c ; }

    /// add constant 
    inline Vector operator+( const Vector& a , const double  b )
    { Vector c { a } ; c += b ; return c ; }

    /// subtract constant 
    inline Vector operator-( const Vector& a , const double  b )
    { Vector c { a } ; c -= b ; return c ; }

    /// multiply vector and matrix 
    inline Vector operator*( const Vector& a , const Matrix& b )
    { return a.multiply ( b ) ; }

    /// scale the vector
    inline Vector operator*( const Vector& a , const double b  )
    { Vector c { a } ; c*= b ; return c ;; }

    /// scale the vector
    inline Vector operator/( const Vector& a , const double b  )
    { Vector c { a } ; c/= b ; return c ;; }
    
    // "right"  forms 
    
    /// scale the vector  from the right 
    inline Vector operator*( const double  b , const Vector& a )
    { Vector c { a } ; c*= b ; return c ;; }    
    /// add constant from rrght 
    inline Vector operator+( const double& b , const Vector& a )
    { Vector c { a } ; c += b ; return c ; }

    // ========================================================================
    // Matrix & permutations 
    // ========================================================================
    
    /// Apply permutation 
    inline Matrix operator*( const Permutation& p ,
			     const Matrix&      a ) { return p.apply    ( a ) ; }
    /// apply permutation 
    inline Matrix operator*( const Matrix&a       ,
			     const Permutation& p ) { return a.multiply ( p ) ; }
    
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
    Permutation PLU
    ( Matrix& A ) ;
    // ========================================================================
    /** perfom LU decomposition  
     *  @param  A   (INOUT)         input matrix 
     *  @param  LU  (UPDATE/OUTPUT) output LU matrix 
     *  @return M-permutation 
     *  @see gsl_linalg_LU_decomp 
     */
    Permutation PLU
    ( const Matrix& A  ,
      Matrix&       LU ) ;
    // ========================================================================
    /** perfom LU decomposition  
     *  @param  A   (INOUT)         input matrix 
     *  @param  L   (UPDATE/OUTPUT) lower triangular matrix 
     *  @param  U   (UPDATE/OUTPUT) upper triangular matrix 
     *  @return M-permutation 
     *  @see gsl_linalg_LU_decomp 
     */
    Permutation PLU
    ( const Matrix& A ,
      Matrix&       L ,
      Matrix&       U ) ;    
    // ========================================================================
    // QR decomposition with column pivoting 
    // ========================================================================
    /** mape QR Decomposion of matrix A : \f$ AP = QR\f$ where 
     *  - A is input                 MxN matrix  
     *  - P is permuutation matrix   NxN 
     *  - Q is orthogonal matrix     MxM 
     *  - R is right triaular matrix MxN 
     *  
     *  @param A  (input) the matrix to decopose 
     *  @param Q  (outpt/update) orthogonal matrix Q 
     *  @param R  (outpt/update) rigth triangular matrix R 
     *  @return permutation P 
     */
    Permutation PQR
    ( const Matrix& A ,
      Matrix&       Q ,
      Matrix&       R ) ;    
    // ======================================================================
    // LQ decomposition
    // ======================================================================
    /** LQ decomposition of matrix A: \f$ A = LQ\f$, where 
     *  - L is lower trapezoidal MxN 
     *  - Q is orthogonal NxN 
     */ 
    void LQ
    ( const Matrix& A ,
      Matrix&       L ,
      Matrix&       Q ) ;
    // ======================================================================
    // QL decomposition
    // ======================================================================
    /** QL decomposition of matrix A: \f$ A = QL\f$, where 
     *  - Q is orthogonal MxM
     *  - L is lower trapezoidal MxN 
     */ 
    void QL
    ( const Matrix& A ,
      Matrix&       Q ,
      Matrix&       L ) ;
    // ======================================================================
    // COD decomposition
    // ======================================================================
    /** COD - Complete Orthogonal Decomposion
     *  \f$ AP = Q R Z^T \f$ 
     *  - A input MxN matrix 
     *  - P is permutation matrix 
     *  - Q is MxM orthogonal matrix 
     *  - Z is NxN orthogonal matrix 
     *  - R is 2x2 block matrix with top-left blobck being right triangular matrix and
     *    other blocks are zeroes 
     */
    Permutation COD
    ( const Matrix& A ,
      Matrix& Q ,
      Matrix& R ,
      Matrix& Z ) ;
    // ========================================================================
    // SVD decomposition
    // ========================================================================
    /** SVD : singular Value Decomposition  \f$ A = U S V^T\f$
     *  - A input MxN matrix 
     *  - K = min ( M , N ) : 
     *  - U MxK orthogonal matrix 
     *  - S KxK Diagonal matrix of singular values 
     *  - V NxK orthogonal matrix 
     *  @param A     (input)  input matrix A 
     *  @param U     (update) orthogonal matrix U 
     *  @param V     (update) orthogonal matrix V 
     *  @param golub (input) use Golub or Jacobi algorithm 
     *  @return vector of singular values 
     * -  Jacobi algorithm is more prrcise  and Golub algorithm is more CPU efficient 
     */
    Vector SVD
    ( const Matrix& A            ,
      Matrix&       U            ,
      Matrix&       V            ,
      const bool    golub = true ) ;
    // ========================================================================

    // ========================================================================
    // Schur decomposiiton of squre matrix
    // ========================================================================
    /** Schur decomposition of square matrix \f$ A = Z T Z^T\f$, where 
     *  - A is inpur MxM (square) matrix
     *  - T is Schur form of matix  
     *  - Z is orthogonam matrix 
     */
    void SCHUR 
    ( const Matrix&  A ,  
      Matrix&        Z , 
      Matrix&        T ) ; 
      
    // ========================================================================
    // Polar decomposition of square matrix 
    // ========================================================================
    /** Polar decompositon of the square matrix A: \f$ A = UP \f$
     *  - U ius orthogonal 
     *  - P is positive semi-definitive 
     */
    void POLAR
    ( const Matrix& A ,
      Matrix      & U ,
      Matrix      & P ) ;
    // ========================================================================    
     
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
    /// 
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
