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
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Span.h"
#include "Ostap/Types.h"
#include "Ostap/Buffer.h"
#include "Ostap/Math.h"
#include "Ostap/StatusCode.h"
// ============================================================================
/** @file Ostap/LinAlg.h
 *  File provides utilities to access to the basic Linear Algebra from GSL
 *
 *  First it provides ilght wrappers for GSL matrices, vectors and permutaitons
 *  with all basic operations included
 *
 *  @see class Ostap::Math::GSL::Matrix
 *  @see class Ostap::Math::GSL::Vectors
 *  @see class Ostap::Math::GSL::Permutation
 *  @see https://www.gnu.org/software/gsl/doc/html/vectors.html
 *
 *  The basic Linear Alegbra functinons are:
 *  - LU decomposition @see https://www.gnu.org/software/gsl/doc/html/linalg.html#lu-decomposition
 *  - PQR: QR decomposition with column pivoting  @see https://www.gnu.org/software/gsl/doc/html/linalg.html#qr-decomposition-with-column-pivoting
 *  - LQ decomposition @see https://www.gnu.org/software/gsl/doc/html/linalg.html#lq-decomposition
 *  - QL decomposition @see https://www.gnu.org/software/gsl/doc/html/linalg.html#ql-decomposition
 *  - COD: Complete Orthogonal Decomposition @see  https://www.gnu.org/software/gsl/doc/html/linalg.html#complete-orthogonal-decomposition
 *  - SVD: Singular Value Decomposition @see  https://www.gnu.org/software/gsl/doc/html/linalg.html#singular-value-decomposition
 *  - LLT: Cholesky decomposition @see https://www.gnu.org/software/gsl/doc/html/linalg.html#cholesky-decomposition
 *  - LDLT: Pivoted Cholesky Decompositon with scale factor @see https://www.gnu.org/software/gsl/doc/html/linalg.html#pivoted-cholesky-decomposition
 *  - D3: Trigiagonal decomposition of symmetric matrices @see https://www.gnu.org/software/gsl/doc/html/linalg.html#tridiagonal-decomposition-of-real-symmetric-matrices
 *  - UHUT: Hessenberg decomposition of matrices @see https://www.gnu.org/software/gsl/doc/html/linalg.html#hessenberg-decomposition-of-real-matrices
 *  - UBVT: Bidiagolalization of general matrices @see https://www.gnu.org/software/gsl/doc/html/linalg.html#bidiagonalization
 *  @see https://www.gnu.org/software/gsl/doc/html/linalg.html
 */
// ============================================================================
namespace Ostap 
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    namespace GSL
    {
      // ======================================================================
      // Forward declarations 
      // ======================================================================
      class Matrix      ;
      class Vector      ;
      class Permutation ;
      // ======================================================================
      // GSL version major 
      std::size_t GSL_version_major () ;
      /// GSL version minor
      std::size_t GSL_version_minor () ;
      /// GSL versionmajor  x 1000 + GSL version minor  
      std::size_t GSL_version_int   () ;
      /// GSL version as string
      std::string GSL_version       () ;
      // ======================================================================
      /** @class Ostap::GSL::Matrix
       *  Internal class to hold GSL-Matrix
       */
      class Matrix
      {
      public : 
        // ====================================================================
        struct Zero {} ;
        struct Id   {} ;
        // ====================================================================
      public : // iterator 
        // ====================================================================
        /** @class const_iterator
         *  Row-major forward iterator for GSL matrix
         */
        class const_iterator 
        {
        public:
          // =================================================================
          using iterator_category = std::forward_iterator_tag ;
          using value_type        = double ;
          using difference_type   = std::ptrdiff_t ;
          using pointer           = const double* ;
          using reference         = double ;
          // =================================================================          
        public:
          // =================================================================
          /// create iterator 
          const_iterator 
          ( const gsl_matrix* m   = nullptr , 
            const std::size_t row = 0       , 
            const std::size_t col = 0       )
            : m_mat ( m   ) 
            , m_row ( row ) 
            , m_col ( col ) {}
          /// dereference iterator 
          inline reference operator* () const 
          { return gsl_matrix_get ( m_mat , m_row , m_col ) ;  }
          /// advance iterator 
          inline const_iterator& operator++ () 
          {
            if ( !m_mat ) { return *this ; }
            ++m_col ;
            if ( m_mat->size2 <= m_col ) 
            {
              m_col = 0 ;
              ++m_row ;
            }
            return *this ;
          }
          /// advance iterator 
          inline const_iterator operator++ ( int ) 
          {
            const_iterator tmp = *this ;
            ++(*this) ;
            return tmp ;
          }
          /// iterator equality 
          inline bool operator== ( const const_iterator& other ) const 
          {
            if ( m_mat != other.m_mat ) return false ;
            if ( m_mat->size1 <= m_row  && other.m_mat->size1 <= other.m_row ) { return true ; } 
            return m_row == other.m_row && m_col == other.m_col ;
          }
          /// iterator non-equality 
          inline bool operator!= ( const const_iterator& other ) const 
          { return !(*this == other) ; }
          // =====================================================================
        private:
          // =====================================================================
          /// the matrix 
          const gsl_matrix* m_mat { nullptr } ; ///<  the matrix
          /// row index 
          std::size_t       m_row { 0 }       ; ///< row index
          /// column index 
          std::size_t       m_col { 0 }       ; ///< column index
          // =====================================================================
        } ;
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
        /// create a diagonal diagonal matrix 
        explicit Matrix ( const Vector&      ) ; 
        /// create a permutation matrix
        explicit Matrix ( const Permutation& ) ;
        // =======================================================================
      public : // create the matrix from shape and continuous buffer 
        // =======================================================================
        /** Create the matrix from the shape and continious buffer:
         *  - to be replaced wth std::span later
         *  
         *  For <code>N1 == N2</code> (square matrix) :
         *  if <code>buffer.size() == N1 * ( N1 + 1 ) / 2 </code> 
         *  the matrix is treated as symmetric
         */
        template <class SCALAR , 
                  typename std::enable_if<Ostap::is_convertible_to_double<SCALAR>::value, bool>::type = true>            
        Matrix
        ( const std::size_t                   N1     , 
          const std::size_t                   N2     , 
          const Ostap::Utils::Buffer<SCALAR>& buffer ) 
          : Matrix ( N1 , N2)
        { this -> fill_impl ( buffer.data () , buffer.size () ) ; }
        // ========================================================================
        /** Create square matrix from the shape and continious buffer:
         *  - to be replaced wth std::span later
         *
         *  If <code>buffer.size() == N * ( N + 1 ) / 2 </code> 
         *  the matrix is treated as symmetric
         */
        template <class SCALAR , 
                  typename std::enable_if<Ostap::is_convertible_to_double<SCALAR>::value, bool>::type = true>            
        Matrix
        ( const std::size_t                   N      ,  
          const Ostap::Utils::Buffer<SCALAR>& buffer ) 
          : Matrix ( N , N , buffer )
        {}
        // =======================================================================
#if defined(OSTAP_HAS_STD_SPAN) && OSTAP_HAS_STD_SPAN
        // ======================================================================
        /** Create the matrix from the shape and continious buffer:
         *  - to be replaced wth std::span later
         *  
         *  For <code>N1 == N2</code> (square matrix) :
         *  if <code>buffer.size() == N1 * ( N1 + 1 ) / 2 </code> 
         *  the matrix is treated as symmetric
         */
        template <class SCALAR, 
                  typename std::enable_if<Ostap::is_convertible_to_double<SCALAR>::value, bool>::type = true>            
        Matrix
        ( const std::size_t                 N1     , 
          const std::size_t                 N2     , 
          const std::span<SCALAR>&          buffer ) 
          : Matrix ( N1 , N2)
        { this -> fill_impl ( buffer.data () , buffer.size () ) ; }
        // ========================================================================
        /** Create square matrix from the shape and continious buffer:
         *  - to be replaced wth std::span later
         *
         *  If <code>buffer.size() == N * ( N + 1 ) / 2 </code> 
         *  the matrix is treated as symmetric
         */
        template <class SCALAR , 
                  typename std::enable_if<Ostap::is_convertible_to_double<SCALAR>::value, bool>::type = true>            
        Matrix
        ( const std::size_t                 N      ,  
          const std::span<SCALAR>&          buffer ) 
          : Matrix ( N , N , buffer )
        {}
        // =======================================================================
#endif  // OSTAP_HAS_STD_SPAN 
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
        /// get the matrix element
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
        // ==============================================================
        /// begin-iterator 
        const_iterator begin () const 
        { return const_iterator ( m_matrix , 0 , 0 ) ; }
        /// end-iterator  
        const_iterator end () const 
        { return const_iterator ( m_matrix , nRows() , 0 ) ; }
        /// begin-iterator 
        const_iterator cbegin () const { return begin() ; }
        /// end-iterator 
        const_iterator cend   () const { return end()   ; }
        // ========================================================================
      public : // rows and colum by value 
        // ========================================================================
        /// get matrix row   (by value)
        Vector  row    ( const std::size_t r ) const ;
        /// get matrix column (by value)
        Vector  column ( const std::size_t r ) const ;  
        // ========================================================================
      public : 
        // ========================================================================
        /// matrices are equal ?
        bool equal  ( const Matrix & right ) const ;
        // ========================================================================
      public:
        // ========================================================================
        /// resize/reset matrix
        Matrix& resize
        ( const std::size_t n1     ,
          const std::size_t n2     ) ;
        /// resize/reset matrix
        Matrix& resize
        ( const std::size_t n1     ,
          const std::size_t n2     ,
          const double      value  ) ;
        /// resize/reset matrix
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
        /// resize to square matrix 
        inline Matrix& resize
        ( const std::size_t n  ) { return resize ( n , n ) ; } 
        /// resize to square matrix 
        inline Matrix& resize
        ( const std::size_t n  ,
          const Zero        z  ) { return resize ( n , n , z  ) ; } 
        /// resize to square matrix 
        inline Matrix& resize
        ( const std::size_t n  ,
          const Id          id ) { return resize ( n , n , id ) ; } 
        // ========================================================================
      public: // simplest math operations   
        // ========================================================================
        /// scale matrix
        Matrix&        imul ( const double  value ) ;
        /// multiply matrices using CBLAS dgemm function 
        Matrix&        imul ( const Matrix& value ) ;
        /// add matrix 
        Matrix&        iadd ( const Matrix& value ) ;
        /// add I*value matrix 
        Matrix&        iadd ( const double  value ) ;   
        /// subtract matrix  
        Matrix&        isub ( const Matrix& value ) ;      
        /// subtract I*value  matrix 
        inline Matrix& isub ( const double  value ) { return iadd ( -value  ) ; }
        /// scale matrix
        inline Matrix& idiv ( const double  value ) { return imul ( 1/value ) ; } 
        // ========================================================================
      public: // multuplication 
        // ========================================================================
        /// multiply  matrices using CBLAS dgemm function 
        Matrix multiply ( const Matrix&      right ) const ;
        /// multiply matrix & vector using CBLAS dgemv function  
        Vector multiply ( const Vector&      right ) const ;
        /// multiply matrix & permutation 
        Matrix multiply ( const Permutation& right ) const ;
        // ========================================================================
      public: // simplest math operations   
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
        /// permute the rows    of the matrix according to permutation 
        Matrix& permute_rows ( const Permutation& p ) ;
        /// permute the columns of the matrix according to permutation 
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
      private :
        // ========================================================================
        /// fill the matrix from continious buffer 
        template <class TYPE>
        void fill_impl
        ( const TYPE*       buffer , 
          const std::size_t size   ) ;
        // ========================================================================
      private:
        // ========================================================================
        /// the  actual pointer to GSL-matrix 
        gsl_matrix* m_matrix { nullptr } ; //! the  actual pointer to GSL-matrix
        // ========================================================================
      } ; //                                    The end of class Ostap::GSL:Matrix 
      // ==========================================================================
      extern template void Matrix::fill_impl<float>       ( const       float* , const std::size_t ) ;
      extern template void Matrix::fill_impl<double>      ( const      double* , const std::size_t ) ;
      extern template void Matrix::fill_impl<long double> ( const long double* , const std::size_t ) ;
      // ==========================================================================
      /** @class Vector
       *  Internal class to  hold GSL-Vector
       */
      class Vector
      {
      public: 
        // ========================================================================
        typedef Ostap::Math::GSL::Matrix::Zero Zero ;
        // ========================================================================
      public :
        // ========================================================================
        /** @class const_iterator
         *  Robust, stride-aware forward iterator for GSL vector 
         *  delegating element access to gsl_vector_get
         */
        class const_iterator 
        {
        public:
          using iterator_category = std::forward_iterator_tag ;
          using value_type        = double ;
          using difference_type   = std::ptrdiff_t ;
          using pointer           = const double* ;
          using reference         = double ;
          /// Construct iterator from GSL vector pointer and logical index
          const_iterator 
          ( const gsl_vector* v     = nullptr , 
            const std::size_t index = 0 )
            : m_vec   ( v     )
            , m_index ( index )
          {}

          /// Dereference operator: fetches value via official GSL API
          reference operator* () const 
          { return gsl_vector_get ( m_vec , m_index ) ; }

          /// Prefix increment operator (++it)
          const_iterator& operator++ () 
          { if ( m_vec ) { ++m_index ; } return *this ; }

          /// Postfix increment operator (it++)
          const_iterator operator++ ( int ) 
          {
            const_iterator tmp { *this } ;
            ++(*this) ;
            return tmp ;
          }
          /// Equality operator handling logical end-of-range comparisons
          bool operator== ( const const_iterator& other ) const 
          {
            // Pointing to the exact same underlying object
            if ( m_vec == other.m_vec ) 
            { return m_index == other.m_index ; }
            // Logical end check for empty or uninitialized vectors
            const std::size_t s1 = ( m_vec       ? m_vec->size       : 0 ) ;
            const std::size_t s2 = ( other.m_vec ? other.m_vec->size : 0 ) ;
            return ( s1 <= m_index  ) && ( s2 <= other.m_index  ) ;
          }
          /// Inequality operator
          bool operator!= ( const const_iterator& other ) const 
          { return !(*this == other) ; }
          // ======================================================================
        private:
          // ======================================================================
          const gsl_vector* m_vec   { nullptr } ;
          std::size_t       m_index { 0 }       ;
          //=======================================================================
        } ;
        // ======================================================================== 
      public: 
        // ========================================================================
        /// allocate vector 
        Vector
        ( const std::size_t N     ) ;
        /// allocate vector & fill with with the value 
        Vector
        ( const std::size_t N     , 
          const double       value ) ;
        /// allocate vector and initiilze it  with zero 
        Vector
        ( const std::size_t N     , 
          const Zero     /* zero */  ) ;
        // =========================================================================
        /// templated contructor from non-empty  sequence of elements 
        template < typename ITERATOR                                                               , 
                   typename CATEGORY = typename std::iterator_traits<ITERATOR>::iterator_category  ,
                   typename VALUE    = typename std::iterator_traits<ITERATOR>::value_type         ,
                   typename          = std::enable_if_t<std::is_base_of_v<std::forward_iterator_tag, CATEGORY>> ,
                   typename          = std::enable_if_t<Ostap::is_convertible_to_double<VALUE>::value   > >
        Vector
        ( ITERATOR begin ,
          ITERATOR end   ) 
        : Vector ( std::distance ( begin , end ) )
        {
          std::size_t index = 0 ;
          for ( ; begin != end ; ++begin , ++index )
          { gsl_vector_set ( this->m_vector , index , static_cast<double> ( *begin ) ) ; }          
        }
        // =========================================================================        
        /// copy constructor 
        Vector
        ( const Vector&  right ) ;
        /// move constructor 
        Vector
        (       Vector&& right ) ;
        /// destructor : free GSL-Vector
        ~Vector () ;
        // ========================================================================
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
        // Iterator accessor methods
        // ========================================================================
        const_iterator begin () const 
        { return const_iterator ( m_vector , 0      ) ; }
        const_iterator end () const 
        { return const_iterator ( m_vector , size() ) ; }
        const_iterator cbegin () const { return begin() ; }
        const_iterator cend   () const { return end()   ; }
        // ========================================================================
      public:
        // ========================================================================
        /// vectors are equal ?
        bool equal ( const Vector& right ) const ;
        // ========================================================================
      public:
        // ========================================================================
        // resize the vector
        Vector& resize
        ( const std::size_t n     ) ; 
        // resize the vector and fill it with defined value
        Vector& resize
        ( const std::size_t n     ,
          const double       value ) ; 
        // resize & mullify the vector 
        Vector& resize
        ( const std::size_t n     ,
          const Zero    /* zero */ ) ;
        // ========================================================================
      public: // some simple math oprations 
        // ========================================================================
        /// add a vector of the same size 
        Vector& iadd        ( const Vector&  value ) ;
        /// add a constant 
        Vector& iadd        ( const double   value ) ;      
        /// add a vector of the same size 
        Vector& isub        ( const Vector&  value ) ;
        /// scale vector 
        Vector& imul        ( const double   value ) ;
        /// multiply by the matrix 
        Vector& imul        ( const Matrix&  value ) ;      
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
        /// cross/tensor product of two vectors 
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
        gsl_vector* m_vector { nullptr } ; //! the  actual pointer to GSL-vector
        // ========================================================================
      };
      // ==========================================================================
      /** @class Permutation
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
        /// get the element 
        inline std::size_t operator() ( const std::size_t i ) const { return get ( i ) ; }
        /// get the element 
        inline std::size_t operator[] ( const std::size_t i ) const { return get ( i ) ; }
        /// the size of the vector 
        inline std::size_t size () const { return m_permutation -> size ; } 
        // ========================================================================
      public: 
        // ========================================================================
        /// resize the permutation 
        Permutation& resize ( const std::size_t n ) ; 
        // ========================================================================
      public : 
        // ========================================================================
        /// swap two permutations 
        void swap ( Permutation& right ) ; // swap two permutations 
        // ========================================================================
      private :
        // ========================================================================
        /// the  actual pointer to GSL-permutation
        gsl_permutation* m_permutation { nullptr } ; //! the  actual pointer to GSL-permutation
        // ========================================================================
      };
      // ==========================================================================

      // ==========================================================================
      /// swap two matrices 
      inline void swap ( Matrix& a      , Matrix&      b ) { a.swap ( b ) ; } 
      /// swap two vectors 
      inline void swap ( Vector& a      , Vector&      b ) { a.swap ( b ) ; } 
      /// swap two permutations 
      inline void swap ( Permutation& a , Permutation& b ) { a.swap ( b ) ; } 
      // =========================================================================
      
      // =========================================================================
      /// equality of two matrices 
      inline bool operator==( const Matrix& a , const Matrix& b )
      { return &a == &b || ( a.nRows() == b.nRows() && a.nCols() == b.nCols() && a.equal ( b ) ) ; }
      /// equality of two vectors 
      inline bool operator==( const Vector& a , const Vector& b )
      { return &a == &b || ( a.size() == b.size() && a.equal ( b ) ) ; }
      
      /// non-equality of two matrices 
      inline bool operator!=( const Matrix& a , const Matrix& b ) { return !( a == b ) ; }
      /// non-equality of two vectors 
      inline bool operator!=( const Vector& a , const Vector& b ) { return !( a == b ) ; }
      
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
      
      // ========================================================================
      // "right"  forms 
      // ========================================================================
      
      /// scale the matrix from the right 
      inline Matrix operator*( const double  b , const Matrix& a )
      { Matrix c { a } ; c*= b ; return c ;; }    
      /// add     b*identity from rrght 
      inline Matrix operator+( const double& b , const Matrix& a )
      { Matrix c { a } ; c += b ; return c ; }
            
      /// add two vectors 
      inline Vector operator+( const Vector& a , const Vector& b )
      { Vector c { a } ; c += b ; return c ; }
      
      /// subtract two vectors 
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
      { Vector c { a } ; c*= b ; return c ; }
      
      /// scale the vector
      inline Vector operator/( const Vector& a , const double b  )
      { Vector c { a } ; c/= b ; return c ; }
      
      // ========================================================================
      // "right"  forms 
      // ========================================================================
      
      /// scale the vector  from the right 
      inline Vector operator*( const double  b , const Vector& a )
      { Vector c { a } ; c *= b ; return c ;; }    
      /// add the constant from right 
      inline Vector operator+( const double& b , const Vector& a )
      { Vector c { a } ; c += b ; return c ; }
      
      // ========================================================================
      // Matrix & permutations 
      // ========================================================================
      
      // ========================================================================
      /// Apply permutation 
      inline Matrix operator*
      ( const Permutation& p ,
        const Matrix&      a ) { return p.apply    ( a ) ; }
      // ========================================================================    
      /// apply permutation 
      inline Matrix operator*
      ( const Matrix&a       ,
        const Permutation& p ) { return a.multiply ( p ) ; }
      // ========================================================================
      
      // ========================================================================
      /** Matrix multiplication:
       *  \f$ C = a^{(T_a)} \times b^{(T_b)}\f$
       *  @param a (INPUT) the first input matrix
       *  @param Ta (INPUT) transpose the first matrix?
       *  @param b (INPUT) the second input matrix
       *  @param Tb (INPUT) transpose the second matrix?
       *  @param c (OUTPUT) the output matrix
       *  @return status code      
       */
      Ostap::StatusCode MM
      ( const Matrix& a  ,
        const bool    Ta ,
        const Matrix& b  ,
        const bool    Tb , 
        Matrix&       c  ) ;

      // =====================================================================
      /** Matrix multiplication of 3 matrices:
       *  \f$ D = a^{(T_a)} \times b^{(T_b)} \times c^{(T_c)}\f$
       *  @param a  (INPUT) the first input matrix
       *  @param Ta (INPUT) transpose the first matrix?
       *  @param b  (INPUT) the second input matrix
       *  @param Tb (INPUT) transpose the second matrix?
       *  @param c  (INPUT) the third input matrix
       *  @param Tc (INPUT) transpose the third matrix?
       *  @param d  (OUTPUT) the output matrix
       *  @return status code      
       */
      Ostap::StatusCode MMM
      ( const Matrix& a  ,
        const bool    Ta ,
        const Matrix& b  ,
        const bool    Tb , 
        const Matrix& c  ,
        const bool    Tc , 
        Matrix&       d  ) ;

      // =====================================================================
      /**  Matrix multiplication of 4 matrices:
       *  \f$ E = a^{(T_a)} \times b^{(T_b)} \times c^{(T_c)} \times d^{(T_d)}\f$
       *  @param a  (INPUT) the first input matrix
       *  @param Ta (INPUT) transpose the first matrix?
       *  @param b  (INPUT) the second input matrix
       *  @param Tb (INPUT) transpose the second matrix?
       *  @param c  (INPUT) the third input matrix
       *  @param Tc (INPUT) transpose the third matrix?
       *  @param d  (INPUT) the fourth input matrix
       *  @param Td (INPUT) transpose the fourth matrix?
       *  @param e  (OUTPUT) the output matrix
       *  @return status code      
       */
      Ostap::StatusCode MMMM
      ( const Matrix& a  ,
        const bool    Ta ,
        const Matrix& b  ,
        const bool    Tb , 
        const Matrix& c  ,
        const bool    Tc , 
        const Matrix& d  ,
        const bool    Td , 
        Matrix&       e  ) ;

      // =====================================================================
      /**Matrix multiplication with a diagonal matrix in the middle:
       *  \f$ C = a^{(T_a)} \times \text{diag}(d) \times b^{(T_b)}\f$
       *
       *  @param a  (INPUT) the first input matrix
       *  @param Ta (INPUT) transpose the first matrix?
       *  @param d  (INPUT) diagonal elements of the middle matrix
       *  @param b  (INPUT) the second input matrix
       *  @param Tb (INPUT) transpose the second matrix?
       *  @param c  (OUTPUT) the output matrix
       *  @return status code      
       */
      Ostap::StatusCode MDM
      ( const Matrix& a  ,
        const bool    Ta ,
        const Vector& d  ,
        const Matrix& b  ,
        const bool    Tb ,
        Matrix&       c  ) ; 

      // ======================================================================
      // DD operations (Diagonal x Diagonal)
      // ======================================================================

      /** Element-wise multiplication of two diagonal matrices: d = a * b
       *  @param a [in]  first diagonal vector
       *  @param b [in]  second diagonal vector
       *  @param d [out] resulting diagonal vector
       *  @return status code (Ostap::StatusCode::SUCCESS on success)
       */
      Ostap::StatusCode DD
      ( const Vector& a , 
        const Vector& b , 
        Vector&       d ) ;

      // ======================================================================
      // DMD operations (Diagonal x Matrix x Diagonal)
      // ======================================================================

      /** Symmetric scaling: r = d * m * d
       *  @param d [in]  diagonal vector
       *  @param m [in]  input matrix
       *  @param r [out] resulting matrix
       *  @return status code
       */
      Ostap::StatusCode DMD 
      ( const Vector& d , 
        const Matrix& m , 
        Matrix&       r ) ;

      /** Asymmetric scaling: r = a * m * b
       *  @param a [in]  left diagonal vector (scales rows)
       *  @param m [in]  input matrix
       *  @param b [in]  right diagonal vector (scales columns)
       *  @param r [out] resulting matrix
       *  @return status code
       */
      Ostap::StatusCode DMD
      ( const Vector& a , 
        const Matrix& m , 
        const Vector& b , 
        Matrix&       r ) ;

      /** Apply permutation transformation to a square matrix: R = P * M * P^T
       *  
       *  @param p [in]  Permutation matrix P
       *  @param m [in]  Square matrix M (N x N)
       *  @param r [out] Resulting permuted matrix (N x N)
       *  @return Status code (Ostap::StatusCode::SUCCESS on success)
       *
       *  @note Safe against argument aliasing (e.g., PMPt(p, m, m))
       */
      Ostap::StatusCode PMP 
      ( const Permutation& P ,
        const Matrix&      m ,
        Matrix&            r ) ;

      /** General asymmetric permutation transformation: R = Pl * M * Pr^T
       *  
       *  @param pl [in]  Left permutation matrix Pl (size matching M.k1())
       *  @param m  [in]  Input matrix M (M x N)
       *  @param pr [in]  Right permutation matrix Pr (size matching M.k2())
       *  @param r  [out] Resulting permuted matrix (M x N)
       *  @return Status code (Ostap::StatusCode::SUCCESS on success)
       *
       *  @note Safe against argument aliasing (e.g., PLMPrT(pl, m, pr, m))
       */
      Ostap::StatusCode PMP 
      ( const Permutation& PL ,
        const Matrix&      m  ,
        const Permutation& PR ,
        Matrix&            r  ) ;

      // ============================================================================
      /** General asymmetric combination: 
       *  \f$ R = P_l \cdot D_1 \cdot M \cdot D_2 \cdot P_r^T \f$
       *  
       *  @param pl [in]  Left permutation matrix Pl
       *  @param d1 [in]  Left diagonal vector D1
       *  @param m  [in]  Input matrix A (M x N)
       *  @param d2 [in]  Right diagonal vector D2
       *  @param pr [in]  Right permutation matrix Pr
       *  @param r  [out] Resulting matrix R (M x N)
       *  @return Status code
       */
      Ostap::StatusCode PDM 
      ( const Permutation& pl ,
        const Vector&      d1 ,
        const Matrix&      m  ,
        const Vector&      d2 ,
        const Permutation& pr ,
        Matrix&            r  ) ;

      /** Symmetric combination of permutation and diagonal scaling: 
       *  \f$ R = P \cdot D \cdot M \cdot D \cdot P^T \f$
       *  
       *  @param p [in]  Permutation matrix P
       *  @param d [in]  Diagonal vector D
       *  @param m [in]  Input square matrix A (N x N)
       *  @param r [out] Resulting matrix R (N x N)
       *  @return Status code (Ostap::StatusCode::SUCCESS on success)
       *
       *  @note Safe against argument aliasing (e.g., PDADPt(p, d, a, a))
       */
      Ostap::StatusCode PDM  
      ( const Permutation& p  ,
        const Vector&      d  ,
        const Matrix&      m  ,
        Matrix&            r  ) ;
      
      // ============================================================================
      /** Compute Moore-Penrose Pseudoinverse using SVD: A^+ = V * Sigma^+ * U^T
       *  @param a     (INPUT)  Input matrix A (m x n)
       *  @param a_pinv(OUTPUT) Pseudoinverse matrix A^+ (n x m)
       *  @param tol   (INPUT)  Tolerance for zeroing small singular values (< 0 for default)
       *  @return status code
       */
      Ostap::StatusCode PINV
      ( const Matrix& a      ,
        Matrix&       a_pinv ,
        double        tol    = - 1 ) ;

      // ========================================================================
      // Linear Algebra 
      // ========================================================================
      
      // ========================================================================
      // Eigen values and eigenvectors 
      // ========================================================================
      
      // ========================================================================
      // LU decomposition with pivoting 
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
       *  @param  A (update) input/update MxN  marix 
       *  @param  P (UPDATE/OUTPUT) permutation 
       *  @return status code 
       *  @see gsl_linalg_LU_decomp 
       */
      Ostap::StatusCode PLU
      ( Matrix&      A  , 
        Permutation& P ) ;        
      // ========================================================================
      /** perfom LU decomposition  
       *  @param  A   (INOUT)         input matrix 
       *  @param  P   (UPDATE/OUTPUT) permutation 
       *  @param  LU  (UPDATE/OUTPUT) output LU matrix 
       *  @return status code 
       *  @see gsl_linalg_LU_decomp 
       */
      Ostap::StatusCode PLU
      ( const Matrix& A  ,
        Permutation&  P  , 
        Matrix&       LU ) ;
      // ========================================================================
      /** perfom LU decomposition  
       *  @param  A   (INOUT)         input matrix 
       *  @param  P   (UPDATE/OUTPUT) permutation 
       *  @param  L   (UPDATE/OUTPUT) lower triangular matrix 
       *  @param  U   (UPDATE/OUTPUT) upper triangular matrix 
       *  @return status code
       *  @see gsl_linalg_LU_decomp 
       */
      Ostap::StatusCode PLU
      ( const Matrix& A ,
        Permutation&  P , 
        Matrix&       L ,
        Matrix&       U ) ;    
      // ========================================================================
      
      // ========================================================================
      // QR decomposition with column pivoting 
      // ========================================================================
      
      // ========================================================================    
      /** make QR Decomposion of matrix A : \f$ AP = QR\f$ where 
       *  - A is input                 MxN matrix  
       *  - P is permutation           N  
       *  - Q is orthogonal matrix     MxM 
       *  - R is right triaular matrix MxN 
       *  
       *  @param A  (input) the matrix to decopose 
       *  @param P  (outpt/update) permutation matrix P
       *  @param Q  (outpt/update) orthogonal matrix Q 
       *  @param R  (outpt/update) rigth triangular matrix R 
       *  @return permutation P 
       */
      Ostap::StatusCode PQR
      ( const Matrix& A ,
        Permutation&  P , 
        Matrix&       Q ,
        Matrix&       R ) ;   
    // ======================================================================== 
      /** make QR Decomposion of matrix A : \f$ AP = QR\f$ where 
       *  - A is input                 MxN matrix  
       *  - P is permutation           N 
       *  - Q is orthogonal matrix     MxM 
       *  - R is right triaular matrix MxN
       *  - r is a condition number of the matrix R 
       *  
       *  @param A  (input) the matrix to decopose 
       *  @param P  (outpt/update) permutation matrix P
       *  @param Q  (outpt/update) orthogonal matrix Q 
       *  @param R  (outpt/update) rigth triangular matrix R 
       *  @param r  (outpt/update) condition number of the matrix R
       *  @return permutation P 
       */
      Ostap::StatusCode PQR
      ( const Matrix& A ,
        Permutation&  P , 
        Matrix&       Q ,
        Matrix&       R , 
        double&       r ) ;    
      // ======================================================================
      
      // ======================================================================
      // LQ & QL decompositions
      // ======================================================================
      
      // ======================================================================
      /** LQ decomposition of matrix A: \f$ A = LQ\f$, where 
       *  - L is lower trapezoidal MxN 
       *  - Q is orthogonal NxN 
       */ 
      Ostap::StatusCode LQ
      ( const Matrix& A ,
        Matrix&       L ,
        Matrix&       Q ) ;
      // ======================================================================
      /** QL decomposition of matrix A: \f$ A = QL\f$, where 
       *  - Q is orthogonal MxM
       *  - L is lower trapezoidal MxN 
       */ 
      Ostap::StatusCode QL
      ( const Matrix& A ,
        Matrix&       Q ,
        Matrix&       L ) ;
      // ======================================================================
      
      // ======================================================================
      // COD decomposition
      // ======================================================================
      
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
      Ostap::StatusCode COD
      ( const Matrix& A ,
        Permutation & P , 
        Matrix&       Q ,
        Matrix&       R ,
        Matrix&       Z ) ;
      // ======================================================================
      
      // ======================================================================
      // Singular values and SVD decomposition
      // ======================================================================
      
      // ======================================================================
      /** SVD : singular Value Decomposition  \f$ A = U S V^T\f$
       *  - A input MxN matrix 
       *  - K = min ( M , N ) : 
       *  - U MxK orthogonal matrix 
       *  - S KxK Diagonal matrix of singular values 
       *  - V NxK orthogonal matrix 
       *  @param A     (input)  input matrix A 
       *  @param S     (update) singular values U 
       *  @param U     (update) orthogonal matrix U 
       *  @param V     (update) orthogonal matrix V 
       *  @param golub (input) use Golub or Jacobi algorithm 
       *  @return vector of singular values 
       * -  Jacobi' algorithm is more precise and Golub' algorithm is more CPU efficient 
       */
      Ostap::StatusCode SVD
      ( const Matrix& A            ,
        Vector&       S            , 
        Matrix&       U            ,
        Matrix&       V            ,
        const bool    golub = true ) ;
      // ======================================================================
      
      // ======================================================================
      /** LLT : Cholesky decomposition of positive definite matrix \f$ A = L L^T\f$, 
       *  Only lower triangular part of the matrix A is used.
       *  @param A (input)  input MxM matrix
       *  @param L (update) lower triangular matrix
       *  @return
       */  
      Ostap::StatusCode LLT
      ( const Matrix& A ,
        Matrix&       L ) ;

      // ======================================================================
      /** LDLT : Cholesky decomposition of positive definite matrix 
       * \f$ PSASP^T = L D L^T\f$, 
       *  Only lower triangular part of the matrix A is used.
       *  @param A (input)  input MxM matrix
       *  @param S (output/update) scale vector/diagonal matrix 
       *  @param P (output/update) permutation 
       *  @param L (utput/update)  lower triangular matrix
       *  @param D (output/update) vector/diagonal matrix   
       *  @return status code
       */  
      Ostap::StatusCode LDLT
      ( const Matrix& A ,
        Vector&       S ,
        Permutation&  P , 
        Matrix&       L , 
        Vector&       D ) ;

      // =======================================================================
      /** D3 : decomposition of symmetric matrix \f$ A = Q D_3 Q^T \f$, where
       *  - \f$ Q \f$ is orthogonal matrix
       *  - \f$ D_2\fF is symmetric  trigiagonal matrix 
       *  @param A (INPUT) input matrix A
       *  @param Q (OUTPUT/UPDATE) orthogonal matrix Q
       *  @param d (OUTPUT/UPDATE) main diagonal of symmetric  matrix \f$ D_3 \f$
       *  @param s (OUTPUT/UPDATE) sub-diagonal of symmetric  matrix \f$ D_3 \f$
       *  @return status code 
       */
      Ostap::StatusCode D3
      ( const Matrix& A ,
        Matrix&       Q ,
        Vector&       d ,
        Vector&       s ) ;
             
      // =======================================================================
      /** D3 : decomposition of symmetric matrix \f$ A = Q D_3 Q^T \f$, where
       *  - \f$ Q \f$ is orthogonal matrix
       *  - \f$ D_2\fF is symmetric  trigiagonal matrix 
       *  @param A (INPUT) input matrix A
       *  @param Q (OUTPUT/UPDATE) orthogonal matrix Q
       *  @param D (OUTPUT/UPDATE) symmetric tridiagonal matrix 
       *  @return status code 
       */
      Ostap::StatusCode D3
      ( const Matrix& A ,
        Matrix&       Q ,
        Matrix&       D ) ;

      // =======================================================================
      /** Hessenberg decomposition of square matrix \f$ A = U H Q^T \f$, where
       *  - \f$ U \f$ is orthogonal 
       *  - \f$ H \f$ is Hessenberg' matrix: \f$ H(i,i)=0 \f$ for \f$ i > j + 1 \f$
       *  @param A (INPUT) input matrix A
       *  @param Q (OUTPUT/UPDATE) orthogonal matrix Q
       *  @param H (OUTPUT/UPDATE) Hessenberg matrix 
       *  @return status code 
       */
      Ostap::StatusCode UHUT
      ( const Matrix& A ,
        Matrix&       Q ,
        Matrix&       H ) ;

      // =======================================================================
      /** Bidiagonalization of of general matrix \f$ A = U B V^T \f$, where
       *  - \f$ A \f$ is \f$ M \times N \f$ matrix
       *  - \f$ U \f$ is \f$ M\times N \f$ orthogonal matrix 
       *  - \f$ B \f$ is \f$ N\times N\f$  square biadiagonal matrix : \f$ B_{i,j} = 0\f$ if \f$ j \ne i,i+1\f$
       *  - \f$ V \f$ is \f$ N\times N \f$ orthogonal matrix 
       *  @param A (INPUT) input matrix A
       *  @param U (OUTPUT/UPDATE) orthogonal matrix U
       *  @param B (OUTPUT/UPDATE) bidiagonal matrix B
       *  @param V (OUTPUT/UPDATE) orthogonal matrix V
       *  @return status code        
       */
      Ostap::StatusCode UBVT
      ( const Matrix& A ,
        Matrix&       U ,
        Matrix&       B , 
        Matrix&       V ) ;
      
      // =======================================================================
      /** Bidiagonalization of of general matrix \f$ A = U B V^T \f$, where
       *  - \f$ A \f$ is \f$ M \times N \f$ matrix
       *  - \f$ U \f$ is \f$ M\times N \f$ orthogonal matrix 
       *  - \f$ B \f$ is \f$ N\times N\f$  square biadiagonal matrix : \f$ B_{i,j} = 0\f$ if \f$ j \ne i,i+1\f$
       *  - \f$ V \f$ is \f$ N\times N \f$ orthogonal matrix 
       *  @param A (INPUT) input matrix A
       *  @param U (OUTPUT/UPDATE) orthogonal matrix U
       *  @param d (OUTPUT/UPDATE) diagonal 
       *  @param s (OUTPUT/UPDATE) super-diagonal  
       *  @param V (OUTPUT/UPDATE) orthogonal matrix V
       *  @return status code        
       */
      Ostap::StatusCode UBVT
      ( const Matrix& A ,
        Matrix&       U ,
        Vector&       d ,
        Vector&       s ,        
        Matrix&       V ) ;      
      
      // ======================================================================
      // Schur' decomposition of square matrix
      // ======================================================================
      
      // ======================================================================
      /** Schur's decomposition of square matrix \f$ A = Z T Z^T\f$, where 
       *  - A is input MxM (square) matrix
       *  - S is Schur' form of matrix  
       *  - Z is orthogonal matrix 
       */
      Ostap::StatusCode SCHUR 
      ( const Matrix&  A ,  
        Matrix&        Z , 
        Matrix&        S ) ; 
      // ======================================================================
      
      // ======================================================================
      // Polar decomposition of square matrix 
      // ======================================================================
      
      // ======================================================================
      /** Polar decompositon of the square matrix A: \f$ A = UP \f$
       *  - U is orthogonal 
       *  - P is positive semi-definitive 
       */
      Ostap::StatusCode POLAR
      ( const Matrix& A ,
        Matrix      & U ,
        Matrix      & P ) ;
      // ======================================================================    
    } //                                  The end of namespace Ostap::Math::GSL
    // ========================================================================
    
    // ========================================================================
    /// get the max element
    double max_element    ( const Ostap::Math::GSL::Matrix& m ) ;
    /// get the min element
    double min_element    ( const Ostap::Math::GSL::Matrix& m ) ;
    /// get the max element
    double max_element    ( const Ostap::Math::GSL::Vector& v ) ;
    /// get the min element
    double min_element    ( const Ostap::Math::GSL::Vector& v ) ;    
    /// get the element with maximal absolute value 
    double maxabs_element ( const Ostap::Math::GSL::Matrix& m ) ;
    /// get the element with maximal absolute value 
    double maxabs_element ( const Ostap::Math::GSL::Vector& v ) ;
    // ========================================================================
    
    // ========================================================================
    /// Is this vector finite    ?
    inline bool isfinite  ( const Ostap::Math::GSL::Vector& v ) { return v.isfinite () ; }
    /// Is this matrix finite    ?
    inline bool isfinite  ( const Ostap::Math::GSL::Matrix& m ) { return m.isfinite () ; }
    /// Is this matrix symmetric ?
    bool        symmetric ( const Ostap::Math::GSL::Matrix& m ) ;
    // =======================================================================
    
    // =======================================================================
    /** Can this matrix be symmetric & positive-definite ?
     *  - Finite 
     *  - Diagonal elements are finite and positive
     *  - Symmetric
     *  - have CholeskyDecomposition
     */
    bool symmetric_positive_definite 
    ( const Ostap::Math::GSL::Matrix& mtrx ) ;

    // ========================================================================
    /** Can this matrix be a covariance matrix?
     *  - Square
     *  - Finite 
     *  - Symmetric 
     *  - Diagonal elements are finite and positive
     *  - Off-diagonal elements are finite and not too large 
     *  - have CholeskyDecomposition
     */
    bool covariance_matrix 
    ( const Ostap::Math::GSL::Matrix& mtrx ) ;

    // =======================================================================
    /// numerical equality of two matrices 
    template <> 
    struct Equal_To<Ostap::Math::GSL::Matrix>
    {
      inline bool operator()
      ( const Ostap::Math::GSL::Matrix& m1 , 
        const Ostap::Math::GSL::Matrix& m2 ) const
      { return m1.nRows() == m2.nRows() && m1.nCols() == m2.nCols() && m1.equal ( m2 ) ; }
    } ;
    
    // =======================================================================
    /// numerical equality of two vectors 
    template <> 
    struct Equal_To<Ostap::Math::GSL::Vector>
    {
      inline bool operator()
      ( const Ostap::Math::GSL::Vector& v1 , 
        const Ostap::Math::GSL::Vector& v2 ) const
      { return v1.size() == v2.size() && v1.equal ( v2 ) ; }
    } ;
    
    // =======================================================================
    /// matrix is numerical zero 
    template <> 
    struct Zero<Ostap::Math::GSL::Matrix>
    {
      inline bool operator()
      ( const Ostap::Math::GSL::Matrix& m ) const
      { return m.iszero () ; }
    } ;

    // =======================================================================
    /// vector is numerical zero 
    template <> 
    struct Zero<Ostap::Math::GSL::Vector>
    {
      inline bool operator()
      ( const Ostap::Math::GSL::Vector& v ) const
      { return v.iszero () ; }
    } ;

    // =======================================================================
  } //                                        The end of namespace Ostap::Math 
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
    ( const Ostap::Math::GSL::Matrix& m ,
      std::ostream&                   s ) ;
    // ========================================================================
    /** print GSL-vector to the stream 
     *  @param v the vector 
     *  @param s the stream 
     *  @return the stream 
     */
    std::ostream& toStream
    ( const Ostap::Math::GSL::Vector& v ,
      std::ostream&                   s ) ;
    // ========================================================================
    /** print GSL-permutation to the stream 
     *  @param p the permutation
     *  @param s the stream 
     *  @return the stream 
     */
    std::ostream& toStream
    ( const Ostap::Math::GSL::Permutation& p ,
      std::ostream&                        s ) ;
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
( std::ostream&                   s ,
  const Ostap::Math::GSL::Matrix& m ) 
{ return Ostap::Utils::toStream ( m , s ) ; }
// ============================================================================
/// print operator 
inline std::ostream& operator<<
( std::ostream&                   s ,
  const Ostap::Math::GSL::Vector& v ) 
{ return Ostap::Utils::toStream ( v , s ) ; }
// ============================================================================
/// print operator 
inline std::ostream& operator<<
( std::ostream&                        s ,
  const Ostap::Math::GSL::Permutation& p ) 
{ return Ostap::Utils::toStream ( p , s ) ; }
// ============================================================================
#endif // OSTAP_GSL_LINALG_H
// ============================================================================
//                                                                      The END 
// ============================================================================
