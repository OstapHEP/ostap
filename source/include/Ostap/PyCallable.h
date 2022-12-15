// ============================================================================
#ifndef OSTAP_PYCALLABLE_H 
#define OSTAP_PYCALLABLE_H 1
// ============================================================================
// Include files
// ============================================================================
// Python
// ============================================================================
struct  _object ;
typedef _object PyObject ;
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Functions 
  {
    // ========================================================================
    /** @class PyCallable Ostap/PyCallable.h
     *  Simple C++   wrapper for the python callable 
     *  @author Vanya Belyaev
     *  @date   2019-09-25
     */
    class PyCallable final 
    {
    public  : 
      // ======================================================================
      /// constructor from the callable object 
      PyCallable ( PyObject* callable , const bool ok ) ;
      /// copy constructor 
      PyCallable ( const PyCallable&  right ) ;
      /// Move constructor 
      PyCallable (       PyCallable&& right ) ;
      /// Destructor 
      ~PyCallable () ;
      // ======================================================================
    public:
      // ======================================================================
      /// the main method:
      double operator () ( const double x ) const { return evaluate ( x ) ; }
      /// the main method 
      double evaluate    ( const double x ) const ;
      // ======================================================================
    private :
      // ======================================================================
      // the callable object 
      PyObject* m_callable   { nullptr } ;
      // the tuple of arguments 
      PyObject* m_arguments  { nullptr } ;
      // ======================================================================
    };
    // ========================================================================
    /** @class PyCallable2 Ostap/PyCallable.h
     *  Simple C++   wrapper for the python callable 
     *  @author Vanya Belyaev
     *  @date   2019-09-25
     */
    class PyCallable2 final 
    {
    public  : 
      // ======================================================================
      /// constructor from the callable object 
      PyCallable2 ( PyObject* callable , const bool ok ) ;
      /// copy constructor 
      PyCallable2 ( const PyCallable2&  right ) ;
      /// Move constructor 
      PyCallable2 (       PyCallable2&& right ) ;
      /// Destructor 
      ~PyCallable2() ;
      // ======================================================================
    public:
      // ======================================================================
      /// the main method:
      double operator () 
        ( const double x , 
          const double y ) const { return evaluate ( x , y ) ; }
      /// the main method 
      double evaluate 
        ( const double x , 
          const double y ) const ;
      // ======================================================================
    private :
      // ======================================================================
      // the callable object 
      PyObject* m_callable   { nullptr } ;
      // the tuple of arguments 
      PyObject* m_arguments  { nullptr } ;
      // ======================================================================
    };
    // ========================================================================
    /** @class PyCallable3 Ostap/PyCallable.h
     *  Simple C++   wrapper for the python callable 
     *  @author Vanya Belyaev
     *  @date   2019-09-25
     */
    class PyCallable3 final 
    {
    public  : 
      // ======================================================================
      /// constructor from the callable object 
      PyCallable3 ( PyObject* callable , const bool ok ) ;
      /// copy constructor 
      PyCallable3 ( const PyCallable3&  right ) ;
      /// Move constructor 
      PyCallable3 (       PyCallable3&& right ) ;
      /// Destructor 
      ~PyCallable3() ;
      // ======================================================================
    public:
      // ======================================================================
      /// the main method:
      double operator () 
        ( const double x , 
          const double y , 
          const double z ) const { return evaluate ( x , y , z ) ; }
      /// the main method 
      double evaluate 
        ( const double x , 
          const double y ,
          const double z ) const ;
      // ======================================================================
    private :
      // ======================================================================
      // the callable object 
      PyObject* m_callable   { nullptr } ;
      // the tuple of arguments 
      PyObject* m_arguments  { nullptr } ;
      // ======================================================================
    };
    // ========================================================================
  } //                                    The end of namespace Ostap::Functions 
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_PYCALLABLE_H
// ============================================================================
