// ============================================================================
#ifndef OSTAP_PARAMS_H 
#define OSTAP_PARAMS_H 1
// ============================================================================
// Include files
// ============================================================================
#include <limits>
#include <string>
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/StatusCode.h"
// ============================================================================
// Forward desclarations 
// ============================================================================ 
// ROOT 
// ============================================================================ 
class TTree ; //  from ROOT 
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    class LegendreSum  ;
    class LegendreSum2 ;
    class LegendreSum3 ;
    class LegendreSum4 ;
    class ChebyshevSum ;
    class Bernstein    ;
    class Bernstein2D  ;
    class Bernstein3D  ;
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
  /** @class  DataParam  Ostap/Params.h
   *  Helper  class to parameterize the data 
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date   2019-07-3
   */
  class DataParam 
  {
  public:
    // =========================================================================
    // 1D LegendreSum
    // =========================================================================
    /** fill Legendre sum with data from the Tree 
     *  @see Ostap::Math::LegendreSum 
     *  @see Ostap::Math::LegendreSum::fill
     *  @param tree       (INPUT)  the input tree 
     *  @param sum        (UPDATE) the parameterization object 
     *  @param expression (INPUT)  expression to be parameterized
     *  @param first      (INPUT)  the first event in Tree 
     *  @param last       (INPUT)  the last  event in Tree 
     *  @code
     *  Tree*  tree = ...
     *  LegendreSum s ( 5 , -1 , 1 ) ;
     *  DataParam::parameterize ( tree , s , "x" ) ;
     *  @endcode
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2019-07-3
     */ 
    static Ostap::StatusCode parameterize 
    ( TTree*                    tree       , 
      Ostap::Math::LegendreSum& sum        , 
      const std::string&        expression , 
      const unsigned long       first      =  0 ,
      const unsigned long       last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** fill Legendre sum with data from the Tree 
     *  @see Ostap::Math::LegendreSum 
     *  @see Ostap::Math::LegendreSum::fill
     *  @param tree       (INPUT)  the input tree 
     *  @param sum        (UPDATE) the parameterization object 
     *  @param expression (INPUT)  expression to be parameterized
     *  @param seelction  (INPUT)  selection/weight to be used 
     *  @param first      (INPUT)  the first event in Tree 
     *  @param last       (INPUT)  the last  event in Tree 
     *  @code
     *  Tree*  tree = ...
     *  LegendreSum s ( 5 , -1 , 1 ) ;
     *  DataParam::parameterize ( tree , s , "x"  , "y>10" ) ;
     *  @endcode
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2019-07-3
     */ 
    static Ostap::StatusCode parameterize 
    ( TTree*                    tree       , 
      Ostap::Math::LegendreSum& sum        , 
      const std::string&        expression , 
      const std::string&        selection  , 
      const unsigned long       first      =  0 ,
      const unsigned long       last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
  public: // Parameterse in terms of Chebyshed polynomials 
    // =========================================================================
    /** fill Chebyshev sum with data from the Tree 
     *  @see Ostap::Math::ChebyshevSum 
     *  @see Ostap::Math::ChebyshevSum::fill
     *  @param tree       (INPUT)  the input tree 
     *  @param sum        (UPDATE) the parameterization object 
     *  @param expression (INPUT)  expression to be parameterized
     *  @param first      (INPUT)  the first event in Tree 
     *  @param last       (INPUT)  the last  event in Tree 
     *  @code
     *  Tree*  tree = ...
     *  ChebyshevSum s ( 5 , -1 , 1 ) ;
     *  DataParam::parameterize ( tree , s , "x" ) ;
     *  @endcode
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2019-07-3
     */ 
    static Ostap::StatusCode parameterize 
    ( TTree*                     tree       , 
      Ostap::Math::ChebyshevSum& sum        , 
      const std::string&         expression , 
      const unsigned long        first      =  0 ,
      const unsigned long        last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** fill Chebyshev sum with data from the Tree 
     *  @see Ostap::Math::ChebyshevSum 
     *  @see Ostap::Math::ChebyshevSum::fill
     *  @param tree       (INPUT)  the input tree 
     *  @param sum        (UPDATE) the parameterization object 
     *  @param expression (INPUT)  expression to be parameterized
     *  @param seelction  (INPUT)  selection/weight to be used 
     *  @param first      (INPUT)  the first event in Tree 
     *  @param last       (INPUT)  the last  event in Tree 
     *  @code
     *  Tree*  tree = ...
     *  ChebyshevSum s ( 5 , -1 , 1 ) ;
     *  DataParam::parameterize ( tree , s , "x"  , "y>10" ) ;
     *  @endcode
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2019-07-3
     */ 
    static Ostap::StatusCode parameterize 
    ( TTree*                     tree       , 
      Ostap::Math::ChebyshevSum& sum        , 
      const std::string&         expression , 
      const std::string&         selection  , 
      const unsigned long        first      =  0 ,
      const unsigned long        last       = std::numeric_limits<unsigned long>::max() ) ;
    // =======================================================================
  public: // Parameterse in terms of Bernstein polynomials 
    // =========================================================================
    /** fill Bernstein sum with data from the Tree 
     *  @see Ostap::Math::Bernstein 
     *  @see Ostap::Math::Bernstein::fill
     *  @param tree       (INPUT)  the input tree 
     *  @param sum        (UPDATE) the parameterization object 
     *  @param expression (INPUT)  expression to be parameterized
     *  @param first      (INPUT)  the first event in Tree 
     *  @param last       (INPUT)  the last  event in Tree 
     *  @code
     *  Tree*  tree = ...
     *  Bernstein s ( 5 , -1 , 1 ) ;
     *  DataParam::parameterize ( tree , s , "x" ) ;
     *  @endcode
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2019-07-3
     */ 
    static Ostap::StatusCode parameterize 
    ( TTree*                     tree       , 
      Ostap::Math::Bernstein&    sum        , 
      const std::string&         expression , 
      const unsigned long        first      =  0 ,
      const unsigned long        last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** fill Bernstein sum with data from the Tree 
     *  @see Ostap::Math::Bernstein 
     *  @see Ostap::Math::Bernstein::fill
     *  @param tree       (INPUT)  the input tree 
     *  @param sum        (UPDATE) the parameterization object 
     *  @param expression (INPUT)  expression to be parameterized
     *  @param seelction  (INPUT)  selection/weight to be used 
     *  @param first      (INPUT)  the first event in Tree 
     *  @param last       (INPUT)  the last  event in Tree 
     *  @code
     *  Tree*  tree = ...
     *  Bernstein s ( 5 , -1 , 1 ) ;
     *  DataParam::parameterize ( tree , s , "x"  , "y>10" ) ;
     *  @endcode
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2019-07-3
     */ 
    static Ostap::StatusCode parameterize 
    ( TTree*                     tree       , 
      Ostap::Math::Bernstein&    sum        , 
      const std::string&         expression , 
      const std::string&         selection  , 
      const unsigned long        first      =  0 ,
      const unsigned long        last       = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
  public:
    // =========================================================================
    /** fill Legendre sum with data from the Tree 
     *  @see Ostap::Math::LegendreSum2 
     *  @see Ostap::Math::LegendreSum2::fill
     *  @param tree        (INPUT)  the input tree 
     *  @param sum         (UPDATE) the parameterization object 
     *  @param xexpression (INPUT)  x-expression to be parameterized
     *  @param yexpression (INPUT)  y-expression to be parameterized
     *  @param first       (INPUT)  the first event in Tree 
     *  @param last        (INPUT)  the last  event in Tree 
     *  @code
     *  Tree*  tree = ...
     *  LegendreSum2 s ( 5 , 3 , -1 , 1 , -2 , 2 ) ;
     *  DataParam::parameterize ( tree , s , "x" , "y/z") ;
     *  @endcode
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2019-07-3
     */ 
    static Ostap::StatusCode parameterize 
    ( TTree*                     tree        , 
      Ostap::Math::LegendreSum2& sum         , 
      const std::string&         xexpression , 
      const std::string&         yexpression , 
      const unsigned long        first       =  0 ,
      const unsigned long        last        = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** fill Legendre sum with data from the Tree 
     *  @see Ostap::Math::LegendreSum2 
     *  @see Ostap::Math::LegendreSum2::fill
     *  @param tree        (INPUT)  the input tree 
     *  @param sum         (UPDATE) the parameterization object 
     *  @param xexpression (INPUT)  x-expression to be parameterized
     *  @param yexpression (INPUT)  y-expression to be parameterized
     *  @param selection   (INPUT)  selection/weight to be used 
     *  @param first       (INPUT)  the first event in Tree 
     *  @param last        (INPUT)  the last  event in Tree 
     *  @code
     *  Tree*  tree = ...
     *  LegendreSum2 s ( 5 , 2 ,  -1 , 1 , -4  , -5 ) ;
     *  DataParam::parameterize ( tree , s , "x"  , "z" , "y>10" ) ;
     *  @endcode
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2019-07-3
     */ 
    static Ostap::StatusCode parameterize 
    ( TTree*                     tree        ,  
      Ostap::Math::LegendreSum2& sum         , 
      const std::string&         xexpression , 
      const std::string&         yexpression , 
      const std::string&         selection   , 
      const unsigned long        first       =  0 ,
      const unsigned long        last        = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** fill BErnstein with data from the Tree 
     *  @see Ostap::Math::Bernstein2D
     *  @see Ostap::Math::Bernstein2D::fill
     *  @param tree        (INPUT)  the input tree 
     *  @param sum         (UPDATE) the parameterization object 
     *  @param xexpression (INPUT)  x-expression to be parameterized
     *  @param yexpression (INPUT)  y-expression to be parameterized
     *  @param first       (INPUT)  the first event in Tree 
     *  @param last        (INPUT)  the last  event in Tree 
     *  @code
     *  Tree*  tree = ...
     *  Bernsteinn2D s ( 5 , 3 , -1 , 1 , -2 , 2 ) ;
     *  DataParam::parameterize ( tree , s , "x" , "y/z") ;
     *  @endcode
     *  @attention it is less CPU efficient than Ostap::Math::LegendreSum2 
     *  @see Ostap::Math::LegendreSum2 
     *  @see Ostap::Math::LegendreSum2::fill
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2022-12-26
     */ 
    static Ostap::StatusCode parameterize 
    ( TTree*                     tree        , 
      Ostap::Math::Bernstein2D&  sum         , 
      const std::string&         xexpression , 
      const std::string&         yexpression , 
      const unsigned long        first       =  0 ,
      const unsigned long        last        = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** fill Bernstein with data from the Tree 
     *  @see Ostap::Math::Bernstein2D
     *  @see Ostap::Math::Bernstein2D::fill
     *  @param tree        (INPUT)  the input tree 
     *  @param sum         (UPDATE) the parameterization object 
     *  @param xexpression (INPUT)  x-expression to be parameterized
     *  @param yexpression (INPUT)  y-expression to be parameterized
     *  @param selection   (INPUT)  selection/weight to be used 
     *  @param first       (INPUT)  the first event in Tree 
     *  @param last        (INPUT)  the last  event in Tree 
     *  @code
     *  Tree*  tree = ...
     *  Bernstein2D s ( 5 , 2 ,  -1 , 1 , -4  , -5 ) ;
     *  DataParam::parameterize ( tree , s , "x"  , "z" , "y>10" ) ;
     *  @endcode
     *  @attention it is less CPU efficient than Ostap::Math::LegendreSum2 
     *  @see Ostap::Math::LegendreSum2 
     *  @see Ostap::Math::LegendreSum2::fill
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2022-12-26
     */ 
    static Ostap::StatusCode parameterize 
    ( TTree*                     tree        ,  
      Ostap::Math::Bernstein2D&  sum         , 
      const std::string&         xexpression , 
      const std::string&         yexpression , 
      const std::string&         selection   , 
      const unsigned long        first       =  0 ,
      const unsigned long        last        = std::numeric_limits<unsigned long>::max() ) ;
    // =======================================================================
  public:  // 3D 
    // ========================================================================
    /** fill Legendre sum with data from the Tree 
     *  @see Ostap::Math::LegendreSum3 
     *  @see Ostap::Math::LegendreSum3::fill
     *  @param tree        (INPUT)  the input tree 
     *  @param sum         (UPDATE) the parameterization object 
     *  @param xexpression (INPUT)  x-expression to be parameterized
     *  @param yexpression (INPUT)  y-expression to be parameterized
     *  @param zexpression (INPUT)  z-expression to be parameterized
     *  @param first       (INPUT)  the first event in Tree 
     *  @param last        (INPUT)  the last  event in Tree 
     *  @code
     *  Tree*  tree = ...
     *  LegendreSum3 s ( 5 , 3 , 2 , -1 , 1 , -2 , 2 , 0 , 4 ) ;
     *  DataParam::parameterize ( tree , s , "x" , "y" , "y/z") ;
     *  @endcode
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2019-07-3
     */ 
    static Ostap::StatusCode parameterize 
    ( TTree*                     tree        , 
      Ostap::Math::LegendreSum3& sum         , 
      const std::string&         xexpression , 
      const std::string&         yexpression , 
      const std::string&         zexpression , 
      const unsigned long        first       =  0 ,
      const unsigned long        last        = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** fill Legendre sum with data from the Tree 
     *  @see Ostap::Math::LegendreSum3 
     *  @see Ostap::Math::LegendreSum3::fill
     *  @param tree        (INPUT)  the input tree 
     *  @param sum         (UPDATE) the parameterization object 
     *  @param xexpression (INPUT)  x-expression to be parameterized
     *  @param yexpression (INPUT)  y-expression to be parameterized
     *  @param zexpression (INPUT)  z-expression to be parameterized
     *  @param selection   (INPUT)  selection/weight to be used 
     *  @param first       (INPUT)  the first event in Tree 
     *  @param last        (INPUT)  the last  event in Tree 
     *  @code
     *  Tree*  tree = ...
     *  LegendreSum3 s ( 5 , 2 , 4 ,  -1 , 1 , -4  , -5 , 0 , 3 ) ;
     *  DataParam::parameterize ( tree , s , "x"  , "y" , "z" , "t>10" ) ;
     *  @endcode
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2019-07-3
     */ 
    static Ostap::StatusCode parameterize 
    ( TTree*                     tree        ,  
      Ostap::Math::LegendreSum3& sum         , 
      const std::string&         xexpression , 
      const std::string&         yexpression , 
      const std::string&         zexpression , 
      const std::string&         selection   , 
      const unsigned long        first       =  0 ,
      const unsigned long        last        = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** fill Bernstein with data from the Tree 
     *  @see Ostap::Math::Bernstein3D
     *  @see Ostap::Math::Bernstein3D::fill
     *  @param tree        (INPUT)  the input tree 
     *  @param sum         (UPDATE) the parameterization object 
     *  @param xexpression (INPUT)  x-expression to be parameterized
     *  @param yexpression (INPUT)  y-expression to be parameterized
     *  @param zexpression (INPUT)  z-expression to be parameterized
     *  @param first       (INPUT)  the first event in Tree 
     *  @param last        (INPUT)  the last  event in Tree 
     *  @code
     *  Tree*  tree = ...
     *  Bernstein3D s ( 5 , 3 , 2 , -1 , 1 , -2 , 2 , 0 , 4 ) ;
     *  DataParam::parameterize ( tree , s , "x" , "y" , "y/z") ;
     *  @endcode
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2019-07-3
     */ 
    static Ostap::StatusCode parameterize 
    ( TTree*                     tree        , 
      Ostap::Math::Bernstein3D&  sum         , 
      const std::string&         xexpression , 
      const std::string&         yexpression , 
      const std::string&         zexpression , 
      const unsigned long        first       =  0 ,
      const unsigned long        last        = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** fill Bernstein with data from the Tree 
     *  @see Ostap::Math::Bernstein3D
     *  @see Ostap::Math::Bernstein3D::fill
     *  @param tree        (INPUT)  the input tree 
     *  @param sum         (UPDATE) the parameterization object 
     *  @param xexpression (INPUT)  x-expression to be parameterized
     *  @param yexpression (INPUT)  y-expression to be parameterized
     *  @param zexpression (INPUT)  z-expression to be parameterized
     *  @param selection   (INPUT)  selection/weight to be used 
     *  @param first       (INPUT)  the first event in Tree 
     *  @param last        (INPUT)  the last  event in Tree 
     *  @code
     *  Tree*  tree = ...
     *  Bernstein3D s ( 5 , 3 , 2 , -1 , 1 , -2 , 2 , 0 , 4 ) ;
     *  DataParam::parameterize ( tree , s , "x"  , "y" , "z" , "t>10" ) ;
     *  @endcode
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2019-07-3
     */ 
    static Ostap::StatusCode parameterize 
    ( TTree*                     tree        ,  
      Ostap::Math::Bernstein3D&  sum         , 
      const std::string&         xexpression , 
      const std::string&         yexpression , 
      const std::string&         zexpression , 
      const std::string&         selection   , 
      const unsigned long        first       =  0 ,
      const unsigned long        last        = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
  public : // 4D 
    // ========================================================================
    /** fill Legendre sum with data from the Tree 
     *  @see Ostap::Math::LegendreSum4 
     *  @see Ostap::Math::LegendreSum4::fill
     *  @param tree        (INPUT)  the input tree 
     *  @param sum         (UPDATE) the parameterization object 
     *  @param xexpression (INPUT)  x-expression to be parameterized
     *  @param yexpression (INPUT)  y-expression to be parameterized
     *  @param zexpression (INPUT)  z-expression to be parameterized
     *  @param uexpression (INPUT)  u-expression to be parameterized
     *  @param first       (INPUT)  the first event in Tree 
     *  @param last        (INPUT)  the last  event in Tree 
     *  @code
     *  Tree*  tree = ...
     *  LegendreSum4 s ( 5 , 3 , 2 , 2 , -1 , 1 , -2 , 2 , 0 , 4 , 0 , 1 ) ;
     *  DataParam::parameterize ( tree , s , "x" , "y" , "z" , "t" ) ;
     *  @endcode
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2019-07-3
     */ 
    static Ostap::StatusCode parameterize 
    ( TTree*                     tree        , 
      Ostap::Math::LegendreSum4& sum         , 
      const std::string&         xexpression , 
      const std::string&         yexpression , 
      const std::string&         zexpression , 
      const std::string&         uexpression , 
      const unsigned long        first       =  0 ,
      const unsigned long        last        = std::numeric_limits<unsigned long>::max() ) ;
    // ========================================================================
    /** fill Legendre sum with data from the Tree 
     *  @see Ostap::Math::LegendreSum4 
     *  @see Ostap::Math::LegendreSum4::fill
     *  @param tree        (INPUT)  the input tree 
     *  @param sum         (UPDATE) the parameterization object 
     *  @param xexpression (INPUT)  x-expression to be parameterized
     *  @param yexpression (INPUT)  y-expression to be parameterized
     *  @param zexpression (INPUT)  z-expression to be parameterized
     *  @param uexpression (INPUT)  u-expression to be parameterized
     *  @param selection   (INPUT)  selection/weight to be used 
     *  @param first       (INPUT)  the first event in Tree 
     *  @param last        (INPUT)  the last  event in Tree 
     *  @code
     *  Tree*  tree = ...
     *  LegendreSum4 s ( 5 , 2 , 4 ,  -1 , 1 , -4  , -5 , 0 , 3 ) ;
     *  DataParam::parameterize ( tree , s , "x"  , "y" , "z" , "t" , "q>10" ) ;
     *  @endcode
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2019-07-3
     */ 
    static Ostap::StatusCode parameterize 
    ( TTree*                     tree        ,  
      Ostap::Math::LegendreSum4& sum         , 
      const std::string&         xexpression , 
      const std::string&         yexpression , 
      const std::string&         zexpression , 
      const std::string&         uexpression , 
      const std::string&         selection   , 
      const unsigned long        first       =  0 ,
      const unsigned long        last        = std::numeric_limits<unsigned long>::max() ) ;
    // ==================================================================================
  } ; //                                      The end of class Ostap::DataParam 
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                     The END 
// ============================================================================
#endif // OSTAP_PARAMS_H
// ============================================================================
