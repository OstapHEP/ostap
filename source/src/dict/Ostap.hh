// ============================================================================
#ifndef OSTAP_OSTAP_HH 
#define OSTAP_OSTAP_HH 1
// ============================================================================
// #include "Python.h"
// ============================================================================
// Include files
// ============================================================================
#include "Ostap/AddBranch.h"
#include "Ostap/AddVars.h"
#include "Ostap/BLOB.h"
#include "Ostap/BSpline.h"
#include "Ostap/Bernstein.h"
#include "Ostap/Bernstein1D.h"
#include "Ostap/Bernstein2D.h"
#include "Ostap/Bernstein3D.h"
#include "Ostap/Binomial.h"
#include "Ostap/BreitWigner.h"
#include "Ostap/Bit.h"
#include "Ostap/ChebyshevApproximation.h"
#include "Ostap/Chi2Fit.h"
#include "Ostap/Chi2Solution.h"
#include "Ostap/Choose.h"
#include "Ostap/Clenshaw.h"
#include "Ostap/Combine.h"
#include "Ostap/Dalitz.h"
#include "Ostap/DalitzIntegrator.h"
#include "Ostap/DataFrameActions.h"
#include "Ostap/DataFrameUtils.h"
#include "Ostap/Digit.h"
#include "Ostap/EigenSystem.h"
#include "Ostap/Exception.h"
#include "Ostap/Error2Exception.h"
#include "Ostap/Formula.h"
#include "Ostap/FitResult.h"
#include "Ostap/FormulaVar.h"
#include "Ostap/Fourier.h"
#include "Ostap/Funcs.h"
#include "Ostap/GenericMatrixTypes.h"
#include "Ostap/GenericVectorTypes.h"
#include "Ostap/GeomFun.h"
#include "Ostap/GetWeight.h"
#include "Ostap/GSL_utils.h"
#include "Ostap/Hesse.h"
#include "Ostap/HFuncs.h"
#include "Ostap/Hash.h"
#include "Ostap/HistoDump.h"
#include "Ostap/HistoHash.h"
#include "Ostap/HistoInterpolation.h"
#include "Ostap/HistoInterpolators.h"
#include "Ostap/HistoMake.h"
#include "Ostap/HistoProject.h"
#include "Ostap/HistoStat.h"
#include "Ostap/KramersKronig.h"
#include "Ostap/Interpolation.h"
#include "Ostap/Interpolants.h"
#include "Ostap/Iterator.h"
#include "Ostap/Line.h"
#include "Ostap/LineTypes.h"
#include "Ostap/Lomont.h"
#include "Ostap/LorentzVectorWithError.h"
#include "Ostap/Kinematics.h"
#include "Ostap/Math.h"
#include "Ostap/MatrixUtils.h"
#include "Ostap/MatrixUtils2.h"
#include "Ostap/MatrixUtilsT.h"
#include "Ostap/MatrixTransforms.h"
#include "Ostap/Models.h"
#include "Ostap/Models2D.h"
#include "Ostap/Models3D.h"
#include "Ostap/Moments.h"
#include "Ostap/MoreMath.h"
#include "Ostap/MoreRooFit.h"
#include "Ostap/MoreVars.h"
#include "Ostap/Mute.h"
#include "Ostap/Notifier.h"
#include "Ostap/NSphere.h"
#include "Ostap/NStatEntity.h"
#include "Ostap/Parameterization.h"
#include "Ostap/Params.h"
#include "Ostap/Peaks.h"
#include "Ostap/PDFs.h"
#include "Ostap/PDFs2D.h"
#include "Ostap/PDFs3D.h"
#include "Ostap/PhaseSpace.h"
#include "Ostap/Piecewise.h"
#include "Ostap/Plane3DTypes.h"
#include "Ostap/Point3DTypes.h"
#include "Ostap/Point3DWithError.h"
#include "Ostap/Polynomials.h"
#include "Ostap/Power.h"
#include "Ostap/Primitives.h"
#include "Ostap/Printable.h"
#include "Ostap/PyCallable.h"   
#include "Ostap/PyFuncs.h"   
#include "Ostap/PyIntegrator.h"   
#include "Ostap/PyIterator.h"
#include "Ostap/PyPdf.h"     
#include "Ostap/PySelector.h"
#include "Ostap/PySelectorWithCuts.h"
#include "Ostap/PyVar.h"     
#include "Ostap/PyBLOB.h"
#include "Ostap/Polarization.h"
#include "Ostap/ProgressBar.h"
#include "Ostap/RootID.h"
#include "Ostap/SFactor.h"
#include "Ostap/StatEntity.h"
#include "Ostap/StatVar.h"
#include "Ostap/StatusCode.h"
#include "Ostap/SVectorWithError.h"
#include "Ostap/SymmetricMatrixTypes.h"
#include "Ostap/Tensors.h"
#include "Ostap/Tee.h"
#include "Ostap/ToStream.h"
#include "Ostap/Topics.h"
#include "Ostap/TreeGetter.h"
#include "Ostap/TypeWrapper.h"
#include "Ostap/Tmva.h"
#include "Ostap/qMath.h"
#include "Ostap/Valid.h"
#include "Ostap/ValueWithError.h"
#include "Ostap/Vector3DTypes.h"
#include "Ostap/Vector3DWithError.h"
#include "Ostap/Vector4DTypes.h"
#include "Ostap/Voigt.h"
#include "Ostap/UStat.h"
#include "Ostap/WStatEntity.h"
#include "Ostap/Workspace.h"
// ============================================================================
#include "RooFormulaVar.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    template <typename aPoint, typename aLine, typename aPlane>
    struct GF
    {
      static bool intersection( const aLine&  line      ,
                                const aPlane& plane     ,
                                aPoint&       intersect ,
                                double&       mu        )
      { return Ostap::Math::intersection<aLine, aPlane, aPoint>(line,
                                                                plane,
                                                                intersect,
                                                                mu); }
      static bool intersection( const aPlane& plane0    ,
                                const aPlane& plane1    ,
                                aLine&        intersect )
      { return Ostap::Math::intersection<aLine, aPlane>(plane0,
                                                        plane1,
                                                        intersect); }
      static  bool intersection( const aPlane& plane0    ,
                                 const aPlane& plane1    ,
                                 const aPlane& plane2    ,
                                 aPoint&       intersect )
      { return Ostap::Math::intersection<aPoint, aPlane>(plane0,
                                                         plane1,
                                                         plane2,
                                                         intersect); }
      static double impactParameter(const aPoint&  point ,
                                    const aLine&   line  )
      { return Ostap::Math::impactParameter<aPoint,aLine>(point, line); }
      
      static double distance( const aLine& line0 ,
                              const aLine& line1 )
      { return Ostap::Math::distance<aLine, aLine>(line0, line1); }
      
      static aPoint project ( const aPlane& plane , 
                              const aPoint& point ) 
      { return plane.ProjectOntoPlane ( point ) ; }
      
      static  bool closestPoints( const aLine& line0 ,
                                  const aLine& line1 ,
                                  aPoint&      p0    ,
                                  aPoint&      p1    )
      { return Ostap::Math::closestPoints<aLine, aLine, aPoint>(line0,
                                                                line1,
                                                                p0,
                                                                p1); }
      static double closestPointParam( const aPoint&  point ,
                                       const aLine&   line  )
      { return Ostap::Math::closestPointParam<aLine, aPoint>(point, line); }
      
      static aPoint closestPoint(const aPoint&  point ,
                                 const aLine& line)
      { return Ostap::Math::closestPoint<aLine, aPoint>(point, line); }
      
      static bool closestPointParams( const aLine& line0 ,
                                      const aLine& line1 ,
                                      double&      mu0   ,
                                      double&      mu1   )
      { return Ostap::Math::closestPointParams<aLine, aLine>(line0,
                                                             line1,
                                                             mu0,
                                                             mu1);
        
      }
      static bool parallel( const aLine& line0 ,
                            const aLine& line1 )
      { return Ostap::Math::parallel<aLine, aLine>(line0, line1); }
      
    };
    typedef GF<XYZPoint, XYZLine, Plane3D> XYZGeomFun;
    // ========================================================================
    class EigenSystems
    {
      // ======================================================================
    public : // eigen values
      // ======================================================================
      // 2x2
      static Ostap::Vector2 eigenValues
      ( const Ostap::SymMatrix2x2& mtrx          ,
        const bool                 sorted = true )
      {
        Ostap::Math::GSL::EigenSystem system ;
        return system.eigenValues ( mtrx , sorted ) ;
      }
      // 3x3
      static Ostap::Vector3 eigenValues
      ( const Ostap::SymMatrix3x3& mtrx          ,
        const bool                 sorted = true )
      {
        Ostap::Math::GSL::EigenSystem system ;
        return system.eigenValues ( mtrx , sorted ) ;
      }
      // 4x4
      static Ostap::Vector4 eigenValues
      ( const Ostap::SymMatrix4x4& mtrx          ,
        const bool                 sorted = true )
      {
        Ostap::Math::GSL::EigenSystem system ;
        return system.eigenValues ( mtrx , sorted ) ;
      }
      // ======================================================================
    public: // eigen vectors
      // ======================================================================
    public: // eigen vectors
      // ======================================================================
      // 2x2
      static StatusCode eigenVectors
      ( const Ostap::SymMatrix2x2&   mtrx          ,
        Ostap::Vector2&              vals          ,
        std::vector<Ostap::Vector2>& vecs          ,
        const bool                   sorted = true )
      {
        Ostap::Math::GSL::EigenSystem system ;
        return system.eigenVectors ( mtrx , vals , vecs , sorted ) ;
      }
      // 3x3
      static StatusCode eigenVectors
      ( const Ostap::SymMatrix3x3&   mtrx          ,
        Ostap::Vector3&              vals          ,
        std::vector<Ostap::Vector3>& vecs          ,
        const bool                   sorted = true )
      {
        Ostap::Math::GSL::EigenSystem system ;
        return system.eigenVectors ( mtrx , vals , vecs , sorted ) ;
      }
      // 4x4
      static StatusCode eigenVectors
      ( const Ostap::SymMatrix4x4&   mtrx          ,
        Ostap::Vector4&              vals          ,
        std::vector<Ostap::Vector4>& vecs          ,
        const bool                   sorted = true )
      {
        Ostap::Math::GSL::EigenSystem system ;
        return system.eigenVectors ( mtrx , vals , vecs , sorted ) ;
      }
      // ======================================================================
    } ;
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================

namespace
{
  // ==========================================================================
  struct Ostap_Instantiations
  {
    Ostap_Instantiations();
    //
    Ostap::XYZLine               __lineXYZ;
    Ostap::XYZVector             __vectXYZ;
    Ostap::Plane3D              __planeXYZ;
    Ostap::XYZPoint             __pointXYZ;
    //
    Ostap::Math::XYZGeomFun __geomFunXYZ;
    //
    Ostap::Math::SVectorWithError<2,double> __sv2 ;
    Ostap::Math::SVectorWithError<3,double> __sv3 ;
    Ostap::Math::SVectorWithError<4,double> __sv4 ;
    Ostap::Math::SVectorWithError<5,double> __sv5 ;
    Ostap::Math::SVectorWithError<6,double> __sv6 ;
    Ostap::Math::SVectorWithError<8,double> __sv8 ;
    //
    std::vector<Ostap::Math::ValueWithError>               _dver1 ;
    std::vector<std::vector<Ostap::Math::ValueWithError> > _dver2 ;
    std::vector<Ostap::WStatEntity>                        _dver3 ;
    //
    std::vector<Ostap::Vector2>  _vct_2 ;
    std::vector<Ostap::Vector3>  _vct_3 ;
    std::vector<Ostap::Vector3>  _vct_4 ;
    //
    Ostap::Math::Chi2Solution<4,2>          __cs11 ;
    Ostap::Math::Chi2Solution<4,2>::DATA    __cs21 ;
    Ostap::Math::Chi2Solution<4,2>::COV2    __cs31 ;
    Ostap::Math::Chi2Solution<4,2>::CMTRX2  __cs41 ;
    Ostap::Math::Chi2Solution<4,2>::COFF    __cs51 ;
    Ostap::Math::Chi2Solution<4,2>::VECT    __cs61 ;

    Ostap::Math::Chi2Solution<6,2>          __cs12 ;
    Ostap::Math::Chi2Solution<6,2>::DATA    __cs22 ;
    Ostap::Math::Chi2Solution<6,2>::COV2    __cs32 ;
    Ostap::Math::Chi2Solution<6,2>::CMTRX2  __cs42 ;
    Ostap::Math::Chi2Solution<6,2>::COFF    __cs52 ;
    Ostap::Math::Chi2Solution<6,2>::VECT    __cs62 ;
    
    Ostap::Math::Equal_To<double>               __eq_1 ;
    Ostap::Math::Equal_To<std::vector<double> > __eq_2 ;
    Ostap::Math::Zero<double>                   __eq_3 ;
    Ostap::Math::Zero<std::vector<double> >     __eq_4 ;
    Ostap::Math::NotZero<double>                __eq_5 ;
    Ostap::Math::NotZero<std::vector<double> >  __eq_6 ;
    //
    Ostap::Math::LessOrEqual<double>     __eq_7 ;
    Ostap::Math::GreaterOrEqual<double>  __eq_8 ;
    //
    // std::function<std::complex<double>(double)>               __ff0 ;
    // std::function<std::complex<double>(double,double)>        __ff1 ;
    // std::function<std::complex<double>(double,double,double)> __ff2 ;
    
  };
  // ==========================================================================
} //                                             The end of anonymous namespace 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_OSTAP_HH
// ============================================================================
