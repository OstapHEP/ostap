#---Create a shared library 
add_library(ostap SHARED src/format.cpp
                         src/gauss.cpp
                         src/AddBranch.cpp
                         src/AddVars.cpp
                         src/BLOB.cpp
                         src/BSpline.cpp
                         src/Bernstein.cpp
                         src/Bernstein1D.cpp
                         src/Bernstein2D.cpp
                         src/Bernstein3D.cpp
                         src/Binomial.cpp
                         src/BreitWigner.cpp
                         src/ChebyshevApproximation.cpp
                         src/Choose.cpp
                         src/Combine.cpp
                         src/Chi2Fit.cpp
                         src/Dalitz.cpp
                         src/DalitzIntegrator.cpp
                         src/DataFrameActions.cpp
                         src/DataFrameUtils.cpp
                         src/EigenSystem.cpp   
                         src/Error2Exception.cpp   
                         src/Exception.cpp
                         src/Faddeeva.cpp 
                         src/Formula.cpp   
                         src/FormulaVar.cpp   
                         src/Fourier.cpp   
                         src/Funcs.cpp   
                         src/GetWeight.cpp 
                         src/GSL_helpers.cpp              
                         src/GSL_sentry.cpp 
                         src/GSL_utils.cpp 
                         src/Hesse.cpp
                         src/HistoDump.cpp
                         src/HistoHash.cpp
                         src/HistoInterpolation.cpp
                         src/HistoInterpolators.cpp
                         src/HistoMake.cpp
                         src/HistoProject.cpp
                         src/HistoStat.cpp
                         src/IFuncs.cpp
                         src/Integrator.cpp
                         src/Interpolation.cpp
                         src/Iterator.cpp
                         src/Kinematics.cpp
                         src/KramersKronig.cpp
                         src/Lomont.cpp
                         src/LorentzVectorWithError.cpp
                         src/Math.cpp
                         src/MatrixUtils.cpp                         
                         src/Models.cpp
                         src/Models2D.cpp
                         src/Models3D.cpp
                         src/Moments.cpp
                         src/MoreMath.cpp
                         src/MoreRooFit.cpp
                         src/MoreVars.cpp
                         src/Mute.cpp
                         src/NStatEntity.cpp
                         src/Notifier.cpp
                         src/Ostap.cpp
                         src/OstapDataFrame.cpp
                         src/P2Quantile.cpp
                         src/Parameterization.cpp
                         src/Params.cpp
                         src/Peaks.cpp
                         src/PDFs.cpp
                         src/PDFs2D.cpp
                         src/PDFs3D.cpp
                         src/PhaseSpace.cpp
                         src/Piecewise.cpp
                         src/Point3DWithError.cpp
                         src/Polynomials.cpp   
                         src/Polarization.cpp
                         src/Primitives.cpp
                         src/Printable.cpp
                         src/PyBLOB.cpp
                         src/PyCallable.cpp 
                         src/PyFuncs.cpp 
                         src/PyPdf.cpp   
                         src/PyIterator.cpp
                         src/PySelector.cpp
                         src/PySelectorWithCuts.cpp
                         src/PyVar.cpp   
                         src/RootID.cpp
                         src/SFactor.cpp
                         src/StatEntity.cpp
                         src/StatVar.cpp
                         src/StatusCode.cpp
                         src/Tee.cpp
                         src/Tensors.cpp
                         src/Topics.cpp
                         src/Tmva.cpp
                         src/TreeGetter.cpp
                         src/UStat.cpp
                         src/Valid.cpp
                         src/ValueWithError.cpp
                         src/Vector3DWithError.cpp
                         src/Voigt.cpp
                         src/Workspace.cpp    
                         src/WStatEntity.cpp    
                         src/nSphere.cpp      
                         src/owens.cpp      
                         src/hcubature.cpp                         
                         src/pcubature.cpp
                        )

target_compile_features ( ostap PUBLIC cxx_constexpr                      )
target_compile_features ( ostap PUBLIC cxx_variadic_templates             )
target_compile_features ( ostap PUBLIC cxx_delegating_constructors        ) 
target_compile_features ( ostap PUBLIC cxx_defaulted_move_initializers    )
target_compile_features ( ostap PUBLIC cxx_decltype                       )
target_compile_features ( ostap PUBLIC cxx_decltype_auto                  )
target_compile_features ( ostap PUBLIC cxx_deleted_functions              )
target_compile_features ( ostap PUBLIC cxx_final                          )
target_compile_features ( ostap PUBLIC cxx_lambdas                        )
target_compile_features ( ostap PUBLIC cxx_inheriting_constructors        )
target_compile_features ( ostap PUBLIC cxx_override                       )
target_compile_features ( ostap PUBLIC cxx_range_for                      )
target_compile_features ( ostap PUBLIC cxx_static_assert                  )
target_compile_features ( ostap PUBLIC cxx_right_angle_brackets           )
target_compile_features ( ostap PUBLIC cxx_nullptr                        )
target_compile_features ( ostap PUBLIC cxx_auto_type                      )
target_compile_features ( ostap PUBLIC cxx_aggregate_default_initializers )
