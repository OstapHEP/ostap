#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// #pragma link C++ namespace ostap;
// #pragma link C++ namespace ostap::math;

#pragma link C++ class ostap::math::ValueWithError;
#pragma link C++ class ostap::math::Equal_To<float>+;

#pragma link C++ function ostap::math::absMin<float>;
#pragma link C++ function ostap::math::pow;
#pragma link C++ function ostap::math::knuth_equal_to_double;
#pragma link C++ function ostap::math::lomont_compare_double;
#pragma link C++ function ostap::math::next_double;


#endif  // __CINT__
