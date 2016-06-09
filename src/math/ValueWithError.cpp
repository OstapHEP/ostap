#include "ostap/math/ValueWithError.h"

#include <cmath>

// ============================================================================
// constructor from the value and covariance
// ============================================================================
ostap::math::ValueWithError::ValueWithError(const double value,
                                            const double covariance)
    : m_value(value), m_cov2(covariance) {
  // need _zero
}

double ostap::math::ValueWithError::error() const {
  // zero
  return 0 <= m_cov2 ? std::sqrt(m_cov2) : -std::sqrt(-m_cov2);
}