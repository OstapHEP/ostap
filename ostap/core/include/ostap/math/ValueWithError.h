// ============================================================================
#ifndef OSTAP_VALUEWITH_ERROR_H
#define OSTAP_VALUEWITH_ERROR_H 1
// ============================================================================

#include <vector>
// ============================================================================
namespace ostap {
namespace math {
// ============================================================================

class ValueWithError {
 public:
  // ======================================================================
  typedef double Value;
  typedef double Covariance;
  typedef std::vector<ValueWithError> Vector;
  // ======================================================================
 public:
  // ======================================================================
  /// constructor from the value and covariance
  ValueWithError(const double value = 0, const double covariance = 0);

 public:  // trivial accessors
  // ======================================================================
  /// get the value
  double value() const { return m_value; }
  /// get the covariance
  double cov2() const { return m_cov2; }
  /// get the covariance
  double covariance() const { return m_cov2; }
  /** get the error
   *  @attention negative erorr is returned for invalid covariance
   *  @return the error estimate
   */
  double error() const;

 private:
  // ======================================================================
  /// the actual value
  double m_value;  //          the actual value
  /// the associated covariance
  double m_cov2;  // the associated covariance
  // ======================================================================
};
// ============================================================================
}
}
#endif  // OSTAP_VALUEWITH_ERROR_H
// ============================================================================
