#pragma once

#include <boost/math/policies/error_handling.hpp>
#include <boost/math/policies/policy.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <cerrno>


namespace sf::math {

// FIXME: Check if long double actually is larger than double.
typedef boost::multiprecision::cpp_bin_float_oct p_value_type;


struct result {
  p_value_type value{};
  int error{};

  operator bool() const { return 0 == error; }
};


inline result beta_incomplete(double aa, double bb, double xx) {
  namespace bm = boost::math;
  namespace bmp = boost::math::policies;

  typedef bmp::policy<bmp::domain_error<bmp::errno_on_error>,
                      bmp::pole_error<bmp::errno_on_error>,
                      bmp::overflow_error<bmp::errno_on_error>,
                      bmp::evaluation_error<bmp::errno_on_error>>
      error_policy;

  typedef p_value_type parameter_type;

  parameter_type const aa_{aa};
  parameter_type const bb_{bb};
  parameter_type const xx_{xx};

  errno = 0; // Clear to make sure.
  auto const res{bm::ibeta(aa_, bb_, xx_, error_policy{})};
  p_value_type const res_(
      res); // For narrowing in case we use some other parameter type.
  return result{res_, errno};
}

}  // namespace sf::math
