#pragma once

#include <boost/math/distributions/binomial.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/policies/policy.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <cerrno>


namespace sf::math::detail {

typedef boost::math::policies::policy<
    boost::math::policies::domain_error<boost::math::policies::errno_on_error>,
    boost::math::policies::pole_error<boost::math::policies::errno_on_error>,
    boost::math::policies::overflow_error<
        boost::math::policies::errno_on_error>,
    boost::math::policies::evaluation_error<
        boost::math::policies::errno_on_error>>
    error_policy;
}


namespace sf::math {

typedef boost::multiprecision::cpp_bin_float_oct p_value_type;


struct result {
  p_value_type value{};
  int error{};

  operator bool() const { return 0 == error; }
};


inline result binomial_cdf_upper_tail(double kk, double pp, double nn) {
  namespace bm = boost::math;

  typedef p_value_type parameter_type;

  parameter_type const nn_{nn};
  parameter_type const pp_{pp};
  parameter_type const kk_{kk};

  errno = 0; // Clear to make sure.
  auto const res{cdf(complement(
      bm::binomial_distribution<p_value_type, detail::error_policy>{nn_, pp_},
      kk_))}; // Use ADL.
  p_value_type const res_(
      res); // For narrowing in case we use some other parameter type.
  return result{res_, errno};
}


inline result beta_incomplete(double aa, double bb, double xx) {
  namespace bm = boost::math;

  typedef p_value_type parameter_type;

  parameter_type const aa_{aa};
  parameter_type const bb_{bb};
  parameter_type const xx_{xx};

  errno = 0; // Clear to make sure.
  auto const res{bm::ibeta(aa_, bb_, xx_, detail::error_policy{})};
  p_value_type const res_(
      res); // For narrowing in case we use some other parameter type.
  return result{res_, errno};
}

}  // namespace sf::math
