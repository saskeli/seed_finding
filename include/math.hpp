#pragma once

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_result.h>


namespace sf::math {

inline double beta_incomplete(double aa, double bb, double xx) {
  gsl_sf_result res;
  int err = gsl_sf_beta_inc_e(aa, bb, xx, &res);
  return err ? 0 : res.val;
}

}  // namespace sf::math
