#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  // check for invalid size
  if ((estimations.size() == 0)
        || (estimations.size() != ground_truth.size())) {
      return rmse;
  }

  // sum of squares
  for (unsigned int i = 0; i < estimations.size(); i++) {
    VectorXd diff = estimations[i] - ground_truth[i];
    diff = diff.array() * diff.array();
    rmse += diff;
  }
  // mean
  rmse /= estimations.size();
  // sqrt of mean
  rmse = rmse.array().sqrt();

  return rmse;
}
