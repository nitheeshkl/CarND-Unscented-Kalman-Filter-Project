#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

#define EPS 0.001 // a small value used to overcome div-by-zero

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  n_x_ = x_.size();
  n_aug_ = n_x_ + 2;
  n_sig_ = 2 * n_aug_ + 1;

  Xsig_pred_ = MatrixXd(n_x_, n_sig_);
  lambda_ = 3 - n_aug_;
  weights_ = VectorXd(n_sig_);

  R_lidar_ = MatrixXd(2, 2);
  R_lidar_ << std_laspx_*std_laspx_, 0,
              0, std_laspy_*std_laspy_;

  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0, std_radrd_*std_radrd_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  // Using CTRV model, state=[px, py, v, phi, phi_dot]
  if (!is_initialized_) {
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_[0], //px
            meas_package.raw_measurements_[1], //py
            0, //v, considered as zero since lidar doesn't provide
            0, //phi, considered as zero since lidar doesn't provide
            0; //phi_dot, considered as zero since lidar doesn't provide
      if (fabs(x_(0)) < EPS && fabs(x_(1)) < EPS) {
        x_(0) = EPS;
        x_(1) = EPS;
      }
    } else { // RADAR measuremenents
      const float rho = meas_package.raw_measurements_[0];
      const float phi = meas_package.raw_measurements_[1];
      const float rho_dot = meas_package.raw_measurements_[2];
      // conversion to state properties
      const float px = rho * cos(phi);
      const float py = rho * sin(phi);
      const float vx = rho_dot * cos(phi);
      const float vy = rho_dot * sin(phi);
      const float v = sqrt(vx*vx + vy*vy);
      x_ << px, py, v, 0, 0;
    }

    // setup weights for sigma points
    weights_(0) = lambda_ / (lambda_ + n_aug_);
    for (int i = 1; i < weights_.size(); i++) {
      weights_(i) = 0.5/(n_aug_ + lambda_);
    }

    time_us_ = meas_package.timestamp_;

    is_initialized_ = true;

    return;
  }

  // get time delta in seconds
  double dt = (meas_package.timestamp_ - time_us_)/1000000.0;
  time_us_ = meas_package.timestamp_;

  // predict
  Prediction(dt);

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double dt) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  const double dt_2 = dt*dt;

  VectorXd x_aug = VectorXd(n_aug_); // augmented state
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_); // augmented state covariance
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_); // sigma points

  // setup matrixes. reuse code from quizes
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;
  P_aug.fill(0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;
  MatrixXd L = P_aug.llt().matrixL(); // sqrt of P_aug

  // generate sigma points
  Xsig_aug.col(0) = x_aug; // first col is same as state
  const double sqrt_lambda_n_aug = sqrt(lambda_ + n_aug_);
  VectorXd sqrt_lambda_n_aug_L;
  for (unsigned int i = 0; i < n_aug_; i++) {
    sqrt_lambda_n_aug_L = sqrt_lambda_n_aug * L.col(i);
    Xsig_aug.col(i+1) = x_aug + sqrt_lambda_n_aug_L; // first half
    Xsig_aug.col(i+1 + n_aug_) = x_aug - sqrt_lambda_n_aug_L; // second half
  }

  // predict sigma points
  for (unsigned int i = 0; i < n_sig_; i++) {
    const double px = Xsig_aug(0, i);
    const double py = Xsig_aug(1, i);
    const double v = Xsig_aug(2, i);
    const double yaw = Xsig_aug(3, i);
    const double yawd = Xsig_aug(4, i);
    const double nu_a = Xsig_aug(5, i);
    const double nu_yawd = Xsig_aug(6, i);

    const double sin_yaw = sin(yaw);
    const double cos_yaw = cos(yaw);
    const double arg = yaw + yawd*dt;

    // predict state
    double px_pred, py_pred;
    // div-by-zero check
    if (fabs(yawd) > EPS) {
      const double v_yawd = v/yawd;
      px_pred = px + v_yawd * (sin(arg) - sin_yaw);
      py_pred = py + v_yawd * (cos_yaw - cos(arg));
    } else {
      const double v_dt = v*dt;
      px_pred = px + v_dt*cos_yaw;
      py_pred = py + v_dt*sin_yaw;
    }

    double v_pred = v;
    double yaw_pred = arg;
    double yawd_pred = yawd;

    // include noise
    px_pred += 0.5 * nu_a * dt_2 * cos_yaw;
    py_pred += 0.5 * nu_a * dt_2 * sin_yaw;
    v_pred += nu_a * dt;
    yaw_pred += 0.5 * nu_yawd * dt_2;
    yawd_pred += nu_yawd  * dt;

    Xsig_pred_(0, i) = px_pred;
    Xsig_pred_(1, i) = py_pred;
    Xsig_pred_(2, i) = v_pred;
    Xsig_pred_(3, i) = yaw_pred;
    Xsig_pred_(4, i) = yawd_pred;
  }

  // predict state mean
  x_ = Xsig_pred_ * weights_;
  // predict state covariance
  P_.fill(0.0);
  for (unsigned int i = 0; i < n_sig_; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = NormalizeAngle(x_diff(3));
    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}

double UKF::NormalizeAngle(double angle) {
  return atan2(sin(angle), cos(angle));
}
