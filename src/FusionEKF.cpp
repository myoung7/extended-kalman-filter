#include "FusionEKF.h"
#include "tools.h"
#include <math.h>
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;

    previous_timestamp_ = 0;
    
    noise_ax = 9;
    noise_ay = 9;

    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    Hj_ = MatrixXd(3, 4);

    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
        0, 0.0225;

    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

    //Laser transformation matrix
    H_laser_ << 1,0,0,0,
                0,1,0,0;
    //Jacobian transformation matrix (used for Radar)
    Hj_ <<  1,0,0,0,
            0,1,0,0,
            0,0,1,0;
    
    VectorXd x = VectorXd(4);
    
    MatrixXd F = MatrixXd(4,4);
    
    MatrixXd P = MatrixXd(4,4);
    
    //Covariance Matrix P
    P <<    1,0,0,0,
            0,1,0,0,
            0,0,1000,0,
            0,0,0,1000;
    
    //Noise process covariance matrix
    MatrixXd Q = MatrixXd(4,4);

    ekf_.Init(x, P, F, R_laser_, Q);
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
    if (!is_initialized_) {
        // first measurement
        cout << "EKF: " << endl;
      
        float px, py, vx, vy = 0;
        previous_timestamp_ = measurement_pack.timestamp_;

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            /**
             Convert radar from polar to cartesian coordinates and initialize state for Radar.
             */
            float rho = measurement_pack.raw_measurements_[0];
            float phi = measurement_pack.raw_measurements_[1];
            float rhodot = measurement_pack.raw_measurements_[2];
        
            px = rho*cos(phi);
            py = rho*sin(phi);
            vx = 0;
            vy = 0;
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            /**
             Initialize state for Laser.
             */
            px = measurement_pack.raw_measurements_[0];
            py = measurement_pack.raw_measurements_[1];
        }
      
        ekf_.x_ << px, py, vx, vy;

    // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }
    
  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

    float dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;
    
    ekf_.F_ <<  1,0,dt,0,
                0,1,0,dt,
                0,0,1,0,
                0,0,0,1;
    
    ekf_.Q_ <<  (noise_ax*pow(dt,4))/4, 0, (noise_ax*pow(dt,3))/2, 0,
                0, (noise_ay*pow(dt,4))/4, 0, (noise_ay*pow(dt,3))/2,
                (noise_ax*pow(dt,3)), 0, (noise_ax*pow(dt,2)), 0,
                0, (noise_ay*pow(dt,3))/2, 0, noise_ay*pow(dt,2);

    ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
      VectorXd z = VectorXd(3);
      
      float z1 = measurement_pack.raw_measurements_[0];
      float z2 = measurement_pack.raw_measurements_[1];
      float z3 = measurement_pack.raw_measurements_[2];
      
      z << z1, z2, z3;
      ekf_.R_ = R_radar_;
      
      try
      {
          ekf_.Update(z, measurement_pack.sensor_type_, Hj_);
      }
      catch(const exception&)
      {
          throw exception();
      }
  } else {
    // Laser updates
      VectorXd z = VectorXd(2);
      
      float z1 = measurement_pack.raw_measurements_[0];
      float z2 = measurement_pack.raw_measurements_[1];
      
      z << z1, z2;
      ekf_.R_ = R_laser_;
      try
      {
          ekf_.Update(z, measurement_pack.sensor_type_, H_laser_);
      } catch(const exception&)
      {
          throw exception();
      }
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
