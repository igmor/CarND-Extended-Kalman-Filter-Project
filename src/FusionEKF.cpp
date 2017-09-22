#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
int iter = 0;

FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

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


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

   float noise_ax = 9.0, noise_ay = 9.0;


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    */
    // first measurement

    previous_timestamp_ = measurement_pack.timestamp_;    

    //state transition matrix
	ekf_.F_ = MatrixXd(4, 4);    
	ekf_.F_ << 1, 0, 1, 0,
			   0, 1, 0, 1,
			   0, 0, 1, 0,
			   0, 0, 0, 1;

    //covariance matrix
	ekf_.P_ = MatrixXd(4, 4);    
	ekf_.P_ << 1, 0, 0, 0,
			   0, 1, 0, 0,
			   0, 0, 1000, 0,
			   0, 0, 0, 1000;

	//measurement covariance
	ekf_.R_ = MatrixXd(2, 2);    
	ekf_.R_ << 0.0225, 0,
			   0, 0.0225;

	//measurement matrix
	H_laser_ << 1, 0, 0, 0,
			    0, 1, 0, 0;

	ekf_.Q_ = MatrixXd(4, 4);        

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        
        cout <<  measurement_pack.raw_measurements_ << endl;

        VectorXd radar_m = VectorXd(4);
        
        float phi = measurement_pack.raw_measurements_[1];
        float ro = measurement_pack.raw_measurements_[0];

        ekf_.x_  = VectorXd(4);        
        ekf_.x_ << ro*cos(phi), ro*sin(phi), 0, 0;

        Hj_ = tools.CalculateJacobian(ekf_.x_);

        ekf_.R_ = R_radar_;
        ekf_.H_ = Hj_;        
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        ekf_.R_ = R_laser_;
        ekf_.H_ = H_laser_;        

        ekf_.x_  = VectorXd(4);
        ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;


  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0, dt,
             0, 0, 1, 0,
             0, 0, 0, 1;

  ekf_.Q_ << dt*dt*dt*dt/4*noise_ax, 0, dt*dt*dt/2*noise_ax, 0,
			 0, dt*dt*dt*dt/4*noise_ay, 0, dt*dt*dt/2*noise_ay,
			 dt*dt*dt/2*noise_ax, 0, dt*dt*noise_ax, 0,
			 0, dt*dt*dt/2*noise_ay, 0, dt*dt*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
      ekf_.R_ = R_radar_;
      ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
      
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);
      
  } else {      
    // Laser updates
      ekf_.R_ = R_laser_;
      ekf_.H_ = H_laser_;        
      
      ekf_.Update(measurement_pack.raw_measurements_);      
  }

  // print the output
  //  cout << iter << ":" << "measurement_pack:" << measurement_pack.raw_measurements_ << endl;
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;

  iter++;
}
