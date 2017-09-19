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

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */


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
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    cout << "EKF:  measurement_pack.raw_measurements_: " <<  measurement_pack.raw_measurements_ << endl;

    ekf_.x_  = VectorXd(4);
	ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    previous_timestamp_ = measurement_pack.timestamp_;    

    cout << "EKF(1) " << endl;

    //state transition matrix
	ekf_.F_ = MatrixXd(4, 4);    
	ekf_.F_ << 1, 0, 1, 0,
			   0, 1, 0, 1,
			   0, 0, 1, 0,
			   0, 0, 0, 1;

    cout << "EKF(2) " << endl;

    //covariance matrix
	ekf_.P_ = MatrixXd(4, 4);    
	ekf_.P_ << 1, 0, 0, 0,
			   0, 1, 0, 0,
			   0, 0, 1000, 0,
			   0, 0, 0, 1000;

    cout << "EKF(3) " << endl;

	//measurement covariance
	ekf_.R_ = MatrixXd(2, 2);    
	ekf_.R_ << 0.0225, 0,
			   0, 0.0225;

    cout << "EKF(4) " << endl;

	//measurement matrix
	ekf_.H_ = MatrixXd(2, 4);        
	ekf_.H_ << 1, 0, 0, 0,
			   0, 1, 0, 0;

	ekf_.Q_ = MatrixXd(4, 4);        

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        ekf_.R_ = R_laser_; //FIXME!!!!!!!!!!!!!!!!!!!!!!!!! not laser
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        ekf_.R_ = R_laser_;
        
      /**
      Initialize state.
      */
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

  cout << "EKF(5) " << endl;

  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0, dt,
             0, 0, 1, 0,
             0, 0, 0, 1;

  cout << "EKF(6) " << endl;

  ekf_.Q_ << dt*dt*dt*dt/4*noise_ax, 0, dt*dt*dt/2*noise_ax, 0,
			 0, dt*dt*dt*dt/4*noise_ay, 0, dt*dt*dt/2*noise_ay,
			 dt*dt*dt/2*noise_ax, 0, dt*dt*noise_ax, 0,
			 0, dt*dt*dt/2*noise_ay, 0, dt*dt*noise_ay;

  cout << "EKF(7) " << endl;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  cout << "EKF(8) " << endl;

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
      //ekf_.Update(measurement_pack.raw_measurements_.transpose());      
  } else {      
    // Laser updates
      ekf_.Update(measurement_pack.raw_measurements_.transpose());      
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
