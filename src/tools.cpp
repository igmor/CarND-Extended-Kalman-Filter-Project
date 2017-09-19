#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
	rmse << 0,0,0,0;

    VectorXd c(4);

	for(int i=0; i < estimations.size(); ++i){
	    c << estimations[i]-ground_truth[i];
	    for(int j=0; j < c.size(); ++j){
	        rmse[j] += c[j]*c[j]/estimations.size();
	    }
	}
    
    return rmse.array().sqrt();
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

	MatrixXd Hj(3,4);

	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//check division by zero
	if (px == 0 && py == 0) {
        std::cout << "ERROR: divison by zero" << std::endl << Hj << std::endl;
	}
	
	//compute the Jacobian matrix
	float f1_px_py = 1/(px*px + py*py);
	float f3_2_px_py = 1/(sqrt(px*px + py*py)*sqrt(px*px + py*py)*sqrt(px*px + py*py));
	
	
    Hj << px*sqrt(f1_px_py), py*sqrt(f1_px_py), 0, 0,
        -py*f1_px_py, px*f1_px_py, 0, 0,
        py*(vx*py-vy*px)*f3_2_px_py, px*(vy*px-vx*py)*f3_2_px_py, px*sqrt(f1_px_py), py*sqrt(f1_px_py);

	return Hj;    
}
