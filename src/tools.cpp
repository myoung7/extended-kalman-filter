#include <iostream>
#include <stdio.h>
#include <math.h>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using namespace std;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
    VectorXd rmse(4);
    rmse << 0,0,0,0;
    
    if(estimations.size() != ground_truth.size() || estimations.size() == 0){
        cout << "Error Tools:CalculateRMSE() - Invalid estimation or ground_truth data" << endl;
        return rmse;
    }
    
    int n=estimations.size();
    
    for(int i=0;i<n;i++)
    {
        VectorXd residual = estimations[i]-ground_truth[i];
        residual = residual.array() * residual.array();
        rmse += residual;
    }
    
    rmse /= n;
    rmse = sqrt(rmse.array());
    
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    
    MatrixXd Hj(3,4);
    
    Hj <<   0,0,0,0,
            0,0,0,0,
            0,0,0,0;
    
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);
    
    //Calculates the denominator value for each partial derivative equation.
    //Used for checking if 0 is in the denominator.
    
    float rho_dpx_den = sqrt(pow(px,2)+pow(py,2));
    float rho_dpy_den = sqrt(pow(px,2)+pow(py,2));
    float phi_dpx_den = pow(px,2)+pow(py,2);
    float phi_dpy_den = pow(px,2)+pow(py,2);
    float rhodot_dpx_den = pow(pow(px,2)+pow(py,2),1.5);
    float rhodot_dpy_den = pow(pow(px,2)+pow(py,2),1.5);
    float rhodot_dvx_den = sqrt(pow(px,2)+pow(py,2));
    float rhodot_dvy_den = sqrt(pow(px,2)+pow(py,2));
    
    //If product_check results in 0, then one of the denominators equated to 0.
    float product_check = rho_dpx_den*rho_dpy_den*phi_dpx_den*phi_dpy_den*rhodot_dpx_den
    *rhodot_dpy_den*rhodot_dvx_den*rhodot_dvy_den;
    
    //Checks if product_check equated to 0. If so, returns an error message.
    if(product_check != 0)
    {
        float rho_dpx = px/rho_dpx_den;
        float rho_dpy = py/rho_dpy_den;
        
        float phi_dpx = -1 * (py/phi_dpx_den);
        float phi_dpy = px/phi_dpy_den;
        
        float rhodot_dpx = (py*(vx*py-vy*px))/rhodot_dpx_den;
        float rhodot_dpy = (px*(vy*px-vx*py))/rhodot_dpy_den;
        float rhodot_dvx = px/rhodot_dvx_den;
        float rhodot_dvy = py/rhodot_dvy_den;
        
        Hj <<   rho_dpx, rho_dpy, 0, 0,
                phi_dpx, phi_dpy, 0, 0,
                rhodot_dpx, rhodot_dpy, rhodot_dvx, rhodot_dvy;
        
    } else{
        cout << "Error Tools::CalculateJacobian() - 0 found in Denominator. Can't divide by 0.";
    }
    
    return Hj;
    
}
