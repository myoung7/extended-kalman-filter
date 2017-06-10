#include "kalman_filter.h"
#include <math.h>
#include <cmath>
#include "tools.h"
#include <stdio.h>
#include "measurement_package.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &R_in, MatrixXd &Q_in) {
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    R_ = R_in;
    Q_ = Q_in;
    
}

//Prediction Step of Kalman Filter.
void KalmanFilter::Predict() {
    x_ = F_*x_;
    P_ = F_*P_*F_.transpose()+Q_;
}

//Measurement (Update) Step of the Kalman Filter. Combines both Kalman Filter and Extended Kalman Filter.
void KalmanFilter::Update(const VectorXd &z, const MeasurementPackage::SensorType &type, MatrixXd &H) {
    
    MatrixXd I(4,4);
    VectorXd y;
    
    if (H.rows() != z.size()){
        printf("Error KalmanFilter::Update - # of H rows does not equal size of z.");
        throw exception();
        return;
    }
    
    //Identity Matrix
    I <<    1,0,0,0,
            0,1,0,0,
            0,0,1,0,
            0,0,0,1;
    
    //If sensor is a Radar, compute the Extended Kalman Filter.
    if(type == MeasurementPackage::RADAR){

        float px = x_(0);
        float py = x_(1);
        float vx = x_(2);
        float vy = x_(3);
        
        float h1 = sqrt(pow(px,2)+pow(py,2));
        float h2 = atan2(py,px);
        float h3 = (px*vx+py*vy)/sqrt(pow(px,2)+pow(py,2));
        
        
        VectorXd h(3);
        h << h1,h2,h3;
    
        //Calculate the Jacobian Matrix
        H = tools.CalculateJacobian(x_);
        
        y = z - h;
        
        //Brings 'phi' value to between -π and π.
        while(y(1) > M_PI){
            y(1) -= 2*M_PI;
        }
        
        while(y(1) < -M_PI){
            y(1) += 2*M_PI;
        }
        
    } else{
        y = z - H*x_;
    }
    
    MatrixXd S = H*P_*H.transpose()+R_;
    MatrixXd K = P_*H.transpose()*S.inverse();
    x_ = x_ + K*y;
    P_ = (I-K*H)*P_;
}


