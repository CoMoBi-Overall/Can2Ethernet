// using class makes the custom function for kinematics and jacobian
//---------------------------- Eigen ---------------------------------//
// #include <stdio.h>
// #include <iostream>
// #include <eigen3/Eigen/Eigen>
// #include <eigen3/Eigen/Dense>
// #include <cmath>
#include "yoon_indy7_func.h"

namespace yoon_indy7_kinematics
{
    // ---- Parameters ---- //  
    // DH-parameter
    double d1=0.3, d4=0.35, d5=0.1865, d6=0.228;
    double a2=0.45;
    double th1, th2, th3, th4, th5, th6;
    Eigen::Matrix4f T01, T12, T23, T34, T45, T56, T06;

    double Ex,Ey,Ez;
    double Eal,Ebe,Ega;
    Eigen::Matrix3f Euler2rot = Eigen::Matrix3f::Zero(3,3);
      
    // ---- Functions ---- //
    
    void Y_indy7_forw(double* th, Eigen::Matrix4f& posture)
    {

        T01 << cos(th[0]), -sin(th[0]),  0,  0,
               sin(th[0]),  cos(th[0]),  0,  0,
                        0,           0,  1, d1,
                        0,           0,  0,  1;

        T12 << -sin(th[1]), -cos(th[1]),  0,  0,
                        0,           0, -1,  0,
                cos(th[1]), -sin(th[1]),  0,  0,
                        0,           0,  0,  1;

        T23 << -sin(th[2]), -cos(th[2]),  0, a2,
                cos(th[2]), -sin(th[2]),  0,  0,
                        0,           0,  1,  0,
                        0,           0,  0,  1;

        T34 << -cos(th[3]),  sin(th[3]),  0,  0,
                        0,           0, -1,-d4,
            -sin(th[3]), -cos(th[3]),  0,  0,
                        0,           0,  0,  1;

        T45 <<  cos(th[4]), -sin(th[4]),  0,  0,
                        0,           0, -1,-d5,
                sin(th[4]),  cos(th[4]),  0,  0,
                        0,           0,  0,  1;

        T56 <<  cos(th[5]), -sin(th[5]),  0,  0,
                        0,           0,  1, d6,
            -sin(th[5]), -cos(th[5]),  0,  0,
                        0,           0,  0,  1;
        
        T06 = T01*T12*T23*T34*T45*T56;
        
        posture = T06;
    }
    void Y_indy7_jacob(double* th, Eigen::MatrixXf& J06)
    {
        th1 = th[0];
        th2 = th[1];
        th3 = th[2];
        th4 = th[3];
        th5 = th[4];
        th6 = th[5];
        J06 << d5*(cos(th1)*cos(th4) - cos(th2 + th3)*sin(th1)*sin(th4)) + d6*(sin(th5)*(cos(th1)*sin(th4) + cos(th2 + th3)*cos(th4)*sin(th1)) + sin(th2 + th3)*cos(th5)*sin(th1)) + sin(th1)*(d4*sin(th2 + th3) + a2*sin(th2)), -cos(th1)*(d4*cos(th2 + th3) + a2*cos(th2) + d6*(cos(th2 + th3)*cos(th5) - sin(th2 + th3)*cos(th4)*sin(th5)) + d5*sin(th2 + th3)*sin(th4)), -cos(th1)*(d4*cos(th2 + th3) + d6*(cos(th2 + th3)*cos(th5) - sin(th2 + th3)*cos(th4)*sin(th5)) + d5*sin(th2 + th3)*sin(th4)), d5*cos(th2 + th3)*cos(th1)*cos(th4) - d5*sin(th1)*sin(th4) + d6*cos(th4)*sin(th1)*sin(th5) + d6*cos(th2 + th3)*cos(th1)*sin(th4)*sin(th5),  d6*(cos(th5)*sin(th1)*sin(th4) + sin(th2 + th3)*cos(th1)*sin(th5) - cos(th2 + th3)*cos(th1)*cos(th4)*cos(th5)),0,
               d5*(cos(th4)*sin(th1) + cos(th2 + th3)*cos(th1)*sin(th4)) + d6*(sin(th5)*(sin(th1)*sin(th4) - cos(th2 + th3)*cos(th1)*cos(th4)) - sin(th2 + th3)*cos(th1)*cos(th5)) - cos(th1)*(d4*sin(th2 + th3) + a2*sin(th2)), -sin(th1)*(d4*cos(th2 + th3) + a2*cos(th2) + d6*(cos(th2 + th3)*cos(th5) - sin(th2 + th3)*cos(th4)*sin(th5)) + d5*sin(th2 + th3)*sin(th4)), -sin(th1)*(d4*cos(th2 + th3) + d6*(cos(th2 + th3)*cos(th5) - sin(th2 + th3)*cos(th4)*sin(th5)) + d5*sin(th2 + th3)*sin(th4)), d5*cos(th1)*sin(th4) + d5*cos(th2 + th3)*cos(th4)*sin(th1) - d6*cos(th1)*cos(th4)*sin(th5) + d6*cos(th2 + th3)*sin(th1)*sin(th4)*sin(th5), -d6*(cos(th1)*cos(th5)*sin(th4) - sin(th2 + th3)*sin(th1)*sin(th5) + cos(th2 + th3)*cos(th4)*cos(th5)*sin(th1)),0,
               0, d5*cos(th2)*cos(th3)*sin(th4) - d4*cos(th2)*sin(th3) - d4*cos(th3)*sin(th2) - a2*sin(th2) - d6*cos(th2)*cos(th5)*sin(th3) - d6*cos(th3)*cos(th5)*sin(th2) - d5*sin(th2)*sin(th3)*sin(th4) - d6*cos(th2)*cos(th3)*cos(th4)*sin(th5) + d6*cos(th4)*sin(th2)*sin(th3)*sin(th5), d5*cos(th2)*cos(th3)*sin(th4) - d4*cos(th3)*sin(th2) - d4*cos(th2)*sin(th3) - d6*cos(th2)*cos(th5)*sin(th3) - d6*cos(th3)*cos(th5)*sin(th2) - d5*sin(th2)*sin(th3)*sin(th4) - d6*cos(th2)*cos(th3)*cos(th4)*sin(th5) + d6*cos(th4)*sin(th2)*sin(th3)*sin(th5), sin(th2 + th3)*(d5*cos(th4) + d6*sin(th4)*sin(th5)), -d6*(cos(th2 + th3)*sin(th5) + sin(th2 + th3)*cos(th4)*cos(th5)), 0,
               0, sin(th1), sin(th1), -sin(th2 + th3)*cos(th1), cos(th4)*sin(th1) + cos(th2 + th3)*cos(th1)*sin(th4),   sin(th5)*(sin(th1)*sin(th4) - cos(th2 + th3)*cos(th1)*cos(th4)) - sin(th2 + th3)*cos(th1)*cos(th5),
               0, -cos(th1), -cos(th1), -sin(th2 + th3)*sin(th1), cos(th2 + th3)*sin(th1)*sin(th4) - cos(th1)*cos(th4), - sin(th5)*(cos(th1)*sin(th4) + cos(th2 + th3)*cos(th4)*sin(th1)) - sin(th2 + th3)*cos(th5)*sin(th1),
               1, 0, 0, cos(th2 + th3), sin(th2 + th3)*sin(th4), cos(th2 + th3)*cos(th5) - sin(th2 + th3)*cos(th4)*sin(th5);
    }
    
    yoon_ori_trans::Quaternion q_out;
    yoon_ori_trans::EulerAngles eul_out;
    Eigen::Matrix3f roty180;
    void Y_indy7_rot2eul(Eigen::Matrix3f& rot, Eigen::Vector3f& out_EulerAnlges)
    {
        // to make the 0 euler angle at initial pos
        roty180 << -1,0,0,
                    0,1,0,
                    0,0,-1;
        rot = roty180*rot;
        q_out = yoon_ori_trans::Rot2Quaternian(rot);
        eul_out = yoon_ori_trans::Qua2EulerAngles(q_out);

        out_EulerAnlges(0) = -eul_out.roll; // x-axis rotation angle
        out_EulerAnlges(1) = eul_out.pitch; // x-axis rotation angle
        out_EulerAnlges(2) = -eul_out.yaw; // x-axis rotation angle
    }

}

namespace yoon_ori_trans
{

    EulerAngles Qua2EulerAngles(Quaternion q)
    {
        EulerAngles angles;

        // roll (x-axis rotation)
        double sinr_cosp = 2*(q.w*q.x+q.y*q.z);
        double cosr_cosp = 1 - 2*(q.x*q.x+q.y*q.y);
        angles.roll = atan2(sinr_cosp,cosr_cosp);
        // std::cout<<sinr_cosp<<", "<< cosr_cosp <<std::endl;
        // if(isnan(angles.roll)) std::cout<<"nan1"<<std::endl;
        // pitch (y-axis rotation);
        double sinp = 2*(q.w*q.y-q.z*q.x);
        if(fabs(sinp) >= 1)
            angles.pitch = copysign(M_PI/2,sinp); // use 90 degrees if out of range
        else
            angles.pitch = asin(sinp);

        // if(isnan(angles.pitch)) std::cout<<"nan2"<<std::endl;
        // yaw (z-axis rotation)
        double siny_cosp = 2*(q.w*q.z+q.x*q.y);
        double cosy_cosp = 1 - 2*(q.y*q.y+q.z*q.z);
        angles.yaw = atan2(siny_cosp,cosy_cosp);

        // if(isnan(angles.yaw)) std::cout<<"nan3"<<std::endl;
        return angles; 


    }

    Quaternion Rot2Quaternian(Eigen::Matrix3f rot)
    {
        // transfer from rotation matrix to quaternion
        Quaternion quat;
        // std::cout << 0.5*pow(1+rot(0,0)+rot(1,1)+rot(2,2),0.5) << std::endl;
        quat.w = 0.5*pow(1+rot(0,0)+rot(1,1)+rot(2,2),0.5);
        // if(isnan(quat.w)) std::cout<<"nan0"<<std::endl;
        quat.x = (rot(2,1)-rot(1,2))/(4*quat.w);
        quat.y = (rot(0,2)-rot(2,0))/(4*quat.w);
        quat.z = (rot(1,0)-rot(0,1))/(4*quat.w);
        // printf("%f,%f,%f,%f \n",quat.w,quat.x,quat.y,quat.z);
        return quat;

    }


}

namespace yoon_indy7_inverse
{
    // ---- Parameters ---- //  
    // For inverse kinematics
    Eigen::MatrixXf Jacobian = Eigen::MatrixXf::Zero(6,6);
    Eigen::MatrixXf inv_Jaco, inv_Jaco_1;
    Eigen::Matrix4f Cur_HTM = Eigen::Matrix4f::Zero(4,4);
    Eigen::VectorXf in_car_divi = Eigen::VectorXf::Zero(6);
    Eigen::VectorXf out_angular_divi = Eigen::VectorXf::Zero(6);
    Eigen::Matrix3f Rot_T06;
    Eigen::Matrix3f Rot_y180;
    yoon_ori_trans::Quaternion rot2qua;
    yoon_ori_trans::EulerAngles qua2euler;

    // ---- Functions ---- //

    void Y_indy7_invk(Eigen::VectorXf& Tar_posture, double* Current_th, double* Out_th)
    {

        // Input : Tar_posture(mm & rad,[x,y,z,roll,pitch,yaw]), Current_th(rad,each joint), Output : Out_th(rad,each joint)
        yoon_indy7_kinematics::Y_indy7_jacob(Current_th, Jacobian); // Jacobian calculation
        inv_Jaco_1 = Jacobian*Jacobian.transpose(); // Moore-Penrose pseudo-inverse 1
        inv_Jaco = Jacobian.transpose()*inv_Jaco_1.inverse(); // Moore-Penrose pseudo-inverse 2


        yoon_indy7_kinematics::Y_indy7_forw(Current_th, Cur_HTM);

        Rot_T06 = Cur_HTM.block(0,0,3,3);
        Rot_y180 << -1,0,0,
                0,1,0,
                0,0,-1;
        Rot_T06 = Rot_y180*Rot_T06;
        rot2qua = yoon_ori_trans::Rot2Quaternian(Rot_T06); //transfer the rotation matrix to quaternion
        qua2euler = yoon_ori_trans::Qua2EulerAngles(rot2qua);       


        in_car_divi(0) = Tar_posture(0)-Cur_HTM(0,3);
        in_car_divi(1) = Tar_posture(1)-Cur_HTM(1,3);
        in_car_divi(2) = Tar_posture(2)-Cur_HTM(2,3);
        in_car_divi(3) = Tar_posture(3)-(-qua2euler.roll);
        in_car_divi(4) = Tar_posture(4)-(qua2euler.pitch);
        in_car_divi(5) = Tar_posture(5)-(-qua2euler.yaw);

        out_angular_divi = inv_Jaco*in_car_divi; // be careful of nan!!!!!!!!!
        // std::cout << out_angular_divi << std::endl;

        for(int i=0;i<6;i++)
        {
            Out_th[i] = Current_th[i] + out_angular_divi(i);

        }

    }
}


namespace yoon_indy7_gravity
{
    
    // DH-parameter
    double d1=0.3, d4=0.35, d5=0.1865, d6=0.228;
    double a2=0.45;
    
    double g = 9.81;

    // Mass
    // double m1 = 11.803; // unit : kg
    // double m2 = 7.993;
    // double m3 = 2.991;
    // double m4 = 2.123;
    // double m5 = 2.229;
    // double m6 = 0.400;

    double m1 = 0.1; // unit : kg
    double m2 = 0.1;
    double m3 = 0.1;
    double m4 = 0.1;
    double m5 = 0.1;
    double m6 = 0.1;
    // Mass center
    // double x_b[6] = {0.00000000, 0.25216149, 0.00000000, 0.00000000, 0.00000000, 0.00000000};
    // double y_b[6] = {0.00000000, 0.00000000, -0.15122139, -0.07574045, 0.09240000, 0.00000000}];
    // double z_b[6] = {0.00000000, 0.16846227, 0.00000000, 0.00000000, 0.00000000, 0.00000000};
    
    double x_b1 = 0.00000000     ,y_b1 = 0.00000000        ,z_b1 = 0.00000000; //unit: m ..
    double x_b2 = 0.25216149     ,y_b2 = 0.00000000        ,z_b2 = 0.16846227; //unit: m ..
    double x_b3 = 0.00000000     ,y_b3 = -0.15122139       ,z_b3 = 0.00000000; //unit: m ..
    double x_b4 = 0.00000000     ,y_b4 = -0.07574045       ,z_b4 = 0.00000000; //unit: m ..
    double x_b5 = 0.00000000     ,y_b5 = 0.09240000        ,z_b5 = 0.00000000; //unit: m ..
    double x_b6 = 0.00000000     ,y_b6 = 0.00000000        ,z_b6 = 0.00000000; //unit: m ..
    
    
    
    double th1,th2,th3,th4,th5,th6;

    void gravity(double* angles,int joint, double& G_torque)
    {
        th1 = angles[0];
        th2 = angles[1];
        th3 = angles[2];
        th4 = angles[3];
        th5 = angles[4];
        th6 = angles[5];

        if(joint == 1){
            G_torque = 0;}
        
        else if (joint == 2){
            G_torque = -g*(m4*z_b4*sin(th2 + th3) - m3*y_b3*sin(th2 + th3) + a2*m3*sin(th2) + a2*m4*sin(th2) + a2*m5*sin(th2) + a2*m6*sin(th2) + m2*y_b2*cos(th2) + m2*x_b2*sin(th2) + d4*m4*sin(th2 + th3) + d4*m5*sin(th2 + th3) + d4*m6*sin(th2 + th3) + m3*x_b3*cos(th2 + th3) - d5*m5*cos(th2 + th3)*sin(th4) - d5*m6*cos(th2 + th3)*sin(th4) + d6*m6*sin(th2 + th3)*cos(th5)
            - m4*x_b4*cos(th2 + th3)*cos(th4) + m4*y_b4*cos(th2 + th3)*sin(th4) + m5*y_b5*sin(th2 + th3)*cos(th5) - m5*z_b5*cos(th2 + th3)*sin(th4) + m6*z_b6*sin(th2 + th3)*cos(th5) + m5*x_b5*sin(th2 + th3)*sin(th5) + d6*m6*cos(th2 + th3)*cos(th4)*sin(th5) - m5*x_b5*cos(th2 + th3)*cos(th4)*cos(th5) + m5*y_b5*cos(th2 + th3)*cos(th4)*sin(th5) + m6*y_b6*cos(th2 + th3)*cos(th6)*sin(th4)
            + m6*z_b6*cos(th2 + th3)*cos(th4)*sin(th5) + m6*x_b6*cos(th2 + th3)*sin(th4)*sin(th6) + m6*x_b6*sin(th2 + th3)*cos(th6)*sin(th5) - m6*y_b6*sin(th2 + th3)*sin(th5)*sin(th6) - m6*x_b6*cos(th2 + th3)*cos(th4)*cos(th5)*cos(th6) + m6*y_b6*cos(th2 + th3)*cos(th4)*cos(th5)*sin(th6));}
        
        else if (joint == 3){
            G_torque = - g*m6*(x_b6*(cos(th6)*(sin(th2 + th3)*sin(th5) - cos(th2 + th3)*cos(th4)*cos(th5)) + cos(th2 + th3)*sin(th4)*sin(th6)) - y_b6*(sin(th6)*(sin(th2 + th3)*sin(th5) - cos(th2 + th3)*cos(th4)*cos(th5)) - cos(th2 + th3)*cos(th6)*sin(th4)) + d4*sin(th2 + th3) + d6*(sin(th2 + th3)*cos(th5) + cos(th2 + th3)*cos(th4)*sin(th5)) 
            + z_b6*(sin(th2 + th3)*cos(th5) + cos(th2 + th3)*cos(th4)*sin(th5)) - d5*cos(th2 + th3)*sin(th4)) - g*m4*(d4*sin(th2 + th3) + z_b4*sin(th2 + th3) - x_b4*cos(th2 + th3)*cos(th4) + y_b4*cos(th2 + th3)*sin(th4)) - g*m3*(x_b3*cos(th2 + th3) - y_b3*sin(th2 + th3)) - g*m5*(d4*sin(th2 + th3) + x_b5*(sin(th2 + th3)*sin(th5) - cos(th2 + th3)*cos(th4)*cos(th5))
            + y_b5*(sin(th2 + th3)*cos(th5) + cos(th2 + th3)*cos(th4)*sin(th5)) - d5*cos(th2 + th3)*sin(th4) - z_b5*cos(th2 + th3)*sin(th4));}
                    
        else if (joint == 4){
            G_torque = g*sin(th2 + th3)*(d5*m5*cos(th4) + d5*m6*cos(th4) - m4*y_b4*cos(th4) + m5*z_b5*cos(th4) - m4*x_b4*sin(th4) + d6*m6*sin(th4)*sin(th5) - m5*x_b5*cos(th5)*sin(th4) - m6*x_b6*cos(th4)*sin(th6) + m5*y_b5*sin(th4)*sin(th5) + m6*z_b6*sin(th4)*sin(th5) - m6*y_b6*cos(th4)*cos(th6) - m6*x_b6*cos(th5)*cos(th6)*sin(th4) + m6*y_b6*cos(th5)*sin(th4)*sin(th6));}
        
        else if (joint == 5){
            G_torque = g*m5*(x_b5*(cos(th2 + th3)*cos(th5) - sin(th2 + th3)*cos(th4)*sin(th5)) - y_b5*(cos(th2 + th3)*sin(th5) + sin(th2 + th3)*cos(th4)*cos(th5))) - g*m6*(d6*(cos(th2 + th3)*sin(th5) + sin(th2 + th3)*cos(th4)*cos(th5)) + z_b6*(cos(th2 + th3)*sin(th5) + sin(th2 + th3)*cos(th4)*cos(th5)) - x_b6*cos(th6)*(cos(th2 + th3)*cos(th5) - sin(th2 + th3)*cos(th4)*sin(th5)) 
            + y_b6*sin(th6)*(cos(th2 + th3)*cos(th5) - sin(th2 + th3)*cos(th4)*sin(th5)));}
        
        else if (joint == 6){
            G_torque = -g*m6*(x_b6*(sin(th6)*(cos(th2 + th3)*sin(th5) + sin(th2 + th3)*cos(th4)*cos(th5)) + sin(th2 + th3)*cos(th6)*sin(th4)) + y_b6*(cos(th6)*(cos(th2 + th3)*sin(th5) + sin(th2 + th3)*cos(th4)*cos(th5)) - sin(th2 + th3)*sin(th4)*sin(th6)));}
        
        else
        {G_torque = 0;}

    }


}





namespace yoon_indy7_traj
{
    void Y_indy7_movj::movj_init(const double init_val, const double final_val)
    {


        this->status = 1;
        this->const_alpha = 1;

        this->init_an = init_val;
        this->final_an = final_val;
        
        this->t_f = sqrt(4*fabs(this->final_an-this->init_an)/this->target_jacc)+const_alpha;
        this->t_b_check = fabs(4*(this->final_an-this->init_an)/pow(this->t_f,2));
        this->y_sign = (this->final_an >= this->init_an) - (this->final_an < this->init_an);


        if(this->t_b_check >= this->target_jacc)
        {
            this->t_b = (this->t_f)/2;
        }
        else
        {
            this->t_b = fabs((this->t_f)/2-sqrt(pow((this->target_jacc)*(this->t_f),2)-4*(this->target_jacc)*fabs(this->final_an-this->init_an))/(2*(this->target_jacc))); 
        }

        this->th_b = this->init_an + (this->y_sign)*(1/2)*(this->target_jacc)*pow(this->t_b,2);
        this->th_fb = (this->t_f - 2*(this->t_b))*this->y_sign*this->target_jacc*this->t_b;
    }

    double Y_indy7_movj::indy7_movj(const int counter, double& out_th)
    {
        // this->t_f = fabs(final_an-init_an)/(this->target_jvel);


        // printf("%f, %f \n",init_an,final_an);


        if(counter*this->sampling_time < this->t_b)
        {
            out_th = this->init_an + this->y_sign*(1/2)*this->target_jacc*pow(counter*this->sampling_time,2);
        }
        else if(counter*this->sampling_time < this->t_f - this->t_b)
        {
            out_th = this->th_b + (counter*this->sampling_time - this->t_b)*this->y_sign*this->target_jacc*this->t_b;
        }
        else
        {
            out_th = this->th_b+this->th_fb + (counter*this->sampling_time - (this->t_f-this->t_b))*this->y_sign*this->target_jacc*this->t_b
             - this->y_sign*(1/2)*this->target_jacc*pow(counter*this->sampling_time-(this->t_f-this->t_b),2);
        }

        if(counter*this->sampling_time >= t_f)
        {
            this->status = -1;
        }
        else
        {
            this->status = 1;
        }
        // printf("%d",status);
        return this->status;

        
    }

    void Y_indy7_movp::movp_init(Eigen::MatrixXf& command_value, double* Present_Jangle)
    {
        
        this->sampling_time = 0.001;
        this->acc_time = 0.1; // defualt : 0.5s
        this->dt = this->sampling_time;
        

        // command_value : x, y, z, alpha, beta, gamma, travel_time(second)
        // Present_pos : x, y, z, alpha, beta, gamma

        // -------- Calculate current posture start (posi ,ori: Euler expression) -------- //

        #if 0 // for the normal case but needed to be added something more
        yoon_indy7_kinematics::Y_indy7_forw(Present_Jangle, FK_T06);
        yoon_indy7_kinematics::Y_forw_Eulered(FK_T06, Vec_FK_T06);
        #endif

        // -------- Calculate current posture end (posi ,ori: Euler expression) -------- //

        std::cout << "asd0"<<std::endl;
        this->Tar_pos = Eigen::MatrixXf::Zero(command_value.rows(),6); // x, y, z, alpha, beta, gamma
        this->Tar_pos << command_value.leftCols(6);
        std::cout << "asd1"<<std::endl;
        this->cmd_travel_time = Eigen::MatrixXf::Zero(command_value.rows(),1);
        // cmd_travel_time << command_value.block(0,6,command_value.rows(),1); // second
        this->cmd_travel_time << command_value.rightCols(1); // second

        this->traj_counter = (int)this->Tar_pos.rows();
        
        for(int i=0; i<this->traj_counter; i++)
        {
            if(i == 0)
            {
                // travel_dist = sqrt(pow(Tar_pos(i,0)-Present_pos(1),2)+pow(Tar_pos(i,1)-Present_pos(0),2));
                // travel_time = abs(travel_dist/Tar_vel(i));
                this->travel_time = this->cmd_travel_time(i);
                this->travel_velX = (this->Tar_pos(i,0)-this->Vec_FK_T06(0))/this->travel_time;
                this->travel_velY = (this->Tar_pos(i,1)-this->Vec_FK_T06(1))/this->travel_time;
                this->travel_velZ = (this->Tar_pos(i,2)-this->Vec_FK_T06(2))/this->travel_time;
                this->travel_velAl = (this->Tar_pos(i,3)-this->Vec_FK_T06(3))/this->travel_time;
                this->travel_velBe = (this->Tar_pos(i,4)-this->Vec_FK_T06(4))/this->travel_time;
                this->travel_velGa = (this->Tar_pos(i,5)-this->Vec_FK_T06(5))/this->travel_time;
                
                this->pre_ini_vel_proX = Eigen::VectorXf::Constant((int)(this->travel_time/this->dt), this->travel_velX);
                this->pre_ini_vel_proY = Eigen::VectorXf::Constant((int)(this->travel_time/this->dt), this->travel_velY);
                this->pre_ini_vel_proZ = Eigen::VectorXf::Constant((int)(this->travel_time/this->dt), this->travel_velZ);
                this->pre_ini_vel_proAl = Eigen::VectorXf::Constant((int)(this->travel_time/this->dt), this->travel_velAl);
                this->pre_ini_vel_proBe = Eigen::VectorXf::Constant((int)(this->travel_time/this->dt), this->travel_velBe);
                this->pre_ini_vel_proGa = Eigen::VectorXf::Constant((int)(this->travel_time/this->dt), this->travel_velGa);

                this->ini_vel_proX = Eigen::VectorXf::Zero(this->pre_ini_vel_proX.size());
                this->ini_vel_proY = Eigen::VectorXf::Zero(this->pre_ini_vel_proY.size());
                this->ini_vel_proZ = Eigen::VectorXf::Zero(this->pre_ini_vel_proZ.size());
                this->ini_vel_proAl = Eigen::VectorXf::Zero(this->pre_ini_vel_proAl.size());
                this->ini_vel_proBe = Eigen::VectorXf::Zero(this->pre_ini_vel_proBe.size());
                this->ini_vel_proGa = Eigen::VectorXf::Zero(this->pre_ini_vel_proGa.size());
            }
            else
            {
                // travel_dist = sqrt(pow(Tar_pos(i,0)-Tar_pos(i-1,0),2)+pow(Tar_pos(i,1)-Tar_pos(i-1,1),2));
                // travel_time = abs(travel_dist/Tar_vel(i));
                this->travel_velX = (this->Tar_pos(i,0)-this->Tar_pos(i-1,0))/this->travel_time;
                this->travel_velY = (this->Tar_pos(i,1)-this->Tar_pos(i-1,1))/this->travel_time;
                this->travel_velZ = (this->Tar_pos(i,2)-this->Tar_pos(i-1,2))/this->travel_time;
                this->travel_velAl = (this->Tar_pos(i,3)-this->Tar_pos(i-1,3))/this->travel_time;
                this->travel_velBe = (this->Tar_pos(i,4)-this->Tar_pos(i-1,4))/this->travel_time;
                this->travel_velGa = (this->Tar_pos(i,5)-this->Tar_pos(i-1,5))/this->travel_time;

                this->pre_ini_vel_proX = Eigen::VectorXf::Constant((int)(this->travel_time/this->dt), this->travel_velX);
                this->pre_ini_vel_proY = Eigen::VectorXf::Constant((int)(this->travel_time/this->dt), this->travel_velY);
                this->pre_ini_vel_proZ = Eigen::VectorXf::Constant((int)(this->travel_time/this->dt), this->travel_velZ);
                this->pre_ini_vel_proAl = Eigen::VectorXf::Constant((int)(this->travel_time/this->dt), this->travel_velAl);
                this->pre_ini_vel_proBe = Eigen::VectorXf::Constant((int)(this->travel_time/this->dt), this->travel_velBe);
                this->pre_ini_vel_proGa = Eigen::VectorXf::Constant((int)(this->travel_time/this->dt), this->travel_velGa);
                                
                this->ini_vel_proX.conservativeResize(this->ini_vel_proX.size() + this->pre_ini_vel_proX.size());
                this->ini_vel_proY.conservativeResize(this->ini_vel_proY.size() + this->pre_ini_vel_proY.size());
                this->ini_vel_proZ.conservativeResize(this->ini_vel_proZ.size() + this->pre_ini_vel_proZ.size());
                this->ini_vel_proAl.conservativeResize(this->ini_vel_proAl.size() + this->pre_ini_vel_proAl.size());
                this->ini_vel_proBe.conservativeResize(this->ini_vel_proBe.size() + this->pre_ini_vel_proBe.size());
                this->ini_vel_proGa.conservativeResize(this->ini_vel_proGa.size() + this->pre_ini_vel_proGa.size());
            }
            std::cout << "hell0" << std::endl;
            this->ini_vel_proX.bottomRows(this->pre_ini_vel_proX.size()) << this->pre_ini_vel_proX;
            this->ini_vel_proY.bottomRows(this->pre_ini_vel_proY.size()) << this->pre_ini_vel_proY;
            this->ini_vel_proZ.bottomRows(this->pre_ini_vel_proZ.size()) << this->pre_ini_vel_proZ;
            this->ini_vel_proAl.bottomRows(this->pre_ini_vel_proAl.size()) << this->pre_ini_vel_proAl;
            this->ini_vel_proBe.bottomRows(this->pre_ini_vel_proBe.size()) << this->pre_ini_vel_proBe;
            this->ini_vel_proGa.bottomRows(this->pre_ini_vel_proGa.size()) << this->pre_ini_vel_proGa;
            std::cout << "hell1" << std::endl;

        }

        // Acceleration & Deceleration part (for position command)
        this->input_velX = Eigen::VectorXf::Zero(2*(int)(this->acc_time/this->dt)+this->ini_vel_proX.size());
        this->input_velY = Eigen::VectorXf::Zero(2*(int)(this->acc_time/this->dt)+this->ini_vel_proY.size());
        this->input_velZ = Eigen::VectorXf::Zero(2*(int)(this->acc_time/this->dt)+this->ini_vel_proZ.size());
        this->input_velAl = Eigen::VectorXf::Zero(2*(int)(this->acc_time/this->dt)+this->ini_vel_proAl.size());
        this->input_velBe = Eigen::VectorXf::Zero(2*(int)(this->acc_time/this->dt)+this->ini_vel_proBe.size());
        this->input_velGa = Eigen::VectorXf::Zero(2*(int)(this->acc_time/this->dt)+this->ini_vel_proGa.size());
        

        std::cout << "hell2-0" << std::endl;
        this->input_velX.topRows((int)(this->acc_time/this->dt)+this->ini_vel_proX.size()) << Eigen::VectorXf::Zero((int)(this->acc_time/this->dt)),this->ini_vel_proX;
        this->input_velY.topRows((int)(this->acc_time/this->dt)+this->ini_vel_proY.size()) << Eigen::VectorXf::Zero((int)(this->acc_time/this->dt)),this->ini_vel_proY;
        this->input_velZ.topRows((int)(this->acc_time/this->dt)+this->ini_vel_proZ.size()) << Eigen::VectorXf::Zero((int)(this->acc_time/this->dt)),this->ini_vel_proZ;
        this->input_velAl.topRows((int)(this->acc_time/this->dt)+this->ini_vel_proAl.size()) << Eigen::VectorXf::Zero((int)(this->acc_time/this->dt)),this->ini_vel_proAl;
        this->input_velBe.topRows((int)(this->acc_time/this->dt)+this->ini_vel_proBe.size()) << Eigen::VectorXf::Zero((int)(this->acc_time/this->dt)),this->ini_vel_proBe;
        this->input_velGa.topRows((int)(this->acc_time/this->dt)+this->ini_vel_proGa.size()) << Eigen::VectorXf::Zero((int)(this->acc_time/this->dt)),this->ini_vel_proGa;
        std::cout << "hell2-1" << std::endl;

        this->output_velX.resize(this->input_velX.size());
        this->output_velY.resize(this->input_velY.size());
        this->output_velZ.resize(this->input_velZ.size());
        this->output_velAl.resize(this->input_velAl.size());
        this->output_velBe.resize(this->input_velBe.size());
        this->output_velGa.resize(this->input_velGa.size());
        // output_velX = Eigen::VectorXf::Zero(input_velX.size());
        // output_velY = Eigen::VectorXf::Zero(input_velX.size());
        std::cout << "hell2-11" << std::endl;
        
        this->output_posX.resize(this->input_velX.size());
        this->output_posY.resize(this->input_velY.size());
        this->output_posZ.resize(this->input_velZ.size());
        this->output_posAl.resize(this->input_velAl.size());
        this->output_posBe.resize(this->input_velBe.size());
        this->output_posGa.resize(this->input_velGa.size());
        // output_posX = Eigen::VectorXf::Zero(input_velX.size());
        // output_posY = Eigen::VectorXf::Zero(input_velY.size());

        std::cout << "hell2-2" << std::endl;
        this->output_posX(0) = this->Vec_FK_T06(0);
        this->output_posY(0) = this->Vec_FK_T06(1);
        this->output_posZ(0) = this->Vec_FK_T06(2);
        this->output_posAl(0) = this->Vec_FK_T06(3);
        this->output_posBe(0) = this->Vec_FK_T06(4);
        this->output_posGa(0) = this->Vec_FK_T06(5);

        std::cout << "hell2-3" << std::endl;
        for(int i=0;i<this->input_velX.size();i++)
        {
            if(i<(int)(this->acc_time/this->dt))
            {
                this->output_velX(i) =  this->input_velX(i);
                this->output_velY(i) =  this->input_velY(i);
                this->output_velZ(i) =  this->input_velZ(i);
                this->output_velAl(i) =  this->input_velAl(i);
                this->output_velBe(i) =  this->input_velBe(i);
                this->output_velGa(i) =  this->input_velGa(i);
                
                // ROS_INFO("hello2");
            }
            else
            {
                this->output_velX(i) =  this->output_velX(i-1)+(1/(this->acc_time/this->dt))*(this->input_velX(i)-this->input_velX(i-(int)(this->acc_time/this->dt)+1));
                this->output_velY(i) =  this->output_velY(i-1)+(1/(this->acc_time/this->dt))*(this->input_velY(i)-this->input_velY(i-(int)(this->acc_time/this->dt)+1));
                this->output_velZ(i) =  this->output_velZ(i-1)+(1/(this->acc_time/this->dt))*(this->input_velZ(i)-this->input_velZ(i-(int)(this->acc_time/this->dt)+1));
                this->output_velAl(i) =  this->output_velAl(i-1)+(1/(this->acc_time/this->dt))*(this->input_velAl(i)-this->input_velAl(i-(int)(this->acc_time/this->dt)+1));
                this->output_velBe(i) =  this->output_velBe(i-1)+(1/(this->acc_time/this->dt))*(this->input_velBe(i)-this->input_velBe(i-(int)(this->acc_time/this->dt)+1));
                this->output_velGa(i) =  this->output_velGa(i-1)+(1/(this->acc_time/this->dt))*(this->input_velGa(i)-this->input_velGa(i-(int)(this->acc_time/this->dt)+1));
                // ROS_INFO("%f",output_velY(i));
                // ROS_INFO("hello3");
            }

            if(i != 0)
            {
                this->output_posX(i) = this->output_posX(i-1) + this->output_velX(i)*this->dt;
                this->output_posY(i) = this->output_posY(i-1) + this->output_velY(i)*this->dt;
                this->output_posZ(i) = this->output_posZ(i-1) + this->output_velZ(i)*this->dt;
                this->output_posAl(i) = this->output_posAl(i-1) + this->output_velAl(i)*this->dt;
                this->output_posBe(i) = this->output_posBe(i-1) + this->output_velBe(i)*this->dt;
                this->output_posGa(i) = this->output_posGa(i-1) + this->output_velGa(i)*this->dt;
            }

            // monitor_msg.data = output_posX(i);
            // monitor_pub.publish(monitor_msg);

        }
        std::cout << "hell3" << std::endl;
        this->output_value = Eigen::MatrixXf::Zero(this->output_posX.rows(),6);
        this->output_value << this->output_posX, this->output_posY, this->output_posZ, this->output_posAl, this->output_posBe, this->output_posGa;
        this->counter = 0;
        this->initialized_bool = true;
        
        
        std::cout << "output_value size" << std::endl;
        std::cout << this->output_value.rows() << std::endl;
        
    }

    double Y_indy7_movp::indy7_movp(double* out_posture)
    {
        if(this->initialized_bool == true)
        {
          if(this->counter < this->output_value.rows())
          {
            for(int i = 0; i<6; i++)
            {
                out_posture[i] = this->output_value(this->counter,i);
            }
            this->counter++;
            // std::cout << this->counter << std::endl;
            return 1;
          }
        }

        else
        {
            this->initialized_bool = false;
            return 0;
            
            std::cout << "finish" << std::endl;
        }
    
    }

}

namespace yoon_indy7_filter
{
    // void VR_Kalman_filter::Kalman_cal(double Kalman_Zin, double& Kalman_Xout)

    float KalmanFilterExample(float input, kalman_filter_data_s* kalman_data)
    {
        float P_minus[4]; /* matrix 2x2 */
        float x_minus[2]; /* vector 2x1 */
        float K_gain[2];  /* matrix 2x1 */
        float temp_help;
        
        /* Prediction Step */
        x_minus[0] = kalman_data->Phi_matrix[0]*kalman_data->x_plus[0] + kalman_data->Phi_matrix[1]*kalman_data->x_plus[1];
        x_minus[1] = kalman_data->Phi_matrix[2]*kalman_data->x_plus[0] + kalman_data->Phi_matrix[3]*kalman_data->x_plus[1];
        P_minus[0] = (kalman_data->Phi_matrix[0]*kalman_data->P_plus[0] + kalman_data->Phi_matrix[1]*kalman_data->P_plus[2])*kalman_data->Phi_matrix[0];
        P_minus[0] += (kalman_data->Phi_matrix[0]*kalman_data->P_plus[1] + kalman_data->Phi_matrix[1]*kalman_data->P_plus[3])*kalman_data->Phi_matrix[1];
        P_minus[0] += kalman_data->Q_matrix[0];
        P_minus[1] = (kalman_data->Phi_matrix[0]*kalman_data->P_plus[0] + kalman_data->Phi_matrix[1]*kalman_data->P_plus[2])*kalman_data->Phi_matrix[2];
        P_minus[1] += (kalman_data->Phi_matrix[0]*kalman_data->P_plus[1] + kalman_data->Phi_matrix[1]*kalman_data->P_plus[3])*kalman_data->Phi_matrix[3];
        P_minus[1] += kalman_data->Q_matrix[1];
        P_minus[2] = (kalman_data->Phi_matrix[2]*kalman_data->P_plus[0] + kalman_data->Phi_matrix[3]*kalman_data->P_plus[2])*kalman_data->Phi_matrix[0];
        P_minus[2] += (kalman_data->Phi_matrix[2]*kalman_data->P_plus[1] + kalman_data->Phi_matrix[3]*kalman_data->P_plus[3])*kalman_data->Phi_matrix[1];
        P_minus[2] += kalman_data->Q_matrix[2];
        P_minus[3] = (kalman_data->Phi_matrix[2]*kalman_data->P_plus[0] + kalman_data->Phi_matrix[3]*kalman_data->P_plus[2])*kalman_data->Phi_matrix[2];
        P_minus[3] += (kalman_data->Phi_matrix[2]*kalman_data->P_plus[1] + kalman_data->Phi_matrix[3]*kalman_data->P_plus[3])*kalman_data->Phi_matrix[3];
        P_minus[3] += kalman_data->Q_matrix[3];
        /* Kalman Gain */
        temp_help = (kalman_data->H_matrix[0]*P_minus[0] + kalman_data->H_matrix[1]*P_minus[2])*kalman_data->H_matrix[0];
        temp_help += (kalman_data->H_matrix[0]*P_minus[1] + kalman_data->H_matrix[1]*P_minus[3])*kalman_data->H_matrix[1];
        temp_help += kalman_data->R_matrix;
        K_gain[0] = (kalman_data->H_matrix[0]*P_minus[0] + kalman_data->H_matrix[1]*P_minus[1])/temp_help; /* temp_help shall be !=0 */
        K_gain[1] = (kalman_data->H_matrix[0]*P_minus[2] + kalman_data->H_matrix[1]*P_minus[3])/temp_help;
        /* Correction Step */
        kalman_data->P_plus[0] = (1.0 - K_gain[0]*kalman_data->H_matrix[0])*P_minus[0] - K_gain[0]*kalman_data->H_matrix[1]*P_minus[2];
        kalman_data->P_plus[1] = (1.0 - K_gain[0]*kalman_data->H_matrix[0])*P_minus[1] - K_gain[0]*kalman_data->H_matrix[1]*P_minus[3];
        kalman_data->P_plus[2] = -K_gain[1]*kalman_data->H_matrix[0]*P_minus[0] + (1.0 - K_gain[1]*kalman_data->H_matrix[1])*P_minus[2];
        kalman_data->P_plus[3] = -K_gain[1]*kalman_data->H_matrix[0]*P_minus[1] + (1.0 - K_gain[1]*kalman_data->H_matrix[1])*P_minus[3];
        kalman_data->x_plus[0] = x_minus[0] + K_gain[0]*(input - x_minus[0]);
        kalman_data->x_plus[1] = x_minus[1] + K_gain[1]*(input - x_minus[0]);
        
        return kalman_data->x_plus[0];
    }

    float KalmanFilter1D(float input, kalman_filter_1D_s* kalman_data1D)
    {
        float x_mi, p_mi, K, x, p;

        x_mi = kalman_data1D->x_pre;
        p_mi = kalman_data1D->p_pre + kalman_data1D->Q;
        
        K = p_mi/(p_mi+kalman_data1D->R);
        x = x_mi + K*(input-x_mi);
        p = (1-K)*p_mi;
        
        kalman_data1D->x_pre = x;
        kalman_data1D->p_pre = p;

        return x;
    }
    float VR_mv_filter(float input, VR_mv_filter_ins* VR_mv_par)
    {
        float output = 0;

        if(VR_mv_par->mv_num <= VR_mv_par->counter)
        {
            VR_mv_par->counter = 0;
        }

        if(VR_mv_par->mv_num > VR_mv_par->counter)
        {
            VR_mv_par->saved_data[VR_mv_par->counter] = input;
            VR_mv_par->counter++;

            for(int i = 0;i<VR_mv_par->mv_num;i++)
            {
                output += VR_mv_par->saved_data[i];
            }
            output /= (float)VR_mv_par->mv_num;
        }

        return output;

    }

    float HPF_2nd(float HPF_in, HPF_2nd_parameter* HPF_2nd_par)
    {
        double HPF_out = 0;

        double w0, T, Q;
        double a0_, a1_, a2_, b0_, b1_, b2_;
        double a1, a2, b0, b1, b2;
        double H0 = 1;
        double sum0;
        
        w0 = 2*(3.141592)*HPF_2nd_par->f_cut;
        T = HPF_2nd_par->ts; // sampling time
        Q = 1/(2*HPF_2nd_par->zeta);
        
        a0_ = 4/(T*T) + 2*w0/(Q*T) + (w0*w0);
        a1_ = -8/(T*T) + 2*(w0*w0);
        a2_ = 4/(T*T) - 2*w0/(Q*T) + (w0*w0);
        
        b0_ = 4*H0/(T*T);
        b1_ = -8*H0/(T*T);
        b2_ = 4*H0/(T*T);
        
        a1 = a1_/a0_;
        a2 = a2_/a0_;
        b0 = b0_/a0_;
        b1 = b1_/a0_;
        b2 = b2_/a0_;
        
        
        sum0 = -a1*HPF_2nd_par->timeZone[1] - a2*HPF_2nd_par->timeZone[0];
        HPF_2nd_par->timeZone[2] = HPF_in + sum0;
        HPF_out = b0*HPF_2nd_par->timeZone[2] + b1*HPF_2nd_par->timeZone[1] + b2*HPF_2nd_par->timeZone[0];

        HPF_2nd_par->timeZone[0] = HPF_2nd_par->timeZone[1];
        HPF_2nd_par->timeZone[1] = HPF_2nd_par->timeZone[2];

        return HPF_out;
    }
    float Yoon_BSF(float BSF_in, Yoon_BSF_parameter* Yoon_BSF_par)
    {
        double BSF_out = 0;

        double w0_peak, Ts, Q;
        double a0_, a1_, a2_, b0_, b1_, b2_;
        double a1, a2, b0, b1, b2;
        double H0 = 1;
        double sum0;

        Ts =  Yoon_BSF_par->ts;
        w0_peak = 2*(3.141592)*Yoon_BSF_par->f_peak;
        Q = Yoon_BSF_par->f_peak/Yoon_BSF_par->bandWidth;
        
        b0_ = H0*4/(Ts*Ts) + H0*(w0_peak*w0_peak);
        b1_ = -2*H0*4/(Ts*Ts) + 2*(w0_peak*w0_peak);
        b2_ = H0*4/(Ts*Ts) + H0*(w0_peak*w0_peak);
        
        a0_ = 4/(Ts*Ts)+2*w0_peak/(Q*Ts)+(w0_peak*w0_peak);
        a1_ = -8/(Ts*Ts)+2*(w0_peak*w0_peak);
        a2_ = 4/(Ts*Ts)-2*w0_peak/(Q*Ts)+(w0_peak*w0_peak);
        
        a1 = a1_/a0_;
        a2 = a2_/a0_;
        b0 = b0_/a0_;
        b1 = b1_/a0_;
        b2 = b2_/a0_;
        
        
        sum0 = -a1*Yoon_BSF_par->timeZone[1] - a2*Yoon_BSF_par->timeZone[0];
        Yoon_BSF_par->timeZone[2] = BSF_in + sum0;
        BSF_out = b0*Yoon_BSF_par->timeZone[2] + b1*Yoon_BSF_par->timeZone[1] + b2*Yoon_BSF_par->timeZone[0];

        Yoon_BSF_par->timeZone[0] = Yoon_BSF_par->timeZone[1];
        Yoon_BSF_par->timeZone[1] = Yoon_BSF_par->timeZone[2];

        return BSF_out;
    }

}

// -------- for FT sensor start ------- //
namespace FT_communication 
{
    void TCP_init(FT_TCP* FT_TCP_struct)
    {
        // ip, port 정의
        char ip[] = "192.168.0.223";
        int port = 4001;

        // 클라이언트 소켓 TCP/IP 프로토콜 생성
        FT_TCP_struct->clnt_sock = socket(PF_INET, SOCK_STREAM, 0);
        if(FT_TCP_struct->clnt_sock == -1) errhandle("socket() ERR!");

        // serv_sock에 bind로 주소 넣기 위한 밑작업
        memset(&FT_TCP_struct->st_serv_addr,0,sizeof(FT_TCP_struct->st_serv_addr));
        FT_TCP_struct->st_serv_addr.sin_family = AF_INET;
        FT_TCP_struct->st_serv_addr.sin_addr.s_addr = inet_addr(ip);
        FT_TCP_struct->st_serv_addr.sin_port = htons(port);

        // connect()으로 서버소켓에 연결요청
        int connret = connect(FT_TCP_struct->clnt_sock,(struct sockaddr*) &FT_TCP_struct->st_serv_addr, sizeof(FT_TCP_struct->st_serv_addr));
        if(connret == -1) errhandle("connect() ERR!");


        // int iResult = send(FT_TCP_struct->clnt_sock, FT_TCP_struct->sendbuf, sizeof(FT_TCP_struct->sendbuf), 0);

        // printf("Bytes Sent: %d\n", iResult);

        FT_TCP_struct->init_flag = true;

    }

    double inter_force[3] = {0,}; // for force value transform
    double inter_moment[3] = {0,}; // for moment value transform
    double inter_posAcc[3] = {0,}; // for linear acceleration value transform
    double inter_angAcc[3] = {0,}; // for angular acceleration value transform

    void TCP_start(FT_TCP* FT_TCP_struct)
    {
        FT_TCP_struct->readstrlen = read(FT_TCP_struct->clnt_sock, (char*)&FT_TCP_struct->recvmsg, sizeof(FT_TCP_struct->recvmsg));
        if(FT_TCP_struct->readstrlen == -1) errhandle("read() ERR!");

        // --------------- Recieved Data Decryption --------------- //
        if (FT_TCP_struct->recvmsg[4] == 0x0A) // if ID is 1 -> motor encoder
        {
            for (int i = 0; i < 2; i++)
            {
                FT_TCP_struct->Count[i] = ((int)FT_TCP_struct->recvmsg[6 + 2 * i] * 256 + (int)FT_TCP_struct->recvmsg[7 + 2 * i]);
                FT_TCP_struct->Deg[i] = FT_TCP_struct->Count[i] * 360 / (1024 *22);
            }

            for (int i= 0; i<8; i++)
            {
                FT_TCP_struct->aq_raw[i] = (int)FT_TCP_struct->recvmsg[6 + i];
            }

            // for (int i = 0; i < 3; i++)
            // {
            //     inter_force[i] = (double)((int)FT_TCP_struct->recvmsg[6 + 2 * i] * 256 + (int)FT_TCP_struct->recvmsg[7 + 2 * i]) / 100 - 300;
            //     // inter_force[i] = FT_TCP_struct->Force_val[i];
            // }

            // // FT sensor transform according to its assembly configuration
            // FT_TCP_struct->Force_val[0] = -inter_force[1] - FT_TCP_struct->init_Force[0]; // -y value -> x value
            // FT_TCP_struct->Force_val[1] = -inter_force[0] - FT_TCP_struct->init_Force[1]; // -x value -> y value
            // FT_TCP_struct->Force_val[2] = -inter_force[2] - FT_TCP_struct->init_Force[2]; // -z value -> z value

        }

        else if (FT_TCP_struct->recvmsg[4] != 0x01) // if ID is not 1 -> something else
        {
            for (int i = 0; i < 2; i++)
            {
                FT_TCP_struct->Count[i] = 1;
            }
        }


    }
    void TCP_send(FT_TCP* FT_TCP_struct)
    {
        // uint8_t Position_Ref_Bit_Calc_1 = 0;
        // uint8_t Position_Ref_Bit_Calc_2 = 0;

        // for (int i = 0; i < 2; i++)
        // {
        //     Position_Ref_Bit_Calc_1 = FT_TCP_struct->Position_Ref[i] / 256;
        //     Position_Ref_Bit_Calc_2 = FT_TCP_struct->Position_Ref[i] % 256;

        //     FT_TCP_struct->sendbuf[6 + 2 * i] = Position_Ref_Bit_Calc_1;
        //     FT_TCP_struct->sendbuf[7 + 2 * i] = Position_Ref_Bit_Calc_2;
        // }
        
        for (int i = 0; i < 5; i++)
        {
            FT_TCP_struct->sendbuf[6 + i] = 0x08;
        }
        

        int iResult = send(FT_TCP_struct->clnt_sock, FT_TCP_struct->sendbuf, sizeof(FT_TCP_struct->sendbuf), 0);
        if(iResult == -1) errhandle("read() ERR!");
        // printf("Bytes Sent: %d\n", iResult);
    }
}
void errhandle(const char *errmsg){ // for FT_communication
  fputs(errmsg, stderr);
  fputc('\n', stderr);
  exit(1);
}

// -------- for FT sensor end ------- //

// -------- for PID tunning start ------- //
namespace yoon_PID_tunning
{
    bool PID_ramp(PID_tunning_parameter* PID_tunning_par)
    {

        //                  **************
        //                 *              *
        //                *                *
        // ***************                  *******************

        if(PID_tunning_par->cmd_counter*PID_tunning_par->sampling_time < PID_tunning_par->start_time) // initial rest
        {
            PID_tunning_par->output_angle = PID_tunning_par->current_angle;
            PID_tunning_par->cmd_counter++;
            return true;
        }
        else if(PID_tunning_par->cmd_counter*PID_tunning_par->sampling_time < PID_tunning_par->start_time + PID_tunning_par->ramp_time) // ramp plus
        {
            PID_tunning_par->output_angle += (PID_tunning_par->ramp_size*PID_tunning_par->sampling_time)/PID_tunning_par->ramp_time;
            PID_tunning_par->cmd_counter++;
            return true;
        }
        else if(PID_tunning_par->cmd_counter*PID_tunning_par->sampling_time < PID_tunning_par->start_time + PID_tunning_par->ramp_time + PID_tunning_par->dwell_time) // dwell
        {
            PID_tunning_par->cmd_counter++;
            return true;
        }
        else if(PID_tunning_par->cmd_counter*PID_tunning_par->sampling_time < PID_tunning_par->start_time + 2*PID_tunning_par->ramp_time + PID_tunning_par->dwell_time) // ramp minus
        {
            PID_tunning_par->output_angle -= (PID_tunning_par->ramp_size*PID_tunning_par->sampling_time)/PID_tunning_par->ramp_time;
            PID_tunning_par->cmd_counter++;
            return true;
        }
        else
        {
            return false;
        }
        
        
    }
}

// -------- for Energy tank start ------- //


namespace Energy_tank
{
    void eg_tank(Energy_tank_par* et_par)
    {
        if(et_par->tank_energy <= et_par->tank_Ulim) et_par->whi = 1;
        else et_par->whi = 0;

        if(et_par->Md_dot <= 0) et_par->gamma = et_par->whi;
        else et_par->gamma = 1;

        et_par->dampE_sum += (et_par->whi)*et_par->Dd*et_par->current_vel*et_par->current_vel;
        et_par->massCE_sum += (et_par->gamma)*0.5*et_par->Md_dot*et_par->current_vel*et_par->current_vel;
        
        et_par->tank_energy = et_par->init_energy + et_par->dampE_sum - et_par->massCE_sum;
    }
}