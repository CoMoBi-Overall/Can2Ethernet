#include <stdio.h>
#include <iostream>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Dense>
#include <cmath>

// for TCP/IP
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <sys/socket.h>


using namespace std;

namespace yoon_indy7_kinematics
{
    // functions
    void Y_indy7_forw(double* th, Eigen::Matrix4f& posture);
    void Y_indy7_jacob(double* th, Eigen::MatrixXf& J06);
    void Y_indy7_rot2eul(Eigen::Matrix3f& rot, Eigen::Vector3f& out_EulerAnlges);
    

}

namespace yoon_ori_trans
{
    struct Quaternion
    {
        double w,x,y,z;
    };
    struct EulerAngles
    {
        double roll, pitch, yaw;
    };

    EulerAngles Qua2EulerAngles(Quaternion q);
    Quaternion Rot2Quaternian(Eigen::Matrix3f rot);
}
namespace yoon_indy7_inverse
{
    // ---- Functions ---- // 
    void Y_indy7_invk(Eigen::VectorXf& Tar_posture, double* Current_th, double* Out_th);
}

namespace yoon_indy7_gravity
{
    void gravity(double* angles,int joint, double& G_torque);
}

namespace yoon_indy7_traj
{
    class Y_indy7_movj
    {

        double sampling_time = 0.001;
        double const_alpha;
        double status;
        
        private:
            double t_b, t_f, t_b_check;
            double y_sign, th_b, th_fb;
            double init_an, final_an;
        public:
            double target_jvel = 0.05;
            double target_jacc = 0.05;
            void movj_init(const double init_val, const double final_val);
            double indy7_movj(const int counter, double& out_th);
    };
    class Y_indy7_movp
    {

        private:

            
        public:
            double acc_time;
            double dt;

            double travel_velX, travel_velY, travel_velZ, travel_velAl, travel_velBe, travel_velGa;
            Eigen::VectorXf pre_ini_vel_proX, pre_ini_vel_proY, pre_ini_vel_proZ, pre_ini_vel_proAl, pre_ini_vel_proBe, pre_ini_vel_proGa;
            Eigen::VectorXf ini_vel_proX, ini_vel_proY, ini_vel_proZ, ini_vel_proAl, ini_vel_proBe, ini_vel_proGa;
            Eigen::VectorXf input_velX, input_velY, input_velZ, input_velAl, input_velBe, input_velGa;
            Eigen::VectorXf output_posX, output_posY, output_posZ, output_posAl, output_posBe, output_posGa;
            Eigen::VectorXf output_velX, output_velY, output_velZ, output_velAl, output_velBe, output_velGa;
            Eigen::MatrixXf Tar_pos, cmd_travel_time;
     
            Eigen::Matrix4f FK_T06 = Eigen::Matrix4f::Zero(4,4);
            Eigen::VectorXf Vec_FK_T06 = Eigen::VectorXf::Zero(6); 
            double travel_time; 
            int traj_counter;
            int counter;
            double sampling_time;
            bool initialized_bool = false;
            Eigen::MatrixXf output_value;  
            
            void movp_init(Eigen::MatrixXf& command_value, double* Present_Jangle);
            double indy7_movp(double* out_posture);
    };
}

/* Export structures */
typedef struct kalman_filter_data
{
    /* Transition matrix: 2x2 */
    float Phi_matrix[4];
    /* Q covariance plant noise matrix: 2x2 */
    float Q_matrix[4];
    /* Sensitivity matrix: 1X2 */
    float H_matrix[2];
    /* Observation noise: R covariance matrix 1x1 */
    float R_matrix;
    /* P plus current covariance matrix 2x2: estimate error */
    float P_plus[4];
    /* x plus current state vector 2x1: value, speed */
    float x_plus[2];
} kalman_filter_data_s;

typedef struct kalman_filter_1D
{
    float x_pre, p_pre;
    float Q, R;

} kalman_filter_1D_s;

typedef struct VR_mv_filter_par
{
    float saved_data[1000]={0,};
    int mv_num;
    int counter = 0;
} VR_mv_filter_ins;

typedef struct HPF_2nd_parameter
{
	double timeZone[3] = {0,};
	double f_cut = 5; // cut-off frequency(Hz)
	double zeta = 0.7; // damping ratio
	double ts = 0.001; // sampling time(s)
	
} HPF_2nd_par;

typedef struct Yoon_BSF_parameter
{
	double timeZone[3] = {0,};
	double f_peak = 5; // stop frequency(Hz)
	double bandWidth = 3; // stop frequency width(Hz)
	double ts = 0.001; // sampling time(s)
	
} Yoon_BSF_par;

namespace yoon_indy7_filter
{
    /* Export function */
    float KalmanFilterExample(float input, kalman_filter_data_s* kalman_data);
    float KalmanFilter1D(float input, kalman_filter_1D_s* kalman_data1D);
    float VR_mv_filter(float input, VR_mv_filter_ins* VR_mv_par);
    float HPF_2nd(float HPF_in, HPF_2nd_parameter* HPF_2nd_par);
    float Yoon_BSF(float BSF_in, Yoon_BSF_parameter* Yoon_BSF_par);
}

// ---------- for TCP/IP --------------- //
#define BUF_SIZE 14
// sockaddr_in 구조체 변수 선언
// struct sockaddr_in st_serv_addr;

typedef struct FT_TCP
{
    // parameters for setting
    char sendbuf[BUF_SIZE] = { 0x04,0x00,0x00,0x00,0x41,0x04,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00}; 
    bool init_flag = false;
    // parameters that do not touch
    int clnt_sock;
    int readstrlen = 0;
    struct sockaddr_in st_serv_addr;

    double CAN_sampling = 0.001;
    double Count[2]; // 1 is the left motor encoder count and 2 is the right motor encoder count
    double Deg[2]; // is the left motor displacement and 2 is the right motor displacement

    uint8_t Position_Ref[2]; // 1 is the left motor reference position and 2 is the right motor reference position
    uint8_t Velocity_Ref[2]; // 1 is the left motor reference velocity and 2 is the right motor reference velocity

    double Force_val[3], Moment_val[3], Pos_acc_val[3], Ang_acc_val[3];
    double aq_raw[8] = {0,};
    double init_Force[3], init_Moment[3];
    double Ang_vel_val[3] = {0,};
    double Ang_pvel_val[3] = {0,};
    unsigned char recvmsg[BUF_SIZE]; // CAN to Ethernet must be set the buffer size 16


} FT_TCP_struct;

void errhandle(const char *errmsg);

namespace FT_communication 
{
    void TCP_init(FT_TCP* FT_TCP_struct);
    void TCP_start(FT_TCP* FT_TCP_struct);
    void TCP_send(FT_TCP* FT_TCP_struct);
}

// ---------- for PID_tunning --------------- //
typedef struct PID_tunning_parameter
{
    double sampling_time = 0.001; // sampling time (s)
    double start_time; // 몇초뒤에 시작할지 (s)
    double ramp_size; // ramp 크기는 얼마인지 (rad 기준)
    double ramp_time; // ramp 를 얼마만에 도달할지 (s)
    double current_angle;   // 현재각도는 얼마인지 (rad)
    double dwell_time; // 도달해서 얼마나 멈춰있을지 (s)
    double cmd_counter = 0; // command counter

    double output_angle; // output_angle(rad)
        
} PID_tunning_par;

namespace yoon_PID_tunning
{
    bool PID_ramp(PID_tunning_parameter* PID_tunning_par);
}

// ---------- for Energy tank --------------- //

typedef struct Energy_tank_par
{
    double init_Md;
    double init_Dd;
    
    double init_energy = 2; // initial energy (unit : J)
    double tank_Ulim = 5; // tank upper limit (unit : J)
    double tank_Llim = 0.1; // tank lower limit (unit : J)

    double tank_energy; // energy stored in tank (unit : J)

    double dampE_sum = 0; // damping energy summation
    double massCE_sum = 0; // mass change energy summation
    double whi; // To prevent overflow of energy
    double gamma; 

    // should be updated at the every loop
    double Md_dot; // admittance control mass derivative
    double Dd; // admittance control damping
    double current_vel; // current vel (m/s or rad/s) 


} et_par;

namespace Energy_tank
{
    void eg_tank(Energy_tank_par* et_par);
}