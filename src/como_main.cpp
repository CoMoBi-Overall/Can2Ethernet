// ------------ basic headers -------------//
#include <stdio.h>
#include <cmath>
#include <iostream>
#include "ros/ros.h"
#include "yoon_indy7_func.h"

// ------------ msg headers -------------//
#include "std_msgs/Float64.h"


// ------------ Xenomai headers -------------//
#include <alchemy/task.h>
#include <errno.h>
#include <sys/mman.h>
#include <signal.h>

#define ms(x) (x*1000000)
RT_TASK control_task;
RT_TASK print_task;
RT_TASK CANsensor_task;
static int run =1;

// ------------ User parameters -------------//
#define Sampling_time 0.01 // second
FT_TCP_struct como_TCP;
double ref_posi[2] = {0,};
double act_posi[2] = {0,};
double act_vel[2] = {0,};
double past_posi[2] = {0,};
double out_torq[2] = {0,};
double Kp[2] = {1, 1};
double Kd[2] = {0, 0};


//---------------------------- Task functions ---------------------------------//
void control(void *arg)
{
    RTIME now, previous; // don't touch
    rt_task_set_periodic(NULL, TM_NOW, (int)ms(Sampling_time*1000)); //1=1000hz, 10=100hz  // don't touch

    previous=rt_timer_read(); // don't touch

    FT_communication::TCP_init(&como_TCP);
    while(1)
    {
        if(!run) // don't touch
        {
            break; // don't touch
        }

        now=rt_timer_read(); // don't touch

        // rt_printf("\n\e[32;1m\t time: %d \n", now-previous);
        // ------------------ start ------------------------//

        act_posi[0] = como_TCP.Count[0];
        act_vel[0] = (act_posi[0]-past_posi[0])/Sampling_time;
        
        out_torq[0] = Kp[0]*(ref_posi[0]-act_posi[0]) - Kd[0]*act_vel[0];

        if(como_TCP.init_flag == true)
		{
            
			FT_communication::TCP_send(&como_TCP);
            // std::cout << "hell" << std::endl;
		}

        
        for(int i=0;i<2;i++) // parameter update
        {
            past_posi[i] = act_posi[i];
        }
        

        // ------------------ start ------------------------//
        previous = now; // don't touch

        
        rt_task_wait_period(NULL);  // don't touch
    }
}


void CAN_communi(void *arg)
{
    // FT_communication::TCP_init(&como_TCP);
	while(1)
	{
		if(como_TCP.init_flag == true)
		{
			FT_communication::TCP_start(&como_TCP);
		}
	}
}


void print(void *arg)
{
    RTIME now, previous; // don't touch
    rt_task_set_periodic(NULL, TM_NOW, (int)ms(100)); //1=1000hz, 10=100hz  // don't touch

    previous=rt_timer_read(); // don't touch
    while(1)
    {
        if(!run) // don't touch
        {
            break; // don't touch
        }

        now=rt_timer_read(); // don't touch
        

        
        #if 1

        // rt_printf("\n\e[32;1m\tcurrent: %ld\n", now);
        // rt_printf("\n\e[32;1m\t Encoder_Count_L: %3.5f \n", como_TCP.Count[0]);
        rt_printf("\n\e[32;1m\t Raw: %d %d %d %d %d %d %d %d \n", como_TCP.aq_raw[0],como_TCP.aq_raw[1],como_TCP.aq_raw[2],como_TCP.aq_raw[3],
        como_TCP.aq_raw[4],como_TCP.aq_raw[5],como_TCP.aq_raw[6],como_TCP.aq_raw[7]);
        // rt_printf("\n\e[32;1m\t Motor_Displacement_L: %3.5f \n", como_TCP.Deg[0]);
        // rt_printf("\n\e[32;1m\t Encoder_Count_R: %3.5f \n", como_TCP.Count[1]);
        // rt_printf("\n\e[32;1m\t Motor_Displacement_R: %3.5f \n", como_TCP.Deg[1]);

        // rt_printf("\n\e[32;1m\t Send Check: %d \n", como_TCP.sendbuf[6]);
        // rt_printf("\n\e[32;1m\t Send Size: %d \n", sizeof(como_TCP.sendbuf));

        // rt_printf("\n\e[32;1m\t Cur_X: %3.5f, Cur_Y: %3.5f, Cur_Z: %3.5f\n", Vec_FK_T06(0)*1000, Vec_FK_T06(1)*1000, Vec_FK_T06(2)*1000);
        // rt_printf("\e[32;1m\t Cur_Al: %3.3f, Cur_Be: %3.3f, Cur_Ga: %3.3f\n", Vec_FK_T06(3)*(180/3.141592), Vec_FK_T06(4)*(180/3.141592), Vec_FK_T06(5)*(180/3.141592));
        
        // rt_printf("\e[32;1m\t Com_X: %3.5f, Com_Y: %3.5f, Com_Z: %3.5f\n", movp_cmd[0], movp_cmd[1], movp_cmd[2]);
        // rt_printf("\e[32;1m\t Com_Al: %3.5f, Com_Be: %3.5f, Com_Ga: %3.5f\n", movp_cmd[3], movp_cmd[4], movp_cmd[5]);
        
        // rt_printf("\e[32;1m\t Com_an1: %3.3f, Com_an2: %3.3f, Com_an3: %3.3f\n", cmd_angle[0]*rad2dgree, cmd_angle[1]*rad2dgree, cmd_angle[2]*rad2dgree);
        // rt_printf("\e[32;1m\t Com_an4: %3.3f, Com_an5: %3.3f, Com_an6: %3.3f\n", cmd_angle[3]*rad2dgree, cmd_angle[4]*rad2dgree, cmd_angle[5]*rad2dgree);
        // // rt_printf("\n\e[32;1m\t Com_an1: %3.3f, Com_an2: %3.3f, Com_an3: %3.3f\n", mon_angle[0]*rad2dgree, mon_angle[1]*rad2dgree, mon_angle[2]*rad2dgree);
        // // rt_printf("\e[32;1m\t Com_an4: %3.3f, Com_an5: %3.3f, Com_an6: %3.3f\n", mon_angle[3]*rad2dgree, mon_angle[4]*rad2dgree, mon_angle[5]*rad2dgree);
        
        // rt_printf("\e[32;1m\t Act_an1: %3.3f, Act_an2: %3.3f, Act_an3: %3.3f\n", act_th[0]*rad2dgree, act_th[1]*rad2dgree, act_th[2]*rad2dgree);
        // rt_printf("\e[32;1m\t Act_an4: %3.3f, Act_an5: %3.3f, Act_an6: %3.3f\n", act_th[3]*rad2dgree, act_th[4]*rad2dgree, act_th[5]*rad2dgree);
               

        // rt_printf("\n");
	    #endif

        previous = now; // don't touch
        rt_task_wait_period(NULL);  // don't touch
    }
}


void signal_handler(int sig)
{
    rt_task_delete(&control_task);
    rt_task_delete(&print_task);
    rt_task_delete(&CANsensor_task);
    
    run = 0;
    exit(1);
}

//---------------------------- Main ---------------------------------//
int main(int argc, char **argv)
{
    ros::init(argc,argv,"como_main");
   

    signal(SIGTERM, signal_handler); //Termination
    signal(SIGINT, signal_handler); //Active

    mlockall(MCL_CURRENT|MCL_FUTURE);
    
    rt_task_create(&control_task, "controling", 0, 90, 0);
    rt_task_start(&control_task, &control, NULL);
    
    rt_task_create(&print_task, "printing", 0, 50, 0);
    rt_task_start(&print_task, &print, NULL);

    rt_task_create(&CANsensor_task, "FTsensing", 0, 99, 0);
    rt_task_start(&CANsensor_task, &CAN_communi, NULL);
    
    pause();

    rt_task_delete(&control_task);
    rt_task_delete(&print_task);
    rt_task_delete(&CANsensor_task);

    return 0;
}