/*
(c) Copyright 2013 - 2016 Xilinx, Inc. All rights reserved. 

This file contains confidential and proprietary information of Xilinx, Inc. and
is protected under U.S. and international copyright and other intellectual
property laws.

DISCLAIMER 
This disclaimer is not a license and does not grant any rights to the materials
distributed herewith. Except as otherwise provided in a valid license issued to
you by Xilinx, and to the maximum extent permitted by applicable law: (1) THESE
MATERIALS ARE MADE AVAILABLE "AS IS" AND WITH ALL FAULTS, AND XILINX HEREBY
DISCLAIMS ALL WARRANTIES AND CONDITIONS, EXPRESS, IMPLIED, OR STATUTORY,
INCLUDING BUT NOT LIMITED TO WARRANTIES OF MERCHANTABILITY, NON-INFRINGEMENT, OR
FITNESS FOR ANY PARTICULAR PURPOSE; and (2) Xilinx shall not be liable (whether
in contract or tort, including negligence, or under any other theory of
liability) for any loss or damage of any kind or nature related to, arising
under or in connection with these materials, including for any direct, or any
indirect, special, incidental, or consequential loss or damage (including loss
of data, profits, goodwill, or any type of loss or damage suffered as a result
of any action brought by a third party) even if such damage or loss was
reasonably foreseeable or Xilinx had been advised of the possibility of the
same.

CRITICAL APPLICATIONS
Xilinx products are not designed or intended to be fail-safe, or for use in any
application requiring fail-safe performance, such as life-support or safety
devices or systems, Class III medical devices, nuclear facilities, applications
related to the deployment of airbags, or any other applications that could lead
to death, personal injury, or severe property or environmental damage
(individually and collectively, "Critical Applications"). Customer assumes the
sole risk and liability of any use of Xilinx products in Critical Applications,
subject only to applicable laws and regulations governing limitations on product
liability.

THIS COPYRIGHT NOTICE AND DISCLAIMER MUST BE RETAINED AS PART OF THIS FILE AT
ALL TIMES. 
*/

#include <iostream>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "sds_lib.h"
#include "mmultadd.h"

#ifndef NUM_TESTS
#define NUM_TESTS 4
#endif





class perf_counter
{
public:
     uint64_t tot, cnt, calls;
     perf_counter() : tot(0), cnt(0), calls(0) {};
     inline void reset() { tot = cnt = calls = 0; }
     inline void start() { cnt = sds_clock_counter(); calls++; };
     inline void stop() { tot += (sds_clock_counter() - cnt); };
     inline uint64_t avg_cpu_cycles() { return ((tot+(calls>>1)) / calls); };
};

void janus_run_golden(REB_PARTICLE_INT_TYPE *p_int_in, REB_PARTICLE_INT_TYPE *p_int, double* p)
{
    for(unsigned int i=0; i<6*N; i++){
        p_int[i] = p_int_in[i];
    }
    double t = 0;
    while (t<2.*M_PI*1e3){ // 1000 year
        const double dt = 0.01;
        // One leapfrog step
        // drift();
        for(unsigned int i=0; i<N; i++){
            p_int[6*i+0] += (REB_PARTICLE_INT_TYPE)(dt/2.*(double)p_int[6*i+3]*scale_vel/scale_pos);
            p_int[6*i+1] += (REB_PARTICLE_INT_TYPE)(dt/2.*(double)p_int[6*i+4]*scale_vel/scale_pos);
            p_int[6*i+2] += (REB_PARTICLE_INT_TYPE)(dt/2.*(double)p_int[6*i+5]*scale_vel/scale_pos);
        }

        // to_double_pos();
        for(unsigned int i=0; i<N; i++){
            p[10*i+0] = ((double)p_int[10*i+0])*scale_pos;
            p[10*i+1] = ((double)p_int[10*i+1])*scale_pos;
            p[10*i+2] = ((double)p_int[10*i+2])*scale_pos;
        }
        // gravity();
        for(unsigned int i=0; i<N; i++){
            p[10*i+6] = 0.;
            p[10*i+7] = 0.;
            p[10*i+8] = 0.;
            for(unsigned int j=0; j<N; j++){
                if (i!=j){
                    const double dx = p[10*i+0] - p[10*j+0];
                    const double dy = p[10*i+1] - p[10*j+1];
                    const double dz = p[10*i+2] - p[10*j+2];
                    const double _r = sqrt(dx*dx + dy*dy + dz*dz);
                    const double prefact = -1/(_r*_r*_r)*p[10*j+9];

                    p[10*i+6]    += prefact*dx;
                    p[10*i+7]    += prefact*dy;
                    p[10*i+8]    += prefact*dz;
                }
            }
        }

        // kick();
        for(unsigned int i=0; i<N; i++){
            p_int[6*i+3] += (REB_PARTICLE_INT_TYPE)(dt*p[10*i+6]/scale_vel);
            p_int[6*i+4] += (REB_PARTICLE_INT_TYPE)(dt*p[10*i+7]/scale_vel);
            p_int[6*i+5] += (REB_PARTICLE_INT_TYPE)(dt*p[10*i+8]/scale_vel);
        }

        // drift();
        for(unsigned int i=0; i<N; i++){
            p_int[6*i+0] += (REB_PARTICLE_INT_TYPE)(dt/2.*(double)p_int[6*i+3]*scale_vel/scale_pos);
            p_int[6*i+1] += (REB_PARTICLE_INT_TYPE)(dt/2.*(double)p_int[6*i+4]*scale_vel/scale_pos);
            p_int[6*i+2] += (REB_PARTICLE_INT_TYPE)(dt/2.*(double)p_int[6*i+5]*scale_vel/scale_pos);
        }

        t += dt;
    }


}

//static int result_check(float *D, float *D_sw)
//{
//     for (int i = 0; i < N * N; i++) {
//          if (D_sw[i] != D[i]) {
//               std::cout << "Mismatch: data index=" << i << "d=" << D_sw[i] 
//                        << ", dout=" << D[i] << std::endl;
//               return 1;
//          }
//     }
//     return 0;
//}

int janus_test(REB_PARTICLE_INT_TYPE* p_int, REB_PARTICLE_INT_TYPE* p_int_out_sw, REB_PARTICLE_INT_TYPE* p_int_out_hw)
{
     std::cout << "Testing " << NUM_TESTS << " iterations of " << N << "-body JANUS ..." << std::endl;

     perf_counter hw_ctr, sw_ctr;
	double p[] = {
			 0.0021709922250528,  0.0057845061154043,  -0.0001290326677066,
			-0.0003084904334499,  0.0003164862379414, 0.0000072860648107,
            0.,0.,0.,
			1.0000000000000000,
			 -0.1529074277548495,  -0.4329649810759809,  -0.0217536815870956,
			1.2130892048755062,  -0.4636664872138580, -0.1492230266727991,
            0.,0.,0.,
			0.0000001660114153,
			 -0.7051385792282048,  0.1305062392874893,  0.0423980407616931,
			-0.2107118903711193,  -1.1628741220859935, -0.0038067721592922,
            0.,0.,0.,
			0.0000024478382878,
			 0.8303864923760965,  0.5551748865431479,  -0.0001556226179998,
			-0.5694403294004744,  0.8300359440285254, -0.0000250486216637,
            0.,0.,0.,
			0.0000030404326480,
			 -1.6007632981663540,  0.4507843866326728,  0.0485350310380760,
			-0.1874661855400607,  -0.7140231189065021, -0.0103688562255236,
            0.,0.,0.,
			0.0000003227156038,
			 -4.5444724195553627,  -2.9811209359531872,  0.1140115745580475,
			0.2354668506120313,  -0.3459544002171689, -0.0038305410200901,
            0.,0.,0.,
			0.0009547919152112,
			 -0.2998316596246585,  -10.0512228718170959,  0.1866942196718307,
			0.3063599906570191,  -0.0107135147677418, -0.0120072161180579,
            0.,0.,0.,
			0.0002858856727222,
			 17.8418531053445939,  8.8433796310403689,  -0.1982994964737093,
			-0.1032131635550300,  0.1941992816066720, 0.0020584917278455,
            0.,0.,0.,
			0.0000436624373583,
			 28.6228992820092181,  -8.7910334836014847,  -0.4786090163574258,
			0.0523633993793736,  0.1755278382196959, -0.0048214129381180,
            0.,0.,0.,
			0.0000515138377263
	};
	for(unsigned int i=0; i<N; i++){
            p_int[6*i+0] = p[10*i+0]/scale_pos;
            p_int[6*i+1] = p[10*i+1]/scale_pos;
            p_int[6*i+2] = p[10*i+2]/scale_pos;
            p_int[6*i+3] = p[10*i+3]/scale_vel;
            p_int[6*i+4] = p[10*i+4]/scale_vel;
            p_int[6*i+5] = p[10*i+5]/scale_vel;
    }

     
     for (int i = 0; i < NUM_TESTS; i++) 
     {
          sw_ctr.start();
          janus_run_golden(p_int,p_int_out_sw,p);
          sw_ctr.stop();
          
          hw_ctr.start();
          janus_run(p_int,p_int_out_hw,p);
          hw_ctr.stop();
     }

     uint64_t sw_cycles = sw_ctr.avg_cpu_cycles();
     uint64_t hw_cycles = hw_ctr.avg_cpu_cycles();
     double speedup = (double) sw_cycles / (double) hw_cycles;

     std::cout << "Average number of CPU cycles running in software: "
               << sw_cycles << std::endl;
     std::cout << "Average number of CPU cycles running in hardware: "
               << hw_cycles << std::endl;
     std::cout << "Speed up: " << speedup << std::endl;

     return 0;
}

int main(int argc, char* argv[]){
     int test_passed = 0;
     REB_PARTICLE_INT_TYPE* p_int = NULL;
     REB_PARTICLE_INT_TYPE* p_int_out_sw = NULL;
     REB_PARTICLE_INT_TYPE* p_int_out_hw = NULL;
     p_int = (REB_PARTICLE_INT_TYPE*) sds_alloc(6*sizeof(REB_PARTICLE_INT_TYPE)*N);
     p_int_out_hw = (REB_PARTICLE_INT_TYPE*) sds_alloc(6*sizeof(REB_PARTICLE_INT_TYPE)*N);
     p_int_out_sw = (REB_PARTICLE_INT_TYPE*) malloc(6*sizeof(REB_PARTICLE_INT_TYPE)*N);

     test_passed = janus_test(p_int, p_int_out_sw, p_int_out_hw);
     
     std::cout << "TEST " << (test_passed ? "FAILED" : "PASSED") << std::endl;

     sds_free(p_int);
     sds_free(p_int_out_hw);
     free(p_int_out_sw);

     return (test_passed ? -1 : 0);
}

