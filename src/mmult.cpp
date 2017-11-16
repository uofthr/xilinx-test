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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mmultadd.h"

/**
 *
 * Design principles to achieve II = 1
 * 1. Stream data into local RAM for inputs (multiple access required)
 * 2. Partition local RAMs into N/2 sub-arrays for fully parallel access (dual-port read)
 * 3. Pipeline the dot-product loop, to fully unroll it
 * 4. Separate multiply-accumulate in inner loop to force two FP operators
 *
 */
//void mmult (float A[N*N], float B[N*N], float C[N*N]) 
//{
//     float Abuf[N][N], Bbuf[N][N];
//#pragma HLS array_partition variable=Abuf block factor=16 dim=2
//#pragma HLS array_partition variable=Bbuf block factor=16 dim=1
//     
//     for(int i=0; i<N; i++) {
//          for(int j=0; j<N; j++) {
//#pragma HLS PIPELINE
//               Abuf[i][j] = A[i * N + j];
//               Bbuf[i][j] = B[i * N + j];
//          }
//     }
//     
//     for (int i = 0; i < N; i++) {
//          for (int j = 0; j < N; j++) {
//#pragma HLS PIPELINE
//               float result = 0;
//               for (int k = 0; k < N; k++) {
//                    float term = Abuf[i][k] * Bbuf[k][j];
//                    result += term;
//               }
//               C[i * N + j] = result;
//          }
//     }
//}
//
void janus_run(REB_PARTICLE_INT_TYPE p_int_in[6*N],REB_PARTICLE_INT_TYPE p_int_out[6*N], double p_mass[10*N])
{
    REB_PARTICLE_INT_TYPE p_int[6*N];
    double p[10*N];
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
            p[10*i+0] = ((double)p_int[6*i+0])*scale_pos;
            p[10*i+1] = ((double)p_int[6*i+1])*scale_pos;
            p[10*i+2] = ((double)p_int[6*i+2])*scale_pos;
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
                    const double prefact = -1/(_r*_r*_r)*p_mass[10*j+9];

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
    for(unsigned int i=0; i<6*N; i++){
        p_int_out[i] = p_int[i];
    }

}
