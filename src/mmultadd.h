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

#ifndef _MMULTADD_H_
#define _MMULTADD_H_

#define N 9
//#define M_PI 3.14159265358979323846264338327950288

// Scale for conversion between FP and INT
const double scale_vel  = 1e-16; 
const double scale_pos  = 1e-16;

// INT datatype used
#define REB_PARTICLE_INT_TYPE int64_t

//struct reb_particle_int {
//    REB_PARTICLE_INT_TYPE x;
//    REB_PARTICLE_INT_TYPE y;
//    REB_PARTICLE_INT_TYPE z;
//    REB_PARTICLE_INT_TYPE vx;
//    REB_PARTICLE_INT_TYPE vy;
//    REB_PARTICLE_INT_TYPE vz;
//};

//struct reb_particle {
//    double x;
//    double y;
//    double z;
//    double vx;
//    double vy;
//    double vz;
//    double ax;
//    double ay;
//    double az;
//    double m;
//};

/**
 * Design principles to achieve best performance
 *
 * 1. Declare secquential access to stream data into accelerators via a hardware FIFO 
 *    interface.  Otherwise, the default RAM interface requires all data to arrive
 *    before starting HLS accelerator
 */
//#pragma SDS data access_pattern(A:SEQUENTIAL, B:SEQUENTIAL, C:SEQUENTIAL)
void janus_run(REB_PARTICLE_INT_TYPE p_int_in[6*N], REB_PARTICLE_INT_TYPE p_int[6*N], double p[10*N]);


#endif /* _MMULTADD_H_ */

