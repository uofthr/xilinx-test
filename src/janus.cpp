#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "janus.h"

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
void janus_run(REB_PARTICLE_INT_TYPE p_int_in[6*N], REB_PARTICLE_INT_TYPE p_int_out[6*N], double p_mass[N],long steps,double dt){
    REB_PARTICLE_INT_TYPE p_int[6*N];
    double p[6*N];  // x y z ax ay az
    for(unsigned int i=0; i<6*N; i++){
        p_int[i] = p_int_in[i];
    }
    for(long i=0;i<steps;i++){
        // One leapfrog step
        // drift();
        for(unsigned int i=0; i<N; i++){
            p_int[6*i+0] += (REB_PARTICLE_INT_TYPE)(dt/2.*(double)p_int[6*i+3]*scale_vel/scale_pos);
            p_int[6*i+1] += (REB_PARTICLE_INT_TYPE)(dt/2.*(double)p_int[6*i+4]*scale_vel/scale_pos);
            p_int[6*i+2] += (REB_PARTICLE_INT_TYPE)(dt/2.*(double)p_int[6*i+5]*scale_vel/scale_pos);
        }

        // to_double_pos();
        for(unsigned int i=0; i<N; i++){
            p[6*i+0] = ((double)p_int[6*i+0])*scale_pos;
            p[6*i+1] = ((double)p_int[6*i+1])*scale_pos;
            p[6*i+2] = ((double)p_int[6*i+2])*scale_pos;
        }
        // gravity();
        for(unsigned int i=0; i<N; i++){
            p[6*i+3] = 0.;
            p[6*i+4] = 0.;
            p[6*i+5] = 0.;
            for(unsigned int j=0; j<N; j++){
                if (i!=j){
                    const double dx = p[6*i+0] - p[6*j+0];
                    const double dy = p[6*i+1] - p[6*j+1];
                    const double dz = p[6*i+2] - p[6*j+2];
                    const double _r = sqrt(dx*dx + dy*dy + dz*dz);
                    const double prefact = -1/(_r*_r*_r)*p_mass[j];

                    p[6*i+3]    += prefact*dx;
                    p[6*i+4]    += prefact*dy;
                    p[6*i+5]    += prefact*dz;
                }
            }
        }

        // kick();
        for(unsigned int i=0; i<N; i++){
            p_int[6*i+3] += (REB_PARTICLE_INT_TYPE)(dt*p[6*i+3]/scale_vel);
            p_int[6*i+4] += (REB_PARTICLE_INT_TYPE)(dt*p[6*i+4]/scale_vel);
            p_int[6*i+5] += (REB_PARTICLE_INT_TYPE)(dt*p[6*i+5]/scale_vel);
        }

        // drift();
        for(unsigned int i=0; i<N; i++){
            p_int[6*i+0] += (REB_PARTICLE_INT_TYPE)(dt/2.*(double)p_int[6*i+3]*scale_vel/scale_pos);
            p_int[6*i+1] += (REB_PARTICLE_INT_TYPE)(dt/2.*(double)p_int[6*i+4]*scale_vel/scale_pos);
            p_int[6*i+2] += (REB_PARTICLE_INT_TYPE)(dt/2.*(double)p_int[6*i+5]*scale_vel/scale_pos);
        }
    }
    for(unsigned int i=0; i<6*N; i++){
        p_int_out[i] = p_int[i];
    }

}
