#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "janus.h"

void janus_run(REB_PARTICLE_INT_TYPE p_int_in[6*N], REB_PARTICLE_INT_TYPE p_int_out[6*N], double p_mass[N],long steps,double dt,double dt12){
    REB_PARTICLE_INT_TYPE p_int[6*N];
    double p[3*N];  // x y z 
    janus_loop:for(unsigned int i=0; i<6*N; i++){
        p_int[i] = p_int_in[i];
    }
    for(long i=0;i<steps;i++){
        // One leapfrog step
        // drift();
        for(unsigned int i=0; i<N; i++){
            //#pragma HLS unroll
            p_int[6*i+0] += (REB_PARTICLE_INT_TYPE)(dt12*(double)p_int[6*i+3]);
            p_int[6*i+1] += (REB_PARTICLE_INT_TYPE)(dt12*(double)p_int[6*i+4]);
            p_int[6*i+2] += (REB_PARTICLE_INT_TYPE)(dt12*(double)p_int[6*i+5]);
        }

        // to_double_pos();
        for(unsigned int i=0; i<N; i++){
            //#pragma HLS unroll
            p[3*i+0] = ((double)p_int[3*i+0])*scale_pos;
            p[3*i+1] = ((double)p_int[3*i+1])*scale_pos;
            p[3*i+2] = ((double)p_int[3*i+2])*scale_pos;
        }
        // gravity();
        for(unsigned int i=0; i<N; i++){
            //#pragma HLS unroll 
            double ax = 0.;
            double ay = 0.;
            double az = 0.;
            for(unsigned int j=0; j<N; j++){
                if (i!=j){
                    const double dx = p[3*j+0] - p[3*i+0];
                    const double dy = p[3*j+1] - p[3*i+1];
                    const double dz = p[3*j+2] - p[3*i+2];
                    const double pre_sqrt = dx*dx + dy*dy + dz*dz;
                    const double _r = sqrt(pre_sqrt);
#pragma HLS RESOURCE variable=_r core=DSqrt
                    const double pre_recip = 1./(pre_sqrt * _r);
#pragma HLS RESOURCE variable=pre_recip core=DRecip
                    const double prefact = -1.*pre_recip*p_mass[j];
#pragma HLS RESOURCE variable=prefact core=dmul

                    ax  += prefact*dx;
                    ay  += prefact*dy;
                    az  += prefact*dz;
                }
            }
            // kick();
            p_int[6*i+3] += (REB_PARTICLE_INT_TYPE)(dt*ax/scale_vel);
            p_int[6*i+4] += (REB_PARTICLE_INT_TYPE)(dt*ay/scale_vel);
            p_int[6*i+5] += (REB_PARTICLE_INT_TYPE)(dt*az/scale_vel);
        }

        // drift();
        for(unsigned int i=0; i<N; i++){
            //#pragma HLS unroll
            p_int[6*i+0] += (REB_PARTICLE_INT_TYPE)(dt12*(double)p_int[6*i+3]);
            p_int[6*i+1] += (REB_PARTICLE_INT_TYPE)(dt12*(double)p_int[6*i+4]);
            p_int[6*i+2] += (REB_PARTICLE_INT_TYPE)(dt12*(double)p_int[6*i+5]);
        }
    }
    for(unsigned int i=0; i<6*N; i++){
        p_int_out[i] = p_int[i];
    }

}
