#ifndef _JANUS_H_
#define _JANUS_H_

// Number of particles
#define N 9

// Scale for conversion between FP and INT
const double scale_vel  = 1e-16; 
const double scale_pos  = 1e-16;
const double scale_vel_pos = 1.;

// INT datatype used
#define REB_PARTICLE_INT_TYPE int64_t

//#pragma SDS data access_pattern(A:SEQUENTIAL, B:SEQUENTIAL, C:SEQUENTIAL)
void janus_run(REB_PARTICLE_INT_TYPE p_int_in[6*N], REB_PARTICLE_INT_TYPE p_int_out[6*N], double p_mass[N],long steps, double dt, double dt12);


#endif /* _JANUS_H_ */


