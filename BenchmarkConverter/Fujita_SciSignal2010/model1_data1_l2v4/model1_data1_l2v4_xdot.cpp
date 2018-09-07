#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "flux.h"

void xdot_model1_data1_l2v4(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot[0] = -1.0*w0 + 1.0*w10 - 1.0*w9;
    xdot[1] = -1.0*w1 + 1.0*w2 - 1.0*w3 + 1.0*w8;
    xdot[2] = 1.0*w1 - 1.0*w2;
    xdot[3] = -1.0*w1 + 1.0*w6;
    xdot[4] = 1.0*w2 - 1.0*w4 + 1.0*w5 - 1.0*w6;
    xdot[5] = -1.0*w4 + 1.0*w7;
    xdot[6] = 1.0*w4 - 1.0*w5;
    xdot[7] = 1.0*w5 - 1.0*w7;
    xdot[8] = 1.0*w0 - 1.0*w8;
}