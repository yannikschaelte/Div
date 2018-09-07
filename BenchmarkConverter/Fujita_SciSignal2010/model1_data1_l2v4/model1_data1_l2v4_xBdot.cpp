#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "adjoint.h"
#include "flux.h"
#include "dwdx.h"

void xBdot_model1_data1_l2v4(realtype *xBdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    xBdot[0] = -1.0*xB8*dwdx0 - (0.0 - 1.0*dwdx0 - 1.0*dwdx1)*xB0;
    xBdot[1] = -1.0*xB2*dwdx2 + 1.0*xB3*dwdx2 - (0.0 - 1.0*dwdx2 - 1.0*dwdx3)*xB1;
    xBdot[2] = -xB1*(0.0 - 1.0*dwdx4 + 1.0*dwdx5) - xB2*(1.0*dwdx4 - 1.0*dwdx5) + 1.0*xB3*dwdx4 - 1.0*xB4*dwdx5;
    xBdot[3] = 1.0*xB1*dwdx6 - 1.0*xB2*dwdx6 + 1.0*xB3*dwdx6;
    xBdot[4] = -1.0*xB3*dwdx8 + 1.0*xB5*dwdx7 - 1.0*xB6*dwdx7 - (-1.0*dwdx7 - 1.0*dwdx8)*xB4;
    xBdot[5] = 1.0*xB4*dwdx9 + 1.0*xB5*dwdx9 - 1.0*xB6*dwdx9;
    xBdot[6] = -xB4*(0.0 - 1.0*dwdx10 + 1.0*dwdx11) + 1.0*xB5*dwdx10 - xB6*(1.0*dwdx10 - 1.0*dwdx11) - 1.0*xB7*dwdx11;
    xBdot[7] = -1.0*xB5*dwdx12 + 1.0*xB7*dwdx12;
    xBdot[8] = 1.0*xB0*dwdx13 - 1.0*xB1*dwdx14 - xB8*(1.0*dwdx13 - 1.0*dwdx14);
}