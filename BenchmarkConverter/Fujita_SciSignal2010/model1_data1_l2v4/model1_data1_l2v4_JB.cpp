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

void JB_model1_data1_l2v4(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 0.0 - 1.0*dwdx0 - 1.0*dwdx1;
    JB[8] = 1.0*dwdx0;
    JB[10] = 0.0 - 1.0*dwdx2 - 1.0*dwdx3;
    JB[11] = 1.0*dwdx2;
    JB[12] = -1.0*dwdx2;
    JB[19] = 0.0 - 1.0*dwdx4 + 1.0*dwdx5;
    JB[20] = 1.0*dwdx4 - 1.0*dwdx5;
    JB[21] = -1.0*dwdx4;
    JB[22] = 1.0*dwdx5;
    JB[28] = -1.0*dwdx6;
    JB[29] = 1.0*dwdx6;
    JB[30] = -1.0*dwdx6;
    JB[39] = 1.0*dwdx8;
    JB[40] = -1.0*dwdx7 - 1.0*dwdx8;
    JB[41] = -1.0*dwdx7;
    JB[42] = 1.0*dwdx7;
    JB[49] = -1.0*dwdx9;
    JB[50] = -1.0*dwdx9;
    JB[51] = 1.0*dwdx9;
    JB[58] = 0.0 - 1.0*dwdx10 + 1.0*dwdx11;
    JB[59] = -1.0*dwdx10;
    JB[60] = 1.0*dwdx10 - 1.0*dwdx11;
    JB[61] = 1.0*dwdx11;
    JB[68] = 1.0*dwdx12;
    JB[70] = -1.0*dwdx12;
    JB[72] = -1.0*dwdx13;
    JB[73] = 1.0*dwdx14;
    JB[80] = 1.0*dwdx13 - 1.0*dwdx14;
}