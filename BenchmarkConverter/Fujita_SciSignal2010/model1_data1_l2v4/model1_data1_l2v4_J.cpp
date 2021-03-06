#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "flux.h"
#include "dwdx.h"

void J_model1_data1_l2v4(realtype *J, const realtype t, const realtype *x, const double *p, const double *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = 0.0 - 1.0*dwdx0 - 1.0*dwdx1;
    J[8] = -1.0*dwdx13;
    J[10] = 0.0 - 1.0*dwdx2 - 1.0*dwdx3;
    J[11] = 0.0 - 1.0*dwdx4 + 1.0*dwdx5;
    J[12] = -1.0*dwdx6;
    J[17] = 1.0*dwdx14;
    J[19] = 1.0*dwdx2;
    J[20] = 1.0*dwdx4 - 1.0*dwdx5;
    J[21] = 1.0*dwdx6;
    J[28] = -1.0*dwdx2;
    J[29] = -1.0*dwdx4;
    J[30] = -1.0*dwdx6;
    J[31] = 1.0*dwdx8;
    J[38] = 1.0*dwdx5;
    J[40] = -1.0*dwdx7 - 1.0*dwdx8;
    J[41] = -1.0*dwdx9;
    J[42] = 0.0 - 1.0*dwdx10 + 1.0*dwdx11;
    J[49] = -1.0*dwdx7;
    J[50] = -1.0*dwdx9;
    J[51] = -1.0*dwdx10;
    J[52] = 1.0*dwdx12;
    J[58] = 1.0*dwdx7;
    J[59] = 1.0*dwdx9;
    J[60] = 1.0*dwdx10 - 1.0*dwdx11;
    J[69] = 1.0*dwdx11;
    J[70] = -1.0*dwdx12;
    J[72] = 1.0*dwdx0;
    J[80] = 1.0*dwdx13 - 1.0*dwdx14;
}