#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "vector.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "flux.h"
#include "dwdx.h"

void Jv_model1_data1_l2v4(realtype *Jv, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *v, const realtype *w, const realtype *dwdx){
    Jv[0] = -1.0*v8*dwdx13 + (0.0 - 1.0*dwdx0 - 1.0*dwdx1)*v0;
    Jv[1] = v2*(0.0 - 1.0*dwdx4 + 1.0*dwdx5) - 1.0*v3*dwdx6 + 1.0*v8*dwdx14 + (0.0 - 1.0*dwdx2 - 1.0*dwdx3)*v1;
    Jv[2] = 1.0*v1*dwdx2 + v2*(1.0*dwdx4 - 1.0*dwdx5) + 1.0*v3*dwdx6;
    Jv[3] = -1.0*v1*dwdx2 - 1.0*v2*dwdx4 - 1.0*v3*dwdx6 + 1.0*v4*dwdx8;
    Jv[4] = 1.0*v2*dwdx5 - 1.0*v5*dwdx9 + v6*(0.0 - 1.0*dwdx10 + 1.0*dwdx11) + (-1.0*dwdx7 - 1.0*dwdx8)*v4;
    Jv[5] = -1.0*v4*dwdx7 - 1.0*v5*dwdx9 - 1.0*v6*dwdx10 + 1.0*v7*dwdx12;
    Jv[6] = 1.0*v4*dwdx7 + 1.0*v5*dwdx9 + v6*(1.0*dwdx10 - 1.0*dwdx11);
    Jv[7] = 1.0*v6*dwdx11 - 1.0*v7*dwdx12;
    Jv[8] = 1.0*v0*dwdx0 + v8*(1.0*dwdx13 - 1.0*dwdx14);
}