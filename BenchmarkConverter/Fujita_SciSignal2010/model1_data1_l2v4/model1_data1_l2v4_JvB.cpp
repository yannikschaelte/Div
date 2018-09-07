#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "vectorB.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "adjoint.h"
#include "flux.h"
#include "dwdx.h"

void JvB_model1_data1_l2v4(realtype *JvB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *vB, const realtype *w, const realtype *dwdx){
    JvB[0] = 1.0*vB8*dwdx0 + (0.0 - 1.0*dwdx0 - 1.0*dwdx1)*vB0;
    JvB[1] = 1.0*vB2*dwdx2 - 1.0*vB3*dwdx2 + (0.0 - 1.0*dwdx2 - 1.0*dwdx3)*vB1;
    JvB[2] = vB1*(0.0 - 1.0*dwdx4 + 1.0*dwdx5) + vB2*(1.0*dwdx4 - 1.0*dwdx5) - 1.0*vB3*dwdx4 + 1.0*vB4*dwdx5;
    JvB[3] = -1.0*vB1*dwdx6 + 1.0*vB2*dwdx6 - 1.0*vB3*dwdx6;
    JvB[4] = 1.0*vB3*dwdx8 - 1.0*vB5*dwdx7 + 1.0*vB6*dwdx7 + (-1.0*dwdx7 - 1.0*dwdx8)*vB4;
    JvB[5] = -1.0*vB4*dwdx9 - 1.0*vB5*dwdx9 + 1.0*vB6*dwdx9;
    JvB[6] = vB4*(0.0 - 1.0*dwdx10 + 1.0*dwdx11) - 1.0*vB5*dwdx10 + vB6*(1.0*dwdx10 - 1.0*dwdx11) + 1.0*vB7*dwdx11;
    JvB[7] = 1.0*vB5*dwdx12 - 1.0*vB7*dwdx12;
    JvB[8] = -1.0*vB0*dwdx13 + 1.0*vB1*dwdx14 + vB8*(1.0*dwdx13 - 1.0*dwdx14);
}