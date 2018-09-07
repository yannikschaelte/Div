#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "flux.h"

void dwdx_model1_data1_l2v4(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dwdx[0] = 1.0*EGF*reaction_1_k1;
    dwdx[1] = 1.0*EGFR_turnover;
    dwdx[2] = 1.0*Akt*reaction_2_k1;
    dwdx[3] = 1.0*reaction_4_k1;
    dwdx[4] = -1.0*reaction_2_k2;
    dwdx[5] = 1.0*reaction_3_k1;
    dwdx[6] = 1.0*pEGFR*reaction_2_k1;
    dwdx[7] = 1.0*S6*reaction_5_k1;
    dwdx[8] = 1.0*reaction_7_k1;
    dwdx[9] = 1.0*pAkt*reaction_5_k1;
    dwdx[10] = -1.0*reaction_5_k2;
    dwdx[11] = 1.0*reaction_6_k1;
    dwdx[12] = 1.0*reaction_8_k1;
    dwdx[13] = -1.0*reaction_1_k2;
    dwdx[14] = 1.0*reaction_9_k1;
}