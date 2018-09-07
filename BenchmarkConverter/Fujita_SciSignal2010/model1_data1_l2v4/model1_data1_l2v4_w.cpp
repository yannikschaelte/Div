#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"

void w_model1_data1_l2v4(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h){
    w[0] = 1.0*(-reaction_1_k2*EGF_EGFR + EGF*EGFR*reaction_1_k1);
    w[1] = 1.0*(-pEGFR_Akt*reaction_2_k2 + Akt*pEGFR*reaction_2_k1);
    w[2] = 1.0*pEGFR_Akt*reaction_3_k1;
    w[3] = 1.0*pEGFR*reaction_4_k1;
    w[4] = 1.0*(-reaction_5_k2*pAkt_S6 + S6*pAkt*reaction_5_k1);
    w[5] = 1.0*reaction_6_k1*pAkt_S6;
    w[6] = 1.0*pAkt*reaction_7_k1;
    w[7] = 1.0*pS6*reaction_8_k1;
    w[8] = 1.0*reaction_9_k1*EGF_EGFR;
    w[9] = 1.0*EGFR*EGFR_turnover;
    w[10] = 68190.0*EGFR_turnover;
}