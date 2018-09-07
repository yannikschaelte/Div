#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "flux.h"

void dwdp_model1_data1_l2v4(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dwdp[0] = 1.0*EGFR;
    dwdp[1] = 68190.0;
    dwdp[2] = 1.0*EGF*EGFR;
    dwdp[3] = -1.0*EGF_EGFR;
    dwdp[4] = 1.0*Akt*pEGFR;
    dwdp[5] = -1.0*pEGFR_Akt;
    dwdp[6] = 1.0*pEGFR_Akt;
    dwdp[7] = 1.0*pEGFR;
    dwdp[8] = 1.0*S6*pAkt;
    dwdp[9] = -1.0*pAkt_S6;
    dwdp[10] = 1.0*pAkt_S6;
    dwdp[11] = 1.0*pAkt;
    dwdp[12] = 1.0*pS6;
    dwdp[13] = 1.0*EGF_EGFR;
}