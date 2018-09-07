#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"

void y_model1_data1_l2v4(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h){
    y[0] = (pEGFR + pEGFR_Akt)*scaleFactor_pEGFR;
    y[1] = scaleFactor_pAkt*(pAkt + pAkt_S6);
    y[2] = pS6*scaleFactor_pS6;
}