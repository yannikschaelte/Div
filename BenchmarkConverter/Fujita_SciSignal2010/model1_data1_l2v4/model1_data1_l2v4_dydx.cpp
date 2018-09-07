#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"

void dydx_model1_data1_l2v4(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h){
    dydx[3] = scaleFactor_pEGFR;
    dydx[6] = scaleFactor_pEGFR;
    dydx[13] = scaleFactor_pAkt;
    dydx[19] = scaleFactor_pAkt;
    dydx[23] = scaleFactor_pS6;
}