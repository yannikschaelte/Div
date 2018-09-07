#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"

void dydp_model1_data1_l2v4(double *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip){
    switch(ip) {
        case 0:
            break;
        case 1:
            break;
        case 2:
            break;
        case 3:
            break;
        case 4:
            break;
        case 5:
            break;
        case 6:
            break;
        case 7:
            break;
        case 8:
            break;
        case 9:
            break;
        case 10:
            break;
        case 11:
            break;
        case 12:
            break;
        case 13:
            break;
        case 14:
            break;
        case 15:
            break;
        case 16:
            dydp[1] = pAkt + pAkt_S6;
            break;
        case 17:
            dydp[0] = pEGFR + pEGFR_Akt;
            break;
        case 18:
            dydp[2] = pS6;
            break;
        case 19:
            break;
        case 20:
            break;
        case 21:
            break;
}
}