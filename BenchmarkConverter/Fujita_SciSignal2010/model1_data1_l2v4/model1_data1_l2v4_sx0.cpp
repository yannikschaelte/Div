#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"

void sx0_model1_data1_l2v4(realtype *sx0, const realtype t,const realtype *x0, const realtype *p, const realtype *k, const int ip){
    switch(ip) {
        case 0:
            break;
        case 1:
            sx0[3] = 1;
            break;
        case 2:
            sx0[0] = 1;
            break;
        case 3:
            sx0[5] = 1;
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
            break;
        case 17:
            break;
        case 18:
            break;
        case 19:
            break;
        case 20:
            break;
        case 21:
            break;
}
}