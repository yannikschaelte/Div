#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "flux.h"
#include "dwdp.h"

void dxdotdp_model1_data1_l2v4(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *dwdp){
    switch(ip) {
        case 0:
            dxdotdp[0] = -1.0*dwdp0 + 1.0*dwdp1;
            break;
        case 1:
            break;
        case 2:
            break;
        case 3:
            break;
        case 4:
            dxdotdp[0] = -1.0*dwdp2;
            dxdotdp[8] = 1.0*dwdp2;
            break;
        case 5:
            dxdotdp[0] = -1.0*dwdp3;
            dxdotdp[8] = 1.0*dwdp3;
            break;
        case 6:
            dxdotdp[1] = -1.0*dwdp4;
            dxdotdp[2] = 1.0*dwdp4;
            dxdotdp[3] = -1.0*dwdp4;
            break;
        case 7:
            dxdotdp[1] = -1.0*dwdp5;
            dxdotdp[2] = 1.0*dwdp5;
            dxdotdp[3] = -1.0*dwdp5;
            break;
        case 8:
            dxdotdp[1] = 1.0*dwdp6;
            dxdotdp[2] = -1.0*dwdp6;
            dxdotdp[4] = 1.0*dwdp6;
            break;
        case 9:
            dxdotdp[1] = -1.0*dwdp7;
            break;
        case 10:
            dxdotdp[4] = -1.0*dwdp8;
            dxdotdp[5] = -1.0*dwdp8;
            dxdotdp[6] = 1.0*dwdp8;
            break;
        case 11:
            dxdotdp[4] = -1.0*dwdp9;
            dxdotdp[5] = -1.0*dwdp9;
            dxdotdp[6] = 1.0*dwdp9;
            break;
        case 12:
            dxdotdp[4] = 1.0*dwdp10;
            dxdotdp[6] = -1.0*dwdp10;
            dxdotdp[7] = 1.0*dwdp10;
            break;
        case 13:
            dxdotdp[3] = 1.0*dwdp11;
            dxdotdp[4] = -1.0*dwdp11;
            break;
        case 14:
            dxdotdp[5] = 1.0*dwdp12;
            dxdotdp[7] = -1.0*dwdp12;
            break;
        case 15:
            dxdotdp[1] = 1.0*dwdp13;
            dxdotdp[8] = -1.0*dwdp13;
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