#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "adjoint.h"
#include "flux.h"
#include "dwdp.h"

void qBdot_model1_data1_l2v4(realtype *qBdot, const int ip, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdp){
    switch(ip) {
        case 0:
            qBdot[0] = -xB0*(-1.0*dwdp0 + 1.0*dwdp1);
            break;
        case 1:
            break;
        case 2:
            break;
        case 3:
            break;
        case 4:
            qBdot[0] = 1.0*xB0*dwdp2 - 1.0*xB8*dwdp2;
            break;
        case 5:
            qBdot[0] = 1.0*xB0*dwdp3 - 1.0*xB8*dwdp3;
            break;
        case 6:
            qBdot[0] = 1.0*xB1*dwdp4 - 1.0*xB2*dwdp4 + 1.0*xB3*dwdp4;
            break;
        case 7:
            qBdot[0] = 1.0*xB1*dwdp5 - 1.0*xB2*dwdp5 + 1.0*xB3*dwdp5;
            break;
        case 8:
            qBdot[0] = -1.0*xB1*dwdp6 + 1.0*xB2*dwdp6 - 1.0*xB4*dwdp6;
            break;
        case 9:
            qBdot[0] = 1.0*xB1*dwdp7;
            break;
        case 10:
            qBdot[0] = 1.0*xB4*dwdp8 + 1.0*xB5*dwdp8 - 1.0*xB6*dwdp8;
            break;
        case 11:
            qBdot[0] = 1.0*xB4*dwdp9 + 1.0*xB5*dwdp9 - 1.0*xB6*dwdp9;
            break;
        case 12:
            qBdot[0] = -1.0*xB4*dwdp10 + 1.0*xB6*dwdp10 - 1.0*xB7*dwdp10;
            break;
        case 13:
            qBdot[0] = -1.0*xB3*dwdp11 + 1.0*xB4*dwdp11;
            break;
        case 14:
            qBdot[0] = -1.0*xB5*dwdp12 + 1.0*xB7*dwdp12;
            break;
        case 15:
            qBdot[0] = -1.0*xB1*dwdp13 + 1.0*xB8*dwdp13;
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