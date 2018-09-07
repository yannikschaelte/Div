#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "species.h"
#include "sensitivity.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "flux.h"
#include "dxdotdp.h"
#include "dwdx.h"
#include "JSparse.h"

void sxdot_model1_data1_l2v4(realtype *sxdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *sx, const realtype *w, const realtype *dwdx, const realtype *J, const realtype *dxdotdp){
    sxdot[0] = dxdotdp0 + J0*sx0 + J25*sx8;
    sxdot[1] = dxdotdp1 + J2*sx1 + J26*sx8 + J5*sx2 + J9*sx3;
    sxdot[2] = dxdotdp2 + J10*sx3 + J3*sx1 + J6*sx2;
    sxdot[3] = dxdotdp3 + J11*sx3 + J12*sx4 + J4*sx1 + J7*sx2;
    sxdot[4] = dxdotdp4 + J13*sx4 + J16*sx5 + J19*sx6 + J8*sx2;
    sxdot[5] = dxdotdp5 + J14*sx4 + J17*sx5 + J20*sx6 + J23*sx7;
    sxdot[6] = dxdotdp6 + J15*sx4 + J18*sx5 + J21*sx6;
    sxdot[7] = dxdotdp7 + J22*sx6 + J24*sx7;
    sxdot[8] = dxdotdp8 + J1*sx0 + J27*sx8;
}