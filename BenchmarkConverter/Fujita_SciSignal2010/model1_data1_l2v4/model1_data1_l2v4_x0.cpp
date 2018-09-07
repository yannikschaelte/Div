#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "parameter.h"
#include "fixed_parameter.h"

void x0_model1_data1_l2v4(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = init_EGFR;
    x0[3] = init_AKT;
    x0[5] = init_S6;
}