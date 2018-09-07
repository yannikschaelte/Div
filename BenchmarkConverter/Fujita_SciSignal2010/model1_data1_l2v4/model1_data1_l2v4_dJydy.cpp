#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "parameter.h"
#include "fixed_parameter.h"
#include "observable.h"
#include "my.h"
#include "sigmay.h"

void dJydy_model1_data1_l2v4(double *dJydy, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my){
    switch(iy) {
        case 0:
            dJydy[0] = 1.0*(-mobservable_pEGFR_tot + observable_pEGFR_tot)/pow(sigmaobservable_pEGFR_tot, 2);
            break;
        case 1:
            dJydy[1] = 1.0*(-mobservable_pAkt_tot + observable_pAkt_tot)/pow(sigmaobservable_pAkt_tot, 2);
            break;
        case 2:
            dJydy[2] = 1.0*(-mobservable_pS6_tot + observable_pS6_tot)/pow(sigmaobservable_pS6_tot, 2);
            break;
}
}