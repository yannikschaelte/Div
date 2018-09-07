#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "parameter.h"
#include "fixed_parameter.h"
#include "observable.h"
#include "my.h"
#include "sigmay.h"

void dJydsigma_model1_data1_l2v4(double *dJydsigma, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = -1.0*pow(-mobservable_pEGFR_tot + observable_pEGFR_tot, 2)/pow(sigmaobservable_pEGFR_tot, 3) + 1.0*pow(sigmaobservable_pEGFR_tot, -1);
            break;
        case 1:
            dJydsigma[1] = -1.0*pow(-mobservable_pAkt_tot + observable_pAkt_tot, 2)/pow(sigmaobservable_pAkt_tot, 3) + 1.0*pow(sigmaobservable_pAkt_tot, -1);
            break;
        case 2:
            dJydsigma[2] = -1.0*pow(-mobservable_pS6_tot + observable_pS6_tot, 2)/pow(sigmaobservable_pS6_tot, 3) + 1.0*pow(sigmaobservable_pS6_tot, -1);
            break;
}
}