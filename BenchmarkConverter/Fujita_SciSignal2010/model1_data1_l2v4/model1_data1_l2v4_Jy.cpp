#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 


#include "parameter.h"
#include "fixed_parameter.h"
#include "observable.h"
#include "my.h"
#include "sigmay.h"

void Jy_model1_data1_l2v4(double *nllh, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my){
    switch(iy) {
        case 0:
            nllh[0] = 0.5*pow(-mobservable_pEGFR_tot + observable_pEGFR_tot, 2)/pow(sigmaobservable_pEGFR_tot, 2) + 0.5*log(2*M_PI*pow(sigmaobservable_pEGFR_tot, 2));
            break;
        case 1:
            nllh[0] = 0.5*pow(-mobservable_pAkt_tot + observable_pAkt_tot, 2)/pow(sigmaobservable_pAkt_tot, 2) + 0.5*log(2*M_PI*pow(sigmaobservable_pAkt_tot, 2));
            break;
        case 2:
            nllh[0] = 0.5*pow(-mobservable_pS6_tot + observable_pS6_tot, 2)/pow(sigmaobservable_pS6_tot, 2) + 0.5*log(2*M_PI*pow(sigmaobservable_pS6_tot, 2));
            break;
}
}