#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 

#include <sundials/sundials_sparse.h>

#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "flux.h"
#include "dwdx.h"

void JSparse_model1_data1_l2v4(SlsMat JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JSparse->indexvals[0] = 0;
    JSparse->indexvals[1] = 8;
    JSparse->indexvals[2] = 1;
    JSparse->indexvals[3] = 2;
    JSparse->indexvals[4] = 3;
    JSparse->indexvals[5] = 1;
    JSparse->indexvals[6] = 2;
    JSparse->indexvals[7] = 3;
    JSparse->indexvals[8] = 4;
    JSparse->indexvals[9] = 1;
    JSparse->indexvals[10] = 2;
    JSparse->indexvals[11] = 3;
    JSparse->indexvals[12] = 3;
    JSparse->indexvals[13] = 4;
    JSparse->indexvals[14] = 5;
    JSparse->indexvals[15] = 6;
    JSparse->indexvals[16] = 4;
    JSparse->indexvals[17] = 5;
    JSparse->indexvals[18] = 6;
    JSparse->indexvals[19] = 4;
    JSparse->indexvals[20] = 5;
    JSparse->indexvals[21] = 6;
    JSparse->indexvals[22] = 7;
    JSparse->indexvals[23] = 5;
    JSparse->indexvals[24] = 7;
    JSparse->indexvals[25] = 0;
    JSparse->indexvals[26] = 1;
    JSparse->indexvals[27] = 8;
    JSparse->indexptrs[0] = 0;
    JSparse->indexptrs[1] = 2;
    JSparse->indexptrs[2] = 5;
    JSparse->indexptrs[3] = 9;
    JSparse->indexptrs[4] = 12;
    JSparse->indexptrs[5] = 16;
    JSparse->indexptrs[6] = 19;
    JSparse->indexptrs[7] = 23;
    JSparse->indexptrs[8] = 25;
    JSparse->indexptrs[9] = 28;
    JSparse->data[0] = 0.0 - 1.0*dwdx0 - 1.0*dwdx1;
    JSparse->data[1] = 1.0*dwdx0;
    JSparse->data[2] = 0.0 - 1.0*dwdx2 - 1.0*dwdx3;
    JSparse->data[3] = 1.0*dwdx2;
    JSparse->data[4] = -1.0*dwdx2;
    JSparse->data[5] = 0.0 - 1.0*dwdx4 + 1.0*dwdx5;
    JSparse->data[6] = 1.0*dwdx4 - 1.0*dwdx5;
    JSparse->data[7] = -1.0*dwdx4;
    JSparse->data[8] = 1.0*dwdx5;
    JSparse->data[9] = -1.0*dwdx6;
    JSparse->data[10] = 1.0*dwdx6;
    JSparse->data[11] = -1.0*dwdx6;
    JSparse->data[12] = 1.0*dwdx8;
    JSparse->data[13] = -1.0*dwdx7 - 1.0*dwdx8;
    JSparse->data[14] = -1.0*dwdx7;
    JSparse->data[15] = 1.0*dwdx7;
    JSparse->data[16] = -1.0*dwdx9;
    JSparse->data[17] = -1.0*dwdx9;
    JSparse->data[18] = 1.0*dwdx9;
    JSparse->data[19] = 0.0 - 1.0*dwdx10 + 1.0*dwdx11;
    JSparse->data[20] = -1.0*dwdx10;
    JSparse->data[21] = 1.0*dwdx10 - 1.0*dwdx11;
    JSparse->data[22] = 1.0*dwdx11;
    JSparse->data[23] = 1.0*dwdx12;
    JSparse->data[24] = -1.0*dwdx12;
    JSparse->data[25] = -1.0*dwdx13;
    JSparse->data[26] = 1.0*dwdx14;
    JSparse->data[27] = 1.0*dwdx13 - 1.0*dwdx14;
}