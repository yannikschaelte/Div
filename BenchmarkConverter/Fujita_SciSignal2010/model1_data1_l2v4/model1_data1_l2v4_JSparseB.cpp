#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
using amici::realtype;
#include <cmath> 

#include <sundials/sundials_sparse.h>

#include "species.h"
#include "parameter.h"
#include "fixed_parameter.h"
#include "adjoint.h"
#include "flux.h"
#include "dwdx.h"

void JSparseB_model1_data1_l2v4(SlsMat JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JSparseB->indexvals[0] = 0;
    JSparseB->indexvals[1] = 8;
    JSparseB->indexvals[2] = 1;
    JSparseB->indexvals[3] = 2;
    JSparseB->indexvals[4] = 3;
    JSparseB->indexvals[5] = 8;
    JSparseB->indexvals[6] = 1;
    JSparseB->indexvals[7] = 2;
    JSparseB->indexvals[8] = 3;
    JSparseB->indexvals[9] = 1;
    JSparseB->indexvals[10] = 2;
    JSparseB->indexvals[11] = 3;
    JSparseB->indexvals[12] = 4;
    JSparseB->indexvals[13] = 2;
    JSparseB->indexvals[14] = 4;
    JSparseB->indexvals[15] = 5;
    JSparseB->indexvals[16] = 6;
    JSparseB->indexvals[17] = 4;
    JSparseB->indexvals[18] = 5;
    JSparseB->indexvals[19] = 6;
    JSparseB->indexvals[20] = 7;
    JSparseB->indexvals[21] = 4;
    JSparseB->indexvals[22] = 5;
    JSparseB->indexvals[23] = 6;
    JSparseB->indexvals[24] = 6;
    JSparseB->indexvals[25] = 7;
    JSparseB->indexvals[26] = 0;
    JSparseB->indexvals[27] = 8;
    JSparseB->indexptrs[0] = 0;
    JSparseB->indexptrs[1] = 2;
    JSparseB->indexptrs[2] = 6;
    JSparseB->indexptrs[3] = 9;
    JSparseB->indexptrs[4] = 13;
    JSparseB->indexptrs[5] = 17;
    JSparseB->indexptrs[6] = 21;
    JSparseB->indexptrs[7] = 24;
    JSparseB->indexptrs[8] = 26;
    JSparseB->indexptrs[9] = 28;
    JSparseB->data[0] = 0.0 - 1.0*dwdx0 - 1.0*dwdx1;
    JSparseB->data[1] = -1.0*dwdx13;
    JSparseB->data[2] = 0.0 - 1.0*dwdx2 - 1.0*dwdx3;
    JSparseB->data[3] = 0.0 - 1.0*dwdx4 + 1.0*dwdx5;
    JSparseB->data[4] = -1.0*dwdx6;
    JSparseB->data[5] = 1.0*dwdx14;
    JSparseB->data[6] = 1.0*dwdx2;
    JSparseB->data[7] = 1.0*dwdx4 - 1.0*dwdx5;
    JSparseB->data[8] = 1.0*dwdx6;
    JSparseB->data[9] = -1.0*dwdx2;
    JSparseB->data[10] = -1.0*dwdx4;
    JSparseB->data[11] = -1.0*dwdx6;
    JSparseB->data[12] = 1.0*dwdx8;
    JSparseB->data[13] = 1.0*dwdx5;
    JSparseB->data[14] = -1.0*dwdx7 - 1.0*dwdx8;
    JSparseB->data[15] = -1.0*dwdx9;
    JSparseB->data[16] = 0.0 - 1.0*dwdx10 + 1.0*dwdx11;
    JSparseB->data[17] = -1.0*dwdx7;
    JSparseB->data[18] = -1.0*dwdx9;
    JSparseB->data[19] = -1.0*dwdx10;
    JSparseB->data[20] = 1.0*dwdx12;
    JSparseB->data[21] = 1.0*dwdx7;
    JSparseB->data[22] = 1.0*dwdx9;
    JSparseB->data[23] = 1.0*dwdx10 - 1.0*dwdx11;
    JSparseB->data[24] = 1.0*dwdx11;
    JSparseB->data[25] = -1.0*dwdx12;
    JSparseB->data[26] = 1.0*dwdx0;
    JSparseB->data[27] = 1.0*dwdx13 - 1.0*dwdx14;
}