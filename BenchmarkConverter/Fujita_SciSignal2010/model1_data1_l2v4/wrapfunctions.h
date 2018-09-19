#ifndef _amici_wrapfunctions_h
#define _amici_wrapfunctions_h
#include "model1_data1_l2v4.h"

/**
 * @brief Wrapper function to instantiate the linked Amici model without knowing the name at compile time.
 * @return
 */
std::unique_ptr<amici::Model> getModel();

#endif /* _amici_wrapfunctions_h */