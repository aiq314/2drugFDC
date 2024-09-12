#ifndef DRUG_SIM_HPP
#define DRUG_SIM_HPP

#include "cipa_t.hpp"
#include "param.hpp"
#include "globals.hpp"
#include "../cellmodels/cellmodel.hpp"

// main simulation function for drug simulation 
void do_drug_sim(const double conc_a, const double conc_b,
std::array<double, 14> ic50_a, std::array<double, 14> ic50_b, 
std::array<double, 6> ddis, std::array<double, 18> cvar,
const param_t *p_param, const unsigned short sample_id, 
Cellmodel *p_cell, cvode_t *p_cvode_t, bool is_firsttime);

#endif
