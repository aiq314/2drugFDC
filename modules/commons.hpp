#ifndef COMMONS_HPP
#define COMMONS_HPP

#include <cstdio>

#include "globals.hpp"
#include "param.hpp"
#include "../cellmodels/cellmodel.hpp"

// custom printf for MPI
// to avoid duplicate printing
void mpi_printf(unsigned short node_id, const char *fmt, ...);
void mpi_fprintf(unsigned short node_id, FILE *stream, const char *fmt, ...);

// CVode setup functions
void init_cvode(cvode_t* p_cvode, Cellmodel* p_cell, double tcurr, bool is_first);
void clean_cvode(cvode_t* p_cvode);

// parameter setup function
void edison_assign_params(int argc, char *argv[], param_t *p_param);

#endif
