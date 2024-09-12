#ifndef GLOBALS_HPP
#define GLOBALS_HPP

#include <array>
#include <vector>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>


// global variable for MPI.
struct mympi
{
  static char host_name[255];
  static int host_name_len;
  static int rank;
  static int size;
};

// CVode struct
// (made by Ali)
typedef struct {
    void* cvode_mem;
    N_Vector states_vec;
    SUNMatrix matrix;
    SUNLinearSolver solver;
} cvode_t;


// data structure for IC50
typedef std::vector< std::array<double,14> > drug_t;

/*==============*/
/* Added by ALI */
/*==============*/
// data structure for DDIs-LA
typedef std::vector< std::array<double,6> > ddis_t;
// data structure for conductance variations 
typedef std::vector< std::array<double,18> > cvar_t;
#endif
