#include "commons.hpp"

#include <cstdarg>
#include <cstdio>
#include <cstring>

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>


void mpi_printf(unsigned short node_id, const char *fmt, ...)
{
  if(mympi::rank == node_id){
    va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
  }
}

void mpi_fprintf(unsigned short node_id, FILE *stream, const char *fmt, ...)
{
  if(mympi::rank == node_id){
    va_list args;
    va_start(args, fmt);
    vfprintf(stream, fmt, args);
    va_end(args);
  }
}

int rhs_fn(realtype time,
    N_Vector y,
    N_Vector ydot,
    void* user_data)
{
  Cellmodel* data = (Cellmodel*)user_data;
  data->computeRates(time,
      data->CONSTANTS,
      N_VGetArrayPointer_Serial(ydot),
      N_VGetArrayPointer_Serial(y),
      data->ALGEBRAIC);
  return 0;
}

void init_cvode(cvode_t* p_cvode, Cellmodel* p_cell, double tcurr, bool is_first)
{
  if (is_first) {
    p_cvode->cvode_mem = CVodeCreate(CV_BDF);
    p_cvode->states_vec = N_VMake_Serial(p_cell->states_size, p_cell->STATES);
    p_cvode->matrix = SUNDenseMatrix(p_cell->states_size, p_cell->states_size);
    p_cvode->solver = SUNLinSol_Dense(p_cvode->states_vec, p_cvode->matrix);
    CVodeInit(p_cvode->cvode_mem, rhs_fn, tcurr, p_cvode->states_vec);
    CVodeSetUserData(p_cvode->cvode_mem, p_cell);
    CVodeSStolerances(p_cvode->cvode_mem, 1.0e-7, 1.0e-7);
    CVodeSetMaxStep(p_cvode->cvode_mem, 0.5);
    CVodeSetLinearSolver(p_cvode->cvode_mem, p_cvode->solver, p_cvode->matrix);
    is_first = false;
  }
  else {
    CVodeReInit(p_cvode->cvode_mem, tcurr, p_cvode->states_vec);
  }
}

void clean_cvode(cvode_t* p_cvode)
{
  SUNMatDestroy(p_cvode->matrix);
  SUNLinSolFree(p_cvode->solver);
  N_VDestroy(p_cvode->states_vec);
  CVodeFree(&(p_cvode->cvode_mem));
}

void edison_assign_params(int argc, char *argv[], param_t *p_param)
{
  bool is_default;
  char buffer[100];
  char key[100];
  char value[100];
  char file_name[255];
  FILE *fp_inputdeck;

  // parameters from arguments
  for (int idx = 1; idx < argc; idx += 2) {
    if (!strcasecmp(argv[idx], "-input_deck"))
      strcpy(file_name, argv[idx + 1]);
    else if (!strcasecmp(argv[idx], "-hill_file_a"))
      strcpy(p_param->hill_file_a, argv[idx + 1]);
    else if (!strcasecmp(argv[idx], "-hill_file_b"))
      strcpy(p_param->hill_file_b, argv[idx + 1]);
    /*==============*/
    /* Added by ALI */
    /*==============*/
    else if (!strcasecmp(argv[idx], "-ddis_file"))
      strcpy(p_param->ddis_file, argv[idx + 1]);
    else if (!strcasecmp(argv[idx], "-cvar_file"))
      strcpy(p_param->cvar_file, argv[idx + 1]);
  }  

  is_default = false;
  fp_inputdeck = fopen( file_name, "r");
  if(fp_inputdeck == NULL){
    fprintf(stderr, "Cannot open input deck file %s!!!\nUse default value as the failsafe.\n", file_name);
    is_default = true;
  }

  // read input_deck line by line
  // and store each line to the buffer
  while ( is_default == false && fgets( buffer, 100, fp_inputdeck ) != NULL ) {
    sscanf( buffer, "%s %*s %s", key, value );
    if (strcasecmp(key, "Celltype") == 0) {
      p_param->celltype = strtod( value, NULL );
    }
    else if (strcasecmp(key, "Is_Dutta") == 0) {
      p_param->is_dutta = strtol( value, NULL, 10 );
    }
    else if (strcasecmp(key, "Is_Using_Output") == 0) {
      p_param->is_using_output = strtol( value, NULL, 10 );
    }
    else if (strcasecmp(key, "Basic_Cycle_Length") == 0) {
      p_param->bcl = strtod( value, NULL );
    }
    else if (strcasecmp(key, "Number_of_Pacing") == 0) {
      p_param->pace_max = strtol( value, NULL, 10 );
    }
    else if (strcasecmp(key, "Time_Step") == 0) {
      p_param->dt = strtod( value, NULL );
    }
    else if (strcasecmp(key, "Write_Step") == 0) {
      p_param->dt_write = strtod( value, NULL );
    }
    else if (strcasecmp(key, "Drug_Name_A") == 0) {
      strcpy( p_param->drug_name_a, value );
    }
    else if (strcasecmp(key, "Drug_Name_B") == 0) {
      strcpy( p_param->drug_name_b, value );
    }  
    else if (strcasecmp(key, "Inet_Vm_Threshold") == 0) {
      p_param->inet_vm_threshold = strtod( value, NULL );
    }
    else if (strcasecmp(key, "Cmax_A") == 0) {
      p_param->cmax_a = strtod( value, NULL );
    }
    else if (strcasecmp(key, "Cmax_B") == 0) {
      p_param->cmax_b = strtod( value, NULL );
    }
    else if (strcasecmp(key, "Radius") == 0) {
      strcat( p_param->radius, value );
    }
    else if (strcasecmp(key, "Angles") == 0) {
      strcat( p_param->angles, value );
    }
    else if (strcasecmp(key, "Is_DDIs_LA") == 0) {
      p_param->is_ddis_la = strtol( value, NULL, 10 );
    }
    else if (strcasecmp(key, "Is_Syntopic") == 0) {
      p_param->is_syntopic = strtol( value, NULL, 10 );
    }
    else if (strcasecmp(key, "Is_Cvar") == 0) {
      p_param->is_cvar = strtol( value, NULL, 10 );
    }

  }

  if( is_default == false ) fclose( fp_inputdeck );
}
