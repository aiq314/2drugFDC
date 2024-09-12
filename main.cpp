#ifdef TOMEK_2019
#include "cellmodels/Tomek_model.hpp"
#else 
#include "cellmodels/Ohara_Rudy_2011.hpp"
#endif
#include "modules/cipa_t.hpp"
#include "modules/commons.hpp"
#include "modules/globals.hpp"

/*==============*/
/* Added by ALI */
/*==============*/
#include "modules/drug_sim.hpp"
#include "modules/dcomb_t.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mpi.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

// POSIX global variable
struct stat st = {0};

// get the IC50 data from file
drug_t get_IC50_data_from_file(const char* file_name);

/*==============*/
/* Added by ALI */
/*==============*/
// get the DDIs-LA data from file
ddis_t get_ddis_data_from_file(const char* file_name);
// get the cvar data from file
cvar_t get_cvar_data_from_file(const char* file_name);

// the main simulation code
std::vector<double> get_data(char *param_data);
int main(int argc, char **argv)
{
  // buffer for writing in snprintf() function
  char buffer[255];
  setvbuf( stdout, NULL, _IONBF, 0 );

  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &mympi::size );
  MPI_Comm_rank( MPI_COMM_WORLD, &mympi::rank );
  MPI_Get_processor_name(mympi::host_name, &mympi::host_name_len);

  // input parameter object
  param_t *p_param;
  p_param = (param_t*)malloc( sizeof( param_t ) );
  p_param->init();
  edison_assign_params(argc,argv,p_param);
  p_param->show_val();

  // cell model and drug induced-related vars
  drug_t ic50_a, ic50_b;
  Cellmodel* p_cell;
  /*==============*/
  /* Added by ALI */
  /*==============*/
  ddis_t ddis_LA;
  std::array<double,6> temp_ddis;
  cvar_t cvar;
  std::array<double,18> temp_cvar;

  // Cvode vars
  cvode_t *ode_solver;
  bool cvode_firsttime;

  /*==============*/
  /* Added by ALI */
  /*==============*/
  // get radius, and angles from input
  std::vector<double> radius;
  std::vector<double> angles;
  radius = get_data(p_param->radius);
  angles = get_data(p_param->angles);
  // extract IC50 data from file and store it into vector
  //snprintf(buffer, sizeof(buffer), 
  //  "./drugs/bepridil/IC50_samples10.csv");
  ic50_a = get_IC50_data_from_file(p_param->hill_file_a);
  ic50_b = get_IC50_data_from_file(p_param->hill_file_b);
  if(ic50_a.size() == 0){
    mpi_printf(0, "Something problem with the IC50_A file!\n");
    MPI_Finalize();
  }
  else if(ic50_b.size() == 0){
    mpi_printf(0, "Something problem with the IC50_B file!\n");
    MPI_Finalize();
  }
  else if(ic50_a.size() != ic50_b.size()){
    mpi_printf(0, "%s\n%s\n", 
                "The size of IC50_A file is not the same as IC50_B file!",
                "Make sure the size of both IC50_A and IC50_B files are the same");
    MPI_Finalize();
  }
  else if(ic50_a.size() > 100000){
    mpi_printf(0, "Too much input! Maximum sample data is 100000!\n");
    MPI_Finalize();
  }
  else if(p_param->pace_max < 750 && p_param->pace_max > 1000){
    mpi_printf(0, "Make sure the maximum pace is around 750 to 1000!\n");
    MPI_Finalize();
  }
  else if(mympi::size > ic50_a.size()){
    mpi_printf(0, "%s\n%s\n", 
                "Overflow of MPI Process!",
                "Make sure MPI Size is less than or equal the number of sample");
    MPI_Finalize();
  }
  else if(p_param->is_ddis_la == true){
    // extract DDIs data from file and store it into vector
    ddis_LA = get_ddis_data_from_file(p_param->ddis_file);
    if(ic50_a.size() != ddis_LA.size()){
      mpi_printf(0, "%s\n%s\n", 
                  "The size of IC50 file is not the same as DDIs file!",
                  "Make sure the size of both IC50 and DDIs files are the same");
      MPI_Finalize();
    }
  }
  if(p_param->is_cvar == true){
    // extract cvar data from file and store it into vector
    cvar = get_cvar_data_from_file(p_param->cvar_file);
    if(ic50_a.size() != cvar.size()){
      mpi_printf(0, "%s\n%s\n", 
                  "The size of IC50 file is not the same as cvar file!",
                  "Make sure the size of both IC50 and cvar files are the same");
      MPI_Finalize();
    }
  }

  // make a directory for each concentration
  // and wait till all directories have been created
  double conc_a, conc_b, pi = 2.0*acos(0.0), temp_angle;
  if(mympi::rank == 0){
    mkdir("result/", 0775);
    for( const auto &rad: radius )
    { // begin concentration loop conc_a
      for( const auto &angle: angles )
      {
        temp_angle = angle*pi/180.0; // angle in radian
        conc_a = p_param->cmax_a*rad*cos(temp_angle);
        conc_b = p_param->cmax_b*rad*sin(temp_angle);
        snprintf( buffer, sizeof(buffer), "result/%.2lf_%.2lf", conc_a, conc_b );
        if(stat("result/", &st) == 0) mkdir(buffer, 0775);
      } // end concentration loop conc_b
    } // end concentration loop conc_a
  }
  MPI_Barrier(MPI_COMM_WORLD);
#ifdef TOMEK_2019
  mpi_printf(0,"Using Tomek cell model\n");
  p_cell = new Tomek_model();
#else
  mpi_printf(0,"Using O'Hara Rudy cell model\n");
  p_cell = new Ohara_Rudy_2011();
#endif
  ode_solver = (cvode_t*) malloc(sizeof(cvode_t));
  cvode_firsttime = true;
 
  MPI_Barrier(MPI_COMM_WORLD);
  double t_begin = MPI_Wtime();
  unsigned short sample_id;
  for( sample_id = mympi::rank;
       sample_id < ic50_a.size();
       sample_id += mympi::size )
  { // begin sample loop drug a and b
    printf("Sample_ID:%d  Count:%d Rank:%d\n", 
          sample_id, ic50_a.size(), mympi::rank );
    for( const auto &rad: radius )
    { // begin radius loop
      for( const auto &angle: angles )
      { // begin angles loop
	if(rad == 0.0 && angle > 0.0){
	  continue;
        }
	else {	
	  // fill the temp_ddis with data if is_ddis_la == true
          if(p_param->is_ddis_la == true){
	    temp_ddis = ddis_LA[sample_id];
  	    mpi_printf(0, "DDIs_LA is used\n");
          }
          if(p_param->is_cvar == true){
	    temp_cvar = cvar[sample_id];
          }
          temp_angle = angle*pi/180.0; // angle in radian
          conc_a = p_param->cmax_a*rad*cos(temp_angle);
          conc_b = p_param->cmax_b*rad*sin(temp_angle);
          mpi_printf(0,"rad = %.2lf, angle = %.2lf \n",rad,angle);
  	  do_drug_sim(conc_a, conc_b, 
		    ic50_a[sample_id], ic50_b[sample_id],
		    temp_ddis, temp_cvar,
	  	    p_param, sample_id,
		    p_cell, ode_solver, cvode_firsttime);
          if(cvode_firsttime == true) cvode_firsttime = false;
	}
      } //end angles loop
    } // end radius loop
  } // end sample loop drug a
  double t_end = MPI_Wtime();
  MPI_Barrier(MPI_COMM_WORLD);

  // memory cleaning and finalize the program
  clean_cvode(ode_solver);
  free(ode_solver);
  free(p_param);
  delete p_cell;

  mpi_printf(0, "Computation Time: %lf minutes\n", (t_end-t_begin)/60.);

  MPI_Finalize();
}

/*=================*/
/* Modified by ALI */
/*=================*/
std::vector<double> get_data(char *param_data)
{
  std::vector<double> data;
  char *token = strtok( param_data, "," );
  while( token != NULL )
  { // begin data tokenizing
      data.push_back( strtod(token, NULL) );
      token = strtok(NULL, ",");
  } // end data tokenizing
  return data;
}

drug_t get_IC50_data_from_file(const char* file_name)
{
  // buffer for writing in snprintf() function
  char buffer[255];
  FILE *fp_drugs;
  drug_t ic50;
  char *token;
  std::array<double,14> temp_array;
  unsigned short idx;

  if( (fp_drugs = fopen(file_name, "r")) == NULL){
    printf("Cannot open file %s in %s at rank %d\n", 
      file_name, mympi::host_name, mympi::rank);
    return ic50;
  }

  fgets(buffer, sizeof(buffer), fp_drugs); // skip header
  while( fgets(buffer, sizeof(buffer), fp_drugs) != NULL )
  { // begin line reading
    token = strtok( buffer, "," );
    idx = 0;
    while( token != NULL )
    { // begin data tokenizing
      temp_array[idx++] = strtod(token, NULL);
      token = strtok(NULL, ",");
    } // end data tokenizing
    ic50.push_back(temp_array);
  } // end line reading

  fclose(fp_drugs);
  return ic50;
}

/*==============*/
/* Added by ALI */
/*==============*/
ddis_t get_ddis_data_from_file(const char* file_name)
{
  // buffer for writing in snprintf() function
  char buffer[255];
  //FILE *fp_drugs;
  //drug_t ic50;
  FILE *fp_ddis;
  ddis_t ddis_LA;
  char *token;
  std::array<double,6> temp_array;
  unsigned short idx;

  if( (fp_ddis = fopen(file_name, "r")) == NULL){
    printf("Cannot open file %s in %s at rank %d\n", 
      file_name, mympi::host_name, mympi::rank);
    return ddis_LA;
  }

  fgets(buffer, sizeof(buffer), fp_ddis); // skip header
  while( fgets(buffer, sizeof(buffer), fp_ddis) != NULL )
  { // begin line reading
    token = strtok( buffer, "," );
    idx = 0;
    while( token != NULL )
    { // begin data tokenizing
      temp_array[idx++] = strtod(token, NULL);
      token = strtok(NULL, ",");
    } // end data tokenizing
    ddis_LA.push_back(temp_array);
  } // end line reading

  fclose(fp_ddis);
  return ddis_LA;
}
cvar_t get_cvar_data_from_file(const char* file_name)
{
  // buffer for writing in snprintf() function
  char buffer[255];
  FILE *fp_cvar;
  cvar_t cvar;
  char *token;
  std::array<double,18> temp_array;
  unsigned short idx;

  if( (fp_cvar = fopen(file_name, "r")) == NULL){
    printf("Cannot open file %s in %s at rank %d\n", 
      file_name, mympi::host_name, mympi::rank);
    return cvar;
  }

  fgets(buffer, sizeof(buffer), fp_cvar); // skip header
  while( fgets(buffer, sizeof(buffer), fp_cvar) != NULL )
  { // begin line reading
    token = strtok( buffer, "," );
    idx = 0;
    while( token != NULL )
    { // begin data tokenizing
      temp_array[idx++] = strtod(token, NULL);
      token = strtok(NULL, ",");
    } // end data tokenizing
    cvar.push_back(temp_array);
  } // end line reading

  fclose(fp_cvar);
  return cvar;
}
