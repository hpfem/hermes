#ifndef __FCT_H
#define __FCT_H

#include "hermes2d.h"
#include "lumped_projection.h"
#include "reg_estimator.h"
using namespace Hermes;
using namespace Hermes::Hermes2D;

class Flux_Correction
{

public:

	Flux_Correction(double theta);
	~Flux_Correction();

	void free();
	
	void init(SpaceSharedPtr<double> new_space);

	UMFPackMatrix<double>* artificialDiffusion( UMFPackMatrix<double>* conv_matrix);
	UMFPackMatrix<double>* massLumping( UMFPackMatrix<double>* mass_matrix);

	void antidiffusiveFlux(UMFPackMatrix<double>* mass_matrix,UMFPackMatrix<double>* lumped_matrix,UMFPackMatrix<double>* conv_matrix,UMFPackMatrix<double>* diffusion,double* u_high, double* u_L, double* u_old,double* flux_scalar,double time_step, Regularity_Estimator* regEst=NULL);

	
	void project_FCT(MeshFunctionSharedPtr<double> sln, double* coeff_vec, double* coeff_vec_2,UMFPackMatrix<double>* mass_matrix,UMFPackMatrix<double>* lumped_matrix, double time_step, OGProjection<double>* ogProjection,	Lumped_Projection* lumpedProjection, Regularity_Estimator* regEst=NULL);
	
	SpaceSharedPtr<double> get_space(){return space;}

protected:
		void lumped_flux_limiter(UMFPackMatrix<double>* mass_matrix,UMFPackMatrix<double>* lumped_matrix, double* u_L, double* u_H,double time_step, int* smooth_dof=NULL);

double theta;
	bool* fct;
	SpaceSharedPtr<double> space;
	AsmList<double>*  al;
	 double* P_plus; 
	 double* P_minus; 
	 double* Q_plus; 
	 double* Q_minus;  
	 double* R_plus; 
	 double* R_minus;



};





#endif
