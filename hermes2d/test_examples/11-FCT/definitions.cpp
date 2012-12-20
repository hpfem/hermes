#include "definitions.h"
//---------------mass-matrix/tau-----------
 CustomWeakFormMassmatrix::CustomWeakFormMassmatrix(double time_step) : WeakForm<double>(1)
  {
		add_matrix_form(new CustomMatrixFormVolMassmatrix(0, 0, time_step));
  }


    template<typename Real, typename Scalar>
    Scalar CustomMatrixFormVolMassmatrix::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, Func<Scalar>  **ext) const 
   {
				Scalar result = Scalar(0); 
				for (int i = 0; i < n; i++)
					result += wt[i] * (u->val[i] * v->val[i])/time_step;
				return result;
    };

   double CustomMatrixFormVolMassmatrix::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, Func<double>  **ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    };

    Ord CustomMatrixFormVolMassmatrix::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, Func<Ord>  **ext) const 
    {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    };

MatrixFormVol<double>* CustomMatrixFormVolMassmatrix::clone() const
{
  return new CustomMatrixFormVolMassmatrix(*this);
}


//---------------Convection-----------

  CustomWeakFormConvection::CustomWeakFormConvection() : WeakForm<double>(1) 
  {
    add_matrix_form(new CustomMatrixFormVolConvection(0, 0));
  };


  template<typename Real, typename Scalar>
  Scalar CustomMatrixFormVolConvection::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                     Func<Real> *v, Geom<Real> *e, Func<Scalar>  **ext) const
  {
		Scalar result = Scalar(0); 
		for (int i = 0; i < n; i++)
			result += -wt[i] * (v->val[i] *(u->dx[i] * (0.5- e->y[i]) + u->dy[i] * (e->x[i]-0.5) ));
		return result;

  };

	double CustomMatrixFormVolConvection::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
               Func<double> *v, Geom<double> *e, Func<double>  **ext) const                  
  {
    return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
  };

 Ord CustomMatrixFormVolConvection::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
          Geom<Ord> *e, Func<Ord>  **ext) const 
 {
    return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
  };

 MatrixFormVol<double>* CustomMatrixFormVolConvection::clone() const
{
	return new CustomMatrixFormVolConvection(*this);
	}



//------------------- Initial condition ----------------

void CustomInitialCondition::derivatives(double x, double y, double& dx, double& dy) const 
{
      
  double radius = 0.;
        //hump
	double x_0 =0.25;
	double y_0= 0.5;	
	radius = (1.0/0.15) * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
	if( radius<= 1.0) 
	{		
		dx = -std::sin(radius*PI)/4.0*(PI/(0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2*x;
		dy = -std::sin(radius*PI)/4.0*(PI/(0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2*y;	
	}
	else
	{			
		//cone
		x_0 = 0.5;
		y_0 = 0.25;
		radius = 1.0/0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
		if((radius< 1.0)&&(x!=x_0)) 
		{ 	
			dx = -(1.0/(0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2*x;
			dy = -(1.0/(0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2*y;	
		}
		else
		{
			dx=0.; dy=0.;
		}	
  
	}
		

};

 double CustomInitialCondition::value(double x, double y) const 
{
       
  double result = 0.0;
	double radius;
        //hump
	double x_0 =0.25;
	double y_0= 0.5;	
	radius = (1.0/0.15) * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
	if( radius<= 1.0) 
	{ 
		 result = (1.0+ std::cos(PI*radius))/4.0;
		return result;	
	}
	//slotted cylinder
	x_0 = 0.5;
	y_0 = 0.75;
	radius = 1.0/0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
	if(radius <= 1) 
	{ 	
		if(fabs((x-x_0))>= 0.025) return 1.0;
		if(y>=0.85) return 1.0;
	}	
	//cone
	x_0 = 0.5;
	y_0 = 0.25;
	radius = 1.0/0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
	if(radius<= 1.0) 
	{ 	
		result = 1.0-radius;
	}	
       return result;
};

 Ord CustomInitialCondition::ord(Ord x, Ord y) const 
 {
      return Ord(10);
	};
 MeshFunction<double>* CustomInitialCondition::clone() const
	{
		return new CustomInitialCondition(this->mesh);

	}



