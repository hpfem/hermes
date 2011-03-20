#ifndef __SOURCE_FILTER_H
#define __SOURCE_FILTER_H

class SourceFilter : public Filter
{
public:
  typedef int (*pos2mat_fn) (double,double);
  
  SourceFilter(MeshFunction* sln1, MeshFunction* sln2, const double fission_xsec[][2], pos2mat_fn get_mat) 
  	: Filter(sln1, sln2), Sf(fission_xsec), get_material(get_mat) {};

  virtual scalar get_pt_value(double x, double y, int item = H2D_FN_VAL_0)
  { error("Not implemented yet"); return 0; }

protected:

  const double(*Sf)[2];
  pos2mat_fn get_material;
  
  virtual void precalculate(int order, int mask);
};

#endif
