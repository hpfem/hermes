//////  Physical parameters.  ////////////////////////////////////////////////

// Reference k_effective reactor eigenvalue for the material properties below
// and geometry from the file 'reactor.mesh'. For this example, it was obtained
// simplistically by a reference calculation on a 3x uniformly refined mesh
// with uniform distribution of polynomial degrees (=4), with convergence
// tolerance set to 5e-11.
const double REF_K_EFF = 1.1409144;

// Reflector properties (0), core properties (1).
const double D[2][N_GROUPS] = {
  {0.0164, 0.0085, 0.00832, 0.00821},
  {0.0235, 0.0121, 0.0119, 0.0116}
};
const double Sa[2][N_GROUPS] = {
  {0.00139, 0.000218, 0.00197, 0.0106},
  {0.00977, 0.162, 0.156, 0.535}
};
const double Sr[2][N_GROUPS] = {
  {1.77139, 0.533218, 3.31197, 0.0106},
  {1.23977, 0.529, 2.436, 0.535}
};
const double Sf[2][N_GROUPS] = {
  {0.0, 0.0, 0.0, 0.0}, 
  {0.00395, 0.0262, 0.0718, 0.346}
};
const double nu[2][N_GROUPS] = {
  {0.0, 0.0, 0.0, 0.0}, 
  {2.49, 2.43, 2.42, 2.42}
};
const double chi[2][N_GROUPS] = {
  {0.0, 0.0, 0.0, 0.0}, 
  {0.9675, 0.03250, 0.0, 0.0}
};
const double Ss[2][N_GROUPS][N_GROUPS] = {
  {
    { 0.0,   0.0,  0.0, 0.0},
    {1.77,   0.0,  0.0, 0.0},
    { 0.0, 0.533,  0.0, 0.0},
    { 0.0,   0.0, 3.31, 0.0}
  },
  {
    { 0.0,   0.0,  0.0, 0.0},
    {1.23,   0.0,  0.0, 0.0},
    { 0.0, 0.367,  0.0, 0.0},
    { 0.0,   0.0, 2.28, 0.0}
  }
};