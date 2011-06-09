class CustomWeakFormPoisson : public WeakForm<double>
{
public:
  CustomWeakFormPoisson(std::string mat_al);
};

CustomWeakFormPoisson::CustomWeakFormPoisson(std::string mat_al) : WeakForm<double>(1)
{
  // Jacobian forms.
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(0, 0, mat_al));
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(0, 0));

  // Residual forms.
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(0));
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(0));
  add_vector_form(new WeakFormsH1::DefaultVectorFormVol<double>(0, HERMES_ANY));
};