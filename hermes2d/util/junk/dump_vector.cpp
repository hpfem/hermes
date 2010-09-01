void dump_element_sln(FILE* f, Solution* sln, int level)
{
  sln->set_quad_order(1, H2D_FN_VAL);
  scalar* v0 = sln->get_fn_values(0);
  scalar* v1 = sln->get_fn_values(1);
  double* x  = sln->get_refmap()->get_phys_x(1);
  double* y  = sln->get_refmap()->get_phys_y(1);

  for (int i = 6; i < 9; i++)
    fprintf(f, "  %g %g %g %g\n", x[i], y[i], v0[i].real(), v1[i].real());

  if (level < 2)
  {
    for (int i = 0; i < 4; i++)
    {
      sln->push_transform(i);
      dump_element_sln(f, sln, level+1);
      sln->pop_transform();
    }
  }
}


extern Quad2D* get_quad_lin();

void dump_vector_sln(char* filename, Solution* sln)
{
  FILE* f = fopen(filename, "w");
  if (f == NULL) error("asadsd");
  fprintf(f, "A=[\n");

  Mesh* mesh = sln->get_space()->get_mesh();
  sln->set_quad_2d(get_quad_lin());

  Element* e;
  for_all_active_elements(e, mesh)
  {
    sln->set_active_element(e);
    dump_element_sln(f, sln, 1);
  }

  fprintf(f, "];\nquiver(A(:,1), A(:,2), A(:,3), A(:,4), 20);\n");
  fclose(f);
}

