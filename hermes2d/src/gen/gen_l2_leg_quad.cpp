// Generates L2 shapeset on quads
// based on product of Legendre polynomials.

#include <stdio.h>
#include <stdlib.h>
#ifdef _MSC_VER
#include <minmax.h>
#endif


void output_fn(int i, int j)
{
  printf
    (
      "static double leg_quad_l%d_l%d(double x, double y)\n"
      "{\n"
      "  return Legendre%d(x) * Legendre%d(y);\n"
      "}\n\n"
      "static double leg_quad_l%d_l%dx(double x, double y)\n"
      "{\n"
      "  return Legendre%dx(x) * Legendre%d(y);\n"
      "}\n\n"
      "static double leg_quad_l%d_l%dy(double x, double y)\n"
      "{\n"
      "  return Legendre%d(x) * Legendre%dx(y);\n"
      "}\n\n"
      "static double leg_quad_l%d_l%dxx(double x, double y)\n"
      "{\n"
      "  return Legendre%dxx(x) * Legendre%d(y);\n"
      "}\n\n"
      "static double leg_quad_l%d_l%dxy(double x, double y)\n"
      "{\n"
      "  return Legendre%dx(x) * Legendre%dx(y);\n"
      "}\n\n"
      "static double leg_quad_l%d_l%dyy(double x, double y)\n"
      "{\n"
      "  return  Legendre%d(x) * Legendre%dxx(y);\n"
      "}\n\n",
      i, j, i, j,
      i, j, i, j,
      i, j, i, j,
      i, j, i, j,
      i, j, i, j,
      i, j, i, j
    );
}


int main(int argc, char* argv[])
{
  int i, j, k, l;

  printf("\n");
  printf("#include \"common.h\"\n"
         "#include \"shapeset.h\"\n"
         "#include \"shapeset_common.h\"\n\n");


  printf("//// quad legendre shapeset /////////////////////////////////////////////////////////////////\n\n\n");

  for (i = 0; i <= 10; i++)
    for (j = 0; j <= 10; j++)
      output_fn(i, j);

/////////////////////////////////////////////////////////////////////////

  printf("static Shapeset::shape_fn_t leg_quad_fn[] = \n"
         "{\n  "
  );
  int r = 0;
  for (i = 0; i <= 10; i++)
  {

    for (j = 0; j <= 10; j++)
    {
      {
        printf("leg_quad_l%d_l%d,   ", i, j);
        r++;
        if (r % 5 == 0) printf("\n  ");
      }
    }
  }
  printf("\n};\n");

  printf("static Shapeset::shape_fn_t leg_quad_fn_dx[] =  \n"
         "{\n  "
  );
  r = 0;
  for (i = 0; i <= 10; i++)
  {

    for (j = 0; j <= 10; j++)
    {
      {
        printf("leg_quad_l%d_l%dx,   ", i, j);
        r++;
        if (r % 5 == 0) printf("\n  ");
      }
    }
  }
  printf("\n};\n");

  printf("static Shapeset::shape_fn_t leg_quad_fn_dy[] =  \n"
         "{\n  "
  );
  r = 0;
  for (i = 0; i <= 10; i++)
  {

    for (j = 0; j <= 10; j++)
    {
      {
        printf("leg_quad_l%d_l%dy,   ", i, j);
        r++;
        if (r % 5 == 0) printf("\n  ");
      }
    }
  }
  printf("\n};\n");

  printf("static Shapeset::shape_fn_t leg_quad_fn_dxx[] =  \n"
         "{\n  "
  );
  r = 0;
  for (i = 0; i <= 10; i++)
  {

    for (j = 0; j <= 10; j++)
    {
      {
        printf("leg_quad_l%d_l%dxx,   ", i, j);
        r++;
        if (r % 5 == 0) printf("\n  ");
      }
    }
  }
  printf("\n};\n");

  printf("static Shapeset::shape_fn_t leg_quad_fn_dxy[] =  \n"
         "{\n  "
  );
  r = 0;
  for (i = 0; i <= 10; i++)
  {

    for (j = 0; j <= 10; j++)
    {
      {
        printf("leg_quad_l%d_l%dxy,   ", i, j);
        r++;
        if (r % 5 == 0) printf("\n  ");
      }
    }
  }
  printf("\n};\n");

  printf("static Shapeset::shape_fn_t leg_quad_fn_dyy[] =  \n"
         "{\n  "
  );
  r = 0;
  for (i = 0; i <= 10; i++)
  {

    for (j = 0; j <= 10; j++)
    {
      {
        printf("leg_quad_l%d_l%dyy,   ", i, j);
        r++;
        if (r % 5 == 0) printf("\n  ");
      }
    }
  }
  printf("\n};\n");

//////////////////////////////////////////////////////////////////////////////////////////////

  printf("Shapeset::shape_fn_t* leg_quad_shape_fn_table[1]     = { leg_quad_fn };\n");

  printf("Shapeset::shape_fn_t* leg_quad_shape_fn_table_dx[1]  = { leg_quad_fn_dx };\n");

  printf("Shapeset::shape_fn_t* leg_quad_shape_fn_table_dy[1]  = { leg_quad_fn_dy };\n");

  printf("Shapeset::shape_fn_t* leg_quad_shape_fn_table_dxx[1] = { leg_quad_fn_dxx };\n");

  printf("Shapeset::shape_fn_t* leg_quad_shape_fn_table_dxy[1] = { leg_quad_fn_dxy };\n");

  printf("Shapeset::shape_fn_t* leg_quad_shape_fn_table_dyy[1] = { leg_quad_fn_dyy };\n");

  printf("\n");
//////////////////////////////////////////////////////////////////////////////////////////////


  for (i = 0; i <= 10; i++)
  {
    for (j = 0; j <= 10; j++)
    {
      printf("static int qb_%d_%d[] = { ", i, j);
      for (k = 0; k <= i; k++)
      {
        for (l = 0; l <= j; l++)
        {
          printf("%d,", k*11 + l);
        }
      }
      printf("};\n");
    }
  }
  printf("\n\n");


  printf("#define NULL16 NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL\n\n");

  printf("int* leg_quad_bubble_indices[] =\n{\n");
  for (i = 0; i <= 10; i++)
  {
    printf("  ");
    for (j = 0; j <= 10; j++)
    {
      printf("qb_%d_%d,  ", i, j);
    }
    printf(" NULL, NULL, NULL, NULL, NULL, NULL16,\n");
  }
  printf("};\n\n");


  printf("#define zero16  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0\n\n");

  printf("int leg_quad_bubble_count[] =\n{\n");
  for (j = 0; j <= 10; j++)
  {
    printf("  ");
    for (i = 0; i <= 10; i++)
    {
      printf("  %d,", (i+1) * (j+1));
    }
    printf("  0,  0,  0,  0,  0, zero16,\n");
  }
  printf("};\n\n");

  printf("int leg_quad_vertex_indices[4] = { -1, -1, -1, -1 };\n\n");
  printf("static int leg_quad_edge_indices_0[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };\n"
         "static int leg_quad_edge_indices_1[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };\n"
         "static int leg_quad_edge_indices_2[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };\n"
         "static int leg_quad_edge_indices_3[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };\n"
  );
  printf("\n\n");

  printf("int* leg_quad_edge_indices[4] = \n"
         "{\n"
         "  leg_quad_edge_indices_0,\n"
         "  leg_quad_edge_indices_1,\n"
         "  leg_quad_edge_indices_2,\n"
         "  leg_quad_edge_indices_3\n"
         "};\n"
 );

 printf("\n\n");

 printf("#define oo make_quad_order\n"
        "#define XX(a,b) oo(a,b), oo(a,b)\n\n");

 #define max(a,b) (((a) > (b)) ? (a) : (b))

 printf("int leg_quad_index_to_order[] = {\n");
  for (i = 0; i <= 10; i++)
  {
    for (j = 0; j <= 10; j++)
    {
       printf("  oo(%d,%d), ", max(i,1), max(j,1));
    }
    printf("\n");
  }
  printf("};\n");

  printf("\n\n");


 return 0;

}

