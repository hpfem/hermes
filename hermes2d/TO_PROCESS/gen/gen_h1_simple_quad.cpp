// This utility generates simple shape functions for quads.

#include <stdio.h>
#include <stdlib.h>
#ifdef _MSC_VER
#include <minmax.h>
#endif


void output_fn(int i, int j)
{
  char c = (i == 0 && j > 1 && (j & 1) || j == 1 && i > 1 && (i & 1)) ? '-' : ' ';

  if (((i == 0 || i == 1) && (j & 1) && (j != 1)) || ((j == 0 || j == 1) && (i & 1) && (i != 1))) {
    printf
    (
      "static double simple_quad_l%d_l%d_0(double x, double y)\n"
      "{\n"
      "  return %c l%d(x) * l%d(y);\n"
      "}\n\n"
      "static double simple_quad_l%d_l%dx_0(double x, double y)\n"
      "{\n"
      "  return %c dl%d(x) * l%d(y);\n"
      "}\n\n"
      "static double simple_quad_l%d_l%dy_0(double x, double y)\n"
      "{\n"
      "  return %c l%d(x) * dl%d(y);\n"
      "}\n\n"
      "static double simple_quad_l%d_l%dxx_0(double x, double y)\n"
      "{\n"
      "  return %c d2l%d(x) * l%d(y);\n"
      "}\n\n"
      "static double simple_quad_l%d_l%dxy_0(double x, double y)\n"
      "{\n"
      "  return %c dl%d(x) * dl%d(y);\n"
      "}\n\n"
      "static double simple_quad_l%d_l%dyy_0(double x, double y)\n"
      "{\n"
      "  return %c l%d(x) * d2l%d(y);\n"
      "}\n\n"
      "static double simple_quad_l%d_l%d_1(double x, double y)\n"
      "{\n"
      "  return -(%c l%d(x) * l%d(y));\n"
      "}\n\n"
      "static double simple_quad_l%d_l%dx_1(double x, double y)\n"
      "{\n"
      "  return -(%c dl%d(x) * l%d(y));\n"
      "}\n\n"
      "static double simple_quad_l%d_l%dy_1(double x, double y)\n"
      "{\n"
      "  return -(%c l%d(x) * dl%d(y));\n"
      "}\n\n"
      "static double simple_quad_l%d_l%dxx_1(double x, double y)\n"
      "{\n"
      "  return -(%c d2l%d(x) * l%d(y));\n"
      "}\n\n"
      "static double simple_quad_l%d_l%dxy_1(double x, double y)\n"
      "{\n"
      "  return -(%c dl%d(x) * dl%d(y));\n"
      "}\n\n"
      "static double simple_quad_l%d_l%dyy_1(double x, double y)\n"
      "{\n"
      "  return -(%c l%d(x) * d2l%d(y));\n"
      "}\n\n",
      i, j, c, i, j,
      i, j, c, i, j,
      i, j, c, i, j,
      i, j, c, i, j,
      i, j, c, i, j,
      i, j, c, i, j,
      i, j, c, i, j,
      i, j, c, i, j,
      i, j, c, i, j,
      i, j, c, i, j,
      i, j, c, i, j,
      i, j, c, i, j
    );
  }
  else {
    printf
    (
      "static double simple_quad_l%d_l%d(double x, double y)\n"
      "{\n"
      "  return %c l%d(x) * l%d(y);\n"
      "}\n\n"
      "static double simple_quad_l%d_l%dx(double x, double y)\n"
      "{\n"
      "  return %c dl%d(x) * l%d(y);\n"
      "}\n\n"
      "static double simple_quad_l%d_l%dy(double x, double y)\n"
      "{\n"
      "  return %c l%d(x) * dl%d(y);\n"
      "}\n\n"
      "static double simple_quad_l%d_l%dxx(double x, double y)\n"
      "{\n"
      "  return %c d2l%d(x) * l%d(y);\n"
      "}\n\n"
      "static double simple_quad_l%d_l%dxy(double x, double y)\n"
      "{\n"
      "  return %c dl%d(x) * dl%d(y);\n"
      "}\n\n"
      "static double simple_quad_l%d_l%dyy(double x, double y)\n"
      "{\n"
      "  return %c l%d(x) * d2l%d(y);\n"
      "}\n\n",
      i, j, c, i, j,
      i, j, c, i, j,
      i, j, c, i, j,
      i, j, c, i, j,
      i, j, c, i, j,
      i, j, c, i, j
    );
  }
}


int main(int argc, char* argv[])
{
  int i, j, k, l;

  printf("\n");
  printf("#include \"h2d_common.h\"\n"
         "#include \"shapeset.h\"\n"
         "#include \"shapeset_common.h\"\n\n");


  printf("//// quad simple lobatto shapeset /////////////////////////////////////////////////////////////////\n\n\n");

  for (i = 0; i <= 10; i++)
    for (j = 0; j <= 10; j++)
      output_fn(i, j);

/////////////////////////////////////////////////////////////////////////

  printf("static Shapeset::shape_fn_t simple_quad_fn[] = \n"
         "{\n  "
  );
  int r = 0;
  for (i = 0; i <= 10; i++)
  {

    for (j = 0; j <= 10; j++)
    {
      if (((i == 0 || i == 1) && (j & 1) && (j != 1)) || ((j == 0 || j == 1) && (i & 1) && (i != 1)))
      {
        printf("simple_quad_l%d_l%d_0, ", i, j);
        r++;
        if (r % 5 == 0) printf("\n  ");
        printf("simple_quad_l%d_l%d_1, ", i, j);
        r++;
        if (r % 5 == 0) printf("\n  ");
      }
      else
      {
        printf("simple_quad_l%d_l%d,   ", i, j);
        r++;
        if (r % 5 == 0) printf("\n  ");
      }
    }
  }
  printf("\n};\n");

  printf("static Shapeset::shape_fn_t simple_quad_fn_dx[] =  \n"
         "{\n  "
  );
  r = 0;
  for (i = 0; i <= 10; i++)
  {

    for (j = 0; j <= 10; j++)
    {
      if (((i == 0 || i == 1) && (j & 1) && (j != 1)) || ((j == 0 || j == 1) && (i & 1) && (i != 1)))
      {
        printf("simple_quad_l%d_l%dx_0, ", i, j);
        r++;
        if (r % 5 == 0) printf("\n  ");
        printf("simple_quad_l%d_l%dx_1, ", i, j);
        r++;
        if (r % 5 == 0) printf("\n  ");
      }
      else
      {
        printf("simple_quad_l%d_l%dx,   ", i, j);
        r++;
        if (r % 5 == 0) printf("\n  ");
      }
    }
  }
  printf("\n};\n");

  printf("static Shapeset::shape_fn_t simple_quad_fn_dy[] =  \n"
         "{\n  "
  );
  r = 0;
  for (i = 0; i <= 10; i++)
  {

    for (j = 0; j <= 10; j++)
    {
      if (((i == 0 || i == 1) && (j & 1) && (j != 1)) || ((j == 0 || j == 1) && (i & 1) && (i != 1)))
      {
        printf("simple_quad_l%d_l%dy_0, ", i, j);
        r++;
        if (r % 5 == 0) printf("\n  ");
        printf("simple_quad_l%d_l%dy_1, ", i, j);
        r++;
        if (r % 5 == 0) printf("\n  ");
      }
      else
      {
        printf("simple_quad_l%d_l%dy,   ", i, j);
        r++;
        if (r % 5 == 0) printf("\n  ");
      }
    }
  }
  printf("\n};\n");

  printf("static Shapeset::shape_fn_t simple_quad_fn_dxx[] =  \n"
         "{\n  "
  );
  r = 0;
  for (i = 0; i <= 10; i++)
  {

    for (j = 0; j <= 10; j++)
    {
      if (((i == 0 || i == 1) && (j & 1) && (j != 1)) || ((j == 0 || j == 1) && (i & 1) && (i != 1)))
      {
        printf("simple_quad_l%d_l%dxx_0, ", i, j);
        r++;
        if (r % 5 == 0) printf("\n  ");
        printf("simple_quad_l%d_l%dxx_1, ", i, j);
        r++;
        if (r % 5 == 0) printf("\n  ");
      }
      else
      {
        printf("simple_quad_l%d_l%dxx,   ", i, j);
        r++;
        if (r % 5 == 0) printf("\n  ");
      }
    }
  }
  printf("\n};\n");

  printf("static Shapeset::shape_fn_t simple_quad_fn_dxy[] =  \n"
         "{\n  "
  );
  r = 0;
  for (i = 0; i <= 10; i++)
  {

    for (j = 0; j <= 10; j++)
    {
      if (((i == 0 || i == 1) && (j & 1) && (j != 1)) || ((j == 0 || j == 1) && (i & 1) && (i != 1)))
      {
        printf("simple_quad_l%d_l%dxy_0, ", i, j);
        r++;
        if (r % 5 == 0) printf("\n  ");
        printf("simple_quad_l%d_l%dxy_1, ", i, j);
        r++;
        if (r % 5 == 0) printf("\n  ");
      }
      else
      {
        printf("simple_quad_l%d_l%dxy,   ", i, j);
        r++;
        if (r % 5 == 0) printf("\n  ");
      }
    }
  }
  printf("\n};\n");

  printf("static Shapeset::shape_fn_t simple_quad_fn_dyy[] =  \n"
         "{\n  "
  );
  r = 0;
  for (i = 0; i <= 10; i++)
  {

    for (j = 0; j <= 10; j++)
    {
      if (((i == 0 || i == 1) && (j & 1) && (j != 1)) || ((j == 0 || j == 1) && (i & 1) && (i != 1)))
      {
        printf("simple_quad_l%d_l%dyy_0, ", i, j);
        r++;
        if (r % 5 == 0) printf("\n  ");
        printf("simple_quad_l%d_l%dyy_1, ", i, j);
        r++;
        if (r % 5 == 0) printf("\n  ");
      }
      else
      {
        printf("simple_quad_l%d_l%dyy,   ", i, j);
        r++;
        if (r % 5 == 0) printf("\n  ");
      }
    }
  }
  printf("\n};\n");

//////////////////////////////////////////////////////////////////////////////////////////////

  printf("Shapeset::shape_fn_t* simple_quad_shape_fn_table[1]     = { simple_quad_fn };\n");

  printf("Shapeset::shape_fn_t* simple_quad_shape_fn_table_dx[1]  = { simple_quad_fn_dx };\n");

  printf("Shapeset::shape_fn_t* simple_quad_shape_fn_table_dy[1]  = { simple_quad_fn_dy };\n");

  printf("Shapeset::shape_fn_t* simple_quad_shape_fn_table_dxx[1] = { simple_quad_fn_dxx };\n");

  printf("Shapeset::shape_fn_t* simple_quad_shape_fn_table_dxy[1] = { simple_quad_fn_dxy };\n");

  printf("Shapeset::shape_fn_t* simple_quad_shape_fn_table_dyy[1] = { simple_quad_fn_dyy };\n");

  printf("\n");
//////////////////////////////////////////////////////////////////////////////////////////////

  int vol[11] = {0, 0, 0, 13, 24, 37, 48, 61, 72, 85, 96};

  for (i = 2; i <= 10; i++)
  {
    for (j = 2; j <= 10; j++)
    {
      printf("static int qb_%d_%d[] = { ", i, j);
      for (k = 2; k <= i; k++)
      {
        for (l = 2; l <= j; l++)
        {
          printf("%d,", 30 + vol[k] + l);
        }
      }
      printf("};\n");
    }
  }
  printf("\n\n");


  printf("#define NULL16 NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL\n\n");

  printf("int* simple_quad_bubble_indices[] =\n{\n");
  printf("  NULL, NULL, NULL,    NULL,    NULL,    NULL,    NULL,    NULL,    NULL,    NULL,    NULL,     NULL, NULL, NULL, NULL, NULL, NULL16,\n");
  printf("  NULL, NULL, NULL,    NULL,    NULL,    NULL,    NULL,    NULL,    NULL,    NULL,    NULL,     NULL, NULL, NULL, NULL, NULL, NULL16,\n");
  for (i = 2; i <= 10; i++)
  {
    printf("  NULL, NULL, ");
    for (j = 2; j <= 10; j++)
    {
      if (j == 10) printf("qb_%d_%d, ", i, j);
      else if (i == 10) printf("qb_%d_%d, ", i, j);
      else printf("qb_%d_%d,  ", i, j);
    }
    printf(" NULL, NULL, NULL, NULL, NULL, NULL16,\n");
  }
  printf("};\n\n");


  printf("#define zero16  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0\n\n");

  printf("int simple_quad_bubble_count[] =\n{\n");
  printf("  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, zero16,\n");
  printf("  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, zero16,\n");
  for (j = 2; j <= 10; j++)
  {
    printf("  0,  0,");
    for (i = 2; i <= 10; i++)
    {
      if ((i-1) * (j-1) > 9) printf(" %d,", (i-1) * (j-1));
      else printf("  %d,", (i-1) * (j-1));
    }
    printf("  0,  0,  0,  0,  0, zero16,\n");
  }
  printf("};\n\n");

  printf("int simple_quad_vertex_indices[4] = { 0, 15, 16, 1 };\n\n");
  printf("static int simple_quad_edge_indices_0[22] =  {  0, 15, 15, 0,  30, 30, 41, 42, 54, 54, 65, 66, 78, 78, 89, 90, 102, 102, 113, 114, 126, 126 };\n"
         "static int simple_quad_edge_indices_1[22] =  { 15, 16, 16, 15, 17, 17, 18, 19, 20, 20, 21, 22, 23, 23, 24, 25, 26,   26,  27,  28,  29,  29 };\n"
         "static int simple_quad_edge_indices_2[22] =  { 16,  1,  1, 16, 31, 31, 43, 44, 55, 55, 67, 68, 79, 79, 91, 92, 103, 103, 115, 116, 127, 127 };\n"
         "static int simple_quad_edge_indices_3[22] =  {  1,  0,  0,  1,  2,  2,  3,  4,  5,  5,  6,  7,  8,  8,  9, 10,  11,  11,  12,  13,  14,  14 };\n"
  );
  printf("\n\n");

  printf("int* simple_quad_edge_indices[4] = \n"
         "{\n"
         "  simple_quad_edge_indices_0,\n"
         "  simple_quad_edge_indices_1,\n"
         "  simple_quad_edge_indices_2,\n"
         "  simple_quad_edge_indices_3\n"
         "};\n"
 );

 printf("\n\n");

 printf("#define oo make_quad_order\n"
        "#define XX(a,b) oo(a,b), oo(a,b)\n\n");

 #define max(a,b) (((a) > (b)) ? (a) : (b))

 printf("int simple_quad_index_to_order[] = {\n");
  for (i = 0; i <= 10; i++)
  {
    for (j = 0; j <= 10; j++)
    {
      if (((i == 0 || i == 1) && (j & 1) && (j != 1)) || ((j == 0 || j == 1) && (i & 1) && (i != 1)))
        printf("  XX(%d,%d), ", max(i,1), max(j,1));
      else {
        if (i == 10) printf("  oo(%d,%d),", max(i,1), max(j,1));
        else printf("  oo(%d,%d), ", max(i,1), max(j,1));
      }
    }
    printf("\n");
  }
  printf("};\n");

  printf("\n\n");


 return 0;

}

