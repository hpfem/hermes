// Generates L2 shapeset on triangles
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
      "static double leg_tri_l%d_l%d(double x, double y)\n"
      "{\n"
      "  return Legendre%d(lambda3(x,y) - lambda2(x,y)) * Legendre%d(lambda2(x,y) - lambda1(x,y));\n"
      "}\n\n"
      "static double leg_tri_l%d_l%dx(double x, double y)\n"
      "{\n"
      "  double l1 = lambda1(x,y), l2 = lambda2(x,y), l3 = lambda3(x,y);\n"
      "  double L1 = Legendre%d(l3 - l2), L1x = Legendre%dx(l3 - l2);\n"
      "  double L2 = Legendre%d(l2 - l1), L2x = Legendre%dx(l2 - l1);\n"
      "  return L1x * (lambda3x(x,y) - lambda2x(x,y)) * L2 + L1 * L2x * (lambda2x(x,y) - lambda1x(x,y));\n"
      "}\n\n"
      "static double leg_tri_l%d_l%dy(double x, double y)\n"
      "{\n"
      "  double l1 = lambda1(x,y), l2 = lambda2(x,y), l3 = lambda3(x,y);\n"
      "  double L1 = Legendre%d(l3 - l2), L1y = Legendre%dx(l3 - l2);\n"
      "  double L2 = Legendre%d(l2 - l1), L2y = Legendre%dx(l2 - l1);\n"
      "  return L1y * (lambda3y(x,y) - lambda2y(x,y)) * L2 + L1 * L2y * (lambda2y(x,y) - lambda1y(x,y));\n"
      "}\n\n",
      i, j, i, j,
      i, j, i, j, i, j,
      i, j, i, j, i, j
    );
}


int main(int argc, char* argv[])
{
  int i, j, k, l;

  printf("\n");
  printf("#include \"h2d_common.h\"\n"
         "#include \"shapeset.h\"\n"
         "#include \"shapeset_l2_all.h\"\n"
         "#include \"shapeset_common.h\"\n\n");


  printf("//// triangle legendre shapeset /////////////////////////////////////////////////////////////////\n\n\n");

  for (i = 0; i <= 10; i++)
  {
    for (j = i; j <= 10 - i; j++)
    {
      output_fn(i, j);
      if (i != j) output_fn(j, i);
    }
  }
/////////////////////////////////////////////////////////////////////////

  printf("static Shapeset::shape_fn_t leg_tri_fn[] = \n"
         "{\n  "  );
  int r = 0;
  for (i = 0; i <= 10; i++)
    for (j = i; j <= 10 - i; j++)
    {
      printf("leg_tri_l%d_l%d,   ", i, j);
      r++;
      if (r % 5 == 0) printf("\n  ");
      if (i != j)
      {
        printf("leg_tri_l%d_l%d,   ", j, i);
        r++;
        if (r % 5 == 0) printf("\n  ");
      }
    }

  printf("\n};\n");

  printf("static Shapeset::shape_fn_t leg_tri_fn_dx[] =  \n"
         "{\n  "  );
  r = 0;
  for (i = 0; i <= 10; i++)
    for (j = i; j <= 10 - i; j++)
    {
      printf("leg_tri_l%d_l%dx,   ", i, j);
      r++;
      if (r % 5 == 0) printf("\n  ");
      if (i != j)
      {
        printf("leg_tri_l%d_l%dx,   ", j, i);
        r++;
        if (r % 5 == 0) printf("\n  ");
      }
    }
  printf("\n};\n");

  printf("static Shapeset::shape_fn_t leg_tri_fn_dy[] =  \n"
         "{\n  "  );
  r = 0;
  for (i = 0; i <= 10; i++)
    for (j = i; j <= 10 - i; j++)
    {
      printf("leg_tri_l%d_l%dy,   ", i, j);
      r++;
      if (r % 5 == 0) printf("\n  ");
      if (i != j)
      {
        printf("leg_tri_l%d_l%dy,   ", j, i);
        r++;
        if (r % 5 == 0) printf("\n  ");
      }
    }
  printf("\n};\n");


//////////////////////////////////////////////////////////////////////////////////////////////

  printf("Shapeset::shape_fn_t* leg_tri_shape_fn_table[1]     = { leg_tri_fn };\n");

  printf("Shapeset::shape_fn_t* leg_tri_shape_fn_table_dx[1]  = { leg_tri_fn_dx };\n");

  printf("Shapeset::shape_fn_t* leg_tri_shape_fn_table_dy[1]  = { leg_tri_fn_dy };\n");

  printf("\n");
//////////////////////////////////////////////////////////////////////////////////////////////

  int tem[11] = {0, 21, 38, 51, 60, 65, 66, 66, 66, 66, 66};

  for (int p = 0; p <= 10; p++)
  {
    printf("static int qb_%d[] = { ", p);
    for (k = 0; k <= p; k++)
    {
      for (l = k; l <= p - k; l++)
      {
        if (l > k) printf("%d, %d, ", tem[k] + 1 + 2*(l - k - 1), tem[k] + 1 + 2*(l - k - 1) + 1);
        else printf("%d, ", tem[k]);
      }
    }
    printf("};\n");
  }
  printf("\n\n");


  printf("int* leg_tri_bubble_indices[11] =\n{\n");
  for (i = 0; i <= 10; i++)
  {
      printf(" qb_%d,  ", i);
  }
  printf("};\n\n");


  printf("int leg_tri_bubble_count[11] =\n{\n");
  for (i = 0; i <= 10; i++)
  {
    printf("  %d,", (i+1) * (i+2) / 2);
  }
  printf("};\n\n");

  printf("int leg_tri_vertex_indices[4] = { -1, -1, -1, -1 };\n\n");
  printf("static int leg_tri_edge_indices_0[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };\n"
         "static int leg_tri_edge_indices_1[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };\n"
         "static int leg_tri_edge_indices_2[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };\n"
         "static int leg_tri_edge_indices_3[22] =  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  };\n"
  );
  printf("\n\n");

  printf("int* leg_tri_edge_indices[4] = \n"
         "{\n"
         "  leg_tri_edge_indices_0,\n"
         "  leg_tri_edge_indices_1,\n"
         "  leg_tri_edge_indices_2,\n"
         "  leg_tri_edge_indices_3\n"
         "};\n"
 );

 printf("\n\n");


 printf("int leg_tri_index_to_order[] = {\n");
 r = 0;
 for (i = 0; i <= 10; i++)
   for (j = i; j <= 10 - i; j++)
 {
   printf("   %d,", i + j);
   r++;
   if (r % 10 == 0) printf("\n");
   if (i != j)
   {
     printf("   %d,", j + i);
     r++;
     if (r % 10 == 0) printf("\n");
   }
 }
 printf("\n};\n");

  printf("\n\n");


 return 0;

}

