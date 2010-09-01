// This utility generates edge shape functions for quads in the
// space Hcurl. Edge functions are the gradients of scalar Lobatto
// edge functions, bubble functions are eigen functions for the
// curl-curl operator.

#include <stdio.h>
#include <stdlib.h>
#ifdef _MSC_VER
#include <minmax.h>
#endif

void whitney(int i)
{
    printf
    (
      "static double gradeigen_quad_p%d_e1_a_0(double x, double y)\n"
      "{\n"
      "  return Legendre%d(x) * l0(y);\n"
      "}\n\n"
       "static double gradeigen_quad_p%d_e1_a_1(double x, double y)\n"
      "{\n"
      "  return -(Legendre%d(x) * l0(y));\n"
      "}\n\n"
     "static double gradeigen_quad_p%d_e1_b(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"

      "static double gradeigen_quad_p%d_e1_ax_0(double x, double y)\n"
      "{\n"
      "  return Legendre%dx(x) * l0(y);\n"
      "}\n\n"
      "static double gradeigen_quad_p%d_e1_ax_1(double x, double y)\n"
      "{\n"
      "  return -(Legendre%dx(x) * l0(y));\n"
      "}\n\n"
      "static double gradeigen_quad_p%d_e1_ay_0(double x, double y)\n"
      "{\n"
      "  return Legendre%d(x) * dl0(y);\n"
      "}\n\n"
       "static double gradeigen_quad_p%d_e1_ay_1(double x, double y)\n"
      "{\n"
      "  return -(Legendre%d(x) * dl0(y));\n"
      "}\n\n"
     "static double gradeigen_quad_p%d_e1_bx(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"
      "static double gradeigen_quad_p%d_e1_by(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n",
     i,i,i, i,i,   i,i, i,i, i,i,i,i,  i,i
    );


    printf
    (
      "static double gradeigen_quad_p%d_e2_a(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"
      "static double gradeigen_quad_p%d_e2_b_0(double x, double y)\n"
      "{\n"
      "  return l1(x) * Legendre%d(y);\n"
      "}\n\n"
      "static double gradeigen_quad_p%d_e2_b_1(double x, double y)\n"
      "{\n"
      "  return -(l1(x) * Legendre%d(y));\n"
      "}\n\n"

      "static double gradeigen_quad_p%d_e2_ax(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"
      "static double gradeigen_quad_p%d_e2_ay(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"
      "static double gradeigen_quad_p%d_e2_bx_0(double x, double y)\n"
      "{\n"
      "  return dl1(x) * Legendre%d(y);\n"
      "}\n\n"
       "static double gradeigen_quad_p%d_e2_bx_1(double x, double y)\n"
      "{\n"
      "  return -(dl1(x) * Legendre%d(y));\n"
      "}\n\n"
      "static double gradeigen_quad_p%d_e2_by_0(double x, double y)\n"
      "{\n"
      "  return l1(x) * Legendre%dx(y);\n"
      "}\n\n"
      "static double gradeigen_quad_p%d_e2_by_1(double x, double y)\n"
      "{\n"
      "  return -(l1(x) * Legendre%dx(y));\n"
      "}\n\n",
    i,i,i,i,i,    i,i, i,i,i,i,i,i,i,i
    );


    printf
    (
      "static double gradeigen_quad_p%d_e3_a_1(double x, double y)\n"
      "{\n"
      "  return Legendre%d(x) * l1(y);\n"
      "}\n\n"
      "static double gradeigen_quad_p%d_e3_a_0(double x, double y)\n"
      "{\n"
      "  return -(Legendre%d(x) * l1(y));\n"
      "}\n\n"
      "static double gradeigen_quad_p%d_e3_b(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"

      "static double gradeigen_quad_p%d_e3_ax_1(double x, double y)\n"
      "{\n"
      "  return Legendre%dx(x) * l1(y);\n"
      "}\n\n"
      "static double gradeigen_quad_p%d_e3_ax_0(double x, double y)\n"
      "{\n"
      "  return -(Legendre%dx(x) * l1(y));\n"
      "}\n\n"
      "static double gradeigen_quad_p%d_e3_ay_1(double x, double y)\n"
      "{\n"
      "  return Legendre%d(x) * dl1(y);\n"
      "}\n\n"
      "static double gradeigen_quad_p%d_e3_ay_0(double x, double y)\n"
      "{\n"
      "  return -(Legendre%d(x) * dl1(y));\n"
      "}\n\n"
      "static double gradeigen_quad_p%d_e3_bx(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"
      "static double gradeigen_quad_p%d_e3_by(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n",
     i,i,i, i,i,i, i,i,i,  i,i, i,i,   i,i
    );

  printf
  (
      "static double gradeigen_quad_p%d_e4_a(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"
      "static double gradeigen_quad_p%d_e4_b_1(double x, double y)\n"
      "{\n"
      "  return l0(x) * Legendre%d(y);\n"
      "}\n\n"
      "static double gradeigen_quad_p%d_e4_b_0(double x, double y)\n"
      "{\n"
      "  return -(l0(x) * Legendre%d(y));\n"
      "}\n\n"

      "static double gradeigen_quad_p%d_e4_ax(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"
      "static double gradeigen_quad_p%d_e4_ay(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"
      "static double gradeigen_quad_p%d_e4_bx_1(double x, double y)\n"
      "{\n"
      "  return dl0(x) * Legendre%d(y);\n"
      "}\n\n"
      "static double gradeigen_quad_p%d_e4_bx_0(double x, double y)\n"
      "{\n"
      "  return -(dl0(x) * Legendre%d(y));\n"
      "}\n\n"
      "static double gradeigen_quad_p%d_e4_by_1(double x, double y)\n"
      "{\n"
      "  return l0(x) * Legendre%dx(y);\n"
      "}\n\n"
      "static double gradeigen_quad_p%d_e4_by_0(double x, double y)\n"
      "{\n"
      "  return -(l0(x) * Legendre%dx(y));\n"
      "}\n\n",   i,i,i,i,i ,   i,i, i,i,i,i,i,i,i,i
  );
}
void edge_fn(int i, int j)
{
 char c1, c2;
 if ((j<2) && (i%2))
 {
  if (j == 1) { c1 = '-'; c2 = ' ';}
  else {c2 = '-'; c1 = ' ';}
  printf
  (
    "static double gradeigen_quad_l%d_l%d_a_0(double x, double y)\n"
    "{\n"
    " return %cdl%d(x) * l%d(y);\n"
    "}\n\n",
    i, j, c1, i, j
  );
  printf
  (
    "static double gradeigen_quad_l%d_l%d_a_1(double x, double y)\n"
    "{\n"
    " return %cdl%d(x) * l%d(y);\n"
    "}\n\n",
    i, j, c2, i, j
  );

  printf
  (
    "static double gradeigen_quad_l%d_l%d_b(double x, double y)\n"
    "{\n"
    " return l%d(x) * dl%d(y);\n"
    "}\n\n",
    i, j, i, j
  );
  printf
  (
    "static double gradeigen_quad_l%d_l%d_ax_0(double x, double y)\n"
    "{\n"
    " return %cd2l%d(x) * l%d(y);\n"
    "}\n\n",
    i, j, c1, i, j
  );
  printf
  (
    "static double gradeigen_quad_l%d_l%d_ax_1(double x, double y)\n"
    "{\n"
    " return %cd2l%d(x) * l%d(y);\n"
    "}\n\n",
    i, j, c2, i, j
  );

  printf
  (
    "static double gradeigen_quad_l%d_l%d_bx(double x, double y)\n"
    "{\n"
    " return dl%d(x) * dl%d(y);\n"
    "}\n\n",
    i, j, i, j
  );
  printf
  (
    "static double gradeigen_quad_l%d_l%d_ay_0(double x, double y)\n"
    "{\n"
    " return %cdl%d(x) * dl%d(y);\n"
    "}\n\n",
    i, j, c1, i, j
  );
  printf
  (
    "static double gradeigen_quad_l%d_l%d_ay_1(double x, double y)\n"
    "{\n"
    " return %cdl%d(x) * dl%d(y);\n"
    "}\n\n",
    i, j, c2, i, j
  );

  printf
  (
    "static double gradeigen_quad_l%d_l%d_by(double x, double y)\n"
    "{\n"
    " return l%d(x) * d2l%d(y);\n"
    "}\n\n",
    i, j, i, j
  );
 }
 else if ((i<2) && (j%2))
 {
  if (i == 0) { c1 = '-'; c2 = ' ';}
  else {c2 = '-'; c1 = ' ';}
  printf
  (
    "static double gradeigen_quad_l%d_l%d_a(double x, double y)\n"
    "{\n"
    " return dl%d(x) * l%d(y);\n"
    "}\n\n",
    i, j, i, j
  );
  printf
  (
    "static double gradeigen_quad_l%d_l%d_b_0(double x, double y)\n"
    "{\n"
    " return %cl%d(x) * dl%d(y);\n"
    "}\n\n",
    i, j, c1, i, j
  );
  printf
  (
    "static double gradeigen_quad_l%d_l%d_b_1(double x, double y)\n"
    "{\n"
    " return %cl%d(x) * dl%d(y);\n"
    "}\n\n",
    i, j, c2, i, j
  );

  printf
  (
    "static double gradeigen_quad_l%d_l%d_ax(double x, double y)\n"
    "{\n"
    " return d2l%d(x) * l%d(y);\n"
    "}\n\n",
    i, j, i, j
  );
  printf
  (
    "static double gradeigen_quad_l%d_l%d_bx_0(double x, double y)\n"
    "{\n"
    " return %cdl%d(x) * dl%d(y);\n"
    "}\n\n",
    i, j, c1, i, j
  );
  printf
  (
    "static double gradeigen_quad_l%d_l%d_bx_1(double x, double y)\n"
    "{\n"
    " return %cdl%d(x) * dl%d(y);\n"
    "}\n\n",
    i, j, c2, i, j
  );
  printf
  (
    "static double gradeigen_quad_l%d_l%d_ay(double x, double y)\n"
    "{\n"
    " return dl%d(x) * dl%d(y);\n"
    "}\n\n",
    i, j, i, j
  );
  printf
  (
    "static double gradeigen_quad_l%d_l%d_by_0(double x, double y)\n"
    "{\n"
    " return %cl%d(x) * d2l%d(y);\n"
    "}\n\n",
    i, j, c1, i, j
  );
  printf
  (
    "static double gradeigen_quad_l%d_l%d_by_1(double x, double y)\n"
    "{\n"
    " return %cl%d(x) * d2l%d(y);\n"
    "}\n\n",
    i, j, c2, i, j
  );

 }
 else
 {
  printf
  (
    "static double gradeigen_quad_l%d_l%d_a(double x, double y)\n"
    "{\n"
    " return dl%d(x) * l%d(y);\n"
    "}\n\n",
    i, j, i, j
  );
  printf
  (
    "static double gradeigen_quad_l%d_l%d_b(double x, double y)\n"
    "{\n"
    " return l%d(x) * dl%d(y);\n"
    "}\n\n",
    i, j, i, j
  );
  printf
  (
    "static double gradeigen_quad_l%d_l%d_ax(double x, double y)\n"
    "{\n"
    " return d2l%d(x) * l%d(y);\n"
    "}\n\n",
    i, j, i, j
  );
  printf
  (
    "static double gradeigen_quad_l%d_l%d_bx(double x, double y)\n"
    "{\n"
    " return dl%d(x) * dl%d(y);\n"
    "}\n\n",
    i, j, i, j
  );
  printf
  (
    "static double gradeigen_quad_l%d_l%d_ay(double x, double y)\n"
    "{\n"
    " return dl%d(x) * dl%d(y);\n"
    "}\n\n",
    i, j, i, j
  );
  printf
  (
    "static double gradeigen_quad_l%d_l%d_by(double x, double y)\n"
    "{\n"
    " return l%d(x) * d2l%d(y);\n"
    "}\n\n",
    i, j, i, j
  );
 }
}

void bubble_fn(int p1, int p2)
{

  printf("/* gradients of Laplace bubbles */\n\n");
  for (int i = 2; i <= p1; i++)
  {
    for (int j = 2; j <= p2; j++)
    {
     printf(
      "static double gradeigen_quad_p%dp%d_bb_%d_%d_a(double x, double y)\n"
      "{\n"
      "  return deigen_laplace_p%d_%d(x) * eigen_laplace_p%d_%d(y);\n"
      "}\n\n"
      "static double gradeigen_quad_p%dp%d_bb_%d_%d_b(double x, double y)\n"
      "{\n"
      "  return eigen_laplace_p%d_%d(x) * deigen_laplace_p%d_%d(y);\n"
      "}\n\n",
      p1,p2,i,j,  p1,i,p2,j,  p1,p2,i,j,  p1,i,p2,j
     );
     printf(
      "static double gradeigen_quad_p%dp%d_bb_%d_%d_ax(double x, double y)\n"
      "{\n"
      "  return d2eigen_laplace_p%d_%d(x) * eigen_laplace_p%d_%d(y);\n"
      "}\n\n"
      "static double gradeigen_quad_p%dp%d_bb_%d_%d_ay(double x, double y)\n"
      "{\n"
      "  return deigen_laplace_p%d_%d(x) * deigen_laplace_p%d_%d(y);\n"
      "}\n\n"
      "static double gradeigen_quad_p%dp%d_bb_%d_%d_bx(double x, double y)\n"
      "{\n"
      "  return deigen_laplace_p%d_%d(x) * deigen_laplace_p%d_%d(y);\n"
      "}\n\n"
      "static double gradeigen_quad_p%dp%d_bb_%d_%d_by(double x, double y)\n"
      "{\n"
      "  return eigen_laplace_p%d_%d(x) * d2eigen_laplace_p%d_%d(y);\n"
      "}\n\n",
      p1,p2,i,j,  p1,i,p2,j, p1,p2,i,j, p1,i,p2,j, p1,p2,i,j,  p1,i,p2,j, p1,p2,i,j,p1,i,p2,j
     );

    }
  }


  printf("/* Antisymetric gradients */\n\n");
  for (int i = 2; i <= p1; i++)
  {
    for (int j = 2; j <= p2; j++)
    {
     printf(
      "static double gradeigen_quad_p%dp%d_an_%d_%d_a(double x, double y)\n"
      "{\n"
      "  return (- eig%d_%d/(eig%d_%d + eig%d_%d)) * deigen_laplace_p%d_%d(x) * eigen_laplace_p%d_%d(y);\n"
      "}\n\n"
      "static double gradeigen_quad_p%dp%d_an_%d_%d_b(double x, double y)\n"
      "{\n"
      "  return (1 - eig%d_%d/(eig%d_%d + eig%d_%d)) * eigen_laplace_p%d_%d(x) * deigen_laplace_p%d_%d(y);\n"
      "}\n\n",
      p1,p2,i,j,  p1,i, p1,i,p2,j,  p1,i,p2,j,  p1,p2,i,j,   p1,i, p1,i,p2,j,   p1,i,p2,j
     );
     printf(
      "static double gradeigen_quad_p%dp%d_an_%d_%d_ax(double x, double y)\n"
      "{\n"
      "  return (- eig%d_%d/(eig%d_%d + eig%d_%d)) * d2eigen_laplace_p%d_%d(x) * eigen_laplace_p%d_%d(y);\n"
      "}\n\n"
      "static double gradeigen_quad_p%dp%d_an_%d_%d_ay(double x, double y)\n"
      "{\n"
      "  return (- eig%d_%d/(eig%d_%d + eig%d_%d)) * deigen_laplace_p%d_%d(x) * deigen_laplace_p%d_%d(y);\n"
      "}\n\n"
      "static double gradeigen_quad_p%dp%d_an_%d_%d_bx(double x, double y)\n"
      "{\n"
      "  return (1 - eig%d_%d/(eig%d_%d + eig%d_%d)) * deigen_laplace_p%d_%d(x) * deigen_laplace_p%d_%d(y);\n"
       "}\n\n"
      "static double gradeigen_quad_p%dp%d_an_%d_%d_by(double x, double y)\n"
      "{\n"
      "  return (1 - eig%d_%d/(eig%d_%d + eig%d_%d)) * eigen_laplace_p%d_%d(x) * d2eigen_laplace_p%d_%d(y);\n"
      "}\n\n",
      p1,p2,i,j,  p1,i, p1,i,p2,j,  p1,i,p2,j,   p1,p2,i,j,  p1,i, p1,i,p2,j,  p1,i,p2,j, p1,p2,i,j,  p1,i, p1,i,p2,j,  p1,i,p2,j,  p1,p2,i,j, p1,i, p1,i,p2,j, p1,i,p2,j
     );

    }
  }

  printf("/*functions with zero-component*/\n\n");
  for (int i = 2; i <= p1; i++)
  {
    printf(
      "static double gradeigen_quad_p%dp%d_bb_%d_0_a(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"
      "static double gradeigen_quad_p%dp%d_bb_%d_0_b(double x, double y)\n"
      "{\n"
      "  return sqrt(1.0/2.0) * eigen_laplace_p%d_%d(x);\n"
      "}\n\n",
      p1,p2,i,    p1,p2,i,  p1,i
    );
     printf(
      "static double gradeigen_quad_p%dp%d_bb_%d_0_ax(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"
      "static double gradeigen_quad_p%dp%d_bb_%d_0_ay(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"
      "static double gradeigen_quad_p%dp%d_bb_%d_0_bx(double x, double y)\n"
      "{\n"
      "  return sqrt(1.0/2.0) * deigen_laplace_p%d_%d(x);\n"
      "}\n\n"
      "static double gradeigen_quad_p%dp%d_bb_%d_0_by(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n",

     p1,p2,i, p1,p2,i,   p1,p2,i,  p1,i,  p1,p2,i
     );
  }
  for (int j = 2; j <= p2; j++)
    {

	 printf(
      "static double gradeigen_quad_p%dp%d_bb_0_%d_a(double x, double y)\n"
      "{\n"
      "  return sqrt(1.0/2.0) * eigen_laplace_p%d_%d(y);\n"
      "}\n\n"
      "static double gradeigen_quad_p%dp%d_bb_0_%d_b(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"

      "}\n\n",
      p1,p2,j,  p2,j,  p1,p2,j
    );
     printf(
      "static double gradeigen_quad_p%dp%d_bb_0_%d_ax(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"
      "static double gradeigen_quad_p%dp%d_bb_0_%d_ay(double x, double y)\n"
      "{\n"
      "  return sqrt(1.0/2.0) * deigen_laplace_p%d_%d(y);\n"
      "}\n\n"

      "static double gradeigen_quad_p%dp%d_bb_0_%d_bx(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"
      "static double gradeigen_quad_p%dp%d_bb_0_%d_by(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n",
     p1,p2,j, p1,p2,j, p2,j,  p1,p2,j,p1,p2,j
     );
    }
}


int main(int argc, char* argv[])
{
  int i, j, k, l;
  printf("/* Whitney fns - constant tangential component */\n\n");
  whitney(0);

  printf("/* Edge fns - gradients of scalar lobatto edge functions */\n\n");
  for (i = 0; i <= 1; i++)
    for (j = 2; j <= 11; j++)
      edge_fn(i, j);
  for (j = 0; j <= 1; j++)
    for (i = 2; i <= 11; i++)
      edge_fn(i, j);

  printf("/* Bubble fns - eigen fns of the curl-curl operator on the square */");

  for (i = 1; i <= 1; i++)
    {
      for (j = 2; j <= 10; j++)
        {
          printf("/* ORDER p1 = %d, p2 = %d */\n\n",i,j);
		        bubble_fn(i,j);
        }
    }
  for (i = 2; i <= 10; i++)
    {
     for (j = 1; j <= 1; j++)
        {
          printf("/* ORDER p1 = %d, p2 = %d */\n\n",i,j);
		        bubble_fn(i,j);
        }
    }
  for (i = 2; i <= 10; i++)
  {
    for (j = 2; j <= 10; j++)
      {
        printf("/* ORDER p1 = %d, p2 = %d */\n\n",i,j);
		    bubble_fn(i,j);
      }
  }

  printf("static Shapeset::shape_fn_t gradeigen_quad_fn_a[] = \n{\n");

  printf("  gradeigen_quad_p0_e1_a_0, gradeigen_quad_p0_e1_a_1, gradeigen_quad_p0_e2_a, gradeigen_quad_p0_e2_a, gradeigen_quad_p0_e3_a_0, gradeigen_quad_p0_e3_a_1, gradeigen_quad_p0_e4_a, gradeigen_quad_p0_e4_a, \n");

  for (j = 2; j <= 11; j++)
  {
    if ((j%2))
      printf("  gradeigen_quad_l%d_l0_a_0, gradeigen_quad_l%d_l0_a_1, gradeigen_quad_l1_l%d_a, gradeigen_quad_l1_l%d_a,  gradeigen_quad_l%d_l1_a_0, gradeigen_quad_l%d_l1_a_1, gradeigen_quad_l0_l%d_a, gradeigen_quad_l0_l%d_a, ", j,j,j,j,j,j,j,j);
    else
      printf("  gradeigen_quad_l%d_l0_a, gradeigen_quad_l%d_l0_a, gradeigen_quad_l1_l%d_a, gradeigen_quad_l1_l%d_a, gradeigen_quad_l%d_l1_a, gradeigen_quad_l%d_l1_a, gradeigen_quad_l0_l%d_a, gradeigen_quad_l0_l%d_a, ", j,j,j,j,j,j,j,j);
    printf("\n");
  }
  printf("\n");
  k = 0;
  for (int p1 = 2; p1 <= 10; p1++)
  {
    for (int p2 = 2; p2 <= 10; p2++)
    {
      for (i = 2; i <= p1; i++)
      {
        for (j = 2; j <= p2; j++)
        {
          printf("  gradeigen_quad_p%dp%d_bb_%d_%d_a, ",p1, p2, i, j);
          k++;
          if (!(k%5)) printf("\n");
          printf("  gradeigen_quad_p%dp%d_an_%d_%d_a, ",p1, p2, i, j);
          k++;
          if (!(k%5)) printf("\n");
        }
      }
      for (i = 2; i <= p1; i++)
      {
          printf("  gradeigen_quad_p%dp%d_bb_%d_0_a, ",p1, p2, i);
          k++;
          if (!(k%5)) printf("\n");
      }
      for (j = 2; j <= p2; j++)
      {
          printf("  gradeigen_quad_p%dp%d_bb_0_%d_a, ",p1, p2, j);
          k++;
          if (!(k%5)) printf("\n");
      }
      printf("\n\n"); k = 0;
    }
  }
  for (int p1 = 1; p1 <= 1; p1++)
  {
    for (int p2 = 2; p2 <= 10; p2++)
    {
      for (i = 1; i <= 1; i++)
      {
        for (j = 2; j <= p2; j++)
        {
          printf("  gradeigen_quad_p%dp%d_bb_0_%d_a,",p1,p2, j);
          k++;
          if (!(k%5)) printf("\n");
		    }
      }
    }
  }
  for (int p2 = 1; p2 <= 1; p2++)
  {
    for (int p1 = 2; p1 <= 10; p1++)
    {
      for (j = 1; j <= 1; j++)
      {
        for (i = 2; i <= p1; i++)
        {
          printf("  gradeigen_quad_p%dp%d_bb_%d_0_a,",p1,p2, i);
          k++;
          if (!(k%5)) printf("\n");
		    }
      }
    }
  }
  printf("};\n\n");




  printf("static Shapeset::shape_fn_t gradeigen_quad_fn_b[] = \n{\n");

  printf("  gradeigen_quad_p0_e1_b, gradeigen_quad_p0_e1_b, gradeigen_quad_p0_e2_b_0, gradeigen_quad_p0_e2_b_1,  gradeigen_quad_p0_e3_b, gradeigen_quad_p0_e3_b, gradeigen_quad_p0_e4_b_0, gradeigen_quad_p0_e4_b_1, \n");

  for (j = 2; j <= 11; j++)
  {
    if (!(j%2))
      printf("  gradeigen_quad_l%d_l0_b, gradeigen_quad_l%d_l0_b, gradeigen_quad_l1_l%d_b, gradeigen_quad_l1_l%d_b, gradeigen_quad_l%d_l1_b, gradeigen_quad_l%d_l1_b, gradeigen_quad_l0_l%d_b,  gradeigen_quad_l0_l%d_b, ", j,j,j,j,j,j,j,j);
    else
      printf("  gradeigen_quad_l%d_l0_b, gradeigen_quad_l%d_l0_b, gradeigen_quad_l1_l%d_b_0, gradeigen_quad_l1_l%d_b_1, gradeigen_quad_l%d_l1_b, gradeigen_quad_l%d_l1_b, gradeigen_quad_l0_l%d_b_0, gradeigen_quad_l0_l%d_b_1, ", j,j,j,j,j,j,j,j);
    printf("\n");
  }
  printf("\n");
  k = 0;
  for (int p1 = 2; p1 <= 10; p1++)
  {
    for (int p2 = 2; p2 <= 10; p2++)
    {
      for (i = 2; i <= p1; i++)
      {
        for (j = 2; j <= p2; j++)
        {
          printf("  gradeigen_quad_p%dp%d_bb_%d_%d_b, ",p1, p2, i, j);
          k++;
          if (!(k%5)) printf("\n");
          printf("  gradeigen_quad_p%dp%d_an_%d_%d_b, ",p1, p2, i, j);
          k++;
          if (!(k%5)) printf("\n");
        }
      }
      for (i = 2; i <= p1; i++)
      {
          printf("  gradeigen_quad_p%dp%d_bb_%d_0_b, ",p1, p2, i);
          k++;
          if (!(k%5)) printf("\n");
      }
      for (j = 2; j <= p2; j++)
      {
          printf("  gradeigen_quad_p%dp%d_bb_0_%d_b, ",p1, p2, j);
          k++;
          if (!(k%5)) printf("\n");
      }
      printf("\n\n"); k = 0;
    }
  }
  for (int p1 = 1; p1 <= 1; p1++)
  {
    for (int p2 = 2; p2 <= 10; p2++)
    {
      for (i = 1; i <= 1; i++)
      {
        for (j = 2; j <= p2; j++)
        {
          printf("  gradeigen_quad_p%dp%d_bb_0_%d_b,",p1,p2, j);
          k++;
          if (!(k%5)) printf("\n");
		    }
      }
    }
  }
  for (int p2 = 1; p2 <= 1; p2++)
  {
    for (int p1 = 2; p1 <= 10; p1++)
    {
      for (j = 1; j <= 1; j++)
      {
        for (i = 2; i <= p1; i++)
        {
          printf("  gradeigen_quad_p%dp%d_bb_%d_0_b,",p1,p2, i);
          k++;
          if (!(k%5)) printf("\n");
		    }
      }
    }
  }
  printf("};\n\n");


  printf("static Shapeset::shape_fn_t gradeigen_quad_fn_ax[] = \n{\n");

  printf("  gradeigen_quad_p0_e1_ax_0, gradeigen_quad_p0_e1_ax_1, gradeigen_quad_p0_e2_ax, gradeigen_quad_p0_e2_ax, gradeigen_quad_p0_e3_ax_0, gradeigen_quad_p0_e3_ax_1, gradeigen_quad_p0_e4_ax, gradeigen_quad_p0_e4_ax, \n");

  for (j = 2; j <= 11; j++)
  {
    if ((j%2))
      printf("  gradeigen_quad_l%d_l0_ax_0, gradeigen_quad_l%d_l0_ax_1, gradeigen_quad_l1_l%d_ax, gradeigen_quad_l1_l%d_ax,  gradeigen_quad_l%d_l1_ax_0, gradeigen_quad_l%d_l1_ax_1, gradeigen_quad_l0_l%d_ax, gradeigen_quad_l0_l%d_ax, ", j,j,j,j,j,j,j,j);
    else
      printf("  gradeigen_quad_l%d_l0_ax, gradeigen_quad_l%d_l0_ax, gradeigen_quad_l1_l%d_ax, gradeigen_quad_l1_l%d_ax, gradeigen_quad_l%d_l1_ax, gradeigen_quad_l%d_l1_ax, gradeigen_quad_l0_l%d_ax, gradeigen_quad_l0_l%d_ax, ", j,j,j,j,j,j,j,j);
    printf("\n");
  }
  printf("\n");
  k = 0;
  for (int p1 = 2; p1 <= 10; p1++)
  {
    for (int p2 = 2; p2 <= 10; p2++)
    {
      for (i = 2; i <= p1; i++)
      {
        for (j = 2; j <= p2; j++)
        {
          printf("  gradeigen_quad_p%dp%d_bb_%d_%d_ax, ",p1, p2, i, j);
          k++;
          if (!(k%5)) printf("\n");
          printf("  gradeigen_quad_p%dp%d_an_%d_%d_ax, ",p1, p2, i, j);
          k++;
          if (!(k%5)) printf("\n");
        }
      }
      for (i = 2; i <= p1; i++)
      {
          printf("  gradeigen_quad_p%dp%d_bb_%d_0_ax, ",p1, p2, i);
          k++;
          if (!(k%5)) printf("\n");
      }
      for (j = 2; j <= p2; j++)
      {
          printf("  gradeigen_quad_p%dp%d_bb_0_%d_ax, ",p1, p2, j);
          k++;
          if (!(k%5)) printf("\n");
      }
      printf("\n\n"); k = 0;
    }
  }
  for (int p1 = 1; p1 <= 1; p1++)
  {
    for (int p2 = 2; p2 <= 10; p2++)
    {
      for (i = 1; i <= 1; i++)
      {
        for (j = 2; j <= p2; j++)
        {
          printf("  gradeigen_quad_p%dp%d_bb_0_%d_ax,",p1,p2, j);
          k++;
          if (!(k%5)) printf("\n");
		    }
      }
    }
  }
  for (int p2 = 1; p2 <= 1; p2++)
  {
    for (int p1 = 2; p1 <= 10; p1++)
    {
      for (j = 1; j <= 1; j++)
      {
        for (i = 2; i <= p1; i++)
        {
          printf("  gradeigen_quad_p%dp%d_bb_%d_0_ax,",p1,p2, i);
          k++;
          if (!(k%5)) printf("\n");
		    }
      }
    }
  }
  printf("};\n\n");




  printf("static Shapeset::shape_fn_t gradeigen_quad_fn_bx[] = \n{\n");

  printf("  gradeigen_quad_p0_e1_bx, gradeigen_quad_p0_e1_bx, gradeigen_quad_p0_e2_bx_0, gradeigen_quad_p0_e2_bx_1,  gradeigen_quad_p0_e3_bx, gradeigen_quad_p0_e3_bx, gradeigen_quad_p0_e4_bx_0, gradeigen_quad_p0_e4_bx_1, \n");

  for (j = 2; j <= 11; j++)
  {
    if (!(j%2))
      printf("  gradeigen_quad_l%d_l0_bx, gradeigen_quad_l%d_l0_bx, gradeigen_quad_l1_l%d_bx, gradeigen_quad_l1_l%d_bx, gradeigen_quad_l%d_l1_bx, gradeigen_quad_l%d_l1_bx, gradeigen_quad_l0_l%d_bx,  gradeigen_quad_l0_l%d_bx, ", j,j,j,j,j,j,j,j);
    else
      printf("  gradeigen_quad_l%d_l0_bx, gradeigen_quad_l%d_l0_bx, gradeigen_quad_l1_l%d_bx_0, gradeigen_quad_l1_l%d_bx_1, gradeigen_quad_l%d_l1_bx, gradeigen_quad_l%d_l1_bx, gradeigen_quad_l0_l%d_bx_0, gradeigen_quad_l0_l%d_bx_1, ", j,j,j,j,j,j,j,j);
    printf("\n");
  }
  printf("\n");
  k = 0;
  for (int p1 = 2; p1 <= 10; p1++)
  {
    for (int p2 = 2; p2 <= 10; p2++)
    {
      for (i = 2; i <= p1; i++)
      {
        for (j = 2; j <= p2; j++)
        {
          printf("  gradeigen_quad_p%dp%d_bb_%d_%d_bx, ",p1, p2, i, j);
          k++;
          if (!(k%5)) printf("\n");
          printf("  gradeigen_quad_p%dp%d_an_%d_%d_bx, ",p1, p2, i, j);
          k++;
          if (!(k%5)) printf("\n");
        }
      }
      for (i = 2; i <= p1; i++)
      {
          printf("  gradeigen_quad_p%dp%d_bb_%d_0_bx, ",p1, p2, i);
          k++;
          if (!(k%5)) printf("\n");
      }
      for (j = 2; j <= p2; j++)
      {
          printf("  gradeigen_quad_p%dp%d_bb_0_%d_bx, ",p1, p2, j);
          k++;
          if (!(k%5)) printf("\n");
      }
      printf("\n\n"); k = 0;
    }
  }
  for (int p1 = 1; p1 <= 1; p1++)
  {
    for (int p2 = 2; p2 <= 10; p2++)
    {
      for (i = 1; i <= 1; i++)
      {
        for (j = 2; j <= p2; j++)
        {
          printf("  gradeigen_quad_p%dp%d_bb_0_%d_bx,",p1,p2, j);
          k++;
          if (!(k%5)) printf("\n");
		    }
      }
    }
  }
  for (int p2 = 1; p2 <= 1; p2++)
  {
    for (int p1 = 2; p1 <= 10; p1++)
    {
      for (j = 1; j <= 1; j++)
      {
        for (i = 2; i <= p1; i++)
        {
          printf("  gradeigen_quad_p%dp%d_bb_%d_0_bx,",p1,p2, i);
          k++;
          if (!(k%5)) printf("\n");
		    }
      }
    }
  }
  printf("};\n\n");


  printf("static Shapeset::shape_fn_t gradeigen_quad_fn_ay[] = \n{\n");

  printf("  gradeigen_quad_p0_e1_ay_0, gradeigen_quad_p0_e1_ay_1, gradeigen_quad_p0_e2_ay, gradeigen_quad_p0_e2_ay, gradeigen_quad_p0_e3_ay_0, gradeigen_quad_p0_e3_ay_1, gradeigen_quad_p0_e4_ay, gradeigen_quad_p0_e4_ay, \n");

  for (j = 2; j <= 11; j++)
  {
    if ((j%2))
      printf("  gradeigen_quad_l%d_l0_ay_0, gradeigen_quad_l%d_l0_ay_1, gradeigen_quad_l1_l%d_ay, gradeigen_quad_l1_l%d_ay,  gradeigen_quad_l%d_l1_ay_0, gradeigen_quad_l%d_l1_ay_1, gradeigen_quad_l0_l%d_ay, gradeigen_quad_l0_l%d_ay, ", j,j,j,j,j,j,j,j);
    else
      printf("  gradeigen_quad_l%d_l0_ay, gradeigen_quad_l%d_l0_ay, gradeigen_quad_l1_l%d_ay, gradeigen_quad_l1_l%d_ay, gradeigen_quad_l%d_l1_ay, gradeigen_quad_l%d_l1_ay, gradeigen_quad_l0_l%d_ay, gradeigen_quad_l0_l%d_ay, ", j,j,j,j,j,j,j,j);
    printf("\n");
  }
  printf("\n");
  k = 0;
  for (int p1 = 2; p1 <= 10; p1++)
  {
    for (int p2 = 2; p2 <= 10; p2++)
    {
      for (i = 2; i <= p1; i++)
      {
        for (j = 2; j <= p2; j++)
        {
          printf("  gradeigen_quad_p%dp%d_bb_%d_%d_ay, ",p1, p2, i, j);
          k++;
          if (!(k%5)) printf("\n");
          printf("  gradeigen_quad_p%dp%d_an_%d_%d_ay, ",p1, p2, i, j);
          k++;
          if (!(k%5)) printf("\n");
        }
      }
      for (i = 2; i <= p1; i++)
      {
          printf("  gradeigen_quad_p%dp%d_bb_%d_0_ay, ",p1, p2, i);
          k++;
          if (!(k%5)) printf("\n");
      }
      for (j = 2; j <= p2; j++)
      {
          printf("  gradeigen_quad_p%dp%d_bb_0_%d_ay, ",p1, p2, j);
          k++;
          if (!(k%5)) printf("\n");
      }
      printf("\n\n"); k = 0;
    }
  }
  for (int p1 = 1; p1 <= 1; p1++)
  {
    for (int p2 = 2; p2 <= 10; p2++)
    {
      for (i = 1; i <= 1; i++)
      {
        for (j = 2; j <= p2; j++)
        {
          printf("  gradeigen_quad_p%dp%d_bb_0_%d_ay,",p1,p2, j);
          k++;
          if (!(k%5)) printf("\n");
		    }
      }
    }
  }
  for (int p2 = 1; p2 <= 1; p2++)
  {
    for (int p1 = 2; p1 <= 10; p1++)
    {
      for (j = 1; j <= 1; j++)
      {
        for (i = 2; i <= p1; i++)
        {
          printf("  gradeigen_quad_p%dp%d_bb_%d_0_ay,",p1,p2, i);
          k++;
          if (!(k%5)) printf("\n");
		    }
      }
    }
  }
  printf("};\n\n");




  printf("static Shapeset::shape_fn_t gradeigen_quad_fn_by[] = \n{\n");

  printf("  gradeigen_quad_p0_e1_by, gradeigen_quad_p0_e1_by, gradeigen_quad_p0_e2_by_0, gradeigen_quad_p0_e2_by_1,  gradeigen_quad_p0_e3_by, gradeigen_quad_p0_e3_by, gradeigen_quad_p0_e4_by_0, gradeigen_quad_p0_e4_by_1, \n");

  for (j = 2; j <= 11; j++)
  {
    if (!(j%2))
      printf("  gradeigen_quad_l%d_l0_by, gradeigen_quad_l%d_l0_by, gradeigen_quad_l1_l%d_by, gradeigen_quad_l1_l%d_by, gradeigen_quad_l%d_l1_by, gradeigen_quad_l%d_l1_by, gradeigen_quad_l0_l%d_by,  gradeigen_quad_l0_l%d_by, ", j,j,j,j,j,j,j,j);
    else
      printf("  gradeigen_quad_l%d_l0_by, gradeigen_quad_l%d_l0_by, gradeigen_quad_l1_l%d_by_0, gradeigen_quad_l1_l%d_by_1, gradeigen_quad_l%d_l1_by, gradeigen_quad_l%d_l1_by, gradeigen_quad_l0_l%d_by_0, gradeigen_quad_l0_l%d_by_1, ", j,j,j,j,j,j,j,j);
    printf("\n");
  }
  printf("\n");
  k = 0;
  for (int p1 = 2; p1 <= 10; p1++)
  {
    for (int p2 = 2; p2 <= 10; p2++)
    {
      for (i = 2; i <= p1; i++)
      {
        for (j = 2; j <= p2; j++)
        {
          printf("  gradeigen_quad_p%dp%d_bb_%d_%d_by, ",p1, p2, i, j);
          k++;
          if (!(k%5)) printf("\n");
          printf("  gradeigen_quad_p%dp%d_an_%d_%d_by, ",p1, p2, i, j);
          k++;
          if (!(k%5)) printf("\n");
        }
      }
      for (i = 2; i <= p1; i++)
      {
          printf("  gradeigen_quad_p%dp%d_bb_%d_0_by, ",p1, p2, i);
          k++;
          if (!(k%5)) printf("\n");
      }
      for (j = 2; j <= p2; j++)
      {
          printf("  gradeigen_quad_p%dp%d_bb_0_%d_by, ",p1, p2, j);
          k++;
          if (!(k%5)) printf("\n");
      }
      printf("\n\n"); k = 0;
    }
  }
  for (int p1 = 1; p1 <= 1; p1++)
  {
    for (int p2 = 2; p2 <= 10; p2++)
    {
      for (i = 1; i <= 1; i++)
      {
        for (j = 2; j <= p2; j++)
        {
          printf("  gradeigen_quad_p%dp%d_bb_0_%d_by,",p1,p2, j);
          k++;
          if (!(k%5)) printf("\n");
		    }
      }
    }
  }
  for (int p2 = 1; p2 <= 1; p2++)
  {
    for (int p1 = 2; p1 <= 10; p1++)
    {
      for (j = 1; j <= 1; j++)
      {
        for (i = 2; i <= p1; i++)
        {
          printf("  gradeigen_quad_p%dp%d_bb_%d_0_by,",p1,p2, i);
          k++;
          if (!(k%5)) printf("\n");
		    }
      }
    }
  }
  printf("};\n\n");

////////////////////////////////////////////////////////////////////////////

  // m = pocet edge fci v tomto shapesetu
  int m = 88;
  int n = m;
  for (i = 2; i <= 10; i++)
  {
    for (j = 2; j <= 10; j++)
    {
      printf("static int qb_%d_%d[] = { ", i, j);
      for (k = 0; k < 2 * ((i-1) * (j-1)) + (i+j-2); k++)
      {
        printf("%d,", n + k);
      }
      n = n + 2 * ((i-1) * (j-1)) + (i+j-2);
      printf("};\n");
    }
  }
  // n = pocet edge fci + pocet bubble druheho a vyssiho radu
  //printf("%d\n",n);
  for (i = 1; i <= 1; i++)
  {
    for (j = 2; j <= 10; j++)
    {
			printf("static int qb_%d_%d[] = { ", i, j);
      for (l = 2; l <= j; l++)
      {
        printf("%d,", n + l-2);
      }
	    n = n + j-1;
      printf("};\n");
    }
  }
  for (j = 1; j <= 1; j++)
  {
    for (i = 2; i <= 10; i++)
    {
      printf("static int qb_%d_%d[] = { ", i, j);
     for (l = 2; l <= i; l++)
     {
        printf("%d,", n + l-2);
     }
	   n = n + i-1;
    printf("};\n");
    }
  }
  printf("\n\n");

  printf("#define NULL16 NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,\n\n");
  printf("static int* gradeigen_quad_bubble_indices[] =\n{\n");
  //printf("  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL, NULL16\n");
  for (i = 1; i <= 1; i++)
  {
    printf("  NULL, ");
    for (j = 2; j <= 10; j++)
      printf("qb_%d_%d, ", i, j);
    printf(" NULL, NULL, NULL, NULL, NULL, NULL, NULL16\n");
  }
  for (i = 2; i <= 10; i++)
  {
    printf("  ");
    for (j = 1; j <= 10; j++)
      printf("qb_%d_%d, ", i, j);
    printf(" NULL, NULL, NULL, NULL, NULL, NULL, NULL16\n");
  }
  printf("};\n\n");

  printf("#define zero16  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,\n\n");
  printf("static int gradeigen_quad_bubble_count[] =\n{\n");
  //printf("  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, zero16\n");
  for (i = 1; i <= 1; i++)
  {
    printf("  0,  ");
    for (j = 2; j <= 10; j++)
    {
      printf("%d,  ", (j-1));
    }
    printf("0,  0,  0,  0,  0, 0, zero16 \n");
  }
  for (i = 2; i <= 10; i++)
  {
    printf("  ");
    for (j = 1; j <= 10; j++)
    {
      printf("%d, ", 2 * ((i-1) * (j-1)) + (i+j-2));
    }
    printf("0,  0,  0,  0,  0, 0, zero16 \n");
  }
  printf("};\n\n");

  printf("static int gradeigen_quad_vertex_indices[4] ={-1, -1, -1, -1};\n\n");

  for (int edge = 0; edge < 4; edge++)
  {
    printf("static int gradeigen_quad_edge_indices_%d[] = { ", edge);
    k = 0;
    for (i = 0; i <= 10; i++)
    {
      printf("%d, %d, ", 2*edge + k, 2*edge + k +1);
      k += 8;
    }
    printf("};\n\n");
  }
  printf("static int* gradeigen_quad_edge_indices[4] =\n{\n");
  for (int i = 0; i < 4; i++)
    printf("  gradeigen_quad_edge_indices_%d,\n",i);
  printf("};\n\n");

  printf("#define oo make_quad_order\n\n");

  printf("static int gradeigen_quad_index_to_order[] = \n{\n");

  printf("  oo(0,1), oo(0,1), oo(1,0), oo(1,0), oo(0,1), oo(0,1), oo(1,0), oo(1,0),\n");

  for (j = 2; j <= 11; j++)
  {
    printf("  oo(%d,1), oo(%d,1), oo(1,%d), oo(1,%d), oo(%d,1), oo(%d,1), oo(1,%d), oo(1,%d),\n", j,j,j,j,j,j,j,j);
  }
  printf("\n");
  k = 0;
  for (int p1 = 2; p1 <= 10; p1++)
  {
    for (int p2 = 2; p2 <= 10; p2++)
    {
      for (i = 2; i <= p1; i++)
      {
        for (j = 2; j <= p2; j++)
        {
          printf("  oo(%d,%d),",p1, p2);
          k++;
          if (!(k%12)) printf("\n");
          printf("  oo(%d,%d),",p1, p2);
          k++;
          if (!(k%12)) printf("\n");
        }
      }
      for (i = 2; i <= p1; i++)
      {
          printf("  oo(%d,%d),",p1, p2);
          k++;
          if (!(k%12)) printf("\n");
      }
      for (j = 2; j <= p2; j++)
      {
          printf("  oo(%d,%d),",p1, p2);
          k++;
          if (!(k%12)) printf("\n");
      }
      printf("\n\n"); k = 0;
    }
  }
  for (int p1 = 1; p1 <= 1; p1++)
  {
    for (int p2 = 2; p2 <= 10; p2++)
    {
      for (i = 1; i <= 1; i++)
      {
        for (j = 2; j <= p2; j++)
        {
          printf("  oo(%d,%d),",p1,p2);
          k++;
          if (!(k%12)) printf("\n");
		    }
      }
    }
  }
  for (int p2 = 1; p2 <= 1; p2++)
  {
    for (int p1 = 2; p1 <= 10; p1++)
    {
      for (j = 1; j <= 1; j++)
      {
        for (i = 2; i <= p1; i++)
        {
          printf("  oo(%d,%d),",p1,p2);
          k++;
          if (!(k%12)) printf("\n");
		    }
      }
    }
  }
  printf("\n};\n\n");


}

