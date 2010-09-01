// This utility generates edge shape functions for triangles in the
// space Hcurl. Edge functions are the gradients of scalar Lobatto
// edge functions, bubble functions are the Legendre-based functions
// from Higher-Order Finite Element Method (Solin, Segeth, Dolezel).

#include <stdio.h>
#include <stdlib.h>
#ifdef _MSC_VER
#include <minmax.h>
#endif

const double len_phi[3] = {1.0, 1.4142135623731, 1.0};
const double len[3] = {1.0, 1.0, 1.0};


void whitney()
{
  printf("\n///////////////////////////////// ORDER 0 //////////////////////////////////\n\n");
  printf("/* Whitney fns - constant tangential component */\n\n");

  for (int edge = 1; edge < 4; edge++)
  {

    printf("// EDGE %d\n\n", edge);
    printf
    (
      "static double gradleg_tri_p0_e%d_a_0(double x, double y)\n"
      "{\n"
      "  return psi0e%d_1(x,y) / %.13f;\n"
      "}\n\n"
       "static double gradleg_tri_p0_e%d_a_1(double x, double y)\n"
      "{\n"
      "  return -(gradleg_tri_p0_e%d_a_0(x,y));\n"
      "}\n\n"
     "static double gradleg_tri_p0_e%d_b_0(double x, double y)\n"
      "{\n"
      "  return psi0e%d_2(x,y) / %.13f;\n"
      "}\n\n"
     "static double gradleg_tri_p0_e%d_b_1(double x, double y)\n"
      "{\n"
      "  return -(gradleg_tri_p0_e%d_b_0(x,y));\n"
      "}\n\n",
      edge, edge, len_phi[edge-1], edge, edge, edge, edge, len_phi[edge-1], edge, edge
    );
    printf
    (
      "static double gradleg_tri_p0_e%d_ax_0(double x, double y)\n"
      "{\n"
      "  return psi0e%dx_1(x,y) / %.13f;\n"
      "}\n\n"
       "static double gradleg_tri_p0_e%d_ax_1(double x, double y)\n"
      "{\n"
      "  return -(gradleg_tri_p0_e%d_ax_0(x,y));\n"
      "}\n\n"
     "static double gradleg_tri_p0_e%d_bx_0(double x, double y)\n"
      "{\n"
      "  return psi0e%dx_2(x,y) / %.13f;\n"
      "}\n\n"
     "static double gradleg_tri_p0_e%d_bx_1(double x, double y)\n"
      "{\n"
      "  return -(gradleg_tri_p0_e%d_bx_0(x,y));\n"
      "}\n\n",
      edge, edge, len_phi[edge-1], edge, edge, edge, edge, len_phi[edge-1], edge, edge
    );
    printf
    (
      "static double gradleg_tri_p0_e%d_ay_0(double x, double y)\n"
      "{\n"
      "  return psi0e%dy_1(x,y) / %.13f;\n"
      "}\n\n"
       "static double gradleg_tri_p0_e%d_ay_1(double x, double y)\n"
      "{\n"
      "  return -(gradleg_tri_p0_e%d_ay_0(x,y));\n"
      "}\n\n"
     "static double gradleg_tri_p0_e%d_by_0(double x, double y)\n"
      "{\n"
      "  return psi0e%dy_2(x,y) / %.13f;\n"
      "}\n\n"
     "static double gradleg_tri_p0_e%d_by_1(double x, double y)\n"
      "{\n"
      "  return -(gradleg_tri_p0_e%d_by_0(x,y));\n"
      "}\n\n",
      edge, edge, len_phi[edge-1], edge, edge, edge, edge, len_phi[edge-1], edge, edge
    );

  }
}


void edge_fn(int o)
{
 char c1, c2, c3, c4;
 printf("/* EDGE FUNCTIONS - order %d*/\n\n", o-1);

 if (!((o-1) % 2))  // orientace
 {
  for (int ed = 0; ed < 3; ed++)
  {
      printf(" /* EDGE %d */\n\n", ed + 1);
      int e3 = ((ed+2) % 3) +1;
      int e2 = ((ed+1) % 3) +1;
      char c1 = ' ';
      char c2 = '-';
      //if (ed == 2) { c1 = '-'; c2 = ' ';}
      printf
      (
        "static double gradleg_tri_p%d_e%d_a_0(double x, double y)\n"
        "{\n"
        " double l%d, l%dx, l%d, l%dx;\n"
        " l%d = lambda%d(x,y); l%dx = lambda%dx(x,y); l%d = lambda%d(x,y); l%dx = lambda%dx(x,y);\n"
        " return %c(l%dx * l%d * phi%d(l%d - l%d) + l%d * l%dx * phi%d(l%d - l%d) + l%d * l%d * phi%dx(l%d - l%d) * (l%dx - l%dx)) / %.13f;\n"
        "}\n\n",
        o-1,ed+1,
        e3,e3,e2,e2,
        e3,e3,e3,e3,e2,e2,e2,e2,
        c1,
        e3,e2,o-2,e3,e2,
        e3,e2,o-2,e3,e2,
        e3,e2,o-2,e3,e2,e3,e2, len[ed]
      );

      printf
      (
        "static double gradleg_tri_p%d_e%d_a_1(double x, double y)\n"
        "{\n"
        " return -(gradleg_tri_p%d_e%d_a_0(x,y));\n"
        "}\n\n",
        o-1, ed+1,  o-1, ed+1
      );

      printf
      (
        "static double gradleg_tri_p%d_e%d_b_0(double x, double y)\n"
        "{\n"
        " double l%d, l%dy, l%d, l%dy;\n"
        " l%d = lambda%d(x,y); l%dy = lambda%dy(x,y); l%d = lambda%d(x,y); l%dy = lambda%dy(x,y);\n"
        " return %c(l%dy * l%d * phi%d(l%d - l%d) + l%d * l%dy * phi%d(l%d - l%d) + l%d * l%d * phi%dx(l%d - l%d) * (l%dy - l%dy)) / %.13f;\n"
        "}\n\n",
        o-1,ed+1,
        e3,e3,e2,e2,
        e3,e3,e3,e3,e2,e2,e2,e2,
        c1,
        e3,e2,o-2,e3,e2,
        e3,e2,o-2,e3,e2,
        e3,e2,o-2,e3,e2,e3,e2, len[ed]
      );

      printf
      (
        "static double gradleg_tri_p%d_e%d_b_1(double x, double y)\n"
        "{\n"
        " return -(gradleg_tri_p%d_e%d_b_0(x,y));\n"
        "}\n\n",
        o-1, ed+1,  o-1, ed+1
      );

     // derivace
      printf
      (
        "static double gradleg_tri_p%d_e%d_ax_0(double x, double y)\n"
        "{\n"
        " double l%d, l%dx, l%d, l%dx;\n"
        " double ker, kerx, kerxx;\n"
        " l%d = lambda%d(x,y); l%dx = lambda%dx(x,y); l%d = lambda%d(x,y); l%dx = lambda%dx(x,y);\n"
        " ker = phi%d(l%d - l%d); kerx = phi%dx(l%d - l%d) * (l%dx - l%dx); kerxx = phi%dxx(l%d - l%d) * sqr(l%dx - l%dx);\n"
        " return %c(2.0 * l%dx * l%dx * ker + 2.0 * l%dx * l%d * kerx + 2.0 * l%d * l%dx * kerx + l%d * l%d * kerxx) / %.13f;\n"
        "}\n\n",
        o-1,ed+1,
        e3,e3,e2,e2,
        e3,e3,e3,e3,e2,e2,e2,e2,
        o-2, e3,e2,
        o-2,e3,e2,e3,e2,
        o-2,e3,e2,e3,e2,
        c1,
        e3,e2,e3,e2,
        e3,e2,e3,e2,len[ed]
      );
      printf
      (
        "static double gradleg_tri_p%d_e%d_ax_1(double x, double y)\n"
        "{\n"
        " return -(gradleg_tri_p%d_e%d_ax_0(x,y));\n"
        "}\n\n",
        o-1, ed+1,  o-1, ed+1
      );

       printf
      (
        "static double gradleg_tri_p%d_e%d_ay_0(double x, double y)\n"
        "{\n"
        " double l%d, l%dx, l%dy, l%d, l%dx, l%dy;\n"
        " double ker, kerx, kery, kerxy;\n"
        " l%d = lambda%d(x,y); l%dx = lambda%dx(x,y); l%dy = lambda%dy(x,y);\n "
        " l%d = lambda%d(x,y); l%dx = lambda%dx(x,y); l%dy = lambda%dy(x,y); \n"
        " ker = phi%d(l%d - l%d); kerx = phi%dx(l%d - l%d) * (l%dx - l%dx); \n"
        " kery = phi%dx(l%d - l%d) * (l%dy - l%dy); kerxy = phi%dxx(l%d - l%d) * (l%dx - l%dx) * (l%dy - l%dy);\n"
        " return %c(l%dx * l%dy * ker + l%dy * l%dx * ker + l%dx * l%d * kery + l%d * l%dx * kery + l%dy * l%d * kerx + l%d * l%dy * kerx +  l%d * l%d * kerxy) / %.13f;\n"
        "}\n\n",
        o-1,ed+1,
        e3,e3,e3,e2,e2,e2,
        e3,e3,e3,e3,e3,e3,
        e2,e2,e2,e2,e2,e2,
        o-2, e3,e2,
        o-2,e3,e2,e3,e2,
        o-2,e3,e2,e3,e2,
        o-2,e3,e2,e3,e2,e3,e2,
        c1,
        e3,e2,e3,e2,e3,e2,e3,e2,
        e3,e2,e3,e2,e3,e2,len[ed]
      );
      printf
      (
        "static double gradleg_tri_p%d_e%d_ay_1(double x, double y)\n"
        "{\n"
        " return -(gradleg_tri_p%d_e%d_ay_0(x,y));\n"
        "}\n\n",
        o-1, ed+1,  o-1, ed+1
      );
       printf
      (
        "static double gradleg_tri_p%d_e%d_bx_0(double x, double y)\n"
        "{\n"
        " double l%d, l%dx, l%dy, l%d, l%dx, l%dy;\n"
        " double ker, kerx, kery, kerxy;\n"
        " l%d = lambda%d(x,y); l%dx = lambda%dx(x,y); l%dy = lambda%dy(x,y);\n "
        " l%d = lambda%d(x,y); l%dx = lambda%dx(x,y); l%dy = lambda%dy(x,y); \n"
        " ker = phi%d(l%d - l%d); kerx = phi%dx(l%d - l%d) * (l%dx - l%dx); \n"
        " kery = phi%dx(l%d - l%d) * (l%dy - l%dy); kerxy = phi%dxx(l%d - l%d) * (l%dx - l%dx) * (l%dy - l%dy);\n"
        " return %c(l%dx * l%dy * ker + l%dy * l%dx * ker + l%dx * l%d * kery + l%d * l%dx * kery + l%dy * l%d * kerx + l%d * l%dy * kerx +  l%d * l%d * kerxy) / %.13f;\n"
        "}\n\n",
        o-1,ed+1,
        e3,e3,e3,e2,e2,e2,
        e3,e3,e3,e3,e3,e3,
        e2,e2,e2,e2,e2,e2,
        o-2, e3,e2,
        o-2,e3,e2,e3,e2,
        o-2,e3,e2,e3,e2,
        o-2,e3,e2,e3,e2,e3,e2,
        c1,
        e3,e2,e3,e2,e3,e2,e3,e2,
        e3,e2,e3,e2,e3,e2,len[ed]
      );
      printf
      (
        "static double gradleg_tri_p%d_e%d_bx_1(double x, double y)\n"
        "{\n"
        " return -(gradleg_tri_p%d_e%d_bx_0(x,y));\n"
        "}\n\n",
        o-1, ed+1,  o-1, ed+1
      );

      printf
      (
        "static double gradleg_tri_p%d_e%d_by_0(double x, double y)\n"
        "{\n"
        " double l%d, l%dy, l%d, l%dy;\n"
        " double ker, kery, keryy;\n"
        " l%d = lambda%d(x,y); l%dy = lambda%dy(x,y); l%d = lambda%d(x,y); l%dy = lambda%dy(x,y);\n"
        " ker = phi%d(l%d - l%d); kery = phi%dx(l%d - l%d) * (l%dy - l%dy); keryy = phi%dxx(l%d - l%d) * sqr(l%dy - l%dy);\n"
        " return %c(2.0 * l%dy * l%dy * ker + 2.0 * l%dy * l%d * kery + 2.0 * l%d * l%dy * kery + l%d * l%d * keryy) / %.13f;\n"
        "}\n\n",
        o-1,ed+1,
        e3,e3,e2,e2,
        e3,e3,e3,e3,e2,e2,e2,e2,
        o-2, e3,e2,
        o-2,e3,e2,e3,e2,
        o-2,e3,e2,e3,e2,
        c1,
        e3,e2,e3,e2,
        e3,e2,e3,e2,len[ed]
      );
      printf
      (
        "static double gradleg_tri_p%d_e%d_by_1(double x, double y)\n"
        "{\n"
        " return -(gradleg_tri_p%d_e%d_by_0(x,y));\n"
        "}\n\n",
        o-1, ed+1,  o-1, ed+1
      );
    }
  }
  else
  {
  for (int ed = 0; ed < 3; ed++)
  {
      printf(" /* EDGE %d */\n\n", ed + 1);
      int e3 = ((ed+2) % 3) +1;
      int e2 = ((ed+1) % 3) +1;
      char c1 = ' ';
      //if (ed == 2) { c1 = '-';}
      printf
      (
        "static double gradleg_tri_p%d_e%d_a(double x, double y)\n"
        "{\n"
        " double l%d, l%dx, l%d, l%dx;\n"
        " l%d = lambda%d(x,y); l%dx = lambda%dx(x,y); l%d = lambda%d(x,y); l%dx = lambda%dx(x,y);\n"
        " return %c(l%dx * l%d * phi%d(l%d - l%d) + l%d * l%dx * phi%d(l%d - l%d) + l%d * l%d * phi%dx(l%d - l%d) * (l%dx - l%dx)) / %.13f;\n"
        "}\n\n",
        o-1,ed+1,
        e3,e3,e2,e2,
        e3,e3,e3,e3,e2,e2,e2,e2,
        c1,
        e3,e2,o-2,e3,e2,
        e3,e2,o-2,e3,e2,
        e3,e2,o-2,e3,e2,e3,e2, len[ed]
      );


      printf
      (
        "static double gradleg_tri_p%d_e%d_b(double x, double y)\n"
        "{\n"
        " double l%d, l%dy, l%d, l%dy;\n"
        " l%d = lambda%d(x,y); l%dy = lambda%dy(x,y); l%d = lambda%d(x,y); l%dy = lambda%dy(x,y);\n"
        " return %c(l%dy * l%d * phi%d(l%d - l%d) + l%d * l%dy * phi%d(l%d - l%d) + l%d * l%d * phi%dx(l%d - l%d) * (l%dy - l%dy)) / %.13f;\n"
        "}\n\n",
        o-1,ed+1,
        e3,e3,e2,e2,
        e3,e3,e3,e3,e2,e2,e2,e2,
        c1,
        e3,e2,o-2,e3,e2,
        e3,e2,o-2,e3,e2,
        e3,e2,o-2,e3,e2,e3,e2, len[ed]
      );


     // derivace
      printf
      (
        "static double gradleg_tri_p%d_e%d_ax(double x, double y)\n"
        "{\n"
        " double l%d, l%dx, l%d, l%dx;\n"
        " double ker, kerx, kerxx;\n"
        " l%d = lambda%d(x,y); l%dx = lambda%dx(x,y); l%d = lambda%d(x,y); l%dx = lambda%dx(x,y);\n"
        " ker = phi%d(l%d - l%d); kerx = phi%dx(l%d - l%d) * (l%dx - l%dx); kerxx = phi%dxx(l%d - l%d) * sqr(l%dx - l%dx);\n"
        " return %c(2.0 * l%dx * l%dx * ker + 2.0 * l%dx * l%d * kerx + 2.0 * l%d * l%dx * kerx + l%d * l%d * kerxx) / %.13f;\n"
        "}\n\n",
        o-1,ed+1,
        e3,e3,e2,e2,
        e3,e3,e3,e3,e2,e2,e2,e2,
        o-2, e3,e2,
        o-2,e3,e2,e3,e2,
        o-2,e3,e2,e3,e2,
        c1,
        e3,e2,e3,e2,
        e3,e2,e3,e2,len[ed]
      );

       printf
      (
        "static double gradleg_tri_p%d_e%d_ay(double x, double y)\n"
        "{\n"
        " double l%d, l%dx, l%dy, l%d, l%dx, l%dy;\n"
        " double ker, kerx, kery, kerxy;\n"
        " l%d = lambda%d(x,y); l%dx = lambda%dx(x,y); l%dy = lambda%dy(x,y);\n "
        " l%d = lambda%d(x,y); l%dx = lambda%dx(x,y); l%dy = lambda%dy(x,y); \n"
        " ker = phi%d(l%d - l%d); kerx = phi%dx(l%d - l%d) * (l%dx - l%dx); \n"
        " kery = phi%dx(l%d - l%d) * (l%dy - l%dy); kerxy = phi%dxx(l%d - l%d) * (l%dx - l%dx) * (l%dy - l%dy);\n"
        " return %c(l%dx * l%dy * ker + l%dy * l%dx * ker + l%dx * l%d * kery + l%d * l%dx * kery + l%dy * l%d * kerx + l%d * l%dy * kerx +  l%d * l%d * kerxy) / %.13f;\n"
        "}\n\n",
        o-1,ed+1,
        e3,e3,e3,e2,e2,e2,
        e3,e3,e3,e3,e3,e3,
        e2,e2,e2,e2,e2,e2,
        o-2, e3,e2,
        o-2,e3,e2,e3,e2,
        o-2,e3,e2,e3,e2,
        o-2,e3,e2,e3,e2,e3,e2,
        c1,
        e3,e2,e3,e2,e3,e2,e3,e2,
        e3,e2,e3,e2,e3,e2,len[ed]
      );
      printf
      (
        "static double gradleg_tri_p%d_e%d_bx(double x, double y)\n"
        "{\n"
        " double l%d, l%dx, l%dy, l%d, l%dx, l%dy;\n"
        " double ker, kerx, kery, kerxy;\n"
        " l%d = lambda%d(x,y); l%dx = lambda%dx(x,y); l%dy = lambda%dy(x,y);\n "
        " l%d = lambda%d(x,y); l%dx = lambda%dx(x,y); l%dy = lambda%dy(x,y); \n"
        " ker = phi%d(l%d - l%d); kerx = phi%dx(l%d - l%d) * (l%dx - l%dx); \n"
        " kery = phi%dx(l%d - l%d) * (l%dy - l%dy); kerxy = phi%dxx(l%d - l%d) * (l%dx - l%dx) * (l%dy - l%dy);\n"
        " return %c(l%dx * l%dy * ker + l%dy * l%dx * ker + l%dx * l%d * kery + l%d * l%dx * kery + l%dy * l%d * kerx + l%d * l%dy * kerx +  l%d * l%d * kerxy) / %.13f;\n"
        "}\n\n",
        o-1,ed+1,
        e3,e3,e3,e2,e2,e2,
        e3,e3,e3,e3,e3,e3,
        e2,e2,e2,e2,e2,e2,
        o-2, e3,e2,
        o-2,e3,e2,e3,e2,
        o-2,e3,e2,e3,e2,
        o-2,e3,e2,e3,e2,e3,e2,
        c1,
        e3,e2,e3,e2,e3,e2,e3,e2,
        e3,e2,e3,e2,e3,e2,len[ed]
      );

      printf
      (
        "static double gradleg_tri_p%d_e%d_by(double x, double y)\n"
        "{\n"
        " double l%d, l%dy, l%d, l%dy;\n"
        " double ker, kery, keryy;\n"
        " l%d = lambda%d(x,y); l%dy = lambda%dy(x,y); l%d = lambda%d(x,y); l%dy = lambda%dy(x,y);\n"
        " ker = phi%d(l%d - l%d); kery = phi%dx(l%d - l%d) * (l%dy - l%dy); keryy = phi%dxx(l%d - l%d) * sqr(l%dy - l%dy);\n"
        " return %c(2.0 * l%dy * l%dy * ker + 2.0 * l%dy * l%d * kery + 2.0 * l%d * l%dy * kery + l%d * l%d * keryy) / %.13f;\n"
        "}\n\n",
        o-1,ed+1,
        e3,e3,e2,e2,
        e3,e3,e3,e3,e2,e2,e2,e2,
        o-2, e3,e2,
        o-2,e3,e2,e3,e2,
        o-2,e3,e2,e3,e2,
        c1,
        e3,e2,e3,e2,
        e3,e2,e3,e2,len[ed]
      );
    }


  }

}

void edge_bubble(int o)
{
  printf("/* BUBBLE */\n\n");
  printf("/* Edge-based BUBBLE - order %d */\n\n", o);
  for (int ed = 0; ed < 3; ed++)
  {
      printf(" // EDGE %d\n", ed + 1);
      int e3 = ((ed+2) % 3) +1;
      int e2 = ((ed+1) % 3) +1;
      printf
      (
        "static double gradleg_tri_p%d_b%d_a(double x, double y)\n"
        "{\n"
        " double l%d, l%d;\n"
        " l%d = lambda%d(x,y); l%d = lambda%d(x,y); \n"
        " return n%d1 * (l%d * l%d * Legendre%d(l%d - l%d));\n"
        "}\n\n",
        o,ed+1,
        e3,e2,
        e3,e3,e2,e2,
        ed+1, e3,e2,o-2,e3,e2
      );
      printf
      (
        "static double gradleg_tri_p%d_b%d_b(double x, double y)\n"
        "{\n"
        " double l%d, l%d;\n"
        " l%d = lambda%d(x,y); l%d = lambda%d(x,y); \n"
        " return n%d2 * (l%d * l%d * Legendre%d(l%d - l%d));\n"
        "}\n\n",
        o,ed+1,
        e3,e2,
        e3,e3,e2,e2,
        ed+1, e3,e2,o-2,e3,e2
      );
      printf
      (
        "static double gradleg_tri_p%d_b%d_ax(double x, double y)\n"
        "{\n"
        " double l%d, l%dx, l%d, l%dx;\n"
        " l%d = lambda%d(x,y); l%dx = lambda%dx(x,y); l%d = lambda%d(x,y); l%dx = lambda%dx(x,y);\n"
        " return n%d1 * (l%dx * l%d * Legendre%d(l%d - l%d) + l%d * l%dx * Legendre%d(l%d - l%d) + l%d * l%d * Legendre%dx(l%d - l%d) * (l%dx - l%dx));\n"
        "}\n\n",
        o,ed+1,
        e3,e3,e2,e2,
        e3,e3,e3,e3,e2,e2,e2,e2,
        ed+1,
        e3,e2,o-2,e3,e2,
        e3,e2,o-2,e3,e2,
        e3,e2,o-2,e3,e2,e3,e2
      );
      printf
      (
        "static double gradleg_tri_p%d_b%d_bx(double x, double y)\n"
        "{\n"
        " double l%d, l%dx, l%d, l%dx;\n"
        " l%d = lambda%d(x,y); l%dx = lambda%dx(x,y); l%d = lambda%d(x,y); l%dx = lambda%dx(x,y);\n"
        " return n%d2 * (l%dx * l%d * Legendre%d(l%d - l%d) + l%d * l%dx * Legendre%d(l%d - l%d) + l%d * l%d * Legendre%dx(l%d - l%d) * (l%dx - l%dx));\n"
        "}\n\n",
        o,ed+1,
        e3,e3,e2,e2,
        e3,e3,e3,e3,e2,e2,e2,e2,
        ed+1,
        e3,e2,o-2,e3,e2,
        e3,e2,o-2,e3,e2,
        e3,e2,o-2,e3,e2,e3,e2
      );
      printf
      (
        "static double gradleg_tri_p%d_b%d_ay(double x, double y)\n"
        "{\n"
        " double l%d, l%dy, l%d, l%dy;\n"
        " l%d = lambda%d(x,y); l%dy = lambda%dy(x,y); l%d = lambda%d(x,y); l%dy = lambda%dy(x,y);\n"
        " return n%d1 * (l%dy * l%d * Legendre%d(l%d - l%d) + l%d * l%dy * Legendre%d(l%d - l%d) + l%d * l%d * Legendre%dx(l%d - l%d) * (l%dy - l%dy));\n"
        "}\n\n",
        o,ed+1,
        e3,e3,e2,e2,
        e3,e3,e3,e3,e2,e2,e2,e2,
        ed+1,
        e3,e2,o-2,e3,e2,
        e3,e2,o-2,e3,e2,
        e3,e2,o-2,e3,e2,e3,e2
      );
      printf
      (
        "static double gradleg_tri_p%d_b%d_by(double x, double y)\n"
        "{\n"
        " double l%d, l%dy, l%d, l%dy;\n"
        " l%d = lambda%d(x,y); l%dy = lambda%dy(x,y); l%d = lambda%d(x,y); l%dy = lambda%dy(x,y);\n"
        " return n%d2 * (l%dy * l%d * Legendre%d(l%d - l%d) + l%d * l%dy * Legendre%d(l%d - l%d) + l%d * l%d * Legendre%dx(l%d - l%d) * (l%dy - l%dy));\n"
        "}\n\n",
        o,ed+1,
        e3,e3,e2,e2,
        e3,e3,e3,e3,e2,e2,e2,e2,
        ed+1,
        e3,e2,o-2,e3,e2,
        e3,e2,o-2,e3,e2,
        e3,e2,o-2,e3,e2,e3,e2
      );
  }


}

void genuine_bubble(int o)
{
  printf("/* Genuine BUBBLE - order %d */\n\n", o);


  for (int n1 = 1; n1 <= o - 2; n1++) {
    for (int n2 = 1; n2 <= o - 2; n2++) {
      if ((n1 + n2) == o - 1)
      {
        printf
        (
          "static double gradleg_tri_b%d_b%d_1_a(double x, double y)\n"
          "{\n"
          " double l1, l2, l3;\n"
          " l1 = lambda1(x,y); l2 = lambda2(x,y); l3 = lambda3(x,y); \n"
          " return l1 * l2 * l3 * Legendre%d(l3 - l2) * Legendre%d(l2 - l1);\n"
          "}\n\n",
          n1, n2,
          n1 -1, n2-1
        );
        printf
        (
          "static double gradleg_tri_b%d_b%d_1_b(double x, double y)\n"
          "{\n"
          " return 0.0;\n"
          "}\n\n",
          n1, n2
        );

        printf
        (
          "static double gradleg_tri_b%d_b%d_1_ax(double x, double y)\n"
          "{\n"
          " double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;\n"
          " l1 = lambda1(x,y); l2 = lambda2(x,y); l3 = lambda3(x,y); \n"
          " l1x = lambda1x(x,y); l2x = lambda2x(x,y); l3x = lambda3x(x,y); \n"
          " L1 = Legendre%d(l3 - l2); L2 = Legendre%d(l2 - l1); \n"
          " L1x = Legendre%dx(l3 - l2) * (l3x - l2x); L2x = Legendre%dx(l2 - l1) * (l2x - l1x); \n"
          " return l1x * l2 *  l3 *  L1 *  L2 + "
          "         l1 * l2x * l3 *  L1 *  L2 + "
          "         l1 * l2 *  l3x * L1 *  L2 + "
          "         l1 * l2 *  l3 *  L1x * L2 + "
          "         l1 * l2 *  l3 *  L1 *  L2x; \n"
          "}\n\n",
          n1, n2,
          n1 -1, n2-1,
          n1 -1, n2-1
        );
        printf
        (
          "static double gradleg_tri_b%d_b%d_1_bx(double x, double y)\n"
          "{\n"
          " return 0.0;\n"
          "}\n\n",
          n1, n2
        );
        printf
        (
          "static double gradleg_tri_b%d_b%d_1_ay(double x, double y)\n"
          "{\n"
          " double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;\n"
          " l1 = lambda1(x,y); l2 = lambda2(x,y); l3 = lambda3(x,y); \n"
          " l1y = lambda1y(x,y); l2y = lambda2y(x,y); l3y = lambda3y(x,y); \n"
          " L1 = Legendre%d(l3 - l2); L2 = Legendre%d(l2 - l1); \n"
          " L1y = Legendre%dx(l3 - l2) * (l3y - l2y); L2y = Legendre%dx(l2 - l1) * (l2y - l1y); \n"
          " return l1y * l2 *  l3 *  L1 *  L2 + "
          "         l1 * l2y * l3 *  L1 *  L2 + "
          "         l1 * l2 *  l3y * L1 *  L2 + "
          "         l1 * l2 *  l3 *  L1y * L2 + "
          "         l1 * l2 *  l3 *  L1 *  L2y; \n"
          "}\n\n",
          n1, n2,
          n1 -1, n2-1,
          n1 -1, n2-1
        );
        printf
        (
          "static double gradleg_tri_b%d_b%d_1_by(double x, double y)\n"
          "{\n"
          " return 0.0;\n"
          "}\n\n",
          n1, n2
        );

        printf
        (
          "static double gradleg_tri_b%d_b%d_2_b(double x, double y)\n"
          "{\n"
          " double l1, l2, l3;\n"
          " l1 = lambda1(x,y); l2 = lambda2(x,y); l3 = lambda3(x,y); \n"
          " return l1 * l2 * l3 * Legendre%d(l3 - l2) * Legendre%d(l2 - l1);\n"
          "}\n\n",
          n1, n2,
          n1 -1, n2-1
        );
        printf
        (
          "static double gradleg_tri_b%d_b%d_2_a(double x, double y)\n"
          "{\n"
          " return 0.0;\n"
          "}\n\n",
          n1, n2
        );

        printf
        (
          "static double gradleg_tri_b%d_b%d_2_bx(double x, double y)\n"
          "{\n"
          " double l1, l2, l3, l1x, l2x, l3x, L1, L2, L1x, L2x;\n"
          " l1 = lambda1(x,y); l2 = lambda2(x,y); l3 = lambda3(x,y); \n"
          " l1x = lambda1x(x,y); l2x = lambda2x(x,y); l3x = lambda3x(x,y); \n"
          " L1 = Legendre%d(l3 - l2); L2 = Legendre%d(l2 - l1); \n"
          " L1x = Legendre%dx(l3 - l2) * (l3x - l2x); L2x = Legendre%dx(l2 - l1) * (l2x - l1x); \n"
          " return l1x * l2 *  l3 *  L1 *  L2 + "
          "         l1 * l2x * l3 *  L1 *  L2 + "
          "         l1 * l2 *  l3x * L1 *  L2 + "
          "         l1 * l2 *  l3 *  L1x * L2 + "
          "         l1 * l2 *  l3 *  L1 *  L2x; \n"
          "}\n\n",
          n1, n2,
          n1 -1, n2-1,
          n1 -1, n2-1
        );
        printf
        (
          "static double gradleg_tri_b%d_b%d_2_ax(double x, double y)\n"
          "{\n"
          " return 0.0;\n"
          "}\n\n",
          n1, n2
        );
        printf
        (
          "static double gradleg_tri_b%d_b%d_2_by(double x, double y)\n"
          "{\n"
          " double l1, l2, l3, l1y, l2y, l3y, L1, L2, L1y, L2y;\n"
          " l1 = lambda1(x,y); l2 = lambda2(x,y); l3 = lambda3(x,y); \n"
          " l1y = lambda1y(x,y); l2y = lambda2y(x,y); l3y = lambda3y(x,y); \n"
          " L1 = Legendre%d(l3 - l2); L2 = Legendre%d(l2 - l1); \n"
          " L1y = Legendre%dx(l3 - l2) * (l3y - l2y); L2y = Legendre%dx(l2 - l1) * (l2y - l1y); \n"
          " return l1y * l2 *  l3 *  L1 *  L2 + "
          "         l1 * l2y * l3 *  L1 *  L2 + "
          "         l1 * l2 *  l3y * L1 *  L2 + "
          "         l1 * l2 *  l3 *  L1y * L2 + "
          "         l1 * l2 *  l3 *  L1 *  L2y; \n"
          "}\n\n",
          n1, n2,
          n1 -1, n2-1,
          n1 -1, n2-1
        );
        printf
        (
          "static double gradleg_tri_b%d_b%d_2_ay(double x, double y)\n"
          "{\n"
          " return 0.0;\n"
          "}\n\n",
          n1, n2
        );
      }
    }
  }
}


int main(int argc, char* argv[])
{
  int i, j, k, l;
  whitney();

  for (int ord = 1; ord <= 10; ord++)
  {
    printf("///////////////////////////////// ORDER %d //////////////////////////////////\n\n", ord);

    edge_fn(ord+1);
    if (ord > 1) edge_bubble(ord);
    if (ord > 2) genuine_bubble(ord);
  }



  printf("static Shapeset::shape_fn_t gradleg_tri_fn_a[] = \n{\n");

  printf(" gradleg_tri_p0_e1_a_0, gradleg_tri_p0_e1_a_1, gradleg_tri_p0_e2_a_0, "
         " gradleg_tri_p0_e2_a_1, gradleg_tri_p0_e3_a_0, gradleg_tri_p0_e3_a_1, \n");

  for (j = 1; j <= 10; j++)
  {
    if (!(j%2))
      printf(
       " gradleg_tri_p%d_e1_a_0, gradleg_tri_p%d_e1_a_1, gradleg_tri_p%d_e2_a_0, "
       " gradleg_tri_p%d_e2_a_1, gradleg_tri_p%d_e3_a_0, gradleg_tri_p%d_e3_a_1, \n", j,j,j,j,j,j);
    else
      printf(
       " gradleg_tri_p%d_e1_a,   gradleg_tri_p%d_e1_a,   gradleg_tri_p%d_e2_a,   "
       " gradleg_tri_p%d_e2_a,   gradleg_tri_p%d_e3_a,   gradleg_tri_p%d_e3_a,   \n",j,j,j,j,j,j);
  }
  printf("\n");

  k = 66;
  int indices1[15][4];
  int indices2[15][15];
  for (int i = 0; i <= 11; i++)
    for (int j = 0; j <= 11; j++)
    {
      indices1[i][j] = 0;
      indices2[i][j] = 0;
    }

  // edge-based
  for (int i = 2; i <= 10; i++){
    for (int ed = 0; ed < 3; ed++)
    {
      printf("  gradleg_tri_p%d_b%d_a, ", i, ed+1);
      indices1[i][ed] = k;
      k++;
    }
  }
  printf("\n\n");
  for (int i = 3; i <= 10; i++){
   for (int n1 = 1; n1 <= i - 2; n1++) {
    for (int n2 = 1; n2 <= i - 2; n2++) {
      if ((n1 + n2) == i - 1)
      {
        printf("  gradleg_tri_b%d_b%d_1_a, gradleg_tri_b%d_b%d_2_a,", n1, n2, n1, n2);
        indices2[n1][n2] = k;
        k = k+2;
      }
     }
   }
  }
  printf("\n\n");
  printf("\n};\n\n");

  printf("static Shapeset::shape_fn_t gradleg_tri_fn_b[] = \n{\n");

  printf(" gradleg_tri_p0_e1_b_0, gradleg_tri_p0_e1_b_1, gradleg_tri_p0_e2_b_0, "
         " gradleg_tri_p0_e2_b_1, gradleg_tri_p0_e3_b_0, gradleg_tri_p0_e3_b_1, \n");

  for (j = 1; j <= 10; j++)
  {
    if (!(j%2))
      printf(
       " gradleg_tri_p%d_e1_b_0, gradleg_tri_p%d_e1_b_1, gradleg_tri_p%d_e2_b_0, "
       " gradleg_tri_p%d_e2_b_1, gradleg_tri_p%d_e3_b_0, gradleg_tri_p%d_e3_b_1, \n", j,j,j,j,j,j);
    else
      printf(
       " gradleg_tri_p%d_e1_b,   gradleg_tri_p%d_e1_b,   gradleg_tri_p%d_e2_b,   "
       " gradleg_tri_p%d_e2_b,   gradleg_tri_p%d_e3_b,   gradleg_tri_p%d_e3_b,   \n",j,j,j,j,j,j);
  }
  printf("\n");


  // edge-based
  for (int i = 2; i <= 10; i++){
    for (int ed = 0; ed < 3; ed++)
    {
      printf("  gradleg_tri_p%d_b%d_b, ", i, ed+1);
    }
  }
  printf("\n\n");
  for (int i = 3; i <= 10; i++){
   for (int n1 = 1; n1 <= i - 2; n1++) {
    for (int n2 = 1; n2 <= i - 2; n2++) {
      if ((n1 + n2) == i - 1)
      {
        printf("  gradleg_tri_b%d_b%d_1_b, gradleg_tri_b%d_b%d_2_b,", n1, n2, n1, n2);
      }
     }
   }
  }
  printf("\n\n");
  printf("\n};\n\n");


  printf("static Shapeset::shape_fn_t gradleg_tri_fn_ax[] = \n{\n");

  printf(" gradleg_tri_p0_e1_ax_0, gradleg_tri_p0_e1_ax_1, gradleg_tri_p0_e2_ax_0, "
         " gradleg_tri_p0_e2_ax_1, gradleg_tri_p0_e3_ax_0, gradleg_tri_p0_e3_ax_1, \n");

  for (j = 1; j <= 10; j++)
  {
    if (!(j%2))
      printf(
       " gradleg_tri_p%d_e1_ax_0, gradleg_tri_p%d_e1_ax_1, gradleg_tri_p%d_e2_ax_0, "
       " gradleg_tri_p%d_e2_ax_1, gradleg_tri_p%d_e3_ax_0, gradleg_tri_p%d_e3_ax_1, \n", j,j,j,j,j,j);
    else
      printf(
       " gradleg_tri_p%d_e1_ax,   gradleg_tri_p%d_e1_ax,   gradleg_tri_p%d_e2_ax,   "
       " gradleg_tri_p%d_e2_ax,   gradleg_tri_p%d_e3_ax,   gradleg_tri_p%d_e3_ax,   \n",j,j,j,j,j,j);
  }
  printf("\n");


  // edge-based
  for (int i = 2; i <= 10; i++){
    for (int ed = 0; ed < 3; ed++)
    {
      printf("  gradleg_tri_p%d_b%d_ax, ", i, ed+1);
    }
  }
  printf("\n\n");
  for (int i = 3; i <= 10; i++){
   for (int n1 = 1; n1 <= i - 2; n1++) {
    for (int n2 = 1; n2 <= i - 2; n2++) {
      if ((n1 + n2) == i - 1)
      {
        printf("  gradleg_tri_b%d_b%d_1_ax, gradleg_tri_b%d_b%d_2_ax,", n1, n2, n1, n2);
      }
     }
   }
  }
  printf("\n\n");
  printf("\n};\n\n");


  printf("static Shapeset::shape_fn_t gradleg_tri_fn_bx[] = \n{\n");

  printf(" gradleg_tri_p0_e1_bx_0, gradleg_tri_p0_e1_bx_1, gradleg_tri_p0_e2_bx_0, "
         " gradleg_tri_p0_e2_bx_1, gradleg_tri_p0_e3_bx_0, gradleg_tri_p0_e3_bx_1, \n");

  for (j = 1; j <= 10; j++)
  {
    if (!(j%2))
      printf(
       " gradleg_tri_p%d_e1_bx_0, gradleg_tri_p%d_e1_bx_1, gradleg_tri_p%d_e2_bx_0, "
       " gradleg_tri_p%d_e2_bx_1, gradleg_tri_p%d_e3_bx_0, gradleg_tri_p%d_e3_bx_1, \n", j,j,j,j,j,j);
    else
      printf(
       " gradleg_tri_p%d_e1_bx,   gradleg_tri_p%d_e1_bx,   gradleg_tri_p%d_e2_bx,   "
       " gradleg_tri_p%d_e2_bx,   gradleg_tri_p%d_e3_bx,   gradleg_tri_p%d_e3_bx,   \n",j,j,j,j,j,j);
  }
  printf("\n");


  // edge-based
  for (int i = 2; i <= 10; i++){
    for (int ed = 0; ed < 3; ed++)
    {
      printf("  gradleg_tri_p%d_b%d_bx, ", i, ed+1);
    }
  }
  printf("\n\n");
  for (int i = 3; i <= 10; i++){
   for (int n1 = 1; n1 <= i - 2; n1++) {
    for (int n2 = 1; n2 <= i - 2; n2++) {
      if ((n1 + n2) == i - 1)
      {
        printf("  gradleg_tri_b%d_b%d_1_bx, gradleg_tri_b%d_b%d_2_bx,", n1, n2, n1, n2);
      }
     }
   }
  }
  printf("\n\n");
  printf("\n};\n\n");


  printf("static Shapeset::shape_fn_t gradleg_tri_fn_ay[] = \n{\n");

  printf(" gradleg_tri_p0_e1_ay_0, gradleg_tri_p0_e1_ay_1, gradleg_tri_p0_e2_ay_0, "
         " gradleg_tri_p0_e2_ay_1, gradleg_tri_p0_e3_ay_0, gradleg_tri_p0_e3_ay_1, \n");

  for (j = 1; j <= 10; j++)
  {
    if (!(j%2))
      printf(
       " gradleg_tri_p%d_e1_ay_0, gradleg_tri_p%d_e1_ay_1, gradleg_tri_p%d_e2_ay_0, "
       " gradleg_tri_p%d_e2_ay_1, gradleg_tri_p%d_e3_ay_0, gradleg_tri_p%d_e3_ay_1, \n", j,j,j,j,j,j);
    else
      printf(
       " gradleg_tri_p%d_e1_ay,   gradleg_tri_p%d_e1_ay,   gradleg_tri_p%d_e2_ay,   "
       " gradleg_tri_p%d_e2_ay,   gradleg_tri_p%d_e3_ay,   gradleg_tri_p%d_e3_ay,   \n",j,j,j,j,j,j);
  }
  printf("\n");


  // edge-based
  for (int i = 2; i <= 10; i++){
    for (int ed = 0; ed < 3; ed++)
    {
      printf("  gradleg_tri_p%d_b%d_ay, ", i, ed+1);
    }
  }
  printf("\n\n");
  for (int i = 3; i <= 10; i++){
   for (int n1 = 1; n1 <= i - 2; n1++) {
    for (int n2 = 1; n2 <= i - 2; n2++) {
      if ((n1 + n2) == i - 1)
      {
        printf("  gradleg_tri_b%d_b%d_1_ay, gradleg_tri_b%d_b%d_2_ay,", n1, n2, n1, n2);
      }
     }
   }
  }
  printf("\n\n");
  printf("\n};\n\n");


  printf("static Shapeset::shape_fn_t gradleg_tri_fn_by[] = \n{\n");

  printf(" gradleg_tri_p0_e1_by_0, gradleg_tri_p0_e1_by_1, gradleg_tri_p0_e2_by_0, "
         " gradleg_tri_p0_e2_by_1, gradleg_tri_p0_e3_by_0, gradleg_tri_p0_e3_by_1, \n");

  for (j = 1; j <= 10; j++)
  {
    if (!(j%2))
      printf(
       " gradleg_tri_p%d_e1_by_0, gradleg_tri_p%d_e1_by_1, gradleg_tri_p%d_e2_by_0, "
       " gradleg_tri_p%d_e2_by_1, gradleg_tri_p%d_e3_by_0, gradleg_tri_p%d_e3_by_1, \n", j,j,j,j,j,j);
    else
      printf(
       " gradleg_tri_p%d_e1_by,   gradleg_tri_p%d_e1_by,   gradleg_tri_p%d_e2_by,   "
       " gradleg_tri_p%d_e2_by,   gradleg_tri_p%d_e3_by,   gradleg_tri_p%d_e3_by,   \n",j,j,j,j,j,j);
  }
  printf("\n");


  // edge-based
  for (int i = 2; i <= 10; i++){
    for (int ed = 0; ed < 3; ed++)
    {
      printf("  gradleg_tri_p%d_b%d_by, ", i, ed+1);
    }
  }
  printf("\n\n");
  for (int i = 3; i <= 10; i++){
   for (int n1 = 1; n1 <= i - 2; n1++) {
    for (int n2 = 1; n2 <= i - 2; n2++) {
      if ((n1 + n2) == i - 1)
      {
        printf("  gradleg_tri_b%d_b%d_1_by, gradleg_tri_b%d_b%d_2_by,", n1, n2, n1, n2);
      }
     }
   }
  }
  printf("\n\n");
  printf("\n};\n\n");

// ////////////////////////////////////////////////////////////////////////////


  printf("static int gradleg_tri_bubble_indices_all_orders[] = { \n");
  for (int o = 2; o <= 10; o++)
  {
    for (int ed = 0; ed < 3; ed++)
    {
      printf(" %d,", indices1[o][ed] );
    }
    for (int n1 = 1; n1 <= o - 2; n1++)
      for (int n2 = 1; n2 <= o - 2; n2++)
        if ((n1 + n2) == o - 1)
        {
          printf(" %d, %d,", indices2[n1][n2], indices2[n1][n2] + 1 );
        }
    printf("\n");
  }
  printf("};\n  ");

  printf("\n\n");

  printf(
  "static int* gradleg_tri_bubble_indices[11] = \n"
  "{\n"
  "  NULL, NULL,\n"
  "  gradleg_tri_bubble_indices_all_orders,\n"
  "  gradleg_tri_bubble_indices_all_orders,\n"
  "  gradleg_tri_bubble_indices_all_orders,\n"
  "  gradleg_tri_bubble_indices_all_orders,\n"
  "  gradleg_tri_bubble_indices_all_orders,\n"
  "  gradleg_tri_bubble_indices_all_orders,\n"
  "  gradleg_tri_bubble_indices_all_orders,\n"
  "  gradleg_tri_bubble_indices_all_orders,\n"
  "  gradleg_tri_bubble_indices_all_orders \n"
  "};\n");

  printf("\n\n");

  printf("static int gradleg_tri_bubble_count[11] = { 0, 0,");
  for (int o = 2; o <=10; o++)
    printf(" %d," , 3*(o-1) + (o-1) * (o-2));
  printf("};\n\n");

  printf("static int gradleg_tri_edge_indices_0[22] =  { ");
  for (int o = 0; o <=10; o++)
    printf(" %d, %d," , 6*o, 6*o +1 );
  printf("};\n");
  printf("static int gradleg_tri_edge_indices_1[22] =  { ");
  for (int o = 0; o <=10; o++)
    printf(" %d, %d," , 6*o + 2, 6*o +3 );
  printf("};\n");
  printf("static int gradleg_tri_edge_indices_2[22] =  { ");
  for (int o = 0; o <=10; o++)
    printf(" %d, %d," , 6*o + 4, 6*o +5 );
  printf("};\n\n");

  printf(
  "static int* gradleg_tri_edge_indices[3] = \n"
  "{\n"
  "  gradleg_tri_edge_indices_0,\n"
  "  gradleg_tri_edge_indices_1,\n"
  "  gradleg_tri_edge_indices_2,\n"
  "};\n");

  printf("\n");


  printf("static int gradleg_tri_vertex_indices[3] = { -1, -1, -1 };\n\n");


  printf("static int gradleg_tri_index_to_order[] = \n{\n");
  for (j = 0; j <= 10; j++)
  {
    printf(" %d, %d, %d, %d, %d, %d,\n", j,j,j,j,j,j);
  }

  for (int i = 2; i <= 10; i++)
    printf(" %d, %d, %d,\n", i,i,i);

  for (int i = 3; i <= 10; i++){
   for (int n1 = 1; n1 <= i - 2; n1++) {
    for (int n2 = 1; n2 <= i - 2; n2++) {
      if ((n1 + n2) == i - 1)
        printf(" %d, %d,",i,i);
    }
   }
   printf("\n");
  }
  printf("\n};\n\n");


}

