
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <GL/glut.h>

#include "hermes2d.h"

#define PI 3.14159265

HcurlShapesetLegendre shapeset;

int shape_n = 0;

static float ramp(float x)
  {
    if (x < -0.25) return 0;
    else if (x < 0) return (x + 0.25) * 4;
    else if (x < 0.25) return 1;
    else if (x < 0.5) return (0.5 - x) * 4;
    else return 0;
  }

float* get_palette_entry(float x)
  {
    static float color[3];
    color[0] = ramp(x - 0.625);
    color[1] = ramp(x - 0.370);
    color[2] = ramp(x - 0.125);
    return color;
  }


void DrawMesh_quad(int steps)
{
  shapeset.set_mode(MODE_QUAD);
  double max = 0.0;
  for (int i = 0; i <= steps; i++)
  {
    for (int j = 0; j <= steps; j++)
    {
      double x = -1.0 + i * 2.0/steps;
      double y = -1.0 + j * 2.0/steps;

      double a = shapeset.get_fn_value(shape_n, x, y, 0);
      double b = shapeset.get_fn_value(shape_n, x, y, 1);
      double length = sqrt(a*a + b*b);
      if (length > max) max = length;
    }
  }

  for (int i = 0; i <= steps; i++)
  {
    for (int j = 0; j <= steps; j++)
    {
      double x = -1.0 + i * 2.0/steps;
      double y = -1.0 + j * 2.0/steps;

      double a = shapeset.get_fn_value(shape_n, x, y, 0);
      double b = shapeset.get_fn_value(shape_n, x, y, 1);
      double length = sqrt(a*a + b*b);

      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
      glTranslated(x,y,0.0);
//      glScaled(length/max, length/max, 1.0);
      glRotated(atan2(b,a) * 180/PI, 0.0, 0.0, 1.0);


      float* color = get_palette_entry(length/max); // range is 0.0 -- 1.0
      glColor3f(color[0], color[1], color[2]);

      glBegin(GL_LINES);
      glVertex2d(0.0,0.0);
      glVertex2d((length/max) * 2.0/steps,0.0);
      glEnd();
      double pomer = 1.5;
      double v = 2.0/(2.0*steps);
      double z = v / pomer;
      if (length/max < 1e-6)
      {
        glBegin(GL_QUADS);
        glVertex2d( z/2,  z/2);
        glVertex2d( z/2, -z/2);
        glVertex2d(-z/2, -z/2);
        glVertex2d(-z/2,  z/2);
        glEnd();
      }
      else
      {
        glBegin(GL_TRIANGLES);
        glVertex2d(((length/max) * 2.0/steps) + (2.0/3.0) * v,  0.0);
        glVertex2d(((length/max) * 2.0/steps) - (1.0/3.0) * v,  z/2);
        glVertex2d(((length/max) * 2.0/steps) - (1.0/3.0) * v, -z/2);
        glEnd();
      }
/*
      glVertex2d(2.0/steps,0.0);
      glVertex2d(4.0/(3.0*steps),0.01);
      glVertex2d(2.0/steps,0.0);
      glVertex2d(4.0/(3.0*steps),-0.01);
      glEnd();
*/
      glLoadIdentity();
    }
  }

}

void DrawMesh_tri(int steps)
{

  double max = 0.0;
  for (int i = 0; i <= steps; i++)
  {
    for (int j = 0; j <= steps - i; j++)
    {
      double x = -1.0 + i * 2.0/steps;
      double y = -1.0 + j * 2.0/steps;

      double a = shapeset.get_fn_value(shape_n, x, y, 0);
      double b = shapeset.get_fn_value(shape_n, x, y, 1);
      double length = sqrt(a*a + b*b);
      if (length > max) max = length;
    }
  }

  for (int i = 0; i <= steps; i++)
  {
    for (int j = 0; j <= steps - i; j++)
    {
      double x = -1.0 + i * 2.0/steps;
      double y = -1.0 + j * 2.0/steps;

      double a = shapeset.get_fn_value(shape_n, x, y, 0);
      double b = shapeset.get_fn_value(shape_n, x, y, 1);
      double length = sqrt(a*a + b*b);

      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
      glTranslated(x,y,0.0);
//      glScaled(length/max, length/max, 1.0);
      glRotated(atan2(b,a) * 180/PI, 0.0, 0.0, 1.0);


      float* color = get_palette_entry(length/max); // range is 0.0 -- 1.0
      glColor3f(color[0], color[1], color[2]);

      glBegin(GL_LINES);
      glVertex2d(0.0,0.0);
      glVertex2d((length/max) * 2.0/steps,0.0);
      glEnd();
      double pomer = 1.5;
      double v = 2.0/(2.0*steps);
      double z = v / pomer;
      if (length/max < 1e-6)
      {
        glBegin(GL_QUADS);
        glVertex2d( z/2,  z/2);
        glVertex2d( z/2, -z/2);
        glVertex2d(-z/2, -z/2);
        glVertex2d(-z/2,  z/2);
        glEnd();
      }
      else
      {
        glBegin(GL_TRIANGLES);
        glVertex2d(((length/max) * 2.0/steps) + (2.0/3.0) * v,  0.0);
        glVertex2d(((length/max) * 2.0/steps) - (1.0/3.0) * v,  z/2);
        glVertex2d(((length/max) * 2.0/steps) - (1.0/3.0) * v, -z/2);
        glEnd();
      }
/*
      glVertex2d(2.0/steps,0.0);
      glVertex2d(4.0/(3.0*steps),0.01);
      glVertex2d(2.0/steps,0.0);
      glVertex2d(4.0/(3.0*steps),-0.01);
      glEnd();
*/
    }
  }

}


void OnDisplay(void)
{
  glClearColor(1.0, 1.0, 1.0, 0.0);
  glClear(GL_COLOR_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glColor3f(0.0, 0.0, 0.0);
  glLineWidth(1.0);
  glBegin(GL_LINE_LOOP);

  glVertex2d(-1, -1);
  glVertex2d( 1, -1);
  glVertex2d(-1,  1);
  glEnd();

/*
  glVertex2d(-1,  1);
  glVertex2d(-1, -1);
  glVertex2d( 1, -1);
  glVertex2d( 1,  1);
  glEnd();
*/

  glLineWidth(1.0);
 // DrawMesh_quad(30);
  DrawMesh_tri(30);

  glutSwapBuffers();
}



void OnReshape(int w, int h)
{
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  double ratio = (double) w / h;
  const double R = 1.1;
  if (w > h)
    glOrtho(-R*ratio, R*ratio, -R, R, -1, 1);
  else
    glOrtho(-R, R, -R/ratio, R/ratio, -1, 1);
}



void OnKeyDown(unsigned char key, int x, int y)
{
  switch (key)
  {
    case 27: exit(0);

    case 'o':
    {
      shape_n--;
      if (shape_n < 0) shape_n = 0;
      break;
    }

    case 'p':
    {
      shape_n++;
      if (shape_n > shapeset.get_max_index()) shape_n = shapeset.get_max_index();
      break;
    }

    default: return;
  }

  char text[200];
  sprintf(text, "ShapeView - vertical order %d, horizontal order %d, number %d", get_v_order(shapeset.get_order(shape_n)), get_h_order(shapeset.get_order(shape_n)), shape_n);
//  sprintf(text, "ShapeView - order %d, number %d", shapeset.get_order(shape_n), shape_n);
  glutSetWindowTitle(text);
  glutPostRedisplay();
}



/*void OnMouseButton(int button, int state, int x, int y)
{
}*/


/*void OnMouseMotion(int x, int y)
{
}*/




int main(int argc, char* argv[])
{
  // initialization
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);

  // opening window of required size at required position
  glutInitWindowPosition(100, 100);
  glutInitWindowSize(800, 600);
  glutCreateWindow("ShapeView");

  // initialization of event handler
  glutDisplayFunc(OnDisplay);
  glutReshapeFunc(OnReshape);
  glutKeyboardFunc(OnKeyDown);
  //glutMouseFunc(OnMouseButton);
  //glutMotionFunc(OnMouseMotion);
  //glutPassiveMotionFunc(OnMouseMotion);

  glutMainLoop();

  return 0;
}
