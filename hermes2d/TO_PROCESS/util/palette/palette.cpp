#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

using namespace std;

template<class T> inline T sqr(T x) { return x*x; }

const int WIDTH = 1000;
const int HEIGHT = 200;
const int EDGE = 20;

#define BASE(n) (EDGE + n*(HEIGHT+EDGE) + HEIGHT+1)


struct POINT
{
  int x, y;
  POINT(int x, int y) : x(x), y(y) {}
};

vector<POINT> points[3];

int n = -1, m = -1;
bool dragging = false;



double get_value(int i, double x)
{
  if (x <= points[i][0].x) return (double) points[i][0].y / HEIGHT;
  int last = points[i].size()-1;
  if (x >= points[i][last].x) return (double) points[i][last].y / HEIGHT;
  int j;
  for (j = 0; points[i][j+1].x < x; j++);
  double d = points[i][j+1].x - points[i][j].x;
  double w1 = (points[i][j+1].x - x) / d;
  double w2 = (x - points[i][j].x) / d;
  return (points[i][j].y * w1 + points[i][j+1].y * w2) / HEIGHT;
}


void on_display(void)
{
  glClearColor(1.0, 1.0, 1.0, 0.0);
  glClear(GL_COLOR_BUFFER_BIT);

  glBegin(GL_LINES);
  glVertex2d(EDGE + WIDTH/2, 0);
  glVertex2d(EDGE + WIDTH/2, 1000);
  glEnd();

  for (int i = 0; i < 3; i++)
  {
    glColor3f(0.96, 0.96, 0.96);
    glBegin(GL_QUADS);
    glVertex2d(EDGE, EDGE + i*(HEIGHT+EDGE));
    glVertex2d(EDGE, BASE(i));
    glVertex2d(EDGE + WIDTH+1, BASE(i));
    glVertex2d(EDGE + WIDTH+1, EDGE + i*(HEIGHT+EDGE));
    glEnd();

    glColor3f(0, 0, 0);
    glBegin(GL_LINE_STRIP);
    glVertex2d(EDGE, EDGE + i*(HEIGHT+EDGE));
    glVertex2d(EDGE, BASE(i));
    glVertex2d(EDGE + WIDTH+1, BASE(i));
    glEnd();

    if (i == 0) glColor3f(1.0, 0.0, 0.0);
    else if (i == 1) glColor3f(0.0, 0.8, 0.0);
    else glColor3f(0.0, 0.0, 1.0);

    glBegin(GL_LINE_STRIP);
    for (int j = 0; j < points[i].size(); j++)
      glVertex2d(EDGE + points[i][j].x, BASE(i) - points[i][j].y);
    glEnd();

    glPointSize(5.0);
    glBegin(GL_POINTS);
    for (int j = 0; j < points[i].size(); j++)
      glVertex2d(EDGE + points[i][j].x, BASE(i) - points[i][j].y);
    glEnd();
  }

  glBegin(GL_QUAD_STRIP);
  for (int i = 0; i <= WIDTH; i += 2)
  {
    glColor3d(get_value(0, i), get_value(1, i), get_value(2, i));
    glVertex2d(EDGE + i, EDGE + 3*(HEIGHT+EDGE));
    glVertex2d(EDGE + i, BASE(3));
  }
  glEnd();

  /*glBegin(GL_QUADS);
  glEnd();*/

  glutSwapBuffers();
}



void on_reshape(int w, int h)
{
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0, w, h-1, -1, -10, 10);
}


void load_points()
{
  int i, j, n, x, y;
  FILE* f = fopen("points.txt", "r");
  if (f == NULL) return;
  for (i = 0; i < 3; i++)
  {
    points[i].clear();
    fscanf(f, "%d", &n);
    for (j = 0; j < n; j++)
    {
      fscanf(f, "%d %d", &x, &y);
      points[i].push_back(POINT(x, y));
    }
  }
  fclose(f);
}


void save_points()
{
  FILE* f = fopen("points.txt", "w");
  for (int i = 0; i < 3; i++)
  {
    fprintf(f, "%d\n", points[i].size());
    for (int j = 0; j < points[i].size(); j++)
      fprintf(f, "%d %d\n", points[i][j].x, points[i][j].y);
    fprintf(f, "\n");
  }
  fclose(f);
}


void save_data()
{
  FILE* f = fopen("view_data.cpp", "w");
  const int num = 256;
  fprintf(f, "\nconst int num_pal_entries = %d;\n\n"
             "const float palette_data[num_pal_entries+1][3] = \n{\n  ", num);
  for (int i = 0; i <= num; i++)
  {
    double t = (double) i / num;
    double x = t * WIDTH;
    fprintf(f, "{ %0.4lf, %0.4lf, %0.4lf }", get_value(0, x), get_value(1, x), get_value(2, x));
    if (i < num) fprintf(f, ", ");
    if (i > 0 && !((i+1)%4)) fprintf(f, "\n  ");
  }
  fprintf(f, "\n};\n\n");
  fclose(f);
}


void on_key_down(unsigned char key, int x, int y)
{
  switch (key)
  {
    case 27:
      save_points();
      exit(1);

    case 's':
      save_points();
      break;

    case 'd':
      save_data();
      break;
  }
}


void on_mouse_motion(int x, int y)
{
  if (!dragging)
  {
    n = -1;
    if (y <= EDGE + HEIGHT + EDGE/2) n = 0;
    else if (y <= 2*(EDGE + HEIGHT) + EDGE/2) n = 1;
    else n = 2;

    m = -1;
    for (int i = 0; i < points[n].size(); i++)
    {
      if (sqrt(sqr(points[n][i].x + EDGE - x) + sqr(BASE(n) - points[n][i].y - y)) <= 5.0)
      {
        m = i; break;
      }
    }
  }
  else
  {
    x -= EDGE;
    y = BASE(n) - y;

    if (x < 0) x = 0;
    else if (x > WIDTH) x = WIDTH;
    if (y < 0) y = 0;
    else if (y > HEIGHT) y = HEIGHT;

    if (m > 0 && m < points[n].size()-1)
    {
      if (x <= points[n][m-1].x) x = points[n][m-1].x + 1;
      if (x >= points[n][m+1].x) x = points[n][m+1].x - 1;
      points[n][m].x = x;
    }
    points[n][m].y = y;
    glutPostRedisplay();
  }
}


void on_mouse_button(int button, int state, int x, int y)
{
  if (button == GLUT_LEFT_BUTTON)
  {
    if (state == GLUT_DOWN)
    {
      if (m != -1)
      {
        dragging = true;
        on_mouse_motion(x, y);
      }
      else
      {
        x -= EDGE;
        y = BASE(n) - y;
        if (x < 0 || x > WIDTH || y < 0 || y > HEIGHT) return;
        int i;
        for (i = 0; points[n][i+1].x < x; i++);
        if (points[n][i].x >= x || points[n][i+1].x <= x) return;
        points[n].insert(points[n].begin() + i+1, POINT(x, y));
        dragging = true;
        m = i+1;
        glutPostRedisplay();
      }
    }
    else
      dragging = false;
  }
  else
  {
    if (state == GLUT_DOWN && m != -1)
    {
      points[n].erase(points[n].begin() + m);
      m = -1;
      glutPostRedisplay();
    }
  }
}


int main(int argc, char* argv[])
{
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);

  glutInitWindowPosition(100, 100);
  glutInitWindowSize(2*EDGE + WIDTH, 5*EDGE + 4*HEIGHT);
  glutCreateWindow("Palette Editor");

  glutDisplayFunc(on_display);
  glutReshapeFunc(on_reshape);
  glutKeyboardFunc(on_key_down);
  glutMouseFunc(on_mouse_button);
  glutMotionFunc(on_mouse_motion);
  glutPassiveMotionFunc(on_mouse_motion);

  for (int i = 0; i < 3; i++)
  {
    points[i].push_back(POINT(0, HEIGHT/2));
    points[i].push_back(POINT(WIDTH, HEIGHT/2));
  }

  load_points();
  glutMainLoop();

  return 0;
}
