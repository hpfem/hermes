// MeshEdit: this is a quick hack which allows moving vertices in a mesh with the mouse.
//           ...don't expect much!

#include "hermes2d.h"
#include <GL/freeglut.h>


static const char* mesh_file_name;

class MeshEdit;
static MeshEdit* instance;


//// MeshEdit //////////////////////////////////////////////////////////////////////////////////////

class MeshEdit : public MeshView
{
public:

  MeshEdit();

  void run(Mesh* mesh)
  {
    this->mesh = mesh;
    show(mesh);
  }

protected:

  Mesh* mesh;
  int foc_node;
  Nurbs* foc_nurbs;
  int nurbs_node;
  Element* nurbs_e;
  bool dragging_node;
  bool hide;
  bool b_lines;

  double untransform_x(double x) { return (x - trans_x - center_x) / scale; }
  double untransform_y(double y) { return (center_y - y - trans_y) / scale; }

  void draw_circle(double x, double y, double r, int n);
  void draw_nurbs(Nurbs* nurbs);
  void hit_test_nodes(int mx, int my);
  void move_vertex_node(Node* node, double x, double y);
  void update_mesh();
  void draw_template();

  virtual void on_create();
  virtual void on_display();
  virtual void on_mouse_move(int x, int y);
  virtual void on_left_mouse_down(int x, int y);
  virtual void on_left_mouse_up(int x, int y);
  virtual void on_key_down(unsigned char key, int x, int y);
  virtual void on_mouse_wheel(int wheel, int dir, int x, int y);

  friend void on_mouse_wheel_stub(int, int, int, int);
};


//// implementation ////////////////////////////////////////////////////////////////////////////////

MeshEdit::MeshEdit()
        : MeshView("MeshEdit")
{
  foc_node = -1;
  foc_nurbs = NULL;
  b_markers = false;
  dragging_node = false;
  hide = false;
  b_ids = false;
  b_lines = false;
}


void on_mouse_wheel_stub(int wheel, int dir, int x, int y)
  {  instance->on_mouse_wheel(wheel, dir, x, y); }

void MeshEdit::on_create()
{
  glutMouseWheelFunc(on_mouse_wheel_stub);
}


void MeshEdit::draw_circle(double x, double y, double r, int n)
{
  for (int i = 0; i < n; i++)
  {
    double a = 2.0 * 3.141592654 / n * i;
    glVertex2d(x + r*cos(a), y + r*sin(a));
  }
}


void MeshEdit::draw_nurbs(Nurbs* nurbs)
{
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glColor4f(0, 0.7, 0, 0.4);
  glBegin(GL_LINE_STRIP);
  for (int i = 0; i < nurbs->np; i++)
    glVertex2d(transform_x(nurbs->pt[i][0]), transform_y(nurbs->pt[i][1]));
  glEnd();
  glDisable(GL_BLEND);

  for (int i = 1; i < nurbs->np-1; i++)
  {
    glBegin(GL_POLYGON);
    draw_circle(transform_x(nurbs->pt[i][0]), transform_y(nurbs->pt[i][1]), 4, 10);
    glEnd();
  }
}


void MeshEdit::on_display()
{
  glPolygonMode(GL_FRONT_AND_BACK, b_lines ? GL_LINE : GL_FILL);
  MeshView::on_display();
  draw_template();
  if (hide) return;
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  // draw nurbs control points
  Element* e;
  for_all_base_elements(e, mesh)
    if (e->cm != NULL)
      for (int i = 0; i < e->nvert; i++)
        if (e->cm->nurbs[i])
          draw_nurbs(e->cm->nurbs[i]);

  // draw nurbs focus
  if (foc_nurbs != NULL)
  {
    double x = transform_x(foc_nurbs->pt[nurbs_node][0]);
    double y = transform_y(foc_nurbs->pt[nurbs_node][1]);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glBegin(GL_POLYGON);
    draw_circle(x, y, 7, 20);
    glEnd();
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    // print weight
    char text[200];
    glColor3f(0, 0, 0);
    sprintf(text, "%g", foc_nurbs->pt[nurbs_node][2]);
    draw_text(x, y-17, text, 0);
  }

  // draw vertex nodes
  glColor3f(1, 0, 0);
  Node* node;
  for_all_vertex_nodes(node, mesh)
  {
    glBegin(GL_POLYGON);
    draw_circle(transform_x(node->x), transform_y(node->y), 4, 10);
    glEnd();
  }

  // draw node focus
  if (foc_node >= 0)
  {
    node = mesh->get_node(foc_node);
    double x = transform_x(node->x);
    double y = transform_y(node->y);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glBegin(GL_POLYGON);
    draw_circle(x, y, 7, 20);
    glEnd();
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    // print weight
    char text[200];
    glColor3f(0, 0, 0);
    sprintf(text, "%g,%g", node->x, node->y);
    draw_text(x, y-17, text, 0);
  }
}


void MeshEdit::hit_test_nodes(int mx, int my)
{
  foc_node = -1;
  foc_nurbs = NULL;

  // hit-test vertex nodes
  Node* node;
  for_all_vertex_nodes(node, mesh)
  {
    if (sqrt(sqr(transform_x(node->x) - mx) + sqr(transform_y(node->y) - my)) < 7.0)
    {
      foc_node = node->id;
      return;
    }
  }

  // hit-test nurbs nodes
  Element* e;
  for_all_base_elements(e, mesh)
    if (e->cm != NULL)
      for (int i = 0; i < e->nvert; i++)
      {
        Nurbs* nurbs = e->cm->nurbs[i];
        if (nurbs != NULL)
        {
          for (int j = 1; j < nurbs->np-1; j++)
            if (sqrt(sqr(transform_x(nurbs->pt[j][0]) - mx) + sqr(transform_y(nurbs->pt[j][1]) - my)) < 7.0)
            {
              foc_nurbs = nurbs;
              nurbs_node = j;
              nurbs_e = e;
              return;
            }
        }
      }
}


void MeshEdit::move_vertex_node(Node* node, double x, double y)
{
  node->x = x;
  node->y = y;

  // handle curved elements
  Element* e;
  for_all_base_elements(e, mesh)
  {
    if (e->cm != NULL)
    {
      for (int i = 0; i < e->nvert; i++)
      {
        if (e->vn[i] == node)
        {
          if (e->cm->nurbs[i])
          {
            e->cm->nurbs[i]->pt[0][0] = x;
            e->cm->nurbs[i]->pt[0][1] = y;
          }
          int j = e->prev_vert(i);
          if (e->cm->nurbs[j])
          {
            int n = e->cm->nurbs[j]->np - 1;
            e->cm->nurbs[j]->pt[n][0] = x;
            e->cm->nurbs[j]->pt[n][1] = y;
          }
        }
      }
      e->cm->update_refmap_coefs(e);
    }
  }
}


void MeshEdit::update_mesh()
{
  Solution sln;
  sln.set_zero(mesh);
  lin.process_solution(&sln, H2D_FN_VAL_0, H2D_EPS_HIGH);
  post_redisplay();
}


void MeshEdit::on_mouse_move(int x, int y)
{
  if (dragging_node)
  {
    if (foc_node >= 0)
    {
      Node* node = mesh->get_node(foc_node);
      move_vertex_node(node, untransform_x(x), untransform_y(y));
      update_mesh();
    }
    else if (foc_nurbs != NULL)
    {
      foc_nurbs->pt[nurbs_node][0] = untransform_x(x);
      foc_nurbs->pt[nurbs_node][1] = untransform_y(y);
      nurbs_e->cm->update_refmap_coefs(nurbs_e);
      update_mesh();
    }
  }
  else
  {
    MeshView::on_mouse_move(x, y);

    if (!hide)
    {
      int old_foc = foc_node;
      Nurbs* old_nurbs = foc_nurbs;
      hit_test_nodes(x, y);
      if (old_foc != foc_node || old_nurbs != foc_nurbs) post_redisplay();
    }
  }
}


void MeshEdit::on_left_mouse_down(int x, int y)
{
  if (foc_node >= 0 || foc_nurbs != NULL)
    dragging_node = true;
  else
    MeshView::on_left_mouse_down(x, y);
}


void MeshEdit::on_left_mouse_up(int x, int y)
{
  dragging_node = false;
  MeshView::on_left_mouse_up(x, y);
}


void MeshEdit::on_key_down(unsigned char key, int x, int y)
{
  switch (key)
  {
    case 's':
      mesh->save(mesh_file_name);
      info("Mesh saved.");
      break;

    case 'n':
      hide = !hide;
      post_redisplay();
      break;

    case 'm':
      break;

    case 'l':
      b_lines = !b_lines;
      post_redisplay();
      break;

    default:
      MeshView::on_key_down(key, x, y);
      break;
  }
}


void MeshEdit::on_mouse_wheel(int wheel, int dir, int x, int y)
{
  if (foc_nurbs != NULL)
  {
    foc_nurbs->pt[nurbs_node][2] += (dir > 0) ? 0.05 : -0.05;
    if (foc_nurbs->pt[nurbs_node][2] < 1e-10)
      foc_nurbs->pt[nurbs_node][2] = 0.0;
    nurbs_e->cm->update_refmap_coefs(nurbs_e);
    update_mesh();
  }
}


//// template //////////////////////////////////////////////////////////////////////////////////////

void MeshEdit::draw_template()
{
  // here you can place an arbitrary drawing, such as an airfoil which you want to fit your mesh to

  glColor3f(0, 0, 0);
  glBegin(GL_LINE_STRIP);
  for (double x = 0.0, dx = 0.0; x <= 1.000000001; dx += 0.00001, x += dx)
    glVertex2d(transform_x(x),
               transform_y((0.12/0.20)*(0.29690*sqrt(x) - 0.12810*x - 0.35160*x*x + 0.28430*x*x*x - 0.10150*x*x*x*x)));
  glEnd();
}


//// main //////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
  if (argc != 2)
  {
     info("Usage: meshedit mesh-filename");
     exit(1);
  }

  verbose_mode = false;
  mesh_file_name = argv[1];

  Mesh mesh;
  mesh.load(mesh_file_name);

  MeshEdit me;
  instance = &me;
  me.run(&mesh);

  View::wait();
  return 0;
}
