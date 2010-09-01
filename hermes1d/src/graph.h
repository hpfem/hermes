// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef __HERMES1D_GRAPH_H
#define __HERMES1D_GRAPH_H

#include <vector>

///  Graph is a utility class storing a simple XY graph (eg., a convergence graph).
///  One or more data rows can be defined by calling add_row(). The actual data
///  are added to the rows by calling add_values(). The resulting graph is saved
///  to a file by calling save().
///
///  Please note that this is a base class that cannot be instantiated.
///  Use MatlabGraph or GnuplotGraph instead.
///
class Graph
{
public:

  Graph(const char* title = NULL, const char* x_axis_name = NULL, const char* y_axis_name = NULL);

  void set_captions(const char* title = NULL, const char* x_axis_name = NULL, const char* y_axis_name = NULL);

  void set_log_x(bool log = true) { logx = log; }
  void set_log_y(bool log = true) { logy = log; }

  void show_legend(bool show = true) { legend = show; }
  void show_grid(bool show = true) { grid = show; }

  int  add_row(const char* name = NULL, const char* color = "k", const char* line = "-", const char* marker = "");
  void set_row_style(int row, const char* color = "k", const char* line = "-", const char* marker = "");

  void add_values(int row, double x, double y);
  void add_values(int row, int n, double* x, double* y);
  void add_values(int row, int n, double2* xy);

  virtual void save(const char* filename) = 0;
  void save_numbered(const char* filename, int number);
  
  // todo: clear

protected:
  
  std::string title, xname, yname;
  bool logx, logy, legend, grid;
  
  struct Values
  {
    double x, y;
  };

  struct Row 
  {
    std::string name, color, line, marker;
    std::vector<Values> data;
  };

  std::vector<Row> rows;

};


///  Outputs a MATLAB graph.
///
class MatlabGraph : public Graph
{
public:
  
  MatlabGraph(const char* title = NULL, const char* x_axis_name = NULL, const char* y_axis_name = NULL)
      : Graph(title, x_axis_name, y_axis_name) {}
        
  virtual void save(const char* filename);

};


///  Outputs a GNUPLOT graph.
///
class GnuplotGraph : public Graph
{
public:
  
  GnuplotGraph(const char* title = NULL, const char* x_axis_name = NULL, const char* y_axis_name = NULL)
       : Graph(title, x_axis_name, y_axis_name) {}
        
  virtual void save(const char* filename);

};


#endif
