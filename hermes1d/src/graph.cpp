// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include <string.h>
#include "common.h"
#include "graph.h"


//// Graph /////////////////////////////////////////////////////////////////////////////////////////

Graph::Graph(const char* title, const char* x_axis_name, const char* y_axis_name)
{
  set_captions(title, x_axis_name, y_axis_name);
  logx = logy = false;
  legend = grid = true;
}


void Graph::set_captions(const char* title, const char* x_axis_name, const char* y_axis_name)
{
  this->title = title ? title : "";
  xname = x_axis_name ? x_axis_name : "";
  yname = y_axis_name ? y_axis_name : "";
}


int Graph::add_row(const char* name, const char* color, const char* line, const char* marker)
{
  Row row;
  if (name == NULL) name = "";
  row.name = name;
  row.color = "k";
  row.line = "-";
  row.marker = "";
  
  rows.push_back(row);
  set_row_style(rows.size()-1, color, line, marker);
  return rows.size()-1;
}

void Graph::set_row_style(int row, const char* color, const char* line, const char* marker)
{
  if (!rows.size()) add_row(NULL);
  rows[row].color  = color;
  rows[row].line   = line;
  rows[row].marker = marker;
}


void Graph::add_values(int row, double x, double y)
{
  if (!rows.size()) add_row(NULL);
  if (row < 0 || row >= rows.size()) error("Invalid row number.");
  Values xy = { x, y };
  rows[row].data.push_back(xy);
}


void Graph::add_values(int row, int n, double* x, double* y)
{
  for (int i = 0; i < n; i++)
    add_values(row, x[i], y[i]);
}


void Graph::add_values(int row, int n, double2* xy)
{
  for (int i = 0; i < n; i++)
    add_values(row, xy[i][0], xy[i][1]);
}


void Graph::save_numbered(const char* filename, int number)
{  
  char buffer[1000];
  sprintf(buffer, filename, number);
  save(buffer);
}


//// MatlabGraph ///////////////////////////////////////////////////////////////////////////////////

void MatlabGraph::save(const char* filename)
{
  int i, j, k;
  
  if (!rows.size()) error("No data rows defined.");
  
  FILE* f = fopen(filename, "w");
  if (f == NULL) error("Error writing to %s", filename);
  
  if (!logx && !logy)
    fprintf(f, "plot(");
  else if (logx && !logy)
    fprintf(f, "semilogx(");
  else if (!logx && logy)
    fprintf(f, "semilogy(");
  else
    fprintf(f, "loglog(");

  for (i = 0; i < rows.size(); i++)
  {
    fprintf(f, "[");
    int rsize = rows[i].data.size();
    for (k = 0; k < 2; k++)
    {
      for (j = 0; j < rsize; j++)
      {
        fprintf(f, "%.14g", k ? rows[i].data[j].y : rows[i].data[j].x);
        if (j < rsize-1) fprintf(f, ", ");
      }
      fprintf(f, (!k ? "], [" : "], '"));
    }
    fprintf(f, "%s%s%s'", rows[i].color.c_str(), rows[i].line.c_str(), rows[i].marker.c_str());
    if (i < rows.size()-1) fprintf(f, ", ");
  }
  fprintf(f, ");\n");
  
  if (title.length()) fprintf(f, "title('%s');\n", title.c_str());
  if (xname.length()) fprintf(f, "xlabel('%s');\n", xname.c_str());
  if (yname.length()) fprintf(f, "ylabel('%s');\n", yname.c_str());

  if (legend && (rows.size() > 1 || rows[0].name.length()))
  {
    fprintf(f, "legend(");
    for (i = 0; i < rows.size(); i++)
    {
      fprintf(f, "'%s'", rows[i].name.c_str());
      if (i < rows.size()-1) fprintf(f, ", ");
    }
    fprintf(f, ");\n");
  }
  else
    fprintf(f, "legend off;\n");
  
  fprintf(f, "grid %s;\n", grid ? "on" : "off");

  info("Graph saved. Run the file '%s' in Matlab.", filename);
  fclose(f);
}


//// GnuplotGraph //////////////////////////////////////////////////////////////////////////////////

static void get_style_types(std::string line, std::string mark, std::string col, int& lt, int& pt, int& ct)
{
  if      (line == "-")  lt = 1; // solid
  else if (line == ":")  lt = 4; // dotted
  else if (line == "-.") lt = 5; // dash dot
  else if (line == "--") lt = 2; // dashed
  else lt = 1;

  if      (mark == ".") pt = 7;  // full circle
  else if (mark == "o") pt = 6;  // empty circle
  else if (mark == "O") pt = 7;  // full circle
  else if (mark == "x") pt = 2;  // cross
  else if (mark == "+") pt = 1;  // cross
  else if (mark == "*") pt = 3;  // star
  else if (mark == "s") pt = 4;  // empty square
  else if (mark == "S") pt = 5;  // full square
  else if (mark == "d") pt = 10; // empty diamond
  else if (mark == "D") pt = 11; // full diamond
  else if (mark == "v") pt = 12; // empty triangle down
  else if (mark == "V") pt = 13; // full triangle down
  else if (mark == "^") pt = 9;  // full triangle up
  else if (mark == "<") pt = 12; // empty triangle down
  else if (mark == ">") pt = 8;  // empty triangle up
  else if (mark == "p") pt = 14; // empty pentagon
  else if (mark == "P") pt = 15; // full pentagon
  else pt = 0;

  if      (col == "k") ct = -1;  // black
  else if (col == "b") ct = 3;   // blue
  else if (col == "g") ct = 2;   // green
  else if (col == "c") ct = 5;   // cyan
  else if (col == "m") ct = 4;   // magenta
  else if (col == "y") ct = 6;   // yellow
  else if (col == "r") ct = 1;   // red
  else ct = -1;
}


void GnuplotGraph::save(const char* filename)
{
  int i, j, k;
  
  if (!rows.size()) error("No data rows defined.");

  FILE* f = fopen(filename, "w");
  if (f == NULL) error("Error writing to %s", filename);
  
  fprintf(f, "set terminal postscript eps enhanced\n");

  int len = strlen(filename);
  char outname[len + 10];
  strcpy(outname, filename);
  char* slash = strrchr(outname, '/');
  if (slash != NULL) strcpy(outname, ++slash);
  char* dot = strrchr(outname, '.');
  if (dot != NULL && dot > outname) *dot = 0;
  strcat(outname, ".eps");
  
  fprintf(f, "set output '%s'\n", outname);  
  
  fprintf(f, "set size 0.8, 0.8\n");

  if (logx && !logy)
    fprintf(f, "set logscale x\n");
  else if (!logx && logy)
    fprintf(f, "set logscale y\n");
  else if (logx && logy)
  {
    fprintf(f, "set logscale x\n");
    fprintf(f, "set logscale y\n");
  }

  if (grid) fprintf(f, "set grid\n");

  if (title.length()) fprintf(f, "set title '%s'\n", title.c_str());
  if (xname.length()) fprintf(f, "set xlabel '%s'\n", xname.c_str());
  if (yname.length()) fprintf(f, "set ylabel '%s'\n", yname.c_str());

  fprintf(f, "plot");
  for (i = 0; i < rows.size(); i++)
  {
    int ct, lt, pt;
    get_style_types(rows[i].line, rows[i].marker, rows[i].color, lt, pt, ct);
 
    if (lt == 0) 
      fprintf(f, " '-' w p pointtype %d title '%s' ", pt, rows[i].name.c_str());
    else if (ct < 0)
      fprintf(f, " '-' w lp linetype %d pointtype %d title '%s' ", lt, pt, rows[i].name.c_str());
    else
      fprintf(f, " '-' w lp linecolor %d linetype %d pointtype %d title '%s' ", ct, lt, pt, rows[i].name.c_str());

    if (i < rows.size() - 1) fprintf(f, ", ");
  }
  fprintf(f,"\n");
  
  for (i = 0; i < rows.size(); i++)
  {
    int rsize = rows[i].data.size();
    for (j = 0; j < rsize; j++)
      fprintf(f, "%.14g  %.14g\n", rows[i].data[j].x, rows[i].data[j].y);
    fprintf(f, "e\n");
  }  

  fprintf(f, "set terminal x11\n");
  printf("Type 'gnuplot %s' to get convergence graph as .eps file.\n", filename);
  fclose(f);
}
