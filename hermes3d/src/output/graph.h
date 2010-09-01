// This file is part of Hermes3D.
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <vector>
#include <string>

///  Abstract class for graphs
///
///  Graph is a utility class storing a simple XY graph (eg. a convergence graph).
///  One or more data rows can be defined by calling add_row(). The actual data
///  are added to the rows by calling add_values(). The resulting graph is saved
///  to a file by calling save().
///
///  Please note that this is a base class that cannot be instantiated.
///  Use MatlabGraph or GnuplotGraph instead.
///
class Graph {
public:
	Graph(const char *title = NULL, const char *x_axis_name = NULL, const char *y_axis_name = NULL);
	virtual ~Graph() { }

	void set_captions(const char *title = NULL, const char *x_axis_name = NULL, const char *y_axis_name = NULL);

	void set_key(const char *key) { this->key = key; }

	void set_log_x(bool log = true) { logx = log; }
	void set_log_y(bool log = true) { logy = log; }

	void enable_legend(bool enable = true) { legend = enable; }
	void enable_grid(bool enable = true) { grid = enable; }

	int add_row(const char *name = NULL, const char *color = "k", const char *line = "-", const char *marker = "");
	void set_row_style(int row, const char *color = "k", const char *line = "-", const char *marker = "");

	void add_value(int row, double x, double y);
	void add_values(int row, int n, double *x, double *y);
	void add_values(int row, int n, double2 *xy);

	virtual void save(const char *filename = NULL) = 0;
	void save_numbered(const char *filename, int number);

	// todo: clear

protected:
	std::string title, xname, yname;
	bool logx, logy, legend, grid;
	const char *key;

	struct Values {
		double x, y;
	};

	struct Row {
		std::string name, color, line, marker;
		std::vector<Values> data;
	};

	std::vector<Row> rows;
};

///  Outputs a MATLAB graph.
///
class MatlabGraph : public Graph {
public:
	MatlabGraph(const char *title = NULL, const char *x_axis_name = NULL, const char *y_axis_name = NULL) :
		Graph(title, x_axis_name, y_axis_name) {
	}

	virtual void save(const char *filename = NULL);
};

///  Outputs a GNUPLOT graph.
///
class GnuplotGraph : public Graph {
public:
	enum Type {
		Line,
		Column
	};

	GnuplotGraph(Type type = Line, const char *title = NULL, const char *x_axis_name = NULL, const char *y_axis_name = NULL) :
		Graph(title, x_axis_name, y_axis_name) {
		this->type = type;
	}

	virtual void save(const char *filename = NULL);

	void set_type(Type type) { this->type = type; }

protected:
	Type type;

	void save_line_graph(FILE *f);
	void save_column_graph(FILE *f);
};

#endif
