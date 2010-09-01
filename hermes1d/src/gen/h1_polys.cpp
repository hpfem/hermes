// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "h1_polys.h"

{% for f in functions %}
static double h1_polys_fn_{{ f.id }}(double x) { return {{ f.expr }}; }
{% endfor %}

shape_fn_t h1_polys_fn_tab_1d[] = {
{% for f in functions %}
h1_polys_fn_{{ f.id }},{% endfor %}
};


{% for f in functions %}
static double h1_polys_der_{{ f.id }}(double x) { return {{ f.expr_diff }}; }
{% endfor %}

shape_fn_t h1_polys_der_tab_1d[] = {
{% for f in functions %}
h1_polys_der_{{ f.id }},{% endfor %}
};

