// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "legendre.h"

int legendre_order_1d[] = {
{% for f in functions %}
{{ f.id }},{% endfor %}
};

{% for f in functions %}
static double legendre_fn_{{ f.id }}(double _x) {
    long double x = _x;
    return {{ f.expr }};
}
{% endfor %}

shape_fn_t legendre_fn_tab_1d[] = {
{% for f in functions %}
legendre_fn_{{ f.id }},{% endfor %}
};


{% for f in functions %}
static double legendre_der_{{ f.id }}(double _x) {
    long double x = _x;
    return {{ f.expr_diff }};
}
{% endfor %}

shape_fn_t legendre_der_tab_1d[] = {
{% for f in functions %}
legendre_der_{{ f.id }},{% endfor %}
};
