// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "lobatto.h"

int lobatto_order_1d[] = {
1, {% for f in functions[1:] %}
{{ f.id }},{% endfor %}
};

{% for f in functions %}
static double lobatto_fn_{{ f.id }}(double _x) {
    long double x = _x;
    return {{ f.expr }};
}
{% endfor %}

shape_fn_t lobatto_fn_tab_1d[] = {
{% for f in functions %}
lobatto_fn_{{ f.id }},{% endfor %}
};


{% for f in functions %}
static double lobatto_der_{{ f.id }}(double _x) {
    long double x = _x;
    return {{ f.expr_diff }};
}
{% endfor %}

shape_fn_t lobatto_der_tab_1d[] = {
{% for f in functions %}
lobatto_der_{{ f.id }},{% endfor %}
};
