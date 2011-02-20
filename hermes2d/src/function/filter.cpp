// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#include "../h2d_common.h"
#include "filter.h"
#include "../mesh/traverse.h"


//// Filter ////////////////////////////////////////////////////////////////////////////////////////
Filter::Filter(Hermes::vector<MeshFunction*> solutions) : MeshFunction()
{
	this->num = solutions.size();
	if(num > 10)
		error("Attempt to create an instance of Filter with more than 10 MeshFunctions.");
	for(int i = 0; i < this->num; i++)
		this->sln[i] = solutions.at(i);
  this->init();
}

void Filter::init(Hermes::vector<MeshFunction*> solutions)
{
	this->num = solutions.size();
	if(num > 10)
		error("Attempt to create an instance of Filter with more than 10 MeshFunctions.");
	for(int i = 0; i < this->num; i++)
		this->sln[i] = solutions.at(i);
  this->init();
}


void Filter::init()
{
  // construct the union mesh, if necessary
  Mesh* meshes[10];
	for(int i = 0; i < this->num; i++)
		meshes[i] = this->sln[i]->get_mesh();
  mesh = meshes[0];
  unimesh = false;

  for (int i = 1; i < num; i++)
	{
    if (meshes[i] == NULL) {
      warn("You may be initializing a Filter with Solution that is missing a Mesh.");
      error("this->meshes[%d] is NULL in Filter::init().", i);
    }
    if (meshes[i]->get_seq() != mesh->get_seq())
    {
			unimesh = true;
			break;
		}
  }

  if (unimesh)
  {
    Traverse trav;
    trav.begin(num, meshes);
    mesh = new Mesh;
    unidata = trav.construct_union_mesh(mesh);
    trav.finish();
  }

  // misc init
  num_components = 1;
  order = 0;

  for(int i = 0; i < 10; i++)
      tables[i] = new LightArray<LightArray<Node*>*>;

  memset(sln_sub, 0, sizeof(sln_sub));
  set_quad_2d(&g_quad_2d_std);
}


Filter::~Filter()
{
  free();
  if (unimesh)
  {
    delete mesh;
    for (int i = 0; i < num; i++)
      ::free(unidata[i]);
    delete [] unidata;
  }
}


void Filter::set_quad_2d(Quad2D* quad_2d)
{
  MeshFunction::set_quad_2d(quad_2d);
  for (int i = 0; i < num; i++)
    sln[i]->set_quad_2d(quad_2d); // nodup
}


void Filter::set_active_element(Element* e)
{
  MeshFunction::set_active_element(e);
  if (!unimesh)
  {
    for (int i = 0; i < num; i++)
      sln[i]->set_active_element(e); // nodup
    memset(sln_sub, 0, sizeof(sln_sub));
  }
  else
  {
    for (int i = 0; i < num; i++) {
      sln[i]->set_active_element(unidata[i][e->id].e);
      sln[i]->set_transform(unidata[i][e->id].idx);
      sln_sub[i] = sln[i]->get_transform();
    }
  }

  if (tables[cur_quad] != NULL) {
    for(unsigned int k = 0; k < tables[cur_quad]->get_size(); k++)
      if(tables[cur_quad]->present(k)) {
        for(unsigned int l = 0; l < tables[cur_quad]->get(k)->get_size(); l++)
          if(tables[cur_quad]->get(k)->present(l))
            ::free(tables[cur_quad]->get(k)->get(l));
        delete tables[cur_quad]->get(k);
      }
    delete tables[cur_quad];
  }

  tables[cur_quad] = new LightArray<LightArray<Node*>*>;

  sub_tables = tables[cur_quad];
  update_nodes_ptr();

  order = 20; // fixme
}


void Filter::free()
{
  for (int i = 0; i < num; i++)
    if (tables[i] != NULL) {
      for(unsigned int k = 0; k < tables[i]->get_size(); k++) 
        if(tables[i]->present(k)) {
          for(unsigned int l = 0; l < tables[i]->get(k)->get_size(); l++)
            if(tables[i]->get(k)->present(l))
              ::free(tables[i]->get(k)->get(l));
          delete tables[i]->get(k);
        }
      delete tables[i];
      tables[i] = NULL;
    }
}


void Filter::reinit()
{
  free();
  init();
}


void Filter::push_transform(int son)
{
  MeshFunction::push_transform(son);
  for (int i = 0; i < num; i++)
  {
    // sln_sub[i] contains the value sln[i]->sub_idx, which the Filter thinks
    // the solution has, or at least had last time. If the filter graph is
    // cyclic, it could happen that some solutions would get pushed the transform
    // more than once. This mechanism prevents it. If the filter sees that the
    // solution already has a different sub_idx than it thinks it should have,
    // it assumes someone else has already pushed the correct transform. This
    // is actually the case for cyclic filter graphs and filters used in multi-mesh
    // computation.

    if (sln[i]->get_transform() == sln_sub[i])
      sln[i]->push_transform(son);
    sln_sub[i] = sln[i]->get_transform();
  }
}


void Filter::pop_transform()
{
  MeshFunction::pop_transform();
  for (int i = 0; i < num; i++)
  {
    if (sln[i]->get_transform() == sln_sub[i])
      sln[i]->pop_transform();
    sln_sub[i] = sln[i]->get_transform();
  }
}


//// SimpleFilter //////////////////////////////////////////////////////////////////////////////////

SimpleFilter::SimpleFilter(void (*filter_fn)(int n, Hermes::vector<scalar*> values, scalar* result),
		Hermes::vector<MeshFunction*> solutions, Hermes::vector<int> items) : filter_fn(filter_fn)
{
	this->num = solutions.size();
	if(num > 10)
		error("Attempt to create an instance of Filter with more than 10 MeshFunctions.");
	if(items.size() != (unsigned) num)
		if(items.size() > 0)
			error("Attempt to create an instance of SimpleFilter with different supplied number of MeshFunctions than the number of types of data used from them.");

	for(int i = 0; i < this->num; i++)
	{
		this->sln[i] = solutions.at(i);
		if(items.size() > 0)
			this->item[i] = items.at(i);
		else
			this->item[i] =	H2D_FN_VAL;
	}
	this->init();
	init_components();
}

void SimpleFilter::init_components()
{
  bool vec1 = false, vec2 = false;
  for (int i = 0; i < num; i++)
  {
    if (sln[i]->get_num_components() > 1) vec1 = true;
    if ((item[i] & H2D_FN_COMPONENT_0) && (item[i] & H2D_FN_COMPONENT_1)) vec2 = true;
    if (sln[i]->get_num_components() == 1) item[i] &= H2D_FN_COMPONENT_0;
  }
  num_components = (vec1 && vec2) ? 2 : 1;
}


void SimpleFilter::precalculate(int order, int mask)
{
  if (mask & (H2D_FN_DX | H2D_FN_DY | H2D_FN_DXX | H2D_FN_DYY | H2D_FN_DXY))
    error("Filter not defined for derivatives.");

  Quad2D* quad = quads[cur_quad];
  int np = quad->get_num_points(order);
  Node* node = new_node(H2D_FN_VAL, np);

  // precalculate all solutions
  for (int i = 0; i < num; i++)
    sln[i]->set_quad_order(order, item[i]);

  for (int j = 0; j < num_components; j++)
  {
    // obtain corresponding tables
    scalar* tab[10];
    for (int i = 0; i < num; i++)
    {
      int a = 0, b = 0, mask = item[i];
      if (mask >= 0x40) { a = 1; mask >>= 6; }
      while (!(mask & 1)) { mask >>= 1; b++; }
      tab[i] = sln[i]->get_values(num_components == 1 ? a : j, b);
      if (tab[i] == NULL) error("Value of 'item%d' is incorrect in filter definition.", i+1);
    }

		Hermes::vector<scalar*> values;
		for(int i = 0; i < this->num; i++)
			values.push_back(tab[i]);

    // apply the filter
		filter_fn(np, values, node->values[j][0]);
  }

  if(nodes->present(order)) {
    assert(nodes->get(order) == cur_node);
    ::free(nodes->get(order));
  }
  nodes->add(node, order);
  cur_node = node;
}

scalar SimpleFilter::get_pt_value(double x, double y, int it)
{
  if (it & (H2D_FN_DX | H2D_FN_DY | H2D_FN_DXX | H2D_FN_DYY | H2D_FN_DXY))
    error("Filter not defined for derivatives.");
  scalar val[10];
  for (int i = 0; i < num; i++)
    val[i] = sln[i]->get_pt_value(x, y, item[i]);

  scalar result;
	Hermes::vector<scalar*> values;
	for(int i = 0; i < this->num; i++)
		values.push_back(&val[i]);

	// apply the filter
	filter_fn(1, values, &result);

	return result;
}

//// DXDYFilter ////////////////////////////////////////////////////////////////////////////////////


DXDYFilter::DXDYFilter(filter_fn_ filter_fn, Hermes::vector<MeshFunction*> solutions) : Filter(solutions), filter_fn(filter_fn)
{
	init_components();
}

void DXDYFilter::init_components()
{
  num_components = sln[0]->get_num_components();
  for (int i = 1; i < num; i++)
    if (sln[i]->get_num_components() != num_components)
      error("Filter: Solutions do not have the same number of components!");
}

void DXDYFilter::init(filter_fn_ filter_fn, Hermes::vector<MeshFunction*> solutions)
{
  Filter::init(solutions);
  this->filter_fn = filter_fn;
  init_components();
}

void DXDYFilter::precalculate(int order, int mask)
{
  Quad2D* quad = quads[cur_quad];
  int np = quad->get_num_points(order);
  Node* node = new_node(H2D_FN_DEFAULT, np);

  // precalculate all solutions
  for (int i = 0; i < num; i++)
    sln[i]->set_quad_order(order, H2D_FN_DEFAULT);

  for (int j = 0; j < num_components; j++)
  {
    // obtain solution tables
    scalar *val[10], *dx[10], *dy[10];
    for (int i = 0; i < num; i++)
    {
      val[i] = sln[i]->get_fn_values(j);
      dx[i]  = sln[i]->get_dx_values(j);
      dy[i]  = sln[i]->get_dy_values(j);
    }

                Hermes::vector<scalar *> values_vector;
                Hermes::vector<scalar *> dx_vector;
                Hermes::vector<scalar *> dy_vector;

		for(int i = 0; i < this->num; i++)
		{
                        values_vector.push_back(val[i]);
                        dx_vector.push_back(dx[i]);
                        dy_vector.push_back(dy[i]);
		}

    // apply the filter
    filter_fn(np, values_vector, dx_vector, dy_vector, node->values[j][0], node->values[j][1], node->values[j][2]);
  }

  if(nodes->present(order)) {
    assert(nodes->get(order) == cur_node);
    ::free(nodes->get(order));
  }
  nodes->add(node, order);
  cur_node = node;
}


//// predefined simple filters /////////////////////////////////////////////////////////////////////

static void magnitude_fn(int n, Hermes::vector<scalar*> values, scalar* result)
{
  for (int i = 0; i < n; i++)
	{
		result[i] = 0;
		for(unsigned int j = 0; j < values.size(); j++)
			result[i] += sqr(values.at(j)[i]);
		result[i] = sqrt(result[i]);
	}
};

MagFilter::MagFilter(Hermes::vector<MeshFunction*> solutions, Hermes::vector<int> items) : SimpleFilter(magnitude_fn, solutions, items)
{
};

MagFilter::MagFilter(MeshFunction* sln1, int item1)
         : SimpleFilter(magnitude_fn, Hermes::vector<MeshFunction*>(sln1, sln1), Hermes::vector<int>(item1 & H2D_FN_COMPONENT_0, item1 & H2D_FN_COMPONENT_1))
{
  if (sln1->get_num_components() < 2)
    error("The single-argument constructor is intended for vector-valued solutions.");
};


static void difference_fn_2(int n, Hermes::vector<scalar*> values, scalar* result)
{
  for (int i = 0; i < n; i++)
		result[i] = values.at(0)[i] - values.at(1)[i];
};

DiffFilter::DiffFilter(Hermes::vector<MeshFunction*> solutions, Hermes::vector<int> items)
          : SimpleFilter(difference_fn_2, solutions, items) {}


static void sum_fn(int n, Hermes::vector<scalar*> values, scalar* result)
{
  for (int i = 0; i < n; i++)
		for (unsigned int j = 0; j < values.size(); j++)
			result[i] += values.at(j)[i];
};

SumFilter::SumFilter(Hermes::vector<MeshFunction*> solutions, Hermes::vector<int> items) : SimpleFilter(sum_fn, solutions, items) {}


static void square_fn_1(int n, Hermes::vector<scalar*> v1, scalar* result)
{
#ifdef H2D_COMPLEX
  for (int i = 0; i < n; i++)
		result[i] = std::norm(v1.at(0)[i]);
#else
  for (int i = 0; i < n; i++)
    result[i] = sqr(v1.at(0)[i]);
#endif
};

SquareFilter::SquareFilter(Hermes::vector<MeshFunction*> solutions, Hermes::vector<int> items)
          : SimpleFilter(square_fn_1, solutions, items)
{
	if (solutions.size() > 1)
    error("SquareFilter only supports one MeshFunction.");
};


static void real_part_fn_1(int n, Hermes::vector<scalar*> v1, scalar* result)
{
#ifndef H2D_COMPLEX
  memcpy(result, v1.at(0), sizeof(scalar) * n);
#else
  for (int i = 0; i < n; i++)
    result[i] = v1.at(0)[i].real();
#endif
};

RealFilter::RealFilter(Hermes::vector<MeshFunction*> solutions, Hermes::vector<int> items)
          : SimpleFilter(real_part_fn_1, solutions, items)
{
	if (solutions.size() > 1)
		error("RealFilter only supports one MeshFunction.");
};


static void imag_part_fn_1(int n, Hermes::vector<scalar*> v1, scalar* result)
{
#ifndef H2D_COMPLEX
  memset(result, 0, sizeof(scalar) * n);
#else
  for (int i = 0; i < n; i++)
    result[i] = v1.at(0)[i].imag();
#endif
};

ImagFilter::ImagFilter(Hermes::vector<MeshFunction*> solutions, Hermes::vector<int> items)
          : SimpleFilter(imag_part_fn_1, solutions, items)
{
	if (solutions.size() > 1)
		error("RealFilter only supports one MeshFunction.");
};


static void abs_fn_1(int n,  Hermes::vector<scalar*> v1, scalar* result)
{
#ifndef H2D_COMPLEX
  for (int i = 0; i < n; i++)
    result[i] = fabs(v1.at(0)[i]);
#else
  for (int i = 0; i < n; i++)
		result[i] = sqrt(sqr(v1.at(0)[i].real()) + sqr(v1.at(0)[i].imag()));
#endif
};

AbsFilter::AbsFilter(Hermes::vector<MeshFunction*> solutions, Hermes::vector<int> items)
          : SimpleFilter(abs_fn_1, solutions, items)
{
		if (solutions.size() > 1)
		error("RealFilter only supports one MeshFunction.");
};


static void angle_fn_1(int n, Hermes::vector<scalar*> v1, scalar* result)
{
#ifndef H2D_COMPLEX
  for (int i = 0; i < n; i++)
    result[i] = 0.0;
#else
  for (int i = 0; i < n; i++)
    result[i] = atan2( v1.at(0)[i].imag(), v1.at(0)[i].real() );
#endif
};

AngleFilter::AngleFilter(Hermes::vector<MeshFunction*> solutions, Hermes::vector<int> items)
  : SimpleFilter(angle_fn_1, solutions, items)
{
	if (solutions.size() > 1)
		error("RealFilter only supports one MeshFunction.");
};


//// VonMisesFilter ////////////////////////////////////////////////////////////////////////////////

#ifndef H2D_COMPLEX
  #define getval(exp) (exp)
#else
  #define getval(exp) (exp.real())
#endif


void VonMisesFilter::precalculate(int order, int mask)
{
  if (mask & (H2D_FN_DX | H2D_FN_DY | H2D_FN_DXX | H2D_FN_DYY | H2D_FN_DXY))
    error("VonMisesFilter not defined for derivatives.");

  Quad2D* quad = quads[cur_quad];
  int np = quad->get_num_points(order);
  Node* node = new_node(H2D_FN_VAL_0, np);

  sln[0]->set_quad_order(order, H2D_FN_VAL | H2D_FN_DX | H2D_FN_DY);
  sln[1]->set_quad_order(order, H2D_FN_DX | H2D_FN_DY);

  scalar *dudx, *dudy, *dvdx, *dvdy;
  sln[0]->get_dx_dy_values(dudx, dudy);
  sln[1]->get_dx_dy_values(dvdx, dvdy);
  scalar *uval = sln[0]->get_fn_values();
  update_refmap();
  double *x = refmap->get_phys_x(order);


  for (int i = 0; i < np; i++)
  {
    // stress tensor
    double tz = lambda*(getval(dudx[i]) + getval(dvdy[i]));
    double tx = tz + 2*mu*getval(dudx[i]);
    double ty = tz + 2*mu*getval(dvdy[i]);
    if (cyl) tz += 2*mu*getval(uval[i]) / x[i];
    double txy = mu*(getval(dudy[i]) + getval(dvdx[i]));

    // Von Mises stress
    node->values[0][0][i] = 1.0/sqrt(2.0) * sqrt(sqr(tx - ty) + sqr(ty - tz) + sqr(tz - tx) + 6*sqr(txy));
  }

  if(nodes->present(order)) {
    assert(nodes->get(order) == cur_node);
    ::free(nodes->get(order));
  }
  nodes->add(node, order);
  cur_node = node;
}


VonMisesFilter::VonMisesFilter(Hermes::vector<MeshFunction*> solutions, double lambda, double mu,
                               int cyl, int item1, int item2)
       : Filter(solutions)
{
  this->mu = mu;
  this->lambda = lambda;
  this->cyl = cyl;
  this->item1 = item1;
  this->item2 = item2;
}

//// LinearFilter //////////////////////////////////////////////////////////////////////////////////


void LinearFilter::precalculate(int order, int mask)
{
  Quad2D* quad = quads[cur_quad];
  int np = quad->get_num_points(order);
  Node* node = new_node(H2D_FN_DEFAULT, np);

  // precalculate all solutions
  for (int i = 0; i < num; i++)
    sln[i]->set_quad_order(order);

  for (int j = 0; j < num_components; j++)
  {
    // obtain solution tables
    scalar *val[4], *dx[4], *dy[4];
    for (int i = 0; i < num; i++)
    {
      val[i] = sln[i]->get_fn_values(j);
      dx[i]  = sln[i]->get_dx_values(j);
      dy[i]  = sln[i]->get_dy_values(j);
    }
    if (num == 2)
      for (int i = 0; i < np; i++)
      {
        node->values[j][0][i] = tau_frac * (val[1][i] - val[0][i]) + val[1][i];
        node->values[j][1][i] = tau_frac * (dx[1][i]  - dx[0][i])  + dx[1][i];
        node->values[j][2][i] = tau_frac * (dy[1][i]  - dy[0][i])  + dy[1][i];
      }
    else
      for (int i = 0; i < np; i++)
      {
        node->values[j][0][i] = val[0][i];
        node->values[j][1][i] = dx[0][i];
        node->values[j][2][i] = dy[0][i];
      }

  }

  if(nodes->present(order)) {
    assert(nodes->get(order) == cur_node);
    ::free(nodes->get(order));
  }
  nodes->add(node, order);
  cur_node = node;
}


LinearFilter::LinearFilter(MeshFunction* old)
          : Filter(Hermes::vector<MeshFunction*>(old))
 {
   init_components();
 }

LinearFilter::LinearFilter(MeshFunction* older, MeshFunction* old, double tau_frac)
          : Filter(Hermes::vector<MeshFunction*>(older, old))
 {
   this->tau_frac = tau_frac;
   init_components();
 }


void LinearFilter::init_components()
{
  num_components = sln[0]->get_num_components();
  for (int i = 1; i < num; i++)
    if (sln[i]->get_num_components() != num_components)
      error("Filter: Solutions do not have the same number of components!");
}


void LinearFilter::set_active_element(Element* e)
{
  Filter::set_active_element(e);

  order = 0;
  for (int i = 0; i < num; i++)
  {
    int o = sln[i]->get_fn_order();
    if (o > order) order = o;
  }
}
