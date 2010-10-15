#include "trans.h"
#define countof(a) (sizeof(a)/sizeof(a[0]))

// taken form RefMap.cpp
extern H1ShapesetJacobi ref_map_shapeset;

// vertical edge
double2 e0_pt[] = {
	{ -1,  1 },
	{ -1, -1 }
};

// horz edge
double2 e1_pt[] = {
	{  1,   -1 },
	{ -1,   -1 }
};

// last edge
double2 e2_pt[] = {
	{   -1,  1 },
	{ -0.5,  0.5 },
	{    0,  0 },
	{  0.5, -0.5 },
	{    1, -1 }
};

// modified get_phys_x,y
// @param[in] pt - point to tranforms
// @return transformed points
double2 *transform_element(Element *e, int np, double2 *pt)
{
	Quad2D *quad_2d = &g_quad_2d_std;

	ref_map_shapeset.set_mode(e->get_mode());

	// transform all x coordinates of the integration points
	double2 *tpt = new double2[np];
	memset(tpt, 0, np * sizeof(double2));

	int indices[70];

	// taken from refmap::set_active_element
	// prepare the shapes and coefficients of the reference map
	int j, k = 0;
	for (unsigned int i = 0; i < e->nvert; i++)
		indices[k++] = ref_map_shapeset.get_vertex_index(i);

	int o = e->cm->order;
	for (unsigned int i = 0; i < e->nvert; i++)
		for (j = 2; j <= o; j++)
			indices[k++] = ref_map_shapeset.get_edge_index(i, 0, j);

	if (e->is_quad()) o = H2D_MAKE_QUAD_ORDER(o, o);
	memcpy(indices + k, ref_map_shapeset.get_bubble_indices(o),
		ref_map_shapeset.get_num_bubbles(o) * sizeof(int));

	// taken from RefMap::get_phys_x
	for (int i = 0; i < e->cm->nc; i++)
	{
		for (j = 0; j < np; j++) {
			double fn = ref_map_shapeset.get_fn_value(indices[i], pt[j][0], pt[j][1], 0);
			tpt[j][0] += e->cm->coeffs[i][0] * fn;
			tpt[j][1] += e->cm->coeffs[i][1] * fn;
		}
	}

	return tpt;
}

HERMES_API double2 *transform(Element *e)
{
	double2 *tpt1 = NULL;
	double2 *tpt2 = NULL;
	double2 *tpt3 = NULL;

/*
	tpt = transform_element(e, countof(e0_pt), e0_pt);
	printf("edge #0\n");
	for (int i = 0; i < countof(e0_pt); i++)
		printf("%lf, %lf\n", tpt[i][0], tpt[i][1]);

	tpt = transform_element(e, countof(e1_pt), e1_pt);
	printf("edge #1\n");
	for (int i = 0; i < countof(e1_pt); i++)
		printf("%lf, %lf\n", tpt[i][0], tpt[i][1]);
*/

	tpt3 = transform_element(e, countof(e2_pt), e2_pt);
	//printf("edge #2\n");
	//for (int i = 0; i < countof(e2_pt); i++)
	//	printf("%lf, %lf\n", tpt3[i][0], tpt3[i][1]);
	return tpt3;
}

HERMES_API void element_polygonal_boundary(Element *e, double2 **tp, int *npoints)
{
	double2 *pt;
    int n;
    int d = 2; // number of points on one side (e.g. d >= 2)

    //*tp = transform_element(e, countof(e2_pt), e2_pt);
    if (e->is_triangle()) {
        if (e->is_curved()) {
            d = 10;
            n = 3*(d-1);
            pt = new double2[n];
            double h = 2.0/(d-1);
            int counter = 0;
            for (int i=0; i < d-1; i++) {
                pt[counter][0] = -1+i*h;
                pt[counter][1] = -1;
                counter++;
            }
            for (int i=0; i < d-1; i++) {
                pt[counter][0] = +1-i*h;
                pt[counter][1] = -1+i*h;
                counter++;
            }
            for (int i=0; i < d-1; i++) {
                pt[counter][0] = -1;
                pt[counter][1] = +1-i*h;
                counter++;
            }
            if (counter != n)
                error("Internal error: counter != n.");
            pt = transform_element(e, n, pt);
        } else {
            n = 3;
            pt = new double2[n];
            for (int i=0; i < 3; i++) {
                pt[i][0] = e->vn[i]->x;
                pt[i][1] = e->vn[i]->y;
            }
        }
    } else if (e->is_quad()) {
        if (e->is_curved()) {
            d = 10;
            n = 4*(d-1);
            pt = new double2[n];
            double h = 2.0/(d-1);
            int counter = 0;
            for (int i=0; i < d-1; i++) {
                pt[counter][0] = -1+i*h;
                pt[counter][1] = -1;
                counter++;
            }
            for (int i=0; i < d-1; i++) {
                pt[counter][0] = +1;
                pt[counter][1] = -1+i*h;
                counter++;
            }
            for (int i=0; i < d-1; i++) {
                pt[counter][0] = +1-i*h;
                pt[counter][1] = +1;
                counter++;
            }
            for (int i=0; i < d-1; i++) {
                pt[counter][0] = -1;
                pt[counter][1] = +1-i*h;
                counter++;
            }
            if (counter != n)
                error("Internal error: counter != n.");
            pt = transform_element(e, n, pt);
        } else {
            n = 4;
            pt = new double2[n];
            for (int i=0; i < 4; i++) {
                pt[i][0] = e->vn[i]->x;
                pt[i][1] = e->vn[i]->y;
            }
        }
    } else
        error("Unsupported element.");
    *tp = pt;
    *npoints = n;
}
