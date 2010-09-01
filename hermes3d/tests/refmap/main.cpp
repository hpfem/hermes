// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

//
// Testing reference mappings
//

#include "config.h"
#include <hermes3d.h>
#include <common/trace.h>
#include <common/error.h>

#define ERROR_SUCCESS								0
#define ERROR_FAILURE								-1


int main(int argc, char *argv[]) {
	int ret = ERROR_SUCCESS;
	set_verbose(false);

	if (argc < 2) error("Not enough parameters");

	Mesh mesh;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(argv[1], &mesh)) error("Loading mesh file '%s'\n", argv[1]);

	RefMap refmap(&mesh);

	FOR_ALL_ACTIVE_ELEMENTS(eid, &mesh) {
		printf("Element #%ld\n", eid);
		Element *e = mesh.elements[eid];
		Quad3D *quad = get_quadrature(e->get_mode());

	    refmap.set_active_element(e);

	    // ref. order
	    order3_t ref_ord = refmap.get_ref_order();
	    printf(" - ref. order: %s\n", ref_ord.str());

	    // inv. ref. order
	    order3_t inv_ref_ord = refmap.get_inv_ref_order();
	    printf(" - inv. ref. order: %s\n", inv_ref_ord.str());

	    // const jacobian
	    bool const_jac = refmap.is_jacobian_const();
	    printf(" - const jacobian: %s\n", const_jac ? "yes" : "no");
	    if (const_jac) {
	    	// TODO: implement me
	    }
	    else {
	    	order3_t o3(2, 2, 2);
	    	int np = quad->get_num_points(o3);

	    	double *jac = refmap.get_jacobian(o3);
	    	printf(" - jacobian: ");
	    	for (int i = 0; i < np; i++) printf(" % lf", jac[i]);
	    	printf("\n");

	    	double3x3 *rm = refmap.get_ref_map(o3);
	    	double3x3 *irm = refmap.get_inv_ref_map(o3);
	    	printf(" - refmap | inv ref map\n");
	    	for (int i = 0; i < np; i++) {
	    		printf("  % lf, % lf, % lf | % lf, % lf, % lf\n", rm[i][0][0], rm[i][0][1], rm[i][0][2], irm[i][0][0], irm[i][0][1], irm[i][0][2]);
	    		printf("  % lf, % lf, % lf | % lf, % lf, % lf\n", rm[i][1][0], rm[i][1][1], rm[i][1][2], irm[i][1][0], irm[i][1][1], irm[i][1][2]);
	    		printf("  % lf, % lf, % lf | % lf, % lf, % lf\n", rm[i][2][0], rm[i][2][1], rm[i][2][2], irm[i][2][0], irm[i][2][1], irm[i][2][2]);
	    		printf("\n");
	    	}

	    	double *x = refmap.get_phys_x(o3);
	    	double *y = refmap.get_phys_y(o3);
	    	double *z = refmap.get_phys_z(o3);
	    	printf(" - phys coords [x, y, z]\n");
	    	for (int i = 0; i < np; i++) {
	    		printf("   % lf, % lf, % lf\n", x[i], y[i], z[i]);
	    	}

	    	// faces
	    	for (int iface = 0; iface < e->get_num_faces(); iface++) {
	    		printf(" * Face %d\n", iface);

	    		bool face_const_jac = refmap.is_face_const_jacobian(iface);
	    	    printf("   - const jacobian: %s\n", face_const_jac ? "yes" : "no");
	    	    if (const_jac) {
	    	    }
	    	    else {
	    	    	order2_t o2(2, 2);
			    	int fnp = quad->get_face_num_points(iface, o2);

	    	    	double *face_jac = refmap.get_face_jacobian(iface, o2);
	    	    	printf("   - jacobian: ");
	    	    	for (int i = 0; i < fnp; i++) printf(" % lf", face_jac[i]);
	    	    	printf("\n");

	    	    	double3x3 *face_rm = refmap.get_face_ref_map(iface, o2);
	    	    	double3x3 *face_irm = refmap.get_face_inv_ref_map(iface, o2);
	    	    	printf("   - refmap | inv ref map\n");
	    	    	for (int i = 0; i < fnp; i++) {
	    	    		printf("    % lf, % lf, % lf | % lf, % lf, % lf\n", face_rm[i][0][0], face_rm[i][0][1], face_rm[i][0][2], face_irm[i][0][0], face_irm[i][0][1], face_irm[i][0][2]);
	    	    		printf("    % lf, % lf, % lf | % lf, % lf, % lf\n", face_rm[i][1][0], face_rm[i][1][1], face_rm[i][1][2], face_irm[i][1][0], face_irm[i][1][1], face_irm[i][1][2]);
	    	    		printf("    % lf, % lf, % lf | % lf, % lf, % lf\n", face_rm[i][2][0], face_rm[i][2][1], face_rm[i][2][2], face_irm[i][2][0], face_irm[i][2][1], face_irm[i][2][2]);
	    	    		printf("\n");
	    	    	}

	    	    	Point3D *face_nm = refmap.get_face_normal(iface, o2);
	    	    	printf("   - outer normals:\n");
	    	    	for (int i = 0; i < fnp; i++) {
	    	    		printf("     [% lf, % lf, % lf]\n", face_nm[i].x, face_nm[i].y, face_nm[i].z);
	    	    	}

	    	    	double *fx = refmap.get_face_phys_x(iface, o2);
	    	    	double *fy = refmap.get_face_phys_y(iface, o2);
	    	    	double *fz = refmap.get_face_phys_z(iface, o2);
	    	    	printf("   - phys coords [x, y, z]\n");
	    	    	for (int i = 0; i < fnp; i++) {
	    	    		printf("     % lf, % lf, % lf\n", fx[i], fy[i], fz[i]);
	    	    	}
	    	    }
	    	}

	    	// TODO: edges
	    }

	}


	return ret;
}
