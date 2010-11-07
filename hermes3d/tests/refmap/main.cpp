#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "config.h"
//#include <getopt.h>
#include <hermes3d.h>

// Testing reference mappings.

int main(int argc, char **args)
{
  // Test variable.
  int success_test = 1;

	if (argc < 2) error("Not enough parameters.");

	// Load the mesh.
  Mesh mesh;
  H3DReader mloader;
  if (!mloader.load(args[1], &mesh)) error("Loading mesh file '%s'.", args[1]);

	RefMap refmap(&mesh);

	FOR_ALL_ACTIVE_ELEMENTS(eid, &mesh) {
		info("Element #%ld", eid);
		Element *e = mesh.elements[eid];
		Quad3D *quad = get_quadrature(e->get_mode());

	    refmap.set_active_element(e);

	    // Reference mapping order.
	    Ord3 ref_ord = refmap.get_ref_order();
	    info(" - ref. order: %s", ref_ord.str());

	    // Inverse reference mapping order.
	    Ord3 inv_ref_ord = refmap.get_inv_ref_order();
	    info(" - inv. ref. order: %s", inv_ref_ord.str());

	    // Constant Jacobian.
      bool const_jac = refmap.is_jacobian_const();
	    info(" - const jacobian: %s", const_jac ? "yes" : "no");
	    if (const_jac) {
	    }
	    else {
	    	Ord3 o3(2, 2, 2);
	    	int np = quad->get_num_points(o3);
        
	    	double *jac = refmap.get_jacobian(np, quad->get_points(o3));
	    	info(" - jacobian: ");
	    	for (int i = 0; i < np; i++) info(" % lf", jac[i]);
	    	info("");

	    	double3x3 *rm = refmap.get_ref_map(np, quad->get_points(o3));
	    	double3x3 *irm = refmap.get_inv_ref_map(np, quad->get_points(o3));
	    	info(" - refmap | inv ref map");
	    	for (int i = 0; i < np; i++) {
	    		info("  % lf, % lf, % lf | % lf, % lf, % lf", rm[i][0][0], rm[i][0][1], rm[i][0][2], irm[i][0][0], irm[i][0][1], irm[i][0][2]);
	    		info("  % lf, % lf, % lf | % lf, % lf, % lf", rm[i][1][0], rm[i][1][1], rm[i][1][2], irm[i][1][0], irm[i][1][1], irm[i][1][2]);
	    		info("  % lf, % lf, % lf | % lf, % lf, % lf", rm[i][2][0], rm[i][2][1], rm[i][2][2], irm[i][2][0], irm[i][2][1], irm[i][2][2]);
	    	}

	    	double *x = refmap.get_phys_x(np, quad->get_points(o3));
	    	double *y = refmap.get_phys_y(np, quad->get_points(o3));
	    	double *z = refmap.get_phys_z(np, quad->get_points(o3));
	    	info(" - phys coords [x, y, z]");
	    	for (int i = 0; i < np; i++) {
	    		info("   % lf, % lf, % lf", x[i], y[i], z[i]);
	    	}

	    	// Faces.
	    	for (int iface = 0; iface < e->get_num_faces(); iface++) {
	    		info(" * Face %d", iface);

	    		bool face_const_jac = refmap.is_jacobian_const();
	    	    info("   - const jacobian: %s", face_const_jac ? "yes" : "no");
	    	    if (const_jac) {
	    	    }
	    	    else {
	    	    	Ord2 o2(2, 2);
			    	int fnp = quad->get_face_num_points(iface, o2);

	    	    	double *face_jac = refmap.get_face_jacobian(iface, o2);
	    	    	info("   - jacobian: ");
	    	    	for (int i = 0; i < fnp; i++) info(" % lf", face_jac[i]);

	    	    	double3x3 *face_rm = refmap.get_face_ref_map(iface, o2);
	    	    	double3x3 *face_irm = refmap.get_face_inv_ref_map(iface, o2);
	    	    	info("   - refmap | inv ref map");
	    	    	for (int i = 0; i < fnp; i++) {
	    	    		info("    % lf, % lf, % lf | % lf, % lf, % lf", face_rm[i][0][0], face_rm[i][0][1], face_rm[i][0][2], face_irm[i][0][0], face_irm[i][0][1], face_irm[i][0][2]);
	    	    		info("    % lf, % lf, % lf | % lf, % lf, % lf", face_rm[i][1][0], face_rm[i][1][1], face_rm[i][1][2], face_irm[i][1][0], face_irm[i][1][1], face_irm[i][1][2]);
	    	    		info("    % lf, % lf, % lf | % lf, % lf, % lf", face_rm[i][2][0], face_rm[i][2][1], face_rm[i][2][2], face_irm[i][2][0], face_irm[i][2][1], face_irm[i][2][2]);
	    	    	}

	    	    	Point3D *face_nm = refmap.get_face_normal(iface, o2);
	    	    	info("   - outer normals:");
	    	    	for (int i = 0; i < fnp; i++) {
	    	    		info("     [% lf, % lf, % lf]", face_nm[i].x, face_nm[i].y, face_nm[i].z);
	    	    	}

	    	    	double *fx = refmap.get_face_phys_x(iface, o2);
	    	    	double *fy = refmap.get_face_phys_y(iface, o2);
	    	    	double *fz = refmap.get_face_phys_z(iface, o2);
	    	    	info("   - phys coords [x, y, z]");
	    	    	for (int i = 0; i < fnp; i++) {
	    	    		info("     % lf, % lf, % lf", fx[i], fy[i], fz[i]);
	    	    	}
	    	    }
	    	}

	    	// TODO: edges
	    }

	}


	if (success_test) {
    info("Success!");
    return ERR_SUCCESS;
  }
  else {
    info("Failure!");
    return ERR_FAILURE;
}
