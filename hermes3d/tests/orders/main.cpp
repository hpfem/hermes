#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "config.h"
//#include <getopt.h>
#include <hermes3d.h>

// Order 2.
int test_order_tri() {
	info("test_order_tri.");

	Ord2 a(1), b(2);

	Ord2 c = a + b;
	if (c.order != 3) return ERR_FAILURE;

	Ord2 d = a;
	d += b;
	if (d.order != 3) return ERR_FAILURE;

	Ord2 e = b * 6;
	if (e.order != 12) return ERR_FAILURE;

	Ord2 m = max(b, e);
	if (m.order != e.order) return ERR_FAILURE;

	return ERR_SUCCESS;
}

int test_order_quad() {
	info("test_order_quad.");

	Ord2 a(1, 2), b(3, 4);

	Ord2 c = a + b;
	if (c.x != 4 || c.y != 6) return ERR_FAILURE;

	Ord2 d = a;
	d += b;
	if (d.x != 4 || d.y != 6) return ERR_FAILURE;

	Ord2 e = b * 6;
	if (e.x != 18 || e.y != 24) return ERR_FAILURE;

	Ord2 f = b * Ord2(2, 3);
	if (f.x != 6 || f.y != 12) return ERR_FAILURE;

	Ord2 m = max(b, e);
	if (m.x != e.x || m.y != e.y) return ERR_FAILURE;

	return ERR_SUCCESS;
}

// Order 3.
int test_order_tetra() {
	info("test_order_tetra.");

	Ord3 a(1), b(2);

	Ord3 c = a + b;
	if (c.order != 3) return ERR_FAILURE;

	Ord3 d = a;
	d += b;
	if (d.order != 3) return ERR_FAILURE;

	Ord3 e = b * 6;
	if (e.order != 12) return ERR_FAILURE;

	Ord3 m = max(b, e);
	if (m.order != e.order) return ERR_FAILURE;

	return ERR_SUCCESS;
}

int test_order_hex() {
	info("test_order_hex.");

	Ord3 a(1, 2, 3), b(3, 4, 2);
	Ord3 x;

	Ord3 c = a + b;
	if (c.x != 4 || c.y != 6 || c.z != 5) return ERR_FAILURE;

	Ord3 d = a;
	d += b;
	if (d.x != 4 || d.y != 6 || c.z != 5) return ERR_FAILURE;

	Ord3 e = b * 6;
	if (e.x != 18 || e.y != 24 || e.z != 12) return ERR_FAILURE;

	Ord3 f = b * Ord3(2, 3, 4);
	if (f.x != 6 || f.y != 12 || f.z != 8) return ERR_FAILURE;

	Ord3 m = max(b, e);
	if (m.x != e.x || m.y != e.y || m.z != e.z) return ERR_FAILURE;

	Ord3 z(2, 1, 4);
	x = a + b + z;
	if (x.x != 6 || x.y != 7 || x.z != 9) return ERR_FAILURE;

	if (c == Ord3(4, 6, 5)) ;
	else return ERR_FAILURE;

	if (c != Ord3(4, 6, 5)) return ERR_FAILURE;

	Ord1 edge_ref_order[] = { 3, 4, 3, 4, 2, 2, 2, 2, 3, 4, 3, 4 };
	for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++) {
		if (b.get_edge_order(iedge) != edge_ref_order[iedge]) return ERR_FAILURE;
	}

	Ord2 face_ref_order[] = {
		Ord2(4, 2), Ord2(4, 2), Ord2(3, 2), Ord2(3, 2), Ord2(3, 4), Ord2(3, 4)
	};
	for (int iface = 0; iface < Hex::NUM_FACES; iface++) {
		if (b.get_face_order(iface) != face_ref_order[iface]) return ERR_FAILURE;
	}

	return ERR_SUCCESS;
}

int main() {
	int ret = ERR_SUCCESS;

  ret = test_order_tri();
  if(ret == ERR_SUCCESS)
    ret = test_order_quad();
  if(ret == ERR_SUCCESS)
    ret = test_order_tetra();
  if(ret == ERR_SUCCESS)
    ret = test_order_hex();

  if (ret == ERR_SUCCESS) {
    info("Success!");
    return ERR_SUCCESS;
  }
  else {
    info("Failure!");
    return ERR_FAILURE;
  }
}
