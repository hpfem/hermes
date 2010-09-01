#!/usr/bin/perl -w
#
# generate shapeset for tertahedron (derivatives)
#
# Usage: $? { dx | dy | dz } > output.cc
#

#
$max_degree = 10;
$num_vertices = 4;
$num_edges = 6;
$num_faces = 4;

# some globals
$shape_fn_type = "shape_fn_t";

@vfn_idx = ( 1, 2, 0, 3 );


@vertex_idx = (
	[ 1, 2 ],
	[ 2, 0 ],
	[ 1, 0 ],
	[ 1, 3 ],
	[ 2, 3 ],
	[ 0, 3 ]
);

@face_idx = (
	[
		[ 1, 2, 3 ],
		[ 2, 3, 1 ],
		[ 3, 1, 2 ],
		[ 1, 3, 2 ],
		[ 2, 1, 3 ],
		[ 3, 2, 1 ],
	],
	[
		[ 2, 0, 3 ],
		[ 0, 3, 2 ],
		[ 3, 2, 0 ],
		[ 2, 3, 0 ],
		[ 0, 2, 3 ],
		[ 3, 0, 2 ],
	],
	[
		[ 1, 0, 3 ],
		[ 0, 3, 1 ],
		[ 3, 1, 0 ],
		[ 1, 3, 0 ],
		[ 0, 1, 3 ],
		[ 3, 0, 1 ],
	],
	[
		[ 1, 2, 0 ],
		[ 2, 0, 1 ],
		[ 0, 1, 2 ],
		[ 1, 0, 2 ],
		[ 2, 1, 0 ],
		[ 0, 2, 1 ],
	]
);

sub out_shape_fn {
	my ($fn_name, $fn) = @_;

	printf "double %s(double x, double y, double z) {\n", $fn_name;
	printf "\treturn ".$fn.";\n";
	printf "}\n";
}

sub get_num_edge_fns {
	my ($order) = @_;
	if ($order >= 2) {
		return $order - 1;
	}
	else {
		return 0;
	}
}

sub get_num_face_fns {
	my ($order) = @_;
	if ($order >= 3) {
		return ($order - 2);
	}
	else {
		return 0;
	}
}

sub get_num_bubble_fns {
	my ($order) = @_;

	if ($order >= 4) {
		return (($order - 3) * ($order - 2)) / 2;
	}
	else {
		return 0;
	}
}

# -- main --

@shape_fns = ();

if (scalar(@ARGV) < 1) {
	print("ERROR: Not enough params.\n");
	exit -1;
}

$dd = $ARGV[0];
if ($dd ne "dx" && $dd ne "dy" && $dd ne "dz") {
	print("ERROR: Specify dx, dy or dz as command line argument.\n");
	exit -1;
}

printf "// Derivatives of the shape functions ($dd)\n";
printf("\n");

# vertex fns
printf "// DEGREE 1\n";
printf "//------------\n";
printf "\n";
printf "// Vertex shape functions, degree 1\n";
printf "\n";

$fn_idx = 0;
$degree = 0;
for ($v = 0; $v < $num_vertices; $v++) {
	$fn_name = sprintf("lobatto_%s_f%d", $dd, $fn_idx);
	$fn = sprintf("lambda%d$dd(x, y, z)", $vfn_idx[$v]);
	out_shape_fn($fn_name, $fn);
	printf("\n");

	push(@shape_fns, $fn_name);

#	printf($shape_fn_type." %s[] = { %s };\n", $fnn, $fn_name);
#	printf("\n");

	$fn_idx++;
}

# higher orders

for ($degree = 2; $degree <= $max_degree; $degree++) {
	# comment
	printf "// DEGREE $degree\n";
	printf "//------------\n";
	printf("\n");

#	$edge_fn_cnt = get_num_edge_fns($degree);

#	print $edge_fn_cnt."\n";

	print "// Edge shape functions, degree $degree\n";
	printf("\n");

#	print join(", ", $vertex_idx[0][1]);

#	exit;

	for ($e = 0; $e < $num_edges; $e++) {
		print "// edge $e\n";

		$fn_name = sprintf("lobatto_%s_f%d", $dd, $fn_idx);

#		dx = 
#			lambda1x(x, y, z) * lambda2(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) + 
#			lambda1(x, y, z) * lambda2x(x, y, z) * phi0(lambda1(x, y, z) - lambda2(x, y, z)) + 
#			lambda1(x, y, z) * lambda2(x, y, z) * phi0x(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1x(x, y, z) - lambda2x(x, y, z));

		$l1 = sprintf("lambda%d", $vertex_idx[$e][0]);
		$l2 = sprintf("lambda%d", $vertex_idx[$e][1]);
		$phi = sprintf("phi%d", $degree - 2);
		$fn = "\n".
			"\t\t$l1$dd(x, y, z) * $l2(x, y, z) * $phi($l1(x, y, z) - $l2(x, y, z)) +\n".
			"\t\t$l1(x, y, z) * $l2$dd(x, y, z) * $phi($l1(x, y, z) - $l2(x, y, z)) +\n".
			"\t\t$l1(x, y, z) * $l2(x, y, z) * ${phi}dx($l1(x, y, z) - $l2(x, y, z)) * ($l1$dd(x, y, z) - $l2$dd(x, y, z))";


		if ($degree % 2 == 0) {
			out_shape_fn($fn_name, $fn);
			printf("\n");

			push(@shape_fns, $fn_name);
		}
		else {
			$fn_name_0 = sprintf("%s_%d", $fn_name, 0);
			out_shape_fn($fn_name_0, $fn);
			printf("\n");

#			dx = -(
#				lambda1x(x, y, z) * lambda2(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) + 
#				lambda1(x, y, z) * lambda2x(x, y, z) * phi1(lambda1(x, y, z) - lambda2(x, y, z)) + 
#				lambda1(x, y, z) * lambda2(x, y, z) * phi1x(lambda1(x, y, z) - lambda2(x, y, z)) * (lambda1x(x, y, z) - lambda2x(x, y, z))
#			);

#			$fn_1 = "-(".
#				"$l1$dd(x, y, z) * $l2(x, y, z) * $phi($l1(x, y, y) - $l2(x, y, z)) +\n".
#				"\t\t$l1(x, y, z) * $l2$dd(x, y, z) * $phi($l1(x, y, y) - $l2(x, y, z)) +\n".
#				"\t\t$l1(x, y, z) * $l2(x, y, z) * ${phi}dx($l1(x, y, y) - $l2(x, y, z)) * ($l1$dd(x, y, z) - $l2$dd(x, y, z))".
#				")";

			$fn_name_1 = sprintf("%s_%d", $fn_name, 1);
			out_shape_fn($fn_name_1, "-(".$fn.")");
			printf("\n");

			push(@shape_fns, $fn_name_0, $fn_name_1);
		}

		$fn_idx++;
	}

	printf("\n");

	# face
	if ($degree >= 3) {
		print "// Face shape functions, degree $degree\n";
		printf("\n");

		$face_cnt = get_num_face_fns($degree)."\n";

		# for all faces
		for ($face = 0; $face < $num_faces; $face++) {
			print "// face $face\n";
	
			# number of face functions on a face
			for ($j = 0; $j < $face_cnt; $j++) {
				$phi[0] = $j;
				$phi[1] = $degree - 3 - $j;

				$fnn = sprintf("lobatto_%s_f%d", $dd, $fn_idx);
#				@fn_names = ();
				for ($fci = 0; $fci < 6; $fci++) {		# 6 possibilities
					$fn_name = sprintf("%s_%d", $fnn, $fci);

#					dx = 
#						lambda1x(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) + 
#						lambda1(x, y, z) * lambda2x(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) + 
#						lambda1(x, y, z) * lambda2(x, y, z) * lambda3x(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) + 
#						lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0x(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2x(x, y, z) - lambda1x(x, y, z)) * phi1(lambda1(x, y, z) - lambda3(x, y, z)) + 
#						lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi1x(lambda1(x, y, z) - lambda3(x, y, z)) * (lambda1x(x, y, z) - lambda3x(x, y, z));

					$l1 = sprintf("lambda%d", $face_idx[$face][$fci][0]);
					$l2 = sprintf("lambda%d", $face_idx[$face][$fci][1]);
					$l3 = sprintf("lambda%d", $face_idx[$face][$fci][2]);
					$phi0 = sprintf("phi%d", $phi[0]);
					$phi1 = sprintf("phi%d", $phi[1]);

					$fn = "\n".
						"\t\t$l1$dd(x, y, z) * $l2(x, y, z) * $l3(x, y, z) * $phi0($l2(x, y, z) - $l1(x, y, z)) * $phi1($l1(x, y, z) - $l3(x, y, z)) +\n".
						"\t\t$l1(x, y, z) * $l2$dd(x, y, z) * $l3(x, y, z) * $phi0($l2(x, y, z) - $l1(x, y, z)) * $phi1($l1(x, y, z) - $l3(x, y, z)) +\n".
						"\t\t$l1(x, y, z) * $l2(x, y, z) * $l3$dd(x, y, z) * $phi0($l2(x, y, z) - $l1(x, y, z)) * $phi1($l1(x, y, z) - $l3(x, y, z)) +\n".
						"\t\t$l1(x, y, z) * $l2(x, y, z) * $l3(x, y, z) * ${phi0}dx($l2(x, y, z) - $l1(x, y, z)) * ($l2$dd(x, y, z) - $l1$dd(x, y, z)) * $phi1($l1(x, y, z) - $l3(x, y, z)) +\n".
						"\t\t$l1(x, y, z) * $l2(x, y, z) * $l3(x, y, z) * $phi0($l2(x, y, z) - $l1(x, y, z)) * ${phi1}dx($l1(x, y, z) - $l3(x, y, z)) * ($l1$dd(x, y, z) - $l3$dd(x, y, z))\n";
					out_shape_fn($fn_name, $fn);
					printf("\n");

					push(@shape_fns, $fn_name);
				}

#				printf($shape_fn_type." %s[] = { %s };\n", $fnn, join(", ", @fn_names));
#				printf("\n");
				$fn_idx++;
			}
			printf("\n");
		}
	}

	# bubble
	if ($degree >= 4) {
#		$bubble_cnt = get_num_bubble_fns($degree);

		print "// Bubble shape functions, degree $degree\n";
		printf("\n");

		$m = $degree - 4;

		@phi = ();
		for ($i = 0; $i <= $m; $i++) {
			$phi[0] = $i;

			for ($j = 0; $j <= $m - $i; $j++) {
				$phi[1] = $j;
				$phi[2] = $m - $i - $j;


				$fn_name = sprintf("lobatto_%s_f%d", $dd, $fn_idx);

#	dx =
#		lambda0x(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z))
#		lambda0(x, y, z) * lambda1x(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z)) +
#		lambda0(x, y, z) * lambda1(x, y, z) * lambda2x(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z))
#		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3x(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z))
#		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0x(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0x(x, y, z) - lambda1x(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z))
#		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0x(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2x(x, y, z) - lambda1x(x, y, z)) * phi2(lambda3(x, y, z) - lambda1(x, y, z))
#		lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * phi0(lambda0(x, y, z) - lambda1(x, y, z)) * phi0(lambda2(x, y, z) - lambda1(x, y, z)) * phi2x(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3x(x, y, z) - lambda1x(x, y, z));

				$phi0 = sprintf("phi%d", $phi[0]);
				$phi1 = sprintf("phi%d", $phi[1]);
				$phi2 = sprintf("phi%d", $phi[2]);

				$fn = "\n". 
					"\t\tlambda0$dd(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * $phi0(lambda0(x, y, z) - lambda1(x, y, z)) * $phi1(lambda2(x, y, z) - lambda1(x, y, z)) * $phi2(lambda3(x, y, z) - lambda1(x, y, z)) +\n".
					"\t\tlambda0(x, y, z) * lambda1$dd(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * $phi0(lambda0(x, y, z) - lambda1(x, y, z)) * $phi1(lambda2(x, y, z) - lambda1(x, y, z)) * $phi2(lambda3(x, y, z) - lambda1(x, y, z)) +\n".
					"\t\tlambda0(x, y, z) * lambda1(x, y, z) * lambda2$dd(x, y, z) * lambda3(x, y, z) * $phi0(lambda0(x, y, z) - lambda1(x, y, z)) * $phi1(lambda2(x, y, z) - lambda1(x, y, z)) * $phi2(lambda3(x, y, z) - lambda1(x, y, z)) +\n".
					"\t\tlambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3$dd(x, y, z) * $phi0(lambda0(x, y, z) - lambda1(x, y, z)) * $phi1(lambda2(x, y, z) - lambda1(x, y, z)) * $phi2(lambda3(x, y, z) - lambda1(x, y, z)) +\n".
					"\t\tlambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * ${phi0}dx(lambda0(x, y, z) - lambda1(x, y, z)) * (lambda0$dd(x, y, z) - lambda1$dd(x, y, z)) * $phi1(lambda2(x, y, z) - lambda1(x, y, z)) * $phi2(lambda3(x, y, z) - lambda1(x, y, z)) +\n".
					"\t\tlambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * $phi0(lambda0(x, y, z) - lambda1(x, y, z)) * ${phi1}dx(lambda2(x, y, z) - lambda1(x, y, z)) * (lambda2$dd(x, y, z) - lambda1$dd(x, y, z)) * $phi2(lambda3(x, y, z) - lambda1(x, y, z)) +\n".
					"\t\tlambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z) * $phi0(lambda0(x, y, z) - lambda1(x, y, z)) * $phi1(lambda2(x, y, z) - lambda1(x, y, z)) * ${phi2}dx(lambda3(x, y, z) - lambda1(x, y, z)) * (lambda3$dd(x, y, z) - lambda1$dd(x, y, z))";

				out_shape_fn($fn_name, $fn);
				printf("\n");

				push(@shape_fns, $fn_name);

#				printf($shape_fn_type." %s[] = { %s };\n", $fnn, $fn_name);
#				printf("\n");

				$fn_idx++;
			}
		}
	}
}

print "\n";
printf($shape_fn_type." lobatto_tetra_%s[] = {\n", $dd);
printf("\t");

$i = 0;
for $fn (@shape_fns) {
	printf("%s, ", $fn);
	$i++;
	if ($i > 6) {
		# simple line wrapping
		printf("\n\t");
		$i = 0;
	}
}
print "};\n";

#$fnn = 0;
#while ($fnn < 286) {
#	@names = ();
#	for ($j = 0; $j < 8; $j++) {
#		push @names, sprintf("lobatto_%s_f%d", $dd, $fnn);
#		$fnn++;
#		last if ($fnn >= 286);
#	}
#	print "\t".join(", ", @names).", \n";
#}
#print "};\n";



