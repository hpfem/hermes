#!/usr/bin/perl -w
#
# generate shapeset for tertahedron
#
# Usage: $? > output.cc
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
		return 1;
#		return $order - 1;
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

# names of shape functions
@shape_fns = ();
@edge_indices = ();
@face_indices = ();
@bubble_indices = ();
@index_to_order = ();

printf "// Tetrahedral shape functions\n";
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
	$fn_name = sprintf("lobatto_f%d", $fn_idx);
#	$fn_name = sprintf("%s_%d", $fnn, 0);
	$fn = sprintf("lambda%d(x, y, z)", $vfn_idx[$v]);
	out_shape_fn($fn_name, $fn);
	printf("\n");

	push(@shape_fns, $fn_name);
	push(@index_to_order, 1);

#	printf($shape_fn_type." %s[] = { %s };\n", $fnn, $fn_name);
#	printf("\n");

	$fn_idx++;
}

# higher orders

#for ($edge = 0; $edge < $num_edges; $edge++) {
#	($f1, $f2) = ($vertex_idx[$edge][0], $vertex_idx[$edge][1]);
#
#	push(@{$edge_indices[$edge]}, ($f1, $f2));		# degree 0 for edges (2 orientations)
#	push(@{$edge_indices[$edge]}, ($f2, $f1));		# degree 1 for edges (2 orientations)
#}

#for ($face = 0; $face < $num_faces; $face++) {
#	push(@{$face_indices[$face]}, (0, 0, 0, 0, 0, 0));		# degree 0 for edges (6 orientations)
#	push(@{$face_indices[$face]}, (0, 0, 0, 0, 0, 0));		# degree 1 for edges (6 orientations)
#	push(@{$face_indices[$face]}, (0, 0, 0, 0, 0, 0));		# degree 2 for edges (6 orientations)
#}


for ($degree = 2; $degree <= $max_degree; $degree++) {
	# comment
	printf "// DEGREE $degree\n";
	printf "//------------\n";
	printf("\n");

	print "// Edge shape functions, degree $degree\n";
	printf("\n");

	for ($edge = 0; $edge < $num_edges; $edge++) {
		print "// edge $edge\n";

		$fn_name = sprintf("lobatto_f%d", $fn_idx);
		$fn = sprintf("lambda%d(x, y, z) * lambda%d(x, y, z) * phi%d(lambda%d(x, y, z) - lambda%d(x, y, z))",
			$vertex_idx[$edge][0], $vertex_idx[$edge][1], $degree - 2, $vertex_idx[$edge][0], $vertex_idx[$edge][1]);

		if ($degree % 2 == 0) {
			out_shape_fn($fn_name, $fn);
			printf("\n");

			push(@{$edge_indices[$edge][0]}, scalar(@shape_fns));
			push(@{$edge_indices[$edge][1]}, scalar(@shape_fns));

			push(@shape_fns, $fn_name);
			push(@index_to_order, $degree);
		}
		else {
			$fn_name_0 = sprintf("%s_%d", $fn_name, 0);
			out_shape_fn($fn_name_0, $fn);
			printf("\n");

			$fn_name_1 = sprintf("%s_%d", $fn_name, 1);
			out_shape_fn($fn_name_1, "-".$fn);
			printf("\n");

			push(@{$edge_indices[$edge][0]}, scalar(@shape_fns));
			push(@{$edge_indices[$edge][1]}, scalar(@shape_fns) + 1);

			push(@shape_fns, $fn_name_0, $fn_name_1);
			push(@index_to_order, $degree, $degree);
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

				$fnn = sprintf("lobatto_f%d", $fn_idx);
#				@fn_names = ();
				for ($fci = 0; $fci < 6; $fci++) {		# 6 possibilities
					$fn_name = sprintf("%s_%d", $fnn, $fci);

					$fn = sprintf("lambda%d(x, y, z) * lambda%d(x, y, z) * lambda%d(x, y, z) * ".
						"phi%d(lambda%d(x, y, z) - lambda%d(x, y, z)) * ".
						"phi%d(lambda%d(x, y, z) - lambda%d(x, y, z))",
						$face_idx[$face][$fci][0], $face_idx[$face][$fci][1], $face_idx[$face][$fci][2],
						$phi[0], $face_idx[$face][$fci][1], $face_idx[$face][$fci][0],
						$phi[1], $face_idx[$face][$fci][0], $face_idx[$face][$fci][2]);
					out_shape_fn($fn_name, $fn);
					printf("\n");

					push(@{$face_indices[$face][$fci]}, scalar(@shape_fns));

					push(@shape_fns, $fn_name);
					push(@index_to_order, $degree);

#					push @fn_names, $fn_name;
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


				$fn_name = sprintf("lobatto_f%d", $fn_idx);
#				@fn_names = ();
#				$fn_name = sprintf("%s_%d", $fnn, 0);

				$fn = sprintf("lambda0(x, y, z) * lambda1(x, y, z) * lambda2(x, y, z) * lambda3(x, y, z)");
				$fn .= " * phi".$phi[0]."(lambda0(x, y, z) - lambda1(x, y, z))";
				$fn .= " * phi".$phi[1]."(lambda2(x, y, z) - lambda1(x, y, z))";
				$fn .= " * phi".$phi[2]."(lambda3(x, y, z) - lambda1(x, y, z))";

				out_shape_fn($fn_name, $fn);
				printf("\n");

				push(@bubble_indices, scalar(@shape_fns));

				push(@shape_fns, $fn_name);
				push(@index_to_order, $degree);

#				printf($shape_fn_type." %s[] = { %s };\n", $fnn, $fn_name);
#				printf("\n");

				$fn_idx++;
			}
		}
	}
}

print "\n";
printf($shape_fn_type." lobatto_tetra_fn[] = {\n");
printf("\t");

$i = 0;
for $fn (@shape_fns) {
	printf("%s, ", $fn);
	$i++;
	if ($i > 7) {
		# simple line wrapping
		printf("\n\t");
		$i = 0;
	}
}
print "};\n";

print "\n// vertices //\n\n";

# vertex indices
print "\n";
printf("static int lobatto_tetra_vertex_indices[] = { 0, 1, 2, 3 };\n");
print "\n";

print "\n// edges //\n\n";

# edge indices
for ($edge = 0; $edge < $num_edges; $edge++) {
	printf("static int lobatto_tetra_edge_indices_%d_0[] = { %s };\n", $edge, join(", ", @{$edge_indices[$edge][0]}));
	printf("static int lobatto_tetra_edge_indices_%d_1[] = { %s };\n", $edge, join(", ", @{$edge_indices[$edge][1]}));
#	$ei[0] = join(", ", @{$edge_indices[$edge][0]});
#	$ei[1] = join(", ", @{$edge_indices[$edge][1]});
	print "\n";
	printf("static int *lobatto_tetra_edge_indices_%d[] = {\n", $edge);
	printf("\tlobatto_tetra_edge_indices_%d_0,\n", $edge);
	printf("\tlobatto_tetra_edge_indices_%d_1\n", $edge);
	printf("};\n");
	print "\n";
}
print "\n";
printf("static int **lobatto_tetra_edge_indices[] = {\n");
for ($face = 0; $face < $num_faces; $face++) {
	printf("\tlobatto_tetra_edge_indices_%d,\n", $face);
}
printf("};\n");
print "\n";


# edge count
@edge_count = ();
for ($order = 0; $order < $max_degree; $order++) {
	$edge_count[$order] = get_num_edge_fns($order);
}
printf("static int lobatto_tetra_edge_count[] = { %s };\n", join(", ", @edge_count));
print "\n";

print "\n// faces //\n\n";

# face indices
for ($face = 0; $face < $num_faces; $face++) {
	for ($ori = 0; $ori < 6; $ori++) {
#		$fi[$ori] = join(", ", @{$face_indices[$face][$ori]});
		printf("static int lobatto_tetra_face_indices_%d_%d[] = {\n", $face, $ori);
		printf("\t%s\n", join(", ", @{$face_indices[$face][$ori]}));
		printf("};\n");
		printf("\n");

		$fi[$ori] = sprintf("lobatto_tetra_face_indices_%d_%d", $face, $ori);
	}

	printf("static int *lobatto_tetra_face_indices_%d[] = {\n", $face);
#	for ($i = 0; $i <= $max_degree; $i++) {
#		printf("\t%s,\n", join(", ", splice(@{$face_indices[$face]}, $i * 6, 6)));
#	}
	printf("\t%s\n", join(",\n\t", @fi));
	printf("};\n");
	printf("\n");
}

print "\n";
printf("static int **lobatto_tetra_face_indices[] = {\n");
for ($face = 0; $face < $num_faces; $face++) {
	printf("\tlobatto_tetra_face_indices_%d,\n", $face);
}
printf("};\n");
print "\n";


# face count
@face_count = ();
for ($order = 0; $order <= $max_degree; $order++) {
	$face_count[$order] = get_num_face_fns($order);
}
printf("static int lobatto_tetra_face_count[] = { %s };\n", join(", ", @face_count));
print "\n";

print "\n// bubbles //\n\n";

# bubble indices
printf("static int lobatto_tetra_bubble_indices_all_orders[] = {\n");
printf("\t%s\n", join(", ", @bubble_indices));
printf("};\n");
print "\n";

printf("static int *lobatto_tetra_bubble_indices[] = {\n");
printf("\tNULL, NULL, NULL, NULL,\n");
for ($order = 4; $order <= $max_degree; $order++) {
	printf("\tlobatto_tetra_bubble_indices_all_orders,\n");
}
printf("};\n");
print "\n";


# bubble count
@bubble_count = ();
for ($order = 0; $order <= $max_degree; $order++) {
	$bubble_count[$order] = get_num_bubble_fns($order);
}
printf("static int lobatto_tetra_bubble_count[] = { %s };\n", join(", ", @bubble_count));
print "\n";


# index to odrer
print "\n";
printf("static int lobatto_tetra_index_to_order[] = {\n");
printf("\t");
$i = 0;
for $fn (@index_to_order) {
	printf("%d, ", $fn);
	$i++;
	if ($i > 30) {
		# simple line wrapping
		printf("\n\t");
		$i = 0;
	}
}
print "};\n";


