#! /usr/bin/env python

from utils import conv_graph
import hydrogen_romanowski as romanowski
import hydrogen_pask_R as pask
import hydrogen_certik_rR as hermes_romanowski
import hydrogen_certik_R
import hydrogen_certik_hp_lowest_R

conv_graph((
    (pask, "pask"),
    (hermes_romanowski, "roman. hermes"),
    (romanowski, "roman. paper"),
    (hydrogen_certik_hp_lowest_R, "hp lowest"),
    ), 0, "R", eigs=3)
#conv_graph2(R_x, R_y, R2_x, R2_y, 0, "rR")
#conv_graph2(R_x, R_y, R2_x, R2_y, 1, "rR")
#conv_graph2(R_x, R_y, R2_x, R2_y, 2, "rR")

#conv_graph2(R3_x, R3_y, R_x, R_y, 0, "R")
#conv_graph2(R3_x, R3_y, R_x, R_y, 1, "R")
#conv_graph2(R3_x, R3_y, R_x, R_y, 2, "R")
