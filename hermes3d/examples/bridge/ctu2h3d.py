#!/usr/bin/env python

if __name__ == "__main__":
    point_start = True
    point_write = True

    h3d_f = open("test.mesh3d", "w")
    pd0_f = open("test.pd0.mesh3d", "w")
    hd0_f = open("test.hd0.mesh3d", "w")

    pd0_id = []

    for line in open("most-sup.top"):

        loop_counter = []
        line = line.strip()
        cols = line.split(' ')

        if(len(cols) == 1):
            item_count = int(cols[0])

            if point_start:
                h3d_f.write("# points\n")
                pd0_f.write("# points with displacement 0\n")
                point_start = False
            else:
                h3d_f.write("# hex\n")
                hd0_f.write("# hex with atleast one point with displacement 0\n")
                point_write = False

            h3d_f.write(str(item_count) + "\n")
        else:
            if point_write:
                h3d_f.write(cols[1] + " " + cols[2] + " " + cols[3] + "\n")

                if (int(cols[-1]) >= 10 and int(cols[-1]) <= 16):
                    pd0_f.write(cols[0] + " " + cols[1] + " " + cols[2] + " " + cols[3] + "\n")
                    pd0_id.append(int(cols[0]))

            else:
                h3d_f.write(cols[2] + " " + cols[3] + " " + cols[4] + " " + cols[5] + " " + cols[6] + " " + cols[7] + " " + cols[8] + " " + cols[9] + "\n")

                str = cols[0] + " " + cols[1] + " "
                hex_with_disp0 = False;
                for i in range(2,10):
                    if int(cols[i]) in pd0_id:
                        str += "|" + cols[i] + "|" + " "
                        hex_with_disp0 = True;
                    else:
                        str += cols[i] + " "

                if hex_with_disp0:
                    hd0_f.write(str + "\n")

    h3d_f.write("\n# prims\n0\n\n# tri\n0\n\n")

    h3d_f.close()
    pd0_f.close()
    hd0_f.close()
