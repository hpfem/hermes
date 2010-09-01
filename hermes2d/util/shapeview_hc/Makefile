include ../../Make.config

HERMES = ../../lib/libhermes2d-$(CC)-g-cplx.a

shapeview_hc: shapeview_hc.o $(HERMES)
	gcc -ggdb -o $@ shapeview_hc.o $(HERMES) -lm -lglut

shapeview_hc.o: shapeview_hc.cpp 
	gcc -c -ggdb -DH2D_COMPLEX shapeview_hc.cpp	 -I ../../src

clean:
	rm -f *.o shapeview_hc
