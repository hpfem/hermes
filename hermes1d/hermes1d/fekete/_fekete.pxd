from numpy cimport ndarray

cpdef double int_f2(ndarray[double, mode="c"] w,
        ndarray[double, mode="c"] values)

cpdef double int_f2_f2(ndarray[double, mode="c"] w,
        ndarray[double, mode="c"] values1,
        ndarray[double, mode="c"] values2,
        )

cpdef get_gauss_points_phys(double a, double b, int n)
