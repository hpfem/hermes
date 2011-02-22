Remote Computing (20-remote-computing)
---------------------

**Git reference:** Tutorial example `20-remote-computing <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/P10-miscellaneous/20-remote-computing>`_. 

This example shows how to save visualization data if you are working 
on a distant computer and cannot use ScalarView, OrderView, or 
related classes directly. There are two basic options:

  * Use the method Solution::save() that saves a complete 
    Solution including Mesh and element orders. Then you can fetch the 
    file and use Solution::load() to restore the Solution
    on your local machine. 
  * Use Linearizer::save_data() that only saves linearized data for direct 
    OpenGL processing. After fetching the file, you can use the methods
    ScalarView::Linearizer::load_data() and ScalarView::show_linearizer_data()
    on your local machine.

The underlying model for computation is the tutorial example 09-timedep. The 
part of the code that is relevant for this example is::

    if (ts % OUTPUT_FREQUENCY == 0) {
      Linearizer lin;
      int item = H2D_FN_VAL_0;
      double eps = HERMES_EPS_NORMAL;
      double max_abs = -1.0;
      MeshFunction* xdisp = NULL; 
      MeshFunction* ydisp = NULL;
      double dmult = 1.0;
      lin.process_solution(&tsln, item, eps, max_abs, xdisp, ydisp, dmult);
      char* filename = new char[100];
      sprintf(filename, "tsln_%d.lin", ts);

      // Save Linearizer data.
      lin.save_data(filename);
      info("Linearizer data saved to file %s.", filename);

      // Save complete Solution.
      sprintf(filename, "tsln_%d.dat", ts);
      bool compress = false;   // Gzip compression not used as it only works on Linux.
      tsln.save(filename, compress);
      info("Complete Solution saved to file %s.", filename);
    }

In the code above, do not worry about the parameters 'xdisp', 'ydisp' and 'dmult'
as these are only used to deform the domain (in linear elasticity problems and such,
see example 08-system).



