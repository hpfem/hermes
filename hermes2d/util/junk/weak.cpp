


/*

  form    := type [ident "," ident ":"] expr
  type    := "vol" | "surf"
  expr    := term | term "+" expr | term "-" expr
  term    := power | power "*" term | power "/" term
  power   := factor | factor "^" expon
  expon   := number | "(" expr ")"
  factor  := number | ident | ident_partial | spvar | "(" expr ")"
  partial := "x" | "y" | "xx" | "yy" | "xy"
  spvar   := "x" | "y"


  TODO: complex unit
  TODO: functions (norm, abs,
  TODO: boolean & ternary expressions?

  operations: ROTx, ADD, SUB, MUL, DIV, NEG, SQR, CUB?, POW
  composites: LAPx, CNVx
  functions:  NRM, ABS
  reductions: SUM

  order table: char[order1][order2]


*/



IDEAS FOR A NEW DISCRETE PROBLEM INTERFACE



  WeakForm      weakform.h
  LinSystem     system.h
  LinSolver     solver.h
  UmfpackSolver solver_umfpack.h
  PardisoSolver solver_pardiso.h
  EigSolver     eigsolver.h


  WeakForm wf(1);
  wf.add_biform(0, 0, biform_0_0, true);
  wf.add_liform(0, liform_0);

  UmfpackSolver solver;

  LinSystem sys(&weakform, &solver);
  sys.set_spaces(&space);
  sys.assemble();
  sys.solve(&sln);

  PardisoSolver solver;
  PetscKSPSolver solver;

  sys.assemble_rhs();
  sys.solve(&sln);


  WeakForm wf(2);
  wf.add_biform(0, 0, "u_x*v_x + u_y*v_y");
  wf.add_biform(1, 1, "u_x*v_x + u_y*v_y");
  wf.add_liform(0, "omega*v");
  wf.add_liform(1, "omega*v");

  WeakForm wf(3);
  wf.def_const("Re", Re);
  wf.def_const(2, "Re", Re, "tau", tau);
  wf.def_ext_fn(2, "xprev", &xprev, "yprev", &yprev);
  wf.def_exa_fn("F", rhs, order, adapt);
  wf.add_biform(0, 0, "(u_x*v_x + u_y*v_y)/Re + u*v/tau");
  wf.add_biform(0, 0, "(xprev*u_x + yprev*u_y)*v");
  wf.add_biform(1, 1, "(u_x*v_x + u_y*v_y)/Re + u*v/tau");
  wf.add_biform(1, 1, "(xprev*u_x + yprev*u_y)*v");
  wf.add_biform(0, 2, "-u*v_x", antisym);
  wf.add_biform(1, 2, "-u*v_y", antisym);
  wf.add_liform(0, "u*v/tau");
  wf.add_liform(1, "u*v/tau");

  wf.def_const(2, "Re", Re, "tau", tau);
  wf.det_base(3, "u1", "u2", "p");
  wf.det_test(3, "v1", "v2", "q");
  wf.def_ext_fn(2, "X", &xprev, "Y", &yprev);
  wf.set_eqn(0, "(u1,v1)/tau + [u1,v1]/Re + (X*u1_x + Y*u1_y, v1) - (p,v1_x) = (xprev,v1)/tau");
  wf.set_eqn(1, "(u1,v1)/tau + [u2,v2]/Re + (X*u2_x + Y*u2_y, v2) - (p,v2_y) = (yprev,v2)/tau");
  wf.set_eqn(2, "(u1_x,q) + (u2_y,q) = 0");

  wf.compile(); // optional

                    +
                 /     \
              *           *
            /   \       /   \
          u_x   v_x   u_y   v_y


   R1 = MULdd(u_x, v_y)
   R2 = MULdd(u_y, v_y)
   R1 = ADD(R1, R2)

   r3 = ADD(r1, r2, np)
   r3 = MUL(r1, r2, np)
   r2 = ROT(rx, ry, ma, mb, np)
   r3 = MULd(r1x, r1y, m1a, m1b, r2, np)
   r3 = MULdd(r1x, r1y, m1a, m1b, r2x, r2y, m2a, m2b, np)
   r3 += LAP(r1x, r1y, r2x, r2y, m00, m01, m10, m11, np)
   r3 += CNV()


   struct { (*op), result, char arg[8] }

  registers:

    0: *u       12:
    1: *u_x
    2: *u_y
    3: *u_xx
    4: *u_yy
    5: *u_xy

    6: *v
    7: *v_x
    8: *v_y
    9: *v_xx
   10: *v_yy
   11: *v_xy

   load instructions?
   program: { op addr, data of specific size } ...
   registers: double* reg[100];  -- all indirect
   - op: load fn val, dx, dy
   - op: declare tmp reg - allocate in buffer

   program: load info, instructions
   - load info: id = -1 base, -2 test, 0 first extern, 1 second extern, ...
                a,b = item


  wf.set_eqn(0, "(r,q)/tau + (X*r_x + Y*r_y, q) + (r*(X_x + Y_y), q) = (R,q)/tau");
  wf.set_eqn(1, "(u1,v1)/tau + (X*u1_x + Y*u1_y, v1) = (X,v1)/tau + (P/R, v1_x)");
  wf.set_eqn(2, "(u2,v2)/tau + (X*u2_x + Y*u2_y, v2) = (Y,v2)/tau + (P/R, v2_y)");
  wf.set_eqn(3, "(e,f)/tau + (X*e_x + Y*e_y, f) + (e*(X_x + Y_y), f) = "
                            " (E,f)/tau - (X*P_x + Y*P_y, f) - (P*(X_x + Y_y), f)");

  surface: " <u,v>_subdom + ... "



  EigSolver
  SlepcEigSolver eigsolver;

  WeakForm mass(1), stiff(1);
  mass.set_eqn (0, "(u,v)");
  stiff.set_eqn(0, "[u,v] = 0");

  EigSystem eig(&mass, &stif, &eigsolver);
  eig.set_spaces(1, &space);
  eig.assemble();
  eig.solve(&sln);


  LSFEM:

  "(u1*v1) + tau * (v1*(u1prev*u1_x + u2prev*u1_y) + u1*(u1prev*v1_x + u2prev*v1_y))"
  "        + tau^2 * (u1prev*u1_x + u2prev*u1_y) * (u1prev*v1_x + u2prev*v1_y)"

  t1 = (u1prev*u1_x + u2prev*u1_y);
  t2 = (u1prev*v1_x + u2prev*v1_y);
  result = u1*v1 + tau*(u1*t2 + v1*t1) + tau*tau*t1*t2;
