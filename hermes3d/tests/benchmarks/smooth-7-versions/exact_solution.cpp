// Exact solution.
double fn(double x, double y, double z)
{
  switch (ANISO_TYPE) {
    case ANISO_X: return sin(x);
    case ANISO_Y: return sin(y);
    case ANISO_Z: return sin(z);

    case ANISO_X | ANISO_Y: return sin(x) * sin(y);
    case ANISO_X | ANISO_Z: return sin(x) * sin(z);
    case ANISO_Y | ANISO_Z: return sin(y) * sin(z);

    case ANISO_X | ANISO_Y | ANISO_Z: return sin(x) * sin(y) * sin(z);

    default:
    return 0;
  }
}

// Needed for calculation of norms and used by visualizer.
double fndd(double x, double y, double z, 
                      double &dx, double &dy, double &dz)
{
  switch (ANISO_TYPE) {
    case ANISO_X:
      dx = cos(x);
      dy = 0;
      dz = 0;
      break;

    case ANISO_Y:
      dx = 0;
      dy = cos(y);
      dz = 0;
      break;

    case ANISO_Z:
      dx = 0;
      dy = 0;
      dz = cos(z);
      break;

    case ANISO_X | ANISO_Y:
      dx = cos(x) * sin(y);
      dy = sin(x) * cos(y);
      dz = 0;
      break;

    case ANISO_X | ANISO_Z:
      dx = cos(x) * sin(z);
      dy = 0;
      dz = sin(x) * cos(z);
      break;

    case ANISO_Y | ANISO_Z:
      dx = 0;
      dy = cos(y) * sin(z);
      dz = sin(y) * cos(z);
      break;

    case ANISO_X | ANISO_Y | ANISO_Z:
      dx = cos(x) * sin(y) * sin(z);
      dy = sin(x) * cos(y) * sin(z);
      dz = sin(x) * sin(y) * cos(z);
      break;
    }

  return fn(x, y, z);
}
