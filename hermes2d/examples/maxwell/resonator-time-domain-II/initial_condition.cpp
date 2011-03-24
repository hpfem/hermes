class InitialConditionWave : public ExactSolutionVector
{
public:
  InitialConditionWave(Mesh* mesh) : ExactSolutionVector(mesh) {};

  virtual scalar2 value (double x, double y) const {
    return scalar2(-cos(x) * sin(y), sin(x) * cos(y));
  }

  virtual void derivatives (double x, double y, scalar2& dx, scalar2& dy) const {
    dx[0] = sin(x) * sin(y);
    dx[1] = cos(x) * cos(y);
    dy[0] = -cos(x) * cos(y);
    dy[1] = -sin(x) * sin(y);
  }

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(20);
  }
};
