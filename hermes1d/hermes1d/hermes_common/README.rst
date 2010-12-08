This directory should eventually go into hermes/hermes_common/, but the problem
is that it's not easy to convince both cmake and cython to compile it properly.

So in the meantime, we just have it here directly, but don't put any dimension
dependent code here (e.g. in particular, no "hermes1d" words should be in
here), so that later on we can refactor this out of H1D.
