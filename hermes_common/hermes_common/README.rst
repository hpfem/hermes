The files in this directory are managed in hermes/hermes_common/hermes_common.

The problem is that it's not easy to convince both cmake and cython to compile
it properly. So in the meantime (until someone figures out a fix), we use
symbolic links in hermes1d and 2d.

Here is how to create such symbolic links from scratch:

cd hermes/hermes1d/hermes1d  # change 1d to 2d for hermes2d
rm -rf hermes_common
mkdir hermes_common
cd hermes_common
lndir ../../../hermes_common/hermes_common/ .
