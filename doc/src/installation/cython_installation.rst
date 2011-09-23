===================================
Installation and Updating of Cython
===================================

How to update and install cython in FEMhub
------------------------------------------

Use these instructions if your FEMhub installation failed with an error "Cython too old".
Firstly, go to the directory where you have all your repositories (such as /home/pavel/repos/). 
Then, in your terminal type::

     git clone https://github.com/cython/cython.git
     cd cython
     femhub --shell
     python setup.py install

Note: If FEMhub is not found, then you need to add the path of your local femhub directory 
to your "PATH" variable in the .bashrc file::

    export PATH=\${PATH}:/home/pavel/repos/femhub/


How to update and install cython (outside of FEMhub)
----------------------------------------------------

Go to the directory where you have all your repositories (such as /home/pavel/repos/). 
Then, in your terminal type::

    git clone https://github.com/cython/cython.git
    cd cython
    python setup.py install --home=~/usr

Add the path of cython to your "PATH" and "PYTHONPATH" variable in your .bashrc file::

    export PYTHONPATH=.:\$HOME/usr/lib/python:\$PYTHONPATH
    export PATH=.:\$HOME/usr/bin:\$PATH



























 







