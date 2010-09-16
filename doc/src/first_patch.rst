=============================
How to Submit the First Patch
=============================

The following is an embarrassingly trivial git primer
whose objective is to show you how to create and send 
your first patch without losing much time and good humor. 
We begin with cloning the Hermes git repository and 
continue through setting 
up the .gitconfig file, creating a new branch, committing 
changes, and generating patches. A good reference for 
further reading is given at the end. 

Clone the Hermes Git Repository
-------------------------------

To clone the repository, type

::

    git clone http://git.hpfem.org/git/hermes.git

This will create a new directory hermes/ with a copy 
of the entire Hermes git repository. Before doing anything 
else, you may want to build Hermes1D, Hermes2D or Hermes3D 
(whiever one you are interested in) to make sure that 
everything works. So change dir to either hermes1d/, hermes2d/
or hermes3d/ and type::

    cmake .
    make

The list of prerequisites and installation instructions 
for various platforms can be found 
`here <http://hpfem.org/hermes/doc/src/intro-2.html>`_.

Create the .gitconfig File
--------------------------

The .gitconfig file can be used to define your identity
for git as well as useful abbreviations. Change dir to your 
home directory. Then adjust and save the following as 
"~/.gitconfig":

::

    [user]
	    name = Pavel Solin
	    email = solin.pavel@gmail.com

    [core]
	    editor = vim

    [color]
	    ui = true
    [color "branch"]
	    current = yellow reverse
	    local = yellow
	    remote = green
    [color "diff"]
	    meta = yellow bold
	    frag = magenta bold
	    old = red bold
	    new = green bold
	    whitespace = red reverse
    [color "status"]
	    added = yellow
	    changed = green
	    untracked = cyan
    [core]
	    whitespace=fix,-indent-with-non-tab,trailing-space,cr-at-eol

    [alias]
	    st = status
	    ci = commit
	    br = branch
	    co = checkout
	    df = diff
	    lg = log -p

Create a Local Branch
---------------------

Change dir back to hermes1d/, hermes2d/ or hermes3d/. 
You can get an overview of existing branches by typing::

    git branch 

This will show you something like this:

  .. image:: hermes2d/img/intro/terminal-git.png
   :align: center
   :width: 600
   :alt: Terminal screenshot

If this is your first time, then you will see
just the master branch with the star next to it,
which tells you that there are no other branches.

If you want to make any changes to the source files, then 
it is a good idea to always create a new branch for it. 
This is done by typing

::

    git co -b test-1

where test-1 is the name of your new local branch. Now you 
can do any changes you like and you do not have to be afraid
of damaging your master branch. HOWEVER, you always must 
commit your changes as described below. 
Unless you commit your changes, git does not 
know that they belong to your local branch. Then you are not 
able to switch branches, and YOU ARE IN TROUBLE!

Make and Commit Your Changes
----------------------------

The best way to get your first patch in 
is to look into the source files in the 
directory src/, find any function that 
is missing a comment, and fix this. 
Say that this file was src/file.cpp.
After adding the comment there, type:

::

    git add src/file.cpp
    git commit

The latter command will invoke a basic text editor 
where you will be asked to enter a one-line comment
describing your changes. If you decide to skip this 
and commit with an empty line, your commit will not 
be accepted. 

Create and Send a Patch
-----------------------

You are almost there! Just type 

::

    git format-patch -1

and a new text file starting with three zeros will be 
created. This is a "patch". The parameter '-1' in there
means that you want only the last commit included in 
the patch. If you typed '-2', git would include the last 
two commits, etc. 

Last, send an email with the patch to the relevant 
mailing list hermes1d@googlegroups.com, hermes2d@googlegroups
or hermes3d@googlegroups.com. Begin the subject 
line with saying "[PATCH] ...", and attach the 
text file with the patch to your email. Someone
will look at your patch shortly.

Change to Master and Update the Repository
------------------------------------------

Before changing to a different branch, type::

    git st

This stands for 'git status'. You will see 
something like this:

  .. image:: hermes2d/img/intro/terminal-git-2.png
   :align: center
   :width: 600
   :alt: Terminal screenshot

The green font tells you that git has the latest 
version of the file. All modified files in red 
need to be added using "git add". It is a good
idea to go through the untracked files too, in case
that you wish to add some of them as well. 
Related to the sample screenshot above, after 
typing 

::

    git add src/intro-2.rst
    git st

you will see

  .. image:: hermes2d/img/intro/terminal-git-3.png
   :align: center
   :width: 600
   :alt: Terminal screenshot

Now you can proceed with "git commit" as described above. 
After the commit, you can switch to the master branch::

    git co master

This brings you to the point where you can 
return to the beginning of this short
tutorial, and start working on a new change.

To update your master to the latest state of
the repository, just type::

    git remote add origin http://git.hpfem.org/git/hermes.git

This tells git where to download the git repository from
(needs to be done just the first time). Then type::

    git pull origin master

Special Note on Sphinx Docs
---------------------------

The Sphinx documentation you are just reading is also 
part of the Hermes git repository and can be developed
in the same way as source files of Hermes. This very 
file can be found in doc/src/intro-2.rst. After 
making any changes to Sphinx docs, type::

    make html

in the doc/ directory. This will launch 
a build process that will take a few seconds. 
After it finishes, type::

    firefox _build/html

This will open a new tab in your Firefox where you will
see something like 

  .. image:: hermes2d/img/intro/firefox.png
   :align: center
   :width: 600
   :alt: Firefox screenshot

Click on the link "index.html" and you should see
the local form of your Sphinx docs that include your 
changes. 

Further Reading and Video
-------------------------

Git is very powerful and we covered just a tiny part of 
it. After the above works for you, please
read more about git in `Pro Git <http://progit.org/book/>`_.

Also watch this `YouTube video <http://www.youtube.com/watch?v=OFkgSjRnay4>`_
by Scott Chacon.

GitHub
------

You should also learn how to upload
your local branch to `GitHub <http://github.com/>`_
instead of sending a patch, since this makes the
work with your changes easier. 

Good luck and let us know if you think 
that this document could be improved!


 







