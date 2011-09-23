===================
Editing Sphinx Docs
===================

The Sphinx documentation that you are reading now is also 
part of the Hermes git repository and can be developed
in the same way as source files of Hermes. This very 
file can be found in doc/src/editing_sphinx.rst. 

Generating and viewing HTML files
---------------------------------

After making your changes, type::

    make html

in the hermes/doc directory. This will launch 
a build process that will take a few seconds. 
After it finishes, type::

    chromium-browser _build/html

or::

    firefox _build/html

(assuming that you use Chromium or Firefox). Other browsers 
can be used as well.

This will open a new tab or window where you will
see something like 

  .. figure:: hermes2d/img/intro/firefox.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: Firefox screenshot

Click on the link "index.html" and you should see
the local form of your Sphinx docs that include your 
changes. 

Generating Latex, PDF and PostScript
------------------------------------

Typing
::

    make latex
 
in hermes/doc will generate the Latex sources of the 
tutorial. The files will be located in _build/latex.
In order to generate a PDF file, issue::

    cd _build/latex
    make all-pdf

A PostScript file can be generated with
::

    make all-ps









