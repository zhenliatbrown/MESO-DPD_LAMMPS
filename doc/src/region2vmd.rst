.. index:: region2vmd

region2cmd command
==================

Syntax
""""""

.. code-block:: LAMMPS

   region2vmd file args

* file = name of VMD script file to write
* args = one or more region IDs may be appended

Examples
""""""""

.. code-block:: LAMMPS

   region2vmd regions.vmd box c1 c2

Description
"""""""""""

.. versionadded:: TBD

Write a `VMD <https:://ks.uiuc.edu/Research/vmd/>`_ Tcl script file with
commands that aim to create a visualization of :doc:`LAMMPS regions
<region>`.  There may be multiple region visualizations stored in a
single file.  Only a limited amount of region styles and settings are
currently supported. See **Restrictions** below.

The visualization is implemented by creating a new (and empty) "VMD
molecule" and then using VMD graphics primitives to represent the region
in VMD.  Each region will be stored in a separate "VMD molecule" with
the name "LAMMPS region <region ID>".

The created file can be loaded into VMD either from the command line
with the -e flag, or from the command prompt with play <script file>, or
from the File menu via "Load VMD visualization state".

----------

Restrictions
""""""""""""

This command is part of the EXTRA-COMMAND package.  It is only enabled
if LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Only the following region styles are currently supported: *block*,
*cylinder*, *cone*, *sphere*.  For region style *cone* one of the two
radii must be zero, since the equivalent VMD graphics primitive does not
support truncated cones.

Moving or rotating regions as well as unions or intersecting regions are
also currently not supported.

Related commands
""""""""""""""""

:doc:`region <region>`

Default
"""""""

none
