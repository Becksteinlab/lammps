.. index:: fix graphics

fix graphics command
====================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID graphics Nevery keyword args ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* graphics = style name of this fix command
* Nevery = update graphics information every this many time steps
* one or more keyword/args pairs may be appended
* keyword = *sphere* or *units*

  .. parsed-literal::

       *sphere* args = type x y z R
         type = an atom type value to select the color of the sphere
         x, y, z = position of the center of the sphere (distance units)
         R = sphere radius (distance units)
         any of x, y, z, R can be a variable (see below)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all graphics 100 sphere 0.0 0.0 15.0 3.0 sphere 0.0 0.0 5.0 1.0
   fix 1 all graphics 1000 sphere v_x v_y 0.0 v_radius

Description
"""""""""""

This fix allows to add arbitrary objects to images rendered with
:doc:`dump image <dump_image>` using the *fix* keyword.

The *Nevery* keyword determines how often the graphics object data is
updated.  This should be the same value as the corresponding *N*
parameter of the :doc:`dump <dump>` image command.  LAMMPS will stop
with an error message if the settings for this fix and the dump command
are not compatible.

The *type* quantity determines the color of the object.  Its represents
an *atom* type and the object will be colored the same as the
corresponding atom type when the *type* coloring scheme is used in the
:doc:`dump image fix <dump_image>` command is used.  The color may also
be that of the atom type's element or just a globally set constant color
for *all* objects of this fix instance, which can be changed using a
:doc:`dump modify fcolor <dump_image>` command.

Available graphics objects are (see above for exact command line syntax):

- *sphere* - a sphere defined by its center location and its radius
     
Most of the quantities defining a graphics object can be specified as an
equal-style :doc:`variable <variable>`, namely *x*, *y*, *z*, or *R* for
a *sphere*.  If any of these values is a variable, it should be
specified as `v_name`, where `name` is the variable name.  In this case,
the variable will be evaluated each *Nevery* timestep, and its value used
to define the indenter geometry.

Note that equal-style variables can specify formulas with various
mathematical functions, and include :doc:`thermo_style <thermo_style>`
command keywords for the simulation box parameters and timestep and
elapsed time.  Thus it is easy to specify graphics object properties
like position, orientation, radius or more that change as a function of
time or span consecutive runs in a continuous fashion.  For the latter,
see the *start* and *stop* keywords of the :doc:`run <run>` command and
the *elaplong* keyword of :doc:`thermo_style custom <thermo_style>` for
details.

For example, if a sphere's x-position is specified as v_x, then this
variable definition will keep it's center at a relative position in the
simulation box, 1/4 of the way from the left edge to the right edge,
even if the box size changes:

.. code-block:: LAMMPS

   variable x equal "xlo + 0.25*lx"

Similarly, either of these variable definitions will move the sphere
from an initial position at 2.5 at a constant velocity of 5:

.. code-block:: LAMMPS

   variable x equal "2.5 + 5*elaplong*dt"
   variable x equal vdisplace(2.5,5)

If a sphere's radius is specified as v_r, then these variable
definitions will grow the size of the sphere at a specified rate.

.. code-block:: LAMMPS

   variable r0 equal 0.0
   variable rate equal 1.0
   variable r equal "v_r0 + step*dt*v_rate"

Dump image info
"""""""""""""""

.. versionadded:: TBD

Fix graphics is designed to be used with the *fix* keyword of :doc:`dump
image <dump_image>`.  The fix will pass geometry information about the
objects listed on the command line to *dump image* so that they are
included in the rendered image.

.. The *fflag1* setting of *dump image fix* has no impact on rendering a

.. The *fflag2* setting allows you to adjust the radius of the rendered


Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.

None of the :doc:`fix_modify <fix_modify>` options apply to this fix.

Restrictions
""""""""""""

none

Related commands
""""""""""""""""

none

Default
"""""""

none
