.. index:: fix graphics/arrows

fix graphics/arrows command
===========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID graphics/arrows Nevery mode keyword args ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* graphics/arrows = style name of this fix command
* Nevery = update graphics information every this many time steps
* mode = one of the following modes *dipole* or *force* or *velocity* or *variable* or *chunk*

  .. parsed-literal::

     *dipole* args = scale radius
       scale = scale factor for the dipole moment to determine the arrow length
       radius = radius for arrows (length units)
     *force* args = scale radius
       scale = scale factor for the force vector to determine the arrow length
       radius = radius for arrows (length units)
     *velocity* args = scale radius
       scale = scale factor for the velocity vector to determine the arrow length
       radius = radius for arrows (length units)
     *variable* args = xval yval zval radius
       xval = x value for arrow vector (may be a variable)
       yval = y value for arrow vector (may be a variable)
       zval = z value for arrow vector (may be a variable)
       radius = radius for arrows (length units)
     *chunk* args = chunk-ID posval vecval scale radius transparency
       chunk-ID = ID of :doc:`compute chunk/atom <compute_chunk_atom>` command
       posval = ID of a per-chunk compute or fix that computes the position vector for the chunk
       vecval = ID of a per-chunk compute or fix that computes the arrow vector for the chunk
       radius = radius for arrows (length units)

* zero or more keyword/value pairs may be appended
* keyword = *autoscale*

  .. parsed-literal::

       *autoscale* value = automatically scale arrows so they have an average length of "value"

..       *units* value = *lattice* or *box*
..         lattice = the geometry is defined in lattice units
..         box = the geometry is defined in simulation box units

Examples
""""""""

.. code-block:: LAMMPS

   fix vec all graphics/arrows 10 velocity 20.0 0.066 autoscale 1.5
   fix vec all graphics/arrows 100 variable v_xnorm v_znorm 0.0 0.066

Description
"""""""""""

.. versionadded:: TBD

This fix allows to add arrows to images rendered with :doc:`dump image
<dump_image>` using the *fix* keyword to represent vector properties
with arrows atoms from all atoms in the fix group.

The *group-ID* sets the group ID of the atoms selected to have the selected
property represented.  This may be a dynamic group.

The *Nevery* keyword determines how often the arrows graphics data is
updated.  This should be the same value as the corresponding *N*
parameter of the :doc:`dump <dump>` image command.  LAMMPS will stop
with an error message if the settings for this fix and the dump command
are not compatible.

There are five keywords available that determine what is shown: *dipole*
will show the per-atom dipole vector, *force* the per-atom force,
*velocity* the per-atom velocity, *variable* a custom vector constructed
from three constants or atom- or equal-style variables. With the *chunk*
keyword per-chunk vectors will be shown.

The *xval*\ , *yval*\ , and *zval*\ , arguments to the *variable* mode
define a custom vector that can be composed of numbers or :doc:`atom- or
equal-style variables <variable>`.  If any of these values is a
variable, it should be specified as *v_name*\ , where "name" is the
variable name.  In this case, the variable will be evaluated each
timestep, and its value used to define the arrow for each atom.  Since
variables can reference :doc:`computes <compute>`, :doc:`fixes <fix>`,
:doc:`custom per-atom properties <fix_property_atom>`, and other
variables, this can be used to construct arrows for almost any per-atom
property available in LAMMPS.

The *scale* quantity determines the length of the arrows.  It should be
chosen so that when multiplied with the per-atom vector quantity the result
is of the same order of magnitude as atom positions, so that the vectors
can be seen well.

The *radius* quantity determines the width of the arrows.

The optional *autoscale* keyword allows to dynamically determine the
*scale* quantity so that the average length of the arrows is set to the
value of the keyword's argument.  The computed scale factor can be
accessed by various :doc:`output commands <Howto_output>` as a global
scalar (see below).

-----------

Dump image info
"""""""""""""""

.. versionadded:: TBD

Fix graphics/arrows is designed to be used with the *fix* keyword of
:doc:`dump image <dump_image>`.  The fix will add spheres based on the
atoms in the fix group across all arrows to *dump image* so that they
are included in the rendered image.

The color of the arrows is by default that of the atoms when using color
styles "type" or "element".  With color style "const" the default value
of "white" can be changed using :doc:`dump_modify fcolor <dump_image>`.
Similarly, the transparency follows the atom type or can be changed for
color style "const" with *dump\_modify ftrans*\ .


The *fflag1* and *fflag2* settings of *dump image fix* are currently ignored.

Restart, fix_modify, output, run start/stop, minimize info
==========================================================

No information about this fix is written to :doc:`binary restart files
<restart>`.

None of the :doc:`fix_modify <fix_modify>` options apply to this fix.

This fix computes a global scalar representing the current scale factor
for displaying the arrows, which can be accessed by various
:doc:`output commands <Howto_output>`.  This is the *autoscale*
keyword argument value divided by the average length of the selected
vector property.  If the *autoscale* keyword is not used, it is the
scale value set by the *fix graphics/arrows* command or 1.0.
The scalar value calculated by this fix is "intensive".

Restrictions
""""""""""""

The *dipole* mode requires the use of :doc:`atom style dipole
<atom_style>` or a hybrid atom style that includes it.

Related commands
""""""""""""""""

:doc:`fix_graphics <fix_graphics>`

Default
"""""""

autoscale is off by default
