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
       radius = radius for arrows
     *force* args = scale radius
       scale = scale factor for the force vector to determine the arrow length
       radius = radius for arrows
     *velocity* args = scale radius
       scale = scale factor for the velocity vector to determine the arrow length
       radius = radius for arrows
     *variable* args = xval yval zval radius
       xval = xvalue for arrow vector
       yval = xvalue for arrow vector
       zval = xvalue for arrow vector
       radius = radius for arrows
     *chunk* args = chunk-ID posval vecval scale radius transparency
       chunk-ID = ID of :doc:`compute chunk/atom <compute_chunk_atom>` command
       posval = ID of a per-chunk compute or fix that computes the position vector for the chunk
       vecval = ID of a per-chunk compute or fix that computes the arrow vector for the chunk
       yval = xvalue for arrow vector
       zval = xvalue for arrow vector
       radius = radius for arrows

Examples
""""""""

.. code-block:: LAMMPS

   fix vec all graphics/arrows 10 velocity 20.0 0.066

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
will show the per-atom dipole vector, *force* the per-atom force, *velocity*
the per-atom velocity, *variable* a custom vector constructed from three
constants or atom-style variables. With the *chunk* keyword per-chunk vectors
will be shown.

The *scale* quantity determines the length of the arrows.  It should be
chosen so that when multiplied with the per-atom vector quantity the result
is of the same order of magnitude as atom positions, so that the vectors
can be seen well.

The *radius* quantity determines the width of the arrows.

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

Restrictions
""""""""""""

None

Related commands
""""""""""""""""

:doc:`fix_graphics <fix_graphics>`

Default
"""""""

none
