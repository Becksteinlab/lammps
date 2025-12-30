.. index:: fix graphics/replica

fix graphics/replica command
============================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID graphics/replica Nevery type keyword args ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* graphics/replica = style name of this fix command
* Nevery = update graphics information every this many time steps
* keyword = *display* or *average*

  .. parsed-literal::

     *display* args = type radius transparency
       type = zero to retain the atom type or an atom type to set the color to that of that atom type
       radius = radius for the atoms
       transparency = transparency setting for the atoms, a value from 0 (invisible) to 1 (fully opaque)
     *average* args = type radius transparency
       type = zero to retain the atom type or an atom type to set the color to that of that atom type
       radius = radius for the atoms or 0 to set the radius to that of the largest distance from the center
       transparency = transparency setting for the atoms, a value from 0 (invisible) to 1 (fully opaque)

Examples
""""""""

.. code-block:: LAMMPS

   fix sf1 water graphics/replica 200 display 0 0.2 1.0 average 3 0.0 0.25

Description
"""""""""""

.. versionadded:: TBD

This fix allows to add spheres to images rendered with :doc:`dump image
<dump_image>` using the *fix* keyword to represent atoms from all
replicas of a multi-replica simulation.

The *group-ID* sets the group ID of the atoms selected to be
represented.  This may be a dynamic group.

The *Nevery* keyword determines how often the replica graphics data is
updated.  This should be the same value as the corresponding *N*
parameter of the :doc:`dump <dump>` image command.  LAMMPS will stop
with an error message if the settings for this fix and the dump command
are not compatible.

There are two keywords available that determine what is shown: *display*
and *average*.  With *display* all atoms from all replica and are in the
fix group will be displayed.  With *average* only the average position
of the atoms with the same atom-ID across all replica will be shown.

The *type* quantity determines the color of the objects.  Its represents
an *atom* type and the atoms will be colored the same as the
corresponding atom type when the *type* coloring scheme is used in the
:doc:`dump image fix <dump_image>` command.  If the value of *type* is
0 then the atom type of the individual atoms is used.

The *radius* quantity determines the radius of the atoms.  The value of
*radius* is 0 then largest distance of an atom to the average position
from all replicas is used.

The *transparency* quantity determines the transparency of the objects.
Its value must be between 0 (invisible) and 1 (fully opaque).

-----------

Dump image info
"""""""""""""""

.. versionadded:: TBD

Fix graphics/replica is designed to be used with the *fix* keyword of
:doc:`dump image <dump_image>`.  The fix will add spheres based on the
atoms in the fix group across all replica to *dump image* so that they
are included in the rendered image.

The *fflag1* setting of *dump image fix* are currently ignored.

and *fflag2* setting of *dump image fix* is used as an adjustment
to the radius of the rendered sphere.  Since the radius is already
determined by this fix, it is recommended to set this flag to 0.0.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.

None of the :doc:`fix_modify <fix_modify>` options apply to this fix.

Restrictions
""""""""""""

This fix is part of the REPLICA package.  It is only only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix graphics <fix_graphics>`, :doc:`fix graphics/arrows <fix_graphics_arrows>`

Default
"""""""

none
