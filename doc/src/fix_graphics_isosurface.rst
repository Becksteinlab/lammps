.. index:: fix graphics/isosurface

fix graphics/isosurface command
===============================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID graphics/isosurface Nevery isovalue radius keyword args ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* graphics/isosurface = style name of this fix command
* Nevery = update graphics information every this many time steps
* isovalue = isovalue for the particle property isosurface selection
* radius = radius describing the spread of the atoms to the density grid (distance units)
* one or more keyword/args pairs may be appended
* keyword = *quality* or *property* or *filename* or *binary* or *pad*

  .. parsed-literal::

     *quality* keyword = isosurface grid resolution setting
        keyword = one of *min*, *low*, *med*, *high*, or *max*
     *property* value = per-atom property used to create the isosurface grid
        value = *none*, *mass*,  c_ID, c_ID[i], f_ID, f_ID[i], v_name
           *none* = 1.0 for all atoms
           *mass* = mass of the atoms
           c_ID = per-atom vector calculated by a compute with ID
           c_ID[I] = Ith column of per-atom array calculated by a compute with ID
           f_ID = per-atom vector calculated by a fix with ID
           f_ID[I] = Ith column of per-atom array calculated by a fix with ID
           v_name = per-atom vector calculated by an atom-style variable with name
     *filename* name = name pattern for output of a sequence of STL format mesh files (must contain a \* character to be replaced by the timestep number)
     *binary* logical = select whether to output a binary STL file (default is text mode)
     *pad* number = pad the timestep in the output file name with zeroes to have this many digits (default is 0)

Examples
""""""""

.. code-block:: LAMMPS

   fix sf1 water graphics/isosurface 200 0.1 2.5 quality high property mass
   fix stl water graphics/isosurface 200 0.01 1.5 filename water-isosurface-*.stl pad 5

Description
"""""""""""

.. versionadded:: TBD

This fix allows to add an isosurface graphics object representing the
triangulated isosurface at a given isovalue on a grid to images rendered
with :doc:`dump image <dump_image>` using the *fix* keyword and
optionally to output the computed mesh as a series of STL format files
for external processing.

The *group-ID* sets the group ID of the atoms selected to be represented
by the isosurface.  This may be a dynamic group.

The *Nevery* keyword determines how often the isosurface graphics data
is updated.  This should be the same value as the corresponding *N*
parameter of the :doc:`dump <dump>` image command.  LAMMPS will stop
with an error message if the settings for this fix and the dump command
are not compatible.

The isosurface objects will be colored by the atom type that is closest
to each isosurface grid cell when the *type* coloring scheme is used in
the :doc:`dump image fix <dump_image>` command.  The color is that of
the atom type's element color instead with the *element* coloring
scheme, or just a globally set constant color for the whole isosurface
with the *const* coloring scheme.  That color can be set with the
*fcolor* keyword of the :doc:`dump modify <dump_image>` command.

The isosurface's transparency setting is fully opaque by default and can
be changed with the *ftrans* keyword of the :doc:`dump modify
<dump_image>` command.

The *isovalue* argument sets the isovalue used to compute the
isosurface.  The optimum value depends on the property on that is being
used and the information that is supposed to be conveyed.  It usually
requires some experimentation in combination with varying the *radius*
setting.

The *radius* argument sets the width of the gaussian distribution
function used to distribute the per-particle data across the grid.  Its
value controls the smoothness of the isosurface and - as mentioned
above - may need some experimentation in combination with the choice of
*isovalue* to achieve the desired output.

The *quality* keyword can have any of these words as argument: "min",
"low", "med", "high", or "max", and selects the grid resolution used
for the isosurface.  The actual grid dimensions depend on the geometry
of the simulation cell.

The optional *property* keyword controls what property is used to set
the values at the grid points for the isosurface.  The default setting
of *none* just uses a value of 1.0, resulting in the data grid
representing a smoothed out number density.  Other possible arguments
are *mass* (for representing the smoothed out mass density) or a
references to a a :doc:`compute <compute>`, a :doc:`fix <fix>`, or a
reference to an atom-style :doc:`variable <variable>`.  The compute or
fix must produce a per-atom vector or array, not a global or local
quantity.  In case the property is a per-atom array, the column must be
selected.

-----------

Dump image info
"""""""""""""""

.. versionadded:: TBD

Fix graphics/isosurface is designed to be used with the *fix* keyword of
:doc:`dump image <dump_image>`.  The fix will construct an isosurface
based on the atom positions, the selected property.  of the atoms in the
fix group and pass the graphics geometry information about it to *dump
image* so that it is included in the rendered image.

The *fflag1* setting of *dump image fix* determines whether the
isosurface will be rendered as a set of connected triangles (1) or as a
mesh of cylinders (2).

If using a mesh of cylinders, the *fflag2* setting determines the
diameter of the cylinders.

The *quality* settings of "min" and "low" work best with the cylinder
mesh setting while the other quality settings are more suitable for a
triangle mesh.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.

None of the :doc:`fix_modify <fix_modify>` options apply to this fix.

Restrictions
""""""""""""

None.

Related commands
""""""""""""""""

:doc:`fix_graphics <fix_graphics>`, :doc:`fix_graphics/arrows <fix_graphics_arrows>`

Defaults
""""""""

quality = medium, property = none, binary = no, pad = 0, filename = none
