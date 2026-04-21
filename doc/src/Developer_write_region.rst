Writing new region styles
^^^^^^^^^^^^^^^^^^^^^^^^^

Region styles define geometric shapes that can be referenced by many
LAMMPS commands (e.g. :doc:`group <group>`, :doc:`fix wall/region
<fix_wall_region>`, :doc:`fix deposit <fix_deposit>`,
:doc:`delete_atoms <delete_atoms>`).  All region styles are derived from
the ``Region`` base class defined in ``src/region.h``.  An overview of
the available methods is given on the :doc:`Modify_region` page.

Unlike most other styles, region styles do not have a dedicated
"EXTRA-REGION" package.  New region styles are placed directly in the
``src/`` directory.  If you feel that your contribution should be added
to a package, please consult with the :doc:`LAMMPS developers
<Intro_authors>` first.  The contributed code needs to support the
:doc:`traditional GNU make build process <Build_make>` **and** the
:doc:`CMake build process <Build_cmake>`.

When adding a new region style to ``src/``:

- Add the ``.cpp`` and ``.h`` filenames (with a leading ``/``) to
  ``src/.gitignore`` so that the copies produced by the traditional
  make build system are not accidentally committed to the repository.
- Add the bare filename (without the leading ``/``) to
  ``src/Purge.list`` so that ``make purge`` removes them correctly.

----

Case study: RegSphere
^^^^^^^^^^^^^^^^^^^^^^

As a concrete example we implement a spherical region.  The complete
implementation can be found in ``src/region_sphere.cpp`` and
``src/region_sphere.h``.  The region accepts four arguments (center
coordinates ``x``, ``y``, ``z`` and ``radius``), each of which can be
either a constant number or an equal-style variable (prefixed with
``v_``).

Header file
"""""""""""

After the standard copyright block, the registration block for a region
style uses ``#ifdef REGION_CLASS`` and ``RegionStyle(name,class)``:

.. code-block:: c++

   #ifdef REGION_CLASS
   // clang-format off
   RegionStyle(sphere,RegSphere);
   // clang-format on
   #else

The block between ``#ifdef REGION_CLASS`` and ``#else`` is included by
the ``Domain`` class in ``domain.cpp`` to build the factory map that
connects the style name ``sphere`` with the class ``RegSphere``.  The
``// clang-format`` comments protect the macro from being reformatted.

The class definition:

.. code-block:: c++

   #ifndef LMP_REGION_SPHERE_H
   #define LMP_REGION_SPHERE_H

   #include "region.h"

   namespace LAMMPS_NS {
   class RegSphere : public Region {
     friend class Region2VMD;
     friend class DumpImage;

    public:
     RegSphere(class LAMMPS *, int, char **);
     ~RegSphere() override;
     void init() override;
     int inside(double, double, double) override;
     int surface_interior(double *, double) override;
     int surface_exterior(double *, double) override;
     void shape_update() override;
     void bbox_update() override;
     void set_velocity_shape() override;
     void velocity_contact_shape(double *, double *) override;

    protected:
     double xc, yc, zc;
     double radius;
     int xstyle, xvar;
     int ystyle, yvar;
     int zstyle, zvar;
     int rstyle, rvar;
     char *xstr, *ystr, *zstr, *rstr;

     void variable_check();
   };

   }    // namespace LAMMPS_NS
   #endif
   #endif

The ``friend`` declarations for ``Region2VMD`` and ``DumpImage`` allow
those classes to access the protected members directly for rendering
purposes.

The three pure virtual methods that every region style must override
are ``inside()``, ``surface_interior()``, and ``surface_exterior()``.
``shape_update()`` and ``bbox_update()`` are optional but required when
the region parameters can be controlled by variables.

Implementation file
"""""""""""""""""""

Constructor
'''''''''''

The constructor parses the region-specific arguments after calling the
base class constructor.  The base class ``Region`` constructor processes
the first two arguments (the region ID and style name) as well as any
trailing keyword/value options (handled by ``options()``).  The
region-specific positional arguments start at index 2 (for ``region
sphere``, these are ``x``, ``y``, ``z``, and ``radius`` at positions
2–5):

.. code-block:: c++

   RegSphere::RegSphere(LAMMPS *lmp, int narg, char **arg) :
     Region(lmp, narg, arg), xstr(nullptr), ystr(nullptr), zstr(nullptr), rstr(nullptr)
   {
     options(narg - 6, &arg[6]);

     if (utils::strmatch(arg[2], "^v_")) {
       xstr = utils::strdup(arg[2] + 2);
       xc = 0.0;
       xstyle = VARIABLE;
       varshape = 1;
     } else {
       xc = xscale * utils::numeric(FLERR, arg[2], false, lmp);
       xstyle = CONSTANT;
     }

     if (utils::strmatch(arg[3], "^v_")) {
       ystr = utils::strdup(arg[3] + 2);
       yc = 0.0;
       ystyle = VARIABLE;
       varshape = 1;
     } else {
       yc = yscale * utils::numeric(FLERR, arg[3], false, lmp);
       ystyle = CONSTANT;
     }

     if (utils::strmatch(arg[4], "^v_")) {
       zstr = utils::strdup(arg[4] + 2);
       zc = 0.0;
       zstyle = VARIABLE;
       varshape = 1;
     } else {
       zc = zscale * utils::numeric(FLERR, arg[4], false, lmp);
       zstyle = CONSTANT;
     }

     if (utils::strmatch(arg[5], "^v_")) {
       rstr = utils::strdup(arg[5] + 2);
       radius = 0.0;
       rstyle = VARIABLE;
       varshape = 1;
     } else {
       radius = xscale * utils::numeric(FLERR, arg[5], false, lmp);
       rstyle = CONSTANT;
     }

     if (varshape) {
       variable_check();
       RegSphere::shape_update();
     }

     if (radius < 0.0) error->all(FLERR, "Illegal region sphere radius: {}", radius);

     if (interior) {
       bboxflag = 1;
       if (dynamic || varshape) {
         RegSphere::bbox_update();
       } else {
         extent_xlo = xc - radius;
         extent_xhi = xc + radius;
         extent_ylo = yc - radius;
         extent_yhi = yc + radius;
         extent_zlo = zc - radius;
         extent_zhi = zc + radius;
       }
     }

     cmax = 1;
     contact = new Contact[cmax];
     tmax = 1;
   }

Key points:

- The constructor initializer list sets all ``char*`` members to
  ``nullptr`` to ensure safe deletion in the destructor even if an
  error is thrown before those members are assigned.
- The base class method ``options()`` must be called with the
  arguments *after* the style-specific positional arguments.  Here
  there are four positional arguments (``x``, ``y``, ``z``,
  ``radius``), so ``options(narg - 6, &arg[6])`` passes any remaining
  keyword/value options.
- Arguments that begin with ``v_`` are treated as variable references.
  ``utils::strmatch(arg[i], "^v_")`` returns true if the string starts
  with ``v_``.  The variable name is stored (without the ``v_`` prefix)
  by duplicating ``arg[i] + 2``.  Setting ``varshape = 1`` tells the
  base class that the region shape can change during a run.
- The ``xscale``, ``yscale``, and ``zscale`` factors (set by the base
  class from the ``units box`` or ``units lattice`` option) are applied
  to convert from input units to internal units.
- If all parameters are constants, the bounding box (``extent_xlo``
  etc.) is set directly.  If any are variable-controlled, ``bbox_update()``
  is called to compute it from the current values.
- ``cmax`` and the ``contact`` array (one contact point per surface
  primitive) and ``tmax`` (maximum number of contact points in one
  region query) must be set here.

The destructor frees the variable name strings and the contact array:

.. code-block:: c++

   RegSphere::~RegSphere()
   {
     if (copymode) return;

     delete[] xstr;
     delete[] ystr;
     delete[] zstr;
     delete[] rstr;
     delete[] contact;
   }

Note the ``if (copymode) return;`` guard, which prevents double-freeing
when the region object is used in Kokkos mode.

The init() method (optional)
'''''''''''''''''''''''''''''

.. code-block:: c++

   void RegSphere::init()
   {
     Region::init();
     if (varshape) variable_check();
   }

The ``init()`` method first calls ``Region::init()`` to perform base
class initialization and then re-checks that the variables used to
control the shape still exist.

The inside() method (required)
'''''''''''''''''''''''''''''''

.. code-block:: c++

   int RegSphere::inside(double x, double y, double z)
   {
     double delx = x - xc;
     double dely = y - yc;
     double delz = z - zc;
     double r = sqrt(delx * delx + dely * dely + delz * delz);

     if (r <= radius) return 1;
     return 0;
   }

This method returns 1 if the point ``(x, y, z)`` is inside or on the
surface of the region and 0 otherwise.  For an exterior region (when
the ``side out`` option is used), the base class inverts the result.

The surface_interior() method (required)
'''''''''''''''''''''''''''''''''''''''''

.. code-block:: c++

   int RegSphere::surface_interior(double *x, double cutoff)
   {
     double delx = x[0] - xc;
     double dely = x[1] - yc;
     double delz = x[2] - zc;
     double r = sqrt(delx * delx + dely * dely + delz * delz);
     if (r > radius || r == 0.0) return 0;

     double delta = radius - r;
     if (delta < cutoff) {
       contact[0].r = delta;
       contact[0].delx = delx * (1.0 - radius / r);
       contact[0].dely = dely * (1.0 - radius / r);
       contact[0].delz = delz * (1.0 - radius / r);
       contact[0].radius = -radius;
       contact[0].iwall = 0;
       contact[0].varflag = 1;
       return 1;
     }
     return 0;
   }

This method determines whether a point is within ``cutoff`` distance of
the *inner* surface of the region (used for fix wall/region with the
``side in`` option).  It returns the number of contact points found
(0 or ``cmax`` at most).  Each contact is recorded in the ``contact``
array with:

- ``contact[0].r``: distance from the point to the surface.
- ``contact[0].delx/y/z``: displacement vector from the nearest
  surface point to the input point ``x``.
- ``contact[0].radius``: radius of curvature of the surface at the
  contact point; negative means the surface curves away from ``x``
  (concave interior).
- ``contact[0].iwall``: index of the surface primitive (0 for a
  sphere with a single surface).
- ``contact[0].varflag``: 1 if the region boundary can move (e.g.,
  variable-controlled parameters or a moving region).

The surface_exterior() method (required)
'''''''''''''''''''''''''''''''''''''''''

.. code-block:: c++

   int RegSphere::surface_exterior(double *x, double cutoff)
   {
     double delx = x[0] - xc;
     double dely = x[1] - yc;
     double delz = x[2] - zc;
     double r = sqrt(delx * delx + dely * dely + delz * delz);
     if (r < radius) return 0;

     double delta = r - radius;
     if (delta < cutoff) {
       contact[0].r = delta;
       contact[0].delx = delx * (1.0 - radius / r);
       contact[0].dely = dely * (1.0 - radius / r);
       contact[0].delz = delz * (1.0 - radius / r);
       contact[0].radius = radius;
       contact[0].iwall = 0;
       contact[0].varflag = 1;
       return 1;
     }
     return 0;
   }

This is the counterpart to ``surface_interior()`` for atoms approaching
from outside.  The radius of curvature is now positive (convex
exterior).

The shape_update() method (optional)
'''''''''''''''''''''''''''''''''''''

The ``shape_update()`` method re-evaluates any variable-controlled
parameters and must be implemented if ``varshape = 1``:

.. code-block:: c++

   void RegSphere::shape_update()
   {
     if (xstyle == VARIABLE) xc = xscale * input->variable->compute_equal(xvar);
     if (ystyle == VARIABLE) yc = yscale * input->variable->compute_equal(yvar);
     if (zstyle == VARIABLE) zc = zscale * input->variable->compute_equal(zvar);

     if (rstyle == VARIABLE) {
       radius = xscale * input->variable->compute_equal(rvar);
       if (radius < 0.0) error->one(FLERR, "Variable evaluation in region gave bad value");
     }
   }

``input->variable->compute_equal(xvar)`` evaluates the equal-style
variable whose index ``xvar`` was obtained by ``variable_check()``.

The bbox_update() method (optional)
'''''''''''''''''''''''''''''''''''''

The ``bbox_update()`` method updates the bounding box when the shape
can change (variable-controlled or moving region):

.. code-block:: c++

   void RegSphere::bbox_update()
   {
     if (varshape || dynamic) {
       extent_xlo = xc - radius;
       extent_xhi = xc + radius;
       extent_ylo = yc - radius;
       extent_yhi = yc + radius;
       extent_zlo = zc - radius;
       extent_zhi = zc + radius;
       if (moveflag) {
         extent_xlo += dx;
         extent_xhi += dx;
         extent_ylo += dy;
         extent_yhi += dy;
         extent_zlo += dz;
         extent_zhi += dz;
       }
     }
   }

The ``moveflag``, ``dx``, ``dy``, ``dz`` variables are set by the base
class when the region has a ``move`` keyword and account for the
time-dependent displacement of a moving region.

The variable_check() helper method
''''''''''''''''''''''''''''''''''''

.. code-block:: c++

   void RegSphere::variable_check()
   {
     if (xstyle == VARIABLE) {
       xvar = input->variable->find(xstr);
       if (xvar < 0) error->all(FLERR, "Variable {} for region sphere does not exist", xstr);
       if (!input->variable->equalstyle(xvar))
         error->all(FLERR, "Variable {} for region sphere is invalid style", xstr);
     }
     // ... (analogous for y, z, radius)
   }

This helper looks up each variable name and verifies that it is an
equal-style variable.  It is called both in the constructor and in
``init()`` so that re-use of the same region across multiple runs
works correctly.

----

Summary of region style requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following table summarizes the methods a new region style must or
may implement.

+---------------------------+----------+------------------------------------------------------+
| Method                    | Required | Purpose                                              |
+===========================+==========+======================================================+
| ``inside()``              | yes      | Determine whether a point is inside the region       |
+---------------------------+----------+------------------------------------------------------+
| ``surface_interior()``    | yes      | Find surface contacts from inside                    |
+---------------------------+----------+------------------------------------------------------+
| ``surface_exterior()``    | yes      | Find surface contacts from outside                   |
+---------------------------+----------+------------------------------------------------------+
| ``shape_update()``        | if vars  | Re-evaluate variable-controlled parameters           |
+---------------------------+----------+------------------------------------------------------+
| ``bbox_update()``         | if vars  | Update the bounding box when shape changes           |
+---------------------------+----------+------------------------------------------------------+
| ``init()``                | no       | Initialization; re-check variable existence          |
+---------------------------+----------+------------------------------------------------------+
| ``set_velocity_shape()``  | no       | Set data needed for wall velocity due to shape change|
+---------------------------+----------+------------------------------------------------------+
| ``velocity_contact_shape()`` | no    | Compute wall velocity contribution from shape change |
+---------------------------+----------+------------------------------------------------------+

The ``set_velocity_shape()`` and ``velocity_contact_shape()`` methods
are only needed for regions used with fix wall/gran/region when the
region boundary can move or change shape.
