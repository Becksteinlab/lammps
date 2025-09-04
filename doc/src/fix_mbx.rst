.. index:: fix mbx

fix mbx command
===============

Syntax
""""""

.. code-block:: LAMMPS

    fix ID group-ID mbx num_mon_types monomer_specification keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* mbx = style name of this fix command
* num_mon_types = number of monomer types
* monomer_specification = for each monomer type, the following
  arguments must be specified:

  .. parsed-literal::

    * monomer_name = name of the monomer type
    * monomer_lower_atom_index = lower atom index of the monomer
      (e.g. 1 for O in water)
    * monomer_upper_atom_index = upper atom index of the monomer
      (e.g. 2 for H in water)
    * monomer_num_atoms = number of atoms in the monomer
    * atom_ids = list of atom IDs in the monomer, in the
      order they appear in the MBX configuration file

* one or more keyword/value pairs may be appended

  .. parsed-literal::
    keyword = *json* or *print/dipoles* or *print/settings*
        *json* arg = name of MBX JSON configuration file
        *print/dipoles* = print dipole moments as part of fix variable output
        *print/settings* = print MBX settings to logfile


Examples
""""""""

.. code-block:: LAMMPS

    # For a system involving water (atom types O=1, H=2)
    processors      * * * map xyz
    pair_style      mbx 9.0
    pair_coeff      * *
    compute         mbx all pair mbx
    fix             mbx_fix all mbx 1 h2o 1 2 3 1 2 2 json mbx.json


    # For a system involving ch4 (atom types C=1, H=2) and
    # water (atom types O=3, H=4)
    processors      * * * map xyz
    pair_style      mbx 9.0
    pair_coeff      * *
    compute         mbx all pair mbx
    fix             mbx_fix all mbx 2 ch4 1 2 5 1 2 2 2 2 h2o 3 4 3 3 4 4 json mbx.json

    # For a system involving water (atom types 0=12, H=13) in a hybrid simulation
    processors * * * map xyz
    pair_style      hybrid/overlay mbx 9.0 lj/cut 9.0 coul/exclude 9.0
    pair_coeff      * * mbx
    pair_coeff      1*11 1*11 coul/exclude
    compute         mbx all pair mbx
    fix             mbx_fix all mbx 2 dp1 1 11 1 1 h2o 12 13 3 12 13 13 json mbx.json

See ``examples/PACKAGES/mbx`` for additional examples of how to use MBX in LAMMPS.


Description
"""""""""""

The MBX (Many-Body eXpansion) software is a C++ library that provides
access to many-body energy (MB-nrg) potential energy functions, such as
the MB-pol water model. Developed over the past decade, these potential
energy functions integrate physics-based and machine-learned many-body
terms trained on electronic structure data calculated at the "gold
standard" coupled-cluster level of theory. :ref:`(Gupta) <Gupta2>`


This fix instructs LAMMPS to call the `MBX library <_mbxwebsite>`_
in order to simulate MB-nrg models such as MB-pol. This fix must be
used in conjunction with the :doc:`pair mbx <pair_mbx>` command.

The MBX library code development is available at
`https://github.com/paesanilab/MBX <https://github.com/paesanilab/MBX>`_.
A detailed discussion of the code can be found in the manuscript :ref:`(Riera) <Riera2>`.


See ``examples/PACKAGES/mbx`` for complete examples of how to use
this fix command.

For hybrid simulations involving MB-nrg and non-MB-nrg molecules in the
same simulation, one can use :doc:`pair_style hybrid/overlay
<pair_hybrid>` to combine the MB-nrg molecules with other pair styles,
such as :doc:`lj/cut <pair_lj>`.  Do note that all electrostatics must be
computed within MBX, so the :doc:`coul/exclude <pair_coul>` pair_style
should usually be applied on the non-MB-nrg molecules.  See
``examples/PACKAGES/mbx`` for a complete hybrid example.


The *num_mon_types* argument specifies the number of different MB-nrg
monomer types in the system.

The *monomer_specification* argument provides the details for each
monomer type.  This information is used by MBX to map the LAMMPS atom
IDs to the corresponding MBX monomer types.  For each monomer type, the
following information must be provided:

* monomer_name = name of the monomer type
  (e.g. h2o for water, ch4 for methane)
* monomer_lower_atom_index = lower atom index of the monomer
  (e.g. 1 for O in water)
* monomer_upper_atom_index = upper atom index of the monomer
  (e.g. 2 for H in water)
* monomer_num_atoms = number of atoms in the monomer
  (e.g. 3 for water, 5 for methane)
* atom_ids = list of LAMMPS atom IDs in the monomer, in the
  order they appear in the MB-nrg potential.
  (e.g. 1 2 2 for water, as it corresponds to O H H)


The *json* argument specifies the name of the MBX JSON configuration
file to use, such as `mbx.json`.  If this file is not provided, the fix
will attempt to use a default configuration.

The *print/dipoles* argument enables the printing of dipole moments as
part of the fix variable output.  This is useful for performing
vibrational spectroscopy calculations such as IR, Raman, and
Sum-Frequency Generation (SFG).

The *print/settings* argument will print the MBX settings to the LAMMPS
logfile at the start of the simulation.  This is used for debugging and
ensuring that the correct settings are being applied.

If you have questions not answered by this documentation, please reach
out to us at `https://groups.google.com/g/mbx-users
<https://groups.google.com/g/mbx-users>`_

Restrictions
""""""""""""

This fix is part of the MBX package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
page for more info. This fix also relies on the presence of :doc:`pair
mbx <pair_mbx>` command.

There can only be one fix mbx command active at a time.

Due to the usage of Partridge and Schwenke charges for MB-pol, all
electrostatic interactions are calculated internally in MBX.  Therefore
one should never calculate coulombic interactions in LAMMPS such as
using `coul/cut` or `coul/long` when also using MBX.



Related commands
""""""""""""""""

:doc:`pair mbx <pair_mbx>`,
:doc:`pair hybrid/overlay <pair_hybrid>`,
:doc:`pair coul/exclude <pair_coul>`

-----------

.. _Riera2:

**(Riera)** M. Riera, C. Knight, E. Bull-Vulpe, X. Zhu, H. Agnew, D. Smith, A. Simmonett, F. Paesani, J. Chem. Phys. 159, 054802 (2023)

.. _Gupta2:

**(Gupta)** S. Gupta, E. Bull-Vulpe, H. Agnew, S. Iyer, X. Zhu, R. Zhou, C. Knight, F. Paesani, J. Chem. Theory Comput. 21, 1938 (2025)

.. _mbxwebsite: https://mbxsimulations.com
