.. index:: pair_style mbx

pair_style mbx command
======================

Syntax
""""""

.. code-block:: LAMMPS

    pair_style mbx cutoff

* cutoff = real-space cutoff for MBX in Angstroms


Examples
""""""""

.. code-block:: LAMMPS

    # For a system involving water (atom types O=1, H=2)
    processors      * * * map xyz
    pair_style      mbx 9.0
    pair_coeff      * * 1 h2o 1 2 2 json mbx.json
    compute         mbx all pair mbx


    # For a system involving ch4 (atom types C=1, H=2) and
    # water (atom types O=3, H=4)
    processors      * * * map xyz
    pair_style      mbx 9.0
    pair_coeff      * * 2 ch4 1 2 2 2 2 h2o 3 4 4 json mbx.json
    compute         mbx all pair mbx

    # For a system involving water (atom types O=12, H=13) in a hybrid simulation
    processors      * * * map xyz
    pair_style      hybrid/overlay mbx 9.0 lj/cut 9.0 coul/exclude 9.0
    pair_coeff      * * mbx 2 dp1 1*11 h2o 12 13 13 json mbx.json
    pair_coeff      1*11 1*11 coul/exclude
    compute         mbx all pair mbx

See ``examples/PACKAGES/mbx`` for full examples of how to use MBX in LAMMPS.


Description
"""""""""""

The MBX (Many-Body eXpansion) software is a C++ library that provides
access to many-body energy (MB-nrg) potential energy functions, such as
the MB-pol water model.  Developed over the past decade, these potential
energy functions integrate physics-based and machine-learned many-body
terms trained on electronic structure data calculated at the "gold
standard" coupled-cluster level of theory. :ref:`(Gupta) <Gupta3>`

This pair_style instructs LAMMPS to call the
`MBX library <https://mbxsimulations.com>`_ in order to simulate
MB-nrg models such as MB-pol. The MBX library code development is available at
`https://github.com/paesanilab/MBX <https://github.com/paesanilab/MBX>`_.
MBX is heavily OpenMP parallelized (OMP), and the OMP_NUM_THREADS
environment variable should be properly set to the number of threads desired.
A detailed discussion of the code structure can be found in the
manuscript :ref:`(Riera) <Riera>`, while a detailed description of the
performance scaling can be found in the manuscript :ref:`(Gupta) <Gupta3>`.

The *cutoff* argument specifies the real-space cutoff for MBX in
Angstroms. For periodic systems, a safe value is 9.0 Angstroms, which
is the cutoff used in the original MB-pol model.  For non-periodic
systems, the cutoff can be set to a large value, such as 100.0
Angstroms, to ensure that all interactions are captured. A larger cutoff
is always safer, but will result in a slower simulation.

For hybrid simulations involving MB-nrg and non-MB-nrg molecules in the
same simulation, one can use :doc:`pair_style hybrid/overlay
<pair_hybrid>` to combine the MB-nrg molecules with other pair styles,
such as :doc:`lj/cut <pair_lj>`. This has been used to simulate
MB-pol water within host frameworks such as metal-organic
frameworks (MOFs) and carbon nanotubes (CNTs).
Do note that all electrostatics must be
computed within MBX, so the :doc:`coul/exclude <pair_coul>` pair_style
should usually be applied on the non-MB-nrg molecules.  See
``examples/PACKAGES/mbx`` for a complete hybrid example.


If you have questions not answered by this documentation, please
reference the MBX website
`mbxsimulations.com <https://mbxsimulations.com>`_ or reach out to us at
`https://groups.google.com/g/mbx-users <https://groups.google.com/g/mbx-users>`_


Pair coeff syntax
"""""""""""""""""

MBX is many-body method, and only a single pair_coeff command is needed
to specify the mapping of LAMMPS atom IDs to MBX monomers. The syntax is as follows:

.. code-block:: LAMMPS

    pair_coeff * * num_mon_types mon_name atom_mapping <mon_name2> <atom_mapping2> ... json mbx.json print/settings

* num_mon_types = number of monomer types in the system
* mon_name = name of the monomer type (e.g. h2o, ch4, etc)
* atom mapping = list of LAMMPS atom types that correspond to the atoms in the monomer
* *json* arg = specifies the name of the MBX json configuration file, such as mbx.json
* print/settings = optionally print MBX settings to logfile


The *num_mon_types* argument specifies the number of different MB-nrg
monomer types in the system.

For each monomer type, the *mon_name* argument specifies the name of
the monomer, such as `h2o` for water or `ch4` for methane. The *atom
mapping* argument specifies then the mapping of LAMMPS atom types to
the atoms in the monomer, such as `1 2 2` for water (O=1, H=2).
For hybrid simulations, the `dp1` (drude particle) monomer
should be used to represent the non-MB-nrg molecules. `dp1` is a
special monomer in MBX in that its *atom_mapping* can be a range of
LAMMPS atom types, such as `1*11` to represent atom types 1 through 11.

.. warning::

    Since the MB-nrg models (e.g. MB-pol) used in MBX include both permanent
    electrostatics and polarization, it is important to ensure that no
    electrostatic interactions are calculated twice.

    When using MBX, **all electrostatics are handled internally by MBX.**
    Therefore, one should never use a coulombic pair style in LAMMPS
    such as `coul/cut` or `coul/long` when also using MBX.

    When performing a hybrid simulation using dp1,
    note that many frameworks (Amber, CHARMM, OPLS, ClayFF etc.) require the usage of
    special_bonds to exclude some bonded coulomb interactions (1-2, 1-3, and/or 1-4).
    Since MBX is handling all
    electrostatics, this should therefore be accounted for using the :doc:`coul/exclude <pair_coul>` command.

    .. code-block:: LAMMPS

        pair_style      hybrid/overlay mbx 9.0 lj/cut 9.0 coul/exclude 9.0
        pair_coeff      * * mbx 2 dp1 1*11 h2o 12 13 13 json mbx.json
        pair_coeff      1*11 1*11 coul/exclude

The *json* argument specifies the name of the MBX JSON configuration
file to use, such as `mbx.json`.  If this file is not provided, the fix
will attempt to use a default configuration.

The *print/settings* argument will print the MBX settings to the LAMMPS
logfile at the start of the simulation.  This is used for debugging and
ensuring that the correct settings are being applied.


Restrictions
""""""""""""

This pair_style is part of the MBX package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Due to the usage of Partridge and Schwenke charges for MB-pol,
all electrostatic interactions are calculated internally in MBX.
Therefore one should never calculate coulombic interactions in
LAMMPS such as using `coul/cut` or `coul/long` when also using MBX.

MBX is primarily tested to work with `units real`. If you encounter
issues with other unit styles, please contact the developers.

Related commands
""""""""""""""""

:doc:`pair hybrid/overlay <pair_hybrid>`,
:doc:`pair coul/exclude <pair_coul>`

-----------

.. _Riera:

**(Riera)** M. Riera, C. Knight, E. Bull-Vulpe, X. Zhu, H. Agnew, D. Smith, A. Simmonett, F. Paesani, J. Chem. Phys. 159, 054802 (2023)

.. _Gupta3:

**(Gupta)** S. Gupta, E. Bull-Vulpe, H. Agnew, S. Iyer, X. Zhu, R. Zhou, C. Knight, F. Paesani, J. Chem. Theory Comput. 21, 1938 (2025)

