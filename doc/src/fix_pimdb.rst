.. index:: fix pimdb/langevin
.. index:: fix pimdb/nvt

fix pimdb/langevin command
=========================

fix pimdb/nvt command
====================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID style keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* style = *pimdb/langevin* or *pimdb/nvt* = style name of this fix command
* zero or more keyword/value pairs may be appended
* keywords for style *pimdb/nvt*

  .. parsed-literal::
       *keywords* = *method* or *fmass* or *sp* or *temp* or *nhc*
       *method* value = *pimd* or *nmpimd*
       *fmass* value = scaling factor on mass
       *sp* value = scaling factor on Planck constant
       *temp* value = temperature (temperature units)
       *nhc* value = Nc = number of chains in Nose-Hoover thermostat

* keywords for style *pimdb/langevin*

  .. parsed-literal::
       *keywords* = *method* or *integrator* or *ensemble* or *fmmode* or *fmass* or *scale* or *temp* or *thermostat* or *tau* or *iso* or *aniso* or *barostat* or *taup* or *fixcom* or *lj*
       *method* value = *pimd*
       *integrator* value = *obabo* or *baoab*
       *fmmode* value = *physical* or *normal*
       *fmass* value = scaling factor on mass
       *temp* value = temperature (temperature unit)
          temperature = target temperature of the thermostat
       *thermostat* values = style seed
          style value = *PILE_L*
          seed = random number generator seed
       *tau* value = thermostat damping parameter (time unit)
       *scale* value = scaling factor of the damping times of non-centroid modes of PILE_L thermostat
       *iso* or *aniso* values = pressure (pressure unit)
         pressure = scalar external pressure of the barostat
       *barostat* value = *BZP* or *MTTK*
       *taup* value = barostat damping parameter (time unit)
       *fixcom* value = *yes* or *no*
       *lj* values = epsilon sigma mass planck mvv2e
          epsilon = energy scale for reduced units (energy units)
          sigma = length scale for reduced units (length units)
          mass = mass scale for reduced units (mass units)
          planck = Planck's constant for other unit style
          mvv2e = mass * velocity^2 to energy conversion factor for other unit style

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all pimdb/nvt method pimd fmass 1.0 sp 1.0 temp 2.0 nhc 4
   fix 1 all pimdb/langevin ensemble npt integrator obabo temp 113.15 thermostat PILE_L 1234 tau 1.0 iso 1.0 barostat BZP taup 1.0

Description
"""""""""""

.. versionchanged:: -

Fix *pimdb/nvt* and fix *pimdb/langevin* were added, inheriting fix *pimd/nvt* and fix *pimd/langevin*, respectively.

These fix commands are based on the fix *pimd/nvt* and fix *pimd/langevin* commands for 
performing quantum molecular dynamics simulations based
on the Feynman path-integral formalism. The key difference is that fix *pimd/nvt* and fix *pimd/langevin* simulate *distinguishable* particles,
while fix *pimdb/nvt* and fix *pimdb/langevin* support simulations of bosons by including exchange effects.
The *pimdb* commands share syntax with the equivalent *pimd* commands. The user is referred to the documentation of the *pimd* commands for a 
detailed syntax description and additional, general capabilities of the commands.

.. note::

   Currently, fix *pimdb/langevin* only supports the "pimd" method, and fix *pimdb/nvt*
   only supports the "pimd" and "nmpimd" methods.

The isomorphism between the partition function of :math:`N` bosonic quantum particles and that of a system of classical ring polymers
at inverse temperature :math:`\beta`
is given by :ref:`(Tuckerman) <book-Tuckerman>`:

.. math::

   Z \propto \int d{\bf q} \cdot \frac{1}{N!} \sum_\sigma \textrm{exp} [ -\beta \left( E^\sigma + V \right) ].

Here, :math:`V` is the potential between different particles at the same imaginary time slice, which is the same for bosons and
distinguishable particles. The sum is over all permutations :math:`\sigma`. Recall that a permutation matches each element :math:`l` in :math:`1, ..., N` to an element :math:`\sigma(l)` in :math:`1, ..., N` without repetitions. The energies :math:`E^\sigma` correspond to the linking of ring polymers of different particles according to the different permutations:

.. math::

   E^\sigma = \frac{\sqrt{Pm^2}}{2\beta \hbar} \sum_{l=1}^N \sum_{j=1}^P \left(\bf{r}_l^j - \bf{r}_l^{j+1}\right)^2,

where :math:`P` is the number of beads and :math:`\bf{r}_l^{P+1}=\bf{r}_{\sigma(l)}^1.` 

Hirshberg et. al. showed that the ring polymer potential 
:math:`-\frac{1}{\beta}\textrm{ln}\left[ \frac{1}{N!} \sum_\sigma \textrm{exp} ^ { -\beta  E^\sigma } \right]`, which scales exponentially with :math:`N`, 
can be replaced by a potential :math:`V^{1,N}` defined through a recurrence relation :ref:`(Hirshberg) <Hirshberg>`:

.. math::

   \textrm{exp} ^ { -\beta  V^{[1,N]} } = \frac{1}{N} \sum_{k=1}^N \textrm{exp} ^ { -\beta \left(  V^{[1,N-k]} + E^{[N-K+1,N} \right)}.

Here, :math:`E^{[N-K+1,N}` is the spring energy of the ring polymer obtained by connecting the beads of particles :math:`N − k + 1, N − k + 2, ..., N` in a cycle.
This potential does not include all :math:`N!` permutations, but samples the same bosonic partition function. The implemented algorithm in LAMMPS for calculating 
the potential is the one developed by Feldman and Hirshberg, which scales like :math:`NP+N^2` :ref:`(Feldman) <Feldman>`. 
The forces are calculated as weighted averages over the representative permutations,
through an algorithm that scales the same as the one for the potential calculation :math:`NP+N^2` :ref:`(Feldman) <Feldman>`.

Output
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Fix *pimdb/nvt* computes a global 4-vector, which can be accessed by
various :doc:`output commands <Howto_output>`.  The three quantities in
the global vector are:

   #. the total spring energy of the quasi-beads,
   #. the current temperature of the classical system of ring polymers,
   #. the current value of the scalar virial estimator for the kinetic
      energy of the quantum system :ref:`(Herman) <Herman>`.
   #. the current value of the scalar primitive estimator for the kinetic
      energy of the quantum system :ref:`(Tuckerman) <book-Tuckerman>`.

The vector values calculated by fix *pimdb/nvt* are "extensive", except for the
temperature, which is "intensive".

Fix *pimdb/langevin* computes a global 6-vector, which
can be accessed by various :doc:`output commands <Howto_output>`. The quantities in the global vector are:

   #. kinetic energy of the normal mode
   #. spring elastic energy of the normal mode
   #. potential energy of the bead
   #. total energy of all beads
   #. primitive kinetic energy estimator
   #. virial energy estimator

The first 3 are different for different log files, and the others are the same for different log files.

----------

.. book-Tuckerman:

**(Tuckerman)** M. Tuckerman, Statistical Mechanics: Theory and Molecular Simulation (Oxford University Press, 2010)

.. Hirshberg:

**(Hirshberg)** B. Hirshberg, V. Rizzi, and M.Parrinello, “Path integral molecular dynamics for bosons,” Proc. Natl. Acad. Sci. U. S. A. 116, 21445 (2019)

.. Feldman:

**(Feldman)** Y. M. Y. Feldman and B. Hirshberg, “Quadratic scaling bosonic path integral molecular dynamics,” J. Chem. Phys. 159, 154107 (2023)

.. _Herman:

**(Herman)** M. F. Herman, E. J. Bruskin, B. J. Berne, J Chem Phys, 76, 5150 (1982).
