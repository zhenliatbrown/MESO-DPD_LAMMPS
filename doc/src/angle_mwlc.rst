.. index:: angle_style mwlc

angle_style mwlc command
==========================

Syntax
""""""

.. code-block:: LAMMPS

   angle_style mwlc

Examples
""""""""

.. code-block:: LAMMPS

   angle_style mwlc
   angle_coeff * 25 1 10 1

Description
"""""""""""

The *mwlc* angle style models a meltable wormlike chain, using a potential that is a canonical-ensemble superposition of
a non-melted and a melted state :ref:`(Farrell) <Farrell>`,

.. math::

    E = -k_{B}T\,\log [q + q^{m}] + E_{0},

where

.. math::
    q = \exp [-k_{1}(1+\cos{\theta})/k_{B}T], \\
    q^{m} = \exp [-(\mu+k_{2}(1+\cos{\theta}))/k_{B}T], \\
    E_{0} = -k_{B}T\,\log [1 + \exp[-\mu/k_{B}T]],

:math:`k_1` is the bending force constant of the non-melted state,
:math:`k_2` is the bending force constant of the melted state,
and :math:`\mu` is the melting energy.

This potential is a continuous version of the two-state potential
introduced by :ref:`(Yan) <Yan>`.

The following coefficients must be defined for each angle type via the
:doc:`angle_coeff <angle_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`k_1` (energy)
* :math:`k_2` (energy)
* :math:`\mu` (energy)
* :math:`T` (temperature)
----------


Restrictions
""""""""""""

This angle style can only be used if LAMMPS was built with the
EXTRA-MOLECULE package.  See the :doc:`Build package <Build_package>` doc page
for more info.

Related commands
""""""""""""""""

:doc:`angle_coeff <angle_coeff>`

Default
"""""""

none

----------

.. _Farrell:

**(Farrell)** `Farrell, Dobnikar, Podgornik, Curk, Phys Rev Lett, 133, 148101 (2024). <https://doi.org/10.1103/PhysRevLett.133.148101>`_

.. _Yan:

**(Yan)** `Yan, Marko, Phys Rev Lett, 93, 108108 (2004). <https://doi.org/10.1103/PhysRevLett.93.108108>`_
