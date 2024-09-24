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
   angle_coeff * 25.0 1.0 10.0

Description
"""""""""""

The *mwlc* angle style models a meltable wormlike chain, using a potential that is a canonical-ensemble superposition of
a non-melted and a melted state :ref:`(Farrell) <Farrell>`,

.. math::

    \beta E = -\log [q + q^{m}],

where

.. math::
    q = \exp [-l_{p}(1-\cos{\theta})/\sigma], \\
    q^{m} = \exp [-\beta\mu-l_{p}^{m}(1-\cos{\theta})/\sigma],

:math:`l_{p}` is the persistence length of the non-melted state,
:math:`l_{p}^{m}` is the persistence length of the melted state,
and :math:`\mu` is the melting energy.

This potential is a continuous version of the two-state potential
introduced by :ref:`(Yan) <Yan>`.

The following coefficients must be defined for each angle type via the
:doc:`angle_coeff <angle_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`l_{p}` (distance)
* :math:`l_{p}^{m}` (distance)
* :math:`\mu` (energy)

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

**(Farrell)** Farrell, Dobnikar, Podgornik, Curk, Phys Rev Lett, in production.

.. _Yan:

**(Yan)** Yan, Marko, Phys Rev Lett, 93, 108108 (2004).
