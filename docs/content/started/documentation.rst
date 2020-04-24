The Mise-à-la-masse method
==========================

Background
----------

The forward operator
--------------------

Green functions generation
--------------------------

Notion of virtual sources
~~~~~~~~~~~~~~~~~~~~~~~~~

Problem linearisation
---------------------

.. math:: Ax=b

Where :math:`\textbf{A}` is a matrix, its columns are the simulated VRTe
R sequences; x is a vector containing the unknown VRTe weights;
:math:`\textbf{b}` is a vector containing the measured sequence of
resistance. Each row in A corresponds to the relative R in the
acquisition sequence, e.g., :math:`\textbf{A}_{1,1}` is the first
resistance extracted from the potential field simulated with injection
at the first VRTe.

= min{x̍̍}

Including Constraints
~~~~~~~~~~~~~~~~~~~~~

The charge conservation is implemented by appending a row of 1’s to
:math:`\textbf{A}` and a corresponding 1 to the vector
:math:`\textbf{b}`. This forces the sum of the VRTe weights to be equal
to 1.

More in including constrainst in Menke W (1989) geophysical data
analysis: discrete inverse theory. International Geophysics Series.
Academic Press, New York.

Problem regularisation
----------------------

Data regularisation
~~~~~~~~~~~~~~~~~~~

Model regularisation
~~~~~~~~~~~~~~~~~~~~

Prior informations
^^^^^^^^^^^^^^^^^^

The initial model :math:`\textbf{m}_{0}` vector is implemented using the
simple misfit between a single source current and the measured data:

.. math::

   \label{eq:PriorObjFct}

       \[f_{1,i}\left(d_{m},\ d_{f,i}\right)=\left\|d_{m}-d_{f,i}\right\|^{2}
   \]

Other approaches can be considered as a 1st attempt to describe regions
of influences without having to go trough a complete inversion (see
`[ssec:OthersApproaches] <#ssec:OthersApproaches>`__.

Spatial regularization
^^^^^^^^^^^^^^^^^^^^^^

-  For the 2d case, since the problem is undetermined a first order
   spatial regularization is added (Menke, 1989). Rows are added to
   express the differences between adjacent VRTe, e.g., the row
   :math:`\left[\begin{matrix}1&-1&\ldots\\\end{matrix}\right]`\ is the
   difference between the first two VRTe weights. The differences are
   added for the entire VRTe grid and set to 0 by adding corresponding
   0’s to b.

-  NEW: a second order spatial regularization with differentiation
   between x and y directions to obtain to different matrices (of the
   same size) such as :math:`D_{x}` and :math:`D_{y}`

-  NEW: for the 3d case, a k-mean with 4 (or more) neighbors sources
   regularization can be used. In that case each source is weighted so
   the sum is equal to 0.

We used a linear solver from Python library, using a least square
inversion which in the current version minimized the following objective
function:

.. math::

   \label{eq:ObjFctFull}
       \[\widetilde{\mathbit{m}}=\ min \left\{\left\|Lr\right\|^{2} + \lambda(\alpha_{s}\left\|m-m0\right\|^{2}+ \alpha_{x}\left\|D_{x}(m-m0)\right\|^{2} + \alpha_{z}\left\|D_{z}(m-m0)\right\|^{2})\right\}
   \]

where :math:`\textbf{m}_{0}` is a reference model to which we believe
the physical property distribution should be close. Often
:math:`\textbf{m}_{0}` is chosen to be a constant average value. In that
case the initial model :math:`\textbf{m}_{0}` vector is implemented
using the simple misfit between a single source current and the measured
data:

Equation `[eq:ObjFctFull] <#eq:ObjFctFull>`__ also contains the
coefficients controlling weight of the relative smallness
:math:`\alpha_{s}`, and the regularization anisotropy wieigth
:math:`\alpha_{x}` and :math:`\alpha_{y}`\ respectively in x and y
directions. Equation `[eq:ObjFctFull] <#eq:ObjFctFull>`__ can be
rewritten as:

.. math::

   \label{eq:ObjFct2}
   \widetilde{m}=\ min\left\{{(Gm-d)}^TW_d(Gm-d)\ +\ \lambda{(m-m_0)}^TW_m(m-m_0)\right\}

Where

.. math::

   \label{eq:}
   W_{d}=L^{T}L

The trade-off between data misfit and solution regularization is
controlled by :math:`\lambda`. The numerical routine includes a “pareto”
functionality wherein regularization and model-to-measurement fit are
traded off while changing the regularization weight. The obtained set of
solutions can be used to construct the “pareto front” (L-curve), which
is a widely accepted way to estimate the optimum regularization weight
:raw-latex:`\cite{hansen1993insect}`.

The solution is further constrained by forcing the linear solver to seek
only positive VRTe weights (i.e., inequality constraint), as the
negative source of current is known to correspond uniquely to the return
electrode. The following equation can be use to solve the inversion
problem:

.. math::

   \label{eq:}
   m={(G^{T}W_{d}G\ +\ \lambda W_{m})}^{-1}(G^{T}W_{d}d\ +\ \lambda W_{m}m_{0})

by solving the system Am=b, with:

.. math::

   \label{eq:}
   A=(G^{T}W_{d}G\ +\ \lambda W_{m})

.. math::

   \label{eq:}
   b=(G^{T}W_{d}d\ +\ \lambda W_{m}m_{0})

Model appraisal
---------------

.. _ssec:OthersApproaches:

Other approaches
----------------

TO WRITE: vraisemblance fct from Binley’s article cite binley paper

---------------