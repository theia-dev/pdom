.. _ref_config:

Configuration
=============

Simulation
----------

To run the simulation :mod:`pdom` needs a couple of information from the user.
For :mod:`pdom`, this data is saved in an ``.ini`` file.
This ensures that the results of a simulation can always be accompanied by the parameters that created them.

A simple config for :mod:`pdom` looks like this:

.. literalinclude:: examples/simple.ini
   :language: ini
   :caption: simple.ini

The parameters are arranged in different sections.
Not all sections need to exist for each simulation.

To create a new ``.ini`` file you can use the configuration tool of :mod:`pdom`.

.. code-block:: shell-session

   $ pdom.config --out 'my_new.ini'

It will guide you through the process by collecting all relevant information.
For the parameters of the initial molecule, it is helpful to look up its chemID in `PubChem <https://pubchem.ncbi.nlm.nih.gov>`_ before you start the process.
This enables the automatic gathering of the molecule parameters from this database.

Most of the time, the automatic process should suffice.
But if you want to build your ``.ini`` file from scratch, take a look into the available :ref:`ref_config_setting`.


Experimental data
-----------------

To extract parameters, we need to compare experimental data to the simulation.
This data needs to be provided in a structured way.
For :program:`pdom` we use a :file:`.json` file.
Depending on the type of fit you want to carry out, the available features differ slightly.

Adsorption-Desorption experiments
+++++++++++++++++++++++++++++++++

Adsorption-Desorption experiments in the dark are analyzed with the single-species model.
As it is common to have multiple repetitions that are based on the same setup.
:program:`pdom` supports multiple time series in its fits.
The initial concentration and time steps can be different between the series.
Below are examples of such a data file.

.. literalinclude:: examples/data_ads_des_multi.json
   :language: json
   :caption: ads_des_multi.json


Degradation experiments
+++++++++++++++++++++++

For degradation experiments all time series have to start with the same initial concentration set in the :file:`.ini` file.
This is the data file from the example :ref:`ref_ex_reac_fit`.

.. literalinclude:: examples/fit_reac/example_reac_fit.json
   :language: json
   :caption: fit_reac_mutli.json

Multi species model
+++++++++++++++++++

To compare multi-species model simulations to experiments, TOC (or NPOC) can be used.
In general, the fit to TOC data is the last step in the experiment analytics.
This fit is limited to a single TOC curve, due to the usually demanding experimental process.
If multiple TOC experiments are available, the results should be average before simulation.
As just one initial concentration is needed, it is taken from the :file:`.ini` file.

This is the data file from the example :ref:`ref_ex_toc_fit`.

.. literalinclude:: examples/toc_fit/example_toc_fit.json
   :language: json
   :caption: fit_toc.json

.. _ref_config_setting:

Configuration settings
----------------------

.. include:: config_setting_details.rst
