Usage
=====

Here you can find a quick summary about :mod:`pdom` command-line tools and the :mod:`pdom` library.
A more in-depth overview can be found in the :ref:`ref_examples` section

CLI pdom
--------

.. argparse::
   :module: pdom.cli
   :func: get_parser_simulation
   :prog: pdom

CLI pdom.config
---------------

.. argparse::
   :module: pdom.cli
   :func: get_parser_config
   :prog: pdom.config


How to create a config is described in the :ref:`ref_config` section.


Library
-------

.. code-block:: python

   import pdom

   simulation = Simulate('simple.ini')
   simulation.run()


See also the source code documentation of :class:`pdom.Simulate`.