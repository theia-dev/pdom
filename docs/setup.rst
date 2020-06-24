Setup
-----

CLI
===

Most features of pdom can be directly accessed from its command-line tools.
The easiest way to install them on your system is via `pipx <https://pipxproject.github.io/pipx/>`_.
So if you just need the CLI, use the following line to install pdom.

.. code-block:: bash

   pipx install pdom

Library
=======

You can also install the :mod:`pdom` library directly through the Python Package Index (`PyPI <https://pypi.org>`_) for use in your own projects.
The use of a `virtual environment <https://docs.python.org/3/tutorial/venv.html>`_ is recommended.

.. code-block:: bash

   pip install pdom


If the stable version of pdom on PyPI is missing a particular function, you can install the latest version directly from the GitHub repository.

.. code-block:: bash

    pip install -U git+https://github.com/theia-dev/pdom.git#egg=pdom


Developing and documentation 
============================

To work with the source, code clone the repository from GitHub and install the requirements.
The source code is accompanied by the documentation and an extensive collection of test cases.

.. code-block:: bash

    git clone https://github.com/theia-dev/pdom.git
    python3 -m venv pdom_env
    . pdom_env/bin/activate
    pip install --upgrade pip
    pip install -r pdom/requirements.txt

Building the documentation locally needs a few extra python packages. 
They can also be installed in the same virtual environment with the following command.

.. code-block:: bash

    pip install -r pdom/docs/requirements.txt

The HTML version of the documentation can then be built with:

.. code-block:: bash

     sphinx-build -b html pdom/docs pdom/docs_html

The tests are located under ``pdom/tests`` and can be started with:

.. code-block:: bash

      python -m unittest discover -s pdom/tests

