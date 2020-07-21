Setup
-----

CLI
===

Most features of pdom can be accessed directly from its command-line tools.
The easiest way to install them on your system is via `pipx <https://pipxproject.github.io/pipx/>`_.
`pipx` is a packagemenager for tools written in python that helps to keep them isolated and up to date.
Some notes on how to get `pipx` running can be found at the end of this section.
If you just need the CLI, use the following line to install pdom.

.. code-block:: shell-session

   $ pipx install pdom



Library
=======

You can also install the :mod:`pdom` library directly through the Python Package Index (`PyPI <https://pypi.org>`_) for use in your own projects.
The use of a `virtual environment <https://docs.python.org/3/tutorial/venv.html>`_ is recommended.

.. code-block:: shell-session

   $ pip install pdom


If the stable version of pdom on PyPI is missing a particular function, you can install the latest development version directly from the GitHub repository.

.. code-block:: shell-session

    $ pip install -U git+https://github.com/theia-dev/pdom.git#egg=pdom


Developing and documentation 
============================

To work with the source code clone the repository from GitHub and install the requirements.
The source code is accompanied by the documentation and an extensive collection of test cases.

.. code-block:: shell-session

    $ git clone https://github.com/theia-dev/pdom.git
    $ python3 -m venv pdom_env
    $ . pdom_env/bin/activate
    (pdom_env) $ pip install --upgrade pip
    (pdom_env) $ pip install -r pdom/requirements.txt

Building the documentation locally needs a few extra python packages. 
They can also be installed in the same virtual environment with the following command.

.. code-block:: shell-session

    (pdom_env) $ pip install -r pdom/docs/requirements.txt

The HTML version of the documentation can then be built:

.. code-block:: shell-session

     (pdom_env) $ sphinx-build -b html pdom/docs pdom/docs_html

The tests are located under ``pdom/tests`` can be started with through:

.. code-block:: shell-session

      (pdom_env) $ python -m unittest discover -s pdom/tests


pipx
====

**Under macOS:**

For macOS the `Homebrew <https://brew.sh>`_ package manager is the easiest way to install pipx.

.. code-block:: shell-session

      $ brew install pipx
      $ pipx ensurepath

**Under Linux:**

For some distributions the python package system `pip` is not installed by default.
On Debin/Ubuntu systems it can be quickly installed.

.. code-block:: shell-session

      $ sudo apt update
      $ sudo apt install python3-pip

Then `pipx` can be added.

.. code-block:: shell-session

     $ python3 -m pip install --user pipx
     $ python3 -m pipx ensurepath

**Under Windows:**

Python is not installed by default under Windows.
You can get the installer from the `Python download page <https://www.python.org/downloads/>`_.
The python package system `pip` is already included in the latest releases.

In the windows commandline `pipx` can then be installed.

.. code-block:: shell-session

     $ python3 -m pip install --user pipx
     $ python3 -m pipx ensurepath

For more information on `pipx` refer to its `documentation <https://pipxproject.github.io/pipx/>`_.