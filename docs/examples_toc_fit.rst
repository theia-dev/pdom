.. _ref_ex_toc_fit:

TOC fit
-------

To fit TOC data with :program:`pdom`, one of the *multi-species* models must be selected:

   * incremental
   * fragmentation
   * excess bonds

For this example, we will use *incremental*.
The model for the desorption constant :math:`k_{\mathrm{ads}}` is always *weak* when a TOC experiment is fitted.

The generation of the config file :download:`example_toc_fit.ini <examples/toc_fit/example_toc_fit.ini>` is again carried out with :program:`pdom.config`.
Lines with require user input are highlighted in yellow.

.. literalinclude:: examples/toc_fit/pdom.config.flow.txt
   :emphasize-lines: 2, 7, 13, 18, 19, 24, 28, 32, 36, 42, 47, 52, 56, 60, 64, 68
   :language: shell-session

After the config is generated, the experimental data set is created.
In this example values published by Houas (2001) :cite:`Houas2001` will be used.

.. literalinclude:: examples/toc_fit/example_toc_fit.json
   :language: json
   :caption: example_reac_fit.json

With both files prepared :program:`pdom` can be started.

.. code-block:: shell-session

    $ pdom example_toc_fit.ini --data example_toc_fit.json
    Start fitting to toc
       Iteration     Total nfev        Cost      Cost reduction    Step norm     Optimality
           0              1         2.6544e-01                                    1.85e+01
           1              2         1.4102e-01      1.24e-01       1.00e-02       8.09e+00
           2              3         4.6277e-02      9.47e-02       1.99e-02       2.48e+00
           3              4         6.3798e-03      3.99e-02       3.24e-02       3.75e-01
           4              5         3.1105e-03      3.27e-03       1.62e-02       3.24e-02
           5              6         3.0716e-03      3.89e-05       2.21e-03       1.61e-04
           6              7         3.0716e-03      9.76e-10       1.16e-05       4.15e-06
           7              8         3.0716e-03      0.00e+00       0.00e+00       4.15e-06
    `xtol` termination condition is satisfied.
    Function evaluations 8, initial cost 2.6544e-01, final cost 3.0716e-03, first-order optimality 4.15e-06.
    Fit finished
        k_ads: 3.000E-09 m/s
        k_des: 6.800E-03 1/s
        k_reac: 9.080E-02 1/s
        beta_0: -2.900E-02 1/s
        beta_1: 5.728E-01 1/s
        error: 1.005E-02
    Results saved in <your_working_dir>/example_toc_fit

The result of the fit is stored under :file:`{<your_working_dir>}/example_toc_fit/fit_toc.json`.

.. literalinclude:: examples/toc_fit/fit_toc.json
   :language: json
   :caption: <your_working_dir>/example_toc_fit/fit_toc.json

In the same folder, you find the raw data files with corresponding units.
The saved :download:`plot <examples/toc_fit/fit_toc.pdf>` shows the TOC development over time compared to the experimental results.

.. image:: examples/toc_fit/fit_toc.png
   :alt: example_toc_fit plot
