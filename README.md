# untangling-complex-dynamical-systems-via-derivative-variable-correlations-srep05030-python

A python implementation of the network reconstruction algorithm from the paper

  _Untangling complex dynamical systems via derivative-variable correlations_, Z. Levnajić & A. Pikovsky, Scientific Reports 2014


Introduction
------------------

This project contains

1. implementation of the general algorithm for network reconstruction described in the paper, and 
2. sample data series and for example 1 from the same paper.

The method is well described in the paper therefore I will not go into details here. The notation in the code closely follows the notation in the paper.

Many thanks to Zoran Levnajic and Marc Gray who helped me to implement their method in python. Marc provided me with source implementation in Fortran and C so if you are more interested in that please leave me a comment and I will help you get in touch with him.


Requirements
---------------
* python 2.7
* [numpy and the scipy stack with matplolib ](https://www.scipy.org/install.html)

And that is it!! I tried to have as few dependencies as possible.

Example 1
---------------

Time series data and adjacency matrix are provided for example 1 of the paper. You can try this example by executing ``srep05030_example1.py''.

Here are the figures that you would obtain:

The time series used for network reconstruction:

![alt text](https://raw.githubusercontent.com/gajduk/untangling-complex-dynamical-systems-via-derivative-variable-correlations-srep05030-python/master/time_series.png)

After calling the reconustrction algorithm you get the following image, where you see the original and reconstructed adjacency matrices and the difference between them

![alt text](https://raw.githubusercontent.com/gajduk/untangling-complex-dynamical-systems-via-derivative-variable-correlations-srep05030-python/master/adjacency_matrix.png)

Not provided in this project are the codes that generate the time series for both example1 and example2 in the paper. These you can find in a [separate repository as matlab scripts](https://github.com/gajduk/egf-ras-mapk-pathway-ode-model/tree/srep05040-data-generation/src/srep05030)


Use the algorithm with your own data
---------------------------------------

If you have your own data and you want to apply this algorithm you will need to make sure that:

Your time points are uniformly spaced e.g. 0 minutes, 5 minutes, 10 minutes, 15 minutes etc.

The dynamics for each entity can be decribed as a *linear combination* of *known functions* i.e. Eq.1 from the paper

dx_i/dt = f(x_i) + Sum over k A_ij h(x_j)

where f(x) and h(x) are known. In biological context `h(x)' can often be modeled as Hill functions (see also example 2 from the paper).

If you don't know the original adjacency matrix and want to see which reconstructed A is best you would need to calculate the T score. Please consult the paper on how to do this, as this is not included in this project.

Closing remarks
-----------------------------

This is meant to serve as a starting point and provide a quick working example for network reconstruction - not a full fledged library for network reconstruction.

I hope this code helps someone!

If the constraint of only one ``h'' which you have to know beforehand is too restrictive I found this paper 

.. _Exact reconstruction of gene regulatory networks suing compressive sensing_, Chang, Gray and Tomlin, BMC Bioinformatics, 2015 

which allows you to specify a set of many different ``h'' functions and it predicts weights for all of them with just one run so check it out. I hope to put together a working example for this paper as well.
