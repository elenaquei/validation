1 Spet 2017

This program validates a branch of periodic orbits of a given ODE system.
A first solution is computed, then the branch is continued by one step ant the segment 
in-between is validated. If the validation succeeds, a new step is performed, if the
validation failed, the segment is recomputed doing a smaller step size (and eventually
adding nodes to the solution).


This library consists of a variety of functions whose aim is to validate a point wise solution, or a segment or a branch of solutions, considering just periodic solutions to polynomial ODEs.

The more important functions you should be concerned with are

continuation (a simple function to help you start your continuation proof, require a non-square problem, including scalar equations, and a single solution as in instance of Xi_vector)

radii_polynomails_cont_new (a more involved function validating a segment between two approximate solutions)


In order to help you approach this two main functions, a couple of helpers functions are there for you:

from_string_to_polynomial_coef is a function that takes as input a string defining a vector field and construct an instance of the class polynomial_coef that corresponds to that string (please view the help function to write the string according to the specified format)

time_series2Xi_vec is a function that takes as input a period and a time series and returns an instance of the class Xi_vector to work with

default_scalar_eq takes as input the constructed Xi_vector and creates an instance of the class scalar_eq corresponding to fixing the first variable of your Xi_vector (i.e. x(0) = const)

Remark: this specific code is built in such a way to start using Intlab, therefore Intlab must be in the path. 



For a first approach to the library, check out the file template_vdP, a script showing you how to use the functions just highlighted to get a point wise validation. 
Furthermore, template_lor_cont shows you how to validate a branch of solution in the simplest way.
A last example provided is example_Hopf, where the code is used to validated a Hopd bifurcation in the simplest way. 