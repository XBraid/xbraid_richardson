Here is a tar file containing all the xbraid code for automatic time adaption
and richardson .

Example ex-07 uses richardson to solve a simple ode. This is good for testing
richardson as there is no spatial errors. (ex-07 -help for example cmd line)

Example ex-02-adaptive uses richardson to adapt the time grid automatically
for whatever problem ex-02 in xbraid was solving. (see help for options)

Example drug_uptake solves a coupled ODE describing drug uptake in the body.
That problem has a lot of spikes. That example uses richardson to adapt the
time grid.

Example ex-02-richardson should use richardson to solve problem. In that case,
when running with richardson, I think we hit spatial discretization error
pretty quickly. See comments from " -help"


===============
Jacob's comment
===============

Looks like ex-02-adaptive sets the option to get the Richardson estimate and 
refine based on that. 

You could just change it to return the Richardson error estimate, rather than
adding it to the solution

Could use his examples (especially ex-07) to regression test Richardson
- or add adaptive capabilities to ex-02 from ex-02-adaptive?


