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
===============

Last time Master pulled in was commit 
7a35d4398f2822bc544bd843779b62a52ccd22de 
(April 13 2017, so a while ago)

The git log for this branch is in the file Branchs_Old_Gitlog


===============
===============

Ben said: 
From what i can tell, the ustore-dev branch added a feature to allow for global
refinement at the end of each v cycle. I have another branch that looks like it
implements the Richardson enhancments and the error estimator. 
 
It's not clear to me if Richardson relies on this ustore-dev branch that allows
for global refinement, or not.  Probably not, because Ben's richardson branch
(this one) just calles SetRFactor().  I don't think any special refinement at
the end of a cycle is needed.


===============
===============
Jacob's comments

Looks like ex-02-adaptive sets the option to get the Richardson estimate and 
refine based on that. 

You could just change it to return the Richardson error estimate, rather than
adding it to the solution

Could use his examples (especially ex-07) to regression test Richardson
- or add adaptive capabilities to ex-02 from ex-02-adaptive?


