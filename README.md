# Self-weight-linear-programming

As per findings of [ArtOfScience](https://github.com/artofscience/SAOR), we implement LP+AML in the standard [DTU Python minimum compliance code](https://www.topopt.mek.dtu.dk/Apps-and-software/Topology-optimization-codes-written-in-Python), adapted for the self-weight problem.

Standard compliance code excecution

`python3 topopt_cholmod.py 180 60 0.4 5.4 3.0 1`

The *raw_input()* command with which it terminates has been replaced by *input()* for Python 3 compatibility. Tested on Ubuntu 20.04 LTS.

Self-weight version of compliance problem solved with LP ([OSQP](https://osqp.org/docs/index.html)) and the adaptive move-limit implied by the [MMA](https://people.kth.se/~krille/mmagcmma.pdf), excecuted with

`python3 topopt_sw_cholmod.py 100 100 0.2 5.4 3.0 1`

Although a clear slow-down is noticed *vs* the much simpler standard minimum compliance OC method, [OSQP](https://osqp.org/docs/index.html) performs decently as an LP solver, given that the author has not managed to make the equivalent SCIPY LP solver perform at all; see *e.g.* 

`python3 topopt_sws_cholmod.py 100 100 0.2 5.4 3.0 1`

