# Self-weight-linear-programming

As per findings of [ArtOfScience](https://github.com/artofscience/SAOR), we implement LP+AML in the standard [DTU Python minimum compliance code](https://www.topopt.mek.dtu.dk/Apps-and-software/Topology-optimization-codes-written-in-Python), adapted for the self-weight problem.

Standard compliance code excecution

`python3 topopt\_cholmod.py 180 60 0.4 5.4 3.0 1`

The `raw_input` command with which it terminates has been replaced by `input` as per Python 3 compatibility.

Self-weight version of compliance problem solved with LP ([OSQP](https://osqp.org/docs/index.html)) and the adaptive move-limit implied by the [MMA](https://people.kth.se/~krille/mmagcmma.pdf), excecuted with

`python3 topopt\_sw\_cholmod.py 100 100 0.2 5.4 3.0 1`
