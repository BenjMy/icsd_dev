Getting Started
===============

The getting-started should show some primary use cases in more detail. The reader will follow a step-by-step procedure to set-up a working prototype

ICSD aims to process Mise-à-la-masse (MALM) datasets for a variety of applications. ICSD has been initially develop for plant root imaging using a python application programming interface

The simpliest processing can be achieved with the python API::

    from icsd3d import Problem
    k = Problem()
    k.createSurvey('test/coverCrop.csv')
    k.invert() # invert measurements
    k.showResults() # display inverted pseudo-section

More examples are available in the Example pseudo-section.