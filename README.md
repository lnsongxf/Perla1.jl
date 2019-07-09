
# Overview

This repository has the code for [A Model of Product Awareness and Industry Life Cycles](http://jesseperla.com/papers/perla_1.pdf)

The [derivation document](https://github.com/jlperla/Perla1.jl/releases/download/dev/perla_1_computational_appendix.pdf) has the complete set of equations implemented for the model, where all equation numbers in the code refer to this document.

## Online Examples (no installation)
Using the Binder technology, the following links will execute the package notebooks in your browser
- [All Notebooks](https://mybinder.org/v2/gh/jlperla/Perla1.jl/master)
- [solving-large-linear-odes.ipynb](https://mybinder.org/v2/gh/jlperla/Perla1.jl/master?filepath=solving_large_linear_odes.ipynb): Methods for solving large systems of equations and markov chains (e.g., 100million by 100million)
- [transition-multiple-cohorts-example-staggered-with-forgetting-experiments.ipynb](https://mybinder.org/v2/gh/jlperla/Perla1.jl/master?filepath=transition_multiple_cohorts_example_staggered_with_forgetting_experiments.ipynb): Transition dynamics example demonstrating multiple cohorts
- [demand-function-with-multiple-cohorts.ipynb](https://mybinder.org/v2/gh/jlperla/Perla1.jl/master?filepath=demand_function_with_multiple_cohorts.ipynb): Solving the demand function for multiple cohorts
    

## Installation for Local Use

1. Follow the instructions to [install Julia and Jupyter](https://lectures.quantecon.org/jl/getting_started.html)

2. Open the Julia REPL (see the documentation above) and then install the package (by entering package mode) with

    ```julia
    ] add https://github.com/jlperla/Perla1.jl.git
    ```

3. There are several ways you can run the notebooks after installation

    Using the built-in Jupyter is straightforward.  In the Julia terminal
    ```julia
    using Perla1, IJulia
    notebook(detached=true, dir=dirname(dirname(pathof(Perla1))))
    ```

    Alternatively, to use a separate Jupyter installation you may have installed with Anaconda,
    ```julia
    using Perla1
    cd(dirname(dirname(pathof(Perla1))))
    ; jupyter lab
    ```
    where the last step runs your `jupyter lab` in the shell.

    **Note** In either case, the first time the `using` it will be very slow

4. The code for individual examples and experiences are in the `.ipynb` files

**NOTE:** When using the notebooks for the first time, it will be very slow as the package and its dependencies are all compiled.
