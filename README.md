# Electrodynamics Research
The default location for my personal research that's more or less ready for outside consumption.  Right now this consists of four main topics, each updated more or less manually as sections are completed.  I've chosen to put this on GitHub as while it's a terrible place to demonstrate the math behind all of this stuff, it's a great place to share code; and, some of what I've done in the code is relevant beyond my personal interests.  Where possible or where I was interested enough to type it up, I've included with each topic a longer form paper describing the math / physics / theory of what's been put in code.  Hopefully some of it makes sense!

The topics so far include:

1. Non-Relativistic Quantum Mechanics,
2. Perfectly Matched Layer for a Standard Wave Equation,
3. Quantum Eigenvectors Arising from Discretized Quantum Mechanics,
4. Fundamental Solutions for a Standard Wave Equation.

Those aren't particularly descriptive titles in and of themselves, but--except for the first--they all have longer papers (in .pdf) in each folder describing my overall thought process.  I have no clue if any of it's correct theoretically, but it works and I got the results that the larger theory predicted so... progress!

## Non-Relativistic Quantum Mechanics (QM)
Over the past few years I've been slowly poking away at figuring out how time-domain versions of non-relatistic quantum mechanics (i.e. Schrodinger's Equation, Pauli's Equation, Electromagnetic Potential Wave Equations) can be shoved into a computer, and to see if we can get some useful results out of it.  This is part of a larger thesis--very much a work in progress--but the completed code I figure I'd upload somewhere just because.

There should be a separate README.md in this folder describing how things can be used, but this contains a few things: 

- General classes that integrate Schrodinger and Pauli's equations in one or two dimensions via a specific, implicit, finite difference method (Thomas' or the tri-diagonal matrix algorithm);
- A couple of test files specific to one or two dimensional potentials; the one-dimensional ones are more complete,
- A whole bunch of standalone electromagnetic scalar potentials, magnetic fields, helper functions, and
- A bunch of common functions that are useful.

Oh, for 2D problems, I've written a mostly optimized version of a tdma solver in C found in `tdma.c`, deep in `/+CommonFunctions/C Source Files/`.  You'll probably have to compile this one your own and place in the `/+CommonFunctions/` folder instead of whatever's in there as `tdma.mexw64`.  To compile it, please use `mex -R2018a 'CFLAGS=-mavx' tdma.c`.

## Perfectly Matched Layer
How do you add decent absorbing boundary conditions so that you can pretend you're simulating real electromagnetic phenomenon except inside of a computer?  How do you do this when you're not solving Maxwell's equations, but wave equations for potentials and not fields?  Well look no further, 'cause some ongoing work is being done all up in this folder.  So far this includes a, "standard," analytic continuation of spatial coordinates into the complex domain, and then discretized and solved using a few different techniques:

- A fully explicit finite difference method using first order equations via an auxiliary differential equation,
- A fully explicit finite difference method using second order equations,
- A semi-implicit finite difference method using first order equations via an auxiliary differential equation.

The nice thing about these methods is that the exact same files should work exactly the same in 3D (albeit quite slow and memory intensive).  Since my larger research requires absolute stability, hopefully a fully implicit method will appear soon.

## Quantum Eigenvectors

## Wave Equation
