# fitTAS
A fortran90 program that reproduces an experimental Transient Absorption Spectrum (TAS) from a linear combination of quantum-chemical calculated spectra

1) Program description

Having the experimental TAS for a given time delay or a specific EADS
(Evolution-Associated Difference Spectra) and the selected calculated
ones for each stationary point of the potential energy surfaces, the
program finds the best linear combination of the calculated spectra to
replicate the experimental one. This is done by generating random sets of
initial coefficients, and then optimizing them by following the gradient
that minimizes the difference with respect to the experimental spectrum.
Multiple initial random spawning maximizes the probability of finding the
global minimum (best fit) in the space spanned by the calculated spectra.

2) How the program works

In first place, the program will ask for the spectral window to be analyzed.
This is important because there may be some sections of the experimental
spectrum that should be avoided in the analysis as they contain spectral
interferences (stimulated emission, solvent absorption or Raman scattering...)
The spectral window is provided introducing the lower and upper wavelengths,
which, for coding reasons, must be integers of 5.

In the second step, the program reads the experimental spectrum, and will
stop if the first value (the lowest wavelength) is larger than the selected
minimum wavelength, or if the last value (the highest wavelength) is lower
than the selected maximum wavelength, at least one point is required outside
each boundary for the next step. If the spectrum is correct, then the program
will give it a proper format for the analysis: it will transform the spectrum
to a uniform stepsize of 5 nm by linear interpolation of the two closest values
for each wavelength, and will normalize it so the maximum intensity, in absolute
value, is equal to one. The program can reproduce the effects of stimulated
emission as it can handle negative intensities.

Then, the program will read sequentially the different theoretical spectra.
Once more, it will check if each spectrum contains enough points for the
formatting process, and will stop if this is not accomplished. These spectra
can be generated from the vertical excitation energies and oscillator strengths
(negative values for fluorescence are tolerated) with external programs like
Gabedit or so.

From an initial set of randomly generated coefficients (conditions), the program
calculates the gradient using a numerical procedure (calculates the error variation
around the coefficient x using x+0.01 and x-0.01), and updates the coefficients
following the gradient direction. The gradient norm is normalized to 0.01, in case
it takes larger values, to prevent excessive step sizes. This provides better
accuracy without affecting the simulation time (in fact, larger norms were tested
and either provided less accurate results or took longer to furnish them). The
optimization stops once the second derivative of the error parameter gets lower
than the convergence criterion, or if the maximum number of steps is reached; both
parameters are provided by the user and can be freely modified (we recommend a
convergence criterion of 10**-4 or 10**-5, and a maximum number of steps of at
least 10**3). Once convergence is reached, the program stores and prints the
optimized results for each set of initial random conditions.

3) Files

3.1 Input files

Only the experimental spectrum and all the theoretical spectra are needed. However,
it is recommended to create an input file to feed the program, so it can be modified
easily if you need to make changes in the parameters.

3.2 Output files

3.2.1 FormattedExpSpectrum.txt, it contains the experimental spectrum with the
format used in the fitting procedure, and only in the analyzed spectral window.

3.2.2 FormattedCalcSpectra.txt, similar to the previous one but for the theoretical
spectra. The header, ignored by gnuplot, identifies each individual spectrum.

3.2.3 FitSummary.txt, a summary that displays the error parameter and the final coefficients
for the linear combinations obtained from each of the initial random conditions.

3.2.4 fitXX.txt, N files which contain the optimized spectra generated from each
initial random coefficients.

3.2.5 fitXX.opt, N files that contain the intermediate coefficients during the
optimization for each initial set.

4) Examples

A set of input and output files is provided alongside the source code so the users
can check whether the compilation was successful and the program works as expected.
In this example, there is a unique linear combination of the calculated absorption spectra
(calc1-3.dat) that produces the best replica of the experimental TAS (exp.dat). If, after
compiling and running the program, you obtain the same output files as provided, that
indicates the program is working properly.

5) How to cite in scientific works

Aside for providing the proper reference to the source code, an additional citation to
the original manuscript in which we described the use of the first version of this code
is greatly appreciated:

E. M. Arpa, M. M. Brister, S. J. Hoehn, C. E. Crespo-Hernández, I. Corral,
J. Phys. Chem. Lett. 2020, 11, 5156-5161.
