## Motivation
Data gained from Molecular Dynamics (MD) and Markov Chain Monte Carlo processes (MD) are highly autocorrelated, thus we can not imply statistic developed for uncorrelated data. The aim of this package is to give simple means to get the mean, the margin of error and effective sample size of correlated data. One should note that some heuristic and educated guesses is needed as you have to find the 'slowest' evolving parameter of your system, check the data for stationarity. In the case of simulation, probably, your system has to be equilibrated.
Thus the package covers most basic cases of 1d stationary autocorrelated process.

There are multiple sampling routines, i.e. binning analysis, but for different reason no unified implementation found. Also every problem is too specific to come up with the universal solution, consider using the package for simple cases.

## Sample to target routine

The routine samples an observable in an infinite loop till one of the targets reached. Then returns the mean, the margins of error and  the effective sample size.

<img src="sample_to_target_scheme.svg">

### Autocorrelation time
Autocorrelation time is found by numerical integration of autocorrelation function.
References:...

### Mean, margins of error and effective sample size

### Sampling targets

## How to use

### Callback

### Other kwargs