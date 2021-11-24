#acf - autocorrelation function
from statsmodels.tsa.stattools import acf
import scipy.stats
import numpy as np
import time
import logging
logger = logging.getLogger("sample_to_target_error")

def get_tau(x, acf_n_lags : int = 200):
    """Get the integrated autocorrelation time i.e. the distance between measurements
    when the data can be considered uncorrelated.

    Args:
        x (array): array of consecutive measurements
        acf_n_lags (int, optional): Number of lags calculated with ACF. Defaults to 200.

    Returns:
        float: integrated autocorrelation time
    """
    acf_arr = acf(x, nlags = acf_n_lags)
    #integrated autocorrelation time
    tau_int =1/2+max(np.cumsum(acf_arr))
    return tau_int

def correlated_data_mean_err(x, tau, ci = 0.95):
    """Returns the mean and the confidence interval of correlated sample.
    That distribution mean lies in sample mean +/- error with a probability of
    conference interval level (95% in most cases).
    Error = z * standard error of the sample

    Margins of the error calculated from standard error of a sample with
    inverse cumulative normal distribution (z-values) if the sample size is big
    (30+) or with a inverse cumulative t-distribution (t-values) for a small
    sample.

    Args:
        x (array): correlated sample
        tau (float): integrated autocorrelation time
        ci (float, optional): Confidence interval level. Defaults to 0.95.

    Returns:
        (float, float): sample mean and the margin of error
    """
    ##References:
    #http://www.stat.yale.edu/Courses/1997-98/101/confint.html
    #https://en.wikipedia.org/wiki/Confidence_interval
    x_mean = np.mean(x)

    #the sample is correlated so the effective size is smaller
    n_eff = np.size(x)/(2*tau)
    #print(f"Effective sample size: {n_eff}")
    if n_eff>30 and ci == 0.95:
        # the most common case when the sample size is big enough
        # we use normal distribution
        # z = ppf(1-(1-0.95)/2) = 1.96
        # where ppf - the percent point function or
        # the inverse cumulative distribution function
        # of normal distribution
        z = 1.96
    else:
        # otherwise we calculate it from t distribution
        # 1-(1-ci)/2 comes from the fact that we are excluding values
        # outside confidence interval from both tails of distribution
        z=scipy.stats.t.ppf(1-(1-ci)/2, n_eff)
    #print(f"z: {z}")
    err = np.std(x)/np.sqrt(n_eff) * z
    #print(f"mean: {x_mean};  error: {err}")
    return x_mean, err

def sample_to_target_error(
        get_data_callback,
        target_error = None,
        target_eff_sample_size = None,
        timeout = 30,
        initial_sample_size = 100,
        tau = None,
        ci = 0.95):
    """The routine sample an autocorrelated observable till one of the exit loop criterion meat.
    Possible criteria sampling time, sampling margins of error and sample size.
    The routine calculate mean value and margins of error in interations
    with ever increasing precision, in the next steps:
        1) get initial sample
        2) calculate autocorr time if not provided
        3) check if any of the criteria to end loop meet
        4) while no criterion meat
            5) double the sample size
            6) update autocorr time
            7) check end loop
        8) returns mean, margin of errors, sample size

    Note using autocorr time to calculate effective sample size.
    ----------------------------------------------------------------------------
    Requirements to callback function:
    use a function of the next signature f(n) -> List[float],
    where len(f(n)) = n
    ----------------------------------------------------------------------------
    For highly autocorrelated data (tau > 150) consider downsampling.
    Try downsample(x, dist) from the package.

    Args:
        get_data_callback ([type]): Sampling function, see requirements above.
        target_error (float, optional): Desired margins of error. Defaults to None.
        target_eff_sample_size (float, optional): Desired effective sample size. Defaults to None.
        timeout (int, optional): Timeout in seconds, in the worst case it takes two times longer. Defaults to 30.
        initial_sample_size (int, optional): Initial sample size. Defaults to 100.
        tau (float, optional): If autocorr time is known provide it to save some time. Defaults to None.
        ci (float, optional): Confidence interval. Defaults to 0.95.

    Returns:
        tuple: Returns mean value, margins of error and effective sample size
    """
    #stop sampling criteria
    def end_loop(elapsed_time, current_error, eff_sample_size):
        if elapsed_time > timeout:
            logger.info('Reached timeout')
            return True
        #check if kwarg provided
        elif target_error is not None:
            if current_error <= target_error:
                logger.info("Reached target error")
                return True
        elif target_eff_sample_size is not None:
            if eff_sample_size >= target_eff_sample_size:
                logger.info("Reached effective sample size")
                return True
        else:
            return False
    #init timer
    start_time = time.time()
    n_samples = initial_sample_size
    #first sample
    x = get_data_callback(n_samples)
    #if non autocorr time provided
    if tau is None: tau = get_tau(x)
    #correlated data mean and err margin
    x_mean, x_err = correlated_data_mean_err(x, tau, ci)
    #sampling loop
    while True:
        #time passed
        logger.info(f'Criteria are not reached')
        logger.info('More data will be collected')
        #append new sample (correlated)
        new_sample = get_data_callback(n_samples)
        x=np.append(x, new_sample)
        #if non autocorr time provided
        if tau is None: tau = get_tau(new_sample)
        #next time take twice the amount of data points
        n_samples = n_samples*2
        #err, sample size, time passed
        x_mean, x_err = correlated_data_mean_err(x, tau, ci)
        n_samples_eff = n_samples/(2*tau)
        elapsed_time = time.time() - start_time
        logger.info(x_err, n_samples_eff, elapsed_time)
        #stop sampling?
        if end_loop(elapsed_time, x_err, n_samples_eff):
            break
    #Done
    logger.info(f'Mean: {x_mean}, err: {x_err}, eff_sample_size: {n_samples/(2*tau)}')
    return x_mean, x_err, n_samples/(2*tau)

#for highly autocorrelated data consider using downsampling prior to sample_to_target routine
def downsample(x, dist):
    import random
    n_chunks = int(len(x)/dist)
    return [random.choice(chunk) for chunk in np.array_split(x,n_chunks)]