import numpy as np
import pandas as pd
import logging
from lightkurve import LightCurve
from importlib.resources import files

logger = logging.getLogger(__name__)



class SLSP:
    """
    Sliding Lomb–Scargle Periodogram (sLSP) analysis class.

    Notes
    -----
    Unit conventions used throughout this class:

    - Time:
        Barycentric Julian Date (BJD), in units of days.

    - Frequency:
        Microhertz (µHz) for Lomb–Scargle periodograms.
        Nanohertz (nHz) for frequency residuals and SNF diagnostics.

    - Amplitude:
        Periodogram power is treated as amplitude.
        Internally stored in ppt (parts-per-thousand) unless otherwise noted.

    - Period:
        Derived from periodogram peak, reported in days.

    All unit conversions are handled explicitly in the code.
    """

    def __init__(self, time, flux, window=200, step=10):
        """
        Initialize the SLSP object and compute the global periodogram.

        Parameters
        ----------
        time : array-like
            Time array (BJD days).
        flux : array-like
            Flux or brightness measurements, and should be normalized by taking the median.
        window : float
            Time length of the sliding window (days).
        step : float
            Step size for window sliding (days).
        """ 
        self.time = np.asarray(time)      # BJD days   
        self.flux = np.asarray(flux)
        self.T = self.time[-1] - self.time[0]

        self.window = window    # days
        self.step = step       # days
        self.corDays = 0
        
        # Internal state variables
        self.t0 = None
        self.tf = None
        self.frequency = None
        self.amplitude = None
        self.df = None  # Full sliding LSP results
        self.snfFile = None

        self._sort_time_series()
        self._build_global_lsp()


    def _sort_time_series(self):
        """Sort the time series and define the time span."""
        idx = np.argsort(self.time)
        self.time = self.time[idx]
        self.flux = self.flux[idx]
        self.t0, self.tf = self.time[0], self.time[-1]


    def _build_global_lsp(self):
        """
        Construct the global Lomb–Scargle periodogram.
    
        Notes
        -----
        - Input light curve:
            time  : BJD (days)
            flux  : relative flux (dimensionless)
    
        - Output periodogram:
            frequency : µHz
            power     : Lomb–Scargle Amplitude (ppt)
    
        """
        
        self.lc = LightCurve(time=self.time, flux=self.flux)
        self.pg = self.lc.to_periodogram(freq_unit="uHz")
        self.frequency = self.pg.frequency.value                # µHz
        self.amplitude = self.pg.power.value * 1000                    # ppt
        

        
    def compute_slsp(self, oversampleFactor=10):
        """
        Compute the sliding Lomb–Scargle periodogram (sLSP).
    
        Parameters
        ----------
        oversampleFactor : int, optional
            Oversampling factor for the Lomb–Scargle periodogram.
    
        Notes
        -----
        - A rectangular time window is used.
        - Window length and step are in days.
        - For each window, a full Lomb–Scargle periodogram is computed.
    
        Output DataFrame (self.df):
        ---------------------------
        time       : window center time (BJD, days)
        frequency  : frequency grid (µHz)
        amplitude  : Lomb–Scargle power scaled to ppt
        """
        
        if self.tf - self.t0 < self.window:
            raise ValueError("Time span shorter than window length.")
        
        times, freqs, amps = [], [], []
        nSteps = int((self.tf - self.t0 - self.window) // self.step)

        for i in range(nSteps):
            t1 = self.t0 + i * self.step
            t2 = t1 + self.window
            tMid = t1 + self.window / 2

            mask = (self.time >= t1) & (self.time < t2)
            if not np.any(mask):
                continue

            try:
                lcPart = LightCurve(self.time[mask], self.flux[mask])
                pgPart = lcPart.to_periodogram(
                    freq_unit="uHz",
                    oversample_factor=oversampleFactor
                )
                times.extend([tMid] * len(pgPart.frequency))
                freqs.extend(pgPart.frequency.value)
                amps.extend(pgPart.power.value*1000)
            except Exception as e:
                logger.warning(f"sLSP failed at step {i}: {e}")

        self.df = pd.DataFrame({
            "time": times,
            "frequency": freqs,
            "amplitude": amps,
        })
        self.referTimes = np.unique(self.df["time"])     # BJD (days)
        

    
    def zoom_frequency(self, freq, deltaF=0.1):
        """
        Zoom into a narrow frequency band and extract the target frequency's amplitude.
    
        Parameters
        ----------
        freq : float
            Target frequency in µHz.
        deltaF : float, optional
            Half-width of the frequency window (µHz).
        """
        self._check_dependency("df", "compute_slsp()")
        self.targetFreq = freq    # µHz  
        self.deltaF = deltaF      # µHz

        mask = (
            (self.df["frequency"] >= freq - deltaF) &
            (self.df["frequency"] <= freq + deltaF)
        )
        self.dfZoom = self.df.loc[mask].copy()
        
        # Ridge extraction: max amplitude at each time
        grouped = self.dfZoom.groupby("time")
        idx = grouped["amplitude"].idxmax()

        self.centerTimes = self.dfZoom.loc[idx, "time"].to_numpy()         # BJD (days)
        self.centerFreqs = self.dfZoom.loc[idx, "frequency"].to_numpy()    # µHz
        
        # Global LSP amplitude near target frequency
        band = (
            (self.frequency >= freq - deltaF) &
            (self.frequency <= freq + deltaF)
        )
        self.targetAmp = np.max(self.amplitude[band])    # ppt
        
        
        

    def analyze_snf(self,freq,deltaF=0.1):
        """
        Evaluate whether a frequency satisfies the Super-Nyquist Frequency (SNF) criteria.
    
        Parameters
        ----------
        freq : float
            Target frequency in µHz.
        deltaF : float, optional
            Frequency window half-width for ridge extraction (µHz).
    
        Notes
        -----
        This method combines:
    
        1. Observational diagnostics:
            - Periodicity of the SNF light curve
            - Signal-to-noise ratio (SNR)
            - Amplitude threshold
    
        2. Theoretical SNF consistency:
            - Comparison with simulated Super-Nyquist frequencies
    
        All frequency residuals used in the SNF criterion are expressed in nHz.
        """

        self._check_dependency("df", "compute_slsp()")
        self.zoom_frequency(freq,deltaF)
        self._load_sim_data()
    
        
        
        """
        Calculates difference between observed ridge and expected Nyquist modulation
        Frequency residuals:
        - (centerFreqs - targetFreq): µHz
        - multiplied by 1000 → nHz
        - snfs: simulated Nyquist frequencies in nHz
        
        The factor of 2 accounts for the SNF reflection symmetry.
        """
        residual = np.abs(
                (self.centerFreqs - self.targetFreq) * 1000.0 - 2.0 * self.snfs
            )    # nHz
        
        self.criterion = np.round(np.nanmean(residual) / 32.0, 5)
            
            
        # Periodicity Analysis of the SNF curve.
        self.snfFitLc = LightCurve(
                                    self.centerTimes,   # BJD (days)
                                    self.centerFreqs    # µHz
                                    ).remove_nans().remove_outliers(sigma=3)
        
        self.snfFitPg = self.snfFitLc.to_periodogram(freq_unit="uHz")
        
        frequency = self.snfFitPg.frequency.value         # µHz
        amplitude = self.snfFitPg.power.value            # µHz

        self.snfNoise = np.nanmedian(amplitude[frequency <= 0.1])   # µHz
        self.snfAmp = amplitude.max()     # µHz
        self.snfSnr = self.snfAmp / self.snfNoise
        self.snfPeriod = self.snfFitPg.period_at_max_power.value * 1e6 / 86400   # days 
        
        self._estimate_errors()
        self._check_snf_state()
        self._print_summary()
        
        
    def _print_summary(self):
        """
        Outputs a human-readable message about the SNF status.
        """
        state = self.snfState
        
        print("\n" + "="*60)
        print("       sLSP4SNFs ANALYSIS REPORT")
        print("="*60)
        
        if state == "Confirmed":
            msg = (
                f"SUCCESS: The target at {self.obj.targetFreq:.2f} uHz is identified as a \n"
                f"Super-Nyquist Frequency (SNF). \n"
            )
        else:
            msg = (
                f"REJECTED: The target at {self.obj.targetFreq:.2f} uHz is NOT an SNF. \n"
            )
        
        print(msg)
        print("="*60 + "\n")
        
    
    def _check_snf_state(self):
        """Determine if the target behaves as a Super-Nyquist Frequency."""
        isPeriodValid = np.abs(self.snfPeriod - 372.0) < (5.0 * self.sigmaPeriod)
        isAmpStrong = (self.snfSnr >= 3.0) and (self.snfAmp > 0.01)
        isTheoryValid = (self.criterion <= 0.3)
        
        self.snfState = "Confirmed" if (isPeriodValid and isAmpStrong) or isTheoryValid else "Rejected"
        
    def _load_sim_data(self):
        """
        Load simulated Nyquist frequencies and match them to center_times.
        """
        data_path = files("sLSP4SNFs").joinpath("data","SNF_simulation_table.csv")
        
        if not data_path.exists():
            raise FileNotFoundError(f"Missing data file at: {data_path}")
        
        df = pd.read_csv(data_path,sep='\s+')
        
        self.snfSimulateTime = df["Time(BJD)"].values
        self.snfSimulateFreq = df["f_ny(nHz)"].values
        
        snfs = []
        for t in self.centerTimes:
            idx = np.argmin(np.abs(self.snfSimulateTime - t))
            snfs.append(self.snfSimulateFreq[idx])

        self.snfs = np.asarray(snfs)
        
        
    def _estimate_errors(self):
        """
        Estimate uncertainties of SNF parameters.
    
        Notes
        -----
        Error estimates follow standard analytical expressions for
        sinusoidal signals in white noise (e.g. Montgomery & O'Donoghue 1999, Zong et al. 2021).
    
        Definitions:
        ------------
        snfAmp        : peak amplitude of SNF periodogram (µHz)
        snfNoise      : median noise level of the periodogram (µHz)
        T             : total duration of the time series (days)
    
        Derived quantities:
        -------------------
        sigmaAmp     : amplitude uncertainty (same unit as snfAmp)
        sigmaFreq    : frequency uncertainty (1/day)
        sigmaPeriod  : period uncertainty (days)
    
        Unit conversions are handled explicitly.
        """
        self.sigmaAmp = self.snfNoise
        self.sigmaFreq = (
            np.sqrt(3) * self.sigmaAmp /
            (self.snfAmp * np.pi * self.T)
            )
        self.sigmaPeriod = self.snfPeriod**2 * self.sigmaFreq 
    
    
    def _check_dependency(self, attr, step):
        """Ensure prerequisite methods are run."""
        if not hasattr(self, attr):
            raise RuntimeError(f"Run `{step}` before this step.")
    


    def plot(self):
        """
        Return a Plot object for visualization.
        """
        from .plotting import Plot
        
        return Plot(self)