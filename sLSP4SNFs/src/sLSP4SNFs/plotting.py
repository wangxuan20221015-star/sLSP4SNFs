import matplotlib.pyplot as plt
from .config import  LABEL_SIZE, set_mpl_style, set_default_colormap

class Plot:
    """Visualization helper for SLSP analysis."""
    
    def __init__(self,slspObject):
        set_mpl_style()
        self.obj = slspObject
        self.custom_cmap = set_default_colormap()
        self.fig = None
        self.axs = {}
        
        
        
    def LSP(self,ax=None):
        """
        Plot the global Lomb–Scargle periodogram around the target frequency.

        Parameters
        ax : matplotlib.axes.Axes, optional
            Axes object to plot on. If None, a new figure and axes are created.
        """
        self._check_status()
        
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 4))
            ax.set_xlabel(r"Frequency ($\mu$Hz)",fontsize=LABEL_SIZE)
            ax.set_ylabel("Amplitude (ppt)",fontsize=LABEL_SIZE)
        
        ax.plot(
            self.obj.frequency - self.obj.targetFreq,
            self.obj.amplitude,
            color='k',
            lw=2
        )
        ax.set_xlim(-self.obj.deltaF ,self.obj.deltaF )
        ax.set_ylim(0.0001,self.obj.targetAmp*1.6)
        
        return ax
        
        
        
    def sLSP(self,ax=None):
        """
        Plot the sliding Lomb–Scargle periodogram (time–frequency map).

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes object to plot on. If None, a new figure and axes are created.
        """
        self._check_status()
        
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 8))
            ax.set_xlabel(r"Frequency ($\mu$Hz)",fontsize=LABEL_SIZE)
            ax.set_ylabel("Time (BJD)",fontsize=LABEL_SIZE)
        
        
        ax.scatter(
            self.obj.dfZoom['frequency'] - self.obj.targetFreq,
            self.obj.dfZoom['time'],
            c=self.obj.dfZoom['amplitude'],
            cmap=self.custom_cmap,
            s=80,
            alpha=0.8
        )
        ax.set_xlim(-self.obj.deltaF,self.obj.deltaF)
        ax.set_ylim(self.obj.referTimes[0],self.obj.referTimes[-1])
        
        return ax
        
    def SNF(self):
        """
        Generate the full SNF diagnostic figure.
        """
        self._check_status()
        
        self._canvas()
        self._plot_left_panels()
        self._plot_right_panels()
        self._plot_theory()
        
        #self._annotate_snf_diagnostics()
        
    def _check_status(self):
        if not hasattr(self.obj, 'snfState'):
            raise RuntimeError("Analysis not complete. Please run analyze_snf() first.")
        
    
    def _canvas(self):
        """Define the layout for the 5-panel SNF diagnostic figure."""
        
        self.fig = plt.figure(figsize=(12, 9))
    
        # Axes layout (left, bottom, width, height)
        axesDef = {
            "ax1": [0.10, 0.72, 0.32, 0.24],
            "ax2": [0.10, 0.10, 0.32, 0.60],
            "ax3": [0.45, 0.76, 0.40, 0.20],
            "ax4": [0.45, 0.56, 0.40, 0.20],
            "ax5": [0.45, 0.10, 0.40, 0.39],
        }
        
        for name, rect in axesDef.items():
            setattr(self, name, self.fig.add_axes(rect))
            
        self.axs = [self.ax1, self.ax2, self.ax3, self.ax4, self.ax5]
        for ax in self.axs:
            for spine in ax.spines.values():
                spine.set_linewidth(1.5)
    
        self._format_axes()
        
        
    def _format_axes(self):
        # Tick configuration
        self.ax1.tick_params(bottom=False, labelbottom=False)
    
        for ax in (self.ax3, self.ax4, self.ax5):
            ax.tick_params(
                left=False,
                labelleft=False,
                right=True,
                labelright=True,
            )
            ax.yaxis.set_label_position("right")
    
        self.ax3.tick_params(bottom=False, labelbottom=False)
    
        self.ax1.set_ylabel("Amplitude (ppt)", fontsize=LABEL_SIZE)
        self.ax2.set_xlabel(r"Frequency ($\mu$Hz)", fontsize=LABEL_SIZE)
        self.ax2.set_ylabel("Time (BJD)", fontsize=LABEL_SIZE)
        self.ax4.set_xlabel("Time (BJD)", fontsize=LABEL_SIZE)
        self.ax5.set_xlabel(r"Frequency ($\mu$Hz)", fontsize=LABEL_SIZE)
        self.ax5.set_ylabel(
            "Amplitude (nHz)",
            fontsize=LABEL_SIZE,
            rotation=270,
            labelpad=35,
        )
    
        
    def _plot_left_panels(self):
        self.LSP(self.ax1)
        self.sLSP(self.ax2)
    
        # Highlight ridge (center frequencies)
        self.ax2.scatter(
            self.obj.centerFreqs - self.obj.targetFreq,
            self.obj.centerTimes,
            c="limegreen",
            s=20,
            alpha=1.0,
            zorder=3,
        )
    
        
    def _plot_right_panels(self):
        
        # Right-top panels: SNF light curve
        self.ax3.scatter(self.obj.centerTimes,
            (self.obj.centerFreqs -self.obj.targetFreq) * 1000,   # nHz
            c="k", 
            s=10, 
            alpha=1
        )
        self.ax3.set_xlim(250, 1500)
        self.ax3.set_ylim(-39, 39)
    
        self.ax4.scatter(self.obj.centerTimes,
            (self.obj.centerFreqs - self.obj.targetFreq) * 1000,   # nHz
            c="k", 
            s=10, 
            alpha=1
        )
    
        
        # Right-bottom panel: SNF periodogram
        self.ax5.plot(
            self.obj.snfFitPg.frequency.value,
            self.obj.snfFitPg.power.value * 1000,
            c="k",
            lw=1.5,
        )
        self.ax5.set_xlim(0, 0.08)
        
    
    def _plot_theory(self):
        
        # Theoretical SNF curve (optional)
        self.ax3.plot(
            self.obj.snfSimulateTime,
            2 * self.obj.snfSimulateFreq,
            c="r",
            lw=2,
            zorder=2,
        )
        
    def _annotate_snf_diagnostics(self):
        """
        Annotate SNF diagnostic metrics on the periodogram panel.
        """
        ax = self.ax5

        
        ax.text(
            0.52, 0.90,
            f"R: {self.obj.criterion:.2f}",
            transform=ax.transAxes,
            fontsize=16,
            color="red",
            ha="left",
        )

        ax.text(
            0.52, 0.82,
            r"A: {:.2f} $\pm$ {:.2f}$\,\mathrm{{nHz}}$".format(
                self.obj.snfAmp * 1000,
                self.obj.sigemaFreq * 1e9 / 86400,
            ),
            transform=ax.transAxes,
            fontsize=16,
            color="red",
            ha="left",
        )

        ax.text(
            0.52, 0.74,
            r"P: {:.2f} $\pm$ {:.2f}$\,\mathrm{{days}}$".format(
                self.obj.snfPeriod,
                self.obj.sigmaPeriod,
            ),
            transform=ax.transAxes,
            fontsize=16,
            color="red",
            ha="left",
        )

        ax.text(
            0.52, 0.66,
            f"SNR: {self.obj.snfSnr:.2f}",
            transform=ax.transAxes,
            fontsize=16,
            color="red",
            ha="left",
        )

        ax.text(
            0.52, 0.58,
            f"SNF: {self.obj.snfState}",
            transform=ax.transAxes,
            fontsize=16,
            color="red",
            ha="left",
        )
        
        ylabel = r"$\mathit{f}-\bar{f}$ (nHz)"
        self.ax3.text(
            1.12,
            -0.4,
            ylabel,
            transform=self.ax3.transAxes,
            fontsize=LABEL_SIZE,
            rotation=270,
            va="center",
        )