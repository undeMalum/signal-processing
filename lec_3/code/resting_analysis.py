import marimo

__generated_with = "0.11.17"
app = marimo.App()


@app.cell
def _():
    # the following line is necessary to have interactive plots. In other case, all plots will be inline (and not interactive)
    # '%matplotlib qt' command supported automatically in marimo
    from pathlib import Path
    import marimo as mo
    import numpy as np
    import matplotlib.pyplot as plt
    import mne
    from mne_bids import BIDSPath, read_raw_bids # you have to install mne_bids library with pip install mne-bids
    from mne.channels import make_standard_montage
    return (
        BIDSPath,
        Path,
        make_standard_montage,
        mne,
        mo,
        np,
        plt,
        read_raw_bids,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Loading the data (in BIDS format)

        The data are from the **"EEG Resting state EEG with closed eyes and open eyes in females from 60 to 80 years old database"** that can be found here [https://openneuro.org/datasets/ds005420/versions/1.0.0].

        More on the data can be found here:
        1. Appelhoff, S., Sanderson, M., Brooks, T., Vliet, M., Quentin, R., Holdgraf, C., Chaumon, M., Mikulan, E., Tavabi, K., Höchenberger, R., Welke, D., Brunner, C., Rockhill, A., Larson, E., Gramfort, A. and Jas, M. (2019). MNE-BIDS: Organizing electrophysiological data into the BIDS format and facilitating their analysis. Journal of Open Source Software 4: (1896). https://doi.org/10.21105/joss.01896

        2. Pernet, C. R., Appelhoff, S., Gorgolewski, K. J., Flandin, G., Phillips, C., Delorme, A., Oostenveld, R. (2019). EEG-BIDS, an extension to the brain imaging data structure for electroencephalography. Scientific Data, 6, 103. https://doi.org/10.1038/s41597-019-0104-81.

        In order to be able to load the EEG recordings of the subjects you need to define the correct path in the ```bids_root``` variable and the corresponding subject number in the ```subject``` variable.

        The result of the following block of code is two data arrays:
        1. ```raw_open``` with the raw EEG recordings while the subject has their eyes open
        2. ```raw_closed``` with the raw EEG recordings while the subject has their eyes closed
        """
    )
    return


@app.cell
def _(BIDSPath, Path, read_raw_bids):
    # Set path to BIDS dataset
    bids_root = Path().cwd() / "lec_3/Resting_State_Database"  # Path to the root of the BIDS dataset
    subject = "1"  # Replace with assigned subject number (without "sub-" prefix)

    bids_path_closed = BIDSPath(
        subject=subject,
        task="oc",  # 'oc' is the task code for eyes closed
        datatype="eeg",
        root=bids_root
    )

    # Load the data using BIDS paths
    raw_closed = read_raw_bids(bids_path_closed, verbose=True)

    # Then load the data into memory
    raw_closed.load_data()
    return bids_path_closed, bids_root, raw_closed, subject


@app.cell
def _(BIDSPath, bids_root, read_raw_bids, subject):
    bids_path_open = BIDSPath(
        subject=subject,
        task="oa",  # 'oa' is the task code for eyes open
        datatype="eeg",
        root=bids_root
    )

    raw_open = read_raw_bids(bids_path_open, verbose=True)

    raw_open.load_data()
    return bids_path_open, raw_open


@app.cell
def _(raw_open):
    # The BIDS reader automatically loads metadata like channel types, sampling frequency, etc.
    print(f"Sampling rate: {raw_open.info['sfreq']} Hz")
    print(f"Number of channels: {len(raw_open.ch_names)}")
    print(f"Channel names: {raw_open.ch_names}")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Now we will rename the channel names into something more typical following the 10-20 international system. In addition we will also rename the EOG channel (this will not be used in our analysis, but feel free to use it if you want to further process the data) into something more standard.""")
    return


@app.cell
def _(make_standard_montage, raw_closed, raw_open):
    # Create a mapping from your channel names to standard 10-20 names
    ch_mapping = {
        'EEG Fp1-A1A2': 'Fp1',
        'EEG Fp2-A1A2': 'Fp2',
        'EEG Fz-A1A2': 'Fz',
        'EEG F3-A1A2': 'F3',
        'EEG F4-A1A2': 'F4',
        'EEG F7-A1A2': 'F7',
        'EEG F8-A1A2': 'F8',
        'EEG Cz-A1A2': 'Cz',
        'EEG C3-A1A2': 'C3',
        'EEG C4-A1A2': 'C4',
        'EEG T3-A1A2': 'T7',  # or 'T7' in newer naming
        'EEG T4-A1A2': 'T8',  # or 'T8' in newer naming
        'EEG Pz-A1A2': 'Pz',
        'EEG P3-A1A2': 'P3',
        'EEG P4-A1A2': 'P4',
        'EEG T5-A1A2': 'P7',  # or 'P7' in newer naming
        'EEG T6-A1A2': 'P8',  # or 'P8' in newer naming
        'EEG O1-A1A2': 'O1',
        'EEG O2-A1A2': 'O2',
    }

    # First, explicitly set the channel type for the EOG channel
    raw_open.set_channel_types({'EEG LOC-ROC': 'eog'})
    raw_closed.set_channel_types({'EEG LOC-ROC': 'eog'})

    # Now rename just the EEG channels
    raw_open.rename_channels(ch_mapping)
    raw_closed.rename_channels(ch_mapping)

    # Set the montage with on_missing='ignore' to ignore the EOG channel
    montage = make_standard_montage('standard_1020')
    raw_open.set_montage(montage, on_missing='ignore')
    raw_closed.set_montage(montage, on_missing='ignore')

    # Plot to verify channel locations
    fig = raw_open.plot_sensors(show_names=True)
    return ch_mapping, fig, montage


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Now we can plot our two data structures, eyes open / eyes closed using interactive plots of MNE.""")
    return


@app.cell
def _(mo, plt, raw_closed):
    raw_closed.plot()
    # plt.gcf() gets the current figure
    mo.mpl.interactive(plt.gcf())

    return


@app.cell
def _(mo, plt, raw_open):
    raw_open.plot()
    mo.mpl.interactive(plt.gcf())
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Signal Preprocessing

        In this part we present the basic workflow of the signal pre-processing, that includes:
        1. Bandpass filtering keeping frequencies in the range of 1-40Hz (since we don't really care for higher frequencies in our analysis).\
        2. Notch filtering the mains frequency (50Hz) for Europe.

        Note that since we are particularly interested in analyzing the brain activity inside the alpha (8-13Hz) band, we could have band-pass filtered our recordings using lower value for the upper-frequency of the filter (i.e. 30Hz)
        """
    )
    return


@app.cell
def _(plt, raw_closed, raw_open):
    # Preprocessing
    # 1. Apply bandpass filter (1-40 Hz)
    raw_open_filtered = raw_open.copy().filter(l_freq=1, h_freq=40, 
                                             method='fir', phase='zero',
                                             fir_design='firwin')

    raw_closed_filtered = raw_closed.copy().filter(l_freq=1, h_freq=40, 
                                                 method='fir', phase='zero',
                                                 fir_design='firwin')

    # 2. Apply notch filter for line noise (50 Hz)
    raw_open_filtered.notch_filter(freqs=50, picks='eeg')
    raw_closed_filtered.notch_filter(freqs=50, picks='eeg')

    # 3. Identify bad channels
    # Option 1: Visual inspection
    raw_open_filtered.plot(duration=30, n_channels=20)  # Look for bad channels
    raw_closed_filtered.plot(duration=30, n_channels=20)

    plt.show()
    return raw_closed_filtered, raw_open_filtered


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Spectral Analysis

        In this part, we perform the spectral analysis of the signals. 

        1. First we plot the Power Spectral Density (PSD) of the signals in both conditions (eyes open/eyes closed) using ```fmin=0Hz``` to ```fmax=48Hz``` for each channel. Finally we plot the **average** PSD using all channels for each condition.
        2. We make a selection of EEG electrodes near the region of interest. Since we are interested for alpha activity which is located at the back of the brain (occipital), we don't want the signal to analyze and plot to be "contaminated" by brain activity not directly related to the eyes open/closed experiment. For this, we consider only the parietal and occipital electrodes and in particular the P3, P4, Pz, O1 and O2 and plot the PSD for each one of these.
        3. Then, we plot the topography of the power as this is distributed on the electrode space (topographic plot) in different basic frequency bands (delta - theta - alpha - beta - gamma). This is very informative, as it can reveal the area that significant PSD is located (we expect this area to be at the back of the head for the eyes closed condition).
        4. Finally we show each channel's PSD relative to its position on the head. This plot is interactive, as we can click on any PSD and have it zoomed in a separate window.
        """
    )
    return


@app.cell
def _(plt, raw_closed_filtered, raw_open_filtered):
    # First we estimate and then plot the mean PSD across all electrodes for both the conditions 
    spectrum_open = raw_open_filtered.compute_psd(fmin=0, fmax = 48)
    spectrum_open.plot(average=True, picks="data", exclude="bads", amplitude=False)

    spectrum_closed = raw_closed_filtered.compute_psd(fmin=0, fmax = 48)
    spectrum_closed.plot(average=True, picks="data", exclude="bads", amplitude=False)

    # Next we define a region of interest on the electrode space and plot the power of each electrode
    channel_sel = ['O1', 'O2', 'P3', 'P4', 'Pz']
    spectrum_open.plot(picks=channel_sel, exclude="bads", amplitude=False)

    spectrum_closed.plot(picks=channel_sel, exclude="bads", amplitude=False)

    # Topographic plots of the distribution of the power on the electrode space for both conditions.
    spectrum_open.plot_topomap()
    spectrum_closed.plot_topomap()

    # Finally we plot all channels' PSD plots in a way relative to their position on the head
    spectrum_open.plot_topo()
    spectrum_closed.plot_topo()


    plt.show()
    return channel_sel, spectrum_closed, spectrum_open


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        1. In this part, we continue estimating the PSD of the signals for both conditions, using specific parameters in the calculation of the PSD, using the Welch method. Then, we plot the PSD of the O1, O2 channels for both conditions on a plot that has highlighted the alpha frequency band (8-13Hz). These plors don't differentiate from the previous that rely solely on MNE's functions and methods, but with a little extra cost programmatically, we can shape the plots in the way we want (to be more informative).

        2. We calculate the mean power for each channel in both cases, but only **inside the alpha band** and we print the power ratio for each channel as

        3. Finally we plot alpha power topography for eyes open, eyes closed as well as the power_ratio for each channel
        """
    )
    return


@app.cell
def _(mne, np, plt, raw_closed_filtered, raw_open_filtered):
    (fmin, fmax) = (0.5, 45)
    (psd_open, freqs) = raw_open_filtered.compute_psd(method='welch', fmin=fmin, fmax=fmax).get_data(return_freqs=True)
    (psd_closed, freqs) = raw_closed_filtered.compute_psd(method='welch', fmin=fmin, fmax=fmax).get_data(return_freqs=True)
    picks = ['O1', 'O2']
    pick_idx = [raw_open_filtered.ch_names.index(ch) for ch in picks]
    plt.figure(figsize=(12, 8))
    for (idx, ch_idx) in enumerate(pick_idx):
        plt.semilogy(freqs, psd_open[ch_idx], label=f'{raw_open_filtered.ch_names[ch_idx]} (Open)', alpha=0.5, linestyle='--')
        plt.semilogy(freqs, psd_closed[ch_idx], label=f'{raw_closed_filtered.ch_names[ch_idx]} (Closed)', alpha=0.5)
    plt.axvspan(8, 13, color='yellow', alpha=0.2, label='Alpha Band (8-13 Hz)')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Power Spectral Density (µV²/Hz)')
    plt.title('Power Spectrum: Eyes Open vs Eyes Closed')
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.xlim(fmin, 30)
    plt.show()
    alpha_band = (8, 13)
    _alpha_idx = np.logical_and(freqs >= alpha_band[0], freqs <= alpha_band[1])
    alpha_power_open = np.mean(psd_open[:, _alpha_idx], axis=1)
    alpha_power_closed = np.mean(psd_closed[:, _alpha_idx], axis=1)
    alpha_ratio = alpha_power_closed / alpha_power_open
    print('Alpha Power (8-13 Hz):')
    for (i, ch_name) in enumerate(raw_open_filtered.ch_names[:-1]):
        print(f'{ch_name}:')
        print(f'  Eyes Open: {alpha_power_open[i]:.2f} µV²/Hz')
        print(f'  Eyes Closed: {alpha_power_closed[i]:.2f} µV²/Hz')
        print(f'  Ratio (Closed/Open): {alpha_ratio[i]:.2f}x')
    mne.viz.plot_topomap(alpha_power_open, raw_open_filtered.info, show=False, cmap='viridis')
    plt.title('Alpha Power - Eyes Open')
    mne.viz.plot_topomap(alpha_power_closed, raw_closed_filtered.info, show=False, cmap='viridis')
    plt.title('Alpha Power - Eyes Closed')
    mne.viz.plot_topomap(alpha_ratio, raw_open_filtered.info, show=False, cmap='viridis')
    plt.title('Alpha Power Ratio (Closed/Open)')
    plt.show()
    return (
        alpha_band,
        alpha_power_closed,
        alpha_power_open,
        alpha_ratio,
        ch_idx,
        ch_name,
        fmax,
        fmin,
        freqs,
        i,
        idx,
        pick_idx,
        picks,
        psd_closed,
        psd_open,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Time Frequency Analysis

        Here we perform further analysis of our signals using Time Frequency analysis techniques. Since the signal that we examine cannot be considered to be stationary, we want to break it in smaller parts and perform the analysis therein.
        For this purpose, we take only the electrode 'O1' (but you are more than welcome to analyze any other electrode you like) and we break it into smaller time windows with duration of 2sec (again feel free to play with the duration and see how the results change).
        Next we perform time-frequency analysis using Wavelet analysis (a common method of multiresolution analysis) inside a range of frequencies 4-30Hz.
        """
    )
    return


@app.cell
def _(mne, np, raw_closed_filtered, raw_open_filtered):
    from mne.time_frequency import tfr_morlet
    ch_name_1 = 'O1'
    ch_idx_1 = raw_open_filtered.ch_names.index(ch_name_1)
    epoch_duration = 2.0
    events_open = mne.make_fixed_length_events(raw_open_filtered, duration=epoch_duration)
    events_closed = mne.make_fixed_length_events(raw_closed_filtered, duration=epoch_duration)
    epochs_open = mne.Epochs(raw_open_filtered, events_open, tmin=0, tmax=epoch_duration, baseline=None, preload=True, reject=None)
    epochs_closed = mne.Epochs(raw_closed_filtered, events_closed, tmin=0, tmax=epoch_duration, baseline=None, preload=True, reject=None)
    freqs_1 = np.arange(4, 30, 1)
    n_cycles = freqs_1 / 2.0
    power_open = epochs_open.compute_tfr(method='morlet', freqs=freqs_1, n_cycles=n_cycles, average=True, return_itc=False, decim=3)
    power_closed = epochs_closed.compute_tfr(method='morlet', freqs=freqs_1, n_cycles=n_cycles, average=True, return_itc=False, decim=3)
    return (
        ch_idx_1,
        ch_name_1,
        epoch_duration,
        epochs_closed,
        epochs_open,
        events_closed,
        events_open,
        freqs_1,
        n_cycles,
        power_closed,
        power_open,
        tfr_morlet,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""After calculating the power of the signals in different frequencies and the different 2sec time intervals, we proceed and plot them. What we would like *ideally* to see is the power to be significantly larger in the alpha band for the eyes closed condition compared to the eyes open condition.""")
    return


@app.cell
def _(ch_idx_1, ch_name_1, plt, power_closed, power_open):
    power_open.plot([ch_idx_1], baseline=None, mode='mean', title=f'Time-Frequency Analysis: {ch_name_1} - Eyes Open', show=False)
    power_closed.plot([ch_idx_1], baseline=None, mode='mean', title=f'Time-Frequency Analysis: {ch_name_1} - Eyes Closed', show=False)
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        Next we want to show a comparative plot of the mean power of the channel we analyze within the frequency band of interest (alpha band 8-13Hz). 

        This plot is different from the previous, as it shows the **mean** power across all epochs (*remember the 2sec blocks we formed earlier?*) for the two different conditions (eyes open/closed). This results in very informative plots as *theoretically* we should see that the power of the occipital channel for the eyes closed condition should be **greater** than the power of the signal in the eyes open condition
        """
    )
    return


@app.cell
def _(ch_idx_1, ch_name_1, freqs_1, np, plt, power_closed, power_open):
    alpha_freqs = np.logical_and(freqs_1 >= 8, freqs_1 <= 13)
    _alpha_idx = np.where(alpha_freqs)[0]
    plt.figure(figsize=(12, 6))
    plt.plot(power_open.times, np.mean(power_open.data[ch_idx_1, _alpha_idx, :], axis=0), label='Eyes Open', color='blue')
    plt.plot(power_closed.times, np.mean(power_closed.data[ch_idx_1, _alpha_idx, :], axis=0), label='Eyes Closed', color='red')
    plt.xlabel('Time (s)')
    plt.ylabel('Alpha Power')
    plt.title(f'Alpha Power Over Time: {ch_name_1}')
    plt.legend()
    plt.grid(True)
    plt.show()
    return (alpha_freqs,)


if __name__ == "__main__":
    app.run()
