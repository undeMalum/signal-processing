import marimo

__generated_with = "0.11.6"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# Definitions""")
    return


@app.cell
def _():
    # '%matplotlib qt' command supported automatically in marimo

    import mne
    import matplotlib.pyplot as plt
    return mne, plt


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# Functions""")
    return


@app.cell
def _(mne):
    def simple_mne_demo():
        """
        Simple demonstration of MNE-Python capabilities using sample data
        """
        print("Loading sample data...")
        # Download and load sample data
        sample_data_folder = mne.datasets.sample.data_path()
        sample_data_raw_file = (sample_data_folder / 'MEG' / 'sample' /
                               'sample_audvis_raw.fif')
        raw = mne.io.read_raw_fif(sample_data_raw_file)

        # Print some basic information about the data
        print("\nDataset Information:")
        print(f"Number of channels: {len(raw.ch_names)}")
        print(f"Sampling frequency: {raw.info['sfreq']} Hz")
        print(f"Recording length: {raw.times.max():.2f} seconds")

        # Plot first 10 seconds of data
        print("\nPlotting first 10 seconds of EEG channels...")
        # Pick only EEG channels for clarity
        raw_eeg = raw.pick_types(meg=False, eeg=True)



        # Plot the first 10 seconds
        # raw_eeg.plot(duration=10, n_channels=5, 
        #             scalings='auto', title='Sample EEG Data')

        # Create a Power Spectral Density plot
        print("\nCreating PSD plot...")
        # raw_eeg.plot_psd(fmax=50, average=True)
        # raw.compute_psd().plot(fmax=20,average=True)
        raw.compute_psd().plot(average=True)
    return (simple_mne_demo,)


@app.cell
def _(simple_mne_demo):
    simple_mne_demo()
    return


@app.cell
def _(mo):
    form = mo.ui.text(label="Welcome aboard").form()
    form
    return (form,)


@app.cell
def _(form):
    form.value
    return


@app.cell
def _(mo):
    import random

    # instead of random.randint, in your notebook you'd use the value of
    # an upstream UI element or other Python object
    dictionary = mo.ui.dictionary({str(i): mo.ui.text() for i in range(random.randint(1, 10))})
    dictionary
    return dictionary, random


@app.cell
def _(dictionary):
    dictionary.value
    return


@app.cell
def _(mo, random):
    # Create a state object that will store the index of the
    # clicked button
    get_state, set_state = mo.state(None)

    # Create an mo.ui.array of buttons - a regular Python list won't work.
    buttons = mo.ui.array(
        [
            mo.ui.button(
                label="button " + str(i), on_change=lambda v, i=i: set_state(i)
            )
            for i in range(random.randint(2, 5))
        ]
    )

    # Put the buttons array into the table
    table = mo.ui.table(
        {
            "Action": ["Action Name"] * len(buttons),
            "Trigger": list(buttons),
        }
    )
    table
    return buttons, get_state, set_state, table


@app.cell
def _(get_state):
    get_state()
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
