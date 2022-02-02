from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns

from .database import Database


class SingleCellAnalysis:
    
    def __init__(self, database: Database):
        self.database = database
        self.events_all_stim_paradigms = {'stimulation_string': list(), 
                                          'amplitude': list(), 
                                          'inter_event_interval': list(), 
                                          'sweep': list()}

    
    def run_analysis(self, global_cell_id: str, show: bool, save: bool):
        filepaths = self.database.list_all_column_values(global_cell_id = global_cell_id, column_name = 'filepath_detected_events')
        filepaths = [Path(elem) for elem in filepaths]
        
        for recording_filepath in filepaths:
            events_single_stim_paradigm = pd.read_csv(recording_filepath)
            columns = list(events_single_stim_paradigm.columns)
            for column_name in columns:
                if 'WaveN' in column_name:
                    sweep_column = column_name
                elif 'PeakT' in column_name:
                    event_time_column = column_name
                elif 'AmpY' in column_name:
                    event_amplitude_column = column_name
                else:
                    continue
        
            stim_string = self.database.cell_recordings_metadata.loc[self.database.cell_recordings_metadata['filepath_detected_events'] == recording_filepath.as_posix(), 'stimulation_string'].values[0]

            mean_inter_event_intervals = list()
            for sweep_id in events_single_stim_paradigm[sweep_column].unique():
                event_times = events_single_stim_paradigm.loc[events_single_stim_paradigm[sweep_column] == sweep_id, event_time_column]
                time_differences = event_times.rolling(2).apply(np.diff)
                mean_inter_event_intervals += list(time_differences.rolling(3, center=True).mean().values)

            self.events_all_stim_paradigms['stimulation_string'] += [stim_string] * events_single_stim_paradigm.shape[0]
            self.events_all_stim_paradigms['sweep'] += list(events_single_stim_paradigm[sweep_column].values)
            self.events_all_stim_paradigms['inter_event_interval'] += mean_inter_event_intervals
            self.events_all_stim_paradigms['amplitude'] += list(events_single_stim_paradigm[event_amplitude_column].values)

        df_all_events = pd.DataFrame(data=self.events_all_stim_paradigms)
        
        parameters_to_compare = ['amplitude', 'inter_event_interval']
        
        stim_paradigms = list(df_all_events['stimulation_string'].unique())
        stim_paradigms.remove('baseline')
        n_multiple_comparisons = len(stim_paradigms)

        fig = plt.figure(figsize=(15, 5*len(stim_paradigms)), facecolor='white')
        gs = fig.add_gridspec(len(stim_paradigms), 2)

        for stim_string in stim_paradigms:
            df_temp = df_all_events.loc[df_all_events['stimulation_string'].isin([stim_string, 'baseline'])]
            for parameter in parameters_to_compare:
                    data_baseline = df_temp.loc[df_temp['stimulation_string'] == 'baseline', parameter].values
                    data_stim_paradigm = df_temp.loc[df_temp['stimulation_string'] == stim_string, parameter].values
                    test_statistic, uncorrected_p_value = stats.ks_2samp(data_baseline, data_stim_paradigm)
                    corrected_p_value = uncorrected_p_value * n_multiple_comparisons
                    if corrected_p_value <= 0.001:
                        stars = '***'
                    elif corrected_p_value <= 0.01:
                        stars = '**'
                    elif corrected_p_value <= 0.05:
                        stars = '*'
                    else:
                        stars = 'n.s.'

                    fig.add_subplot(gs[stim_paradigms.index(stim_string), parameters_to_compare.index(parameter)])
                    sns.ecdfplot(data = df_temp, x = parameter, hue = 'stimulation_string', hue_order = ['baseline', stim_string])
                    plt.title(f'{stim_string}: {stars} ({round(corrected_p_value, 3)})')
        plt.suptitle(f'comparison of cell: {global_cell_id}')
        plt.tight_layout()
        plt.show()

        