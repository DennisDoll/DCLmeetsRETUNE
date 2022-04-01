from pathlib import Path
from typing import List, Dict, Optional

import pandas as pd
from pandas import DataFrame
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns

from .database import Database


class BoxplotAnalysis:
    
    def __init__(self, database: Database, df_to_use: DataFrame, recording_type: str, plot_title: str):
        self.database = database
        self.df = df_to_use
        self.recording_type = recording_type
        self.plot_title = plot_title
        
    
    def run_analysis(self, group_column: str, group_id: str, show: bool, save: bool):
        self.group_column = group_column
        self.group_id = group_id
        data = self.df.loc[self.df[group_column] == group_id]
        self.plot_boxplots(data = data, show = show, save = save)
        
    
    def plot_boxplots(self, data: DataFrame, show: bool, save: bool):
        columns = list(data.columns)
        relevant_columns = columns[columns.index('pharmacology') + 1 : columns.index('filepath_main_excel_sheet')]
        n_cols = 5
        if len(relevant_columns) % n_cols == 0:
            n_rows = int(len(relevant_columns) / n_cols)
        else:
            n_rows = int(len(relevant_columns) / n_cols) + 1
            
        fig = plt.figure(figsize=(16, 5*n_rows), facecolor='white')
        gs = fig.add_gridspec(n_rows, n_cols)
        
        for column_name in relevant_columns:
            row_idx = int(relevant_columns.index(column_name) / n_cols)
            col_idx = relevant_columns.index(column_name) % n_cols
            ax = fig.add_subplot(gs[row_idx, col_idx])
                
            sns.boxplot(data = data, x = 'stimulation_string', y = column_name, showfliers = False, palette = 'viridis', ax=ax)
            sns.stripplot(data = data, x = 'stimulation_string', y = column_name, color = 'black', size=7, ax=ax)
            
            if relevant_columns.index(column_name) < len(relevant_columns) - n_cols:
                plt.xlabel('')
                plt.xticks([])
            else:
                plt.xlabel('')
                plt.xticks(rotation = 90)
        plt.suptitle(self.plot_title, y=1.0, fontsize=12)
        plt.tight_layout()

        if save:
            if len(list(data['global_cell_id'])) == 1:
                directory = self.database.subdirectories.single_cell_analyses.as_posix()
                filename = f'{list(data["global_cell_id"])[0]}_{self.recording_type}'
            else:
                directory = self.database.subdirectories.group_analyses.as_posix()
                filename = f'{self.group_id}_in_{self.group_column}_{self.recording_type}'
            plt.savefig(f'{directory}/{filename}.png', dpi=300)
        
        if show:    
            plt.show()
        else:
            plt.close()
            
        
        


class CDFAnalysis:
    
    def __init__(self, database: Database, df_to_use: DataFrame, recording_type: str, plot_title: Optional[str]=None):
        self.database = database
        self.df = df_to_use
        self.recording_type = recording_type
        self.plot_title = plot_title
        self.events_all_stim_paradigms = {'global_cell_id': list(),
                                          'stimulation_string': list(), 
                                          'amplitude': list(), 
                                          'inter_event_interval': list(), 
                                          'sweep': list()}
    
    def get_stars_string(self, p: float) -> str:
        if p <= 0.001:
            stars_string = '***'
        elif p <= 0.01:
            stars_string = '**'
        elif p <= 0.05:
            stars_string = '*'
        else:
            stars_string = 'n.s.'
        return stars_string 
    
    def get_data_for_cumulative_distributions(self, filepath: Path):
            events_single_stim_paradigm = pd.read_csv(filepath)
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
            
            stim_string = self.df.loc[self.df['filepath_detected_events'] == filepath.as_posix(), 'stimulation_string'].values[0]
            global_cell_id = self.df.loc[self.df['filepath_detected_events'] == filepath.as_posix(), 'global_cell_id'].values[0]

            mean_inter_event_intervals = list()
            for sweep_id in events_single_stim_paradigm[sweep_column].unique():
                event_times = events_single_stim_paradigm.loc[events_single_stim_paradigm[sweep_column] == sweep_id, event_time_column]
                time_differences = event_times.rolling(2).apply(np.diff)
                mean_inter_event_intervals += list(time_differences.rolling(3, center=True).mean().values)

            self.events_all_stim_paradigms['global_cell_id'] += [global_cell_id] * events_single_stim_paradigm.shape[0]
            self.events_all_stim_paradigms['stimulation_string'] += [stim_string] * events_single_stim_paradigm.shape[0]
            self.events_all_stim_paradigms['sweep'] += list(events_single_stim_paradigm[sweep_column].values)
            self.events_all_stim_paradigms['inter_event_interval'] += mean_inter_event_intervals
            self.events_all_stim_paradigms['amplitude'] += list(events_single_stim_paradigm[event_amplitude_column].values)
    
    def get_all_recording_filepaths(self):
        self.global_cell_ids = list(self.df.loc[self.df[self.group_column] == self.group_id, 'global_cell_id'].unique())
        filepaths = list()
        for global_cell_id in self.global_cell_ids:
            tmp_filepaths = self.database.list_all_column_values(global_cell_id = global_cell_id, recording_type = self.recording_type, column_name = 'filepath_detected_events')
            tmp_filepaths = [Path(elem) for elem in tmp_filepaths]
            filepaths = filepaths + tmp_filepaths
        return filepaths

    def run_analysis(self, group_column: str, group_id: str, show: bool, save: bool):
        self.group_column = group_column
        self.group_id = group_id
        filepaths = self.get_all_recording_filepaths()
        for recording_filepath in filepaths:
            self.get_data_for_cumulative_distributions(filepath = recording_filepath)
        df_all_events = pd.DataFrame(data=self.events_all_stim_paradigms)
        self.plot_cdfs(data = df_all_events,
                       columns_for_cdfs = ['amplitude', 'inter_event_interval'], 
                       show = show,
                       save = save)
    
    def plot_cdfs(self, data: pd.DataFrame, columns_for_cdfs: List, show: bool, save: bool):
        stim_paradigms = list(data['stimulation_string'].unique())
        stim_paradigms.remove('baseline')
        n_multiple_comparisons = len(stim_paradigms)

        fig = plt.figure(figsize=(12, 3*len(stim_paradigms)), facecolor='white')
        gs = fig.add_gridspec(len(stim_paradigms), 2)

        for stim_string in stim_paradigms:
            df_temp = data.loc[data['stimulation_string'].isin([stim_string, 'baseline'])]
            for column_name in columns_for_cdfs:
                    data_baseline = df_temp.loc[df_temp['stimulation_string'] == 'baseline', column_name].values
                    n_baseline_cells = df_temp.loc[df_temp['stimulation_string'] == 'baseline', 'global_cell_id'].unique().shape[0]
                    data_stim_paradigm = df_temp.loc[df_temp['stimulation_string'] == stim_string, column_name].values
                    n_stim_cells = df_temp.loc[df_temp['stimulation_string'] == stim_string, 'global_cell_id'].unique().shape[0]
                    test_statistic, uncorrected_p_value = stats.ks_2samp(data_baseline, data_stim_paradigm)
                    corrected_p_value = uncorrected_p_value * n_multiple_comparisons
                    stars_string = self.get_stars_string(p = corrected_p_value)

                    ax = fig.add_subplot(gs[stim_paradigms.index(stim_string), columns_for_cdfs.index(column_name)])
                    sns.ecdfplot(data = df_temp, x = column_name, palette='viridis',
                                 hue = 'stimulation_string', hue_order = ['baseline', stim_string])
                    plt.text(0.5, 0.5, f'{stars_string} ({round(corrected_p_value, 3)})', 
                             horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
                    ax.spines['right'].set_visible(False)
                    ax.spines['top'].set_visible(False)
                    if columns_for_cdfs.index(column_name) < len(columns_for_cdfs) - 1:
                        ax.get_legend().remove()
                    else:
                        ax.get_legend().remove()
                        if len(self.global_cell_ids) == 1:
                            plt.legend(labels=['baseline', stim_string], loc='center left', bbox_to_anchor=(1, 0.5))
                        else:
                            plt.legend(labels=[f'baseline ({n_baseline_cells})', 
                                               f'{stim_string} ({n_stim_cells})'], 
                                               loc='center left', bbox_to_anchor=(1, 0.5))
                    plt.title(f'{stim_string} - {column_name}', pad = 10, fontsize=10)
        if len(self.global_cell_ids) == 1:
            plt.suptitle(f'Single cell analyses of global_cell_id #{self.global_cell_ids[0]}', fontsize=12)
        else:
            plt.suptitle(self.plot_title, y=1.0, fontsize=12)
        plt.tight_layout()
        
        if save:
            if len(self.global_cell_ids) == 1:
                directory = self.database.subdirectories.single_cell_analyses.as_posix()
                filename = f'{self.global_cell_ids[0]}_{self.recording_type}'
            else:
                directory = self.database.subdirectories.group_analyses.as_posix()
                filename = f'{self.group_id}_in_{self.group_column}_{self.recording_type}'
            for parameter in columns_for_cdfs:
                filename = filename + f'_{parameter}'
            plt.savefig(f'{directory}/{filename}.png', dpi=300)
        
        if show:    
            plt.show()
        else:
            plt.close()