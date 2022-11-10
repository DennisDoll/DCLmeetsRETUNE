from typing import List, Dict, Optional, Union, Tuple
from abc import ABC, abstractmethod

from pathlib import Path
from pandas import DataFrame
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pingouin as pg
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
            
        
class BaseCDFAnalysis(ABC):

    def __init__(self, database: Database, df_to_use: DataFrame, recording_type: str, plot_title: Optional[str]=None) -> None:
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


    def get_data_for_cumulative_distributions(self, filepath: Path) -> None:
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

            mean_inter_event_intervals = []
            for sweep_id in events_single_stim_paradigm[sweep_column].unique():
                event_times = events_single_stim_paradigm.loc[events_single_stim_paradigm[sweep_column] == sweep_id, event_time_column]
                time_differences = event_times.rolling(2).apply(np.diff)
                mean_inter_event_intervals += list(time_differences.rolling(3, center=True).mean().values)

            self.events_all_stim_paradigms['global_cell_id'] += [global_cell_id] * events_single_stim_paradigm.shape[0]
            self.events_all_stim_paradigms['stimulation_string'] += [stim_string] * events_single_stim_paradigm.shape[0]
            self.events_all_stim_paradigms['sweep'] += list(events_single_stim_paradigm[sweep_column].values)
            self.events_all_stim_paradigms['inter_event_interval'] += mean_inter_event_intervals
            self.events_all_stim_paradigms['amplitude'] += list(events_single_stim_paradigm[event_amplitude_column].values)

    """
    def get_all_recording_filepaths(self):
        # kept for backwards compatibility
        self.global_cell_ids = list(self.df.loc[self.df[self.group_column] == self.group_id, 'global_cell_id'].unique())
        filepaths = list()
        for global_cell_id in self.global_cell_ids:
            tmp_filepaths = self.database.list_all_column_values(global_cell_id = global_cell_id, recording_type = self.recording_type, column_name = 'filepath_detected_events')
            tmp_filepaths = [Path(elem) for elem in tmp_filepaths]
            filepaths = filepaths + tmp_filepaths
        return filepaths
    """
    
    def get_all_recording_filepaths_from_df(self):
        self.global_cell_ids = list(self.df['global_cell_id'].values)
        filepaths_as_strings = self.df['filepath_detected_events'].values
        filepaths = [Path(elem) for elem in filepaths_as_strings]
        return filepaths
    
    @abstractmethod
    def run_analysis(self, group_column: str, group_id: str, show: bool, save: bool) -> None:
        pass



class CDFAnalysis(BaseCDFAnalysis):


    def run_analysis(self, group_column: str, group_id: str, show: bool, save: bool):
        self.group_column = group_column
        self.group_id = group_id
        #filepaths = self.get_all_recording_filepaths()
        filepaths = self.get_all_recording_filepaths_from_df()
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



class MeanComparisonOfCDFs(BaseCDFAnalysis):
    
    @property
    def measurements(self):
        return ['amplitude', 'inter_event_interval']
    
    def run_analysis(self, group_column: str, group_id: str, percentile: int, show: bool, save: bool):
        self.group_column = group_column
        self.group_id = group_id
        self.percentile = percentile
        #filepaths = self.get_all_recording_filepaths()
        filepaths = self.get_all_recording_filepaths_from_df()
        for recording_filepath in filepaths:
            self.get_data_for_cumulative_distributions(filepath = recording_filepath)
        df_all_events = pd.DataFrame(data=self.events_all_stim_paradigms)
        self.percentile_data_per_stim_string = self.get_data_of_specific_percentile(df_all_events = df_all_events, percentile = percentile)
        self.plot_percentile_data(percentile_datasets = self.percentile_data_per_stim_string, show = show, save = save)
        
    
    def get_data_for_export(self) -> Tuple[List[pd.DataFrame], List[str]]:
        dfs, tab_names = [], []
        d = self.percentile_data_per_stim_string.copy()
        for stimulation_paradigm in d.keys():
            for analyzed_feature, data_per_condition in d[stimulation_paradigm].items():
                dcl_stats_n_plots_export_scheme = {'data': [], 'group_id': []}
                for condition, data_values in data_per_condition.items():
                    if condition == 'global_cell_id':
                        continue
                    else:
                        dcl_stats_n_plots_export_scheme['data'] += data_values
                        dcl_stats_n_plots_export_scheme['group_id'] += [condition] * len(data_values)
                dfs.append(pd.DataFrame(data=dcl_stats_n_plots_export_scheme))
                if analyzed_feature == 'inter_event_interval':
                    feature_name_for_tab = 'IEI'
                else:
                    feature_name_for_tab = analyzed_feature
                tab_names.append(f'{feature_name_for_tab}_at_{stimulation_paradigm}')        
        return dfs, tab_names
    
    
    def get_data_of_specific_percentile(self, df_all_events: pd.DataFrame, percentile: Union[int, str]) -> Dict:
        percentile_data_per_stim_string = dict()
        stim_paradigms = list(df_all_events['stimulation_string'].unique())
        stim_paradigms.remove('baseline')
        for stim_string in stim_paradigms:
            percentile_data_per_stim_string[stim_string] = dict()
            for measurement in self.measurements:
                percentile_data_per_stim_string[stim_string][measurement] = {'baseline': list(),
                                                                             'stimulated': list(), 
                                                                             'global_cell_id': list()}
                for global_cell_id in df_all_events.loc[df_all_events['stimulation_string'] == stim_string, 'global_cell_id'].unique():
                    if type(percentile) == int:
                        baseline = np.nanpercentile(df_all_events.loc[(df_all_events['stimulation_string'] == 'baseline') &
                                                                      (df_all_events['global_cell_id'] == global_cell_id), measurement].values,
                                                    percentile)
                        stimulated = np.nanpercentile(df_all_events.loc[(df_all_events['stimulation_string'] == stim_string) &
                                                                        (df_all_events['global_cell_id'] == global_cell_id), measurement].values,
                                                      percentile)
                    else: # percentile == mean:
                        baseline = np.nanmean(df_all_events.loc[(df_all_events['stimulation_string'] == 'baseline') &
                                                             (df_all_events['global_cell_id'] == global_cell_id), measurement].values)
                        stimulated = np.nanmean(df_all_events.loc[(df_all_events['stimulation_string'] == stim_string) &
                                                               (df_all_events['global_cell_id'] == global_cell_id), measurement].values)
                    percentile_data_per_stim_string[stim_string][measurement]['baseline'] += [baseline]
                    percentile_data_per_stim_string[stim_string][measurement]['stimulated'] += [stimulated]
                    percentile_data_per_stim_string[stim_string][measurement]['global_cell_id'] += [global_cell_id]
        return percentile_data_per_stim_string
    
    
    def plot_percentile_data(self, percentile_datasets: Dict, show: bool, save: bool, jitter: float=0.1) -> None:
        stimulation_strings = percentile_datasets.keys()
        n_rows = len(stimulation_strings)
        n_columns = len(self.measurements)
        
        fig = plt.figure(figsize=(6*n_columns, 3*n_rows), facecolor='white')
        gs = fig.add_gridspec(n_rows, n_columns)
        
        for stim_string, row_index in zip(stimulation_strings, range(n_rows)):
            for measurement, column_index in zip(self.measurements, range(n_columns)):
                df_temp = pd.DataFrame(data = percentile_datasets[stim_string][measurement])
                stats = pg.wilcoxon(x=df_temp['baseline'], y=df_temp['stimulated'])
                pval = stats['p-val']['Wilcoxon']
                stars_string = self.get_stars_string(p = pval)
                n_cells = df_temp.shape[0]
                df_x_jitter = pd.DataFrame(np.random.normal(loc=0, scale=jitter, size=df_temp.values.shape), columns=df_temp.columns)
                df_x_jitter += np.arange(len(df_temp.columns))
                colormixer = plt.cm.viridis(np.linspace(0, 1, len(df_temp.columns) + 2))
                colormixer = colormixer[1 : -1]
                ax = fig.add_subplot(gs[row_index, column_index])
                for column, color in zip(df_temp.columns, colormixer):
                    ax.plot(df_x_jitter[column], df_temp[column], 'o', color=color, alpha=.60, zorder=1, ms=8, mew=1)
                ax.set_xticks(range(2))
                ax.set_xticklabels(['baseline', 'stimulated'])
                ax.set_xlim(-0.5,2-0.5)
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)               
                for idx in df_temp.index:
                    ax.plot(df_x_jitter.loc[idx,['baseline','stimulated']], df_temp.loc[idx,['baseline','stimulated']], color = 'grey', linewidth = 0.5, linestyle = '--', zorder=-1)
                plt.text(0.5, 0.85, f'{stars_string} ({round(pval, 3)})', 
                         horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
                plt.title(f'{stim_string} - {measurement}', pad = 10, fontsize=10)
        if self.percentile == 'mean':
            plt.suptitle(f'{self.plot_title} at {self.percentile} values' , y=1.0, fontsize=12)
        else:
            plt.suptitle(f'{self.plot_title} at {self.percentile}th percentile' , y=1.0, fontsize=12)
        plt.tight_layout()

        if save:
            directory = self.database.subdirectories.group_analyses.as_posix()
            if self.percentile == 'mean':
                filename = f'{self.group_id}_in_{self.group_column}_{self.recording_type}_at_{self.percentile}'
            else:
                filename = f'{self.group_id}_in_{self.group_column}_{self.recording_type}_at_{self.percentile}th_percentile'
            plt.savefig(f'{directory}/{filename}.png', dpi=300)
        
        if show:    
            plt.show()
        else:
            plt.close()
        