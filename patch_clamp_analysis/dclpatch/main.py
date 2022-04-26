from pathlib import Path
import os
import pandas as pd
from pandas import DataFrame

from .database import Database, listdir_nohidden
from .analysis import CDFAnalysis, MeanComparisonOfCDFs, BoxplotAnalysis
from typing import List, Dict, Optional


CRITICAL_COLUMNS = ['mouse_line', 'sex', 'brain_region', 'cell_type', 'recording_type', 'internal_solution', 'pharmacology']


class PatchProject:
    
    def __init__(self, root_dir: Path):
        self.database = Database(root_dir = root_dir)
        
    def save_database(self):
        self.database.save_all()
        
    def load_database(self):
        self.database.load_all()
    
    def add_cell_to_database(self, path_to_cell_recordings_dir: Path, overwrite: bool=False):
        self.database.add_new_cell_recording(cell_recordings_dir = path_to_cell_recordings_dir, overwrite = overwrite)
        
    def remove_cell_from_database(self):
        pass
    
    def add_all_cells_to_database(self, path_with_cell_recordings_as_subdirs: Path, overwrite: bool=False):
        subdirectories = list()
        elements = listdir_nohidden(path = path_with_cell_recordings_as_subdirs)
        for elem in elements:
            subdir = f'{path_with_cell_recordings_as_subdirs.as_posix()}/{elem}'
            if os.path.isdir(subdir):
                subdirectories.append(Path(subdir))
        for cell_recordings_dir in subdirectories:
            path_to_cell_recordings_dir = Path(cell_recordings_dir)
            try:
                self.add_cell_to_database(path_to_cell_recordings_dir = path_to_cell_recordings_dir, overwrite = overwrite)
            except:
                user_warning_line1 = 'Warning! When I tried to load the data of the following path:\n'
                user_warning_line2 = f'{path_to_cell_recordings_dir.as_posix()}\n'
                user_warning_line3 = 'For now, I skipped these data to continue. Try loading that dataset individually to see more details about the Error.\n'
                user_warning_line4 = '_________________________________________'
                user_warning = user_warning_line1 + user_warning_line2 + user_warning_line3 + user_warning_line4
                print(user_warning)
        
    def compare_on_single_cell_level(self, global_cell_id: str, analysis_type: str, recording_type: str, 
                                     show: bool=True, save: bool=False, export: bool=False) -> None:
        self.check_for_valid_input_to_plotting_methods(analysis_type = analysis_type, recording_type = recording_type)
        df_to_use = self.select_corresponding_dataframe(recording_type = recording_type)
        if analysis_type == 'CDF':
            if recording_type == 'current_clamp':
                raise ValueError('CDF analysis is only possible for IPSPs or EPSPs.')
            single_cell_analysis = CDFAnalysis(database = self.database, df_to_use = df_to_use, recording_type = recording_type)
        elif analysis_type == 'Boxplot':
            raise ValueError('Boxplot analysis is only possible in the "compare_within_group()" method and not on single cell level.')
        single_cell_analysis.run_analysis(group_column = 'global_cell_id', group_id = global_cell_id, show = show, save = save)
        if export:
            filename = f'{analysis_type}_single_cell_analysis_{recording_type}_{global_cell_id}.xlsx'
            filepath = self.database.subdirectories.exported_excel_sheets.joinpath(filename)
            df_to_use.to_excel(filepath)
    
    
    def compare_within_group(self, group_column: str, group_id: str, analysis_type: str, recording_type: str, include: Optional[Dict]=None, 
                             exclude: Optional[Dict]=None, show: bool=True, save: bool=False, export: bool=False):
        self.check_for_valid_input_to_plotting_methods(analysis_type = analysis_type, recording_type = recording_type)
        df_to_use = self.select_corresponding_dataframe(recording_type = recording_type)
        if type(include) == dict:
            df_to_use = self.apply_inclusion_criteria(df = df_to_use, inclusion_criteria = include)
        if type(exclude) == dict:
            df_to_use = self.apply_exclusion_criteria(df = df_to_use, exclusion_criteria = exclude)
        df_to_use = df_to_use.loc[df_to_use[group_column] == group_id]
        self.check_for_multiple_groups_in_critical_columns(df = df_to_use)
        plot_title = self.create_plot_title(group_column = group_column, group_id = group_id, recording_type = recording_type, include = include, exclude = exclude)
        if 'CDF' in analysis_type: #CDF@50
            if recording_type == 'current_clamp':
                raise ValueError('CDF analysis is only possible for IPSPs or EPSPs.')
            if '@' in analysis_type:
                percentile = int(analysis_type[analysis_type.find('@') + 1:])
                group_analysis = MeanComparisonOfCDFs(database = self.database, df_to_use = df_to_use, recording_type = recording_type, plot_title = plot_title)
                group_analysis.run_analysis(group_column = group_column, group_id = group_id, percentile = percentile, show = show, save = save)
                if export:
                    percentile_dataset = group_analysis.get_data_for_export()
                    df_to_use = pd.DataFrame(data = percentile_dataset)
            elif analysis_type == 'CDF':
                group_analysis = CDFAnalysis(database = self.database, df_to_use = df_to_use, recording_type = recording_type, plot_title = plot_title)
                group_analysis.run_analysis(group_column = group_column, group_id = group_id, show = show, save = save)
        elif analysis_type == 'Boxplot':
            group_analysis = BoxplotAnalysis(database = self.database, df_to_use = df_to_use, recording_type = recording_type, plot_title = plot_title)
            group_analysis.run_analysis(group_column = group_column, group_id = group_id, show = show, save = save)
        if export:
            filename = f'{analysis_type}_group_analysis_{recording_type}_{group_column}_{group_id}.xlsx'
            filepath = self.database.subdirectories.exported_excel_sheets.joinpath(filename)
            #df_to_use.to_excel(filepath)
            return percentile_dataset
    
    
    def compare_between_groups(self):
        pass
    
    
    def create_plot_title(self, group_column: str, group_id: str, recording_type: str, include: Optional[Dict], exclude: Optional[Dict]) -> str:
        if (include == None) & (exclude == None):
            extra_criteria = 'no further criteria'
        if (include != None) & (exclude == None):
            extra_criteria = 'additional inclusion criteria'
        elif (include == None) & (exclude != None):
            extra_criteria = 'additional exclusion criteria'
        else:
            extra_criteria = 'addition in- and exclusion criteria' 
        return f'Group analyses of all {recording_type} recordings of cells matching:\n{group_id} in {group_column} and {extra_criteria}'
    
    
    def check_for_multiple_groups_in_critical_columns(self, df: DataFrame) -> None:
        for critical_column in CRITICAL_COLUMNS:
            if len(list(df[critical_column].unique())) > 1:
                print(f'Warning: Please be aware that I detected inconsistent values for the column "{critical_column}":\n')
                print(f'--> {list(df[critical_column].unique())}\n')
                print('You might consider adding all but one to the exclusion criteria!')
                print('_______________________________________')         
        
    
    def select_corresponding_dataframe(self, recording_type: str) -> DataFrame:
        if recording_type == 'current_clamp':
            df = self.database.current_clamp_recordings
        elif recording_type == 'IPSPs':
            df = self.database.inhibitory_voltage_clamp_recordings
        elif recording_type == 'EPSPs':
            df = self.database.excitatory_voltage_clamp_recordings
        return df.copy()
    
    
    def check_for_valid_input_to_plotting_methods(self, analysis_type: str, recording_type: str) -> None:
        if analysis_type not in ['CDF', 'Boxplot']:
            if 'CDF@' not in analysis_type:
                raise ValueError(f'The "analysis_type" attribute you provide has to be one of the following: "CDF" or "Boxplot".')
            else:
                try:
                    probability_to_compare = int(analysis_type[analysis_type.find('@') + 1:])
                except ValueError:
                    error_message_line0 = 'The probability you tried to pass is not in a valid format.\n'
                    error_message_line1 = 'For instance, if you want to compare the groups at a probability of 50, use: "analysis_type" = "CDF@50"'
                    raise ValueError(error_message_line0 + error_message_line1)
        if recording_type not in ['current_clamp', 'IPSPs', 'EPSPs']:
            raise ValueError(f'The "recording_type" attribute you provide has to be one of the following: "current_clamp", "IPSPs", or "EPSPs".')
        
    
    def apply_inclusion_criteria(self, df: DataFrame, inclusion_criteria: Dict) -> DataFrame:
        chunks_of_df_to_include = []
        for key, value in inclusion_criteria.items():
            if key not in df.columns:
                raise ValueError(f'{key} has to be in columns of the DataFrame. Please find a list of all possible columns: {list(df.columns)}')
            else:
                if type(value) == str:
                    chunks_of_df_to_include.append(df.loc[df[key] == value].copy())
                else:
                    raise ValueError(f'{value} is not a valid input! Only strings (e.g. "vlPAG") are possible, simply add a column multiple times with different values if needed.')
        df_after_inclusions = pd.concat(chunks_of_df_to_include, ignore_index = True)
        df_after_inclusions = df_after_inclusions.drop_duplicates(ignore_index = True)
        return df_after_inclusions

    
    def apply_exclusion_criteria(self, df: DataFrame, exclusion_criteria: Dict) -> DataFrame:
        df_after_exclusions = df.copy()
        for key, value in exclusion_criteria.items():
            if key not in df.columns:
                raise ValueError(f'{key} has to be in columns of the DataFrame. Please find a list of all possible columns: {list(df.columns)}')
            else:
                if type(value) == str:
                    df_after_exclusions = df_after_exclusions.loc[df_after_exclusions[key] != value]
                else:
                    raise ValueError(f'{value} is not a valid input! Only strings (e.g. "vlPAG") are possible, simply add a column multiple times with different values if needed.')
        return df_after_exclusions