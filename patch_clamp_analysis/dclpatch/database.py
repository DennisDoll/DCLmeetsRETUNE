from pathlib import Path
from abc import ABC, abstractmethod
from typing import List, Dict
from pandas import DataFrame
import pickle
import os
import datetime

import pandas as pd
import numpy as np

SORTED_COLUMNS = ['global_cell_id', 'date', 'session_cell_id', 'mouse_line', 'sex', 'age', 'brain_region', 'cell_type', 'recording_type',
                  'internal_solution', 'stimulation_in_percent', 'pharmacology',
                  'stimulation_string', 'stimulation_frequency-Hz', 'stimulation_duration-ms', 'filepath_main_excel_sheet', 'filepath_detected_events']

#'Recordings and cell properties' for reasons of backwards compatibility    
VOLTAGE_CLAMP_IDENTIFIERS = ['Voltage-clamp', 'voltage-clamp', 'Voltage_clamp', 'voltage_clamp', 'Recordings and cell properties'] 
CURRENT_CLAMP_IDENTIFIERS = ['Current-clamp', 'current-clamp', 'Current_clamp', 'current_clamp']


def listdir_nohidden(path: Path) -> List:
    return [f for f in os.listdir(path.as_posix()) if f.startswith('.') == False]


class Database:
    
    """
    The database is supposed to hold all information about all recorded cells that were added to it.
    These information exclude the raw data (only contain filepaths to raw data),
    but include general metadata like celltype, stimulus type, brain region, pharamacological treatment.
    Potentially, some intermediate data could be added (like how the cell reacted upon stimulation, see dashboard).
    Overall, the database should allow the selection of cells based on mix-and-match criteria 
    for further (statistical) analyses and plotting.
    """
    
    def __init__(self, root_dir: Path):
        self.root_dir = root_dir
        self.subdirectories = Subdirectories(root_dir = root_dir)
        self.global_cell_id = 0
        self.cells_in_database = list()


    def list_all_column_values(self, global_cell_id: str, recording_type: str, column_name: str) -> List:
        if recording_type not in ['current_clamp', 'IPSPs', 'EPSPs']:
            raise ValueError('The attribute "recording_type" has to be one of the following: "current_clamp", "IPSPs", or "EPSPs"')
        if recording_type == 'current_clamp':
            attribute = 'current_clamp_recordings'
        elif recording_type == 'IPSPs':
            attribute = 'inhibitory_voltage_clamp_recordings'
        elif recording_type == 'EPSPs':
            attribute = 'excitatory_voltage_clamp_recordings'
            
        if hasattr(self, attribute):
            df_to_use = getattr(self, attribute)
            if global_cell_id in df_to_use['global_cell_id'].unique():
                values = list(df_to_use.loc[df_to_use['global_cell_id'] == global_cell_id, column_name].values)
                return values
            else:
                raise ValueError(f'{global_cell_id} is not in the database yet!')
        else:
            print('No entries in the database yet')
                                 
    
    def add_new_cell_recording(self, cell_recordings_dir: Path, overwrite: bool):
        increase_global_cell_id_count = True
        if overwrite:
            if cell_recordings_dir.name in self.cells_in_database:
                date = cell_recordings_dir.name[:10]
                session_cell_id = cell_recordings_dir.name[cell_recordings_dir.name.rfind('_')+1:]
                global_cell_id = self.cell_recordings_metadata.loc[(self.cell_recordings_metadata['date'] == date) &
                                                                   (self.cell_recordings_metadata['session_cell_id'] == session_cell_id),
                                                                   'global_cell_id'].iloc[0]
                increase_global_cell_id_count = False
                self.cell_recordings_metadata.drop(self.cell_recordings_metadata.loc[self.cell_recordings_metadata['global_cell_id'] == global_cell_id].index, inplace = True)
                self.cells_in_database.remove(cell_recordings_dir.name)
            else:
                global_cell_id = self.global_cell_id
        else:
            global_cell_id = self.global_cell_id

        if cell_recordings_dir.name not in self.cells_in_database:
            new_recording = CellRecording(cell_recordings_dir = cell_recordings_dir, global_cell_id = global_cell_id)
            recording_overview = new_recording.create_recordings_overview()
            if hasattr(self, 'cell_recordings_metadata') == False:
                for column_name in list(recording_overview.columns):
                    if column_name not in SORTED_COLUMNS:
                        print(f'Warning: {column_name} not defined in SORTED_COLUMNS.')
                self.cell_recordings_metadata = recording_overview[SORTED_COLUMNS]
            else:
                for column_name in list(recording_overview.columns):
                    if column_name not in list(self.cell_recordings_metadata.columns):
                        print(f'Warning: {column_name} represents a new column. Consider updating of all previously added cell recordings.')
                self.cell_recordings_metadata = pd.concat([self.cell_recordings_metadata, recording_overview], ignore_index = True)
            self.cell_recordings_metadata.sort_values(by=['global_cell_id'], inplace = True, ignore_index = True)
            
            for recording_type_key, dataframe in new_recording.recording_type_specific_datasets.items():
                if type(dataframe) == DataFrame:
                    if hasattr(self, recording_type_key) == False:
                        setattr(self, recording_type_key, dataframe)
                    else:
                        # check if all columns are matching?
                        current_df = getattr(self, recording_type_key)
                        updated_df = pd.concat([current_df, dataframe], ignore_index = True)
                        updated_df.sort_values(by = ['global_cell_id'], inplace = True, ignore_index = True)
                        setattr(self, recording_type_key, updated_df)
                        

            self.cells_in_database.append(cell_recordings_dir) #simply add .name to get cell_id (probably requires adaptations at multiple other sites too) 
            if increase_global_cell_id_count:
                self.global_cell_id += 1
        
            # Trigger update of mix-and-match categories?
        else:
            note_line1 = f'Note: The recordings in "{cell_recordings_dir.as_posix()}" were already added to the database.\n'
            note_line2 =  '      Consider passing "overwrite = True" when calling this function in order to update the information.'
            print(note_line1 + note_line2)
            
            
    def save_all(self):
        self.save_all_attributes_as_pickle()
        self.save_cell_recordings_metadata_as_csv()

        
    def save_all_attributes_as_pickle(self):
        project_summary = self.__dict__.copy()
        filepath = f'{self.subdirectories.project_summary.as_posix()}/{datetime.datetime.now().strftime("%Y_%m_%d")}_dclpatch_project_summary.p'        
        with open(filepath, 'wb') as io:
            pickle.dump(project_summary, io)

    def save_cell_recordings_metadata_as_csv(self):
        filepath = f'{self.subdirectories.project_summary.as_posix()}/{datetime.datetime.now().strftime("%Y_%m_%d")}_dclpatch_database_overview.csv'  
        self.cell_recordings_metadata.to_csv(filepath)
    
    
    def load_all(self):
        result_files = [fname for fname in listdir_nohidden(self.subdirectories.project_summary) if fname.endswith('dclpatch_project_summary.p')]
        result_files.sort(reverse = True)
        with open(f'{self.subdirectories.project_summary.as_posix()}/{result_files[0]}', 'rb') as io:
            project_summary = pickle.load(io)

        for key, value in project_summary.items():
            if hasattr(self, key) == False:
                setattr(self, key, value)
        
    


class CellRecording:
    
    def __init__(self, cell_recordings_dir: Path, global_cell_id: int) -> None:
        self.cell_recordings_dir = cell_recordings_dir
        self.global_cell_id = str(global_cell_id).zfill(4)
        self.recording_type_specific_datasets = {'excitatory_voltage_clamp_recordings': None,
                                                 'inhibitory_voltage_clamp_recordings': None,
                                                 'current_clamp_recordings': None}
        self.extract_all_data_from_directory()
        
        
    def extract_all_data_from_directory(self) -> None:
        self.general_metadata_df = self.create_general_metadata_df()
        recording_type_data_extractors = self.prepare_all_recording_type_data_extractors()
        for data_extractor in recording_type_data_extractors:
            extracted_df = data_extractor.create_recording_specific_dataframe(cell_recording_obj = self)
            self.recording_type_specific_datasets[data_extractor.recording_type_key] = extracted_df
                
       
    def create_general_metadata_df(self) -> DataFrame:
        metadata_df = pd.read_excel(self.cell_recordings_dir.joinpath(f'{self.cell_recordings_dir.name}.xlsx'), sheet_name = 'General information')
        general_metadata = {'date': self.cell_recordings_dir.name[:10],
                            'session_cell_id': self.cell_recordings_dir.name[self.cell_recordings_dir.name.rfind('_') + 1:],
                            'mouse_line': metadata_df['Animal line'][0],
                            'brain_region': metadata_df['Region'][0],
                            'cell_type': metadata_df['Type'][0],
                            'sex': metadata_df['Sex'][0],
                            'age': metadata_df['Age (days)'][0],
                            'stimulation_in_percent': metadata_df['Stimulation (%)'][0],
                            'internal_solution': metadata_df['Internal used'][0]}
        return pd.DataFrame(general_metadata, index=[0])
    
    
    def prepare_all_recording_type_data_extractors(self) -> List: #returns a list of DataExtractor objects
        main_excel_sheet = pd.read_excel(self.cell_recordings_dir.joinpath(f'{self.cell_recordings_dir.name}.xlsx'), sheet_name = None)
        recording_type_data_extractors = list()
        for tab_name in main_excel_sheet.keys():
            if tab_name in VOLTAGE_CLAMP_IDENTIFIERS:
                recording_type_data_extractors.append(ExtractVoltageClampData(df = main_excel_sheet[tab_name]))
            elif tab_name in CURRENT_CLAMP_IDENTIFIERS:
                recording_type_data_extractors.append(ExtractCurrentClampData(df = main_excel_sheet[tab_name]))
            elif tab_name == 'General information':
                continue
            else:
                raise ValueError(f'"{tab_name}" is not a valid identifier of recording type for {self.cell_recordings_dir.as_posix()}')
        return recording_type_data_extractors
                
    
    def create_recordings_overview(self) -> DataFrame:
        recording_specific_dataframes = []
        for recording_type, value in self.recording_type_specific_datasets.items():
            if type(value) == DataFrame:
                relevant_columns1 = ['global_cell_id', 'recording_type','pharmacology', 'filepath_main_excel_sheet', 'filepath_detected_events']
                relevant_columns2 = ['stimulation_string', 'stimulation_frequency-Hz', 'stimulation_duration-ms']
                df_tmp = value[relevant_columns1 + relevant_columns2].copy()
                recording_specific_dataframes.append(df_tmp)
        overview_df = pd.concat(recording_specific_dataframes)
        overview_df = overview_df.reset_index(drop = True)
        shape_adapted_general_metadata_df = pd.concat([self.general_metadata_df]*overview_df.shape[0], ignore_index = True)
        overview_df = pd.concat([shape_adapted_general_metadata_df, overview_df], axis=1)
        return overview_df


    
class DataExtractor(ABC):
    
    def __init__(self, df: DataFrame) -> None:
        self.input_df = df

    
    @property
    @abstractmethod
    def recording_method(self) -> str:
        # e.g. 'voltage_clamp' or 'current_clamp'
        pass

    
    def create_recording_specific_dataframe(self, cell_recording_obj: CellRecording) -> DataFrame:
        input_df = self.reshape_df(df = self.input_df.copy())
        self.check_for_baseline_duplicates(df = input_df, cell_recording_obj = cell_recording_obj)
        columns = list(input_df.columns)
        relevant_data_columns = columns[columns.index('Pharmacology') + 1 :]
        all_final_single_stimulation_dfs = []
        for row_index in range(input_df.shape[0]):
            single_stimulation_df = input_df.iloc[[row_index]].copy()
            remaining_data_df = single_stimulation_df[relevant_data_columns].copy()
            remaining_data_df = remaining_data_df.reset_index(drop=True)
            sorted_single_stimulation_df = self.create_sorted_single_stimulation_df(single_stimulation_df = single_stimulation_df)
            all_final_single_stimulation_dfs.append(pd.concat([sorted_single_stimulation_df, remaining_data_df], axis=1))
        recording_specific_dataframe = pd.concat(all_final_single_stimulation_dfs)
        final_df = pd.concat([cell_recording_obj.general_metadata_df, recording_specific_dataframe], axis=1)
        final_df.insert(0, 'global_cell_id', [cell_recording_obj.global_cell_id]*final_df.shape[0])
        final_df = self.add_corresponding_filepaths(df_without_filepaths = final_df, cell_recordings_dir = cell_recording_obj.cell_recordings_dir)
        return final_df


    def reshape_df(self, df: DataFrame) -> DataFrame:
        df.columns = df.iloc[1]
        df = df.iloc[2:, 1:]
        df = df.loc[df['Type of protocol'].notnull()]
        return df
    
    
    def check_for_baseline_duplicates(self, df: DataFrame, cell_recording_obj: CellRecording) -> None:
        listed_protocols = list(df['Type of protocol'].values)
        if listed_protocols.count('Bsl') > 1:
            filepath = cell_recording_obj.cell_recordings_dir.joinpath(f'{cell_recording_obj.cell_recordings_dir.name}.xlsx')
            error_message_line1 = 'There are multiple "Bsl" identifiers in the following excel sheet:\n'
            error_message_line2 = f'{filepath}\n'
            error_message_line3 = 'If this is a control baseline, please use "Bsl2" as identifier.\n'
            error_message_line4 = 'If these baselines come from IPSP and EPSP recordings, please split these data in two different tabs.'
            error_message = error_message_line1 + error_message_line2 + error_message_line3 + error_message_line4
            raise ValueError(error_message)
    
    
    def create_sorted_single_stimulation_df(self, single_stimulation_df: DataFrame) -> DataFrame:
            if type(single_stimulation_df['Time since cutting'].iloc[0]) == datetime.time:
                time_since_cutting = single_stimulation_df['Time since cutting'].iloc[0].hour * 60 + single_stimulation_df['Time since cutting'].iloc[0].minute
            else:
                time_since_cutting = 'not_available'
                
            if self.recording_method == 'voltage_clamp':
                if single_stimulation_df['Vh (mV)'].iloc[0] == -70:
                    recording_type = 'excitatory_voltage_clamp_recordings'
                elif single_stimulation_df['Vh (mV)'].iloc[0] == 0:
                    recording_type = 'inhibitory_voltage_clamp_recordings'
                else:
                    raise ValueError(f'{single_stimulation_df["Vh (mV)"].iloc[0]} is not a valid input for holding potential!')
            else:
                recording_type = 'current_clamp_recordings'
            self.recording_type_key = recording_type

            if pd.isna(single_stimulation_df['Pharmacology'].iloc[0]):
                pharmacology = 'none'
            else:
                pharmacology = single_stimulation_df['Pharmacology'].iloc[0]

            if single_stimulation_df['Type of protocol'].iloc[0] == 'Bsl':
                stimulation_frequency = np.NaN
                stimulation_duration = np.NaN
                stimulation_string = 'baseline'
            elif single_stimulation_df['Type of protocol'].iloc[0] == 'Bsl2':
                stimulation_frequency = np.NaN
                stimulation_duration = np.NaN
                stimulation_string = 'control_baseline'
            elif 'Hz' in single_stimulation_df['Type of protocol'].iloc[0]:
                stimulation_frequency = single_stimulation_df['Type of protocol'].iloc[0].replace(' Hz', '')
                stimulation_duration = int(single_stimulation_df['Timing of stim (s)'].iloc[0] * 1000)
                stimulation_string = f'{stimulation_frequency}-Hz_for_{stimulation_duration}-ms'
            sorted_data = {'recording_type': [recording_type],
                           'stimulation_string': [stimulation_string],
                           'stimulation_frequency-Hz': [stimulation_frequency],
                           'stimulation_duration-ms': [stimulation_duration],
                           'pharmacology': [pharmacology],
                           'time_since_cutting-M': [time_since_cutting]}
            return pd.DataFrame(data = sorted_data)
        
        
    def add_corresponding_filepaths(self, df_without_filepaths: DataFrame, cell_recordings_dir: Path) -> DataFrame:
        filepath_main_excel_sheet = cell_recordings_dir.joinpath(f'{cell_recordings_dir.name}.xlsx').as_posix()
        df_with_filepaths = df_without_filepaths.copy()
        df_with_filepaths['filepath_main_excel_sheet'] = [filepath_main_excel_sheet]*df_with_filepaths.shape[0]
        df_with_filepaths['filepath_detected_events'] = 'not_available'
        if self.recording_method == 'voltage_clamp':
            filepaths_of_detected_events = list()
            for elem in cell_recordings_dir.iterdir():
                if 'datapoints' not in elem.name:
                    if elem.name.endswith('.csv'):
                        if elem.name.startswith('.') == False:
                            filepaths_of_detected_events.append(elem)
            if len(filepaths_of_detected_events) > 0:
                stimulation_strings_from_filenames = self.get_stimulation_strings_from_filenames(filepaths_of_detected_events = filepaths_of_detected_events)
                df_with_filepaths = self.add_event_filepaths(df_without_event_filepaths = df_with_filepaths, strings_to_compare = stimulation_strings_from_filenames)
        else:
            df_with_filepaths
        return df_with_filepaths
    
    def add_event_filepaths(self, df_without_event_filepaths: DataFrame, strings_to_compare: Dict) -> DataFrame:
        df_strings_to_compare = pd.DataFrame(data=strings_to_compare)
        df = df_without_event_filepaths.copy()
        for stimulation_string, recording_type in zip(df['stimulation_string'].values, df['recording_type'].values):
            postsynaptic_current_type = recording_type[:recording_type.find('_')]
            try:
                filepath_detected_events = df_strings_to_compare.loc[(df_strings_to_compare['stimulation_string'] == stimulation_string) &
                                                                     (df_strings_to_compare['postsynaptic_current_type'] == postsynaptic_current_type),
                                                                     'filepath_detected_events'].values[0]
                df.loc[(df['stimulation_string'] == stimulation_string) & (df['recording_type'] == recording_type), 'filepath_detected_events'] = filepath_detected_events
            except IndexError:
                cell_identifier = f'{df["date"].iloc[0]}_{df["session_cell_id"].iloc[0]}'
                raise IndexError(f'Could not find a file with recorded events for cell: {cell_identifier} '
                                 f'for {recording_type} with stimulation type: {stimulation_string}. Please check '
                                 'the corresponding files!')
        return df
        

    def add_cross_validated_event_filepaths(self, df_without_event_filepaths: DataFrame, strings_to_compare: Dict) -> DataFrame:
        df = df_without_event_filepaths.copy()
        for expected_stimulation_string in df['stimulation_string'].values:
            if expected_stimulation_string in strings_to_compare['stimulation_string']:
                idx = strings_to_compare['stimulation_string'].index(expected_stimulation_string)
                df.loc[df['stimulation_string'] == expected_stimulation_string, 'filepath_detected_events'] = strings_to_compare['filepath_detected_events'][idx]
            else:
                cell_identifier = f'{df["date"].iloc[0]}_{df["session_cell_id"].iloc[0]}'
                raise ValueError(f'Could not find a file with recorded events for cell: {cell_identifier} and stimulation type: {expected_stimulation_string}')
        return df
                

    def get_stimulation_strings_from_filenames(self, filepaths_of_detected_events: List) -> Dict:
        stimulation_paradigms = {'stimulation_string': [],
                                 'stimulation_frequency-Hz': [],
                                 'stimulation_duration-ms': [], 
                                 'filepath_detected_events': [],
                                 'postsynaptic_current_type': []}
        for filepath in filepaths_of_detected_events:
            filename = filepath.name
            filename = filename.replace('.csv', '')
            if filename.endswith('_i'):
                postsynaptic_current_type = 'inhibitory'
                filename = filename.replace('_i', '')
            elif filename.endswith('_I'):
                postsynaptic_current_type = 'inhibitory'
                filename = filename.replace('_I', '')
            else:
                postsynaptic_current_type = 'excitatory'
            yyyy_mm_dd_ = filename[:11]
            filename = filename.replace(yyyy_mm_dd_, '')
            cell_id = filename[:filename.index('_')]
            filename = filename.replace(cell_id + '_', '')
            if 'Hz' in filename:
                stimulation_frequency = int(filename[:filename.index('Hz')])
                stimulation_duration = filename[filename.index('Hz') + 2 :]
                if 'ms' not in stimulation_duration:
                    stimulation_duration = int(stimulation_duration[:stimulation_duration.find('s')]) * 1000
                else:
                    stimulation_duration = int(stimulation_duration[:stimulation_duration.find('ms')])
                stimulation_paradigm = f'{stimulation_frequency}-Hz_for_{stimulation_duration}-ms'
            elif filename.endswith('Bsl'):
                stimulation_frequency = np.NaN
                stimulation_duration = np.NaN
                stimulation_paradigm = 'baseline'
            elif filename.endswith('Bsl2'):
                stimulation_frequency = np.NaN
                stimulation_duration = np.NaN
                stimulation_paradigm = 'control_baseline'            
            else:
                print(f'Warning: stimulation paradigm could not be identified for: {filepath.name}')
                stimulation_frequency = np.NaN
                stimulation_duration = np.NaN
                stimulation_paradigm = 'unknown'
            stimulation_paradigms['stimulation_string'].append(stimulation_paradigm)
            stimulation_paradigms['stimulation_frequency-Hz'].append(stimulation_frequency)
            stimulation_paradigms['stimulation_duration-ms'].append(stimulation_duration)
            stimulation_paradigms['filepath_detected_events'].append(filepath.as_posix())
            stimulation_paradigms['postsynaptic_current_type'].append(postsynaptic_current_type)
        return stimulation_paradigms

    
class ExtractVoltageClampData(DataExtractor):

    @property
    def recording_method(self) -> str:
        return 'voltage_clamp' 



class ExtractCurrentClampData(DataExtractor):

    @property
    def recording_method(self) -> str:
        return 'current_clamp' 
    
    
"""      
class CellRecording:
    
    def create_general_metadata_df(self, path_to_recordings_dir: Path) -> DataFrame:
        metadata_df = pd.read_excel(path_to_recordings_dir.joinpath(f'{path_to_recordings_dir.name}.xlsx'),
                                    sheet_name = 'General information')
        general_metadata = {'date': path_to_recordings_dir.name[:10],
                            'session_cell_id': path_to_recordings_dir.name[path_to_recordings_dir.name.rfind('_') + 1:],
                            'mouse_line': metadata_df['Animal line'][0],
                            'brain_region': metadata_df['Region'][0],
                            'cell_type': metadata_df['Type'][0],
                            'sex': metadata_df['Sex'][0],
                            'age': metadata_df['Age (days)'][0],
                            'stimulation_in_percent': metadata_df['Stimulation (%)'][0],
                            'internal_patch_solution': metadata_df['Internal used'][0]}
        return pd.DataFrame(general_metadata, index=[0])
    
    
    def create_stimulation_paradigms_df(self, path_to_recordings_dir: Path) -> DataFrame:
        filepaths_stimulation_recordings = list()
        for elem in path_to_recordings_dir.iterdir():
            if 'datapoints' not in elem.name:
                if elem.name.endswith('.csv'):
                    filepaths_stimulation_recordings.append(elem)
        stimulation_paradigms = {'stimulation_string': list(),
                                 'stimulation_frequency-Hz': list(),
                                 'stimulation_duration-ms': list(), 
                                 'filepath_detected_events': list()}
        for filepath in filepaths_stimulation_recordings:
            filename = filepath.name
            filename = filename.replace('.csv', '')
            yyyy_mm_dd_ = filename[:11]
            filename = filename.replace(yyyy_mm_dd_, '')
            cell_id = filename[:filename.index('_')]
            filename = filename.replace(cell_id + '_', '')
            if 'Hz' in filename:
                stimulation_frequency = int(filename[:filename.index('Hz')])
                stimulation_duration = filename[filename.index('Hz') + 2 :]
                if 'ms' not in stimulation_duration:
                    stimulation_duration = int(stimulation_duration[:stimulation_duration.find('s')]) * 1000
                else:
                    stimulation_duration = int(stimulation_duration[:stimulation_duration.find('ms')])
                stimulation_paradigm = f'{stimulation_frequency}-Hz_for_{stimulation_duration}-ms'          
            elif 'Bsl' in filename:
                stimulation_frequency = np.NaN
                stimulation_duration = np.NaN
                stimulation_paradigm = 'baseline'
            else:
                print(f'Warning: stimulation paradigm could not be identified for: {filepath.name}')
                stimulation_frequency = np.NaN
                stimulation_duration = np.NaN
                stimulation_paradigm = 'unknown'
            stimulation_paradigms['stimulation_string'].append(stimulation_paradigm)
            stimulation_paradigms['stimulation_frequency-Hz'].append(stimulation_frequency)
            stimulation_paradigms['stimulation_duration-ms'].append(stimulation_duration)
            stimulation_paradigms['filepath_detected_events'].append(filepath.as_posix())
        return pd.DataFrame(stimulation_paradigms)
    
    
    def create_recordings_overview(self, cell_recordings_dir: Path, global_cell_id: int) -> DataFrame:
        general_metadata = self.create_general_metadata_df(cell_recordings_dir)
        stimulation_paradigms = self.create_stimulation_paradigms_df(cell_recordings_dir)
        stimulations_adjusted_general_metadata = pd.concat([general_metadata]*stimulation_paradigms.shape[0], ignore_index=True)
        recordings_overview = pd.concat([stimulations_adjusted_general_metadata, stimulation_paradigms], axis=1)
        recordings_overview['global_cell_id'] = str(global_cell_id).zfill(4)
        return recordings_overview
"""

class Subdirectories:
    
    def __init__(self, root_dir: Path):
        self.root_dir = root_dir
        self.create_missing_subdirectories()
        
    def create_missing_subdirectories(self):
        subdirectories = listdir_nohidden(self.root_dir)
        
        try: self.single_cell_analyses = self.root_dir.joinpath([elem for elem in subdirectories if 'single_cell' in elem][0])
        except:
            self.single_cell_analyses = self.root_dir.joinpath('single_cell_analysis')
            os.mkdir(self.single_cell_analyses.as_posix()) 
        
        try: self.group_analyses = self.root_dir.joinpath([elem for elem in subdirectories if 'group' in elem][0])
        except:
            self.group_analyses = self.root_dir.joinpath('group_analysis')
            os.mkdir(self.group_analyses.as_posix()) 
            
        try: self.project_summary = self.root_dir.joinpath([elem for elem in subdirectories if 'project_summary' in elem][0])
        except:
            self.project_summary = self.root_dir.joinpath('project_summary')
            os.mkdir(self.project_summary.as_posix()) 
        
        try: self.exported_excel_sheets = self.root_dir.joinpath([elem for elem in subdirectories if 'exported_excel_sheets' in elem][0])
        except:
            self.exported_excel_sheets = self.root_dir.joinpath('exported_excel_sheets')
            os.mkdir(self.exported_excel_sheets.as_posix()) 