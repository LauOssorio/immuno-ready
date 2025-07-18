
import pandas as pd

def fill_group_II_status(data_frame):
    """
    Label peptides as shared between MHC class I and II based on sequence containment.

    For each peptide in the DataFrame, checks whether:
    - A class II peptide contains any class I peptide
    - A class I peptide is contained in any class II peptide

    Adds a new column 'mhc_status' with the result.
    """
    data_frame['mhc_status'] = None

    # Get unique peptide sequences from each class
    class_I_peps = data_frame[data_frame['MHC Restriction - Class'] == 'I']['Epitope - Name'].unique()
    class_II_peps = data_frame[data_frame['MHC Restriction - Class'] == 'II']['Epitope - Name'].unique()

    for idx, row in data_frame.iterrows():
        peptide = row['Epitope - Name']
        mhc_class = row['MHC Restriction - Class']

        if pd.isna(peptide) or pd.isna(mhc_class):
            continue

        if mhc_class == 'II':
            # Check if this class II peptide contains any class I peptide
            if any(class_I_pep in peptide for class_I_pep in class_I_peps):
                data_frame.at[idx, 'mhc_status'] = 'peptide shared in MHC I and II'
            else:
                data_frame.at[idx, 'mhc_status'] = 'peptide not shared'

        elif mhc_class == 'I':
            # Check if this class I peptide is contained in any class II peptide
            if any(peptide in class_II_pep for class_II_pep in class_II_peps):
                data_frame.at[idx, 'mhc_status'] = 'peptide shared in MHC I and II'
            else:
                data_frame.at[idx, 'mhc_status'] = 'peptide not shared'

    return data_frame
