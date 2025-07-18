
# def check_in_group_II(class_II_peptide, class_I_seqs):
#     return 'peptide shared in MHC I and II' if any(class_I_pep in class_II_peptide for class_I_pep in class_I_seqs) else 'peptide not shared'

# def fill_group_II_status(data_frame):
#     data_frame['mhc_status'] = None

#     class_I_seqs = data_frame[data_frame['MHC Restriction - Class'] == 'I']['Epitope - Name'].unique()

#     # Only label MHC Class II peptides
#     data_frame.loc[data_frame['MHC Restriction - Class'] == 'II', 'mhc_status'] = \
#         data_frame[data_frame['MHC Restriction - Class'] == 'II']['Epitope - Name']\
#             .apply(lambda x: check_in_group_II(x, class_I_seqs))

#     #data_frame = data_frame[data_frame['MHC Restriction - Class'] == 'II']

#     return data_frame

# import re

# def fill_group_II_status(data_frame):
#     data_frame['mhc_status'] = None

#     # Get unique MHC class I peptide sequences
#     class_I_seqs = data_frame[data_frame['MHC Restriction - Class'] == 'I']['Epitope - Name'].unique()

#     # Build regex pattern from class I peptides (escaped for safety)
#     class_I_pattern = '|'.join(re.escape(p) for p in class_I_seqs)

#     # Filter MHC class II entries
#     mask_class_II = data_frame['MHC Restriction - Class'] == 'II'

#     # Vectorized substring search: does any class I peptide appear in each class II peptide?
#     data_frame.loc[mask_class_II, 'mhc_status'] = data_frame.loc[mask_class_II, 'Epitope - Name'] \
#         .str.contains(class_I_pattern) \
#         .map(lambda x: 'peptide shared in MHC I and II' if x else 'peptide not shared')

#     return data_frame

# import re

# def fill_group_II_status(data_frame):
#     data_frame['mhc_status'] = None

#     # Get unique peptide sequences
#     class_I_seqs = data_frame[data_frame['MHC Restriction - Class'] == 'I']['Epitope - Name'].unique()
#     class_II_seqs = data_frame[data_frame['MHC Restriction - Class'] == 'II']['Epitope - Name'].unique()

#     # Build regex patterns
#     class_I_pattern = '|'.join(re.escape(p) for p in class_I_seqs)
#     class_II_pattern = '|'.join(re.escape(p) for p in class_II_seqs)

#     # Masks
#     mask_I = data_frame['MHC Restriction - Class'] == 'I'
#     mask_II = data_frame['MHC Restriction - Class'] == 'II'

#     # Class II peptides: do they contain any class I peptide?
#     data_frame.loc[mask_II, 'mhc_status'] = data_frame.loc[mask_II, 'Epitope - Name'] \
#         .str.contains(class_I_pattern) \
#         .map(lambda x: 'peptide shared in MHC I and II' if x else 'peptide not shared')

# # 🔁 Reverse check: for each class I peptide, check if it's in any class II peptide
#     def check_if_in_class_II(peptide, class_II_list):
#         return any(peptide in longer_pep for longer_pep in class_II_list)

#     mask_I = data_frame['MHC Restriction - Class'] == 'I'
#     data_frame.loc[mask_I, 'mhc_status'] = data_frame.loc[mask_I, 'Epitope - Name'] \
#         .apply(lambda pep: 'peptide shared in MHC I and II' if check_if_in_class_II(pep, class_II_seqs) else 'peptide not shared')
import pandas as pd

def fill_group_II_status(data_frame):
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



if __name__ == "__main__":
    # Load the IEDB data and clean it
   fill_group_II_status(data_frame=None)  # Replace None with actual DataFrame when running
