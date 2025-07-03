

def check_in_group_II(seq, group_II_seqs):
    return 'peptide found in MHC II' if any(seq in s for s in group_II_seqs) else 'peptide only in MHC I'

def fill_group_II_status(data_frame):
    data_frame['mhc_status'] = None

    class_II_seqs = data_frame[data_frame['MHC Restriction - Class'] == 'II']['Epitope - Name'].unique()

    data_frame.loc[data_frame['MHC Restriction - Class'] == 'I', 'mhc_status'] = data_frame[data_frame['MHC Restriction - Class'] == 'I']['Epitope - Name'].apply(
        lambda x: check_in_group_II(x, class_II_seqs)
    )

    data_frame = data_frame[data_frame["MHC Restriction - Class"] == "I"]

    return data_frame


if __name__ == "__main__":
    # Load the IEDB data and clean it
   fill_group_II_status(data_frame=None)  # Replace None with actual DataFrame when running
