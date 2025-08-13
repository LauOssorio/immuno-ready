import numpy as np


def create_target_features(data_frame):
    """
    Creates binary and continuous immunogenicity target features for peptides.

    This function performs the following:
    - Maps immunization conditions to a binary target (`peptide_base`)
    - Assigns strength scores to immune responses (`peptide_strength`)
    - Aggregates peptide-level responses across duplicate entries
    - Computes a final immunogenicity strength target (`target_strength`)

    Parameters:
    -----------
    data_frame : pd.DataFrame
        Input dataset with epitope metadata and assay results.

    Returns:
    --------
    pd.DataFrame
        DataFrame with a new 'target_strength' column and non-informative
        rows removed. Intermediate columns used in calculations are dropped.
    """

    # ## Hard coding the decisions on immunogenicity
    # conditions_immuno = [
    #     data_frame["1st in vivo Process - Process Type"] == 'Occurrence of infectious disease',
    #     data_frame["1st in vivo Process - Process Type"] =='Occurrence of allergy',
    #     data_frame["1st in vivo Process - Process Type"] =='Exposure with existing immune reactivity without evidence for disease',
    #     data_frame["1st in vivo Process - Process Type"] == 'Unknown',
    #     data_frame["1st in vivo Process - Process Type"] == 'Occurrence of autoimmune disease',
    #     data_frame["1st in vivo Process - Process Type"] =='Environmental exposure to endemic/ubiquitous agent without evidence for disease',
    #     data_frame["1st in vivo Process - Process Type"] =='Prophylactic vaccination',
    #     data_frame["1st in vivo Process - Process Type"] =='Administration in vivo',
    #     data_frame["1st in vivo Process - Process Type"] =='Occurrence of disease',
    #     data_frame["1st in vivo Process - Process Type"] =='Exposure without evidence for disease',
    #     #data_frame["1st in vivo Process - Process Type"] =='Occurrence of cancer',
    #     data_frame["1st in vivo Process - Process Type"] =='Documented exposure without evidence for disease',
    #     data_frame["1st in vivo Process - Process Type"] =='Transplant/transfusion',
    #     data_frame["1st in vivo Process - Process Type"] =='Vaccination',
    #     data_frame["1st in vivo Process - Process Type"] =='Therapeutic vaccination',
    #     data_frame["1st in vivo Process - Process Type"] =='Administration in vivo to cause disease',
    #     data_frame["1st in vivo Process - Process Type"] == np.nan,
    #     data_frame["1st in vivo Process - Process Type"] =='Administration in vivo to prevent or reduce disease',
    #     data_frame["1st in vivo Process - Process Type"] == 'None'
    #     ]

    # choices_immuno = [1,
    #     1,
    #     1,
    #     "unknown",
    #     1,
    #     1,
    #     1,
    #     1,
    #     1,
    #     1,
    #     #1,
    #     1,
    #     1,
    #     1,
    #     1,
    #     1,
    #     "unknown",
    #     1,
    #     0]

    # choices_immuno = [str(x) for x in choices_immuno]

    # data_frame["peptide_base"] = np.select(conditions_immuno, choices_immuno, default = "other")

    # # filter the unkown out
    # data_frame = data_frame[data_frame["peptide_base"] != "unknown"]
    # data_frame = data_frame[data_frame["peptide_base"] != "other"]

    # # Convert to integer the column peptide safety
    # data_frame['peptide_base'] = data_frame['peptide_base'].fillna(0).astype(int)

    # # if a peptide has multiple entries, keep the immunogenic (1)
    # peptide_base_sum = data_frame.groupby("Epitope - Name")["peptide_base"].sum()
    # data_frame["peptide_base_sum"] = data_frame["Epitope - Name"].map(peptide_base_sum)
    # data_frame['peptide_base'] = data_frame['peptide_base_sum'].clip(upper=1)

    #     # Immunity strength

    conditions_immuno = [data_frame["Assay - Qualitative Measurement"] == 'Positive',
                        data_frame["Assay - Qualitative Measurement"] =='Positive-Low',
                        data_frame["Assay - Qualitative Measurement"] =='Positive-Intermediate',
                        data_frame["Assay - Qualitative Measurement"] =='Positive-High',
                        data_frame["Assay - Qualitative Measurement"] =='Assay - Qualitative Measurement',
                        data_frame["Assay - Qualitative Measurement"] =='Negative'
                        ]


    choices_immuno = [
        1,
        1,
        1,
        1,
        1,
        0
    ]
    data_frame["target_strength"] = np.select(conditions_immuno, choices_immuno, default = "other")

    return data_frame


    # data_frame["peptide_strength"] = np.select(conditions_strength, choices_strength, default = np.nan)

    # data_frame["averaged_strength"] = (
    #     data_frame.groupby("Epitope - Name")["peptide_strength"]
    #     .transform("mean")
    # )



    # data_frame["target_strength"] = data_frame["peptide_base"]

    # return data_frame.drop(columns=[
    #     "Assay - Qualitative Measurement",
    #     "1st in vivo Process - Process Type",
    #     "peptide_base",
    #     "peptide_strength",
    #     "averaged_strength",
    #     "peptide_base_sum"
    # ])
