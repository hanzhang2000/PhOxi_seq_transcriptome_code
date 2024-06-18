#!/usr/bin/env python3

import argparse
import pandas as pd

def filter_bamreadcount_data(Condition_sequencing_depth, Condition_multiallelic_deletion_sig, 
                             VAF_treated, VAF_untreated, VAF_treated_versus_untreated, 
                             VAF_treated_THUMPD3_kd_versus_control,df_treated,df_untreated,
                             df_treated_THUMPD3_kd,df_untreated_THUMPD3_kd):
    
    # (1) control

    # 1. filtering sites with sequencing depth >= 50
    df_treated = df_treated[df_treated['depth'] >= Condition_sequencing_depth] # [Condition 1]


    # 2 filtering sites with multiallelic or deletion signature

    if Condition_multiallelic_deletion_sig == 'Yes': # [Condition 2]

        # Create a boolean mask for rows where the values in the first two columns occur more than once
        mask1 = df_treated.duplicated(subset=[df_treated.columns[0], 
                                              df_treated.columns[1]], keep=False)

        # Create a boolean mask for rows where the 'base' column contains the character '-'
        mask2 = df_treated['base'].str.contains('-')

        # Combine the two masks using logical OR
        final_mask = mask1 | mask2

        # Apply the final mask to retain the desired rows
        df_treated = df_treated[final_mask]


    # Create a new DataFrame to store the combined values
    df_treated_vaf = df_treated.groupby(['chrom', 'position', 'ref']).agg({'vaf': 'sum', 'depth': 'max'}).reset_index()

    # ---------------------------------

    # 3. filtering sites with variant frequency in treated samples >= 10%

    df_treated_vaf = df_treated_vaf[df_treated_vaf['vaf'] >= VAF_treated] # [Condition 3]

    print(df_treated_vaf.head(5))


    # filtering untreated samples
    df_untreated = df_untreated[df_untreated['depth'] >= Condition_sequencing_depth] # [Condition 1]

    filtered_treated_position = df_treated_vaf[['chrom', 'position']]

    df_untreated_matched = pd.merge(filtered_treated_position, df_untreated, on=['chrom', 'position'], how='inner')
    column_names =['chrom', 'position', 'ref', 'vaf', 'depth']
    df_untreated_vaf = pd.DataFrame(columns = column_names, dtype=object)

    # Loop over the rows of DataFrame 1
    for index, row in df_treated_vaf.iterrows():

        chrom_value = row['chrom']
        position_value = row['position']
        ref_value = row['ref']

        matching_rows = df_untreated_matched[(df_untreated_matched['chrom'] == chrom_value) 
                                            & (df_untreated_matched['position'] == position_value)]
        
        if len(matching_rows) > 0:
            depth_value = matching_rows['depth'].tolist()

            matching_rows_filtered = matching_rows[matching_rows['ref'] != matching_rows['base']]
            vaf_values = matching_rows_filtered['vaf'].tolist()

            if len(vaf_values) > 0:
                combined_vaf = sum(vaf_values)
            else:
                combined_vaf = 0

            new_row = {'chrom': chrom_value, 'position': position_value, 'ref':ref_value, 
                       'vaf': combined_vaf, 'depth': depth_value[0]}

            df_untreated_vaf = df_untreated_vaf.append(new_row, ignore_index=True)
    
    print(df_untreated_vaf.head(5))
    

    # merge treated and untreated samples on the same dataframe
    df_combined_all_control = pd.merge(df_treated_vaf, df_untreated_vaf, on=['chrom', 'position'], suffixes=('_treated', '_untreated'))

    # ---------------------------------

    # 5. Compare treated versus
    # boolean mask:  VAF treated - VAF untreated >= VAF_treated_versus_untreated
    mask1 = df_combined_all_control['vaf_treated'] >= (df_combined_all_control['vaf_untreated'] + VAF_treated_versus_untreated) 
    mask2 = df_combined_all_control['vaf_untreated'] <= VAF_untreated
    final_mask = mask1 & mask2

    df_combined_all_control = df_combined_all_control[final_mask]
    
    print(df_combined_all_control.head(5))


    # (2) THUMPD3-kd


    # 1. filtering sites with sequencing depth >= 50
    df_treated_THUMPD3_kd = df_treated_THUMPD3_kd[df_treated_THUMPD3_kd['depth'] >= Condition_sequencing_depth] # [Condition 1]

    # Create a new DataFrame to store the combined values
    df_treated_THUMPD3_kd_vaf = df_treated_THUMPD3_kd.groupby(['chrom', 'position', 'ref']).agg({'vaf': 'sum', 'depth': 'max'}).reset_index()

    print(df_treated_THUMPD3_kd_vaf.head(5))

    df_untreated_THUMPD3_kd = df_untreated_THUMPD3_kd[df_untreated_THUMPD3_kd['depth'] >= Condition_sequencing_depth] # [Condition 1]
    
    filtered_treated_position = df_treated_THUMPD3_kd_vaf[['chrom', 'position']]
    df_untreated_matched = pd.merge(filtered_treated_position, df_untreated_THUMPD3_kd, on=['chrom', 'position'], how='inner')

    column_names =['chrom', 'position', 'ref', 'vaf', 'depth']
    df_untreated_THUMPD3_kd_vaf = pd.DataFrame(columns = column_names, dtype=object)

    # Loop over the rows of DataFrame 1
    for index, row in df_treated_THUMPD3_kd_vaf.iterrows():

        chrom_value = row['chrom']
        position_value = row['position']
        ref_value = row['ref']

        matching_rows = df_untreated_matched[(df_untreated_matched['chrom'] == chrom_value) 
                                            & (df_untreated_matched['position'] == position_value)]
        
        if len(matching_rows) > 0:
            depth_value = matching_rows['depth'].tolist()

            matching_rows_filtered = matching_rows[matching_rows['ref'] != matching_rows['base']]
            vaf_values = matching_rows_filtered['vaf'].tolist()

            if len(vaf_values) > 0:
                combined_vaf = sum(vaf_values)
            else:
                combined_vaf = 0

            new_row = {'chrom': chrom_value, 'position': position_value, 'ref':ref_value, 
                       'vaf': combined_vaf, 'depth': depth_value[0]}

            df_untreated_THUMPD3_kd_vaf = df_untreated_THUMPD3_kd_vaf.append(new_row, ignore_index=True)

    print(df_untreated_THUMPD3_kd_vaf.head(5))

    # 3. variant frequency in untreated samples <= 5%
    df_untreated_THUMPD3_kd_vaf = df_untreated_THUMPD3_kd_vaf[df_untreated_THUMPD3_kd_vaf['vaf'] <= VAF_untreated] # [Condition 4]

    # ---------------------------------

    # merge treated and untreated samples on the same dataframe
    df_combined_all_THUMPD3_kd = pd.merge(df_treated_THUMPD3_kd_vaf, df_untreated_THUMPD3_kd_vaf, on=['chrom', 'position'], suffixes=('_treated', '_untreated'))
    
    print(df_combined_all_THUMPD3_kd.head(5))
    
    # (3) Control versus THUMPD3-kd

    # Merge filtered dataframes of control and THUMPD3-kd data
    df_combined_all = df_combined_all_control.merge(df_combined_all_THUMPD3_kd, on=['chrom', 'position'], suffixes=('_control', '_THUMPD3_kd'))
    
    print(df_combined_all.head(5))

    # Create a boolean mask based on the condition that sequencing error of control_treated - THUMPD3-kd treated >= 5%
    mask1 = df_combined_all['vaf_treated_control'] >= (df_combined_all['vaf_treated_THUMPD3_kd'] + VAF_treated_THUMPD3_kd_versus_control) # [Condition 6]

    # Use the boolean mask to filter the dataframe
    df_combined_all = df_combined_all[mask1]

    df_combined_all = df_combined_all.drop_duplicates()
    
    print(df_combined_all.head(5))

    return df_combined_all, df_combined_all_control, df_combined_all_THUMPD3_kd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--control_treated', type=str, help='Path to control treated CSV file')
    parser.add_argument('--control_untreated', type=str, help='Path to control untreated CSV file')
    parser.add_argument('--THUMPD3_kd_treated', type=str, help='Path to THUMPD3-kd treated CSV file')
    parser.add_argument('--THUMPD3_kd_untreated', type=str, help='Path to THUMPD3-kd untreated CSV file')
    parser.add_argument('--output_combined_all', type=str, help='Path to save combined all CSV file')
    parser.add_argument('--output_combined_all_control', type=str, help='Path to save combined all control CSV file')
    parser.add_argument('--output_combined_all_THUMPD3_kd', type=str, help='Path to save combined all THUMPD3-kd CSV file')
    parser.add_argument('--Condition_sequencing_depth', type=int, default=50, help='Condition for sequencing depth')
    parser.add_argument('--Condition_multiallelic_deletion_sig', type=str, default='Yes', help='Condition for multiallelic or deletion signature')
    parser.add_argument('--VAF_treated', type=float, default=0.1, help='Variant frequency in treated samples')
    parser.add_argument('--VAF_untreated', type=float, default=0.05, help='Variant frequency in untreated samples')
    parser.add_argument('--VAF_treated_versus_untreated', type=float, default=0.1, help='Variant frequency in treated samples versus untreated samples')
    parser.add_argument('--VAF_treated_THUMPD3_kd_versus_control', type=float, default=0.05, help='Variant frequency in treated THUMPD3-kd samples versus control samples')

    args = parser.parse_args()

# Load the CSV files
df_treated = pd.read_csv(args.control_treated, sep='\t')
df_untreated = pd.read_csv(args.control_untreated, sep='\t')
df_treated_THUMPD3_kd = pd.read_csv(args.THUMPD3_kd_treated, sep='\t')
df_untreated_THUMPD3_kd = pd.read_csv(args.THUMPD3_kd_untreated, sep='\t')

# in the way code and data are, only THUMPD3-kd treated samples needs to be filtered
mask = df_treated_THUMPD3_kd.apply(lambda row: row['ref'] == row['base'], axis=1)
df_treated_THUMPD3_kd = df_treated_THUMPD3_kd[~mask]

# Call the main function
df_combined_all, df_combined_all_control, df_combined_all_THUMPD3_kd = filter_bamreadcount_data(
    args.Condition_sequencing_depth,
    args.Condition_multiallelic_deletion_sig,
    args.VAF_treated,
    args.VAF_untreated,
    args.VAF_treated_versus_untreated,
    args.VAF_treated_THUMPD3_kd_versus_control,
    df_treated,
    df_untreated,
    df_treated_THUMPD3_kd,
    df_untreated_THUMPD3_kd
)

# Save the output to CSV files
df_combined_all.to_csv(args.output_combined_all, index=False)
df_combined_all_control.to_csv(args.output_combined_all_control, index=False)
df_combined_all_THUMPD3_kd.to_csv(args.output_combined_all_THUMPD3_kd, index=False)



