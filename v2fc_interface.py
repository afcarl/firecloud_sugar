import os
import pandas as pd
from firecloud import api as firecloud_api
import datetime
import glob

"""
 _____ ___ ____  _____ ____ _     ___  _   _ ____    ____  _   _  ____    _    ____  
|  ___|_ _|  _ \| ____/ ___| |   / _ \| | | |  _ \  / ___|| | | |/ ___|  / \  |  _ \ 
| |_   | || |_) |  _|| |   | |  | | | | | | | | | | \___ \| | | | |  _  / _ \ | |_) |
|  _|  | ||  _ <| |__| |___| |__| |_| | |_| | |_| |  ___) | |_| | |_| |/ ___ \|  _ < 
|_|   |___|_| \_\_____\____|_____\___/ \___/|____/  |____/ \___/ \____/_/   \_\_| \_\
"""

########################################################
# Sample functions
########################################################
def get_samples(paths_to_batches_info, google_bucket_id, sublist=None):
    """Compile samples from multiple batches
    Args: Self-explanatory
        - paths_to_samples_info: .xlsx file containing paths to files containing sample_info
        - sublist: list of tsca_ids to only compile data from certain batches. If None, compile data from all batches.
    Returns: 
        - df with samples from all batches
    """
    paths_to_samples = pd.read_excel(paths_to_batches_info, index_col=0)
    df_list = []

    for tsca_id, paths in paths_to_samples.iterrows():
        if sublist is not None and tsca_id not in sublist:
            continue
        # Make data Firecloud-compatible
        batch_data = prepare_batch_samples_for_metadata_export(paths.path_to_samples_info, tsca_id, google_bucket_id)
        df_list.append(batch_data)

    all_samples = pd.concat(df_list, axis=0)
    return all_samples

def get_samples_with_cohort(latest_tsca_id, paths_to_batches_info, google_bucket_id):
    """Retrieve df with all samples, including the cohort code they belong to
    """
    # Retrieve list of samples with corresponding cohort from bsp.broadinstitute.org
    samples_with_cohort = pd.read_excel('cohort_files/bsp_latest_all_samples_%s.xls'%latest_tsca_id)

    # All samples, without cohort data
    all_samples = get_samples(paths_to_batches_info, google_bucket_id)

    # Add cohort data to all samples
    data = pd.merge(all_samples, samples_with_cohort[['Sample ID', 'Collection']], \
                     left_on='bsp_sample_id_validation', \
                     right_on='Sample ID', \
                     how='inner') \
                    .drop(['Sample ID'], axis=1)

    # FC doesn't accept cohort names with non-alphanumeric characters, so use cohort codes instead
    # Load dictionary of {long cohort name : short cohort code}
    cohort_formatted_names = pd.read_table('cohort_files/cohort_names_dictionary.txt', \
                                           header=None, names=['Collection', 'cohort_code'])

    # Add cohort codes to data
    data = pd.merge(data, cohort_formatted_names, on='Collection', how='outer')
    return data

def prepare_batch_samples_for_metadata_export(path, tsca_id, google_bucket_id):
    """Prepare the file to export samples metadata to firecloud
    Args:
        path_id: path to file ending in {}.import_samples.txt
        tsca_id: TSCAXX
        google_bucket_id: id of google bucket ('gs://google_bucket_id')
    Returns:
        pd.DF of data ready for export
    Saves:
        ./{tsca_id}/fc_upload_samples_tsca_{tsca_id}.txt
    """
    # export raw data
    data = pd.read_table(path)
    
    # Rename columns to match firecloud requirements
    data = data.rename(columns={'sample_id':'entity:sample_id', 'individual_id':'participant_id'})
    # Locations of BAM files in google bucket
    path_in_bucket_full = "gs://%s/seq_data/%s" % (google_bucket_id, tsca_id)
    # Extract bam filename
    data['bam_filename'] = data.apply(lambda row: row['clean_bam_file_capture'].split('/')[-1], axis=1)
    # Create bai filename (change extension on .bam file)
    data['bai_filename'] = data.apply(lambda row: "%s%s" %(row['bam_filename'][:-3], 'bai'), axis=1)
    # Change BAM path from xchip to Google cloud
    data['clean_bam_file_capture'] = \
        data.apply( lambda row: "%s/%s/%s" \
                   %(path_in_bucket_full, row['external_id_validation'], row['bam_filename']), axis=1)
    # Add location of .bai file 
    data['clean_bai_file_capture'] = \
        data.apply( lambda row: "%s/%s/%s" \
                   %(path_in_bucket_full, row['external_id_validation'], row['bai_filename']), axis=1)
    # Add TSCA ID
    data['tsca_id'] = tsca_id

    return data

def save_and_upload_samples(data, namespace, workspace, tsca_id):
    """Create FC uploading file and upload patients to FC
    Args:
        - data: participants df
    Writes: 
        - {tsca_id}/fc_upload_patients_tsca_{tsca_id}.txt
    """
    os.system('mkdir -p %s'%tsca_id)
    filename = '%s/fc_upload_samples_tsca_%s.txt' % (tsca_id, tsca_id)
    data.to_csv(filename, sep='\t', index=False)
    res = upload_entities_from_tsv(namespace, workspace, filename)
    return res

########################################################
# Pairs list
########################################################
def create_pairs_list(all_samples):
    """Creates DF with pairs for firecloud metadata export.
    Args:
        - all_samples: all samples.
    """
    dfs = []
    # Find match normals for tumor samples only
    tumor_samples = all_samples[all_samples.sample_type=="Tumor"]
    i = 0
    for index, row in tumor_samples.iterrows():
        # Find all samples from same individual (same individual_id, different sample_id)
        patient_samples = all_samples[ (all_samples['participant_id'] == row['participant_id']) \
                                          & (all_samples['entity:sample_id'] != row['entity:sample_id']) ]

        # NOTE: If more than one match tumor tissue or match normal found, select one at random.
        # The match normal is used to compute allelic fractions in Mutect2, so for now we ignore the conditions it was grown in.

        ######## Match normal: Add match normal
        match_normal = patient_samples[ patient_samples['sample_type'] == "Normal"]
        #   > No match normal found
        if match_normal.empty: 
            control_sample_id = "NA"
            control_sample_tsca_id = "NA"
        #   > Match normal found
        elif match_normal.shape[0] > 0:
            match_normal = match_normal.iloc[0]
            control_sample_id = match_normal['entity:sample_id']
            control_sample_tsca_id = match_normal['tsca_id']
        
        # Create DF with Tumor/Normal pair set
        pair_id = "%s_%s_TN" % (row['entity:sample_id'], control_sample_id)
        df_dict = {'entity:pair_id': pair_id, 'case_sample_id': row['entity:sample_id'], \
                    'control_sample_id': control_sample_id, 'participant_id': row['participant_id'], 'match_type': 'tumor_normal', \
                    'case_sample_tsca_id': row['tsca_id'], 'control_sample_tsca_id': control_sample_tsca_id}
        dfs.append(pd.DataFrame(df_dict, index=[i], columns=df_dict.keys()))
        i+=1
        
        ######## Tumor tissue: Add primary tumor tissue
        match_primary_tumor = patient_samples[ patient_samples['external_id_validation'] \
                                              .str.contains('primary|prim|tissue|tiss|Primary|Tissue') ]
        #    > No primary tumor tissue found
        if match_primary_tumor.empty:
            control_sample_id = "NA"
            control_sample_tsca_id = "NA"
        #    > Sample itself is a primary tumor tissue
        elif any(substring in row['external_id_validation'] for substring in ['primary', 'prim', 'tissue', 'tiss', 'Primary', 'Tissue']):
            control_sample_id = "NA"
            control_sample_tsca_id = "NA"
        #    > Tumor tissue found
        elif match_primary_tumor.shape[0] > 0:
            match_primary_tumor = match_primary_tumor.iloc[0]
            control_sample_id = match_primary_tumor['entity:sample_id']
            control_sample_tsca_id = match_primary_tumor['tsca_id']
        
        # Create DF with Tumor/Primary pair set
        pair_id = "%s_%s_TP" % (row['entity:sample_id'], control_sample_id)
        df_dict = {'entity:pair_id': pair_id, 'case_sample_id': row['entity:sample_id'], \
                    'control_sample_id': control_sample_id, 'participant_id': row['participant_id'], 'match_type': 'tumor_primary', \
                    'case_sample_tsca_id': row['tsca_id'], 'control_sample_tsca_id': control_sample_tsca_id}
        dfs.append(pd.DataFrame(df_dict, index=[i], columns=df_dict.keys()))
        i+=1
   
    return pd.concat(dfs, axis=0)

def save_and_upload_pairs(namespace, workspace, pairs, blacklist=[]):
    """Updates pairs to firecloud. 
    NOTE: All pairs need to be updated with every new batch,
    as it may contain match normals or primary matches for previous batches.
    Args:
        - pairs: df of all pairs, as created by create_pairs_list
    Returns: 
        - res: json response from http request
    Creates: 
        - ./Pairs/fc_upload_pairs.txt file
    """
    pairs = pairs[ ~pairs['case_sample_id'].isin(blacklist)]
    os.system('mkdir -p Pairs')
    filename = './Pairs/fc_upload_pairs.txt'
    pairs.to_csv(filename, '\t', index=False)
    res = upload_entities_from_tsv(namespace, workspace, filename)
    return res

###############################################
# Participant functions
###############################################
def prepare_participants_for_metadata_export(paths_to_batches_info, tsca_id):
    """Create participant entities DF for Firecloud. 
    Participants need to exist before you can upload their respective samples
    Args:
        path_id: path to file ending in {}.import_samples.txt
    """    
    paths_to_samples = pd.read_excel(paths_to_batches_info, index_col=0)
    path_to_samples_info  = paths_to_samples.loc[tsca_id, 'path_to_samples_info']

    raw = pd.read_table(path_to_samples_info)
    print( "%d Participants in this batch" % raw['individual_id'].unique().shape[0] )
    # Data to upload
    data = pd.DataFrame(raw.individual_id.drop_duplicates()).rename(columns={'individual_id':'entity:participant_id'})
    return data

def save_and_upload_participants(data, namespace, workspace, tsca_id):
    """Create FC uploading file and upload patients to FC
    Args:
        - data: participants df
    Writes: 
        - {tsca_id}/fc_upload_patients_tsca_{tsca_id}.txt
    """
    os.system('mkdir -p %s'%tsca_id)
    filename = './%s/fc_upload_patients_%s.txt' % (tsca_id, tsca_id)
    data.to_csv(filename, '\t', index=False)
    res = upload_entities_from_tsv(namespace, workspace, filename)
    return res

########################################################
# Sample set functions
########################################################
def prepare_batch_sample_set_for_metadata_export(path, tsca_id):
    """Create dfs and write files with sample_set metadata for Firecloud.
    Three sample sets are created: all samples, normal samples, and tumor samples in batch
    A sample for a given batch 
    Args:
        - path: path to file ending in {}.import_samples.txt
        - tsca_id: batch tsca id
    Returns:
        - (all_samples, tumor_samples, normal_samples)
    Writes:
        - files for all_samples, tumor_samples, normal_samples
    """
    raw = pd.read_table(path)
    print( "%d Samples in this batch" % raw.shape[0] )

    # Create dfs to upload
    all_samples = pd.concat([pd.DataFrame(index=raw.index, columns=['membership:sample_set_id'], data=tsca_id), \
                      raw[ ['sample_id', 'sample_type'] ]], axis=1)


    tumors  = all_samples.loc[ all_samples['sample_type'] == "Tumor", ['membership:sample_set_id', 'sample_id'] ]
    tumors.loc[: , 'membership:sample_set_id'] = "%s_T"%tsca_id
    
    normals = all_samples.loc[ all_samples['sample_type'] == "Normal", ['membership:sample_set_id', 'sample_id'] ]
    normals.loc[: , 'membership:sample_set_id'] = "%s_N"%tsca_id

    all_samples = all_samples.drop('sample_type', axis=1)
    return (all_samples, tumors, normals)

def filter_existing_samples(df, sample_id_colname, remote_samples):
    """If sample set is uploaded more than once, samples will be duplicated, which is undesirable behavior.
    Therefore, need to download samples existing remotely (in FC), and remove them from sample set to be uploaded.
    Args:
        - df: dataframe to filter
        - sample_id_colname: name of column containing sample_id
    """
    remote_sample_ids = remote_samples['entity:sample_id'].tolist()
    df_clean = df[~df[sample_id_colname].isin(remote_sample_ids)]
    return df_clean

def save_and_upload_batch_sample_sets(batch_samples, batch_tumors, batch_normals, tsca_id, namespace, workspace):
    """Create FC uploading file and upload patients to FC
    """
    # Save to file
    os.system('mkdir -p %s'%tsca_id)
    batch_samples_filename = './%s/fc_upload_sample_set_tsca_%s.txt' % (tsca_id, tsca_id)
    batch_tumors_filename = './%s/fc_upload_sample_set_tsca_%s_tumors.txt' % (tsca_id, tsca_id)
    batch_normals_filename = './%s/fc_upload_sample_set_tsca_%s_normals.txt' % (tsca_id, tsca_id)
    
    batch_samples.to_csv(batch_samples_filename , sep="\t", index=False )
    batch_tumors.to_csv(batch_tumors_filename , sep="\t", index=False )
    batch_normals.to_csv(batch_normals_filename , sep="\t", index=False )

    r1 = upload_entities_from_tsv(namespace, workspace, batch_samples_filename)
    r2 = upload_entities_from_tsv(namespace, workspace, batch_tumors_filename)
    r3 = upload_entities_from_tsv(namespace, workspace, batch_normals_filename)
    return (r1, r2, r3)

########################################################
# Cohort functions
########################################################
def prepare_cohorts_for_metadata_export(paths_to_batches_info, google_bucket_id, latest_tsca_id, namespace, workspace, blacklist=[]):
    """Creates sample sets corresponding to cohorts for Firecloud export.
    Args:
        - tsca_id: latest tsca_id
        - paths_to_batches_info: file containing paths to all paths_to_batch_samples_info files (usually "paths_to_batches_info.xlsx")
    Returns: 
        - metadata with samples and cohorts they belong to for FC upload
    """
    # Retrieve list of samples with corresponding cohort from bsp.broadinstitute.org
    samples_with_cohort = pd.read_excel('cohort_files/bsp_latest_all_samples_%s.xls'%latest_tsca_id)

    # All samples, without cohort data
    all_samples = get_samples(paths_to_batches_info, google_bucket_id)

    # Add cohort data to all samples
    data = pd.merge(all_samples, samples_with_cohort[['Sample ID', 'Collection']], \
                     left_on='bsp_sample_id_validation', \
                     right_on='Sample ID', \
                     how='inner') \
                    .drop(['Sample ID'], axis=1)

    # FC doesn't accept cohort names with non-alphanumeric characters, so use cohort codes instead
    # Load dictionary of {long cohort name : short cohort code}
    cohort_formatted_names = pd.read_table('cohort_files/cohort_names_dictionary.txt', \
                                           header=None, names=['Collection', 'cohort_code'])

    # Add cohort codes to data
    data = pd.merge(data, cohort_formatted_names, on='Collection', how='outer')

    # Prepare for FC export format
    data = data.rename(columns={'cohort_code': 'membership:sample_set_id', 'entity:sample_id': 'sample_id'})
    data_clean = data[['membership:sample_set_id', 'sample_id']]

    # Remove blacklist
    data_clean = data_clean[ ~data_clean['sample_id'].isin(blacklist)]

    return data_clean

def save_and_upload_cohorts(data, latest_tsca_id, namespace, workspace):
    """Save and upload cohort metadata to FC
    """
    # Write FC import file 
    filename = 'cohort_files/fc_upload_sample_set_cohorts_%s.txt'%latest_tsca_id
    data.to_csv(filename, index=False, sep="\t")
    res = upload_entities_from_tsv(namespace, workspace, filename)
    return res

def prepare_cohort_pairsets_for_metadata_exports(latest_tsca_id, pairs, all_samples, blacklist=[]):
    """Create DF with cohort pairsets, used to create cohort reports for SNVs.
    Args:
        - pairs: pair list
        - all_samples: all samples with cohort
    Returns: 
        - DF with [cohort_code, pair_id]
    """
    # Create list of pairs
    clean_pairs_list = pairs[ ~pairs['case_sample_id'].isin(blacklist)]
    # Add cohorts to pairs
    pairs_with_cohort = pd.merge(clean_pairs_list, all_samples[['entity:sample_id', 'cohort_code']], \
             left_on='case_sample_id', right_on='entity:sample_id')
    # Prepare DF for FC export
    pairs_with_cohort_clean = pairs_with_cohort[['cohort_code', 'entity:pair_id']] \
        .rename(columns={'entity:pair_id': 'pair_id', 'cohort_code': 'membership:pair_set_id'})

    return pairs_with_cohort_clean

def save_and_upload_cohort_pairsets(namespace, workspace, pairsets):
    os.system('mkdir -p Pairs')
    filename = './Pairs/fc_upload_cohort_pairsets.txt'
    pairsets.to_csv(filename, '\t', index=False)
    res = upload_entities_from_tsv(namespace, workspace, filename)
    return res

def save_and_upload_cohort_all_samples(all_samples, name, namespace, workspace, blacklist=[]):
    """Create and upload cohort all samples
    Args: 
        - Self-explanatory
    """
    df = all_samples[['entity:sample_id']].rename(columns={'entity:sample_id': 'sample_id'})
    df['membership:sample_set_id'] = name

    # Re-arrange columns
    cols = ['membership:sample_set_id', 'sample_id']
    df = df[cols]

    # Blacklist
    df = df[ ~df['sample_id'].isin(blacklist) ]
    df.to_csv('all_samples/fc_upload_%s.txt'%name, index=None, sep="\t")
    res = upload_entities_from_tsv(namespace, workspace, 'all_samples/fc_upload_%s.txt'%name)
    return res

def save_and_upload_cohort_all_tumors(all_samples, name, namespace, workspace, blacklist=[]):
    """Create and upload cohort of all tumor samples across all batches
    Args: 
        - Self-explanatory
        - name: cohort name (usually Cum_Tumors_{LATEST_TSCA_ID}_all)
        - paths_to_samples_info: .xlsx file containing paths to files containing sample_info
    """
    tumor_samples = all_samples[all_samples.sample_type == "Tumor"]

    # Prepare column names
    df = tumor_samples[['entity:sample_id']].rename(columns={'entity:sample_id': 'sample_id'})
    df['membership:sample_set_id'] = name

    # Re-arrange columns
    cols = ['membership:sample_set_id', 'sample_id']
    df = df[cols]

    # Blacklist
    df = df[ ~df['sample_id'].isin(blacklist) ]
    df.to_csv('tumor_samples/fc_upload_%s.txt'%name, index=None, sep="\t")
    res = upload_entities_from_tsv(namespace, workspace, 'tumor_samples/fc_upload_%s.txt'%name)
    return res

########################################################
# PoN Functions
########################################################
def create_panel_of_normals_from_batch(batch_id, paths_to_batches_info, N=20):
    """Create panel of normals with samples from a given batch. Add N random samples from other batches.
    Args:
        batch_id: (int) batch_id of batch in question
        paths_to_samples_info: .xlsx file containing paths to files containing sample_info
        N: (int) number of samples in panel of normals
    Returns: 
        pon: pon df
        name: pon name
    """

    paths_to_samples = pd.read_excel(paths_to_batches_info, index_col=0)
    path_to_samples_info  = paths_to_samples.loc[batch_id, 'path_to_samples_info']
    all_samples = pd.read_table(path_to_samples_info)
    normal_samples = all_samples[all_samples['sample_type'] == "Normal"]

    # Iterate over batches except the one in question (to select N random normal samples from them)
    dfs = []
    for i, row in paths_to_samples.loc[~paths_to_samples.index.isin([batch_id])].iterrows():
        df_tmp = pd.read_table(row['path_to_samples_info'])
        dfs.append(df_tmp)
    df = pd.concat(dfs, axis=0)

    other_normal_samples = df[df['sample_type'] == "Normal"].sample(N)
    
    name = 'PoN_%s_plus_%s_random' %(batch_id, N)
    pon = pd.concat([normal_samples, other_normal_samples], axis=0)
    pon['membership:sample_set_id'] = name
    pon = pon[['membership:sample_set_id', 'sample_id']]

    return pon, name

def upload_pon(pon_df, pon_name, namespace, workspace):
    """Upload PoN to FC
    Args:
        - pon_df: dataframe with normal samples in PoN
        - pon_name: name of PoN
    """
    os.system('mkdir -p PoNs')
    filename = './PoNs/fc_upload_PoN_%s.txt' % (pon_name)
    pon_df.to_csv(filename, '\t', index=False)
    res = upload_entities_from_tsv(namespace, workspace, 'PoNs/fc_upload_PoN_%s.txt'%pon_name)
    return res

################################################
# Helper Functions
###############################################
def upload_entities_from_tsv(namespace, workspace, entities_tsv_file):
    """Upload entities from tsv file
    Args: 
        Self-explanatory
        entities_tsv_file: path to tsv file
    Returns: 
        HTTP Response
    """
    res = firecloud_api.upload_entities_tsv(namespace, workspace, entities_tsv=entities_tsv_file)
    return res

def delete_pair(namespace, workspace, pair_id):
    """Delete pair from workspace/namespace
    """
    body = [{"entityType": "pair", "entityName": pair_id}]
    res = firecloud_api.delete_entities(namespace, workspace, body)
    return res

def delete_pair_set(namespace, workspace, pair_set_id):
    """Delete pair set from workspace/namespace
    Args: 
        Self-explanatory
    Returns: 
        HTTP Response
    """
    body = [{"entityType": "pair_set", "entityName": pair_set_id}]
    res = firecloud_api.delete_entities(namespace, workspace, body)
    return res

def delete_sample(namespace, workspace, sample_id):
    """Delete sample from workspace/namespace
    Args: 
        Self-explanatory
    Returns: 
        HTTP Response
    """
    body = [{"entityType": "sample", "entityName": sample_id}]
    res = firecloud_api.delete_entities(namespace, workspace, body)
    return res

def delete_sample_set(namespace, workspace, sample_set_id):
    """Delete sample set from workspace/namespace
    Args: 
        Self-explanatory
    Returns: 
        HTTP Response
    """
    body = [{"entityType": "sample_set", "entityName": sample_set_id}]
    res = firecloud_api.delete_entities(namespace, workspace, body)
    return res

def delete_participant(namespace, workspace, participant_id):
    """Delete participant from workspace/namespace
    Args: 
        Self-explanatory
    Returns: 
        HTTP Response
    """
    body = [{"entityType": "participant", "entityName": participant_id}]
    res = firecloud_api.delete_entities(namespace, workspace, body)
    return res

def download_remote_samples(namespace, workspace):
    """Download remote samples from Firecloud
    Writes:
        - remote_samples.txt: samples in FC
    """
    res = firecloud_api.get_entities_tsv(namespace, workspace, "sample")
    with open('remote_files/remote_samples.txt', 'w') as outfile:
        outfile.write(res.text)
    return

################################################
### MAIN
###############################################
def upload_data_new_batch(tsca_id, latest_tsca_id, path_to_samples_info, namespace, workspace, google_bucket_id):
    """Upload data for new batch (sample set)
    This function only needs to be called once for a new batch. However, if it is run multiple times, it's ok as FC won't 
    create any duplicates.
    Args:
        - Self-explanatory
        - path_to_samples_info: path to *.import_samples.txt file
    """
    participants = prepare_participants_for_metadata_export("paths_to_batches_info.xlsx", tsca_id)
    r1 = save_and_upload_participants(participants, namespace, workspace, tsca_id)

    batch_samples = prepare_batch_samples_for_metadata_export(path_to_samples_info, tsca_id, google_bucket_id)
    r2 = save_and_upload_samples(batch_samples, namespace, workspace, tsca_id)

    all_samples = get_samples_with_cohort(latest_tsca_id, "paths_to_batches_info.xlsx", google_bucket_id)
    pairs = create_pairs_list(all_samples)
    r3 = save_and_upload_pairs(namespace, workspace, pairs)

    # Sample set
    batch_sample_set, batch_tumor_set, batch_normal_set = prepare_batch_sample_set_for_metadata_export(path_to_samples_info, tsca_id)
    r4a, r4b, r4c = save_and_upload_batch_sample_sets(batch_sample_set, batch_tumor_set, batch_normal_set, tsca_id, namespace, workspace)
    
    # Sample set
    pon, name = create_panel_of_normals_from_batch(tsca_id, "paths_to_batches_info.xlsx")
    r5 = upload_pon(pon, name, namespace, workspace)

    return (r1, r2, r3, r4a, r4b, r4c, r5)

def update_cohorts(latest_tsca_id, namespace, workspace, google_bucket_id):
    """Update cohorts (sample sets that span multiple batches)
    This function needs to be called once for a new batch. Before updating a cohort sample set, we need to remove samples that 
    already belong to that cohort remotely, because if we don't they will be duplicated.
    """
    # Necessary to be up to date
    download_remote_samples(namespace, workspace)
    # Pre-requisites
    all_samples = get_samples_with_cohort(latest_tsca_id, "paths_to_batches_info.xlsx", google_bucket_id)
    pairs = create_pairs_list(all_samples)

    # DF of remote samples
    remote_samples = pd.read_table('remote_files/remote_samples.txt')
    # DF of remote [sample < > sample set ]
    remote_sample_sets = pd.read_table('remote_files/sample_set_membership_%s.tsv'%latest_tsca_id)
    # DF of remote [pair < > pair set]
    remote_pair_sets = pd.read_table('remote_files/pair_set_membership_%s.tsv'%latest_tsca_id)

    #### UPDATE COHORT SAMPLE SETS
    # DF of [samples < > sample set] to be updated
    cohorts = prepare_cohorts_for_metadata_export("paths_to_batches_info.xlsx", google_bucket_id, latest_tsca_id, namespace, workspace, blacklist=[])
    # Remove the samples that already belong to the cohort in FC 
    sample_ids_in_cohort = remote_sample_sets['sample'].tolist()
    cohorts_clean = cohorts[~cohorts['sample_id'].isin(sample_ids_in_cohort)]
    r1 = save_and_upload_cohorts(cohorts_clean, latest_tsca_id, namespace, workspace)

    #### UPDATE COHORT PAIR SETS
    # Retrieve cohort pairsets
    cohort_pairsets = prepare_cohort_pairsets_for_metadata_exports(latest_tsca_id, pairs, all_samples, blacklist=[])
    # Remove the pairs that already belong to the cohort in FC
    pair_ids_in_cohort = remote_pair_sets['pair'].tolist()
    cohort_pairsets_clean = cohort_pairsets[~cohort_pairsets['pair_id'].isin(pair_ids_in_cohort)]
    r2 = save_and_upload_cohort_pairsets(namespace, workspace, cohort_pairsets_clean)

    # Remove samples that already exist in FC
    remote_sample_ids = remote_samples['entity:sample_id'].tolist()
    all_samples_clean = all_samples[~all_samples['entity:sample_id'].isin(remote_sample_ids)]
    r3 = save_and_upload_cohort_all_samples(all_samples_clean, "Cum_%s_all"%latest_tsca_id, namespace, workspace, blacklist=[])
    r4 = save_and_upload_cohort_all_tumors(all_samples_clean, "Cum_Tumors_%s_all"%latest_tsca_id, namespace, workspace, blacklist=[])

    return (r1, r2, r3, r4)


