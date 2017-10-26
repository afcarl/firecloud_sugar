import os
import pandas as pd
from firecloud import api as firecloud_api
import datetime
import glob

################################################
### Interfacing with FC entities
################################################

blacklist = ["DW039-Tumor-SM-DB2IF"]

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

def delete_pair(namespace, workspace, pair_id):
    """Delete pair from workspace/namespace
    """
    body = [{"entityType": "pair", "entityName": pair_id}]
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

def delete_workspace_config(namespace, workspace, cnamespace, config):
    """Delete workspace configuration
    Args: 
        Self-explanatory
    Returns: 
        HTTP Response
    """
    res = firecloud_api.delete_workspace_config(namespace, workspace, cnamespace, config)
    return res

def delete_entity_attributes(namespace, workspace, entity_type, entity_name, attrs):
    """Delete entity attributes
    Args: 
    - attrs: list of attributes to delete
    """
    attr_update = [{"op": "RemoveAttribute", "attributeName":  attr} for attr in attrs]
    res = firecloud_api.update_entity(namespace, workspace, entity_type, entity_name, attr_update)
    return res

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

################################################
### Get functions
################################################
def get_sample_by_id(sample_id):
    """Get sample given its sample_id
    Args:
        - sample_id
    """
    all_samples = fc_interface.get_samples_multiple_batches(paths_to_samples_info, google_bucket_id)
    return all_samples[all_samples['entity:sample_id'] == sample_id]

################################################
### PoN Functions
################################################
def create_panel_of_normals(paths, N, name):
    """Create panel of normals sample set for Firecloud from multiple TSCA batches.
    Randomly select N samples from samples present in files listed in paths.
    Args:
        paths: (list) paths to file ending in {}.import_samples.txt
        N: (int) number of samples in panel of normals
        name: (string) name of Panel of Normals
    """
    dfs = [ pd.read_table(paths[0]) ]
    for i, path in enumerate(paths[1:]):
        df_to_concat = pd.read_table(path)
        dfs.append(df_to_concat)
    df = pd.concat(dfs, axis=0)
    # Shuffle samples to pick from all batches
    df = df.sample(frac=1).reset_index(drop=True)
    normals = df[df.sample_type=="Normal"][:N]['sample_id']
    if N==-1: print ("Creating panel of %d normals" %normals.shape[0])
    else: print ("Creating panel of %d normals" %N)
    
    # Compile data
    data = pd.concat([pd.DataFrame(index=normals.index, columns=['membership:sample_set_id'], data=name), \
                        normals], axis=1)

    return data

def upload_cohort_all_tumors(paths_to_samples_info, google_bucket_id, name, namespace, workspace, blacklist=[]):
    """Create and upload cohort of all tumor samples across all batches
    Args: 
        - Self-explanatory
        - name: cohort name (usually Cum_Tumors_{LATEST_TSCA_ID}_all)
        - paths_to_samples_info: .xlsx file containing paths to files containing sample_info
    """
    all_samples = get_samples_multiple_batches(paths_to_samples_info, google_bucket_id)
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

def upload_cohort_all_samples(paths_to_samples_info, google_bucket_id, name, namespace, workspace, blacklist=[]):
    """Create and upload cohort all samples
    Args: 
        - Self-explanatory
    """
    all_samples = get_samples_multiple_batches(paths_to_samples_info, google_bucket_id)
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

def create_panel_of_normals_from_batch(batch_id, paths_to_samples_info, N=20):
    """Create panel of normals with samples from a given batch. Add N random samples from other batches.
    Args:
        batch_id: (int) batch_id of batch in question
        paths_to_samples_info: .xlsx file containing paths to files containing sample_info
        N: (int) number of samples in panel of normals
    Returns: 
        pon: pon df
        name: pon name
    """

    paths_to_samples = pd.read_excel(paths_to_samples_info, index_col=0)
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

    res = fc_interface.upload_entities_from_tsv(namespace, workspace, cohorts_sample_set_metadata)

def upload_cohorts(namespace, workspace, tsca_id):
    """Upload cohort metadata to FC
    Args: 
        - self-explanatory
    """
    cohorts_sample_set_metadata = "cohort_files/fc_upload_sample_set_cohorts_%s.txt"%tsca_id
    res = upload_entities_from_tsv(namespace, workspace, cohorts_sample_set_metadata)
    return res

################################################
### Samples Functions
################################################
def get_samples_multiple_batches(paths_to_samples_info, google_bucket_id, sublist=None):
    """Compile samples from multiple batches
    Args: Self-explanatory
        - paths_to_samples_info: .xlsx file containing paths to files containing sample_info
        - sublist: list of tsca_ids to only compile data from certain batches. If None, compile data from all batches.
    Returns: 
        - df with samples from all batches
    """
    paths_to_samples = pd.read_excel(paths_to_samples_info, index_col=0)
    df_list = []

    for tsca_id, paths in paths_to_samples.iterrows():
        if sublist is not None and tsca_id not in sublist:
            continue
        # Make data Firecloud-compatible
        batch_data = prepare_batch_samples_for_metadata_export(paths.path_to_samples_info, tsca_id, google_bucket_id)
        df_list.append(batch_data)

    all_samples = pd.concat(df_list, axis=0)
    return all_samples

def merge_walkupseq_files(paths_to_batches_info):
    """Merge all walkupseq files in the walkupseq directory
    """
    paths_to_batches_info = glob.glob('walkupseq_files/*sample_info*')

    dfs = []
    for f in paths_to_batches_info:
        tmp = pd.read_table(f, encoding='latin1')
        dfs.append(tmp)

    df = pd.concat(dfs, axis=0)
    df.to_csv('walkupseq_files/walkupseq_all_combined.txt', sep="\t", index=None)
    return df

def download_remote_samples(namespace, workspace):
    """Download remote samples from Firecloud
    Writes:
        - remote_samples.txt: samples in FC
    """
    res = firecloud_api.get_entities_tsv(namespace, workspace, "sample")
    with open('remote_samples.txt', 'w') as outfile:
        outfile.write(res.text)
    return
################################################
### Pairs Functions
################################################s
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
            match_normal = match_normal.sample(n=1)
            control_sample_id = match_normal['entity:sample_id'].item()
            control_sample_tsca_id = match_normal['tsca_id'].item()
        
        # Create DF with Tumor/Normal pair set
        pair_id = "%s_%s_TN" % (row['entity:sample_id'], control_sample_id)
        df_dict = {'entity:pair_id': pair_id, 'case_sample_id': row['entity:sample_id'], \
                    'control_sample_id': control_sample_id, 'participant_id': row['participant_id'], 'match_type': 'tumor_normal', \
                    'case_sample_tsca_id': row['tsca_id'], 'control_sample_tsca_id': control_sample_tsca_id}
        dfs.append(pd.DataFrame(df_dict, index=[i], columns=df_dict.keys()))
        i+=1
        
        ######## Tumor tissue: Add primary tumor tissue
        match_primary_tumor = patient_samples[ patient_samples['external_id_validation'] \
                                              .str.contains('primary|prim|tissue|tiss') ]
        #    > No primary tumor tissue found
        if match_primary_tumor.empty: 
            control_sample_id = "NA"
            control_sample_tsca_id = "NA"
        #    > Tumor tissue found
        elif match_primary_tumor.shape[0] > 0:
            match_primary_tumor = match_primary_tumor.sample(n=1)
            control_sample_id = match_primary_tumor['entity:sample_id'].item()
            control_sample_tsca_id = match_primary_tumor['tsca_id'].item()
        
        # Create DF with Tumor/Primary pair set
        pair_id = "%s_%s_TP" % (row['entity:sample_id'], control_sample_id)
        df_dict = {'entity:pair_id': pair_id, 'case_sample_id': row['entity:sample_id'], \
                    'control_sample_id': control_sample_id, 'participant_id': row['participant_id'], 'match_type': 'tumor_primary', \
                    'case_sample_tsca_id': row['tsca_id'], 'control_sample_tsca_id': control_sample_tsca_id}
        dfs.append(pd.DataFrame(df_dict, index=[i], columns=df_dict.keys()))
        i+=1
   
    return pd.concat(dfs, axis=0)

def upload_pairs(namespace, workspace, pairs):
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
    os.system('mkdir -p Pairs')
    filename = './Pairs/fc_upload_pairs.txt'
    pairs.to_csv(filename, '\t', index=False)
    res = upload_entities_from_tsv(namespace, workspace, filename)
    return res

def update_pair_attrs(namespace, workspace, pairs, attrs):
    """Update the FC (remote) pair attributes, listed in @attrs, present in the @pairs dataframe.
    The purpose of this function is to "push" to Firecloud changes to pairs made locally.
    Args:
        - pairs: df of pairs to update
        - attrs: list of attributes to update
    """
    for idx, row in pairs.iterrows():
        attr_update = [{"op": "AddUpdateAttribute", "attributeName":  attr, \
                        "addUpdateAttribute": row[attr]} for attr in attrs]
        res = firecloud_api.update_entity(namespace, workspace, "pair", row['entity:pair_id'], attr_update)
    return

def create_pairsets(all_samples, pairs):
    """Creates pair sets for FC export, both tumor-normal and tumor-primary pairs
    Args:
        - all_samples: all samples
        - pairs: pairs df as created by create_pairs_list
    Returns: 
        - Tuple of dfs: (tumor-normal, tumor-primary)
    """
    tn_pairs = pairs[pairs['match_type'] == "tumor_normal"]
    tp_pairs = pairs[pairs['match_type'] == "tumor_primary"]

    tn_pairsets = pd.merge(tn_pairs, all_samples[['entity:sample_id', 'tsca_id']], \
                            left_on='case_sample_id', right_on='entity:sample_id', \
                            how='inner')[['tsca_id', 'entity:pair_id']] \
                            .rename(columns={'tsca_id': 'membership:pair_set_id', 'entity:pair_id': 'pair_id'})

    tp_pairsets = pd.merge(tp_pairs, all_samples[['entity:sample_id', 'tsca_id']], \
                            left_on='case_sample_id', right_on='entity:sample_id', \
                            how='inner')[['tsca_id', 'entity:pair_id']] \
                            .rename(columns={'tsca_id': 'membership:pair_set_id', 'entity:pair_id': 'pair_id'})

    # Append _TN/_TP to the end of the tumor-normal/tumor-primary pair set ids
    tn_pairsets['membership:pair_set_id'] = tn_pairsets['membership:pair_set_id'].apply(lambda x: "%s_TN"%x)
    tp_pairsets['membership:pair_set_id'] = tp_pairsets['membership:pair_set_id'].apply(lambda x: "%s_TP"%x)
    
    return (tn_pairsets, tp_pairsets)

def upload_pairsets(namespace, workspace, pairsets, pairset_type):
    os.system('mkdir -p Pairs')
    filename = './Pairs/fc_upload_pairsets_%s.txt'%pairset_type
    pairsets.to_csv(filename, '\t', index=False)
    res = upload_entities_from_tsv(namespace, workspace, filename)
    return res

################################################
### Preparations for FC export
# These functions create files used to export data to FC
################################################
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

    data.to_csv('%s/fc_upload_samples_tsca_%s.txt' % (tsca_id, tsca_id), sep='\t', index=False)
    return data

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
    # Save to file
    os.system('mkdir -p %s'%tsca_id)
    all_samples.to_csv( './%s/fc_upload_sample_set_tsca_%s.txt' % (tsca_id, tsca_id), sep="\t", index=False )
    tumors.to_csv( './%s/fc_upload_sample_set_tsca_%s_tumors.txt' % (tsca_id, tsca_id), sep="\t", index=False )
    normals.to_csv( './%s/fc_upload_sample_set_tsca_%s_normals.txt' % (tsca_id, tsca_id), sep="\t", index=False )
    return (all_samples, tumors, normals)

def prepare_patients_for_metadata_export(path_to_samples_info, tsca_id):
    """Create participant entities file for Firecloud. 
    Patients need to exist before you can upload their respective samples
    Args:
        path_id: path to file ending in {}.import_samples.txt
        tsca_id: tsca id
    Saves: 
        ./tsca_id/fc_upload_patients_tsca_{tsca_id}.csv:
            contains patient ids in tsca batch
    """    
    raw = pd.read_table(path_to_samples_info)
    print( "%d Participants in this batch" % raw['individual_id'].unique().shape[0] )
    # Data to upload
    data = pd.DataFrame(raw.individual_id.drop_duplicates()).rename(columns={'individual_id':'entity:participant_id'})
    os.system('mkdir -p %s'%tsca_id)
    filename = './%s/fc_upload_patients_tsca_%s.txt' % (tsca_id, tsca_id)
    data.to_csv(filename, '\t', index=False)
    return data

def prepare_cohorts_for_metadata_export(paths_to_batches_info, google_bucket_id, blacklist, tsca_id, namespace, workspace):
    """Creates sample sets corresponding to cohorts for Firecloud export.
    Args:
        - tsca_id: latest tsca_id
        - paths_to_batches_info: file containing paths to all paths_to_batch_samples_info files (usually paths_to_batches_info.xlsx)
    Writes:
        - in cohort_files: fc_upload_sample_set_{cohort_name}.txt
    Returns: 
        - data_clean 
    """
    # Retrieve list of samples with corresponding cohort from bsp.broadinstitute.org
    samples_with_cohort = pd.read_excel('cohort_files/bsp_latest_all_samples_%s.xls'%tsca_id)

    # All samples, without cohort data
    all_samples = get_samples_multiple_batches(paths_to_batches_info, google_bucket_id)

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

    # Samples in new batch
    samples_new_batch = pd.read_table("%s/fc_upload_sample_set_tsca_%s.txt"%(tsca_id, tsca_id))
    
    # Only upload new samples. As of 10/17, Firecloud doesn't filter duplicate entries
    data_clean = data_clean[ data_clean['sample_id'].isin(samples_new_batch.sample_id.tolist()) ]

    # Write FC import file 
    data_clean.to_csv('cohort_files/fc_upload_sample_set_cohorts_%s.txt'%tsca_id, index=False, sep="\t")
    
    return data_clean

def prepare_all_data_for_metadata_export(tsca_id, path_to_batch_samples_info, paths_to_batches_info, google_bucket_id):    
    """Prepare all batch metadata for uploading to Firecloud
    Args:
        - path_to_batch_samples_info: path to info on batch samples
    """
    prepare_patients_for_metadata_export(path_to_batch_samples_info, tsca_id)
    prepare_batch_sample_set_for_metadata_export(path_to_batch_samples_info, tsca_id)
    prepare_batch_samples_for_metadata_export(path_to_batch_samples_info, tsca_id, google_bucket_id)
    prepare_cohorts_for_metadata_export(paths_to_batches_info, google_bucket_id)
    return

################################################
### Prepare metadata
# These functions prepare the walkupseq and BSP metadata
################################################
def update_walkup_sheets(tsca_id):
    """The walkupseq sheets contain experimental data (collaborator, tumor type, tissue site, etc.).
    Compile all the sheets onto one walkupseq_all_combined_TSCAXX.txt
    Args: 
    - tsca_id: latest TSCA id added
    Writes: 
    - walkupseq_files: walkupseq_all_combined_TSCA{XX}.txt: combined data from all walkupseq files
    """
    files = glob.glob('/Users/mimoun/production_fc/metadata_exports/walkupseq_files/tsca*')
    dfs = []
    for f in files:
        dfs.append( pd.read_table(f, comment='#', encoding='latin1'))
    combined = pd.concat(dfs, axis=0)
    fname = "/Users/mimoun/production_fc/metadata_exports/walkupseq_files/walkupseq_all_combined_%s.txt"%tsca_id
    combined.to_csv(fname, sep="\t", index=None)
    return combined

################################################
### Export batch metadata to FC
################################################
def export_batch_metadata_to_fc(tsca_id, namespace, workspace):
    """Export metadata into Firecloud
    Args: 
        - tsca_id
    """
    patient_metadata            = "%s/fc_upload_patients_tsca_%s.txt" % (tsca_id, tsca_id)
    sample_set_metadata         = "%s/fc_upload_sample_set_tsca_%s.txt" % (tsca_id, tsca_id)
    tumors_sample_set_metadata  = "%s/fc_upload_sample_set_tsca_%s_tumors.txt" % (tsca_id, tsca_id)
    normals_sample_set_metadata = "%s/fc_upload_sample_set_tsca_%s_normals.txt" % (tsca_id, tsca_id)
    samples_metadata            = "%s/fc_upload_samples_tsca_%s.txt" % (tsca_id, tsca_id)
    pon_metadata                = "PoNs/fc_upload_PoN_sample_set_tsca_%s_PoN.txt" %(tsca_id)
    # cohorts_sample_set_metadata = "cohort_files/fc_upload_sample_set_cohorts.txt"

    # Upload metadata
    r1 = upload_entities_from_tsv(namespace, workspace, patient_metadata)
    r2 = upload_entities_from_tsv(namespace, workspace, samples_metadata)
    r3 = upload_entities_from_tsv(namespace, workspace, sample_set_metadata)
    r4 = upload_entities_from_tsv(namespace, workspace, tumors_sample_set_metadata)
    r5 = upload_entities_from_tsv(namespace, workspace, normals_sample_set_metadata)
    # r6 = upload_entities_from_tsv(namespace, workspace, cohorts_sample_set_metadata)
    return (r1, r2, r3, r4, r5) #, r6)

def update_batch_metadata(tsca_id, path_to_batch_samples_info, paths_to_batches_info, namespace, workspace, google_bucket_id):
    """Upload to Firecloud all the metadata necessary to run TSCA pipeline on new batch
    Args: 
        - tsca_id: tsca_id of the batch to run TSCA pipeline on
        - path_to_batch_samples_info: ends in *.export_samples.txt
    """
    # Prepare all metadata for batch
    prepare_all_data_for_metadata_export(tsca_id, path_to_batch_samples_info, paths_to_batches_info, google_bucket_id)
    # Upload to Firecloud
    export_batch_metadata_to_fc(tsca_id, namespace, workspace)
    return

################################################
### Download method configs
################################################
def download_method_configs(namespace, workspace):
    """Download the latest method configs being used in firecloud, and write them to file
    """
    configs = firecloud_api.list_workspace_configs(namespace, workspace).json()
    configs_list = []
    for config in configs:
        df_list = [ config['name'], \
                   config['methodRepoMethod']['methodName'], \
                   config['methodRepoMethod']['methodVersion'], \
                   config['rootEntityType']]
        configs_list.append(df_list)

    method_configs = pd.DataFrame(data=configs_list, columns=['method_config', 'method_name', 'snapshot', 'entity'])
    timestamp = datetime.datetime.now().strftime("%Y_%m_%d_%H:%M")
    method_configs.to_csv('method_configs/%s_method_configs.txt'%timestamp, index=False, sep="\t")
    method_configs.to_csv('method_configs/latest_method_configs.txt', index=False, sep="\t")
    return

def download_remote_wdls(namespace, workspace):
    """Update WDL scripts in wdls/production directory, to the ones currently being used in Firecloud
    """
    download_method_configs(namespace, workspace)
    method_configs = pd.read_table('method_configs/latest_method_configs.txt')
    for idx, method in method_configs.iterrows():
        res = firecloud_api.get_repository_method('tsca', method['method_name'], method['snapshot'])
        print("Updating WDL for %s:%s"%(method['method_name'], method['snapshot']))
        if res.status_code == 200:
            text_file = open("../wdls/production/%s.wdl"%method['method_name'], "w")
            text_file.write(res.json()['payload'])
            text_file.close()
    return

def main():
    """Run for every new batch.
    It is also necessary to upload metadata for all previous batches, because this batch can contain 
    match normals for samples in previous batches.
    """
    print ("For new batch: please add batch_id and path_to_batch_info to paths_to_samples_info.xlsx file")
    s = input('Have you updated the paths_to_samples_info.xlsx file? (Y/N)')
    if s == "N":
        print("Please do so before proceeding...")
        return
    path_to_all_samples_info = "paths_to_samples_info.xlsx"
    batches_info = pd.read_excel(path_to_all_samples_info)
#     for idx, batch in batches_info.iterrows():
#         update_batch_metadata(batch.tscaid, batch.path_to_samples_info, path_to_all_samples_info)