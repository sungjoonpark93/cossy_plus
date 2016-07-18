__author__ = 'SungJoonPark'

import pandas as pd
import numpy as np
import code.network.smoothing as smoothing
import filter

#preprocessing TCGA data came from ICGC
#preprocessed df means gene x patient dataframe.



def preprocess_exp(input_filename=None,output_filename=None,exp_type='seq'):
    #make patient x gene matrix from raw TCGA exp file from ICGG and output to file
    print('start preprocessing expression data\n' +"type:"+exp_type+'\n'+"input filename:" + input_filename)
    df = pd.read_csv(input_filename,sep='\t')
    if exp_type=='seq':
        preprocessed_df = df.pivot_table(index=['gene_id'],columns=['submitted_sample_id'],values= 'normalized_read_count',aggfunc=np.mean)
        preprocessed_df=preprocessed_df.fillna(0)
        if "?" in preprocessed_df.index:
            preprocessed_df = preprocessed_df.drop('?')
    elif exp_type=='array':
        preprocessed_df = df.pivot_table(index=['gene_id'],columns=['submitted_sample_id'],values= 'normalized_expression_value',aggfunc=np.mean)
        #there could be nan value in array data.Then, fill with 0
        preprocessed_df=preprocessed_df.fillna(0)
        if "?" in preprocessed_df.index:
            preprocessed_df = preprocessed_df.drop('?')

    preprocessed_df.to_csv(output_filename)


def preprocess_mut(input_filename=None,output_filename=None, score_type='binary'):
    """
    smoothingnetwork_filename is only worth when is_smoohting='True'
    smoothingnetwork_filename should be sif format.
    smoothing network takes a quite long time
    """
    print('start preprocessing mutation\n' +"score type:"+score_type+'\n'+"input filename:" + input_filename)

    df = pd.read_csv(input_filename,sep='\t',low_memory=False)
    #there is no gene symbol in ICGC mutation file, only Ensemble ID, so wee need to make gene symbol from EnsembleID
    mapping_file = "Q:/COSSY+/data/mapping_file/HGNC_ApprovedSymbol_EnsembleID.txt"
    map_df = pd.read_csv(mapping_file,sep='\t')
    map_dict = map_df.set_index('Ensembl ID(supplied by Ensembl)')['Approved Symbol'].to_dict()
    #filter the ensemble id into the one in the mapping file.
    df = df[df['gene_affected'].isin(map_dict.keys())]
    df['Gene Symbol'] = [map_dict[ensemble_id] for ensemble_id in df['gene_affected']]

    #function for score the mutation type.
    def get_score(consequence_type):
        if consequence_type in ['disruptive_inframe_deletion','disruptive_inframe_insertion','exon_variant','frameshift_variant','inframe_deletion','inframe_insertion','missense_variant','stop_gained']:
            score=1
            return score
        else:
            score=0
            return score
    #make score field and give the mutation score
    df['score'] = map(get_score,df['consequence_type'])

    if score_type == 'binary':
        #Even though patient has multiple mutations on the gene, only counting mutation once.
        tumor_sample_preprocessed_df = df.pivot_table(values='score',index='Gene Symbol',columns='submitted_sample_id',aggfunc=np.max).fillna(0)
    elif score_type =='freq':
        #If patient has multiple mutations on the gene, we sum the count, difference is aggfunc.
        tumor_sample_preprocessed_df = df.pivot_table(values='score',index='Gene Symbol',columns='submitted_sample_id',aggfunc=np.sum).fillna(0)

    #We need to concat normal sample matrix in which value is all 0
    normal_sample_list = list(set(df['submitted_matched_sample_id']))
    zero_matrix = np.zeros((len(tumor_sample_preprocessed_df.index),len(normal_sample_list)))
    normal_sample_preprocessed_df = pd.DataFrame(zero_matrix,index=tumor_sample_preprocessed_df.index,columns=normal_sample_list)
    preprocessed_df = pd.concat([tumor_sample_preprocessed_df,normal_sample_preprocessed_df],axis=1)


    preprocessed_df.to_csv(output_filename)


if __name__ == '__main__':
    pass
    #---exp---#
    # input_filename = "Q:/COSSY+/data/raw_data/TCGA/ICGC/release21/Breast Cancer/exp_array.BRCA-US.tsv/exp_array.BRCA-US.tsv"
    # output_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Breast Cancer/exp/exp_array.BRCA-US.tsv_preprocessed.csv"
    # preprocess_exp(input_filename=input_filename, output_filename=output_filename,exp_type='array')


     #---mut---#
    # input_filename = "Q:/COSSY+/data/raw_data/TCGA/ICGC/release21/Breast Cancer/simple_somatic_mutation.open.BRCA-US.tsv/simple_somatic_mutation.open.BRCA-US.tsv"
    # output_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Breast Cancer/mut/simple_somatic_mutation.open.BRCA-US.tsv_preprocessed.csv"
    # preprocess_mut(input_filename=input_filename, output_filename=output_filename,score_type='binary')



