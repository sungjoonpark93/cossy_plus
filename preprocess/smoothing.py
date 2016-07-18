__author__ = 'SungJoonPark'
__author__ = 'SungJoonPark'

'''
Using preprocessed file, smoothing, or propagating the mutation matrix with PPI.
'''

import pandas as pd
import numpy as np
import time
import code.network.preprocess as network_preprocess
import code.preprocess.filter as filter


def network_smoothing(mut_preproceesed_filename=None,output_filename =None,network_file=None,alpha=0.7):
    '''
    preprocessed_filename:input preprocessed mutation filename
    network_file ; sif formatted network file
    alpha : random walk parameter. Weight on restart part
    '''
    print "network smoothing start\n" + "prprocessed mutation filename:"+mut_preprocessed_filename+"\n"+"network filename for smoothing:"+network_file+"\n"+"output filename"+output_filename
    s_time = time.time()
    n_iter=0

    mut_preprocessed_df = pd.read_csv(mut_preproceesed_filename,index_col=0)
    network_df = network_preprocess.get_network(network_file)
    mut_preprocessed_df,network_df = filter.match_gene(mut_preprocessed_df=mut_preprocessed_df,network_df=network_df)
    patient_mutation_df = mut_preprocessed_df.T

    F_0 = patient_mutation_df.copy(deep=True)
    F = patient_mutation_df.copy(deep=True)

    while True:
        n_iter = n_iter + 1

        F_t = F.copy(deep=True)
        F = (alpha*np.matmul(F,network_df)) + ((1-alpha)*F_0)

        if np.sum(np.asmatrix(np.square((F-F_t)))) < (1/1000000.0):
            print "network smoothing end. ", str(n_iter), "iterations, " ,str((time.time()-s_time) / 60.0) ," minutes elapsed\n"
            break


    (F.T).to_csv(output_filename)


def network_smoothing_with_exp(mut_preprocessed_filename=None,exp_preprocessed_filename=None,output_filenmae=None,network_file=None,alpha=0.7,reverse=False):
    '''
    mut_preprocessed_filename:input preprocessed mutation filename
    exp_preprocessed_filename : input preprocessed exp filename
    network_file ; sif formatted network file
    alpha : random walk parameter. Weight on restart part
    reverse : TRUE if you want to reverse mut_preprocessed matrix. 0 to 1 and 1 to 0.
    '''

    mut_preprocessed_df = pd.read_csv(mut_preprocessed_filename,index_col=0)
    if reverse==True:
        mut_preprocessed_df = mut_preprocessed_df.replace([0,1],[1,0])
    exp_preprocessed_df = pd.read_csv(exp_preprocessed_filename,index_col=0)
    print "network smoothing with exp start\n" + "prprocessed mutation filename:"+mut_preprocessed_filename+"\n"+"preprocessed exp filename:"+exp_preprocessed_filename+"\n"+"network filename for smoothing:"+network_file+"\n"+"output filename"+output_filename
    s_time = time.time()
    n_iter=0
    network_df = network_preprocess.get_network(network_file)
    mut_preprocessed_df,exp_preprocessed_df=filter.match_sample(mut_preprocessed_df=mut_preprocessed_df,exp_preprocessed_df=exp_preprocessed_df,conserve_exp_normal_sample=True)
    mut_preprocessed_df,network_df,exp_preprocessed_df = filter.match_gene(mut_preprocessed_df=mut_preprocessed_df,network_df=network_df,is_exp=True,exp_preprocessed_df=exp_preprocessed_df)

    patient_mutation_df = mut_preprocessed_df.T
    patient_exp_df = exp_preprocessed_df.T


    if (patient_mutation_df.shape != patient_exp_df.shape):
        raise Exception("mutation and exp data shape are not matched")

    F_0 = patient_mutation_df.copy(deep=True)
    F = patient_mutation_df.copy(deep=True)
    while True:
        n_iter = n_iter + 1

        F_t = F.copy(deep=True)
        F = np.multiply(((alpha*np.matmul(F,network_df)) + ((1-alpha)*F_0)), patient_exp_df/np.max(np.max(np.abs(patient_exp_df))))

        if np.sum(np.asmatrix(np.square((F-F_t)))) < (1/1000000.0):
            print "network smoothing end. ", str(n_iter), "iterations, " ,str((time.time()-s_time) / 60.0) ," minutes elapsed\n"
            break

    (F.T).to_csv(output_filenmae)




if __name__ == '__main__':
    pass
     #---mut with smoothing(reverse Fasle)---#
    reactome_filename = "Q:/COSSY+/data/network/reactome/PathwayCommons.8.reactome.BINARY_SIF_DRUG_deleted.hgnc.txt.sif"
    string_filename = "Q:/COSSY+/data/network/string/string_fix_excel_problem_tab_seperator.sif"
    mut_preprocessed_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Breast Cancer/mut/simple_somatic_mutation.open.BRCA-US.tsv_smoothing_string_preprocessed.csv"
    exp_preprocessed_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Breast Cancer/exp/exp_seq.BRCA-US.tsv_preprocessed.csv"
    output_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Breast Cancer/mut_with_exp/BRCA_mut_with_exp_seq_alpha0.7_string.csv"
    network_smoothing_with_exp(mut_preprocessed_filename=mut_preprocessed_filename,exp_preprocessed_filename=exp_preprocessed_filename,output_filenmae=output_filename,network_file=string_filename)
    output_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Breast Cancer/mut_with_exp/BRCA_mut_with_exp_seq_alpha0.7_reactome.csv"
    network_smoothing_with_exp(mut_preprocessed_filename=mut_preprocessed_filename,exp_preprocessed_filename=exp_preprocessed_filename,output_filenmae=output_filename,network_file=reactome_filename)


    mut_preprocessed_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Colon Adenocarcinoma/mut/simple_somatic_mutation.open.COAD-US.tsv_smoothing_string_preprocessed.csv"
    exp_preprocessed_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Colon Adenocarcinoma/exp/exp_seq.COAD-US.tsv_preprocessed.csv"
    output_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Colon Adenocarcinoma/mut_with_exp/COAD_mut_with_exp_seq_alpha0.7_string.csv"
    network_smoothing_with_exp(mut_preprocessed_filename=mut_preprocessed_filename,exp_preprocessed_filename=exp_preprocessed_filename,output_filenmae=output_filename,network_file=string_filename)
    output_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Colon Adenocarcinoma/mut_with_exp/COAD_mut_with_exp_seq_alpha0.7_reactome.csv"
    network_smoothing_with_exp(mut_preprocessed_filename=mut_preprocessed_filename,exp_preprocessed_filename=exp_preprocessed_filename,output_filenmae=output_filename,network_file=reactome_filename)

    mut_preprocessed_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Gastric Adenocarcinoma/mut/simple_somatic_mutation.open.STAD-US.tsv_smoothing_string_preprocessed.csv"
    exp_preprocessed_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Gastric Adenocarcinoma/exp/exp_seq.STAD-US.tsv_preprocessed.csv"
    output_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Gastric Adenocarcinoma/mut_with_exp/STAD_mut_with_exp_seq_alpha0.7_string.csv"
    network_smoothing_with_exp(mut_preprocessed_filename=mut_preprocessed_filename,exp_preprocessed_filename=exp_preprocessed_filename,output_filenmae=output_filename,network_file=string_filename)
    output_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Gastric Adenocarcinoma/mut_with_exp/STAD_mut_with_exp_seq_alpha0.7_reactome.csv"
    network_smoothing_with_exp(mut_preprocessed_filename=mut_preprocessed_filename,exp_preprocessed_filename=exp_preprocessed_filename,output_filenmae=output_filename,network_file=reactome_filename)


    mut_preprocessed_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Lung Squamous Cell Carcinoma/mut/simple_somatic_mutation.open.LUSC-US.tsv_smoothing_string_preprocessed.csv"
    exp_preprocessed_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Lung Squamous Cell Carcinoma/exp/exp_seq.LUSC-US.tsv_preprocessed.csv"
    output_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Lung Squamous Cell Carcinoma/mut_with_exp/LUSC_mut_with_exp_seq_alpha0.7_string.csv"
    network_smoothing_with_exp(mut_preprocessed_filename=mut_preprocessed_filename,exp_preprocessed_filename=exp_preprocessed_filename,output_filenmae=output_filename,network_file=string_filename)
    output_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Lung Squamous Cell Carcinoma/mut_with_exp/LUSC_mut_with_exp_seq_alpha0.7_reactome.csv"
    network_smoothing_with_exp(mut_preprocessed_filename=mut_preprocessed_filename,exp_preprocessed_filename=exp_preprocessed_filename,output_filenmae=output_filename,network_file=reactome_filename)


    mut_preprocessed_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Prostate Adenocarcinoma/mut/simple_somatic_mutation.open.PRAD-US.tsv_smoothing_string_preprocessed.csv"
    exp_preprocessed_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Prostate Adenocarcinoma/exp/exp_seq.PRAD-US.tsv_preprocessed.csv"
    output_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Prostate Adenocarcinoma/mut_with_exp/PRAD_mut_with_exp_seq_alpha0.7_string.csv"
    network_smoothing_with_exp(mut_preprocessed_filename=mut_preprocessed_filename,exp_preprocessed_filename=exp_preprocessed_filename,output_filenmae=output_filename,network_file=string_filename)
    output_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Prostate Adenocarcinoma/mut_with_exp/PRAD_mut_with_exp_seq_alpha0.7_reactome.csv"
    network_smoothing_with_exp(mut_preprocessed_filename=mut_preprocessed_filename,exp_preprocessed_filename=exp_preprocessed_filename,output_filenmae=output_filename,network_file=reactome_filename)


     #---mut with smoothing(reverse True)---#
    reactome_filename = "Q:/COSSY+/data/network/reactome/PathwayCommons.8.reactome.BINARY_SIF_DRUG_deleted.hgnc.txt.sif"
    string_filename = "Q:/COSSY+/data/network/string/string_fix_excel_problem_tab_seperator.sif"
    mut_preprocessed_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Breast Cancer/mut/simple_somatic_mutation.open.BRCA-US.tsv_smoothing_string_preprocessed.csv"
    exp_preprocessed_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Breast Cancer/exp/exp_seq.BRCA-US.tsv_preprocessed.csv"
    output_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Breast Cancer/mut_with_exp/BRCA_mut_with_exp_seq_reversed_alpha0.7_string.csv"
    network_smoothing_with_exp(mut_preprocessed_filename=mut_preprocessed_filename,exp_preprocessed_filename=exp_preprocessed_filename,output_filenmae=output_filename,network_file=string_filename,reverse=True)
    output_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Breast Cancer/mut_with_exp/BRCA_mut_with_exp_seq_reversed_alpha0.7_reactome.csv"
    network_smoothing_with_exp(mut_preprocessed_filename=mut_preprocessed_filename,exp_preprocessed_filename=exp_preprocessed_filename,output_filenmae=output_filename,network_file=reactome_filename,reverse=True)


    mut_preprocessed_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Colon Adenocarcinoma/mut/simple_somatic_mutation.open.COAD-US.tsv_smoothing_string_preprocessed.csv"
    exp_preprocessed_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Colon Adenocarcinoma/exp/exp_seq.COAD-US.tsv_preprocessed.csv"
    output_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Colon Adenocarcinoma/mut_with_exp/COAD_mut_with_exp_seq_reversed_alpha0.7_string.csv"
    network_smoothing_with_exp(mut_preprocessed_filename=mut_preprocessed_filename,exp_preprocessed_filename=exp_preprocessed_filename,output_filenmae=output_filename,network_file=string_filename,reverse=True)
    output_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Colon Adenocarcinoma/mut_with_exp/COAD_mut_with_exp_seq_reversed_alpha0.7_reactome.csv"
    network_smoothing_with_exp(mut_preprocessed_filename=mut_preprocessed_filename,exp_preprocessed_filename=exp_preprocessed_filename,output_filenmae=output_filename,network_file=reactome_filename,reverse=True)



    mut_preprocessed_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Gastric Adenocarcinoma/mut/simple_somatic_mutation.open.STAD-US.tsv_smoothing_string_preprocessed.csv"
    exp_preprocessed_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Gastric Adenocarcinoma/exp/exp_seq.STAD-US.tsv_preprocessed.csv"
    output_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Gastric Adenocarcinoma/mut_with_exp/STAD_mut_with_exp_seq_reversed_alpha0.7_string.csv"
    network_smoothing_with_exp(mut_preprocessed_filename=mut_preprocessed_filename,exp_preprocessed_filename=exp_preprocessed_filename,output_filenmae=output_filename,network_file=string_filename,reverse=True)
    output_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Gastric Adenocarcinoma/mut_with_exp/STAD_mut_with_exp_seq_reversed_alpha0.7_reactome.csv"
    network_smoothing_with_exp(mut_preprocessed_filename=mut_preprocessed_filename,exp_preprocessed_filename=exp_preprocessed_filename,output_filenmae=output_filename,network_file=reactome_filename,reverse=True)



    mut_preprocessed_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Lung Squamous Cell Carcinoma/mut/simple_somatic_mutation.open.LUSC-US.tsv_smoothing_string_preprocessed.csv"
    exp_preprocessed_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Lung Squamous Cell Carcinoma/exp/exp_seq.LUSC-US.tsv_preprocessed.csv"
    output_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Lung Squamous Cell Carcinoma/mut_with_exp/LUSC_mut_with_exp_seq_reversed_alpha0.7_string.csv"
    network_smoothing_with_exp(mut_preprocessed_filename=mut_preprocessed_filename,exp_preprocessed_filename=exp_preprocessed_filename,output_filenmae=output_filename,network_file=string_filename,reverse=True)
    output_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Lung Squamous Cell Carcinoma/mut_with_exp/LUSC_mut_with_exp_seq_reversed_alpha0.7_reactome.csv"
    network_smoothing_with_exp(mut_preprocessed_filename=mut_preprocessed_filename,exp_preprocessed_filename=exp_preprocessed_filename,output_filenmae=output_filename,network_file=reactome_filename,reverse=True)



    mut_preprocessed_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Prostate Adenocarcinoma/mut/simple_somatic_mutation.open.PRAD-US.tsv_smoothing_string_preprocessed.csv"
    exp_preprocessed_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Prostate Adenocarcinoma/exp/exp_seq.PRAD-US.tsv_preprocessed.csv"
    output_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Prostate Adenocarcinoma/mut_with_exp/PRAD_mut_with_exp_seq_reversed_alpha0.7_string.csv"
    network_smoothing_with_exp(mut_preprocessed_filename=mut_preprocessed_filename,exp_preprocessed_filename=exp_preprocessed_filename,output_filenmae=output_filename,network_file=string_filename,reverse=True)
    output_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Prostate Adenocarcinoma/mut_with_exp/PRAD_mut_with_exp_seq_reversed_alpha0.7_reactome.csv"
    network_smoothing_with_exp(mut_preprocessed_filename=mut_preprocessed_filename,exp_preprocessed_filename=exp_preprocessed_filename,output_filenmae=output_filename,network_file=reactome_filename,reverse=True)


