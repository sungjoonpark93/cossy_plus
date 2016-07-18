__author__ = 'SungJoonPark'
import pandas as pd
import os
import filter

def make_gct_file(df,output_filename=None):
    with open(output_filename,'w') as w:
        w.write("#1.2\n")
        num_row = len(df.index)
        num_column = len(df.columns)
        w.write(str(num_row)+"\t"+str(num_column)+"\n")
        w.write('Name'+'\t'+'Description'+'\t'+'\t'.join(df.columns)+"\n")
        for index in range(len(df.index)):
            gene_name = df.iloc[index,:].name
            value_list = map(lambda x:str(x) ,df.iloc[index,:].values.tolist())
            w.write(gene_name+"\t"+"---"+"\t"+"\t".join(value_list)+"\n")

def make_cls_file(ordered_df, tumor_normal_divide_index,output_filename=None):
    with open(output_filename,'w') as w:
        w.write(str(len(ordered_df.columns))+" "+'2'+' '+'1'+"\n")
        w.write(str('#'+'\n'))
        w.write(' '.join(['tumor']*len((ordered_df.iloc[:,0:tumor_normal_divide_index]).columns) + ['normal']*len((ordered_df.iloc[:,tumor_normal_divide_index:]).columns)))
        w.write('\n')

def make_kid_file(gct_filename=None,network='string',output_filename=None):
    if network =='string':
        HGNC_STRING_mapping_file = "Q:/COSSY+/data/mapping_file/HGNC_Entrez_STRING_fix_excel_problem.csv"
        mapping_df = pd.read_csv(HGNC_STRING_mapping_file)
        mapping_dict= mapping_df.set_index('HGNC_Symbol')['STRING_Locus_ID'].to_dict()
        hgnc_list=pd.read_csv(gct_filename, sep='\t',skiprows=2,index_col=0).index.tolist()
        #print hgnc_list
        with open(output_filename,'w') as w:
            for hgnc in hgnc_list:
                if hgnc not in mapping_dict:
                    w.write('-')
                    w.write('\n')
                elif hgnc in mapping_dict:
                    if '.' in str(mapping_dict[hgnc]):
                        w.write(mapping_dict[hgnc])
                        w.write('\n')
                    else:
                        w.write('-')
                        w.write('\n')

    elif network =='kegg':
    #for kegg since kegg id is equal to entrez, we use HGNC_ENTREZ mapping file
        HGNC_Entrez_mapping_file = "Q:/COSSY+/data/mapping_file/HGNC_Entrez_fix_excel_problem.csv"
        HGNC_kegg_id_dict = pd.read_csv(HGNC_Entrez_mapping_file).set_index('Approved Symbol')['Entrez Gene ID(supplied by NCBI)']
        hgnc_list=pd.read_csv(gct_filename, sep='\t',skiprows=2,index_col=0).index.tolist()
        with open(output_filename,'w') as w:
            for kegg_id in HGNC_kegg_id_dict[hgnc_list].values.tolist():
                if str(kegg_id) =='nan':
                    w.write('-')
                    w.write('\n')
                else:
                    w.write(str(int(kegg_id)))
                    w.write('\n')

    else:
        raise Exception("uncorrected network")


#make gctfile, clsfile,kidfile
def make_cossy_input(preproceed_filename,kid_type='kegg'):
    print(preproceed_filename)

    df = pd.read_csv(preproceed_filename,index_col=0)
    print('odering the column to make tumor samples come first..')
    ordered_df, tumor_normal_divide_index = filter.order_sample(df)
    print('end odering..\n')

    print('start making gctfile..')
    gct_filename = os.path.splitext(preproceed_filename)[0]+'.gct'
    make_gct_file(ordered_df,output_filename=gct_filename)
    print('end making gctfile..\n')

    print('start making clsfile..')
    cls_filename = os.path.splitext(preproceed_filename)[0]+'.cls'
    make_cls_file(ordered_df,tumor_normal_divide_index,output_filename=cls_filename)
    print('end making clsfile..\n')

    print('start making kidfile..')
    kid_filename = os.path.splitext(preproceed_filename)[0]+'-'+kid_type+'.kid'
    make_kid_file(gct_filename=gct_filename,network=kid_type,output_filename=kid_filename)
    print('end making kidfile..\n')

if __name__ =='__main__':
    pass
    #exp
    # preprocessed_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Breast Cancer/exp/exp_array.BRCA-US.tsv_preprocessed.csv"
    # make_cossy_input(preprocessed_filename)
    # preprocessed_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Breast Cancer/exp/exp_seq.BRCA-US.tsv_preprocessed.csv"
    # make_cossy_input(preprocessed_filename)


    #mut
    # preprocessed_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Breast Cancer/mut/simple_somatic_mutation.open.BRCA-US.tsv_smoothing_reactome_preprocessed.csv"
    # make_cossy_input(preprocessed_filename)
    # preprocessed_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Breast Cancer/mut/simple_somatic_mutation.open.BRCA-US.tsv_smoothing_string_preprocessed.csv"
    # make_cossy_input(preprocessed_filename)


    #mut with exp(reverse True)
    # preprocessed_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Breast Cancer/mut_with_exp/BRCA_mut_with_exp_seq_reversed_alpha0.7_reactome.csv"
    # make_cossy_input(preprocessed_filename)
    # preprocessed_filename ="Q:/COSSY+/data/preprocessed/TCGA/ICGC/release21/Breast Cancer/mut_with_exp/BRCA_mut_with_exp_seq_reversed_alpha0.7_string.csv"
    # make_cossy_input(preprocessed_filename)

