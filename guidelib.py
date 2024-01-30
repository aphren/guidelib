#!/usr/bin/env python

## Author: Andrew Hren
## Email: Andrew.Hren@colorado.edu

from Bio import SeqIO
from Bio import SeqFeature
import pandas
import time
from collections import Counter
from ast import literal_eval

#clear the existing file and add a header line



def indexer(PAM,sequence,component):      
#Runs thru the full genbank sequence and finds all instances of the passed PAM
#Runs until the index goes past its limit, becoming -1.
    
    indices = []
    index_count = 0
    new_index = 0
    working_index = 0
    while working_index != -1:
        working_index = sequence.find(PAM,new_index)
        if working_index != -1:
            indices.append(working_index)
            index_count += 1
            new_index = working_index + 1
    return [indices, index_count]




def find_spacer(component,PAM_index,PAM):
#Takes in the full genome component data, the PAM index of interest, 
#and the PAM that was indexed (GG or CC), which determines the location
#from which you pull the spacer    
#Pams are indexed asymmetrically on the two strands, 
#so the number adjustments are used to select the correct spacer location     
    
    sequence = component.seq
    
    #if component.annotations['topology'] == 'circular':
    edge_sequence = sequence[-21:] + sequence[:21] 
    
    
    fwd_pams = ['GG','AG']
    rev_pams = ['CC','CT']

    
    if PAM in fwd_pams and PAM_index > 20:
        spacer = sequence[PAM_index-21:PAM_index-1]
    elif PAM in fwd_pams and PAM_index <= 20:
        spacer = edge_sequence[PAM_index:PAM_index+20]
        
    if PAM in rev_pams and PAM_index < len(sequence)-20:
        spacer = sequence[PAM_index+3:PAM_index+23].reverse_complement()
        
    elif PAM in rev_pams and PAM_index >= len(sequence)-20:
        new_index = PAM_index-len(sequence[:-21])
        spacer = edge_sequence[new_index+3:new_index+23].reverse_complement()
    
    return spacer


#next(it,None)
#skips first loop iteration for debugging


def create_list(PAM_type,filename):
#calls the indexer and find_spacer functions to gather all the PAM sites
#and their relevant spacers.
#This info is combined with the strand and genome location, then 
#returned as a pandas dataframe.
    
    data = {}
    PAM_count = 0
    master_index = []
    master_component = []
    master_strand = []
    master_spacer = []
    
    it = iter(enumerate(SeqIO.parse(filename,'genbank')))
    
    for index, record in it:
        data[str(record.id)] = record
        #save file data to a dictionary
        component = data[str(record.id)]
        print('Moving to component ' + component.id)
    
        print('Length: ' + str(len(component.seq)))
        #print(sequence)
        sequence = component.seq
        
        if PAM_type == 'NGG':        
            PAMs = ['GG','CC']
        if PAM_type == 'NAG':
            PAMs = ['AG','CT']

        
        for PAM in PAMs:
            func_output = indexer(PAM,sequence,component)
            index_list = func_output[0]
            index_count = func_output[1]
            for x in index_list:
                master_index.append(x)
                
            
            PAM_count += index_count
            #spacer_list = []
            
            #Link genome component to PAM
            genome_component = [component.id]*index_count
            for x in genome_component:
                master_component.append(x)
            #Link source strand to PAM
            if PAM == 'GG' or PAM == 'AG':
                strand = ['+']*index_count
                for x in strand:
                    master_strand.append(x)
                    
            else:
                strand = ['-']*index_count
                for x in strand:
                    master_strand.append(x)
            for val in index_list:
                spacer = find_spacer(component,val,PAM)
                master_spacer.append(str(spacer))


    index_series = pandas.Series(master_index)
    all_pams = index_series.to_frame(name='PAM_index')
    all_pams['PAM_type'] = pandas.Series([PAM_type]*len(all_pams))
    all_pams['Spacer'] = pandas.Series(master_spacer)
    all_pams['Strand'] = pandas.Series(master_strand)
    all_pams['Component'] = pandas.Series(master_component)
    all_pams['index'] = all_pams.index
    
    
    
    
    return all_pams


def gene_targeting_filter(data_df,filename):
#Selects all NGG pams which target within a CDS and adds the locus tag
#returns both the noncoding (best) guides and coding (filler) guides    
    
    print('gene targeting filter')
    data = {}
    feature_index = []
    cod_feature_index = []
    
    noncoding_df_fin = pandas.DataFrame(data=None,columns=['PAM_index','PAM_type','Spacer','Strand','Component','index'])
    coding_df_fin = pandas.DataFrame(data=None,columns=['PAM_index','PAM_type','Spacer','Strand','Component','index'])
    data_df = data_df[data_df['PAM_type'] =='NGG']

    it = iter(enumerate(SeqIO.parse(filename,'genbank')))
    
    #next(it,None)
    
    for index, record in it:
        data[str(record.id)] = record
        #save file data to a dictionary
      
        component = data[str(record.id)]
        #this loop goes through all CDS's and creates a list of all base indices
        #that are within a CDS, essentially cutting out intergeneic regions
        
    
        loop_data = data_df[data_df['Component'] ==  str(component.id)]
        for feature in component.features:
            if feature.type == 'CDS' or feature.type == 'tRNA':
                
                
                #location adjustments extend search frame to ensure all
                #guides which cut inside the CDS are picked up
                if str(feature.location)[-2] == '+':
                    cds_df = loop_data[loop_data['PAM_index'].isin(range(feature.location.start-5,feature.location.end-6))]
                if str(feature.location)[-2] == '-':
                    cds_df = loop_data[loop_data['PAM_index'].isin(range(feature.location.start+5,feature.location.end+4))]
                noncoding_df = cds_df[cds_df['Strand'] != str(feature.location)[-2]]
                coding_df =  cds_df[cds_df['Strand'] == str(feature.location)[-2]]           
                
                #appending dataframes is slower than generating a list and 
                #later using .isin(), but there are duplicate entries thanks to
                #overlapping CDS's. By doing it this slower way, I can maintain
                #the order and add a column with the corresopnding locus_tags
                #without worrying about the order being switched 
                
                noncoding_df_fin = noncoding_df_fin.append(noncoding_df)
                coding_df_fin = coding_df_fin.append(coding_df)
                
                #cod_feature_index.extend(feature.qualifiers['locus_tag'] for i in range(len(coding_df)))
                #feature_index.extend(feature.qualifiers['locus_tag'] for i in range(len(noncoding_df)))
                
                cod_feature_index.extend(feature.qualifiers['db_xref'][0] for i in range(len(coding_df)))
                feature_index.extend(feature.qualifiers['db_xref'][0] for i in range(len(noncoding_df)))
                 
    noncoding_df_fin['locus_tag'] = feature_index    
    coding_df_fin['locus_tag'] = cod_feature_index
    
    #TEMP TEMP TEMP TEMP
    #noncoding_df_fin.to_csv('noncoding_df_fin.csv')
    #coding_df_fin.to_csv('coding_df_fin.csv')
    
    print('done')
    
    return noncoding_df_fin, coding_df_fin


######



def seed_comparison(all_pams,strand1df,strand2df):
    # Identifies the uniqueness or "seed score" of each guide 
    
    
    print('seed comparison')
    data = all_pams.drop_duplicates(subset='Spacer',keep='first')
    guidesA = strand1df
    guidesB = strand2df
    # guides1 = strand1df.drop_duplicates(subset='Spacer',keep='first')
    # guides2 = strand2df.drop_duplicates(subset='Spacer',keep='first')

    temp = pandas.concat([guidesA, guidesB], keys=[0, 1])
    temp2 = temp.drop_duplicates(subset='Spacer', keep=False)
    guides1, guides2 = temp2.xs(0), temp2.xs(1)
    guides1 = guides1.set_index(guides1['index'])
    guides2 = guides2.set_index(guides2['index'])

    #Check that concatenation is correct.
    #guidesA.to_csv('guidesA.csv', index=False)
    #guides1.to_csv('guides1.csv', index=False)

    data = data.assign(score='na')
    spacers = dict(zip(data.PAM_index,data.Spacer))
    
    seed_lengths = [i for i in range(4,21)] #10,11,12, etc
    seed_lib = {}
    
    start = time.perf_counter()   
    
    #creates a dictionary of dictionaries, each with spacers truncated to 
    #a certain seed truncation
    new_lib = {}
    for seed in seed_lengths:
        lib = 'res{0}'.format(seed)
        seed_lib[lib] = {} 
        new_lib[lib] = {}
        for key in spacers:
            seed_lib[lib][key] = spacers[key][-seed:]
        
        ## Counts instances of each truncated spacer
        ## If it only occurs once, it is unique within
        ## that cutoff.
        cnt = Counter()
        for idx in seed_lib[lib].values():
            cnt.update([idx])
        for pamindex,seq in seed_lib[lib].items():
            if cnt[seq] == 1:
                new_lib[lib][pamindex] = seq
        
    
       
    finish = time.perf_counter()
    print(f'It took {finish-start} seconds')
    
    scores = {}
    seed_lib = dict(reversed(list(seed_lib.items())))
    
    
    ## Starting from the shortest truncation, assign 
    ## a uniqueness score to each guide in a particular
    ## cutoff. If it is already in the dictionary, ignore.
    p = 0
    for dictionary in new_lib:
        trunc_size = seed_lengths[p]
        for key, value in new_lib[dictionary].items():
            if key in scores:
                pass
            else:
                scores[key] = trunc_size
        p += 1
        
        
    scoredf = pandas.DataFrame.from_dict(scores,orient='index')
    scoredf.reset_index(inplace=True)
    scoredf.columns = ['PAM_index','score']
    final1 = pandas.merge(guides1,scoredf,on='PAM_index',how='left')
    final2 = pandas.merge(guides2,scoredf,on='PAM_index',how='left')
        
    return final1,final2
    
def fix_format(row):
    locus = str(row.locus_tag)
    return locus
    

def top_guides(guides_per_locus,noncoding_scores, coding_scores):
    ## Picks the top guides for each locus (# = guides_per_locus)
    
    print('top guide selection')
    
    data = noncoding_scores
    
    data2 = coding_scores
    
    
    #removes guides that target two genes
    #double_hits = data[data['locus_tag'].apply(lambda x: len(str(x)) > 20)].index
    #data = data.drop(index=double_hits)
    
    #double_hits2 = data2[data2['locus_tag'].apply(lambda x: len(str(x)) > 20)].index
    #data2 = data2.drop(index=double_hits2)
    
    
    loci = list(set(data['locus_tag']))
    library_list = pandas.DataFrame(data=None,columns=['index','PAM_index','PAM_type','Component','Strand','Spacer','score','locus_tag'])
    
    for idx, locus in enumerate(loci):
        gene_guides = data[data['locus_tag'] == locus]
        gene_guides = gene_guides.sort_values(by='score')
        gene_guides['target'] = 'coding'
        top_picks = gene_guides.head(guides_per_locus)
        library_list = pandas.concat([top_picks,library_list],sort=True)
        
        
        if len(top_picks) != guides_per_locus:
            coding_guides = data2[data2['locus_tag'] == locus]
            coding_guides = coding_guides.sort_values(by='score')
            remain = guides_per_locus - len(top_picks)
            fill_picks = coding_guides.head(remain)
            fill_picks['target'] = 'noncoding'
            library_list = pandas.concat([fill_picks,library_list],sort=True)
            
    return library_list
    


########
    



#######

def run(filename):
    start = time.perf_counter()
    print('NGG PAMs')
    NGG_pams                    = create_list('NGG',filename)    
    print('NAG PAMs now')
    NAG_pams                    = create_list('NAG',filename)
    all_pams                    = pandas.concat([NGG_pams,NAG_pams]) 
    noncoding_df,coding_df      = gene_targeting_filter(all_pams,filename)
    
    finish = time.perf_counter()
    print(f'It took {finish-start} seconds')
    start = time.perf_counter()      
           
    
    noncoding_scores            = seed_comparison(all_pams,noncoding_df)
    
    finish = time.perf_counter()
    print(f'It took {finish-start} seconds')
    start = time.perf_counter()   
    
    coding_scores               = seed_comparison(all_pams,coding_df)
    
    finish = time.perf_counter()
    print(f'It took {finish-start} seconds')
    start = time.perf_counter()   
    
    noncoding_list,coding_list  = top_guides(noncoding_scores,coding_scores)

    
    finish = time.perf_counter()
    print(f'It took {finish-start} seconds')







filename = 'RAST_annotated_GB03.gbk'
start = time.perf_counter()
print('NGG PAMs')
NGG_pams                    = create_list('NGG',filename)    
print('NAG PAMs now')
NAG_pams                    = create_list('NAG',filename)
all_pams                    = pandas.concat([NGG_pams,NAG_pams])

       
noncoding_df,coding_df      = gene_targeting_filter(all_pams,filename)

finish = time.perf_counter()
print(f'It took {finish-start} seconds')


non_scores, cod_scores      = seed_comparison(all_pams,noncoding_df,coding_df)
non_scores['locus_tag']     =    non_scores.apply(fix_format,axis=1)
cod_scores['locus_tag']     =    cod_scores.apply(fix_format,axis=1)

library_list = top_guides(10,non_scores,cod_scores)
library_list.to_csv('GBO3_targeting.csv')
#run(filename)



            

