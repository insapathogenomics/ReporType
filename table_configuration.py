# -*- coding: utf-8 -*-

import pandas as pd
import os
import csv
import numpy as np
import sys


pd.options.mode.chained_assignment = None  

def results_path(path):
    FILES=[]
    for p,_,files in os.walk(os.path.abspath(path)):
        for file in files:
            file=os.path.join(p, file)
            FILES.append (file)
    return (FILES)


def types (list,coverage,id,gene,acce):
    SUB_I=[]
    id_type,cov_type,gene_type,acce_type,id_sub,cov_sub,gene_sub,acce_sub=(None,None,None,None,None,None,None,None)
    for i in range(len(list)):
        if list[i] == 'species' or list[i] == 'type' :
            id_type=gene[i]+'-'+str(id[i])
            cov_type=gene[i]+'-'+str(coverage[i])
            acce_type=gene[i]+'-'+str(acce[i])
            gene_type=gene[i]
        if list[i] == 'subtype' or list[i] == 'genus' or list[i] == 'lineage':
            SUB_I.append(int(i))
            ID=[]
            COV=[]
            GENE=[]
            ACCE=[]
            for x in SUB_I:
                id_sueb=str(gene[x])+'-'+str(id[x])
                cov_sueb=str(gene[x])+'-'+str(coverage[x])
                acce_sueb=str(gene[x])+'-'+str(acce[x])
                gene_sueb=str(gene[x])
                ID.append(id_sueb)
                COV.append(cov_sueb)
                GENE.append(gene_sueb)
                ACCE.append(acce_sueb)
            GENE=sorted(GENE)
            ID=sorted(ID)
            COV=sorted(COV)
            ACCE=sorted(ACCE)
            join=' : '
            id_sub=join.join(ID)
            cov_sub=join.join(COV)
            gene_sub=join.join(GENE)
            acce_sub=join.join(ACCE) 
    return id_type,cov_type,gene_type,acce_type,id_sub,cov_sub,gene_sub,acce_sub

def check_only_type(database,TYPES):
    
    for i in range(len(TYPES)):
        if database == TYPES[i]:
            unique=True
        else:
            unique=False
    if database=="ompA" or database=="refseqid" or database=="Lp" or database=="HCV" or database=="HTVL-1" or database=="lp":
        unique=True
    if database=="influenza" or database=="coronavirus":
        unique=False
    return(unique)

def check_multi(multi,file_in_name):
    if file_in_name in multi:
        mult=True
    elif "all" in multi:
        mult=True
    else:
        mult=False
    return (mult) 

FILES=sys.argv[1:-4]


if len(FILES)==0:
    print("  ERROR: Some error occured during the analysis, verify your samples  ")

args=sys.argv[1:]


file_name = sys.argv[-4]

database_name=sys.argv[-3]

multi=sys.argv[-2]

sort=sys.argv[-1]


def transform_in_list(obj):
   LIST=[]
   for line in obj:
       line=str(line)
       LIST.append(line)
   return(LIST)

for file in FILES:
    file_in_name=str(file)
    file_in_name=file_in_name.strip()
    file_in_name=file_in_name.split("/")[-1]
    file_in_name=file_in_name.replace(".tab","")
    mult = check_multi(multi,file_in_name)
    file = pd.read_table(file,header=0,usecols=['#FILE','SEQUENCE','GENE','%COVERAGE','%IDENTITY','DATABASE','ACCESSION'])
    n=len(file)
    if n == 0:
        file=file.drop(columns=["SEQUENCE"])
        file= pd.DataFrame({'#FILE':[file_in_name],'GENE':["0 genes found (check if there is any warning during the analysis)"],'%COVERAGE':["---"],'%IDENTITY':["---"],'DATABASE':["---"],'ACCESSION':["---"]})

        file.to_csv(file_name,sep="\t",header=False, index=False,mode='a')

    else:
        file=file.sort_values(by=[sort], ascending=False)

        if mult == True:
            file=file.drop_duplicates(subset=["GENE","SEQUENCE"])
            file["#FILE"]=file["SEQUENCE"] 
            
        else:
            file=file.drop_duplicates(subset="GENE")
            FILE=transform_in_list(file["#FILE"])
            for f in FILE:
                FILE=[]
                f=f.replace(".fasta","")
                FILE.append(f)
                n=len(file["#FILE"])
                file["#FILE"]=FILE*n
        file=file.drop(columns=["SEQUENCE"])
        type=file["DATABASE"]
        gene_1=file["GENE"]
        database_1=file["DATABASE"]
        
        TYPES=[]
        DBS=[]
        for line in type:
            line=str(line)
            line.strip()
            database=line.split('_')[0]
            type=line.split('_')[-1]
            if type !='DATABASE':
                TYPES.append(type)
            if database!='DATABASE':
                DBS.append(database)
        database=DBS[0]
        n=len(database_1)
        unique= check_only_type(database,TYPES)
        for index, row in file.iterrows():
            if index < len(gene_1):
                if row["DATABASE"] == gene_1[index]:
                     
                     file["DATABASE"] = [database_name]*n
            file.to_csv(file_name,sep="\t",header=False, index=False,mode='a')
            
            break
      
        else:
        
            if unique == True:
                file.to_csv(file_name,sep="\t",header=False, index=False,mode='a')


            else:
                grup=file.groupby("#FILE")
                for f_file, columns in grup:
                    type=columns["DATABASE"]
                    TYPES=[]
                    DBS=[]
                    for line in type:
                        line=line.strip()
                        database=line.split('_')[0]
                        type=line.split('_')[-1]
                        if type !='DATABASE':
                            TYPES.append(type)
                        if database!='DATABASE':
                            DBS.append(database)
                    database=DBS[0]
                    n=len(columns)
                    columns["DATABASE"]=[database]*n
                    coverage=columns["%COVERAGE"]
                    coverage=transform_in_list(coverage)
                    id=columns['%IDENTITY']
                    id=transform_in_list(id)
                    gene=columns['GENE']
                    gene=transform_in_list(gene)
                    acce=columns['ACCESSION']
                    acce=transform_in_list(acce)
                    id_type,cov_type,gene_type,acce_type,id_sub,cov_sub,gene_sub,acce_sub=types(TYPES,coverage,id,gene,acce)
                    if gene_type is None:
                        columns['%IDENTITY']=[str(id_sub)]*n
                        columns['GENE']=[str(gene_sub)]*n
                        columns['%COVERAGE']=[str(cov_sub)]*n
                        columns['ACCESSION']=[str(acce_sub)]*n
                    elif gene_sub is None:
                        columns['%IDENTITY']=[str(id_type)]*n
                        columns['GENE']=[str(gene_type)]*n
                        columns['%COVERAGE']=[str(cov_type)]*n
                        columns['ACCESSION']=[str(acce_type)]*n   
                    else:
                        if database =="influenza":
                            ID_SUB_B=[]
                            ID_SUB_A=[]
                            COV_SUB_B=[]
                            COV_SUB_A=[]
                            ACCE_SUB_B=[]
                            ACCE_SUB_A=[]
                            GENE_SUB_A=[]
                            GENE_SUB_B=[]
                            cov_sub=cov_sub.split(':')
                            gene_sub=gene_sub.split(':')
                            id_sub=id_sub.split(':')
                            acce_sub=acce_sub.split(':')
                            for id in id_sub:
                                if "Yamagata" in id or "Victoria" in id:
                                    ID_SUB_B.append(id)
                                if "H" in id or "N" in id:
                                    ID_SUB_A.append(id)                            
                            id_sub_b=':'.join(ID_SUB_B).replace(" ","")
                            id_sub_a=':'.join(ID_SUB_A).replace(" ","")
                            for id in cov_sub:
                                if "Yamagata" in id or "Victoria" in id:
                                    COV_SUB_B.append(id)
                                if "H" in id or "N" in id:
                                    COV_SUB_A.append(id)                            
                            cov_sub_b=':'.join(ID_SUB_B).replace(" ","")
                            cov_sub_a=':'.join(ID_SUB_A).replace(" ","")
                            for id in acce_sub:
                                if "Yamagata" in id or "Victoria" in id:
                                    ACCE_SUB_B.append(id)
                                if "H" in id or "N" in id:
                                    ACCE_SUB_A.append(id)                            
                            acce_sub_b=':'.join(ACCE_SUB_B).replace(" ","")
                            acce_sub_a=':'.join(ACCE_SUB_A).replace(" ","")
                            for id in gene_sub:
                                if "Yamagata" in id or "Victoria" in id:
                                    GENE_SUB_B.append(id)
                                if "H" in id or "N" in id:
                                    GENE_SUB_A.append(id)                            
                            gene_sub_b=', other hits:'.join(GENE_SUB_B).replace(" ","")
                            gene_sub_a=''.join(GENE_SUB_A).replace(" ","")
                            if gene_type=="A":
                                columns['%IDENTITY']=[str(id_type)+':'+str(id_sub_a)]*n
                                columns['GENE']=[str(gene_type)+'-'+str(gene_sub_a)]*n
                                columns['%COVERAGE']=[str(cov_type)+':'+str(cov_sub_a)]*n
                                columns['ACCESSION']=[str(acce_type)+':'+str(acce_sub_a)]*n 
                                
                            if gene_type=="B":
                                columns['%IDENTITY']=[str(id_type)+':'+str(id_sub_b)]*n
                                columns['GENE']=[str(gene_type)+'-'+str(gene_sub_b)]*n
                                columns['%COVERAGE']=[str(cov_type)+':'+str(cov_sub_b)]*n
                                columns['ACCESSION']=[str(acce_type)+':'+str(acce_sub_b)]*n 
                        
                            
                        else:
                            join=" : "
                            columns['%IDENTITY']=[str(id_type)+str(join)+str(id_sub)]*n
                            columns['GENE']=[str(gene_type)+'-'+str(gene_sub)]*n
                            columns['%COVERAGE']=[str(cov_type)+str(join)+str(cov_sub)]*n
                            columns['ACCESSION']=[str(acce_type)+str(join)+str(acce_sub)]*n 
                    file=columns.drop_duplicates(subset="#FILE")
                
                    file.to_csv(file_name,sep="\t",header=False, index=False,mode='a')
                
                



final_file=pd.read_csv(file_name, sep='\t',names=['SAMPLE','HIT','COVERAGE(%)','IDENTITY(%)','DATABASE','ACCESSION'])


if final_file['DATABASE'][0] != "influenza":

    final_file = final_file.groupby('SAMPLE').agg(lambda x: ','.join(x.astype(str))).reset_index()

    duplicated = final_file.duplicated('SAMPLE')




    final_file = final_file.drop_duplicates(subset='SAMPLE')
    

    final_file['COVERAGE(%)'] = final_file.apply(lambda row: f"{row['HIT'].split(',')[0]}-{row['COVERAGE(%)'].split(',')[0]}; other hits: {', '.join(f'{a}-{b}' for a, b in zip(row['HIT'].split(',')[1:], row['COVERAGE(%)'].split(',')[1:]))}" if len(row['HIT'].split(',')) > 1 else f"{row['HIT']}-{row['COVERAGE(%)']}", axis=1)
    

    final_file['ACCESSION'] = final_file.apply(lambda row: f"{row['HIT'].split(',')[0]}-{row['ACCESSION'].split(',')[0]}; other hits: {', '.join(f'{a}-{b}' for a, b in zip(row['HIT'].split(',')[1:], row['ACCESSION'].split(',')[1:]))}" if len(row['HIT'].split(',')) > 1 else f"{row['HIT']}-{row['ACCESSION']}", axis=1)
    
    
    final_file['IDENTITY(%)'] = final_file.apply(lambda row: f"{row['HIT'].split(',')[0]}-{row['IDENTITY(%)'].split(',')[0]}; other hits: {', '.join(f'{a}-{b}' for a, b in zip(row['HIT'].split(',')[1:], row['IDENTITY(%)'].split(',')[1:]))}" if len(row['HIT'].split(',')) > 1 else f"{row['HIT']}-{row['IDENTITY(%)']}", axis=1)
    
    

    final_file['HIT'] = final_file.apply(lambda row: f"{str(row['HIT']).split(',')[0]}; other hits: {', '.join(str(x) for x in str(row['HIT']).split(',')[1:])}" if len(str(row['HIT']).split(',')) > 1 else str(row['HIT']), axis=1)

    final_file['DATABASE']=final_file.apply(lambda row:f"{row['DATABASE'].split(',')[0]}", axis=1)
    
    

    mask = ~final_file['HIT'].str.contains('other hits')

    cov_values = final_file.loc[mask, 'COVERAGE(%)'].str.split('-', expand=True)

    final_file.loc[mask, 'COVERAGE(%)'] = cov_values.apply(lambda row: row[row.last_valid_index()], axis=1)



    id_values = final_file.loc[mask, 'IDENTITY(%)'].str.split('-', expand=True)
    final_file.loc[mask, 'IDENTITY(%)'] = id_values.apply(lambda row: row[row.last_valid_index()], axis=1)
    



final_file.to_csv(file_name,header=True, index=False,mode='w', sep='\t')