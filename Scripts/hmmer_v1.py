#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 20:25:28 2019

@author: ane
"""
import os
import subprocess
from Bio import SearchIO
import json
import csv
#%%

def hmmprofile_list(hmmfile):
    """Obtains hmm profile file names from input csv file and returns them
    as a list.
    """
    f=open(hmmfile,"r")
    lines=f.readlines()
    hmmprofiles=[]
    for x in lines:
        hmmprofiles.append(x.split(',')[0])
    hmmprofiles=hmmprofiles[1:]
    return hmmprofiles

#%%
hmmprofile_list("hmm_filenames.csv")

#%%

def databases_list(dbfile):
    """Obtains database file names from input csv file and returns them
    as a list.
    """
    f=open(dbfile,"r")
    lines=f.readlines()
    databases=[]
    for x in lines:
        databases.append(x.split(',')[0])
    databases=databases[1:]
    return databases

#%%
databases_list("targets_fafiles.csv")

#%%

def call_hmmsearch(hmmsearch_args):

    """Calls hmmersearch with given parameters
    Args:
        hmmsearch_args (string): command line arguments, including query file
        and hmmer database file.
    """
    subprocess.call(hmmsearch_args)
    
#%%

def run_hmmsearch(hmmprofile, eval, database):

    hmmsearch_args = ['hmmsearch', '-o', 'output.txt', '-E', eval, hmmprofile, database]

    call_hmmsearch(hmmsearch_args)
    

#%%
def main_function(hmmfile, dbfile, json_outputfile, csv_outputfile, eval):
    POI = hmmprofile_list(hmmfile)
    dbs = databases_list(dbfile)
	
    output = {}
    with open(csv_outputfile, 'w') as myfile:
        wr = csv.writer(myfile)
        wr.writerow(['id', 'strain','accession', 'description', 'e-value'])
	
        for id in POI:
            output[id] = []
            
            for strain in dbs:            
        
                run_hmmsearch(id, eval, strain)
                results = SearchIO.parse('output.txt', 'hmmer3-text')
                
                a = {}
                a[strain] = {}
                a[strain]['accession'] = []
                a[strain]['description'] = []
                a[strain]['e-value'] = []
        
                for res in results:
                    hits = res.hits
                    for hit in hits:
                        a[strain]['accession'].append(hit.id)
                        a[strain]['description'].append(hit.description)
                        a[strain]['e-value'].append(hit.evalue)
                        
                        wr.writerow([id, strain, hit.id, hit.description, hit.evalue])
                            
                output[id].append(a)
                        
    json.dump(output, open(json_outputfile, "w"))
        
    return output

#%%
    
main_function("hmm_filenames.csv", "targets_fafiles.csv", "results_hmmer_enteric.json", "results_hmmer_enteric.csv", "0.001")
