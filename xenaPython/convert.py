from os.path import join, isfile, isdir
import os, sys
import datetime, json
import scanpy as sc
import cellbrowser.cellbrowser as cb

def buildsjson_scRNA_geneExp(output, cohort, label = None):
    fout = open(output +'.json', 'w')
    J = {}
    J['type'] ='genomicMatrix'
    J['dataSubtype'] = 'scRNA-seq gene expression'
    if label:
        J['label'] = label
    else:
        J['label'] = os.path.basename(output)
    J["colNormalization"] = True
    J['cohort'] = cohort
    J['version'] = datetime.date.today().isoformat()
    json.dump(J, fout, indent = 4)
    fout.close()

def buildsjson_scRNA_meta(output, cohort, label = None):
    fout = open(output +'.json', 'w')
    J = {}
    J['type'] ='clinicalMatrix'
    J['dataSubtype'] = 'phenotype'
    if label:
        J['label'] = label
    else:
        J['label'] = os.path.basename(output)
    J['cohort'] = cohort
    J['version'] = datetime.date.today().isoformat()
    json.dump(J, fout, indent = 4)
    fout.close()

def adataToXena(adata, path, studyName):
    """
    Given a scanpy object, write dataset to a dataset directory under path.
    """

    if not isdir(path):
        os.mkdir(path)

    # build expression data file
    expfile = 'exprMatrix.tsv'
    matName = join(path, expfile)
    if isfile(matName):
        overwrite  = input("%s already exists. Overwriting existing files? Yes or No: " % matName)
        if overwrite.upper() == "YES":
            cb.anndataMatrixToTsv(adata, matName)
    else:
        cb.anndataMatrixToTsv(adata, matName)
    
    # build expression data .json file
    buildsjson_scRNA_geneExp(matName, studyName)

    # build meta data (phenotype data) file
    metafile = 'meta.tsv'
    metaName = join(path, metafile)
    adata.obs.to_csv(metaName, sep='\t')

    # build meta data .json file
    buildsjson_scRNA_meta(metaName, studyName)

    
def loomToXena(matrixFname, path, studyName):
    """
    Given a loom file, write dataset to a dataset directory under path.
    """
    adata = sc.read(matrixFname, first_column_names=True)
    adataToXena(adata, path, studyName)

def h5adToXena(matrixFname, path, studyName):
    """
    Given a h5ad file, write dataset to a dataset directory under path.
    """
    adata = sc.read(matrixFname, first_column_names=True)
    adataToXena(adata, path, studyName)
