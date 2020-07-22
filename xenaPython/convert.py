from os.path import join, isfile, isdir
import os, sys
import datetime, json
import scanpy as sc

def buildsjson_scRNA_geneExp(output, cohort, label = None):
    fout = open(output +'.json', 'w')
    J = {}
    J['type'] ='genomicMatrix'
    J['dataSubtype'] = 'gene expression'
    if label:
        J['label'] = label
    else:
        J['label'] = os.path.basename(output)
    J["colNormalization"] = True
    J['cohort'] = cohort
    J['version'] = datetime.date.today().isoformat()
    json.dump(J, fout, indent = 4)
    fout.close()

def buildsjson_phenotype(output, cohort, label = None):
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

def anndataMatrixToTsv(adata, matFname):
    """
    write adata expression matrix to .tsv file"
    """
    import pandas as pd
    import scipy.sparse

    mat = adata.X
    var = adata.var

    # Transposing matrix, as scanpy has the samples on the rows
    mat = mat.T

    if scipy.sparse.issparse(mat):
        mat = mat.tocsr() # makes writing to a file ten times faster, thanks Alex Wolf!

    ofh = open(matFname, "w")

    sampleNames = adata.obs.index.tolist()
    ofh.write("gene\t")
    ofh.write("\t".join(sampleNames))
    ofh.write("\n")

    genes = var.index.tolist()

    print("Writing %d genes in total" % len(genes))
    for i, geneName in enumerate(genes):
        if i % 2000==0:
            print("Wrote %d genes" % i)
        ofh.write(geneName)
        ofh.write("\t")
        if scipy.sparse.issparse(mat):
            row = mat.getrow(i).todense()
        else:
            row = mat[i,:]

        row.tofile(ofh, sep="\t", format="%.7g")
        ofh.write("\n")

    ofh.close()

def adataToXena(adata, path, studyName):
    """
    Given a scanpy anndata (adata) object, write dataset to a dataset directory under path.
    """

    if not isdir(path):
        os.makedirs(path)

    # build expression data file
    expfile = 'exprMatrix.tsv'
    matName = join(path, expfile)
    if isfile(matName):
        overwrite  = input("%s already exists. Overwriting existing files? Yes or No: " % matName)
        if overwrite.upper() == "YES":
            anndataMatrixToTsv(adata, matName)
    else:
        anndataMatrixToTsv(adata, matName)
    
    # build expression data .json file
    buildsjson_scRNA_geneExp(matName, studyName)

    # build meta data (phenotype data) file
    metafile = 'meta.tsv'
    metaName = join(path, metafile)
    adata.obs.to_csv(metaName, sep='\t')

    # build meta data .json file
    buildsjson_phenotype(metaName, studyName)

    
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

def visiumToXena(visiumDataDir, count_file, path, studyName):
    """
    Given a visium spacerane output data directory, write dataset to a dataset directory under path.
    """
    # https://scanpy.readthedocs.io/en/stable/api/scanpy.read_visium.html
    adata = sc.read_visium(visiumDataDir, count_file=count_file)
    adataToXena(adata, path, studyName)

    import numpy, pandas
    data = pandas.DataFrame(adata.obsm["X_spatial"], columns=['X', 'Y'])
    data = data.set_index(adata.obs.index)
    spatial_coord_file = 'spatial.tsv'
    label = "spatial XY coordinate"

    data.to_csv(join(path, spatial_coord_file), sep='\t')
    buildsjson_phenotype(join(path, spatial_coord_file), studyName, label)
