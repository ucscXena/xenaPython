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

def anndataMatrixToTsv(adata, matFname, transpose = True):
    """
    write adata expression matrix to .tsv file"
    """
    import pandas as pd
    import scipy.sparse

    mat = adata.X
    var = adata.var
    obs = adata.obs

    # Transposing matrix, has the samples on the rows: scanpy
    # Do not transpose, has the cells on the rows: starfish
    if (transpose):
        mat = mat.T
    if scipy.sparse.issparse(mat):
        mat = mat.tocsr() # makes writing to a file ten times faster, thanks Alex Wolf!

    ofh = open(matFname, "w")

    if (transpose):
        sampleNames = obs.index.tolist()
    else:
        sampleNames = var.index.tolist()
    ofh.write("gene\t")
    ofh.write("\t".join(sampleNames))
    ofh.write("\n")

    if (transpose):
        genes = var.index.tolist()
    else:
        genes = obs.genes.tolist()
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

def adataToXena(adata, path, studyName, transpose = True):
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
            anndataMatrixToTsv(adata, matName, transpose)
    else:
        anndataMatrixToTsv(adata, matName, transpose)
    
    # build expression data .json file
    buildsjson_scRNA_geneExp(matName, studyName)

    # build meta data (phenotype data) file
    metafile = 'meta.tsv'
    metaName = join(path, metafile)
    if (transpose):
        adata.obs.to_csv(metaName, sep='\t')
    else:
        adata.var.to_csv(metaName, sep='\t')

    # build meta data .json file
    buildsjson_phenotype(metaName, studyName)

    # spatial coordinate
    if adata.obsm and 'X_spatial' in adata.obsm:
        import numpy, pandas
        data = pandas.DataFrame(adata.obsm["X_spatial"], columns=['X', 'Y'])
        data = data.set_index(adata.obs.index)
        spatial_coord_file = 'spatial.tsv'
        label = "spatial XY coordinate"

        data.to_csv(join(path, spatial_coord_file), sep='\t')
        buildsjson_phenotype(join(path, spatial_coord_file), studyName, label)

def scanpyLoomToXena(matrixFname, path, studyName, transpose = True):
    """
    Given a scanpy loom file, write dataset to a dataset directory under path.
    Transposing matrix needed, as scanpy has the samples on the rows
    """
    loomToXena(matrixFname, path, studyName, transpose)

def starfishLoomToXena(matrixFname, path, studyName, transpose = False):
    """
    Given a starfish loom file, write dataset to a dataset directory under path.
    Transposing matrix not needed, as starfish has the cells on the rows
    """
    loomToXena(matrixFname, path, studyName, transpose)

def loomToXena(matrixFname, path, studyName, transpose = True):
    """
    Given a loom file, write dataset to a dataset directory under path.
    """
    adata = sc.read(matrixFname, first_column_names=True)
    adataToXena(adata, path, studyName, transpose)

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

