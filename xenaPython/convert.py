from os.path import join, isfile, isdir
import os, sys
import datetime, json
import scanpy as sc
import pandas as pd

def dim_name(mapName, dim):
    return mapName + '_' + str(dim+1)

def buildsjson_scRNA_geneExp(output, cohort, label = None, metaPara = None):
    fout = open(output +'.json', 'w')
    J = {}
    J['type'] ='genomicMatrix'
    J['dataSubtype'] = 'gene expression'
    if label:
        J['label'] = label
    else:
        J['label'] = os.path.basename(output)
    if metaPara:
        J.update(metaPara)
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

def buildsjson_map (output, map_type, map_meta, cohort, label = None):
    fout = open(output +'.json', 'w')
    J = {}
    J['type'] ='clinicalMatrix'
    J['dataSubtype'] = map_type
    if label:
        J['label'] = label
    else:
        J['label'] = os.path.basename(output)
    J['cohort'] = cohort
    J['version'] = datetime.date.today().isoformat()

    J['map'] = map_meta
    json.dump(J, fout, indent = 4)
    fout.close()

def buildsjson_cluster(output, cluster_meta, cohort, label = None):
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
    J['cellCluster'] = cluster_meta
    json.dump(J, fout, indent = 4)
    fout.close()

def anndataMatrixToTsv(adata, matFname, transpose = True, geneColumn = "var.index"):
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
        if geneColumn == "var.index":
            genes = var.index.tolist()
        else:
            genes = var[geneColumn].tolist()
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


def adataToXena(adata, path, studyName, transpose = True, metaPara = None, geneColumn = "var.index"):
    """
    Given an anndata (adata) object, write dataset to a dataset directory under path.
    """

    if not isdir(path):
        os.makedirs(path)

    # build expression data file
    expfile = 'exprMatrix.tsv'
    matName = join(path, expfile)
    if isfile(matName):
        overwrite  = input("%s already exists. Overwriting existing files? Yes or No: " % matName)
        if overwrite.upper() == "YES":
            anndataMatrixToTsv(adata, matName, transpose = transpose, geneColumn = geneColumn)
    else:
        anndataMatrixToTsv(adata, matName, transpose = transpose, geneColumn = geneColumn)
    
    # build expression data .json file
    buildsjson_scRNA_geneExp(matName, studyName, metaPara = metaPara)

    # build cell meta data (phenotype data) file, without the cluster columns
    metafile = 'meta.tsv'
    metaName = join(path, metafile)
    if (transpose):
        adata.obs.loc[:, ~adata.obs.columns.isin(['leiden', 'louvain'])].to_csv(metaName, sep='\t')
    else:
        adata.var.loc[:, ~adata.var.columns.isin(['leiden', 'louvain'])].to_csv(metaName, sep='\t')

    # build cell meta data .json file
    buildsjson_phenotype(metaName, studyName, label="cell metadata")

    # build maps and the associated matadata
    adataToMap(adata, path, studyName)

    # build cluster and associated metadata
    assayDataset = expfile
    adataToCluster(adata, path, studyName, assayDataset)  

def adataToMap(adata, path, studyName):
    # tsne, umap, spatial coordinates
    if adata.obsm is not None:
        import numpy

        for map in adata.obsm.keys():
            cols =[]
            if map == 'X_umap':
                mapName = "umap"
                dataSubType = 'embedding'
                label = 'umap'
                map_file = 'umap.tsv'
            elif map == 'X_tsne':
                mapName = "tsne"
                dataSubType = 'embedding'
                label = 'tsne'
                map_file =  'tsne.tsv'
            elif map == 'X_spatial':
                mapName = 'spatial_map'
                dataSubType = 'spatial'
                label = 'spatial map'
                map_file =  'spatial_map.tsv'
            elif map == 'spatial': # visium
                mapName = 'spatial_map'
                dataSubType = 'spatial'
                label = 'spatial map'
                map_file =  'spatial_map.tsv'
            else:
                print("unrecognized or ignored map:", map)
                continue

            row,col = adata.obsm[map].shape
            col = min(col, 3)

            for i in range (0, col):
                colName = dim_name(mapName, i)
                cols.append(colName)

            df = pd.DataFrame(adata.obsm[map][:,range(col)], columns=cols)


            df = df.set_index(adata.obs.index)
            df_meta = [{
                'label': label,
                'type': dataSubType,
                'dimension':cols
            }]

            df.to_csv(join(path, map_file), sep='\t')
            map_type = dataSubType
            buildsjson_map(join(path, map_file), map_type, df_meta, studyName, label)


def adataToCluster (adata, path, studyName, assayDataset):
    df = pd.DataFrame()
    cluster_file = 'cluster.tsv'
    label = 'cell clusters'
    df_meta = []

    for cluster in adata.obs.keys():
        if cluster == 'leiden':
            df['leiden'] = adata.obs['leiden']
            feature = 'leiden'
            assay = 'scanpy leiden'
            assayDataset = assayDataset
            assayDatasetLabel = 'scRNA-seq'
            assayParameter = 'default parameter'

        elif cluster == 'louvain':
            df['louvain'] = adata.obs['louvain']
            feature = 'louvain'
            assay = 'scanpy louvain'
            assayDataset = assayDataset
            assayDatasetLabel = 'scRNA-seq'
            assayParameter = 'default parameter'
        else:
            continue

        df_meta.append({
            'feature': feature,
            'assay': assay,
            'assayParameter': assayParameter,
            'assayDataset': {
                'host': './',
                'name': assayDataset
            },
            'assayDatasetLabel': assayDatasetLabel
        })

    if len(df.columns) >0:
        df.to_csv(join(path, cluster_file), sep='\t')
        buildsjson_cluster(join(path, cluster_file), df_meta, studyName, label)


def starfishExpressionMatrixToXena(mat, path, studyName):
    """
    Given a starfish ExpressionMatrix object (mat), write dataset to a dataset directory under path.
    """

    # build expression data file
    expfile = 'exprMatrix.tsv'
    matName = join(path, expfile)

    if isfile(matName):
        overwrite  = input("%s already exists. Overwriting existing files? Yes or No: " % matName)
        if overwrite.upper() == "YES":
            mat.to_pandas().transpose().to_csv(matName, sep='\t')
    else:
        mat.to_pandas().transpose().to_csv(matName, sep='\t')

    # build expression data .json file
    buildsjson_scRNA_geneExp(matName, studyName)

    # build meta data (phenotype data) file
    metafile = 'meta.tsv'
    metaName = join(path, metafile)

    cells = mat.cells.data.tolist()
    features = mat.cells.coords

    ofh = open(metaName, "w")
    ofh.write("\t")
    ofh.write("\t".join(features))
    ofh.write("\n")
    for i, cell in enumerate(cells):
        ofh.write(str(cell))
        for k in features:
            ofh.write("\t" + str(features[k].values[i]))
        ofh.write("\n")
    ofh.close()

    # build meta data .json file
    buildsjson_phenotype(metaName, studyName)


def scanpyLoomToXena(matrixFname, outputpath, studyName, transpose = True):
    """
    Given a scanpy loom file, write dataset to a dataset directory under path.
    Transposing matrix needed, as scanpy has the samples on the rows
    """
    loomToXena(matrixFname, outputpath, studyName, transpose = transpose)

def starfishLoomToXena(matrixFname, outputpath, studyName, transpose = False):
    """
    Given a starfish loom file, write dataset to a dataset directory under path.
    Transposing matrix not needed, as starfish has the cells on the rows
    """
    loomToXena(matrixFname, outputpath, studyName, transpose = transpose)

def loomToXena(matrixFname, outputpath, studyName, transpose = True):
    """
    Given a loom file, write dataset to a dataset directory under path.
    """
    adata = sc.read(matrixFname, first_column_names=True)
    adataToXena(adata, outputpath, studyName, transpose = transpose)

def h5adToXena(h5adFname, outputpath, studyName):
    """
    Given a h5ad file, write dataset to a dataset directory under path.
    """
    adata = sc.read_h5ad(h5adFname)
    adataToXena(adata, outputpath, studyName)

def basic_analysis(adata, normalization = True):
    # normalize_total_count (or intensity), log1p, pca, 3D umap (dense) and clustering (leiden, louvain)

    n_components = 3

    if (normalization):
        sc.pp.filter_cells(adata, min_genes=0)
        sc.pp.filter_cells(adata, min_counts=0)
        sc.pp.normalize_total(adata, inplace=True)
        sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata)

    #PCA
    sc.tl.pca(adata, svd_solver='arpack')

    # UMAP 3D
    import umap
    # run umap in dense mode  https://www.nature.com/articles/s41587-020-00801-7
    dens_lambda = 1 # default = 2
    embedding = umap.UMAP(densmap=True, n_components = n_components, dens_lambda= dens_lambda).fit(adata.obsm['X_pca'])
    adata.obsm['X_umap'] = embedding.embedding_

    # clustering
    sc.pp.neighbors(adata, n_pcs=20)
    sc.tl.leiden(adata)
    sc.tl.louvain(adata)
    return adata

def visiumToXena(visiumDataDir, count_file, outputpath, studyName):
    """
    Given a visium spaceranger output data directory, write dataset to a dataset directory under path.
    """
    # https://scanpy.readthedocs.io/en/stable/api/scanpy.read_visium.html

    if not (count_file.endswith("filtered_feature_bc_matrix.h5")):
        print (count_file, "is not the correct format")
        sys.exit()
    else:
        print (count_file)

    adata = sc.read_visium(visiumDataDir, count_file = count_file)
    
    adata = basic_analysis(adata)

    metaPara = {}
    metaPara['unit'] = "LogNorm(count+1)"
    metaPara['wrangling_procedure'] = "download filtered_feature_bc_matrix.h5, normalize count data using scanpy sc.pp.normalize_total(adata), then sc.pp.log1p(adata)"
    adataToXena(adata, outputpath, studyName, metaPara = metaPara)

def vizgenToXena(vizgenDataDir, outputpath, studyName):
    """
    Given a vizgen output data directory, write dataset to a dataset directory under path.
    """
    # https://f.hubspotusercontent40.net/hubfs/9150442/Vizgen%20MERFISH%20Mouse%20Receptor%20Map%20File%20Descriptions%20.pdf?__hstc=30510752.65b077e2f6b41ba4f2e0c44a2103598e.1631299341334.1631299341334.1631299341334.1&__hssc=30510752.3.1631299341334&__hsfp=3105977984&hsCtaTracking=f0a4edb5-afb5-4b5c-b3fe-5b73da111821%7Ce87c6069-24f9-4538-a9a1-e54304c082b2

    for file in os.listdir(vizgenDataDir):
        import re
        exp_pattern = 'cell_by_gene.*csv$'
        meta_pattern = 'cell_metadata.*csv$'
        if re.search(exp_pattern, file):
            count_file = file
            print (count_file)
        if re.search(meta_pattern, file):
            meta_file = file
            print(meta_file)

    adata = sc.read_csv(count_file, first_column_names = True)
    meta_cell = pd.read_csv(meta_file, index_col=0)
    adata.obs = meta_cell

    adata = basic_analysis(adata)

    metaPara = {}
    metaPara['unit'] = "log(count+1)"
    metaPara['wrangling_procedure'] = "download cell_by_gene.csv, normalize count data using scanpy sc.pp.normalize_total(adata), then sc.pp.log1p(adata)"
    adataToXena(adata, outputpath, studyName, metaPara = metaPara)

def sprmToXena(sprmDataDir, outputpath, studyName):
    """
    Given a SPRM output data directory, write dataset to a dataset directory under path.
    """
    # https://github.com/hubmapconsortium/sprm
    # https://github.com/hubmapconsortium/codex-pipeline
    # https://view.commonwl.org/workflows/github.com/hubmapconsortium/codex-pipeline/blob/f3d6e97408b1c542641b313c1ea8d3115d72e3f8/pipeline.cwl

    for file in os.listdir(sprmDataDir):
        import re
        exp_pattern = 'cell_channel_mean.csv$'
        meta_pattern = 'cell_centers.csv$'
        if re.search(exp_pattern, file):
            count_file = file
            print (count_file)
        if re.search(meta_pattern, file):
            meta_file = file
            print(meta_file)

    import pandas as pd
    adata = sc.read_csv(count_file, first_column_names = True)
    meta_cell = pd.read_csv(meta_file, index_col=0, names=['y','x'])  # the hubmap current output from sprm (2021-09) might have the header x and y swapped
    meta_cell.index = meta_cell.index.astype(str)
    meta_cell = meta_cell.filter(items = list(adata.obs_names), axis=0)
    adata.obs = meta_cell

    adata = basic_analysis(adata)

    metaPara = {}
    metaPara['dataSubtype'] = 'protein expression'
    metaPara['unit'] = "log(intensity+1)"
    metaPara['wrangling_procedure'] = "download cell_channel_mean.csv, normalize cell mean intensity data using scanpy sc.pp.normalize_total(adata), then sc.pp.log1p(adata)"
    adataToXena(adata, outputpath, studyName, metaPara = metaPara)
