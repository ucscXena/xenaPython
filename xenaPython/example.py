import xenaPython as xena

# use case 1
# query six samples and three identifers expression
hub = "https://toil.xenahubs.net"
dataset = "tcga_RSEM_gene_tpm"
samples = ["TCGA-02-0047-01","TCGA-02-0055-01","TCGA-02-2483-01","TCGA-02-2485-01","TCGA-04-1331-01","TCGA-04-1332-01"]
probes = ["ENSG00000000003.14", "ENSG00000000005.5", "ENSG00000000419.12"]
print(xena.xenaAPI.Probes_values (hub, dataset, samples, probes))

# use case 2
# query six samples and three genes expression, when the dataset you want to query has a identifier-to-gene mapping
hub = "https://toil.xenahubs.net"
dataset = "tcga_RSEM_gene_tpm"
samples = ["TCGA-02-0047-01","TCGA-02-0055-01","TCGA-02-2483-01","TCGA-02-2485-01","TCGA-04-1331-01","TCGA-04-1332-01"]
genes =["TP53", "RB1", "PIK3CA"]
print(xena.xenaAPI.Genes_values (hub, dataset, samples, genes))

# use case 3
# if the dataset does not have id-to-gene mapping, but the dataset used gene names as its identifier,
# you can query gene expression like use case 1, use case 2 will not work
hub = "https://toil.xenahubs.net"
dataset = "tcga_RSEM_Hugo_norm_count"
samples = ["TCGA-02-0047-01","TCGA-02-0055-01","TCGA-02-2483-01","TCGA-02-2485-01","TCGA-04-1331-01","TCGA-04-1332-01"]
probes =["TP53", "RB1", "PIK3CA"]
print(xena.xenaAPI.Probes_values (hub, dataset, samples, probes))

# use case 4
# find out the samples in a dataset
hub = "https://tcga.xenahubs.net"
dataset = "TCGA.BLCA.sampleMap/HiSeqV2"
print(xena.xenaAPI.dataset_samples (hub, dataset))

# use case 5
# find out the identifiers in a dataset
hub = "https://tcga.xenahubs.net"
dataset = "TCGA.BLCA.sampleMap/HiSeqV2"
print(xena.xenaAPI.dataset_fields (hub, dataset))

# find out hub id, dataset id
# use xena browser datasets tool:  https://xenabrowser.net/datapages/
