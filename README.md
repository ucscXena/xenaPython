# xenaPython
Python API for Xena Hub

---------

#### Requirement
    support python2


#### Installation
    pip install xenaPython


#### Upgrade Installation
    pip install --upgrade xenaPython


#### Usage

    >>import xenaPython as xena
    >>xena.xenaAPI


#### Examples

##### 1: Query four samples and three identifers expression
    import xenaPython as xena

    hub = "https://toil.xenahubs.net"
    dataset = "tcga_RSEM_gene_tpm"
    samples = ["TCGA-02-0047-01","TCGA-02-0055-01","TCGA-02-2483-01","TCGA-02-2485-01"]
    probes = ['ENSG00000282740.1', 'ENSG00000000005.5', 'ENSG00000000419.12']
    print xena.xenaAPI.Probes_values (hub, dataset, samples, probes)

##### 2: Query four samples and three genes expression, when the dataset you want to query has a identifier-to-gene mapping
    hub = "https://toil.xenahubs.net"
    dataset = "tcga_RSEM_gene_tpm"
    samples = ["TCGA-02-0047-01","TCGA-02-0055-01","TCGA-02-2483-01","TCGA-02-2485-01"]
    genes =["TP53", "RB1", "PIK3CA"]
    print xena.xenaAPI.Genes_values (hub, dataset, samples, genes)

##### 3: If the dataset does not have id-to-gene mapping, but the dataset used gene names as its identifier, you can query gene expression like example 1, example 2 will not work.
    hub = "https://toil.xenahubs.net"
    dataset = "tcga_RSEM_Hugo_norm_count"
    samples = ["TCGA-02-0047-01","TCGA-02-0055-01","TCGA-02-2483-01","TCGA-02-2485-01"]
    probes =["TP53", "RB1", "PIK3CA"]
    print xena.xenaAPI.Probes_values (hub, dataset, samples, probes)

##### 4: Find out the samples in a dataset
    hub = "https://tcga.xenahubs.net"
    dataset = "TCGA.BLCA.sampleMap/HiSeqV2"
    print xena.xenaAPI.dataset_samples (hub, dataset)

##### 5: Find out the identifiers in a dataset
    hub = "https://tcga.xenahubs.net"
    dataset = "TCGA.BLCA.sampleMap/HiSeqV2"
    print xena.xenaAPI.dataset_fields (hub, dataset)

##### 6. Find out hub id, dataset id
    use xena browser datasets tool:  https://xenabrowser.net/datapages/



#### Contact
     http://xena.ucsc.edu/
     https://groups.google.com/forum/#!forum/ucsc-cancer-genomics-browser
     genome-cancer@soe.ucsc.edu


#### Acknowledgement
[@GiriB](https://github.com/GiriB)

