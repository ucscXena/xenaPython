# xenaPython
Python API for Xena Hub

---------

#### Requirement
    support python2 python3


#### Installation
    pip install xenaPython


#### Upgrade Installation
    pip install --upgrade xenaPython


#### Usage
    >>> import xenaPython as Xena

#### Examples

##### 1: Query four samples and three identifers expression
    import xenaPython as xena

    hub = "https://toil.xenahubs.net"
    dataset = "tcga_RSEM_gene_tpm"
    samples = ["TCGA-02-0047-01","TCGA-02-0055-01","TCGA-02-2483-01","TCGA-02-2485-01"]
    probes = ['ENSG00000282740.1', 'ENSG00000000005.5', 'ENSG00000000419.12']
    [position, [ENSG00000282740_1, ENSG00000000005_5, ENSG00000000419_12]] = xena.dataset_probe_values(hub, dataset, samples, probes)
    ENSG00000282740_1
    
##### 2: Query four samples and three genes expression, when the dataset you want to query has a identifier-to-gene mapping (i.e. xena probeMap)
    hub = "https://toil.xenahubs.net"
    dataset = "tcga_RSEM_gene_tpm"
    samples = ["TCGA-02-0047-01","TCGA-02-0055-01","TCGA-02-2483-01","TCGA-02-2485-01"]
    genes =["TP53", "RB1", "PIK3CA"]
    xena.dataset_gene_probe_avg(hub, dataset, samples, genes)

##### 3: If the dataset does not have id-to-gene mapping, but the dataset used gene names as its identifier, you can query gene expression like example 1, example 2 will not work.
    hub = "https://toil.xenahubs.net"
    dataset = "tcga_RSEM_Hugo_norm_count"
    samples = ["TCGA-02-0047-01","TCGA-02-0055-01","TCGA-02-2483-01","TCGA-02-2485-01"]
    probes =["TP53", "RB1", "PIK3CA"]
    [position, [TP53, RB1, PIK3CA]] = xena.dataset_probe_values (hub, dataset, samples, probes)
    TP53

##### 4: Find out the samples in a dataset
    hub = "https://tcga.xenahubs.net"
    dataset = "TCGA.BLCA.sampleMap/HiSeqV2"
    xena.dataset_samples (hub, dataset, 10)
    xena.dataset_samples (hub, dataset, None)

##### 5: Find out the identifiers in a dataset
    hub = "https://tcga.xenahubs.net"
    dataset = "TCGA.BLCA.sampleMap/HiSeqV2"
    xena.dataset_field (hub, dataset)

##### 6. Find out the number of idnetifiers in a dataset
    hub = "https://tcga.xenahubs.net"
    dataset = "TCGA.BLCA.sampleMap/HiSeqV2"
    xena.dataset_field_n (hub, dataset)
    
##### 7. Find out hub id, dataset id
    use xena browser datasets tool:  https://xenabrowser.net/datapages/

#### Help
    >>> import xenaPython
    >>> help(xenaPython)
    
Help on package xenaPython:

NAME

    xenaPython - Methods for querying data from UCSC Xena hubs

DESCRIPTION

    Data rows are associated with "sample" IDs.
    Sample IDs are unique within a "cohort".
    A "dataset" is a particular assay of a cohort, e.g. gene expression.
    Datasets have associated metadata, specifying their data type and cohort.
    
    There are three primary data types: dense matrix (samples by probes),
    sparse (sample, position, variant), and segmented (sample, position, value).
    
    
    Dense matrices can be genotypic or phenotypic. Phenotypic matrices have
    associated field metadata (descriptive names, codes, etc.).
    
    Genotypic matricies may have an associated probeMap, which maps probes to
    genomic locations. If a matrix has hugo probeMap, the probes themselves
    are gene names. Otherwise, a probeMap is used to map a gene location to a
    set of probes.

FUNCTIONS

    all_cohorts(host, exclude)
    
    all_datasets(host)
    
    all_datasets_n(host)
        Count the number datasets with non-null cohort
    
    all_field_metadata(host, dataset)
        Metadata for all dataset fields (phenotypic datasets)
    
    cohort_samples(host, cohort, limit)
        All samples in cohort
    
    cohort_summary(host, exclude)
        Count datasets per-cohort, excluding the given dataset types
        
        xena.cohort_summary(xena.PUBLIC_HUBS["pancanAtlasHub"], ["probeMap"])
    
    dataset_fetch(host, dataset, samples, probes)
        Probe values for give samples
    
    dataset_field(host, dataset)
        All field (probe) names in dataset
    
    dataset_field_examples(host, dataset, count)
        Field names in dataset, up to <count>
    
    dataset_field_n(host, dataset)
        Number of fields in dataset
    
    dataset_gene_probe_avg(host, dataset, samples, genes)
        Probe average, per-gene, for given samples
    
    dataset_gene_probes_values(host, dataset, samples, genes)
        Probe values in gene, and probe genomic positions, for given samples
    
    dataset_list(host, cohorts)
        Dataset metadata for datasets in the given cohorts
    
    dataset_metadata(host, dataset)
        Dataset metadata
    
    dataset_probe_signature(host, dataset, samples, probes, weights)
        Computed probe signature for given samples and weight array
    
    dataset_probe_values(host, dataset, samples, probes)
        Probe values for given samples, and probe genomic positions
        
        host = xena.PUBLIC_HUBS["pancanAtlasHub"]
        dataset = "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"
        samples = xena.dataset_samples(host, dataset, None)
        [position, [foxm1, tp53]] = xena.dataset_probe_values(host, dataset, samples, ["FOXM1", "TP53"])
    
    dataset_samples(host, dataset, limit)
        All samples in dataset (optional limit)
        
        samples = xena.dataset_samples(xena.PUBLIC_HUBS["pancanAtlasHub"], "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena", None)
    
    dataset_samples_n_dense_matrix(host, dataset)
        All samples in dataset (faster, for dense matrix dataset only)

    datasets_null_rows(host)
    
    feature_list(host, dataset)
        Dataset field names and long titles (phenotypic datasets)
    
    field_codes(host, dataset, fields)
        Codes for categorical fields
    
    field_metadata(host, dataset, fields)
        Metadata for given fields (phenotypic datasets)
    
    gene_transcripts(host, dataset, gene)
        Gene transcripts
    
    match_fields(host, dataset, names)
        Find fields matching names (must be lower-case)
    
    probe_count(host, dataset)
    
    probemap_list(host)
        Find probemaps
    
    ref_gene_exons(host, dataset, genes)
        Gene model
    
    ref_gene_position(host, dataset, gene)
        Gene position from gene model
    
    ref_gene_range(host, dataset, chr, start, end)
        Gene models overlapping range
    
    segment_data_examples(host, dataset, count)
        Initial segmented data rows, with limit
    
    segmented_data_range(host, dataset, samples, chr, start, end)
        Segmented (copy number) data overlapping range
    
    sparse_data(host, dataset, samples, genes)
        Sparse (mutation) data rows for genes
    
    sparse_data_examples(host, dataset, count)
        Initial sparse data rows, with limit
    
    sparse_data_match_field(host, field, dataset, names)
        Genes in sparse (mutation) dataset matching given names
    
    sparse_data_match_field_slow(host, field, dataset, names)
        Genes in sparse (mutation) dataset matching given names, case-insensitive (names must be lower-case)
    
    sparse_data_match_partial_field(host, field, dataset, names, limit)
        Partial match genes in sparse (mutation) dataset
    
    sparse_data_range(host, dataset, samples, chr, start, end)
        Sparse (mutation) data rows overlapping the given range, for the given samples
    
    transcript_expression(host, transcripts, studyA, subtypeA, studyB, subtypeB, dataset)

    
DATA

    LOCAL_HUB = 'https://local.xena.ucsc.edu:7223'
    PUBLIC_HUBS = {'gdcHub': 'https://gdc.xenahubs.net', 'icgcHub': 'https...
    excludeType = ['probeMap', 'probemap', 'genePredExt']
    
#### Contact
     http://xena.ucsc.edu/
     https://groups.google.com/forum/#!forum/ucsc-cancer-genomics-browser
     genome-cancer@soe.ucsc.edu


