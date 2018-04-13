"""Methods for querying data from UCSC Xena hubs

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
"""
import json
from glob import glob
import re
import os.path
import sys

from . import xenaAPI
from . import xenaQuery

PUBLIC_HUBS = {
    'publicHub': 'https://ucscpublic.xenahubs.net',
    'tcgaHub': 'https://tcga.xenahubs.net',
    'icgcHub': 'https://icgc.xenahubs.net',
    'toilHub': 'https://toil.xenahubs.net',
    'pcawgHub': 'https://pcawg.xenahubs.net',
    'singlecellHub': 'https://singlecell.xenahubs.net',
    'pancanAtlasHub': 'https://pancanatlas.xenahubs.net',
    'treehouseHub': 'https://treehouse.xenahubs.net',
    'gdcHub': "https://gdc.xenahubs.net"
}

LOCAL_HUB = 'https://local.xena.ucsc.edu:7223'

#
# Dynamically build methods from queries/*.xq
#

def _to_snake(name):
    "convert camel case to snake case"
    s = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s).lower()

_MODULE_NAME = globals()['__name__']
_MODULE = sys.modules[_MODULE_NAME]

_DOCS = {
    'cohort_samples': 'All samples in cohort',
    'all_datasets_n': 'Count the number datasets with non-null cohort',
    'all_field_metadata': 'Metadata for all dataset fields (phenotypic datasets)',
    'cohort_summary': 'Count datasets per-cohort, excluding the given dataset types\n\nxena.cohort_summary(xena.PUBLIC_HUBS["pancanAtlasHub"], ["probeMap"])',
    'dataset_fetch': 'Probe values for give samples',
    'dataset_field': 'All field (probe) names in dataset',
    'dataset_field_examples': 'Field names in dataset, up to <count>',
    'dataset_field_n': 'Number of fields in dataset',
    'dataset_gene_probe_avg': 'Probe average, per-gene, for given samples',
    'dataset_gene_probes_values': 'Probe values in gene, and probe genomic positions, for given samples',
    'dataset_list': 'Dataset metadata for datasets in the given cohorts',
    'dataset_metadata': 'Dataset metadata',
    'dataset_probe_signature': 'Computed probe signature for given samples and weight array',
    'dataset_probe_values': """Probe values for given samples, and probe genomic positions

host = xena.PUBLIC_HUBS["pancanAtlasHub"]
dataset = "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"
samples = xena.dataset_samples(host, dataset, None)
[position, [foxm1, tp53]] = xena.dataset_probe_values(host, dataset, samples, ["FOXM1", "TP53"])""",
    'dataset_samples': 'All samples in dataset (optional limit)\n\nsamples = xena.dataset_samples(xena.PUBLIC_HUBS["pancanAtlasHub"], "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena", None)',
    'dataset_samples_n_dense_matrix': 'All samples in dataset (faster, for dense matrix dataset only)',
    'feature_list': 'Dataset field names and long titles (phenotypic datasets)',
    'field_codes': 'Codes for categorical fields',
    'field_metadata': 'Metadata for given fields (phenotypic datasets)',
    'gene_transcripts': 'Gene transcripts',
    'match_fields': 'Find fields matching names (must be lower-case)',
    'ref_gene_exons': 'Gene model',
    'ref_gene_position': 'Gene position from gene model',
    'ref_gene_range': 'Gene models overlapping range',
    'segment_data_examples': 'Initial segmented data rows, with limit',
    'segmented_data_range': 'Segmented (copy number) data overlapping range',
    'sparse_data': 'Sparse (mutation) data rows for genes',
    'sparse_data_examples': 'Initial sparse data rows, with limit',
    'sparse_data_match_field': 'Genes in sparse (mutation) dataset matching given names',
    'sparse_data_match_field_slow': 'Genes in sparse (mutation) dataset matching given names, case-insensitive (names must be lower-case)',
    'sparse_data_match_partial_field': 'Partial match genes in sparse (mutation) dataset',
    'sparse_data_range': 'Sparse (mutation) data rows overlapping the given range, for the given samples'
}

QUERIES = {}
def _create_methods():
    "inject query methods into global namespace"
    xq = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'queries', '*.xq')
    for method in glob(xq):
        with open(method) as x:
            query = x.read()
        name = _to_snake(os.path.splitext(os.path.basename(method))[0])
        # This really should be edn parse, instead of regex.
        params = re.split(
            r"\s+",
            re.sub(r"^[^[]+[[]([^]]*)[]].*$", r"\1", query, flags=re.DOTALL))

        QUERIES[name] = query
        call = eval('lambda %s: json.loads(xenaQuery.post(host, xenaQuery.call(QUERIES["%s"], [%s])))' % (
            ', '.join(['host'] + params), name, ', '.join(params)))

        if name in _DOCS:
            setattr(call, '__doc__', _DOCS[name])
        setattr(call, '__name__', name)
        setattr(call, '__module__', _MODULE_NAME)
        setattr(_MODULE, name, call)

_create_methods()
