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
    'treehouseHub': 'https://xena.treehouse.gi.ucsc.edu',
    'gdcHub': "https://gdc.xenahubs.net",
    'reference': "https://reference.xenahubs.net"
}

LOCAL_HUB = 'https://local.xena.ucsc.edu:7223'

excludeType = ["probeMap", "probemap", "genePredExt"]

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
    'probemap_list': 'Find probemaps',
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

def _filehash(file):
    import hashlib
    sha_hash = hashlib.sha1()
    with open(file,"rb") as f:
        # Read and update hash string value in blocks of 4K
        for byte_block in iter(lambda: f.read(4096), b""):
            sha_hash.update(byte_block)
    return sha_hash.hexdigest()

def _name_to_path(name):
    import os
    home = os.environ['HOME']
    return os.path.join(home, 'xena', 'files', name)

def _file_loaded(host, name):
    status = dataset_status(host, name)
    try:
        return status[0]['status'] == 'loaded'
    except:
        return False

def _hashes_match(host, name):
    hashes = dataset_sources(host, name)
    return all(_filehash(_name_to_path(file['name'])) == file['hash'] for file in hashes)

# XXX need timeout and error handling
# should this do the file copy, too?
def load_file(name, host = "http://127.0.0.1:7222"):
    import requests
    import time
    requests.post(host + '/update/', data = {'file': _name_to_path(name)})
    while not _file_loaded(host, name) or not _hashes_match(host, name):
        time.sleep(20)

ipython_instance = None

def open_browser(url = 'heatmap/', columns = []):
    import json, urllib
    if len(columns) > 0:
        withhost = [{'name': c['name'], 'fields': c['fields'], 'host': 'notebook:'}
                for c in columns]
        cols = '?columns=' + urllib.parse.quote(json.dumps(withhost), safe='~()*!.\'')
    else:
        cols = ''
    ipython_instance.run_cell("%%%%javascript\nwindow.open(window.xenabrowser.url + '/%s%s')" % (url, cols))

# notebook support
def load_ipython_extension(ipython):
    "jupyter support method"
    global ipython_instance
    ipython_instance = ipython
    import os
    try:
        host = os.environ["XENA_BROWSER"]
    except KeyError:
        host = 'https://xenabrowser.net'
    dir = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(dir, "jupyter.py"), "r") as f:
        init_py_comm = f.read()
    with open(os.path.join(dir, "jupyter.js"), "r") as f:
        init_js_comm = f.read()
    ipython.run_cell('%%javascript\nwindow.xenabrowser = {url: "' + host + '"};')
    ipython.run_cell('%%javascript\n' + init_js_comm)
    ipython.run_cell('import ' + globals()['__name__'])
    ipython.run_cell(init_py_comm, False, True, False)
    os.system(os.path.join(dir, 'startxena.sh'))
    print("loading")

def unload_ipython_extension(ipython):
    "jupyter support method"
    print("unloading")
