"""
Utilities for xena queries.

A basic query example.
Queries are scheme expressions.

>>> import xenaPython.xenaQuery as xena
>>> xena.post("https://ucscpublic.xenahubs.net", "(+ 1 2)")
'3.0'

>>> xena.post("https://ucscpublic.xenahubs.net", "(let [x 2 y (+ x 3)] (* x y))")
'10.0'

Find all the identifiers of a dataset
>>> probes = xena.dataset_field("https://tcga.xenahubs.net", "TCGA.PANCAN.sampleMap/HiSeqV2")
>>> print probes[-10:]
'[u'ZXDB', u'ZXDC', u'ZYG11A', u'ZYG11B', u'ZYX', u'ZZEF1', u'ZZZ3', u'psiTPTE22', u'sampleID', u'tAKR']'

Find all the samples of a dataset
>>> samples = xena.dataset_samples("https://tcga.xenahubs.net", "TCGA.PANCAN.sampleMap/HiSeqV2")
>>> print samples[:5]
'[u'TCGA-S9-A7J2-01', u'TCGA-G3-A3CH-11', u'TCGA-EK-A2RE-01', u'TCGA-44-6778-01', u'TCGA-VM-A8C8-01']'

Find value matrix of a particular set of samples and identifiers (probes)
>>> values = xena.dataset_probe_values("https://tcga.xenahubs.net", "TCGA.PANCAN.sampleMap/HiSeqV2",["TCGA-44-6778-01","TCGA-44-6778-01"],["TP53"])
>>> print values
[[10.4169, 10.4169]]

Find value matrix of a particular set of samples and genes (identifier to gene mapping exists)
>>> values = xena.dataset_gene_values("https://tcga.xenahubs.net", "TCGA.PANCAN.sampleMap/HiSeqV2",["TCGA-44-6778-01","TCGA-44-6778-01"],["TP53"])
>>> print values

Looking up sample ids for the TCGA LGG cohort.

>>> r = xena.post("https://genome-cancer.ucsc.edu/proj/public/xena",
                  xena.patient_to_sample_query("TCGA.LGG.sampleMap",
                                               ["TCGA-CS-4938",
                                                "TCGA-HT-7693",
                                                "TCGA-CS-6665",
                                                "TCGA-S9-A7J2",
                                                "TCGA-FG-A6J3"]))
'{"TCGA.LGG.sampleMap":["TCGA-CS-4938-01","TCGA-CS-6665-01","TCGA-FG-A6J3-01","TCGA-HT-7693-01","TCGA-S9-A7J2-01"]}'

>>> r = xena.post("https://genome-cancer.ucsc.edu/proj/public/xena",
                  xena.find_sample_by_field_query("TCGA.LGG.sampleMap",
                                                    "_PATIENT",
                                                    ["TCGA-CS-4938",
                                                     "TCGA-HT-7693",
                                                     "TCGA-CS-6665",
                                                     "TCGA-S9-A7J2",
                                                     "TCGA-FG-A6J3"]))
'{"TCGA.LGG.sampleMap":["TCGA-CS-4938-01","TCGA-CS-6665-01","TCGA-FG-A6J3-01","TCGA-HT-7693-01","TCGA-S9-A7J2-01"]}'
>>> import json
>>> json.loads(r)
{u'TCGA.LGG.sampleMap': [u'TCGA-CS-4938-01', u'TCGA-CS-6665-01', u'TCGA-FG-A6J3-01', u'TCGA-HT-7693-01', u'TCGA-S9-A7J2-01']}
"""


import re
from functools import reduce

try:
    basestring
except NameError:
    basestring = (str, bytes)

def compose1(f, g):
    def composed(*args, **kwargs):
        return f(g(*args, **kwargs))
    return composed

# funcitonal composition, e.g.
# compose(f, g)(a, ...) == f(g(a, ...))
compose = lambda *funcs: reduce(compose1, funcs)

def quote(s):
    return '"' + s + '"'

def array_fmt(l):
    return '[' + ', '.join((quote(s) for s in l)) + ']'

def strip_first_url_dir(path):
    return re.sub(r'^[^/]*', '', path)

# proj/<proj>/xena/<proj>/<path>
# download/<proj>/xena/<path>
def name_to_url(base_url, name):
    return base_url.replace('/proj/', '/download/') + strip_first_url_dir(name)

# The strategy here is
#   o Do table scan on code to find codes matching field values
#   o Do IN query on unpack(field, x) to find rows matching codes
#   o Project to unpack(sample, x) to get sampleID code
#   o Join with code to get sampleID values
#
# Note the :limit on the table scan. This makes the table scan exit after we've
# found enough values, rather than continuing to the end. We can do this because
# enumerated values are unique. An alternative would be to index all the enumerated
# values in the db.


import json
try:
    from urllib.parse import urlparse, urlencode
    from urllib.request import urlopen, Request
    from urllib.error import HTTPError
except ImportError:
    from urlparse import urlparse
    from urllib import urlencode
    from urllib2 import urlopen, Request, HTTPError
    
headers = { 'Content-Type' : "text/plain" }

def post(url, query):
    """POST a xena data query to the given url."""
    req = Request(url + '/data/', query.encode(), headers)
    response = urlopen(req)
    result = response.read().decode('utf-8')
    return result

def quote(s):
    "quote a string value"
    return 'nil' if s is None else '"' + s + '"' # should escape "


def arrayfmt(l):
    "format an array"
    return '[' + ' '.join([marshall_param(x) for x in l]) + ']'

def marshall_param(p):
    "format a parameter"
    if isinstance(p, basestring):
        return quote(p)

    if isinstance(p, list):
        return arrayfmt(p)

    if p is None:
        return 'nil'
    return str(p)

def call(query_fn, params):
    "marshall parameters and build the lisp call form"
    return '(%s %s)' % (
        query_fn, ' '.join(map(marshall_param, params)))
