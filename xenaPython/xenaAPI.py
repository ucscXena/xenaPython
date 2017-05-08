import xenaQuery as xena

def Gene_values (hub, dataset, samples, gene):
    values = xena.dataset_gene_values (hub, dataset, samples, [gene])
    return values[0]["scores"][0]

def Genes_values (hub, dataset, samples, genes):
    values = [x["scores"][0] for x in xena.dataset_gene_values (hub, dataset, samples, genes)]
    return values

def Probe_values (hub, dataset, samples, probe):
    values = xena.dataset_probe_values (hub, dataset, samples, [probe])
    return values[0]

def Probes_values (hub, dataset, samples, probes):
    values = xena.dataset_probe_values (hub, dataset, samples, probes)
    return values

def dataset_samples (hub,dataset):
    return xena.dataset_samples(hub, dataset)

def dataset_fields (hub, dataset):
    return xena.dataset_field (hub, dataset)
