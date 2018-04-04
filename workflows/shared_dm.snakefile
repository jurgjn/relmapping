l_stage = ['emb', 'l1', 'l2', 'l3', 'l4', 'ya']
l_rep = ['rep1', 'rep2', 'pep1', 'pep2']

def techreps_collapse(l_label):
    def collapse_(s): return s if s[-1].isdigit() else s[:-1]
    return list(sorted(set(map(collapse_, l_label))))

def techreps_retrieve(sample, config_):
    return list(sorted(v for k, v in config_.items() if k.startswith(sample)))
