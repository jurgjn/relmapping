l_stage = ['emb', 'l1', 'l2', 'l3', 'l4', 'ya']
l_rep = ['rep1', 'rep2', 'pep1', 'pep2']

def techreps_collapse(l_label_, include_raw=False):
    def collapse_(s): return s if s[-1].isdigit() else s[:-1]
    l_label = list(l_label_)
    l_collapsed = list(sorted(set(map(collapse_, l_label))))
    if include_raw:
        return list(sorted(l_collapsed + l_label))
    else:
        return l_collapsed

def techreps_retrieve(sample, config_):
    return list(sorted(v for k, v in config_.items() if k.startswith(sample)))
