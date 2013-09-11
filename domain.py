import re
from matrix import ContingencyMatrix

class DomainSetBuilder():
    ''' 
    Given arguments relative to domain analysis, this class enables the ability
    to analyze a consensus sequence and extract out domains which satisfy
    user-provided arguments. 
    '''
    def __init__(self, consensus, win, max_gap, is_strip, is_enum=False):
        self.consensus = consensus # input consensus sequence object
        self.is_enumerate = is_enum # enumerate window size and max #/gaps
        self.win = win # sliding-window size
        self.max_gap = max_gap # maximum #/gaps in each domain
        self.is_strip = is_strip
        
    def build(self):
        wins = [self.win]
        if self.is_enumerate: # enumeration yields dynamic window and gaps
            wins = list(range(1, self.win + 1)) # window sizes from 1 .. win
            
        # iterate over window and gaps (if enumeration) and identify domains
        domains = {} # key => domain, value => count (abundance)
        for w in wins:
            if w > self.consensus.get_length():
                raise IOError('-win must be less than consensus length')
            else: # for each window, pull-out the respective domain
                for idx in range(len(self.consensus.seq)):
                    sub_str = self.consensus.seq[idx: idx + w] # reference
                    if self.is_strip: # if gaps are found, remove them
                        sub_str = sub_str.replace('-', '')
                    num_gap = sub_str.count('-') # count number of gaps
                    if len(sub_str) == w: # domain must equal sliding window
                        # only-gapped sequences are ignored
                        if not re.match('^-+$', sub_str) and num_gap <= self.max_gap:
                            if sub_str not in domains:
                                domains[sub_str] = 0
                            domains[sub_str] += 1 # increment its abundance
        return domains # return dictionary of domains and their abundances
    
class DomainAbundanceBuilder():
    ''' 
    This class enables the ability to build contingency matrices so as to
    enable the ability to contrast abundance of an item across different 
    groups. This then leads to the fact that a domain, D, could be abundant in
    1 of 3 possible scenarios: 
    1) D is present in query, Q, and not baseline, B.
    2) D is present in baseline, B, and not query, Q.
    3) D is present in both baseline, B, and query, Q.
    Thus, for all scenarios, the contingency matrix must model group-specific
    abundances so that truly over-represented domains can be identified.
    '''
    def __init__(self, query, baseline):
        self.query = query # query domain abundances
        self.baseline = baseline # baseline domain abundances
    
    def _query_exclusive_domains(self):
        names = set()
        for q in self.query:
            if q not in self.baseline:
                names.add(q)
        return names # references domains only found in the query
    
    def _baseline_exclusive_domains(self):
        names = set()
        for b in self.baseline:
            if b not in self.query:
                names.add(b)
        return names # references domains only found in the baseline
    
    def _intersection_domains(self):
        # get intersection of domains present in both query and baseline sets
        return set(self.query.keys()).intersection(self.baseline.keys())
    
    def build(self):
        # Suppose we set the following: 
        # Query (G) and Baseline (! G)
        # Domain (i) and all-other domains (! i) 
        size_query = sum(self.query.values()) # same as n(G)
        size_baseline = sum(self.baseline.values()) # same as n(! G)
        matrices = []
        
        # building matrices in query and baseline
        for d in self._intersection_domains():
            i_and_G = self.query[d] # domain count in query (idx: 0, 0)
            i_and_not_G = self.baseline[d] # domain count in query (idx: 0, 1)
            not_i_and_G = size_query - i_and_G # not-domain in query (idx: 1, 0)
            not_i_and_not_G = size_baseline - i_and_not_G # not-domain in baseline (idx: 1, 1)
            cm = ContingencyMatrix(node = d, i_g = i_and_G, 
                                   i_not_g = i_and_not_G, not_i_g=not_i_and_G, 
                                   not_i_not_g = not_i_and_not_G)
            matrices.append(cm)
            
        # building matrices in baseline only
        for d in self._baseline_exclusive_domains():
            i_and_G = 0.01 # domain count in query (idx: 0, 0)
            i_and_not_G = self.baseline[d] # domain count in query (idx: 0, 1)
            not_i_and_G = size_query - i_and_G # not-domain in query (idx: 1, 0)
            not_i_and_not_G = size_baseline - i_and_not_G # not-domain in baseline (idx: 1, 1)
            cm = ContingencyMatrix(node = d, i_g = i_and_G, 
                                   i_not_g = i_and_not_G, not_i_g=not_i_and_G, 
                                   not_i_not_g = not_i_and_not_G)
            matrices.append(cm)
        
        # building matrices in query only
        for d in self._query_exclusive_domains():
            i_and_G = self.query[d] # domain count in query (idx: 0, 0)
            i_and_not_G = 0.01 # domain count in query (idx: 0, 1)
            not_i_and_G = size_query - i_and_G # not-domain in query (idx: 1, 0)
            not_i_and_not_G = size_baseline - i_and_not_G # not-domain in baseline (idx: 1, 1)
            cm = ContingencyMatrix(node = d, i_g = i_and_G, 
                                   i_not_g = i_and_not_G, not_i_g=not_i_and_G, 
                                   not_i_not_g = not_i_and_not_G)
            matrices.append(cm)
        
        return matrices
    
class DomainPrettyPrinter():
    ''' 
    Prints statistically-significant domains after computing a hypergeometric
    p-value representing the abundance for it.
    @param domains: List of identified domains
    @param pval: User-provided p-value cutoff
    '''
    def __init__(self, domains, pval, out):
        self.domains = domains
        self.pval = pval
        self.fname = out
    
    def display(self):
        ''' 
        Prints statistically-significant domains to the screen.
        '''
        handle = open(self.fname, 'w')
        handle.write('Domain\tp-value\n') # write header
        handle.flush()
        for i in self.domains:
            dom_pval = i.get_hypergeometric_prob() # compute domain p-value
            if dom_pval <= self.pval: # domain must be less than p-value cutoff
                handle.write(i.name + '\t' + str(round(dom_pval, 6)) + '\n')
                handle.flush()
        handle.close()
        