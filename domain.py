import re

class DomainSetBuilder():
    ''' 
    Given arguments relative to domain analysis, this class enables the ability
    to analyze a consensus sequence and extract out domains which satisfy
    user-provided arguments. 
    '''
    def __init__(self, consensus, win, max_gap, is_enum=False):
        self.consensus = consensus # input consensus sequence object
        self.is_enumerate = is_enum # enumerate window size and max #/gaps
        self.win = win # sliding-window size
        self.max_gap = max_gap # maximum #/gaps in each domain
        
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
                    num_gap = sub_str.count('-') # count number of gaps
                    if len(sub_str) == w: # domain must equal sliding window
                        # only-gapped sequences are ignored
                        if not re.match('^-+$', sub_str) and num_gap <= self.max_gap:
                            if sub_str not in domains:
                                domains[sub_str] = 0
                            domains[sub_str] += 1 # increment its abundance
        return domains # return dictionary of domains and their abundances