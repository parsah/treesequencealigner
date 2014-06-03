'''
Provides the ability to model and represent various concrete objects that
are to be used throughout application runtime.
'''

default_nodetypes = {'A':'A','C':'C','T':'T'}

# Reads in a file containing node types (A,C,T) and a string of specific characters of that type
# Format for a given line should be: "<node-type>:<character string>", for example: "C:BRPD", or the simple case: "C:C"
def parse_node_types(fname):
    node_types = {}
    for line in open(fname):
        line = line.strip()
        if len(line) == 0:
            break
        else:
            vals = line.split(':')
            node_types[vals[0]] = vals[1]
    return node_types

def generate_identity_matrix(nodetypes=default_nodetypes):
    submat = {}
    ac_types = ['A','C']
    for char1 in nodetypes.keys():
        nodetype1 = nodetypes[char1]
        for char2 in nodetypes.keys():
            nodetype2 = nodetypes[char2]
            if nodetype1 == nodetype2 or nodetype1 in ac_types and nodetype2 in ac_types:
                submat[(char1,char2)] = 1
                submat[(char2,char1)] = 1
            else:
                submat[(char1,char2)] = -40
                submat[(char2,char1)] = -40
    return submat


class NeuriteSequence():
    ''' 
    A NeuriteSequence object is simply a FASTA object but solely references
    neuronal sequences.
    @param sequence: input sequence string.
    @param header: header for the respective sequence.
    '''
    def __init__(self, name, seq):
        self.seq = seq
        self.name = name
        
    def get_length(self):
        ''' 
        Returns the length of the neurite-sequence object.
        @return: integer referencing sequence length.
        '''
        return len(self.seq)
    
    def get_name(self):
        '''
        Returns the name of the neurite-sequence object.
        @return: string referencing sequence object.
        '''
        return self.name
    
    def get_sequence(self):
        ''' 
        Returns the current sequence as a string.
        @return: string referencing the current neuronal sequence.
        '''
        return self.seq
    
    def __str__(self):
        return '{ ' + self.get_name() + '; ' + self.get_sequence() +\
            '; ' + str(self.get_length()) + ' bases }'
    
    def __repr__(self):
        return self.__str__()
    
class ConsensusSequence(NeuriteSequence):
    ''' 
    Encapsulates a consensus sequence so it can be parsed and analyzed using 
    domain-analysis mode.
    '''
    def __init__(self, name, seq, threshold=None, score=None):
        super(ConsensusSequence, self).__init__(name, seq)
        self.num_gaps = 0
        self.threshold = threshold
        self.score = score
        self.ungapped_consensus = seq.replace('-','')
        
tree_sequence_types = ['complete_tree','incomplete_tree','multiple_trees_complete','multiple_trees_incomplete']

def tree_sequence_type(seq,node_types=default_nodetypes):
    a_nodes = 0
    #AStack.append(-1) # Begin the AStack with a sentinal value
    complete_trees = 0
    tree_complete = True
    # Visit each character in the sequence
    for index in range(0,len(seq)):
        tree_complete = False
        # If the current node is an A-node
        if seq[index] in node_types['A']:
            # push it onto the top of the stack
            a_nodes += 1
 
        # If the current node is a T-node
        elif seq[index] in node_types['T']:
            if a_nodes == 0:
                complete_trees += 1
                tree_complete = True
            else:
                a_nodes -= 1

    if complete_trees == 0:
        return 'incomplete_tree'
    elif tree_complete:
        if complete_trees == 1:
            return 'complete_tree'
        else: # more than one complete tree
            return 'multiple_trees_complete'
    else:
        return 'multiple_trees_incomplete'

class MultipleSequenceAlignment():
    '''
    An object class containing a complete multiple sequence alignment, including component
    sequences with MSA gaps, composite sequence, and optionally one or more consensus
    sequences.
    '''

    def __init__(self, composite, alignments, consensuses=[]):
        self.composite = composite
        self.alignments = alignments
        self.consensuses = consensuses

    def add_consensus(consensus):
        self.consensuses.append(consensus)
