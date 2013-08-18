'''
Provides the ability to model and represent various concrete objects that
are to be used throughout application runtime.
'''

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
    def __init__(self, name, seq):
        super(ConsensusSequence, self).__init__(name, seq)
        self.num_gaps = 0
        