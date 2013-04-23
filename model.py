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
    def __init__(self, name, sequence):
        self.seq = sequence
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
        return '{ ' + self.get_name() + '; ' + self.get_sequence() + ' }'
    
    def __repr__(self):
        return self.__str__()