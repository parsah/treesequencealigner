'''
A high-level module which performs functionality having to do with writing 
results produced from analysis
'''
from xml.dom import minidom # used for XML writing

class XMLBuildWriter():
    ''' 
    Writes XML output files, otherwise known as a 'build'
    '''
    def __init__(self, fname, msa_obj):
        self.handle = open(fname, 'w')
        self.doc = minidom.Document()
        self.root = None
        self.msa_obj = msa_obj
        self._build_elements() # set the core collection of elements
        self._build_alignments() # set all alignments making up this analysis
        
    def _build_elements(self):
        self.root = self.doc.createElement('msa') # root element
        self.root.setAttribute('software', 'pasta')
        
        # set the consensus node
        consensus_elem = self.doc.createElement('consensus')
        consensus_seq = self.doc.createElement('sequence') # consensus sequence
        consensus_len = self.doc.createElement('length') # consensus length
        consensus_thresh = self.doc.createElement('threshold') # consensus threshold
    
        # set text values for respective elements
        text_consensus_thresh = self.doc.createTextNode(str(self.msa_obj.threshold))
        text_consensus_seq = self.doc.createTextNode(self.msa_obj.consensus.get_sequence())
        text_consensus_len = self.doc.createTextNode(str(len(self.msa_obj.consensus.get_sequence())))
        
        # append element textual information to the respective element
        consensus_seq.appendChild(text_consensus_seq)
        consensus_len.appendChild(text_consensus_len)
        consensus_thresh.appendChild(text_consensus_thresh)
        consensus_elem.appendChild(consensus_seq)
        consensus_elem.appendChild(consensus_len)
        consensus_elem.appendChild(consensus_thresh)
        self.root.appendChild(consensus_elem)
    
    def _build_alignments(self):
        # add each query sequence to the document
        aligns_elem = self.doc.createElement('alignments')
        aligns_elem.setAttribute('n', str(len(self.msa_obj.alignments))) # number of alignments
        for i, seq in enumerate(self.msa_obj.alignments):
            text_align_seq = self.doc.createTextNode(seq)
            an_align_elem = self.doc.createElement('sequence')
            an_align_elem.appendChild(text_align_seq)
            aligns_elem.appendChild(an_align_elem)
        
        # write XML tree hierarchy to the user output file
        self.root.appendChild(aligns_elem)
        self.doc.appendChild(self.root) # add root element to the document 
    
    def write(self):
        self.doc.writexml(self.handle, addindent="  ", newl='\n')
        self.handle.close()
