'''
A high-level module which performs functionality having to do with writing 
results produced from analysis
'''
from xml.dom import minidom # used for XML writing
from xml.etree import ElementTree # used for XML reading
from sequence import ConsensusSequence, MultipleSequenceAlignment

def status_message(s, state='OK'):
    print(s + ' [' + state + ']')

class XMLBuildReader():
    ''' 
    Reads a user-provided XML build file and extracts its consensus sequence.
    '''
    def __init__(self, fname):
        self.fname = fname
        
    def parse(self):
        tree = ElementTree.parse(self.fname)
        composite = tree.find('composite/sequence').text
        alignments = []
        for elem in tree.iterfind('alignments/sequence'):
            alignments.append(elem.text)

        consensuses = []
        for elem in tree.iterfind('consensus'):
            cons_seq = elem.find('sequence').text
            threshold = float(elem.find('threshold').text)
            score = float(elem.find('score').text)
            consensus = ConsensusSequence(self.fname+'_T'+str(threshold),cons_seq,threshold,score)
            consensuses.append(consensus)
            
        return MultipleSequenceAlignment(composite,alignments,consensuses)
       # seq = tree.find('consensus/sequence').text # use XPath to get consensus
       # return ConsensusSequence(name = self.fname, seq = seq)

class XMLBuildWriter():
    ''' 
    Writes XML output files, otherwise known as a 'build'
    '''
    def __init__(self, fname, msa_obj, consensus_obj):
        self.handle = open(fname, 'w')
        self.doc = minidom.Document()
        self.root = None
        self.msa_obj = msa_obj
        self.consensus_obj = consensus_obj
        self._build_elements() # set the core collection of elements
        self._build_alignments() # set all alignments making up this analysis
        
    def _build_elements(self):
        self.root = self.doc.createElement('msa') # root element
        self.root.setAttribute('software', 'pasta')

        # set the COMPOSITE node
        composite_elem = self.doc.createElement('composite')
        composite_seq = self.doc.createElement('sequence') # composite sequence
        composite_len = self.doc.createElement('length') # composite length
        composite_score = self.doc.createElement('score') # composite score
    
        # set text values for respective elements
        text_composite_seq = self.doc.createTextNode(self.msa_obj.composite.get_sequence())
        text_composite_len = self.doc.createTextNode(str(len(self.msa_obj.composite.get_sequence())))
        text_composite_score = self.doc.createTextNode(str(round(self.msa_obj.composite_score,3)))
        
        # append element textual information to the respective element
        composite_seq.appendChild(text_composite_seq)
        composite_len.appendChild(text_composite_len)
        composite_score.appendChild(text_composite_score)
        composite_elem.appendChild(composite_seq)
        composite_elem.appendChild(composite_len)
        composite_elem.appendChild(composite_score)
        self.root.appendChild(composite_elem)
        
        # set the CONSENSUS node
        consensus_elem = self.doc.createElement('consensus')
        consensus_seq = self.doc.createElement('sequence') # consensus sequence
        consensus_seq_cond = self.doc.createElement('condensed-sequence') # composite sequence stripped of '-'
        consensus_len = self.doc.createElement('length') # consensus length
        consensus_thresh = self.doc.createElement('threshold') # consensus threshold
        consensus_score = self.doc.createElement('score') # consensus score
    
        # set text values for respective elements
        text_consensus_thresh = self.doc.createTextNode(str(self.consensus_obj.get_threshold_percent()))
        text_consensus_seq = self.doc.createTextNode(self.consensus_obj.raw_consensus)
        text_consensus_seq_cond = self.doc.createTextNode(self.consensus_obj.consensus)
        text_consensus_len = self.doc.createTextNode(str(len(self.consensus_obj.consensus)))
        text_consensus_score = self.doc.createTextNode(str(round(self.consensus_obj.score,3)))
        
        # append element textual information to the respective element
        consensus_seq.appendChild(text_consensus_seq)
        consensus_seq_cond.appendChild(text_consensus_seq_cond)
        consensus_len.appendChild(text_consensus_len)
        consensus_thresh.appendChild(text_consensus_thresh)
        consensus_score.appendChild(text_consensus_score)

        consensus_elem.appendChild(consensus_seq)
        consensus_elem.appendChild(consensus_seq_cond)
        consensus_elem.appendChild(consensus_len)
        consensus_elem.appendChild(consensus_thresh)
        consensus_elem.appendChild(consensus_score)
        self.root.appendChild(consensus_elem)
    
    def _build_alignments(self):
        # add each query sequence to the document
        aligns_elem = self.doc.createElement('alignments')
        aligns_elem.setAttribute('n', str(len(self.msa_obj.alns))) # number of alignments
        for _, seq in enumerate(self.msa_obj.alns):
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
