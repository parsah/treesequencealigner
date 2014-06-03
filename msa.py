from pairwise import NeedlemanWunsch, PositionWeightedMatcher
import math
from sequence import NeuriteSequence
import sequence
from collections import Counter
# from random import shuffle
from tree import TreeLogicFactory
from buffer import XMLBuildWriter

class MultipleSequenceDriver():
    ''' 
    A high-level class to perform multiple sequence pairwise.
    If iterate is an integer the align phase will run multiple times with a pwm based on the
    prior run.
    '''
    def __init__(self, queries, input_state):
        self.queries = queries
        self.costs = input_state.get_penalties() # set costs to core
        self.submat = input_state.get_submatrix() # set submatrix to core
        self.iterate = input_state.get_args()['iterate']
        self.composite = None # initially, no composite exists
        self.composite_score = 0
        self.alignment_file = input_state.alignment_file
        self.alns = []
        self.pwm = None
        self.node_types = input_state.node_types
        self.consensus_check_percent = .4

    def build_composite(self):
        ''' 
        Takes the first 2 sequences, aligns them and renders the pairwise
        alignment a composite (formerly "pre-consensus"). Every sequence thereof is then 
        subsequently pairwise-aligned to this composite.
        '''
        
        queries = self.queries
        # get the first two input sequences
        s0 = queries[0]
        s1 = queries[1]
        # pass them both into the tree--based Needleman--Wunsch algorithm.
        nw = NeedlemanWunsch(s1=s0, s2=s1, costs=self.costs, submat=self.submat, 
                                       node_types=self.node_types)
        first_align, second_align = nw.prettify()[1]
        # feed respective alignments into an analysis class and get consensus.
        composite = TreeLogicFactory(str1=first_align, 
                                       str2=second_align).get_alignment()
        # since the first two sequences have been aligned, focus on all others.
        for i in range(2, len(queries)):
            curr_seq = queries[i]
            print(curr_seq)
            nw = NeedlemanWunsch(s1=composite, s2=curr_seq, 
                                         costs=self.costs, submat=self.submat, 
                                         node_types=self.node_types)
            align_sA, align_sB = nw.prettify()[1]
            composite = TreeLogicFactory(str1=align_sA, 
                                       str2=align_sB).get_alignment()
            self.composite = composite # of type NeuriteSequence

    def align(self):
        ''' 
        A composite alignment is a manifestation of pairwise alignments
        between sequence pairs. In other words, sequence (i) and (i+1) are
        aligned and the resultant alignment (composite) is saved. Next,
        sequence (i+2) is aligned against this composite and the resultant
        alignment updates the composite. This logic is repeated for all
        queries. This function re-aligns the input queries back onto the
        composite so that the actual consensus sequence can be derived.
        If the 'iterate' parameter is set, the last part will repeat using
        a pwm produced using the previous run.
        '''
        if not self.composite:
            raise IOError('Composite alignment required.')
        if self.alignment_file:
            align_handle = open(self.alignment_file, 'w')
            
        name_lengths = []
        iterate = self.iterate
        if iterate < 1:
            iterate_count = 1 # Iterate until %change is within the iterate threshold given in iterate

        # Initialize pwm with equal weights per position
        pwm = PositionWeightedMatrix(self.composite,list([{'total':1,char:1} for char in self.composite.seq]))

        continue_iter = True
        iter_count = 0
        change = 1
        consensus_check_score = 0
        prev_consensus_score = 0
        # Iternate until change in consensus strength is sufficiently small or some number of iterations have been run
        while iterate < 1 and change > iterate or iterate >= 1 and iter_count < iterate:
            iter_count += 1

#        for curr_it in range(iterate_count):
            self.alns = [] # New alignments at each iteration, final iteration gives final alignment
            # Align each query with composite
            for curr_seq in self.queries:
                # Setting 'consensus=2' tells NW that s2 is the consensus and will prevent gaps from appearing in s1 alignemtn
                #nw = NeedlemanWunsch(s1=curr_seq, s2=self.composite, 
                pw_matcher = pairwise.PositionWeightedMatcher(sequence=curr_seq, pwm=pwm, 
                                    costs=self.costs, node_types=self.node_types)
                # sequence sA is the query alignment while sB is the composite.
                # we only need sA because sB does not change; all sequences are
                # mapped to this and therefore the resultant alignment is desired.

                align_sA, _ = pw_matcher.prettify()[1]
                self.alns.append(align_sA)

                if iterate is None or iter_count == iterate-1:
                    name_lengths.append(len(curr_seq.name))

            # Get the PWM given the alignments for the next iteration
            pwm = self.get_pwm()

            # Calculate change between previous and current alignment (even if not required)
            prev_consensus_score = consensus_check_score
            consensus_check_obj = self.build_consensus(threshold_ratio=self.consensus_check_percent)
            consensus_check_score = consensus_check_obj.score
            change = consensus_check_score - prev_consensus_score
#            if iterate < 1:

            #print('\n'+consensus_check_obj.raw_consensus)
            #print(str(consensus_check_score))
            #for aln in self.alns:
            #    print(aln)
        
        # Write alignments to file
        if self.alignment_file:
            total_space = max(12,max(name_lengths))+1
            align_handle.write(('Composite'+' '*(total_space-12))+self.composite.seq+'\n') # write header
            for index,curr_seq in enumerate(self.queries):
                align_handle.write(curr_seq.name+(' '*(total_space-len(curr_seq.name)))+(''.join(self.alns[index]))+'\n') # write header
            align_handle.close()

    # Generate a position weighted matrix of the form: array of {character:weight}
    def get_pwm(self):        
        # Determine composite score now that alignments exist, and generate pwm
        col_counts = []
        composite_count = 0
        height = len(self.queries)

        for col_num in range(len(self.composite.seq)): # process each column
            char_counts = self.enumerate_column(col_num)
            col_counts.append(char_counts)
            if '-' in char_counts:
                count_base = height - char_counts['-']
            else:
                count_base = height
            composite_count += count_base
        self.composite_score = float(composite_count)/height/len(self.composite.seq)

        composite_node_types = NeuriteSequence('PWM Node Types',''.join(list([self.node_types[char] for char in self.composite.seq])))
        self.pwm = PositionWeightedMatrix(composite_node_types,col_counts)

        return self.pwm

    def enumerate_column(self, num):
        '''
        Parses a specific column and counts the number of times a specific
        value is found in that respective column.
        '''

        char_counts = dict(Counter(''.join(aln[num] for aln in self.alns)))

#        a_column = [] # stores values for a single column
#        for row_num in range(len(self.alns)):            
#            for col_num in range(len(self.alns[row_num])):
#                if col_num == num:
#                    a_column.append(self.alns[row_num][col_num])
#        char_counts = dict(Counter(a_column))
        return char_counts

    def build_consensus(self, threshold_ratio, threshold_type='percent'):
        ''' 
        Iterates through each column, extracts respective values, counts its
        abundance and determines the consensus character given various neural 
        logic clauses. 
        Currently this might break if an expanded character set (of node types beyond A,C,T) is used
        '''
        if self.composite is None:
            self.build_composite()

        if len(self.alns) == 0:
            self.align()

        height = len(self.alns) # number of entries comprising alignment
        width = len(self.composite.seq) # all alignments are the same length

        if threshold_type == 'sqrt':
            threshold_count = threshold_ratio*math.sqrt(height) # consensus threshold
        else: # Default to percent
            threshold_count = threshold_ratio*height # consensus threshold

        consensus = None # actual generated consensus sequence 
        consensus_counts = []
        consensus_score = 0

        # store consensus as list (initially) so characters can be easily
        # placed easily; akin to array indexing.
        consensus = ['-'] * width 
        consensus_counts = [0] * width
        count_sum = 0 # for consensus score
        consensus_length = 0

        for col_num in range(width): # process each column
            #char_counts = self.enumerate_column(col_num)
            char_counts = self.pwm.pwm[col_num].copy()
            del char_counts['total'] # delete 'total'
            if '-' in char_counts:
                count_base = height - char_counts['-']
            else:
                count_base = height
            consensus_counts[col_num] = count_base

            if count_base >= threshold_count:
                consensus_length += 1
                count_sum += count_base
            ############################################################
            #### Logic for when only 1 base is exclusive to a column ###
            ############################################################
            if len(char_counts) == 1:
                chars = list(char_counts.keys()) # get the bases without counts
                char = chars[0] # get the key for that character
                consensus[col_num] = char # assign consensus character
                continue
            ############################################################
            #### Logic for when only 1 base and dashes are in column ###
            ############################################################
            if len(char_counts) == 2 and '-' in char_counts:
                # since we do not know what the known base is (either A, T, C),
                # subtract 1 from the dash count and use that to guide whether
                # the known base will become part of the consensus.
                if count_base >= threshold_count:
                    del char_counts['-'] # delete dash; remaining is character
                    char = list(char_counts.keys())[0]
                    consensus[col_num] = char # assign consensus character
                    continue
            #####################################################
            #### Logic for when 2 bases are found in a column ###
            #####################################################
            # This case should never happen. We could put in an error condition here...
            #if 'C' in char_counts and 'T' in char_counts:
            #    # we get the counts for both C and T, and contrast their
            #    # respective scores; assigning to the consensus whichever is
            #    # not only the largest score but also exceed the threshold.
            #    count_C, count_T = char_counts['C'], char_counts['T']
            #    if count_C >= count_T and count_C >= self.threshold:
            #        consensus[col_num] = 'C' # select C over T
            #    elif count_T >= count_C and count_T >= self.threshold:
            #        consensus[col_num] = 'T' # select T over C
            # choosing the better of A or C
            if 'C' in char_counts and 'A' in char_counts:
                # we get the counts for both C and A, and contrast their
                # respective scores; assigning to the consensus whichever is
                # not only the largest score but also exceed the threshold.
                if char_counts['A'] >= threshold_count:
                    consensus[col_num] = 'A' # There are enough A's
                elif char_counts['A'] + char_counts['C'] >= threshold_count:
                    consensus[col_num] = 'C' # There aren't enough A's, but there are enough characters

        consensus = NeuriteSequence(name='consensus', seq=''.join(consensus))
        consensus_score = float(count_sum)/height/consensus_length
        consensus_counts = consensus_counts
        #consensus_object = Consensus(consensus,consensus_score,consensus_counts)
        consensus_object = Consensus(consensus,threshold_ratio,threshold_type,height,consensus_score)
        return consensus_object

class PositionWeightedMatrix():
    '''
    Contains a string of node types (A,C,T) along with an array of dictionaries giving the weight of a position and
    an associated character.
    '''
    def __init__(self, node_type_sequence, pwm):
        self.node_type_sequence = node_type_sequence
        self.pwm = pwm
        # Calculate total weight of each column
        for pos in range(len(self.pwm)):
            self.pwm[pos]['total'] = sum(self.pwm[pos][key] for key in self.pwm[pos] if not key=='total' and not key=='-')

class Consensus():
    '''
    Contains a gapped and ungapped version of the consensus given an MSA and threshold values
    '''
    def __init__(self,consensus,threshold_ratio,threshold_type,height,score):
        self.raw_consensus = consensus.seq
        self.consensus = consensus.seq.replace('-','')
        self.threshold_ratio = threshold_ratio
        self.threshold_type = threshold_type
        self.height = height
        self.score = score
    
    def get_threshold_count(self):
        if self.threshold_type == 'sqrt':
            threshold_count = self.threshold_ratio*math.sqrt(self.height) # consensus threshold
        else: # Default to percent
            threshold_count = self.threshold_ratio*self.height # consensus threshold
        return threshold_count

    def get_threshold_percent(self):
        if self.threshold_type == 'sqrt':
            threshold_percent = self.threshold_ratio*math.sqrt(self.height)/self.height
        else: # Default to percent
            threshold_percent = self.threshold_ratio
        return threshold_percent


class ConsensusFilterFactory():
    ''' 
    When multiple sequence alignment is complete, alignments per query against
    the composite are produced. Each alignment must then be parsed given
    neural logic rules, enabling derivation of a sound consensus sequence.
    '''
    
#    def __init__(self, alignments, composite, threshold_ratio, threshold_type):
#        self.alns = alignments # matrix representing alignments
#        self.threshold = threshold_ratio
    def __init__(self, msa_object, consensus_object):
        self.msa_obj = msa_object
        self.consensus_obj = consensus_object
            
    def write(self, fname):
        ''' 
        Saves MSA analysis as an XML file.
        @param fname: Output filename
        '''
        w = XMLBuildWriter(fname, self.msa_obj, self.consensus_obj)
        w.write() # write the XML data-structure to the disk
    

