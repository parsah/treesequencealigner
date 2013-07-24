import concurrent.futures
import sys
import os
import core
import math
from Bio.SubsMat import MatrixInfo
from Bio import SeqIO
from model import NeuriteSequence
from collections import Counter
from random import shuffle
from core import stdout

version = 0.1

class TreeIndexLogic():
    ''' 
    A class which captures tree-specific logic given two characters.
    @param char1: First character
    @param char2: Second character
    '''
    def __init__(self, char1, char2):
        self.char1 = char1
        self.char2 = char2
        
    def get(self):
        if (self.char1 == 'A' and self.char2 == 'C') or (self.char2 == 'A' and self.char1 == 'C'):
            return 'A'
        if (self.char1 == 'A' and self.char2 == '-') or (self.char2 == 'A' and self.char1 == '-'):
            return 'A'
        if (self.char1 == 'T' and self.char2 == '-') or (self.char2 == 'T' and self.char1 == '-'):
            return 'T'
        if (self.char1 == 'C' and self.char2 == '-') or (self.char2 == 'C' and self.char1 == '-'):
            return 'C'

class TreeLogicFactory():
    '''
    Parses and processes the consensus string to ultimately yield a single
    string which encapsulate the pairwise alignment.
    '''
    def __init__(self, str1, str2):
        self.str1 = str1
        self.str2 = str2
        
    def get_alignment(self):
        ''' 
        Simple function to merge two strings and produce a consensus.
        @return: NeuriteSequence object representing the consensus sequence.
        '''
        consensus = ''
        for idx, char1 in enumerate(self.str1):
            char2 = self.str2[idx]
            if char1 == self.str2[idx]:
                consensus += char1
            else:
                # Apply neuronal logic given two specific characters.
                consensus += TreeIndexLogic(char1, char2).get()
        return NeuriteSequence(name='alignment', sequence=consensus)

class MultipleSequenceDriver():
    ''' 
    A high-level class to perform multiple sequence alignment.
    '''
    def __init__(self, queries, input_state):
        self.queries = queries
        self.costs = input_state.get_penalties() # set costs to core
        self.submat = input_state.get_submatrix() # set submatrix to core
        self.preconsensus = None # initially, no pre-consensus exists
        self.alignment_file = input_state.alignment_file
        self.random_order = False

    def set_randomize(self, bool):
        self.random_order = bool

    def build_preconsensus(self):
        ''' 
        Takes the first 2 sequences, aligns them and renders the pairwise
        alignment a pre-consensus. Every sequence thereof is then subsequently
        pairwise-aligned to this pre-consensus.
        '''
        
        queries = self.queries
        if self.random_order:
            # Randomize the order of the sequences
            shuffle(queries)

        stdout('--- Multiple sequence alignment mode ---')
        # get the first two input sequences
        s0 = queries[0]
        s1 = queries[1]
        # pass them both into the tree--based Needleman--Wunsch algorithm.
        nw = core.NeedlemanWunsch(s1=s0, s2=s1, costs=self.costs, submat=self.submat, nodeTypes=core.default_nodetypes)
        first_align, second_align = nw.prettify()[1]
        # feed respective alignments into an analysis class and get consensus.
        consensus = TreeLogicFactory(str1=first_align, 
                                       str2=second_align).get_alignment()
        # since the first two sequences have been aligned, focus on all others.
        for i in range(2, len(queries)):
            curr_seq = queries[i]
            nw = core.NeedlemanWunsch(s1=consensus, s2=curr_seq, 
                                         costs=self.costs, submat=self.submat, 
                                         nodeTypes=core.default_nodetypes)
            align_sA, align_sB = nw.prettify()[1]
            consensus = TreeLogicFactory(str1=align_sA, 
                                       str2=align_sB).get_alignment()
            self.preconsensus = consensus # of type NeuriteSequence.

    def align(self):
        ''' 
        A pre-consensus alignment is a manifestation of pairwise alignments
        between sequence pairs. In other words, sequence (i) and (i+1) are
        aligned and the resultant alignment (pre-consensus) is saved. Next,
        sequence (i+2) is aligned against this pre-consensus and the resultant
        alignment updates the pre-consensus. This logic is repeated for all
        queries. This function re-aligns the input queries back onto the
        pre-consensus so that the actual consensus sequence can be derived.
        '''
        if not self.preconsensus:
            raise IOError('Pre-consensus alignment required.')
        if self.alignment_file:
            align_handle = open(self.alignment_file, 'w')
            
        alignments = [] # references alignments against the pre-consensus
        name_lengths = []
        for curr_seq in self.queries:
            # Setting 'consensus=2' tells NW that s2 is the consensus and will prevent gaps from appearing in s1 alignemtn
            nw = core.NeedlemanWunsch(s1=curr_seq, s2=self.preconsensus, 
                                    costs=self.costs, submat=self.submat, 
                                    nodeTypes=core.default_nodetypes,consensus=2)
            # sequence sA is the query alignment while sB is the pre-consensus.
            # we only need sA because sB does not change; all sequences are
            # mapped to this and therefore the resultant alignment is desired.
            align_sA, align_sB = nw.prettify()[1]
            alignments.append(list(align_sA))
            print(align_sA)
            name_lengths.append(len(curr_seq.name))
            
        # Write alignments to file
        if self.alignment_file:
            total_space = max(12,max(name_lengths))+1
            align_handle.write(('Preconsensus'+' '*(total_space-12))+self.preconsensus.seq+'\n') # write header
            for index,curr_seq in enumerate(self.queries):
                align_handle.write(curr_seq.name+(' '*(total_space-len(curr_seq.name)))+(''.join(alignments[index]))+'\n') # write header
            align_handle.close()
            
        return alignments # return the matrix of all alignments

class ConsensusFilterFactory():
    ''' 
    When multiple sequence alignment is complete, alignments per query against
    the pre-consensus are produced. Each alignment must then be parsed given
    neural logic rules, enabling derivation of a sound consensus sequence.
    '''
    
    def __init__(self, alignments, threshold_ratio, threshold_type):
        self.alignments = alignments # matrix representing alignments
        if threshold_type == 'sqrt':
            self.threshold_count = threshold_ratio*math.sqrt(len(alignments)) # consensus threshold
        else: # Default to percent
            self.threshold_count = threshold_ratio*len(alignments) # consensus threshold
        self.height = len(alignments) # number of entries comprising alignment
        self.width = len(alignments[0]) # all alignments are the same length
        
    def enumerate_column(self, num):
        '''
        Parses a specific column and counts the number of times a specific
        value is found in that respective column.
        '''
        a_column = [] # stores values for a single column
        for row_num in range(len(self.alignments)):
            for col_num in range(len(self.alignments[row_num])):
                if col_num == num:
                    a_column.append(self.alignments[row_num][col_num])
        char_counts = dict(Counter(a_column))
        # set the counts as a fraction so it can be easily computed in relation
        # to the threshold. EDIT: using count threshold instead to enable other types of thresholding
        #for char in char_counts:
        #    frac = char_counts[char]/float(self.height)
        #    char_counts[char] = round(frac, 2)
        return char_counts
    
    def build_consensus(self):
        ''' 
        Iterates through each column, extracts respective values, counts its
        abundance and determines the consensus character given various neural 
        logic clauses. 
        '''
        # store consensus as list (initially) so characters can be easily
        # placed easily; akin to array indexing.
        consensus = ['-'] * self.width 
        for col_num in range(self.width): # process each column
            char_counts = self.enumerate_column(col_num)
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
                count_base = self.height - char_counts['-']
                if count_base >= self.threshold_count:
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
                if char_counts['A'] >= self.threshold_count:
                    consensus[col_num] = 'A' # There are enough A's
                elif char_counts['A'] + char_counts['C'] >= self.threshold_count:
                    consensus[col_num] = 'C' # There aren't enough A's, but there are enough characters
        
        print()
        print('consensus:')
        print(''.join(consensus))


class PairwiseDriver():
    '''
    Executes the pairwise application
    '''
    def __init__(self, targets, queries, input_state):
        self.targets = targets # core operates given targets and queries
        self.queries = queries
        self.costs = input_state.get_penalties() # set costs to core
        self.submat = input_state.get_submatrix() # set submatrix to core
        self.num_complete = 0# len(self.priorCompletions) # for how many sequences have been aligned
        #stdout(str(self.num_complete)+" complete of "+str(len(targets)))

        # Set openMode to append if some targets have already been run and completed
        if self.num_complete > 0:
            openMode = 'a'
        else:
            openMode = 'w'
            
        self.scorehandle = open(input_state.get_args()['o'], openMode) # output file
        self.alignhandle = None
        if input_state.get_args()['a'] is not None:
            self.alignhandle = open(input_state.get_args()['a'], openMode) # alignments file
            
        self.num_workers = input_state.get_args()['n']
        # Get node type lists
        if input_state.get_args()['nodeTypes'] is None:
            self.nodeTypes = core.default_nodetypes
        else:
            self.nodeTypes = core.parse_nodetypes(input_state.get_args()['nodeTypes'])

    # Initialize the core given query sequences and input arguments
    def start(self):
        stdout('--- Pairwise alignment mode ---')
        executor = concurrent.futures.ProcessPoolExecutor(self.num_workers)
        try:
            for target in self.targets: # per fasta, create a concurrent job, f.
                f = executor.submit(_aligner, target, queries, self.costs, self.submat, self.nodeTypes)
                f.add_done_callback(self._callback)
            executor.shutdown()
            self.close_output_buffers()
            stdout('** Analysis Complete **')
        except KeyboardInterrupt:
            executor.shutdown()

    # Close all I/O buffers such as file handles
    def close_output_buffers(self):
        if self.alignhandle is not None:
            self.alignhandle.close()
        self.scorehandle.close()

    # Get the headers, i.e. top-most row for the score matrix
    def _create_header(self, results):
        h = '\t' +'\t'.join([h[-1] for h in results])
        self.scorehandle.write(h + '\n')
        self.scorehandle.flush()
        
    # Callback function once a thread is complete
    def _callback(self, return_val):
        res = return_val.result() # get result once thread is complete
        target, results = res
        results = sorted(results, key=lambda x: x[-1]) # sort by query (last item)
        if self.num_complete == 0: # for the first result, write headers
            self._create_header(results)

        # save scores to the alignment matrix
        scores = '\t'.join([str(s[0]) for s in results])
        self.scorehandle.write(target + '\t' + scores + '\n')
        self.scorehandle.flush()

        # also save actual alignment string
        if self.alignhandle is not None:
            for r in results:
                if r[1] is not None:
                    align_target, align_query = r[1]
                    out_str = target + '\t' + r[-1] +'\t'+ align_target +'\t'+ align_query
                    self.alignhandle.write(out_str + '\n')
                    self.alignhandle.flush()
        self.num_complete += 1
        stdout(' --> ' + target + ' [OK] '+str(self.num_complete) +
            ' of '+ str(len(self.targets))) # print-stdout progress

# Maps each query sequence against a set of targets (itself)
def _aligner(target, queries, costs, submat, nodeTypes):
    results = [] # K => target, V => aligned queries 
    # get the gap and substitution matrix
    for query in queries:
        NW = core.NeedlemanWunsch(target, query, costs, submat, nodeTypes)
        output = NW.prettify()
        results.append(output)
    return target.name, results

if __name__ == '__main__':
    try:
        from parameter import ArgumentValidator, CommandLineParser, InputWrapperState
        args = CommandLineParser().parse_args()
        ArgumentValidator(args) # test all arguments are correct

        input_state = InputWrapperState(args)
        input_state.assign_matrix() # parse in-built or custom matrix
        targets = input_state.parse_fasta(input_state.fname) # next, parse fasta file
        if input_state.fname2 is None:
            queries = targets
        else:
            queries = input_state.parse_fasta(input_state.fname2)
        if args['mode'] == 'local':
            driver = PairwiseDriver(targets, queries, input_state)
            driver.start() # start only pairwise alignment
            # TODO there is no 'consensus' mode; it is part of msa mode. 
#         elif args['mode'] == 'consensus' and args['a'] is not None:
#             alignments = parse_alignments(args['a'])
#             consensus_fact = ConsensusFilterFactory(alignments, args['t'], args['threshold_type'])
#             consensus_fact.build_consensus()
        else: # else, start multiple-sequence alignment (MSA)
            driver = MultipleSequenceDriver(queries, input_state)
            driver.build_preconsensus()
            # map queries back onto consensus and build a filtered consensus
            alignments = driver.align()
            consensus_fact = ConsensusFilterFactory(alignments, args['t'], args['threshold_type'])
            consensus_fact.build_consensus()
    except (IOError, KeyboardInterrupt, IndexError) as e:
        stdout(str(e)+'\n')
