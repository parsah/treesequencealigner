import argparse
import platform
import concurrent.futures
import sys
import os
import factory
from Bio.SubsMat import MatrixInfo
from Bio import SeqIO
from model import NeuriteSequence

# Validates user-provided command-line arguments
class ArgumentValidator():
    def __init__(self, args):
        self.args = args
        self.check_python()
        self.check_biopython()
        self.check_args()

    # Check python 3.3 is available
    def check_python(self):
        if platform.python_version_tuple()[0:2] >= ('3', '0'): # >= 3.0
            out('Python v. '+ str(platform.python_version()) + ' found [OK]')
        else:
            raise RuntimeError('Python 3.2+ recommended')

    # BioPython must be installed
    def check_biopython(self):
        try:
            import Bio
            out('BioPython v.' + Bio.__version__ + ' found [OK]')
        except ImportError:
            raise ImportError('Make sure BioPython is installed')

    # Checks user-provided arguments are valid
    def check_args(self):
        return all([self.test_num_workers(), self.test_mutual_matrices(),
                self.test_valid_matrix()])

    # Test either a custom matrix or in-built matrix is selected
    def test_mutual_matrices(self):
        if self.args['custom'] and self.args['matrix']:
            raise IOError('Either a custom or in-built matrix must be selected')
        else:
            return True

    # Test a valid substitution matrix is selected
    def test_valid_matrix(self):
        all_matrices = set(MatrixInfo.available_matrices) # al sub. matrices
        if self.args['matrix'] in all_matrices or self.args['custom']:
            return True
        else:
            err = 'An in-built (URL below) or custom matrix is required.\n'+\
            'http://biopython.org/DIST/docs/api/Bio.SubsMat.MatrixInfo-module.html'
            raise IOError(err)

    # Test a valid number of workers are provided
    def test_num_workers(self):
        if self.args['n'] >= 1:
            return True
        else:
            raise IOError('>= 1 worker processes must be provided')

# Helper-class to parse input arguments
class CommandLineParser():
    def __init__(self):
        desc = 'Script to execute exhaustive brute-force pairwise alignment'
        u='%(prog)s [options]' # command-line usage
        self.parser = argparse.ArgumentParser(description=desc, add_help=False, usage=u)
        self._init_params()

    # Create parameters to be used throughout the application
    def _init_params(self):
        param_reqd = self.parser.add_argument_group('Required Parameters')
        param_opts = self.parser.add_argument_group('Optional Parameters')
        #param_costs = self.parser.add_argument_group('Specific Costs')

        # Specify required arguments
        param_reqd.add_argument('-f', metavar='FILE', required=True,
                    help='Input fasta file [na]')
        param_reqd.add_argument('--mode', metavar='MODE', required=True,
                    help='Either pairwise (local) or multiple-sequence (msa) [req]')
        
        # Specify optional arguments
        param_opts.add_argument('-f2', metavar='FILE', required=False, default=None,
                    help='Second input fasta file [na]')
        param_opts.add_argument('--gap', metavar='INT', default=-8, type=int,
                    help='Gap extension penalty [-8]')
        param_opts.add_argument('--gapopen', metavar='INT', default=0, type=int,
                    help='Gap open penalty (in addition to, not instead of, extension penalty) [0]')
        param_opts.add_argument('-custom', metavar='FILE', default=None,
                    help='Custom substitution matrix [na]')
        param_opts.add_argument('-nodeTypes', metavar='FILE', default=None,
                    help='Node Type Specifications [na]')
        param_opts.add_argument('-matrix', metavar='STR', default=None,
                    help='Matrix name; see Biopython MatrixInfo for all matrices [na]')
        param_opts.add_argument('-n', metavar='INT', default=2, type=int,
                    help='Number of worker processes [2]')
        param_opts.add_argument('-o', metavar='FILE', default='scores.tab', 
                    help='File to write/append output [scores.tab]')
        param_opts.add_argument('-a', metavar='FILE', default='', 
                    help='File to write/append alignments [none]')
        param_opts.add_argument('-s', metavar='STR', default='alignment', 
                    help='Type of score to write to output file [alignment]\n\talignment,gaps,excess_gaps,short_normalized,long_normalized')
        param_opts.add_argument('--forceQuery', action='store_const', const=True, default=False)
        param_opts.add_argument('-h','--help', action='help',
                    help='Show this help screen and exit')

    # Get the arguments for each parameter
    def parse_args(self):
        return vars(self.parser.parse_args()) # parse arguments

# An InputStateWrapper wraps data used as input, eg. fasta file,
# custom matrix, user-provided arguments, etc.
class InputWrapperState():
    def __init__(self, args):
        self.args = args # reference user-provided arguments
        self.subsmat = None # references data for substitution matrix
        self.fname = args['f'] # input filename
        self.fname2 = args['f2'] # input filename
        self.score_type = args['s']

    # Get arguments
    def get_args(self):
        return self.args

    # Get arguments relative to penalties
    def get_penalties(self):
        cost_ids = ('gap','gapopen') # all possible costs, might in the future include a separate gap open and gap extension cost
        return {k: self.args[k] for k in cost_ids} # get costs per penality

    # Get the substitution matrix which will be used
    def get_submatrix(self):
        return self.subsmat

    # Get the type of score to be used in the output file
    def get_scoretype(self):
        return self.score_type
        
    # Trivial function to parse a fasta file
    def parse_fasta(self, fname):
        parsed_fasta = list(SeqIO.parse(fname, 'fasta')) # easy indexing
        queries = [] # references list of parsed sequences
        for i in parsed_fasta: # cast as a neuronal sequence; easy modeling.
            s = NeuriteSequence(sequence=str(i.seq), name=i.name)
            queries.append(s)
        out(str(len(queries)) + ' queries parsed [OK]')
        return queries # return set of fasta entries

    # Trivial function to write parameter arguments to a file 
    def write_args(self):
        outhandle = open('param_args.tab', 'w')
        outhandle.write('Parameter\tValue\n') # write header
        outhandle.flush()
        for i in sorted(self.args):
            outhandle.write(i+'\t'+str(self.args[i])+'\n') # write args
            outhandle.flush()
        outhandle.close()
        out('') # write new line

    # Set the desired matrix the user wishes to add
    def assign_matrix(self):
        if not self.args['matrix']: # if no in-built matrix, parse custom matrix
            self.subsmat = self.__parse_custom_matrix() # custom matrix
        else:
            self.subsmat = self.__parse_inbuilt_matrix()

    # Get the user-provided in-built matrix
    def __parse_inbuilt_matrix(self):
        return getattr(MatrixInfo, args['matrix']) # get substitution matrix

    # Function to parse custom scoring matrix.
    def __parse_custom_matrix(self):
        """ 
        A schema is organized such that you have 3 columns: A, B, C.
        Columns A and B represents the query and target base, respectively.
        Column C is the real-number assigned to the mismatch given A and B.
        This file must be tab-delimited.
        Example:
        A    T    -5
        A    C    2
        A    G    -2
        A    A    8
        T    A    2 
        ... 
        etc. 
        """
        submat = {} # the substitution matrix
        for line in open(self.args['custom']): # parse custom matrix file
            line = line.strip()
            if len(line) == 0: # if an empty line, terminate analysis
                break
            else: # pertains to parsing of the custom substition matrix
                line = line.split('\t')
                if len(line) < 3:
                    raise IndexError('Custom matrix must have 3x columns')
                else: # set the Query (A) and Target (B) scores
                    a, b, score = line[0], line[1], float(line[2])
                    submat[(a, b)] = score
        return submat

# Helper-function to write a string
def out(s):
    sys.stdout.write(s+'\n')

# Function to parse a preexisting version of the file that will be written to, returning a list of sequences already aligned
def parse_output(fname):
    alreadyDone = []
    if os.path.isfile(fname):
        for line in open(fname):
            if len(line) == 0:
                break
            else:
                sequenceName = line.split('\t',1)[0]
                if len(sequenceName) > 0:
                    alreadyDone.append(sequenceName)
    return alreadyDone

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
            return 'C'
        if (self.char1 == 'A' and self.char2 == '-') or (self.char2 == 'A' and self.char1 == '-'):
            return 'A'
        if (self.char1 == 'T' and self.char2 == '-') or (self.char2 == 'T' and self.char1 == '-'):
            return 'T'
        if (self.char1 == 'C' and self.char2 == '-') or (self.char2 == 'C' and self.char1 == '-'):
            return 'C'

# Parses and processes the consensus string to ultimately yield a single
# string which encapsulate the pairwise alignment.
class ConsensusProcessor():
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

# Executes multiple-sequence alignment
class MultipleSequenceDriver():
    def __init__(self, queries, input_state):
        self.queries = queries
        self.costs = input_state.get_penalties() # set costs to factory
        self.submat = input_state.get_submatrix() # set submatrix to factory
        self.score_type = input_state.get_scoretype()

    # As of now, MSA alignment is a sequential algorithm.
    def align(self):
        out('--- Multiple sequence alignment mode ---')
        # get the first two input sequences
        s0 = self.queries[0]
        s1 = self.queries[1]
        # pass them both into the tree--based Needleman--Wunsch algorithm.
        nw = factory.NeedlemanWunsch(s1=s0, s2=s1, costs=self.costs, submat=self.submat, nodeTypes=factory.default_nodetypes())
        first_align, second_align = nw.prettify()[1]
        # feed respective alignments into an analysis class and get consensus.
        consensus = ConsensusProcessor(str1=first_align, 
                                       str2=second_align).get_alignment()
        # since the first two sequences have been aligned, focus on all others.
        for i in range(2, len(self.queries)):
            curr_seq = self.queries[i]
            nw = factory.NeedlemanWunsch(s1=consensus, s2=curr_seq, 
                                         costs=self.costs, submat=self.submat, 
                                         nodeTypes=factory.default_nodetypes())
            align_sA, align_sB = nw.prettify()[1]
            consensus = ConsensusProcessor(str1=align_sA, 
                                       str2=align_sB).get_alignment()
        aln_string = NeuriteSequence(name='consensus', sequence=consensus)
        return aln_string

# Executes the pairwise application
class PairwiseDriver():
    def __init__(self, targets, queries, input_state):
        self.targets = targets # factory operates given targets and queries
        self.queries = queries
        self.costs = input_state.get_penalties() # set costs to factory
        self.submat = input_state.get_submatrix() # set submatrix to factory
        self.score_type = input_state.get_scoretype()
        self.forceQuery = input_state.get_args()['forceQuery']

        # Get sequences already completed and remove from queries
        self.priorCompletions = parse_output(input_state.get_args()['o'])
        self.num_complete = len(self.priorCompletions) # for how many sequences have been aligned
        #out(str(self.num_complete)+" complete of "+str(len(targets)))

        # Set openMode to append if some targets have already been run and completed
        if self.num_complete > 0:
            openMode = 'a'
        else:
            openMode = 'w'
            
        self.scorehandle = open(input_state.get_args()['o'], openMode) # output file
        self.alignhandle = None
        if input_state.get_args()['a'] is not '':
            self.alignhandle = open(input_state.get_args()['a'], openMode) # alignments file
            
        self.num_workers = input_state.get_args()['n']
        # Get node type lists
        if input_state.get_args()['nodeTypes'] is None:
            self.nodeTypes = factory.default_nodetypes()
        else:
            self.nodeTypes = factory.parse_nodetypes(input_state.get_args()['nodeTypes'])

    # Initialize the factory given query sequences and input arguments
    def start(self):
        out('--- Pairwise alignment mode ---')
        executor = concurrent.futures.ProcessPoolExecutor(self.num_workers)
        if self.forceQuery:
            queryCompletions = []
        else:
            queryCompletions = self.priorCompletions
        try:
            for target in self.targets: # per fasta, create a concurrent job, f.
                if target.name not in self.priorCompletions:
                    f = executor.submit(pairwise_mapper, target, queries, self.costs, self.submat, self.nodeTypes, queryCompletions)
                    f.add_done_callback(self._callback)
            executor.shutdown()
            self.close_output_buffers()
            out('** Analysis Complete **')
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

    # Determines which score to use and call the appropriate function
    def calc_score(self, result):
        if result[0] is None:
            return None
        elif self.score_type is 'alignment':
            return result[0]
        elif self.score_type is 'gaps':
            return self.count_gaps(result)
        elif self.score_type is 'excess_gaps':
            return self.count_excess_gaps(result)
        elif self.score_type is 'short_normalized':
            return self.calc_short_normalized(result)
        elif self.score_type is 'long_normalized':
            return self.calc_long_normalized(result)

    # Just gets the total number of gaps in the alignment
    def count_gaps(self, result):
        al1 = result[1][0]
        al2 = result[1][1]
        return al1.count('-')+al2.count('-')
        
    # Determine the number of gaps in excess of those required from the length difference between the two sequences
    def count_excess_gaps(self, result):
        al1 = result[1][0]
        al2 = result[1][1]
        gapsIn1 = al1.count('-')
        gapsIn2 = al2.count('-')
        seq1Len = len(al1)-gapsIn1
        seq2Len = len(al2)-gapsIn2
        return gapsIn1+gapsIn2-abs(seq1Len-seq2Len)
        
    # Divides the score by the smaller sequence length
    def calc_short_normalized(self, result):
        seq1Len = len(self.al1)-result[1][0].count('-')
        seq2Len = len(self.al2)-result[1][1].count('-')
        return result[0]/min(seq1Len,seq2Len)
    
    # Divides the score by the larger sequence length    
    def calc_long_normalized(self, result):
        seq1Len = len(self.al1)-result[1][0].count('-')
        seq2Len = len(self.al2)-result[1][1].count('-')
        return result[0]/max(seq1Len,seq2Len)
        
    # Callback function once a thread is complete
    def _callback(self, return_val):
        res = return_val.result() # get result once thread is complete
        target, results = res
        results = sorted(results, key=lambda x: x[-1]) # sort by query (last item)
        if self.num_complete == 0: # for the first result, write headers
            self._create_header(results)

        # save scores to the alignment matrix
        scores = '\t'.join([str(self.calc_score(s)) for s in results])
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
        out(' --> ' + target + ' [OK] '+str(self.num_complete) +
            ' of '+ str(len(self.targets))) # print-out progress

# Maps each query sequence against a set of targets (itself)
def pairwise_mapper(target, queries, costs, submat, nodeTypes, priorCompletions):
    results = [] # K => target, V => aligned queries 
    # get the gap and substitution matrix
    for query in queries:
        # Doesn't run the current query if it has already been run as a target (avoid duplicating effort)
        if query.name not in priorCompletions:
            NW = factory.NeedlemanWunsch(target, query, costs, submat, nodeTypes)
            output = NW.prettify()
            results.append(output)
        else:
            #print(query.name+' already completed')
            results.append([None,None,query.name])
    return target.name, results

if __name__ == '__main__':
    try:
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
        else: # else, start multiple-sequence alignment (MSA)
            driver = MultipleSequenceDriver(queries, input_state)
            aln_string = driver.align()
            print(aln_string.seq)
            for i in queries:
                print(i)

    except (IOError, KeyboardInterrupt, IndexError) as e:
        out(str(e)+'\n')