import argparse
import platform
from Bio.SubsMat import MatrixInfo
from Bio import SeqIO
from sequence import NeuriteSequence
import sequence
from random import shuffle

class DomainArgumentValidator():
    ''' 
    Validates that user-provided arguments relative to domain analysis are
    sound and of valid ranges.
    '''
    def __init__(self, args):
        self.args = args
        self.check_args()
    
    def test_laplace(self):
        if 0 <= self.args['l'] <= 1.0:
            return True
        else:
            raise IOError('Laplace correction must be a fraction [error]')
    
    def test_window_size(self):
        if self.args['win'] > 0:
            return True
        else:
            raise IOError('Sliding window must be positive [error]')
    
    def test_max_gaps(self):
        if self.args['max_g'] > 0:
            return True
        else:
            raise IOError('Max #/gaps must be positive [error]')
    
    def check_args(self):
        return all([self.test_laplace(), self.test_max_gaps(), 
                    self.test_window_size()])

# Validates user-provided command-line arguments relative to sequence alignment
class AlignmentArgumentValidator():
    def __init__(self, args):
        self.args = args
        self.check_python()
        self.check_args()

    # Check python 3.0+ is available
    def check_python(self):
        if not platform.python_version_tuple()[0:2] >= ('3', '0'): # >= 3.0
            raise RuntimeError('Python 3.0+ recommended')

    # Checks user-provided arguments are valid
    def check_args(self):
        return all([self.test_num_workers(),
                self.test_valid_matrix(), self.test_threshold()])

    # Test a valid substitution matrix is selected
    def test_valid_matrix(self):
        #all_matrices = set(MatrixInfo.available_matrices) # al sub. matrices
        if self.args['custom']:
            return True
        else:
            err = 'A custom matrix is required.'
            raise IOError(err)

    # Test a valid number of workers are provided
    def test_num_workers(self):
        if self.args['n'] >= 1:
            return True
        else:
            raise IOError('>= 1 worker processes must be provided')
    
    # Test that a valid consensus threshold is provided
    def test_threshold(self):
        if self.args['thresh'] >= 0:
            return True
        else:
            raise IOError('Consensus threshold (thresh) must be positive')

# Helper-class to parse input arguments
class AlignmentCommandParser():
    def __init__(self):
        desc = 'Pairwise and multiple sequence alignment (MSA) of neurite sequences.'
        u='%(prog)s [options]' # command-line usage
        self.parser = argparse.ArgumentParser(description=desc, add_help=False, usage=u)
        self._init_params()

    # Create parameters to be used throughout the application
    def _init_params(self):
        param_aln = self.parser.add_argument_group('Required alignment parameters')
        param_local = self.parser.add_argument_group('Pairwise (local) parameters')
        param_msa = self.parser.add_argument_group('MSA parameters')
        param_opts = self.parser.add_argument_group('Optional parameters')

        # Specify required arguments
        param_aln.add_argument('-f', metavar='FILE',
                    help='Input FASTA file [na]')
        
        param_aln.add_argument('-mode', metavar='MODE',
                    choices=['local', 'msa', 'domain'],
                    help='Analysis mode {local, msa} [msa]')
        
#        param_aln.add_argument('-matrix', metavar='STR', default=None,
#                    help='Matrix; see Biopython MatrixInfo [na]')
        
        # Local-alignment specific parameters
        param_local.add_argument('-gap', metavar='INT', default=-1, type=int,
                    help='Gap extension penalty [-1]')
        
        param_local.add_argument('-gapopen', metavar='INT', default=0, type=int,
                    help='Gap open penalty [0]')
        
        # MSA specific parameters
        param_msa.add_argument('-thresh', metavar='FLOAT', default=0.7, type=float, 
                    help='Consensus threshold [0.7]')
        
        param_msa.add_argument('-type', metavar='STR', default='percent', 
                    choices=['percent', 'sqrt'],
                    help='Threshold type {percent, sqrt} [percent]')
        
        param_msa.add_argument('-build', metavar='FILE', default='alns.xml', 
                    type=str, help='Output file of consensus & alignments [./alns.xml]')

        param_msa.add_argument('-iterate', metavar='FLOAT', default=1, type=float,
                    help='Number of MSA iterations (using a PWM) or threshold for change in 40% composite score to continue iterating [1]')
        
        param_opts.add_argument('-subsample', metavar='FLOAT', default=1, type=float, 
                    help='Subsample of data, taking the first n sequences. Value treated as proportion of total if (0,1] and explicit number for [2,N].')

        param_opts.add_argument('-subsample_start', metavar='FLOAT', default=0, type=float, 
                    help='If taking a subsample, subsample_start determines from which sequence to start the subset (as proportion or explicit number, indexed from 0). If size and start of subsample lead to an exhaustion of the available sequences an error will be be displayed.')
        
        param_opts.add_argument('--random_subset', action='store_const', const=True, default=False,
                    help='Subset is a random sample of data (shuffle sequences before taking subset).')

        param_opts.add_argument('--random_order', action='store_const', const=True, default=False,
                    help='Order of sequences (in subsample of data if using --subset) is shuffled.')

        # General, optional arguments
        param_opts.add_argument('-f2', metavar='FILE', default=None,
                    help='Second input FASTA file [na]')
        
        param_opts.add_argument('-custom', metavar='FILE', default=None,
                    help='Custom substitution matrix [na]')
        
        param_opts.add_argument('-node_types', metavar='FILE', default=None,
                    help='Node Type Specifications [na]')
        
        param_opts.add_argument('-n', metavar='INT', default=2, type=int,
                    help='Number of worker processes [2]')
        
        param_opts.add_argument('-o', metavar='FILE', default='./scores.tab', 
                    help='File to write/append output [./scores.tab]')
        
        param_opts.add_argument('-a', metavar='FILE', default=None, 
                    help='File to write/append alignments [na]')
        
        param_opts.add_argument('-h','--help', action='help',
                    help='Show this help screen and exit')

    # Get the arguments for each parameter
    def parse_args(self):
        return vars(self.parser.parse_args()) # parse arguments

class DomainCommandParser():
    def __init__(self):
        desc = 'Domain analysis of neurite multiple sequence alignments.'
        u='%(prog)s [options]' # command-line usage
        self.parser = argparse.ArgumentParser(description=desc, add_help=False, usage=u)
        self._init_params()

    def _init_params(self):
        param_reqd = self.parser.add_argument_group('Required parameters')
        param_msa = self.parser.add_argument_group('MSA parameters')
        param_domain = self.parser.add_argument_group('Domain parameters')
        param_opts = self.parser.add_argument_group('Optional parameters')

        # Domain-specific parameters
        param_reqd.add_argument('-mode', metavar='MODE',
                    choices=['single','multiple'],
                    help='Analysis mode {single, multiple} [single]')

        param_reqd.add_argument('-query', metavar='BUILD', 
                    help='Set query consensus build (if mode=single) or fasta file (if mode=multiple) [na]')
        
        param_reqd.add_argument('-baseline', metavar='BUILD',
                    help='Set baseline consensus build (if mode=single) or fasta file (if mode=multiple) [na]')
        
        param_domain.add_argument('-max-g', metavar='INT', default=3, type=int, 
                    help='Maximum #/gaps in domain [3]')
        
        param_domain.add_argument('-win', metavar='INT', default=7, type=int, 
                    help='Sliding window length; characters [7]')
        
        param_domain.add_argument('-l', metavar='FLOAT', default=0.4, type=float, 
                    help='LaPlace correction cutoff [0.4]')
        
        param_domain.add_argument('-o', metavar='FILE', default='./doms.txt', type=str, 
                    help='Output file for saving identified domains [./doms.txt]')
        
        param_domain.add_argument('-p', metavar='FLOAT', default=0.05, type=float, 
                    help='Hypergeometric p-value cutoff [0.05]')
        
        param_domain.add_argument('--ipf', action='store_const', const=True, 
                    default=False, help = 'IPF normalize [false]')
        
        param_domain.add_argument('--enumerate', action='store_const', const=True, 
                    default=False, help = 'Use many window & gap cutoffs [false]')

        param_domain.add_argument('-minwin', metavar='INT', default=1, type=int, 
                    help = 'While using enumerate, the minimum window size [1]')
        
        param_domain.add_argument('--strip', action='store_const', const=True,
                    default=False, help = 'Strip gaps within sliding window [false]')

        param_domain.add_argument('-runs', metavar='INT', default=10, type=int, 
                    help='Number of times to run MSA if mode=multiple [10]')

        param_domain.add_argument('--cluster', action='store_const', const=True,
                    default=False, help='Cluster domains before running statistical analysis [false]')

#        param_domain.add_argument('-max_cluster_dist', metavar='INT', default=0, type=int, 
#                    help='Maximum distance between seed domain and any domain in the cluster [0]')

        param_domain.add_argument('-thresh', metavar='FLOAT', default=0.7, type=float, 
                    help='Consensus threshold [0.7]')
        
        param_domain.add_argument('-type', metavar='STR', default='percent', 
                    choices=['percent', 'sqrt'],
                    help='Threshold type {percent, sqrt} [percent]')

        # MSA specific parameters
        param_msa.add_argument('-node_types', metavar='FILE', default=None,
                    help='Node Type Specifications [na]')

        param_msa.add_argument('-custom', metavar='FILE', default=None,
                    help='Custom substitution matrix [na]')

        param_msa.add_argument('-gap', metavar='INT', default=-1, type=int,
                    help='Gap extension penalty [-1]')
        
        param_msa.add_argument('-gapopen', metavar='INT', default=0, type=int,
                    help='Gap open penalty [0]')

        param_msa.add_argument('-iterate', metavar='FLOAT', default=1, type=float,
                    help='Number of MSA iterations (using a PWM) or threshold for change in 40% composite score to continue iterating [1]')
        
        param_opts.add_argument('--overlap', action='store_const', const=True, default=False,
                    help='Allow query and baseline sets to contain the same sequences. If false, remove overlapping sequences from baseline set [False]')

        param_opts.add_argument('--disjoint_subset', action='store_const', const=True, default=False,
                    help='Ensure that baseline subset contains no sequences from query subset [False]')
        
        param_opts.add_argument('-subsample', metavar='FLOAT', default=1, type=float, 
                    help='Subsample of data, taking the first n sequences. Value treated as proportion of total if (0,1] and explicit number for [2,N]')

        param_opts.add_argument('-subsample_start', metavar='FLOAT', default=0, type=float, 
                    help='If taking a subsample, subsample_start determines from which sequence to start the subset (as proportion or explicit number, indexed from 0). If size and start of subsample lead to an exhaustion of the available sequences an error will be be displayed.')
        
        param_opts.add_argument('--random_subset', action='store_const', const=True, default=False,
                    help='Subset is a random sample of data (shuffle sequences before taking subset) [False]')

        param_opts.add_argument('--random_order', action='store_const', const=True, default=False,
                    help='Order of sequences (in subsample of data if using --subset) is shuffled [False]')

        param_opts.add_argument('-n', metavar='INT', default=2, type=int,
                    help='Number of worker processes [2]')

        param_opts.add_argument('-h','--help', action='help',
                    help='Show this help screen and exit')
    
    def parse_args(self):
        return vars(self.parser.parse_args()) # parse arguments

class ConsensusStatsCommandParser():
    def __init__(self):
        desc = 'Statistics of neurite multiple sequence alignments.'
        u='%(prog)s [options]' # command-line usage
        self.parser = argparse.ArgumentParser(description=desc, add_help=False, usage=u)
        self._init_params()

    def _init_params(self):
        param_reqd = self.parser.add_argument_group('Required parameters')
        param_opts = self.parser.add_argument_group('Optional parameters')

        # Domain-specific parameters
#        param_reqd.add_argument('-mode', metavar='MODE',
#                    choices=['?','?'],
#                    help='Analysis mode {?}')

        param_reqd.add_argument('-build', metavar='BUILD', 
                    help='XML build file containing composite, consensus, and MSA [na]')

        param_reqd.add_argument('-newick', metavar='BUILD', default=None,
                    help='Newick file for output build file [None]')
        
        param_opts.add_argument('-h','--help', action='help',
                    help='Show this help screen and exit')
    
    def parse_args(self):
        return vars(self.parser.parse_args()) # parse arguments

# An InputStateWrapper wraps data used as input, eg. FASTA file,
# custom matrix, user-provided arguments, etc.
class InputWrapperState():
    def __init__(self, args):
        self.args = args # reference user-provided arguments
        self.subsmat = None # references data for substitution matrix
        self.fname = args['f'] # input filename
        self.fname2 = args['f2'] # input filename
        self.alignment_file = args['a']
        if 'node_types' in args.keys() and args['node_types'] is not None:
            self.node_types = sequence.parse_node_types(args['node_types'])
        else:
            self.node_types = sequence.default_nodetypes
        self.assign_matrix()

    # Get arguments
    def get_args(self):
        return self.args

    # Get node types
    def get_node_types(self):
        return self.node_types       

    # Get arguments relative to penalties
    def get_penalties(self):
        cost_ids = ('gap','gapopen') # all possible costs, might in the future include a separate gap open and gap extension cost
        return {k: self.args[k] for k in cost_ids} # get costs per penality

    # Get the substitution matrix which will be used
    def get_submatrix(self):
        return self.subsmat
        
    # Trivial function to parse a fasta file
    def parse_fasta(self, fname):
        parsed_fasta = list(SeqIO.parse(fname, 'fasta')) # easy indexing
        queries = [] # references list of parsed sequences
        for i in parsed_fasta: # cast as a neuronal sequence; easy modeling.
            s = NeuriteSequence(seq=str(i.seq), name=i.name)
            queries.append(s)
        print(str(len(queries)) + ' queries parsed [OK]')
        return queries # return set of fasta entries

    # Set the desired matrix the user wishes to add
    def assign_matrix(self):
        if 'custom' not in self.args or not self.args['custom']: # if no custom matrix, use identity
            self.subsmat = sequence.generate_identity_matrix(self.node_types)
        else:
            self.subsmat = self.__parse_custom_matrix() # custom matrix

    # Get the user-provided in-built matrix
    #def __parse_inbuilt_matrix(self):
    #    return getattr(MatrixInfo, self.args['matrix']) # get substitution matrix

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

def generate_sequence_set(sequences,num_seqs=0,random_set=False,random_order=True,seq_start=0,disallowed=[]):
    sequence_set = sequences
    if len(disallowed) > 0:
        # Remove any disallowed sequences from the sequence set
        sequence_set = list([sequence for sequence in sequences if sequence.name not in list([dis_seq.name for dis_seq in disallowed])])

#    print("Num_seqs: "+str(num_seqs)+"; RandomSet: "+str(random_set)+"; RandomOrder: "+str(random_order))

    if num_seqs != 1:
        if num_seqs < 1:
            num_seqs = int(len(sequences)*num_seqs)
        else:
            num_seqs = int(min(len(sequences),num_seqs))

        if seq_start < 1:
            start = int(seq_start*len(sequences))
        else:
            start = int(seq_start)

        if random_set:
            ids = [x for x in range(len(sequences))]
            shuffle(ids)
            ids = sorted(ids[0:num_seqs])
            sequence_set = list([sequences[i] for i in ids])
        else:
            sequence_set = list([sequences[i] for i in range(start,start+num_seqs)])

    if random_order:
        shuffle(sequence_set)

    return sequence_set
