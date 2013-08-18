''' 
Performs both local and multiple sequence alignment (MSA) given a list of
user-provided input sequences as a FASTA file.
'''

from parameter import AlignmentArgumentValidator, AlignmentCommandParser, InputWrapperState
from pairwise import PairwiseDriver
from msa import MultipleSequenceDriver, ConsensusFilterFactory
        
version = 0.2

def run_local(targets, queries, input_state):
    driver = PairwiseDriver(targets, queries, input_state)
    driver.start() # start only pairwise alignment

def run_msa(queries, input_state):
    driver = MultipleSequenceDriver(queries, input_state)
    driver.build_preconsensus()
    driver.align() # align sequences and build resultant consensus
    consensus_fact = ConsensusFilterFactory(driver.alns, args['thresh'], args['type'])
    consensus_fact.build_consensus()
    consensus_fact.write(fname=args['build']) # write analysis
    
if __name__ == '__main__':
    try:
        args = AlignmentCommandParser().parse_args()
        AlignmentArgumentValidator(args) # test all arguments are correct

        input_state = InputWrapperState(args)
        input_state.assign_matrix() # parse in-built or custom matrix
        targets = input_state.parse_fasta(input_state.fname) # next, parse fasta file
        if input_state.fname2 is None:
            queries = targets
        else:
            queries = input_state.parse_fasta(input_state.fname2)
        
        if args['mode'] == 'local':
            run_local(targets, queries, input_state)
        elif args['mode'] == 'msa': # start multiple-sequence alignment (MSA)
            run_msa(queries, input_state)
            
    except (IOError, KeyboardInterrupt, IndexError) as e:
        print(str(e)+'\n')
