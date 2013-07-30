version = 0.1

if __name__ == '__main__':
    try:
        from parameter import ArgumentValidator, CommandLineParser, InputWrapperState
        from pairwise import PairwiseDriver
        from msa import MultipleSequenceDriver, ConsensusFilterFactory
        
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
        elif args['mode'] == 'msa': # start multiple-sequence alignment (MSA)
            driver = MultipleSequenceDriver(queries, input_state)
            driver.build_preconsensus()
            # map queries back onto consensus and build a filtered consensus
            alignments = driver.align()
            consensus_fact = ConsensusFilterFactory(alignments, args['thresh'], 
                                                    args['threshold_type'])
            consensus_fact.build_consensus()
    except (IOError, KeyboardInterrupt, IndexError) as e:
        print(str(e)+'\n')
