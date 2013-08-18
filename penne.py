''' 
Enables the ability to analyze multiple sequence alignment (MSA) builds and
identify over-represented domains within such analyses.
'''

from parameter import DomainArgumentValidator, DomainCommandParser
from domain import DomainSetBuilder
from buffer import XMLBuildReader, status_message

version = 0.1
    
if __name__ == '__main__':
    try:
        args = DomainCommandParser().parse_args()
        DomainArgumentValidator(args) # test all arguments are correct
        print('penne - v.' + str(version) + '\n=============')
        cons_query = XMLBuildReader(args['query']).parse()
        cons_baseline = XMLBuildReader(args['baseline']).parse()
        
        # next, yield domains for both query and baseline datasets. 
        doms_query = DomainSetBuilder(win=args['win'], max_gap=args['max_g'], 
                         is_enum=args['enumerate'], consensus=cons_query)
        doms_baseline = DomainSetBuilder(win=args['win'], max_gap=args['max_g'], 
                         is_enum=args['enumerate'], consensus=cons_baseline)
        counts_query = doms_query.build() # build abundance counts
        counts_baseline = doms_baseline.build()
        status_message('Abundance counting', 'OK')
        
    except (IOError, KeyboardInterrupt, IndexError) as e:
        print(str(e))
