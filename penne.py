''' 
Enables the ability to analyze multiple sequence alignment (MSA) builds and
identify over-represented domains within such analyses.
'''

from parameter import DomainArgumentValidator, DomainCommandParser

version = 0.1
    
if __name__ == '__main__':
    try:
        args = DomainCommandParser().parse_args()
        DomainArgumentValidator(args) # test all arguments are correct
            
    except (IOError, KeyboardInterrupt, IndexError) as e:
        print(str(e))
