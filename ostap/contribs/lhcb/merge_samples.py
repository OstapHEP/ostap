#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file
# Helper function to "merge" sub-samples for LHCb Run 1&2 data 
# =============================================================================
""" Helper function to "merge" sub-samples for LHCb Run 1&2 data 
"""
# =============================================================================\
__all__ = (
    'merge_samples' , ## Helper function to "merge" sub-samples for LHcb data 
    )
# =============================================================================
import ostap.contribs.lhcb.data 
import ostap.trees.data_utils
import ostap.io.files 
# =============================================================================
from ostap.logger.logger import getLogger
logger = getLogger ( 'ostap.contrib.lhcb.merge_samples')
# =============================================================================
## Helper function to "merge" sub-samples for LHcb data
#  @code
#  data   = ... ## dictionary for file/tree/data collections
#  merged = merge_samples ( data ) 
#  @endcodce 
def merge_samples ( data ) :
    """ Helper function to "merge" sub-samples for LHcb data 
    >>> data = ... ## dictionary for file/tree/data collections
    >>> merged = merge_samples ( data ) 
    """
    ## samples by year 
    if not '2011'      in data and '2011u' in data and '2011d' in data :
        data [ '2011'     ] = data [ '2011u'    ] + data [ '2011d'    ]
        
    if not '2012'      in data and '2012u' in data and '2012d' in data :
        data [ '2012'     ] = data [ '2012u'    ] + data [ '2012d'    ]

    if not '2015'      in data and '2015u' in data and '2015d' in data :
        data [ '2015'     ] = data [ '2015u'    ] + data [ '2015d'    ]

    if not '2016'      in data and '2016u' in data and '2016d' in data :
        data [ '2016'     ] = data [ '2016u'    ] + data [ '2016d'    ]

    if not '2017'      in data and '2017u' in data and '2017d' in data :
        data [ '2017'     ] = data [ '2017u'    ] + data [ '2017d'    ]

    if not '2018'      in data and '2018u' in data and '2018d' in data :
        data [ '2018'     ] = data [ '2018u'    ] + data [ '2018d'    ]
                
    if not 'Run1'      in data and '2011'  in data and '2012'  in data :
        data [ 'Run1'     ] = data [ '2011'     ] + data [ '2012'     ]

    if not 'Run1Up'    in data and '2011u' in data and '2012u' in data :
        data [ 'Run1Up'   ] = data [ '2011u'    ] + data [ '2012u'    ]

    if not 'Run1Down'  in data and '2011d' in data and '2012d' in data :
        data [ 'Run1Down' ] = data [ '2011d'    ] + data [ '2012d'    ]

    if not 'Run2'      in data and ( '2015' in data and
                                     '2016' in data and
                                     '2017' in data and
                                     '2018' in data ) :
        data [ 'Run2' ] = data [ '2015' ] + data [ '2016' ] + data [ '2017' ] + data [ '2018' ]

    if not 'Run2Up'    in data and ( '2015u' in data and
                                     '2016u' in data and
                                     '2017u' in data and
                                     '2018u' in data ) :
        data [ 'Run2Up' ] = data [ '2015u' ] + data [ '2016u' ] + data [ '2017u' ] + data [ '2018u' ]

    if not 'Run2Down'  in data and ( '2015d' in data and
                                     '2016d' in data and
                                     '2017d' in data and
                                     '2018d' in data ) :
        data [ 'Run2Down' ] = data [ '2015d' ] + data [ '2016d' ] + data [ '2017d' ] + data [ '2018d' ]
                              
    if not 'MagUp'    in data and 'Run1Up' in data and 'Run2Up' in data :
        data [ 'MagUp'    ] = data [ 'Run1Up'   ] + data [ 'Run2Up'   ]

    if not 'MagDown'  in data and 'Run1Down' in data and 'Run2Down' in data :
        data [ 'MagDown'  ] = data [ 'Run1Down' ] + data [ 'Run2Down' ]

    if not 'All'      in data and 'Run1' in data and 'Run2' in data :
        data [ 'All'      ] = data [ 'Run1'     ] + data [ 'Run2'     ]

    if not 'ALL'      in data and 'Run1' in data and 'Run2' in data :
        data [ 'ALL'      ] = data [ 'Run1'     ] + data [ 'Run2'     ]
        
    return data 

# =============================================================================
if '__main__' == __name__ :

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
##                                                                      The END 
# =============================================================================

