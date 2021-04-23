'''
Created on 9 Oct 2018
Updated 26 Feb 2021

@author: thomasgumbricht

The smap package and module links to geomimage.assets for hadling all DAAC data holdings
'''

# Standard library imports


# Package application imports

from geoimagine.assets import AccessOnlineData

class ProcessSmap(AccessOnlineData):
    '''class for SMAP specific processing
    '''   
    
    def __init__(self, pp, session):
        ''' Expects the pp processing instructions and an open postgres database session
        '''
        
        # Initiate the package for Online data access
        AccessOnlineData.__init__(self)
        
        self.session = session
                
        self.pp = pp  
        
        self.verbose = self.pp.process.verbose 
        
        self.session._SetVerbosity(self.verbose)

        print ('        ProcessSmap',self.pp.process.processid) 

        # Direct to sub-processes
        
        if self.pp.process.processid.lower() == 'searchsmapproducts':
            
            self._SearchOnlineProducts()
            
        elif self.pp.process.processid.lower() == 'smapsearchtodb':
            
            self._SearchToDB()
            
        elif self.pp.process.processid.lower() == 'downloadsmapdaac':
            
            self._DownLoadProduct()
            
        elif self.pp.process.processid.lower() == 'extractsmaphdf':
            
            self._ExtractHdf()
            
        elif self.pp.process.processid.lower() == 'checksmap':
            
            self._CheckSmap()
            
        else:
            
            exitstr = 'Exiting, processid %(p)s missing in ProcessSmap' %{'p':self.pp.process.processid}
            
            exit(exitstr)
            
            