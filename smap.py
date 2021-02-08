'''
Created on 9 Oct 2018

@author: thomasgumbricht
'''
import urllib.request
from html.parser import HTMLParser
import os
from sys import exit
from shutil import move
import support.karttur_dt as mj_dt

from params import Composition, LayerCommon, RegionLayer, VectorLayer, RasterLayer
#import geoimagine.smap.hdf5_2_geotiff as hdf5_2_geotiff
from smap import hdf5_2_geotiff 

class SmapComposition:
    '''
    class for sentinel compositions
    '''
    def __init__(self, compD):  
        for key in compD:
            if '_' in compD[key]:
                exitstr = 'the "%s" parameter can not contain underscore (_): %s ' %(key, compD[key])
                exit(exitstr) 
            setattr(self, key, compD[key])
        if not hasattr(self, 'content'):
            exitstr = 'All SMAP compositions must contain a content'
            exit(exitstr)
            
class SmapTile(LayerCommon):
    '''Class for sentinel tiles'''
    def __init__(self, smapid, composition, locusD, datumD, filepath, FN): 
        """The constructor expects an instance of the composition class."""
        LayerCommon.__init__(self)
        self.smapid = smapid
        self.comp = composition
        
        self.locus = locusD['locus']
        self.locuspath = locusD['path']

        self.path = filepath
        self.FN = FN

        self.datum = lambda: None
        for key, value in datumD.items():
            setattr(self.datum, key, value)
        if self.datum.acqdate:
            self._SetDOY()
            self._SetAcqdateDOY()
        self._SetPath()
        self._SetQuery()
        
    def _SetPath(self):
        """Sets the complete path to sentinel tiles"""
        
        self.FP = os.path.join('/Volumes',self.path.volume, self.comp.system, self.comp.source, self.comp.division, self.comp.content, self.locuspath, self.datum.acqdatestr)
        self.FPN = os.path.join(self.FP,self.FN)
        if ' ' in self.FPN:
            exitstr = 'EXITING smap FPN contains space %s' %(self.FPN)
            exit(exitstr)
            
    def _SetQuery(self):
        '''
        '''
        self.query = {'smapid':self.smapid, 'smapfilename':self.FN,'source':self.comp.source,'product':self.comp.product,
                 'version':self.comp.version,'acqdate':self.datum.acqdate, 'doy':self.datum.doy, 'content':self.comp.content}

class ProcessSmap:
    '''class for SMAP specific processing
    '''   
    
    def __init__(self, pp, session):
        '''
        '''
        
        self.session = session
                
        self.pp = pp  
        
        self.verbose = self.pp.process.verbose 
        
        self.session._SetVerbosity(self.verbose)

        if self.verbose > 0:

            print ('        ProcessSmap',self.pp.process.processid) 

        #Direct to SMAP sub-processes
        if self.pp.process.processid.lower() == 'searchsmapproducts':
            
            self._SearchSmapProducts()
            
        elif self.pp.process.processid.lower() == 'smapsearchtodb':
            
            self._SmapSearchToDb()
            
        elif self.pp.process.processid.lower() == 'downloadsmapdaac':
            
            self._DownLoadSmapDaac()
            
        elif self.pp.process.processid.lower() == 'extractsmaphdf':
            
            self._ExtractSmapHdf()
            
        elif self.pp.process.processid.lower() == 'checksmap':
            
            self._CheckSmap()
            
        else:
            exitstr = 'Exiting, processid %(p)s missing in ProcessSmap' %{'p':self.pp.process.processid}
            exit(exitstr)
            
            
    def _SmapPath(self):
        ''' Set paths to SMAP products
        '''
        
        # get the product, note the conflict between kartturs and smap namining convention
        self.product = self.pp.process.parameters.product.replace('-','_')
        
        self.version =self.pp.process.parameters.version
        
        #check that the version is correctly stated
        
        if not len(self.version) == 3:
            
            exit('The smap version must be 3 digits, e.g. "005" or "006"')
            
        if not self.version.isdigit():
            
            exit('The smap version must be 3 digits, e.g. "005" or "006"')
            
        self.prodPath ='%s.%s' %(self.product,self.version)
        
        if self.verbose > 1:
            
            print ('    Processing SMAP product.', self.prodPath)
            
        #create the localpath where the search data (html) will be saved
        self.localPath = os.path.join('/volumes',self.pp.dstPath.volume,'DAAC-SMAP',self.prodPath)
        
        if not os.path.exists(self.localPath):
            
            os.makedirs(self.localPath)
        
 
    def _SearchSmapProducts(self):
        '''IMPORTANT the user credentials must be in a hidden file in users home directory called ".netrc"
        '''
        
        # Set the path for this SMAP product
        self._SmapPath()
        
        #Set todays date 
        today = mj_dt.Today()
        
        #Set the serverurl, pand the SMAP roduct and version to search for 
        self.serverurl = self.pp.process.parameters.serverurl
           
        #Set the sensorpath on the server   
        sensorurl = 'SMAP'
                    
        # Check the existence of aready processed data
        doneFP = os.path.join(self.localPath,'done')
        
        done = True
            
        if not os.path.exists(doneFP):
            
            done = False
            
        #change to the local directory
        
        cmd ='cd %s;' %(self.localPath)
        
        os.system(cmd)
        
        #Loop over the dates defined in process
        
        for datum in self.pp.srcPeriod.datumD:
            
            if self.verbose > 1:
            
                print ('         date:', datum)
            
            # Search for the data
            
            if self.pp.srcPeriod.datumD[datum]['acqdate'] > today:
                
                #skip all dates later than today (there can be no images from the future)
                continue
            
            # Convert date to pointed string used on the server
            dateStr = mj_dt.DateToStrPointDate( self.pp.srcPeriod.datumD[datum]['acqdate'] )
            
            htmlFN = '%s.html' %(dateStr)
            
            # Set the local file name path
            localFPN = os.path.join(self.localPath,htmlFN)
            
            # Define the complete url to the SMAP data
            url = os.path.join(self.serverurl,sensorurl,self.prodPath,dateStr)

            if os.path.exists(localFPN) and not self.pp.process.overwrite:
                
                continue
            
            if done and os.path.exists( os.path.join(doneFP,htmlFN) ) and not self.pp.process.overwrite:
                
                continue
            
            indexFPN = os.path.join(self.localPath,'index.html')
    
            cmd ='cd %s;' %(self.localPath)
            
            #Run the wget command including definition of the cookie needed for accessing the server
            
            #cmd ='%(cmd)s /usr/local/bin/wget -L --load-cookies --spider --no-parent ~/.cookies --save-cookies ~/.cookies %(url)s' %{'cmd':cmd, 'url':url}
            
            # Updated in Feb 2021
            cmd ='%(cmd)s /usr/local/bin/wget -L --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --auth-no-challenge=on --keep-session-cookies --content-disposition %(url)s' %{'cmd':cmd, 'url':url}

            # The online recommended code for single url
            # wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --auth-no-challenge=on --keep-session-cookies --content-disposition <url>
            #print (cmd)
            
            os.system(cmd)
            
            #rename the retrieved index.html to the designated filename
            
            cmd = 'mv %s %s' %(indexFPN, localFPN)
            
            os.system(cmd)
          
    def _SmapSearchToDb(self):
        '''Load search holdings to local db
            Does not utilize the layer class but take parameters directly from xml
        '''
        
        # Set the path for this SMAP product
        self._SmapPath()
        
        #Set todays date
        today = mj_dt.Today()
                 
        #Create a sub-folder called done, when the search results are transferred to the db the html will be moved into the done folder
        doneFP = os.path.join(self.localPath,'done')
            
        if not os.path.exists(doneFP):
            
            os.makedirs(doneFP)
        
        #Loop over the dates
        for datum in self.pp.srcPeriod.datumD:
            
            if self.pp.srcPeriod.datumD[datum]['acqdate'] > today:
                #skip all dates later than today (there can be no images from the future)
                continue
                        
            # Convert date to pointed string used on the server
            dateStr = mj_dt.DateToStrPointDate( self.pp.srcPeriod.datumD[datum]['acqdate'] )
            
            htmlFN = '%s.html' %(dateStr)
            
            # Set the local path
            localFPN = os.path.join(self.localPath,htmlFN)
            
            #Create a sub-folder called done, when the search results are transferred to the db the html will be moved into the done folder
            dstFPN = os.path.join(doneFP,htmlFN)
                       
            if os.path.exists(localFPN):  
                  
                self._ReadSMAPhtml(self.session,localFPN)
                
                move(localFPN, dstFPN)
    
    def _IniTileDownload(self,statusD):
        '''
        '''
        self.dlL = []
        
        # Create a temp folder to which the download will be directed, only when the download is complete will the data be moved in place
        if not os.path.exists(self.tempFP):
        
            os.makedirs(self.tempFP)
            
        # If asscript, the whole downloading will be written as a shell script
        if self.pp.process.parameters.asscript:
            
            shFN = 'download_%(prod)s.sh' %{'prod':self.pp.process.parameters.product}
            
            shFP = os.path.join(self.tempFP, 'script')
            
            if not os.path.exists(shFP):
            
                os.makedirs(shFP)
            
            self.downloadShFPN = os.path.join(shFP,shFN)
            
            self.downloadScriptF = open(self.downloadShFPN,'w')
           
            cmd = 'mkdir -p %(fp)s;\n' %{'fp':shFP}
            
            self.downloadScriptF.write(cmd)
        
        #Get the tiles
        tiles = self.session._SelectSmapData(self.pp.srcPeriod, self.pp.process.parameters, statusD)
        
        return tiles
    
    def _CheckSmap(self):
        #Set the expected layers and parameters for filling the db
        queryD = {}
        queryD['product'] = {'val':self.pp.process.parameters.product, 'op':'=' }
        queryD['retrieve'] = {'val':'Y', 'op':'=' }
        self.paramL = ['source', 'product', 'content', 'layerid', 'prefix', 'suffix', 'celltype', 'dataunit', 'scalefac', 'offsetadd', 'cellnull', 'measure', 'retrieve', 'hdffolder', 'hdfgrid']
        self.compL = ['source', 'product', 'content', 'layerid', 'prefix', 'suffix', 'celltype', 'dataunit', 'scalefac', 'offsetadd', 'cellnull', 'measure']
        self.extractL = self.session._SelectSMAPTemplate( queryD, self.paramL )

        #from geoimagine.support.modis import DisentangleModisTileName as convFN 
        #First loop over the src folder structure to find all tiles at this position
        #Construct a dummy tile, to get the FP
        smapid = 'smapid' 
        hdfFN = '*.%(hdr)s' % {'hdr':self.pp.srcPath.hdrfiletype}

        product = self.pp.process.parameters.product
        version = self.pp.process.parameters.version
        source = '%(p)s.%(v)s' %{'p':product,'v':version}

        acqdate = mj_dt.Today()

        tile = (hdfFN, smapid, source, product, version, 'original', acqdate)

        smapTile = self._ConstructDaacTile(tile,self.pp.srcPath)
        datepath = os.path.split(smapTile.FPN)[0]
        locuspath = os.path.split(datepath)[0]

        for root, directories, filenames in os.walk(locuspath):
            for filename in filenames:
                
                if filename.endswith(self.pp.srcPath.hdrfiletype):

                    queryD = {'smapfilename':filename}
                    paramL = ['smapid', 'smapfilename', 'source', 'product', 'version', 'acqdate']
                    tile = self.session._SelectSingleSMAPDaacTile(queryD,paramL)
                    smapid, smapfilename, source, product, version, acqdate = tile
                    tile = (smapfilename, smapid, source, product, version, 'original', acqdate)
                    smapTile = self._ConstructDaacTile(tile,self.pp.srcPath)

                    #Replace is needed for adjusting between SMAP and Karttur default naming conventions
                    source = source.replace('-E','_E')
                    source = source.replace('-S','_S')
                    if os.path.exists(smapTile.FPN): 
                        self.session._InsertSmapData(smapTile.query)
                        statusD = {'smapid': smapid,'column':'downloaded', 'status': 'Y'}
                        self.session._UpdateSmapStatus(statusD)
 
                        #Only tiles found on file are checked, should it be updated
                        self._SearchExtractLayers(acqdate)

                        if self.nrExploded == len(self.extractL):
                            statusD = {'smapid': smapid,'column':'organized', 'status': 'Y'}
                            self.session._UpdateSmapStatus(statusD)
                            statusD = {'smapid': smapid,'column':'exploded', 'status': 'Y'}
                            self.session._UpdateSmapStatus(statusD)
                    else:
                        pass
                        #This should not happen  
                                                 
    def _SearchExtractLayers(self,acqdate):
        '''Search for extracted layers for specific SMAP tile
        '''
        self.nrExploded = 0
        # self.explodeD is not used
        self.explodeD = {}
        for extcomp in self.extractL:
            paramD = dict(zip(self.paramL,extcomp))
            compD = dict(zip(self.compL,extcomp))
                        
            comp = Composition(compD, self.process.system.dstsystem, self.process.system.dstdivision)
            #Set the datum
            
            acqdatestr = mj_dt.DateToStrDate(acqdate)

            datumD = {'acqdatestr': acqdatestr, 'acqdate':acqdate}

            #Construct the locus dictionary
            locusD = {'locus':'global','path':'global'}
            filepath = lambda: None
            filepath.volume = self.pp.dstPath.volume; filepath.hdrfiletype = self.pp.dstPath.hdrfiletype
            
            #Create a standard raster layer
            layer = RasterLayer(comp, locusD, datumD, filepath)

            if not layer._Exists() or self.process.overwrite:
                self.explodeD[paramD['layerid']] = {'layer':layer,'params':paramD}
            elif layer._Exists():
                self.session._InsertLayer(layer,self.process.overwrite,self.process.delete)
                self.nrExploded += 1                       
                                                
    def _DownLoadSmapDaac(self):
        '''
        '''
        
        #create a temp folder to which the download will be directed, only when the download is complete will the data be moved in place  
        self.tempFP = os.path.join('/Volumes',self.pp.dstPath.volume, 'smap', 'temp')
        
        statusD = {}
        
        # TGTODO downloaded must be in xml, defaulted to N and not obligatory
        statusD['downloaded'] = self.pp.process.parameters.downloaded
        
        #tiles = self.session._SelectSmapData(self.process.srcperiod, self.pp.process.parameters, statusD)
        tiles = self._IniTileDownload(statusD)
        
        for tile in tiles:
            
            self._AddDownload(tile,self.pp.dstPath)  
                     
        self._AccessSMAP()
        
        if self.pp.process.parameters.asscript:
            
            self.downloadScriptF.close()
            
            infostr = 'SMAP dnwnload script file:{n    %s' %(self.downloadShFPN)
            
            print (infostr)
            
    def _AddDownload(self,tile,sdpath):
        '''
        '''
        
        smapTile = self._ConstructDaacTile(tile,sdpath)
        
        smapfilename, smapid, source, product, version, content, acqdate = tile
        
        # replace the hyphen to an underscore 
        source = source.replace('-E','_E')
        
        if os.path.exists(smapTile.FPN): 
            
            self.session._InsertSmapData(smapTile.query)
            
            statusD = {'smapid': smapid,'column':'downloaded', 'status': 'Y'}
            
            self.session._UpdateSmapStatus(statusD)
            
        else:

            if self.pp.process.parameters.asscript:
                
                cmd = 'mkdir -p %(FP)s;\n' %{'FP':smapTile.FP}
                
                self.downloadScriptF.write(cmd)
                
            datedir = mj_dt.DateToStrPointDate(acqdate)
            
            localTempFPN = os.path.join(self.tempFP,smapTile.FN)
            
            self.dlL.append({'query':smapTile.query,'productversion':source,'datedir':datedir,'fn':smapfilename,'dstFPN':smapTile.FPN,'tempFPN':localTempFPN,'smapid':smapid})
  
    def _ConstructDaacTile(self,tile,sdpath):
        '''
        '''
        smapfilename, smapid, source, product, version, content, acqdate = tile
        #construct the composition
        compD = {'source':source, 'product':product, 'version':version, 'content':content, 'system':'smap', 'division':'region'}
        #Invoke the composition
        comp = SmapComposition(compD)
        #Set the datum
        datumD = {'acqdatestr': mj_dt.DateToStrDate(acqdate), 'acqdate':acqdate}
        #Set the filename
        FN = smapfilename
        #Set the locus         
        loc = 'global'
        #Set the locuspath
        locusPath = 'global'
        #Construct the locus dictionary
        locusD = {'locus':loc, 'path':locusPath}
        #Invoke and return a SentinelTile             
        return SmapTile(smapid, comp, locusD, datumD, sdpath, FN)
    
    def _ConstructSmapLayer(self,compD,acqdate,compFormatD):
        '''
        '''
        comp = Composition(compD, self.process.system.dstsystem, self.process.system.dstdivision)
        comp._Update(compFormatD)
        datumD = {'acqdatestr': mj_dt.DateToStrDate(acqdate), 'acqdate':acqdate}
            
        #Set the locus         
        loc = 'global'
            
        #Set the locuspath
        locusPath = 'global'
            
        #Construct the locus dictionary
        locusD = {'locus':loc, 'path':locusPath}
            
        filepath = lambda: None
        filepath.volume = self.pp.dstPath.volume; filepath.hdrfiletype = self.pp.dstPath.hdr
            
        #Create a standard reaster layer
        bandR = RasterLayer(comp, locusD, datumD, filepath)

        return bandR
                    
    def _ReadSMAPhtml(self,session,srcFPN):
        '''
        '''
        
        queryD = self._ParseSmapWgetHTML(srcFPN)
        
        session._InsertSmapData(queryD)
        
    def _ParseSmapWgetHTML(self, FPN):
        '''
        '''
        
        tmpFP = os.path.split(FPN)[0]
        
        tmpFP = os.path.split(tmpFP)[0]
        
        tmpFP = os.path.join(tmpFP,'tmpcsv')
        
        if not os.path.exists(tmpFP):
        
            os.makedirs(tmpFP)

        FPN = 'file://%(fpn)s' %{'fpn':FPN}
        
        req = urllib.request.Request(FPN)
        
        with urllib.request.urlopen(req) as response:
            
            html = response.read()
            
        parser = MjHTMLParser()

        parser.queryD = {}
        
        parser.feed(str(html)) 
        
        return (parser.queryD)

    def _AccessSMAP(self):   
        '''This is similar to _AccessMODIS
        '''
        serverurl = self.pp.process.parameters.serverurl
        for tile in self.dlL:
            remotepath = os.path.join(serverurl,'SMAP',tile['productversion'],tile['datedir'])
            url = os.path.join(remotepath,tile['fn']) 

            home = os.path.expanduser("~")
            cookieFPN = os.path.join(home,'.urs_cookies')
            cmd = "curl -n -L -c %(c)s -b %(c)s  %(r)s --output %(l)s;" %{'u':self.pp.process.parameters.remoteuser, 'c':cookieFPN, 'r':url, 'l':tile['tempFPN']}
            cmd = "%(cmd)s mv %(output)s %(dstFPN)s;" %{'cmd':cmd,'output':tile['tempFPN'], 'dstFPN':tile['dstFPN']}
            if self.pp.process.parameters.asscript:
                cmdL = cmd.split(';')
                for c in cmdL:
                    if len(c) > 1:
                        writeln = '%(c)s;\n' %{'c':c}
                        self.downloadScriptF.write(writeln)
            else:
                os.system(cmd)
                statusD = {'smapid': tile['smapid'],'column':'downloaded', 'status': 'Y'}
                self.session._UpdateSmapStatus(statusD)            

    def _ExtractSmapHdf(self):
        '''Extract the SMAP hdf file
        '''
        #Set asscript to True, this will create a shell file for downloading all missing tiles, if any
        #self.pp.process.parameters.asscript = True
        self.tempFP = os.path.join('/Volumes',self.pp.srcPath.volume, 'smap', 'temp')
        if self.pp.process.parameters.asscript:
            
            shFP = os.path.join(self.tempFP, 'script')
            if not os.path.exists(shFP):
                os.makedirs(shFP)
            shFN = 'explode_%(prod)s.sh' %{'prod':self.pp.process.parameters.product}
            explodeShFPN = os.path.join(shFP,shFN)
            shFN = 'download_%(prod)s.sh' %{'prod':self.pp.process.parameters.product}
            downloadShFPN = os.path.join(shFP,shFN)
            self.explodeScriptF = open(explodeShFPN,'w')
            self.downloadScriptF = open(downloadShFPN,'w')
            
        #Get the tiles
        statusD = {}
        statusD['downloaded'] = 'Y'
        if not self.process.overwrite and self.pp.process.parameters.exploded:
            statusD['exploded'] = 'Y'
            
        tiles = self._IniTileDownload(statusD)

        #Search template for layers to extract
        #Get the layers to extract for this product + version
        self.paramL = ['source', 'product', 'content', 'layerid', 'prefix', 'suffix', 'celltype', 'dataunit', 'cellnull', 'scalefac', 'measure', 'offsetadd', 'region', 'fileext', 'hdffolder', 'hdfgrid']
        queryD = {'source': '%(p)s.%(v)s' %{'p':self.pp.process.parameters.product, 'v':self.pp.process.parameters.version},'retrieve':'Y'}
        self.extractLayerL = self.session._SelectTemplateLayersOnSource(queryD, self.paramL)

        if len(self.extractLayerL) == 0:
            exitstr = 'No layers to exract for smap', queryD
            exit(exitstr)
        missingFlag = False

        for tile in tiles:
            #Construct the smap tile
            smapTile = self._ConstructDaacTile(tile,self.pp.srcPath)
            smapfilename, smapid, source, product, version, content, acqdate = tile
            if not smapTile._Exists():
                warnstr = ('warning the smaptile missing: %s' %(smapTile.FPN))
                print (warnstr)
                self._AddDownload(tile,self.pp.srcPath)
                missingFlag = True
                continue 
            nrExploded = self._ExplodeH5(smapTile, acqdate, product)
            print ('    smap.h5, nrexploded',   smapfilename,nrExploded)

            if nrExploded == len(self.extractLayerL):    
                statusD = {'smapid': smapid,'column':'organized', 'status': 'Y'}
                self.session._UpdateSmapStatus(statusD)
                statusD = {'smapid': smapid,'column':'exploded', 'status': 'Y'}
                self.session._UpdateSmapStatus(statusD)
        #Write the missing tiles to the access shell script
        self._AccessSMAP()
        if self.pp.process.parameters.asscript:
            self.explodeScriptF.close()
            self.downloadScriptF.close()
            printstr = 'To explode tiles you can run the scriptfile %(fpn)s' %{'fpn':explodeShFPN}
            print (printstr)
            if missingFlag:
                printstr = 'To download missing tiles you can run the scriptfile %(fpn)s' %{'fpn':self.downloadShFPN}
                print (printstr)
            
            
           
                   
    def _ExplodeH5(self, smapTile, acqdate, product):
        #  
        nrExploded = 0 
        for extraclayer in self.extractLayerL:
            extractD = dict(zip(self.paramL,extraclayer))

            dstLayer = self._ConstructLayer(extractD,acqdate)

            if dstLayer._Exists():
                if self.process.overwrite:

                    os.remove(dstLayer.FPN)
                else:
                    nrExploded += 1
                    self.session._InsertLayer(dstLayer, self.process.overwrite, self.process.delete)
                    continue
            self._Hdf5_2_geotiff(extractD,product,smapTile,dstLayer)

            if os.path.isfile(dstLayer.FPN):
                nrExploded += 1
                self.session._InsertLayer(dstLayer, self.process.overwrite, self.process.delete)
            
            '''
            #The giving of fixed coordinates is not good
            -17349514.3530680164694786,-7296524.6913595553487539 : 17349514.3530680164694786,7296524.6913595534861088
            cmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdal_translate '

            cmd = '%(cmd)s -a_ullr -17349514.353 7296524.691 17349514.353 -7296524.691 ' %{'cmd':cmd}


            #SET proj to EASE GRID 2 (epsg:6033)
            cmd = '%(cmd)s -a_srs "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" ' %{'cmd':cmd}
            
            cmd = '%(cmd)s  -a_nodata -9999  ' %{'cmd':cmd}
            cmd = '%(cmd)s HDF5:"%(hdf)s"://%(folder)s/%(grid)s %(dst)s' %{'cmd':cmd,
                    'hdf':smapTile.FPN,'folder':extractD['hdffolder'],'grid':extractD['hdfgrid'], 
                    'dst':dstLayer.FPN}
            print (cmd)

            if self.pp.process.parameters.asscript:
                cmd = '%(cmd)s;\n' %{'cmd':cmd}
                self.explodeScriptF.write(cmd)
            else:   
                os.system(cmd)
                #register band
                if os.path.isfile(dstLayer.FPN):
                    nrExploded += 1
                    self.session._InsertLayer(dstLayer, self.process.overwrite, self.process.delete)
            '''
        return nrExploded
    
    def _Hdf5_2_geotiff(self,extractD,product,smapTile,dstLayer):
        #Some products have more than one hdffolder and thus more than one lon/lat pair 
        if product == 'SPL3SMP' and extractD['hdffolder'] == 'Soil_Moisture_Retrieval_Data_PM':
            queryD = {'hdfgrid':'longitude_pm', 'product':extractD['product'], 'source':extractD['source'], 'hdffolder':extractD['hdffolder']} 
        elif product == 'SPL3SMP-E' and extractD['hdffolder'] == 'Soil_Moisture_Retrieval_Data_PM':
            queryD = {'hdfgrid':'longitude_pm', 'product':extractD['product'], 'source':extractD['source'], 'hdffolder':extractD['hdffolder']} 

        else:
            queryD = {'hdfgrid':'longitude', 'product':extractD['product'], 'source':extractD['source'], 'hdffolder':extractD['hdffolder']} 
        lonLayer = self.session._SelectTemplateLayersOnGrid(queryD, self.paramL)

        if lonLayer == None:
            exitstr = 'No lon/lat data found for SMAP extraction',queryD
            exit(exitstr)

        lonD = dict(zip(self.paramL,lonLayer))
        extractD['longrid'] = lonD['hdfgrid']
        if product == 'SPL3SMP' and extractD['hdffolder'] == 'Soil_Moisture_Retrieval_Data_PM':
            queryD = {'hdfgrid':'latitude_pm', 'product':extractD['product'], 'source':extractD['source'], 'hdffolder':extractD['hdffolder']} 
        elif product == 'SPL3SMP-E' and extractD['hdffolder'] == 'Soil_Moisture_Retrieval_Data_PM':
            queryD = {'hdfgrid':'latitude_pm', 'product':extractD['product'], 'source':extractD['source'], 'hdffolder':extractD['hdffolder']} 

        else:
            queryD = {'hdfgrid':'latitude', 'product':extractD['product'], 'source':extractD['source'], 'hdffolder':extractD['hdffolder']} 
        latLayer = self.session._SelectTemplateLayersOnGrid(queryD, self.paramL)
        latD = dict(zip(self.paramL,latLayer))

        extractD['latgrid'] = latD['hdfgrid']
        extractD['lonlatfolder'] = latD['hdffolder']
        print ('retrieving, ',smapTile.FPN,dstLayer.FPN)
        hdf5_2_geotiff.Retrieve(smapTile.FPN, extractD, dstLayer.FPN)
  
    def _ConstructLayer(self,extractD,acqdate):
        '''
        '''
        compD = extractD
        comp = Composition(compD, 'smap', 'region')
        datumD = {'acqdatestr': mj_dt.DateToStrDate(acqdate), 'acqdate':acqdate}
            
        #Set the locus         
        loc = extractD['region']
        
        #Set the locuspath
        locusPath = extractD['region']
        
        #Construct the locus dictionary
        locusD = {'locus':loc,  'path':locusPath}
        
        filepath = lambda: None
        filepath.volume = self.pp.dstPath.volume; filepath.hdrfiletype = extractD['fileext']
        
        #Create a standard raster layer
        return RasterLayer(comp, locusD, datumD, filepath)
      
class MjHTMLParser(HTMLParser):
            
    def handle_starttag(self, tag, attrs):
        # Only parse the 'anchor' tag.
        if tag == "a":
            # Check the list of defined attributes.
            for name, value in attrs:
                # If href is defined, print it.
                if name == "href" and 'SMAP' in value:
                    if value[0:6] == '/SMAP/':
                        source = value.split('/')[2]
                        product,version = source.split('.')
                        self.queryD['source'] = source.replace('_','-')
                        self.queryD['product'] = product.replace('_','-')
                        self.queryD['version'] = version
                    elif value[0:4] == 'SMAP' and os.path.splitext(value)[1] == '.h5':
                        smapfilename = value
                        fnParts = value.split('_')
                        if len(fnParts) == 8 and '_E_' in smapfilename:
                            sensor, level, type, code, enhanced,  acqdatestr, Rcode, vext = value.split('_')
                        elif len(fnParts) == 7:
                            sensor, level, type, code,acqdatestr, Rcode, vext = value.split('_')
                        else:
                            errorstringnotthere
                        acqdate = mj_dt.yyyymmddDate(acqdatestr)
                        self.queryD['smapid'] = os.path.splitext(smapfilename)[0]
                        self.queryD['smapfilename'] = smapfilename
                        self.queryD['acqdate'] = acqdate
                        self.queryD['doy'] = mj_dt.DateToDOY(acqdate)
                    elif value[0:4] == 'SMAP' and os.path.splitext(value)[1] == '.xml':
                        metafilename = value
                    elif value[0:4] == 'SMAP' and os.path.splitext(value)[1] == '.qa':
                        qafilename = value
