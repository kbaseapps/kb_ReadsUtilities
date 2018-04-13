# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import sys
import shutil
import hashlib
import subprocess
import requests
import re
import traceback
import uuid
from datetime import datetime
from pprint import pprint, pformat
import numpy as np
import gzip
import random

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
#from biokbase.workspace.client import Workspace as workspaceService                                                              
from Workspace.WorkspaceClient import Workspace as workspaceService
from requests_toolbelt import MultipartEncoder  # added                                                                           
from biokbase.AbstractHandle.Client import AbstractHandle as HandleService

# SDK Utils                                                                                                                       
from ReadsUtils.ReadsUtilsClient import ReadsUtils
from SetAPI.SetAPIServiceClient import SetAPI
from KBaseReport.KBaseReportClient import KBaseReport

# silence whining                                                                                                                 
import requests
requests.packages.urllib3.disable_warnings()

#END_HEADER


class kb_ReadsUtilities:
    '''
    Module Name:
    kb_ReadsUtilities

    Module Description:
    ** A KBase module: kb_ReadsUtilities
**
** This module contains basic utility Apps for manipulating Reads Libraries
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "1.0.2"
    GIT_URL = "https://github.com/kbaseapps/kb_ReadsUtilities"
    GIT_COMMIT_HASH = "7caaaa1fd572e5d774a5aaa194150cb9f3fef2f5"

    #BEGIN_CLASS_HEADER
    workspaceURL = None
    shockURL = None
    handleURL = None
    serviceWizardsURL = None
    callbackURL = None
    scratch = None


    # target is a list for collecting log messages                                                                                
    def log(self, target, message):
        # we should do something better here...                                                                                   
        if target is not None:
            target.append(message)
        print(message)
        sys.stdout.flush()

    def get_single_end_read_library(self, ws_data, ws_info, forward):
        pass

    def get_feature_set_seqs(self, ws_data, ws_info):
        pass

    def get_genome_feature_seqs(self, ws_data, ws_info):
        pass

    def get_genome_set_feature_seqs(self, ws_data, ws_info):
        pass


    # Helper script borrowed from the transform service, logger removed                                                           
    #                                                                                                                             
    def upload_file_to_shock(self,
                             console,  # DEBUG                                                                                    
                             shock_service_url = None,
                             filePath = None,
                             ssl_verify = True,
                             token = None):
        """                                                                                                                       
        Use HTTP multi-part POST to save a file to a SHOCK instance.                                                              
        """
        self.log(console,"UPLOADING FILE "+filePath+" TO SHOCK")

        if token is None:
            raise Exception("Authentication token required!")

        #build the header                                                                                                         
        header = dict()
        header["Authorization"] = "Oauth {0}".format(token)
        if filePath is None:
            raise Exception("No file given for upload to SHOCK!")

        dataFile = open(os.path.abspath(filePath), 'rb')
        m = MultipartEncoder(fields={'upload': (os.path.split(filePath)[-1], dataFile)})
        header['Content-Type'] = m.content_type

        #logger.info("Sending {0} to {1}".format(filePath,shock_service_url))                                                     
        try:
            response = requests.post(shock_service_url + "/node", headers=header, data=m, allow_redirects=True, verify=ssl_verify)
            dataFile.close()
        except:
            dataFile.close()
            raise
        if not response.ok:
            response.raise_for_status()
        result = response.json()
        if result['error']:
            raise Exception(result['error'][0])
        else:
            return result["data"]


    def upload_SingleEndLibrary_to_shock_and_ws (self,
                                                 ctx,
                                                 console,  # DEBUG                                                                
                                                 workspace_name,
                                                 obj_name,
                                                 file_path,
                                                 provenance,
                                                 sequencing_tech):

        self.log(console,'UPLOADING FILE '+file_path+' TO '+workspace_name+'/'+obj_name)

        # 1) upload files to shock                                                                                                
        token = ctx['token']
        forward_shock_file = self.upload_file_to_shock(
            console,  # DEBUG                                                                                                     
            shock_service_url = self.shockURL,
            filePath = file_path,
            token = token
            )
        #pprint(forward_shock_file)                                                                                               
        self.log(console,'SHOCK UPLOAD DONE')

        # 2) create handle                                                                                                        
        self.log(console,'GETTING HANDLE')
        hs = HandleService(url=self.handleURL, token=token)
        forward_handle = hs.persist_handle({
                                        'id' : forward_shock_file['id'],
                                        'type' : 'shock',
                                        'url' : self.shockURL,
                                        'file_name': forward_shock_file['file']['name'],
                                        'remote_md5': forward_shock_file['file']['checksum']['md5']})

        # 3) save to WS                                                                                                           
        self.log(console,'SAVING TO WORKSPACE')
        single_end_library = {
            'lib': {
                'file': {
                    'hid':forward_handle,
                    'file_name': forward_shock_file['file']['name'],
                    'id': forward_shock_file['id'],
                    'url': self.shockURL,
                    'type':'shock',
                    'remote_md5':forward_shock_file['file']['checksum']['md5']
                },
                'encoding':'UTF8',
                'type':'fasta',
                'size':forward_shock_file['file']['size']
            },
            'sequencing_tech':sequencing_tech
        }
        self.log(console,'GETTING WORKSPACE SERVICE OBJECT')
        ws = workspaceService(self.workspaceURL, token=ctx['token'])
        self.log(console,'SAVE OPERATION...')
        new_obj_info = ws.save_objects({
                        'workspace':workspace_name,
                        'objects':[
                            {
                                'type':'KBaseFile.SingleEndLibrary',
                                'data':single_end_library,
                                'name':obj_name,
                                'meta':{},
                                'provenance':provenance
                            }]
                        })
        self.log(console,'SAVED TO WORKSPACE')

        return new_obj_info[0]

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.handleURL = config['handle-service-url']
        self.serviceWizardURL = config['service-wizard-url']

        #self.callbackURL = os.environ['SDK_CALLBACK_URL'] if os.environ['SDK_CALLBACK_URL'] != None else 'https://kbase.us/services/njs_wrapper'  # DEBUG                                                                                                         
        self.callbackURL = os.environ.get('SDK_CALLBACK_URL')
#        if self.callbackURL == None:                                                                                             
#            self.callbackURL = os.environ['SDK_CALLBACK_URL']                                                                    
        if self.callbackURL == None:
            raise ValueError ("SDK_CALLBACK_URL not set in environment")

        self.scratch = os.path.abspath(config['scratch'])
        # HACK!! temporary hack for issue where megahit fails on mac because of silent named pipe error                           
        #self.host_scratch = self.scratch                                                                                         
        #self.scratch = os.path.join('/kb','module','local_scratch')                                                              
        # end hack                                                                                                                
        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)
        #END_CONSTRUCTOR
        pass


    def KButil_FASTQ_to_FASTA(self, ctx, params):
        """
        :param params: instance of type "KButil_FASTQ_to_FASTA_Params"
           (KButil_FASTQ_to_FASTA() ** ** Method for Converting a FASTQ
           SingleEndLibrary to a FASTA SingleEndLibrary) -> structure:
           parameter "workspace_name" of type "workspace_name" (** The
           workspace object refs are of form: ** **    objects =
           ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_ref" of type "data_obj_ref", parameter
           "output_name" of type "data_obj_name"
        :returns: instance of type "KButil_FASTQ_to_FASTA_Output" ->
           structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN KButil_FASTQ_to_FASTA
        console = []
        invalid_msgs = []
        self.log(console,'Running KButil_FASTQ_to_FASTA with params=')
        self.log(console, "\n"+pformat(params))
        report = ''
#        report = 'Running KButil_FASTQ_to_FASTA with params='
#        report += "\n"+pformat(params)


        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
        if 'input_ref' not in params:
            raise ValueError('input_ref parameter is required')
        if 'output_name' not in params:
            raise ValueError('output_name parameter is required')
        input_ref = params['input_ref']


        # Download Reads
        #
        sequencing_tech = 'N/A'
        self.log (console, "DOWNLOADING READS")  # DEBUG
        try:
            readsUtils_Client = ReadsUtils (url=self.callbackURL, token=ctx['token'])  # SDK local
        except Exception as e:
            raise ValueError('Unable to get ReadsUtils Client' +"\n" + str(e))
        try:
            readsLibrary = readsUtils_Client.download_reads ({'read_libraries': [input_ref],
                                                             'interleaved': 'false'
                                                             })
        except Exception as e:
            raise ValueError('Unable to download read library sequences from workspace: (' + str(input_ref) +")\n" + str(e))

        forward_reads_file_path = readsLibrary['files'][input_ref]['files']['fwd']
        sequencing_tech     = readsLibrary['files'][input_ref]['sequencing_tech']


        #### Create the file to upload
        ##
        output_file_name   = params['output_name']+'.fna'
        output_file_path  = os.path.join(self.scratch,output_file_name)
        input_file_handle  = open(forward_reads_file_path, 'r', -1)
        output_file_handle = open(output_file_path, 'w', -1)
        self.log(console, 'PROCESSING reads file: '+str(forward_reads_file_path))

        seq_cnt = 0
        header = None
        last_header = None
        last_seq_buf = None
        rec_line_i = -1
        for line in input_file_handle:
            rec_line_i += 1
            if rec_line_i == 3:
                rec_line_i = -1
            elif rec_line_i == 0:
                if not line.startswith('@'):
                    raise ValueError ("badly formatted rec line: '"+line+"'")
                seq_cnt += 1
                header = line[1:]
                if last_header != None:
                    output_file_handle.write('>'+last_header)
                    output_file_handle.write(last_seq_buf)
                last_seq_buf = None
                last_header = header
            elif rec_line_i == 1:
                last_seq_buf = line
        if last_header != None:
            output_file_handle.write('>'+last_header)
            output_file_handle.write(last_seq_buf)

        input_file_handle.close()
        output_file_handle.close()
        

        # load the method provenance from the context object
        #
        self.log(console,"SETTING PROVENANCE")  # DEBUG
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
        provenance[0]['input_ws_objects'].append(params['input_ref'])
        provenance[0]['service'] = 'kb_ReadsUtilities'
        provenance[0]['method'] = 'KButil_FASTQ_to_FASTA'


        # Upload results
        #
        if len(invalid_msgs) == 0:
            self.log(console,"UPLOADING RESULTS")  # DEBUG
            self.upload_SingleEndLibrary_to_shock_and_ws (ctx,
                                                      console,  # DEBUG
                                                      params['workspace_name'],
                                                      params['output_name'],
                                                      output_file_path,
                                                      provenance,
                                                      sequencing_tech
                                                      )

        # build output report object
        #
        self.log(console,"BUILDING REPORT")  # DEBUG
        if len(invalid_msgs) == 0:
            report += 'sequences in library:  '+str(seq_cnt)+"\n"
            reportObj = {
                'objects_created':[{'ref':params['workspace_name']+'/'+params['output_name'], 'description':'KButil_FASTQ_to_FASTA'}],
                'text_message':report
                }
        else:
            report += "FAILURE:\n\n"+"\n".join(invalid_msgs)+"\n"
            reportObj = {
                'objects_created':[],
                'text_message':report
                }

        reportName = 'kbutil_fastq_to_fasta_report_'+str(uuid.uuid4())
        ws = workspaceService(self.workspaceURL, token=ctx['token'])
        report_obj_info = ws.save_objects({
#                'id':info[6],
                'workspace':params['workspace_name'],
                'objects':[
                    {
                        'type':'KBaseReport.Report',
                        'data':reportObj,
                        'name':reportName,
                        'meta':{},
                        'hidden':1,
                        'provenance':provenance
                    }
                ]
            })[0]

        self.log(console,"BUILDING RETURN OBJECT")
#        returnVal = { 'output_report_name': reportName,
#                      'output_report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
#                      'output_filtered_ref': params['workspace_name']+'/'+params['output_filtered_name']
#                      }
        returnVal = { 'report_name': reportName,
                      'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
                      }
        self.log(console,"KButil_FASTQ_to_FASTA DONE")
        #END KButil_FASTQ_to_FASTA

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method KButil_FASTQ_to_FASTA return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def KButil_Split_Reads(self, ctx, params):
        """
        :param params: instance of type "KButil_Split_Reads_Params"
           (KButil_Split_Reads() ** **  Method for spliting a ReadsLibrary
           into evenly sized ReadsLibraries) -> structure: parameter
           "workspace_name" of type "workspace_name" (** The workspace object
           refs are of form: ** **    objects = ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_ref" of type "data_obj_ref", parameter
           "output_name" of type "data_obj_name", parameter "split_num" of
           Long, parameter "desc" of String
        :returns: instance of type "KButil_Split_Reads_Output" -> structure:
           parameter "report_name" of type "data_obj_name", parameter
           "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN KButil_Split_Reads
        console = []
        report = ''
        self.log(console, 'Running KButil_Split_Reads() with parameters: ')
        self.log(console, "\n"+pformat(params))
        
        token = ctx['token']
        wsClient = workspaceService(self.workspaceURL, token=token)
        headers = {'Authorization': 'OAuth '+token}
        env = os.environ.copy()
        env['KB_AUTH_TOKEN'] = token
        
        #SERVICE_VER = 'dev'  # DEBUG
        SERVICE_VER = 'release'

        # param checks
        required_params = ['workspace_name',
                           'input_ref', 
                           'output_name'
                           ]
        for required_param in required_params:
            if required_param not in params or params[required_param] == None:
                raise ValueError ("Must define required param: '"+required_param+"'")
            
        # and param defaults
        defaults = { 'split_num': 10
                   }
        for arg in defaults.keys():
            if arg not in params or params[arg] == None or params[arg] == '':
                params[arg] = defaults[arg]

        # load provenance
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        provenance[0]['input_ws_objects']=[str(params['input_ref'])]


        # Determine whether read library is of correct type
        #
        try:
            # object_info tuple
            [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)
            
            input_reads_ref = params['input_ref']

            input_reads_obj_info = wsClient.get_object_info_new ({'objects':[{'ref':input_reads_ref}]})[0]
            input_reads_obj_type = input_reads_obj_info[TYPE_I]
            #input_reads_obj_version = input_reads_obj_info[VERSION_I]  # this is object version, not type version

        except Exception as e:
            raise ValueError('Unable to get read library object info from workspace: (' + str(input_reads_ref) +')' + str(e))

        input_reads_obj_type = re.sub ('-[0-9]+\.[0-9]+$', "", input_reads_obj_type)  # remove trailing version

        acceptable_types = ["KBaseFile.PairedEndLibrary", "KBaseFile.SingleEndLibrary"]
        if input_reads_obj_type not in acceptable_types:
            raise ValueError ("Input reads of type: '"+input_reads_obj_type+"'.  Must be one of "+", ".join(acceptable_types))


        # Download Reads
        #
        self.log (console, "DOWNLOADING READS")  # DEBUG
        try:
            readsUtils_Client = ReadsUtils (url=self.callbackURL, token=ctx['token'])  # SDK local
        except Exception as e:
            raise ValueError('Unable to get ReadsUtils Client' +"\n" + str(e))
        try:
            readsLibrary = readsUtils_Client.download_reads ({'read_libraries': [input_reads_ref],
                                                             'interleaved': 'false'
                                                             })
        except Exception as e:
            raise ValueError('Unable to download read library sequences from workspace: (' + str(input_reads_ref) +")\n" + str(e))


        # Paired End
        #
        if input_reads_obj_type == "KBaseFile.PairedEndLibrary":

            # Download reads Libs to FASTQ files
            input_fwd_file_path = readsLibrary['files'][input_reads_ref]['files']['fwd']
            input_rev_file_path = readsLibrary['files'][input_reads_ref]['files']['rev']
            sequencing_tech     = readsLibrary['files'][input_reads_ref]['sequencing_tech']

            # file paths
            input_fwd_path = re.sub ("\.fastq$", "", input_fwd_file_path)
            input_fwd_path = re.sub ("\.FASTQ$", "", input_fwd_path)
            input_rev_path = re.sub ("\.fastq$", "", input_rev_file_path)
            input_rev_path = re.sub ("\.FASTQ$", "", input_rev_path)
            output_fwd_paired_file_path_base   = input_fwd_path+"_fwd_paired"
            output_rev_paired_file_path_base   = input_rev_path+"_rev_paired"
            output_fwd_unpaired_file_path = input_fwd_path+"_fwd_unpaired.fastq"
            output_rev_unpaired_file_path = input_rev_path+"_rev_unpaired.fastq"
            # set up for file io
            total_paired_reads = 0
            total_paired_reads_by_set = []
            total_unpaired_fwd_reads = 0
            total_unpaired_rev_reads = 0
            fwd_ids = dict()
            paired_lib_i = dict()
            unpaired_buf_size = 0
            paired_buf_size = 100000
            recs_beep_n = 1000000

            # read fwd file to get fwd ids
#            rec_cnt = 0  # DEBUG
            self.log (console, "GETTING IDS")  # DEBUG
            with open (input_fwd_file_path, 'r', 0) as input_reads_file_handle:
                rec_line_i = -1
                for line in input_reads_file_handle:
                    rec_line_i += 1
                    if rec_line_i == 3:
                        rec_line_i = -1
                    elif rec_line_i == 0:
                        if not line.startswith('@'):
                            raise ValueError ("badly formatted rec line: '"+line+"'")
                        read_id = line.rstrip('\n')
                        read_id = re.sub ("[ \t]+.*$", "", read_id)
                        read_id = re.sub ("[\/\.\_\-\:\;][012lrLRfrFR53]\'*$", "", read_id)
                        fwd_ids[read_id] = True
                        # DEBUG
#                        if rec_cnt % 100 == 0:
#                            self.log(console,"read_id: '"+str(read_id)+"'")
#                        rec_cnt += 1 

            # determine paired and unpaired rev, split paired rev
            #   write unpaired rev, and store lib_i for paired
            self.log (console, "WRITING REV SPLIT PAIRED")  # DEBUG
            paired_output_reads_file_handles = []
            for lib_i in range(params['split_num']):
                paired_output_reads_file_handles.append(open (output_rev_paired_file_path_base+"-"+str(lib_i)+".fastq", 'w', paired_buf_size))
                total_paired_reads_by_set.append(0)

            rec_buf = []
            unpaired_rev_buf = []
            last_read_id = None
            paired_cnt = 0
            capture_type_paired = False

            with open (input_rev_file_path, 'r', 0) as input_reads_file_handle:
                rec_line_i = -1
                for line in input_reads_file_handle:
                    rec_line_i += 1
                    if rec_line_i == 3:
                        rec_line_i = -1
                    elif rec_line_i == 0:
                        if not line.startswith('@'):
                            raise ValueError ("badly formatted rec line: '"+line+"'")
                        if last_read_id != None:
                            if capture_type_paired:
                                lib_i = paired_cnt % params['split_num']
                                total_paired_reads_by_set[lib_i] += 1
                                paired_lib_i[last_read_id] = lib_i
                                paired_output_reads_file_handles[lib_i].writelines(rec_buf)
                                paired_cnt += 1
                                if paired_cnt != 0 and paired_cnt % recs_beep_n == 0:
                                    self.log(console,"\t"+str(paired_cnt)+" recs processed")
                            else:
                                unpaired_rev_buf.extend(rec_buf)
                            rec_buf = []
                        read_id = line.rstrip('\n')
                        read_id = re.sub ("[ \t]+.*$", "", read_id)
                        read_id = re.sub ("[\/\.\_\-\:\;][012lrLRfrFR53]\'*$", "", read_id)
                        last_read_id = read_id
                        try:
#                            self.log(console,"CHECKING: '"+str(read_id)+"'") # DEBUG
                            found = fwd_ids[read_id]
#                            self.log(console,"FOUND PAIR: '"+str(read_id)+"'") # DEBUG
                            total_paired_reads += 1
                            capture_type_paired = True
                        except:
                            total_unpaired_rev_reads += 1
                            capture_type_paired = False
                    rec_buf.append(line)
                # last record
                if len(rec_buf) > 0:
                    if capture_type_paired:
                        lib_i = paired_cnt % params['split_num']
                        total_paired_reads_by_set[lib_i] += 1
                        paired_lib_i[last_read_id] = lib_i
                        paired_output_reads_file_handles[lib_i].writelines(rec_buf)
                        paired_cnt += 1
                        if paired_cnt != 0 and paired_cnt % recs_beep_n == 0:
                            self.log(console,"\t"+str(paired_cnt)+" recs processed")
                    else:
                        unpaired_rev_buf.extend(rec_buf)
                    rec_buf = []

            for output_handle in paired_output_reads_file_handles:
                output_handle.close()

            self.log(console,"\t"+str(paired_cnt)+" recs processed")
            self.log (console, "WRITING REV UNPAIRED")  # DEBUG
            output_reads_file_handle = open (output_rev_unpaired_file_path, 'w', 0)
            output_reads_file_handle.writelines(unpaired_rev_buf)
            output_reads_file_handle.close()


            # split fwd paired and write unpaired fwd
            self.log (console, "WRITING FWD SPLIT PAIRED")  # DEBUG
            paired_output_reads_file_handles = []
            for lib_i in range(params['split_num']):
                paired_output_reads_file_handles.append(open (output_fwd_paired_file_path_base+"-"+str(lib_i)+".fastq", 'w', paired_buf_size))

            rec_buf = []
            unpaired_fwd_buf = []
            last_read_id = None
            paired_cnt = 0
            capture_type_paired = False

            with open (input_fwd_file_path, 'r', 0) as input_reads_file_handle:
                rec_line_i = -1
                for line in input_reads_file_handle:
                    rec_line_i += 1
                    if rec_line_i == 3:
                        rec_line_i = -1
                    elif rec_line_i == 0:
                        if not line.startswith('@'):
                            raise ValueError ("badly formatted rec line: '"+line+"'")
                        if last_read_id != None:
                            if capture_type_paired:
                                lib_i = paired_lib_i[last_read_id]
                                paired_output_reads_file_handles[lib_i].writelines(rec_buf)
                                paired_cnt += 1
                                if paired_cnt != 0 and paired_cnt % recs_beep_n == 0:
                                    self.log(console,"\t"+str(paired_cnt)+" recs processed")
                            else:
                                unpaired_fwd_buf.extend(rec_buf)
                            rec_buf = []
                        read_id = line.rstrip('\n')
                        read_id = re.sub ("[ \t]+.*$", "", read_id)
                        read_id = re.sub ("[\/\.\_\-\:\;][012lrLRfrFR53]\'*$", "", read_id)
                        last_read_id = read_id
                        try:
                            found = paired_lib_i[read_id]
                            capture_type_paired = True
                        except:
                            total_unpaired_fwd_reads += 1
                            capture_type_paired = False
                    rec_buf.append(line)
                # last rec
                if len(rec_buf) > 0:
                    if capture_type_paired:
                        lib_i = paired_lib_i[last_read_id]
                        paired_output_reads_file_handles[lib_i].writelines(rec_buf)
                        paired_cnt += 1
                        if paired_cnt != 0 and paired_cnt % recs_beep_n == 0:
                            self.log(console,"\t"+str(paired_cnt)+" recs processed")
                    else:
                        unpaired_fwd_buf.extend(rec_buf)
                    rec_buf = []

            for output_handle in paired_output_reads_file_handles:
                output_handle.close()

            self.log(console,"\t"+str(paired_cnt)+" recs processed")
            self.log (console, "WRITING FWD UNPAIRED")  # DEBUG
            output_reads_file_handle = open (output_fwd_unpaired_file_path, 'w', 0)
            output_reads_file_handle.writelines(unpaired_fwd_buf)
            output_reads_file_handle.close()


            # store report
            #
            report += "TOTAL PAIRED READS: "+str(total_paired_reads)+"\n"
            report += "TOTAL UNPAIRED FWD READS: "+str(total_unpaired_fwd_reads)+"\n"
            report += "TOTAL UNPAIRED REV READS: "+str(total_unpaired_rev_reads)+"\n"
            report += "\n"
            for lib_i in range(params['split_num']):
                report += "PAIRED READS IN SET "+str(lib_i)+": "+str(total_paired_reads_by_set[lib_i])+"\n"


            # upload paired reads
            #
            self.log (console, "UPLOAD PAIRED READS LIBS")  # DEBUG
            paired_obj_refs = []
            for lib_i in range(params['split_num']):
                output_fwd_paired_file_path = output_fwd_paired_file_path_base+"-"+str(lib_i)+".fastq"
                output_rev_paired_file_path = output_rev_paired_file_path_base+"-"+str(lib_i)+".fastq"
                if not os.path.isfile (output_fwd_paired_file_path) \
                     or os.path.getsize (output_fwd_paired_file_path) == 0 \
                   or not os.path.isfile (output_rev_paired_file_path) \
                     or os.path.getsize (output_rev_paired_file_path) == 0:
                    
                    raise ValueError ("failed to create paired output")
                else:
                    output_obj_name = params['output_name']+'_paired-'+str(lib_i)
                    self.log(console, 'Uploading paired reads: '+output_obj_name)
                    paired_obj_refs.append (readsUtils_Client.upload_reads ({ 'wsname': str(params['workspace_name']),
                                                                              'name': output_obj_name,
                                                                              # remove sequencing_tech when source_reads_ref working
                                                                              #'sequencing_tech': sequencing_tech,
                                                                              'source_reads_ref': input_reads_ref,
                                                                              'fwd_file': output_fwd_paired_file_path,
                                                                              'rev_file': output_rev_paired_file_path
                                                                              })['obj_ref'])
                    

            # upload reads forward unpaired
            self.log (console, "UPLOAD UNPAIRED FWD READS LIB")  # DEBUG
            unpaired_fwd_ref = None
            if os.path.isfile (output_fwd_unpaired_file_path) \
                and os.path.getsize (output_fwd_unpaired_file_path) != 0:

                output_obj_name = params['output_name']+'_unpaired-fwd'
                self.log(console, '\nUploading trimmed unpaired forward reads: '+output_obj_name)
                unpaired_fwd_ref = readsUtils_Client.upload_reads ({ 'wsname': str(params['workspace_name']),
                                                                     'name': output_obj_name,
                                                                     # remove sequencing_tech when source_reads_ref working
                                                                     #'sequencing_tech': sequencing_tech,
                                                                     'source_reads_ref': input_reads_ref,
                                                                     'fwd_file': output_fwd_unpaired_file_path
                                                                     })['obj_ref']
                

            # upload reads reverse unpaired
            self.log (console, "UPLOAD UNPAIRED REV READS LIB")  # DEBUG
            unpaired_rev_ref = None
            if os.path.isfile (output_rev_unpaired_file_path) \
                and os.path.getsize (output_rev_unpaired_file_path) != 0:

                output_obj_name = params['output_name']+'_unpaired-rev'
                self.log(console, '\nUploading trimmed unpaired reverse reads: '+output_obj_name)
                unpaired_rev_ref = readsUtils_Client.upload_reads ({ 'wsname': str(params['workspace_name']),
                                                                     'name': output_obj_name,
                                                                     # remove sequencing_tech when source_reads_ref working
                                                                     #'sequencing_tech': sequencing_tech,
                                                                     'source_reads_ref': input_reads_ref,
                                                                     'fwd_file': output_rev_unpaired_file_path
                                                                     })['obj_ref']
                

        # SingleEndLibrary
        #
        elif input_reads_obj_type == "KBaseFile.SingleEndLibrary":
            self.log(console, "Downloading Single End reads file...")

            # Download reads Libs to FASTQ files
            input_fwd_file_path = readsLibrary['files'][input_reads_ref]['files']['fwd']
            sequencing_tech     = readsLibrary['files'][input_reads_ref]['sequencing_tech']

            # file paths
            input_fwd_path = re.sub ("\.fastq$", "", input_fwd_file_path)
            input_fwd_path = re.sub ("\.FASTQ$", "", input_fwd_path)
            output_fwd_paired_file_path_base   = input_fwd_path+"_fwd_paired"

            # set up for file io
            total_paired_reads = 0
            total_paired_reads_by_set = []
            paired_buf_size = 1000000

            # split reads
            self.log (console, "WRITING SPLIT SINGLE END READS")  # DEBUG
            paired_output_reads_file_handles = []
            for lib_i in range(params['split_num']):
                paired_output_reads_file_handles.append(open (output_fwd_paired_file_path_base+"-"+str(lib_i)+".fastq", 'w', paired_buf_size))
                total_paired_reads_by_set.append(0)

            rec_buf = []
            last_read_id = None
            paired_cnt = 0
            recs_beep_n = 1000000
            with open (input_fwd_file_path, 'r', 0) as input_reads_file_handle:
                rec_line_i = -1
                for line in input_reads_file_handle:
                    rec_line_i += 1
                    if rec_line_i == 3:
                        rec_line_i = -1
                    elif rec_line_i == 0:
                        if not line.startswith('@'):
                            raise ValueError ("badly formatted rec line: '"+line+"'")
                        if last_read_id != None:
                            total_paired_reads += 1
                            lib_i = paired_cnt % params['split_num']
                            total_paired_reads_by_set[lib_i] += 1
                            paired_output_reads_file_handles[lib_i].writelines(rec_buf)
                            paired_cnt += 1
                            if paired_cnt != 0 and paired_cnt % recs_beep_n == 0:
                                self.log(console,"\t"+str(paired_cnt)+" recs processed")
                            rec_buf = []
                        read_id = line.rstrip('\n')
                        read_id = re.sub ("[ \t]+.*$", "", read_id)
                        read_id = re.sub ("[\/\.\_\-\:\;][012lrLRfrFR53]\'*$", "", read_id)
                        last_read_id = read_id
                    rec_buf.append(line)
                # last rec
                if len(rec_buf) > 0:
                    total_paired_reads += 1
                    lib_i = paired_cnt % params['split_num']
                    total_paired_reads_by_set[lib_i] += 1
                    paired_output_reads_file_handles[lib_i].writelines(rec_buf)
                    paired_cnt += 1
                    if paired_cnt != 0 and paired_cnt % recs_beep_n == 0:
                        self.log(console,"\t"+str(paired_cnt)+" recs processed")
                    rec_buf = []

            for output_handle in paired_output_reads_file_handles:
                output_handle.close()


            # store report
            #
            report += "TOTAL READS: "+str(total_paired_reads)+"\n"
            for lib_i in range(params['split_num']):
                report += "SINGLE END READS IN SET "+str(lib_i)+": "+str(total_paired_reads_by_set[lib_i])+"\n"


            # upload reads
            #
            self.log (console, "UPLOADING SPLIT SINGLE END READS")  # DEBUG
            unpaired_fwd_ref = None
            unpaired_rev_ref = None
            paired_obj_refs = []
            for lib_i in range(params['split_num']):
                output_fwd_paired_file_path = output_fwd_paired_file_path_base+"-"+str(lib_i)+".fastq"
                if not os.path.isfile (output_fwd_paired_file_path) \
                        or os.path.getsize (output_fwd_paired_file_path) == 0:
                    
                    raise ValueError ("failed to create paired output")
                else:
                    output_obj_name = params['output_name']+'-'+str(lib_i)
                    self.log(console, 'Uploading paired reads: '+output_obj_name)
                    paired_obj_refs.append( readsUtils_Client.upload_reads ({ 'wsname': str(params['workspace_name']),
                                                                              'name': output_obj_name,
                                                                              # remove sequencing_tech when source_reads_ref working
                                                                              #'sequencing_tech': sequencing_tech,
                                                                              'source_reads_ref': input_reads_ref,
                                                                              'fwd_file': output_fwd_paired_file_path
                                                                              })['obj_ref'])
                                            
        else:
            raise ValueError ("unknown ReadLibrary type as input: "+str(input_reads_obj_type))


        # save output readsSet
        #
        self.log (console, "SAVING READSSET")  # DEBUG
        items = []
        for lib_i,lib_ref in enumerate(paired_obj_refs):
            label = params['output_name']+'-'+str(lib_i)
            items.append({'ref': lib_ref,
                          'label': label
                          #'data_attachment': ,
                          #'info':
                              })
        description = params['desc']
        output_readsSet_obj = { 'description': params['desc'],
                                'items': items
                                }
        output_readsSet_name = str(params['output_name'])
        setAPI_Client = SetAPI (url=self.serviceWizardURL, token=ctx['token'])  # for dynamic service
        readsSet_ref = setAPI_Client.save_reads_set_v1 ({'workspace_name': params['workspace_name'],
                                                         'output_object_name': output_readsSet_name,
                                                         'data': output_readsSet_obj
                                                         })['set_ref']
                              

        # build report
        #
        self.log (console, "SAVING REPORT")  # DEBUG        
        reportObj = {'objects_created':[], 
                     'text_message': report}

        reportObj['objects_created'].append({'ref':readsSet_ref,
                                             'description':params['desc']})

        if unpaired_fwd_ref != None:
            reportObj['objects_created'].append({'ref':unpaired_fwd_ref,
                                                 'description':params['desc']+" unpaired fwd reads"})

        if unpaired_rev_ref != None:
            reportObj['objects_created'].append({'ref':unpaired_rev_ref,
                                                 'description':params['desc']+" unpaired rev reads"})


        # save report object
        #
        report = KBaseReport(self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
        report_info = report.create({'report':reportObj, 'workspace_name':params['workspace_name']})

        returnVal = { 'report_name': report_info['name'], 'report_ref': report_info['ref'] }
        #END KButil_Split_Reads

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method KButil_Split_Reads return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def KButil_Random_Subsample_Reads(self, ctx, params):
        """
        :param params: instance of type
           "KButil_Random_Subsample_Reads_Params" -> structure: parameter
           "workspace_name" of type "workspace_name" (** The workspace object
           refs are of form: ** **    objects = ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_ref" of type "data_obj_ref", parameter
           "output_name" of type "data_obj_name", parameter
           "subsample_fraction" of type "Fractionate_Options"
           (KButil_Random_Subsample_Reads() ** **  Method for random
           subsampling of reads library) -> structure: parameter "split_num"
           of Long, parameter "reads_num" of Long, parameter "reads_perc" of
           Double, parameter "desc" of String, parameter "seed" of Long
        :returns: instance of type "KButil_Random_Subsample_Reads_Output" ->
           structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN KButil_Random_Subsample_Reads
        console = []
        invalid_msgs = []
        self.log(console, 'Running KButil_Random_Subsample_Reads() with parameters: ')
        self.log(console, "\n"+pformat(params))
        report = ''
        
        token = ctx['token']
        wsClient = workspaceService(self.workspaceURL, token=token)
        headers = {'Authorization': 'OAuth '+token}
        env = os.environ.copy()
        env['KB_AUTH_TOKEN'] = token
        
        #SERVICE_VER = 'dev'  # DEBUG
        SERVICE_VER = 'release'

        # init randomizer
        if 'seed' in params and params['seed'] != None:
            random.seed(params['seed'])
        else:
            random.seed()

        # param checks
        required_params = ['workspace_name',
                           'input_ref', 
                           'output_name'
                           ]
        for required_param in required_params:
            if required_param not in params or params[required_param] == None:
                raise ValueError ("Must define required param: '"+required_param+"'")
            
#        # and param defaults
#        defaults = { 'split_num': 10
#                   }
#        for arg in defaults.keys():
#            if arg not in params or params[arg] == None or params[arg] == '':
#                params[arg] = defaults[arg]

        if 'subsample_fraction' not in params or params['subsample_fraction'] == None:
            raise ValueError ("Missing subsample_fraction params")
        if 'split_num' not in params['subsample_fraction'] or params['subsample_fraction']['split_num'] == None or params['subsample_fraction']['split_num'] < 0:
            raise ValueError ("Missing split_num")

        # use split_num to create reads_perc if neither reads_num or reads_perc defined
        use_reads_num  = False
        use_reads_perc = False
        if ('reads_num' in params['subsample_fraction'] and params['subsample_fraction']['reads_num'] != None and params['subsample_fraction']['reads_num'] > 0):
            self.log (console, "Ignoring reads_perc and just using reads_num: "+str(params['subsample_fraction']['reads_num']))
            use_reads_num  = True
            
        elif ('reads_perc' in params['subsample_fraction'] and params['subsample_fraction']['reads_perc'] != None and params['subsample_fraction']['reads_perc'] > 0 and params['subsample_fraction']['reads_perc'] <= 100):
            self.log (console, "Ignoring reads_num and just using reads_perc: "+str(params['subsample_fraction']['reads_perc']))
            use_reads_perc = True

        elif ('reads_num' not in params['subsample_fraction'] or params['subsample_fraction']['reads_num'] == None or params['subsample_fraction']['reads_num'] <= 0) \
                and ('reads_perc' not in params['subsample_fraction'] or params['subsample_fraction']['reads_perc'] == None or params['subsample_fraction']['reads_perc'] <= 0):

            params['subsample_fraction']['reads_perc'] = int(100.0 * 1.0/params['subsample_fraction']['split_num'])
            use_reads_perc = True

        else:
            raise ValueError ("Badly configured subsample_fraction params")
            

        # load provenance
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        provenance[0]['input_ws_objects']=[str(params['input_ref'])]


        # Determine whether read library is of correct type
        #
        try:
            # object_info tuple
            [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)
            
            input_reads_ref = params['input_ref']
            input_reads_obj_info = wsClient.get_object_info_new ({'objects':[{'ref':input_reads_ref}]})[0]
            input_reads_obj_name = input_reads_obj_info[NAME_I]
            input_reads_obj_type = input_reads_obj_info[TYPE_I]
            input_reads_obj_type = re.sub ('-[0-9]+\.[0-9]+$', "", input_reads_obj_type)  # remove trailing version
            #input_reads_obj_version = input_reads_obj_info[VERSION_I]  # this is object version, not type version

        except Exception as e:
            raise ValueError('Unable to get read library object info from workspace: (' + str(input_reads_ref) +')' + str(e))

        PE_types = ["KBaseFile.PairedEndLibrary"]
        SE_types = ["KBaseFile.SingleEndLibrary"]
        acceptable_types = PE_types + SE_types
        if input_reads_obj_type not in acceptable_types:
            raise ValueError ("Input reads of type: '"+input_reads_obj_type+"'.  Must be one of "+", ".join(acceptable_types))


        # divide reads_num by 2 if paired end library because only counting pairs
        #
        if input_reads_obj_type in PE_types and 'reads_num' in params['subsample_fraction'] and params['subsample_fraction']['reads_num'] != None and params['subsample_fraction']['reads_num'] != '':
            params['subsample_fraction']['reads_num'] = int (params['subsample_fraction']['reads_num']/2.0 + 0.5)


        # Download Reads
        #
        self.log (console, "DOWNLOADING READS")  # DEBUG
        try:
            readsUtils_Client = ReadsUtils (url=self.callbackURL, token=ctx['token'])  # SDK local
        except Exception as e:
            raise ValueError('Unable to get ReadsUtils Client' +"\n" + str(e))
        try:
            readsLibrary = readsUtils_Client.download_reads ({'read_libraries': [input_reads_ref],
                                                             'interleaved': 'false'
                                                             })
        except Exception as e:
            raise ValueError('Unable to download read library sequences from workspace: (' + str(input_reads_ref) +")\n" + str(e))


        # Paired End
        #
        if input_reads_obj_type == "KBaseFile.PairedEndLibrary":

            # Download reads Libs to FASTQ files
            input_fwd_file_path = readsLibrary['files'][input_reads_ref]['files']['fwd']
            input_rev_file_path = readsLibrary['files'][input_reads_ref]['files']['rev']
            sequencing_tech     = readsLibrary['files'][input_reads_ref]['sequencing_tech']

            # file paths
            input_fwd_path = re.sub ("\.fastq$", "", input_fwd_file_path)
            input_fwd_path = re.sub ("\.FASTQ$", "", input_fwd_path)
            input_rev_path = re.sub ("\.fastq$", "", input_rev_file_path)
            input_rev_path = re.sub ("\.FASTQ$", "", input_rev_path)
            output_fwd_paired_file_path_base   = input_fwd_path+"_fwd_paired"
            output_rev_paired_file_path_base   = input_rev_path+"_rev_paired"
            # set up for file io
            total_paired_reads = 0
            total_unpaired_fwd_reads = 0
            total_unpaired_rev_reads = 0
            total_paired_reads_by_set = []
            fwd_ids = dict()
            paired_ids = dict()
            paired_ids_list = []
            paired_lib_i = dict()
            paired_buf_size = 100000
            recs_beep_n = 1000000

            # read fwd file to get fwd ids
#            rec_cnt = 0  # DEBUG
            self.log (console, "GETTING IDS")  # DEBUG
            with open (input_fwd_file_path, 'r', 0) as input_reads_file_handle:
                rec_line_i = -1
                for line in input_reads_file_handle:
                    rec_line_i += 1
                    if rec_line_i == 3:
                        rec_line_i = -1
                    elif rec_line_i == 0:
                        if not line.startswith('@'):
                            raise ValueError ("badly formatted rec line: '"+line+"'")
                        read_id = line.rstrip('\n')
                        read_id = re.sub ("[ \t]+.*$", "", read_id)
                        read_id = re.sub ("[\/\.\_\-\:\;][012lrLRfrFR53]\'*$", "", read_id)
                        fwd_ids[read_id] = True

                        # DEBUG
#                        if rec_cnt % 100 == 0:
#                            self.log(console,"read_id: '"+str(read_id)+"'")
#                        rec_cnt += 1 


            # read reverse to determine paired
            self.log (console, "DETERMINING PAIRED IDS")  # DEBUG
            with open (input_rev_file_path, 'r', 0) as input_reads_file_handle:
                rec_line_i = -1
                for line in input_reads_file_handle:
                    rec_line_i += 1
                    if rec_line_i == 3:
                        rec_line_i = -1
                    elif rec_line_i == 0:
                        if not line.startswith('@'):
                            raise ValueError ("badly formatted rec line: '"+line+"'")
                        read_id = line.rstrip('\n')
                        read_id = re.sub ("[ \t]+.*$", "", read_id)
                        read_id = re.sub ("[\/\.\_\-\:\;][012lrLRfrFR53]\'*$", "", read_id)
                        if fwd_ids[read_id]:
                            paired_ids[read_id] = True
                            paired_ids_list.append(read_id)

                        # DEBUG
#                        if rec_cnt % 100 == 0:
#                            self.log(console,"read_id: '"+str(read_id)+"'")
#                        rec_cnt += 1 
            total_paired_reads = len(paired_ids_list)
            self.log (console, "TOTAL PAIRED READS CNT: "+str(total_paired_reads))  # DEBUG


            # Determine sublibrary sizes
            if use_reads_num:
                reads_per_lib = params['subsample_fraction']['reads_num']
                if reads_per_lib > total_paired_reads // params['subsample_fraction']['split_num']:
                    raise ValueError ("must specify reads_num <= total_paired_reads_cnt / split_num.  You have reads_num:"+str(params['subsample_fraction']['reads_num'])+" > total_paired_reads_cnt:"+str(total_paired_reads)+" / split_num:"+str(params['subsample_fraction']['split_num'])+".  Instead try reads_num <= "+str(total_paired_reads // params['subsample_fraction']['split_num']))
            elif use_reads_perc:
                reads_per_lib = int ((params['subsample_fraction']['reads_perc']/100.0) * total_paired_reads)
                if reads_per_lib > total_paired_reads // params['subsample_fraction']['split_num']:
                    raise ValueError ("must specify reads_perc <= 1 / split_num.  You have reads_perc:"+str(params['subsample_fraction']['reads_perc'])+" > 1 / split_num:"+str(params['subsample_fraction']['split_num'])+".  Instead try reads_perc <= "+ str(int(100 * 1/params['subsample_fraction']['split_num'])))
            else:
                raise ValueError ("error in logic reads_num vs. reads_perc logic")

            
            # Determine random membership in each sublibrary
            self.log (console, "GETTING RANDOM SUBSAMPLES")  # DEBUG

            for i,read_id in enumerate(random.sample (paired_ids_list, reads_per_lib * params['subsample_fraction']['split_num'])):
                lib_i = i % params['subsample_fraction']['split_num']
                paired_lib_i[read_id] = lib_i


            # split fwd paired
            self.log (console, "WRITING FWD SPLIT PAIRED")  # DEBUG
            paired_output_reads_file_handles = []
            for lib_i in range(params['subsample_fraction']['split_num']):
                paired_output_reads_file_handles.append(open (output_fwd_paired_file_path_base+"-"+str(lib_i)+".fastq", 'w', paired_buf_size))
                total_paired_reads_by_set.append(0)

            rec_buf = []
            last_read_id = None
            paired_cnt = 0
            capture_type_paired = False

            with open (input_fwd_file_path, 'r', 0) as input_reads_file_handle:
                rec_line_i = -1
                for line in input_reads_file_handle:
                    rec_line_i += 1
                    if rec_line_i == 3:
                        rec_line_i = -1
                    elif rec_line_i == 0:
                        if not line.startswith('@'):
                            raise ValueError ("badly formatted rec line: '"+line+"'")
                        if last_read_id != None:
                            if capture_type_paired:
                                lib_i = paired_lib_i[last_read_id]
                                paired_output_reads_file_handles[lib_i].writelines(rec_buf)
                                paired_cnt += 1
                                total_paired_reads_by_set[lib_i] += 1
                                if paired_cnt != 0 and paired_cnt % recs_beep_n == 0:
                                    self.log(console,"\t"+str(paired_cnt)+" recs processed")
                            else:
                                #unpaired_fwd_buf.extend(rec_buf)
                                pass
                            rec_buf = []
                        read_id = line.rstrip('\n')
                        read_id = re.sub ("[ \t]+.*$", "", read_id)
                        read_id = re.sub ("[\/\.\_\-\:\;][012lrLRfrFR53]\'*$", "", read_id)
                        last_read_id = read_id
                        try:
                            found = paired_lib_i[read_id]
                            capture_type_paired = True
                        except:
                            total_unpaired_fwd_reads += 1
                            capture_type_paired = False
                    rec_buf.append(line)
                # last rec
                if len(rec_buf) > 0:
                    if capture_type_paired:
                        lib_i = paired_lib_i[last_read_id]
                        paired_output_reads_file_handles[lib_i].writelines(rec_buf)
                        paired_cnt += 1
                        if paired_cnt != 0 and paired_cnt % recs_beep_n == 0:
                            self.log(console,"\t"+str(paired_cnt)+" recs processed")
                    else:
                        #unpaired_fwd_buf.extend(rec_buf)
                        pass
                    rec_buf = []

            for output_handle in paired_output_reads_file_handles:
                output_handle.close()

            self.log(console,"\t"+str(paired_cnt)+" FWD recs processed")


            # split rev paired
            self.log (console, "WRITING REV SPLIT PAIRED")  # DEBUG
            paired_output_reads_file_handles = []
            for lib_i in range(params['subsample_fraction']['split_num']):
                paired_output_reads_file_handles.append(open (output_rev_paired_file_path_base+"-"+str(lib_i)+".fastq", 'w', paired_buf_size))

            rec_buf = []
            last_read_id = None
            paired_cnt = 0
            capture_type_paired = False

            with open (input_rev_file_path, 'r', 0) as input_reads_file_handle:
                rec_line_i = -1
                for line in input_reads_file_handle:
                    rec_line_i += 1
                    if rec_line_i == 3:
                        rec_line_i = -1
                    elif rec_line_i == 0:
                        if not line.startswith('@'):
                            raise ValueError ("badly formatted rec line: '"+line+"'")
                        if last_read_id != None:
                            if capture_type_paired:
                                lib_i = paired_lib_i[last_read_id]
                                paired_output_reads_file_handles[lib_i].writelines(rec_buf)
                                paired_cnt += 1
                                if paired_cnt != 0 and paired_cnt % recs_beep_n == 0:
                                    self.log(console,"\t"+str(paired_cnt)+" recs processed")
                            else:
                                #unpaired_fwd_buf.extend(rec_buf)
                                pass
                            rec_buf = []
                        read_id = line.rstrip('\n')
                        read_id = re.sub ("[ \t]+.*$", "", read_id)
                        read_id = re.sub ("[\/\.\_\-\:\;][012lrLRfrFR53]\'*$", "", read_id)
                        last_read_id = read_id
                        try:
                            found = paired_lib_i[read_id]
                            capture_type_paired = True
                        except:
                            total_unpaired_rev_reads += 1
                            capture_type_paired = False
                    rec_buf.append(line)
                # last rec
                if len(rec_buf) > 0:
                    if capture_type_paired:
                        lib_i = paired_lib_i[last_read_id]
                        paired_output_reads_file_handles[lib_i].writelines(rec_buf)
                        paired_cnt += 1
                        if paired_cnt != 0 and paired_cnt % recs_beep_n == 0:
                            self.log(console,"\t"+str(paired_cnt)+" recs processed")
                    else:
                        #unpaired_fwd_buf.extend(rec_buf)
                        pass
                    rec_buf = []

            for output_handle in paired_output_reads_file_handles:
                output_handle.close()

            self.log(console,"\t"+str(paired_cnt)+" REV recs processed")


            # store report
            #
            report += "TOTAL PAIRED READS: "+str(total_paired_reads)+"\n"
            report += "TOTAL UNPAIRED FWD READS (discarded): "+str(total_unpaired_fwd_reads)+"\n"
            report += "TOTAL UNPAIRED REV READS (discarded): "+str(total_unpaired_rev_reads)+"\n"
            report += "\n"
            for lib_i in range(params['subsample_fraction']['split_num']):
                report += "PAIRED READS IN SET "+str(lib_i)+": "+str(total_paired_reads_by_set[lib_i])+"\n"


            # upload paired reads
            #
            self.log (console, "UPLOAD PAIRED READS LIBS")  # DEBUG
            paired_obj_refs = []
            for lib_i in range(params['subsample_fraction']['split_num']):
                output_fwd_paired_file_path = output_fwd_paired_file_path_base+"-"+str(lib_i)+".fastq"
                output_rev_paired_file_path = output_rev_paired_file_path_base+"-"+str(lib_i)+".fastq"
                if not os.path.isfile (output_fwd_paired_file_path) \
                     or os.path.getsize (output_fwd_paired_file_path) == 0 \
                   or not os.path.isfile (output_rev_paired_file_path) \
                     or os.path.getsize (output_rev_paired_file_path) == 0:
                    
                    raise ValueError ("failed to create paired output")
                else:
                    output_obj_name = params['output_name']+'_paired-'+str(lib_i)
                    self.log(console, 'Uploading paired reads: '+output_obj_name)
                    paired_obj_refs.append (readsUtils_Client.upload_reads ({ 'wsname': str(params['workspace_name']),
                                                                              'name': output_obj_name,
                                                                              # remove sequencing_tech when source_reads_ref is working
                                                                              #'sequencing_tech': sequencing_tech,
                                                                              'source_reads_ref': input_reads_ref,
                                                                              'fwd_file': output_fwd_paired_file_path,
                                                                              'rev_file': output_rev_paired_file_path
                                                                              })['obj_ref'])
                    
                

        # SingleEndLibrary
        #
        elif input_reads_obj_type == "KBaseFile.SingleEndLibrary":
            self.log(console, "Downloading Single End reads file...")

            # Download reads Libs to FASTQ files
            input_fwd_file_path = readsLibrary['files'][input_reads_ref]['files']['fwd']
            sequencing_tech     = readsLibrary['files'][input_reads_ref]['sequencing_tech']

            # file paths
            input_fwd_path = re.sub ("\.fastq$", "", input_fwd_file_path)
            input_fwd_path = re.sub ("\.FASTQ$", "", input_fwd_path)
            output_fwd_paired_file_path_base   = input_fwd_path+"_fwd_paired"

            # get "paired" ids
            self.log (console, "DETERMINING IDS")  # DEBUG
            paired_ids = dict()
            paired_ids_list = []
            paired_lib_i = dict()
            paired_buf_size = 100000
            recs_beep_n = 1000000

            with open (input_fwd_file_path, 'r', 0) as input_reads_file_handle:
                rec_line_i = -1
                for line in input_reads_file_handle:
                    rec_line_i += 1
                    if rec_line_i == 3:
                        rec_line_i = -1
                    elif rec_line_i == 0:
                        if not line.startswith('@'):
                            raise ValueError ("badly formatted rec line: '"+line+"'")
                        read_id = line.rstrip('\n')
                        read_id = re.sub ("[ \t]+.*$", "", read_id)
                        read_id = re.sub ("[\/\.\_\-\:\;][012lrLRfrFR53]\'*$", "", read_id)
                        if read_id in paired_ids:
                            self.log (console, "WARNING: repeat read_id "+read_id+" in reads library "+input_reads_obj_name)
                        paired_ids[read_id] = True
                        paired_ids_list.append(read_id)
                        # DEBUG
#                        if rec_cnt % 100 == 0:
#                            self.log(console,"read_id: '"+str(read_id)+"'")
#                        rec_cnt += 1 
            total_paired_reads = len(paired_ids_list)
            self.log (console, "TOTAL READS CNT: "+str(total_paired_reads))  # DEBUG


            # Determine sublibrary sizes
            if use_reads_num:
                reads_per_lib = params['subsample_fraction']['reads_num']
                if reads_per_lib > total_paired_reads // params['subsample_fraction']['split_num']:
                    raise ValueError ("must specify reads_num <= total_reads_cnt / split_num.  You have reads_num:"+str(params['subsample_fraction']['reads_num'])+" > total_reads_cnt:"+str(total_paired_reads)+" / split_num:"+str(params['subsample_fraction']['split_num'])+".  Instead try reads_num <= "+str(total_paired_reads // params['subsample_fraction']['split_num']))
            elif use_reads_perc:
                reads_per_lib = int ((params['subsample_fraction']['reads_perc']/100.0) * total_paired_reads)
                if reads_per_lib > total_paired_reads // params['subsample_fraction']['split_num']:
                    raise ValueError ("must specify reads_perc <= 1 / split_num.  You have reads_perc:"+str(params['subsample_fraction']['reads_perc'])+" > 1 / split_num:"+str(params['subsample_fraction']['split_num'])+".  Instead try reads_perc <= "+ str(int(100 * 1/params['subsample_fraction']['split_num'])))
            else:
                raise ValueError ("error in logic reads_num vs. reads_perc logic")

            
            # Determine random membership in each sublibrary
            self.log (console, "GETTING RANDOM SUBSAMPLES")  # DEBUG

            for i,read_id in enumerate(random.sample (paired_ids_list, reads_per_lib * params['subsample_fraction']['split_num'])):
                lib_i = i % params['subsample_fraction']['split_num']
                paired_lib_i[read_id] = lib_i


            # set up for file io
            total_paired_reads = 0
            total_paired_reads_by_set = []
            paired_buf_size = 1000000


            # split reads
            self.log (console, "WRITING SPLIT SINGLE END READS")  # DEBUG
            paired_output_reads_file_handles = []
            for lib_i in range(params['subsample_fraction']['split_num']):
                paired_output_reads_file_handles.append(open (output_fwd_paired_file_path_base+"-"+str(lib_i)+".fastq", 'w', paired_buf_size))
                total_paired_reads_by_set.append(0)

            rec_buf = []
            last_read_id = None
            paired_cnt = 0
            recs_beep_n = 1000000
            with open (input_fwd_file_path, 'r', 0) as input_reads_file_handle:
                rec_line_i = -1
                for line in input_reads_file_handle:
                    rec_line_i += 1
                    if rec_line_i == 3:
                        rec_line_i = -1
                    elif rec_line_i == 0:
                        if not line.startswith('@'):
                            raise ValueError ("badly formatted rec line: '"+line+"'")
                        total_paired_reads += 1
                        if last_read_id != None:
                            try:
                                lib_i = paired_lib_i[last_read_id]
                                total_paired_reads_by_set[lib_i] += 1
                                paired_output_reads_file_handles[lib_i].writelines(rec_buf)
                                paired_cnt += 1
                            except:
                                pass
                            if paired_cnt != 0 and paired_cnt % recs_beep_n == 0:
                                self.log(console,"\t"+str(paired_cnt)+" recs processed")
                            rec_buf = []
                        read_id = line.rstrip('\n')
                        read_id = re.sub ("[ \t]+.*$", "", read_id)
                        read_id = re.sub ("[\/\.\_\-\:\;][012lrLRfrFR53]\'*$", "", read_id)
                        last_read_id = read_id
                    rec_buf.append(line)
                # last rec
                if len(rec_buf) > 0:
                    if last_read_id != None:
                        try:
                            lib_i = paired_lib_i[last_read_id]
                            total_paired_reads_by_set[lib_i] += 1
                            paired_output_reads_file_handles[lib_i].writelines(rec_buf)
                            paired_cnt += 1
                        except:
                            pass
                    if paired_cnt != 0 and paired_cnt % recs_beep_n == 0:
                        self.log(console,"\t"+str(paired_cnt)+" recs processed")
                    rec_buf = []

            for output_handle in paired_output_reads_file_handles:
                output_handle.close()


            # store report
            #
            report += "TOTAL READS: "+str(total_paired_reads)+"\n"
            for lib_i in range(params['subsample_fraction']['split_num']):
                report += "SINGLE END READS IN SET "+str(lib_i)+": "+str(total_paired_reads_by_set[lib_i])+"\n"


            # upload reads
            #
            self.log (console, "UPLOADING SPLIT SINGLE END READS")  # DEBUG
            paired_obj_refs = []
            for lib_i in range(params['subsample_fraction']['split_num']):
                output_fwd_paired_file_path = output_fwd_paired_file_path_base+"-"+str(lib_i)+".fastq"
                if not os.path.isfile (output_fwd_paired_file_path) \
                        or os.path.getsize (output_fwd_paired_file_path) == 0:
                    
                    raise ValueError ("failed to create single end library output")
                else:
                    output_obj_name = params['output_name']+'-'+str(lib_i)
                    self.log(console, 'Uploading single end reads: '+output_obj_name)
                    paired_obj_refs.append( readsUtils_Client.upload_reads ({ 'wsname': str(params['workspace_name']),
                                                                              'name': output_obj_name,
                                                                              # remove sequencing_tech when source_reads_ref is working
                                                                              #'sequencing_tech': sequencing_tech,
                                                                              'source_reads_ref': input_reads_ref,
                                                                              'fwd_file': output_fwd_paired_file_path
                                                                              })['obj_ref'])
                                            
        else:
            raise ValueError ("unknown ReadLibrary type as input: "+str(input_reads_obj_type))


        # save output readsSet
        #
        self.log (console, "SAVING READSSET")  # DEBUG
        items = []
        for lib_i,lib_ref in enumerate(paired_obj_refs):
            label = params['output_name']+'-'+str(lib_i)
            items.append({'ref': lib_ref,
                          'label': label
                          #'data_attachment': ,
                          #'info':
                              })
        description = params['desc']
        output_readsSet_obj = { 'description': params['desc'],
                                'items': items
                                }
        output_readsSet_name = str(params['output_name'])
        setAPI_Client = SetAPI (url=self.serviceWizardURL, token=ctx['token'])  # for dynamic service
        readsSet_ref = setAPI_Client.save_reads_set_v1 ({'workspace_name': params['workspace_name'],
                                                         'output_object_name': output_readsSet_name,
                                                         'data': output_readsSet_obj
                                                         })['set_ref']
                              

        # build report
        #
        self.log (console, "SAVING REPORT")  # DEBUG        
        reportObj = {'objects_created':[], 
                     'text_message': report}

        reportObj['objects_created'].append({'ref':readsSet_ref,
                                             'description':params['desc']})


        # save report object
        #
        report = KBaseReport(self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
        report_info = report.create({'report':reportObj, 'workspace_name':params['workspace_name']})

        returnVal = { 'report_name': report_info['name'], 'report_ref': report_info['ref'] }
        #END KButil_Random_Subsample_Reads

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method KButil_Random_Subsample_Reads return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def KButil_Merge_ReadsSet_to_OneLibrary(self, ctx, params):
        """
        :param params: instance of type
           "KButil_Merge_ReadsSet_to_OneLibrary_Params"
           (KButil_Merge_ReadsSet_to_OneLibrary() ** **  Method for merging a
           ReadsSet into one library) -> structure: parameter
           "workspace_name" of type "workspace_name" (** The workspace object
           refs are of form: ** **    objects = ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_ref" of type "data_obj_ref", parameter
           "output_name" of type "data_obj_name", parameter "desc" of String
        :returns: instance of type
           "KButil_Merge_ReadsSet_to_OneLibrary_Output" -> structure:
           parameter "report_name" of type "data_obj_name", parameter
           "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN KButil_Merge_ReadsSet_to_OneLibrary
        console = []
        report = ''
        self.log(console, 'Running KButil_Merge_ReadsSet_to_OneLibrary with parameters: ')
        self.log(console, "\n"+pformat(params))

        token = ctx['token']
        wsClient = workspaceService(self.workspaceURL, token=token)
        headers = {'Authorization': 'OAuth '+token}
        env = os.environ.copy()
        env['KB_AUTH_TOKEN'] = token

        #SERVICE_VER = 'dev'  # DEBUG
        SERVICE_VER = 'release'

        # param checks
        required_params = ['workspace_name',
                           'input_ref', 
                           'output_name'
                           ]
        for required_param in required_params:
            if required_param not in params or params[required_param] == None:
                raise ValueError ("Must define required param: '"+required_param+"'")
            
        # load provenance
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects']=[str(params['input_ref'])]

        # Determine whether read library or read set is input object
        #
        try:
            # object_info tuple
            [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)

            input_reads_obj_info = wsClient.get_object_info_new ({'objects':[{'ref':params['input_ref']}]})[0]
            input_reads_obj_type = input_reads_obj_info[TYPE_I]
            input_reads_obj_type = re.sub ('-[0-9]+\.[0-9]+$', "", input_reads_obj_type)  # remove trailing version
            #input_reads_obj_version = input_reads_obj_info[VERSION_I]  # this is object version, not type version

        except Exception as e:
            raise ValueError('Unable to get read library object from workspace: (' + str(params['input_ref']) +')' + str(e))

        acceptable_types = ["KBaseSets.ReadsSet"]
        if input_reads_obj_type not in acceptable_types:
            raise ValueError ("Input reads of type: '"+input_reads_obj_type+"'.  Must be one of "+", ".join(acceptable_types))


        # get set
        #
        readsSet_ref_list = []
        readsSet_names_list = []
        try:
            setAPI_Client = SetAPI (url=self.serviceWizardURL, token=ctx['token'])  # for dynamic service
            input_readsSet_obj = setAPI_Client.get_reads_set_v1 ({'ref':params['input_ref'],'include_item_info':1})
        except Exception as e:
            raise ValueError('SetAPI FAILURE: Unable to get read library set object from workspace: (' + str(params['input_ref'])+")\n" + str(e))

        for readsLibrary_obj in input_readsSet_obj['data']['items']:
            readsSet_ref_list.append(readsLibrary_obj['ref'])
            NAME_I = 1
            readsSet_names_list.append(readsLibrary_obj['info'][NAME_I])


        # check type of readsLibrary memebers of set
        #
        report = ''
        read_library_type = None
        for input_reads_library_ref in readsSet_ref_list:

            # make sure library types are consistent
            try:
                # object_info tuple
                [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)

                input_reads_obj_info = wsClient.get_object_info_new ({'objects':[{'ref':input_reads_library_ref}]})[0]
                input_reads_obj_type = input_reads_obj_info[TYPE_I]
                input_reads_obj_type = re.sub ('-[0-9]+\.[0-9]+$', "", input_reads_obj_type)  # remove trailing version

            except Exception as e:
                raise ValueError('Unable to get read library object from workspace: (' + str(params['input_ref']) +')' + str(e))

            acceptable_types = ["KBaseFile.PairedEndLibrary", "KBaseFile.SingleEndLibrary"]
            if input_reads_obj_type not in acceptable_types:
                raise ValueError ("Input reads of type: '"+input_reads_obj_type+"'.  Must be one of "+", ".join(acceptable_types))

            if read_library_type == None:
                read_library_type = input_reads_obj_type
            elif input_reads_obj_type != read_library_type:
                raise ValueError ("incompatible read library types in ReadsSet "+params['input_ref'])
            
        # combine read libraries
        #
        self.log (console, "CREATING COMBINED INPUT FASTQ FILES")

        # make dir
        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()*1000)
        input_dir = os.path.join(self.scratch,'input.'+str(timestamp))
        if not os.path.exists(input_dir):
            os.makedirs(input_dir)

        # connect to ReadsUtils Client
        try:
            readsUtils_Client = ReadsUtils (url=self.callbackURL, token=ctx['token'])  # SDK local
        except:
            raise ValueError("Unable to get readsUtils_Client\n" + str(e))

        # start combined file
        read_buf_size  = 65536
        write_buf_size = 65536
        combined_input_fwd_path = os.path.join (input_dir, 'input_reads_fwd.fastq')
        combined_input_rev_path = os.path.join (input_dir, 'input_reads_rev.fastq')
        combined_input_fwd_handle = open (combined_input_fwd_path, 'w', write_buf_size)
        combined_input_rev_handle = open (combined_input_rev_path, 'w', write_buf_size)


        # add libraries, one at a time
        sequencing_tech = None
        for this_input_reads_ref in readsSet_ref_list:
            self.log (console, "DOWNLOADING FASTQ FILES FOR ReadsSet member: "+str(this_input_reads_ref))
            try:
                readsLibrary = readsUtils_Client.download_reads ({'read_libraries': [this_input_reads_ref],
                                                                  'interleaved': 'false'
                                                                  })
            except Exception as e:
                raise ValueError('Unable to get reads object from workspace: (' + this_input_reads_ref +")\n" + str(e))

            this_input_fwd_path = readsLibrary['files'][this_input_reads_ref]['files']['fwd']

            if input_reads_obj_type == "KBaseFile.PairedEndLibrary":
                this_input_rev_path = readsLibrary['files'][this_input_reads_ref]['files']['rev']

            this_sequencing_tech = readsLibrary['files'][this_input_reads_ref]['sequencing_tech']
            if sequencing_tech == None:
                sequencing_tech = this_sequencing_tech
            elif this_sequencing_tech != sequencing_tech:
                sequencing_tech = 'N/A'


            # append fwd
            self.log (console, "APPENDING FASTQ FILES FOR ReadsSet member: "+str(this_input_reads_ref))
            this_input_path = this_input_fwd_path
            cat_file_handle = combined_input_fwd_handle
            with open (this_input_path, 'r', read_buf_size) as this_input_handle:
                while True:
                    read_data = this_input_handle.read(read_buf_size)
                    if read_data:
                        cat_file_handle.write(read_data)
                    else:
                        break
            os.remove (this_input_path)  # create space since we no longer need the piece file

            # append rev
            if input_reads_obj_type == "KBaseFile.PairedEndLibrary":
                this_input_path = this_input_rev_path
                cat_file_handle = combined_input_rev_handle
                with open (this_input_path, 'r', read_buf_size) as this_input_handle:
                    while True:
                        read_data = this_input_handle.read(read_buf_size)
                        if read_data:
                            cat_file_handle.write(read_data)
                        else:
                            break
                os.remove (this_input_path)  # create space since we no longer need the piece file

        combined_input_fwd_handle.close()
        combined_input_rev_handle.close()


        # upload reads
        #
        self.log (console, "UPLOADING MERGED READS LIB")  # DEBUG
        if not os.path.isfile (combined_input_fwd_path) \
                or os.path.getsize (combined_input_fwd_path) == 0:
            raise ValueError ("failed to create fwd read library output")
        if input_reads_obj_type == "KBaseFile.PairedEndLibrary":
            if not os.path.isfile (combined_input_rev_path) \
                or os.path.getsize (combined_input_rev_path) == 0:
                    
                raise ValueError ("failed to create rev read library output")

        output_obj_name = params['output_name']
        self.log(console, 'Uploading reads library: '+output_obj_name)

        if input_reads_obj_type == "KBaseFile.PairedEndLibrary":
            reads_library_ref = readsUtils_Client.upload_reads ({ 'wsname': str(params['workspace_name']),
                                                                  'name': output_obj_name,
                                                                  # remove sequencing_tech when source_reads_ref is working
                                                                  #'sequencing_tech': sequencing_tech,
                                                                  'source_reads_ref': readsSet_ref_list[0],
                                                                  'single_genome': 0,
                                                                  'fwd_file': combined_input_fwd_path,
                                                                  'rev_file': combined_input_rev_path
                                                                  })['obj_ref']
        else:
            reads_library_ref = readsUtils_Client.upload_reads ({ 'wsname': str(params['workspace_name']),
                                                                  'name': output_obj_name,
                                                                  # remove sequencing_tech when source_reads_ref is working
                                                                  #'sequencing_tech': sequencing_tech,
                                                                  'source_reads_ref': readsSet_ref_list[0],
                                                                  'single_genome': 0,
                                                                  'fwd_file': combined_input_fwd_path,
                                                                  })['obj_ref']
            

        # build report message
        report += "NUM READS LIBRARIES COMBINED INTO ONE READS LIBRARY: " + str(len(readsSet_ref_list))+"\n"

        # build report
        #
        self.log (console, "SAVING REPORT")  # DEBUG        
        reportObj = {'objects_created':[], 
                     'text_message': report}

        reportObj['objects_created'].append({'ref':reads_library_ref,
                                             'description':params['desc']})


        # save report object
        #
        report = KBaseReport(self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
        report_info = report.create({'report':reportObj, 'workspace_name':params['workspace_name']})

        returnVal = { 'report_name': report_info['name'], 'report_ref': report_info['ref'] }
        #END KButil_Merge_ReadsSet_to_OneLibrary

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method KButil_Merge_ReadsSet_to_OneLibrary return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def KButil_Merge_MultipleReadsLibs_to_OneLibrary(self, ctx, params):
        """
        :param params: instance of type
           "KButil_Merge_MultipleReadsLibs_to_OneLibrary_Params"
           (KButil_Merge_MultipleReadsLibs_to_OneLibrary() ** **  Method for
           merging ReadsLibs into one library) -> structure: parameter
           "workspace_name" of type "workspace_name" (** The workspace object
           refs are of form: ** **    objects = ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_refs" of type "data_obj_ref", parameter
           "output_name" of type "data_obj_name", parameter "desc" of String
        :returns: instance of type
           "KButil_Merge_MultipleReadsLibs_to_OneLibrary_Output" ->
           structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN KButil_Merge_MultipleReadsLibs_to_OneLibrary
        console = []
        invalid_msgs = []
        report = ''
        self.log(console, 'Running KButil_Merge_MultipleReadsLibs_to_OneLibrary with parameters: ')
        self.log(console, "\n"+pformat(params))

        token = ctx['token']
        wsClient = workspaceService(self.workspaceURL, token=token)
        headers = {'Authorization': 'OAuth '+token}
        env = os.environ.copy()
        env['KB_AUTH_TOKEN'] = token

        #SERVICE_VER = 'dev'  # DEBUG
        SERVICE_VER = 'release'

        # param checks
        required_params = ['workspace_name',
                           'input_refs', 
                           'output_name'
                           ]
        for required_param in required_params:
            if required_param not in params or params[required_param] == None:
                raise ValueError ("Must define required param: '"+required_param+"'")
            
        # clean input_refs
        clean_input_refs = []
        for ref in params['input_refs']:
            if ref != None and ref != '':
                clean_input_refs.append(ref)
        params['input_refs'] = clean_input_refs

        if len(params['input_refs']) < 2:
            self.log(console,"Must provide at least two ReadsLibs or ReadsSets")
            self.log(invalid_msgs,"Must provide at least two ReadsLibs or ReadsSets")

        # load provenance
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects']=[]
        for input_ref in params['input_refs']:
            provenance[0]['input_ws_objects'].append(input_ref)

        # get set
        #
        readsSet_ref_list = []
        readsSet_names_list = []
        for reads_ref in params['input_refs']:
            try:
                # object_info tuple
                [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)
                
                input_reads_obj_info = wsClient.get_object_info_new ({'objects':[{'ref':reads_ref}]})[0]
                input_reads_obj_type = input_reads_obj_info[TYPE_I]
                input_reads_obj_type = re.sub ('-[0-9]+\.[0-9]+$', "", input_reads_obj_type)  # remove trailing version
                #input_reads_obj_version = input_reads_obj_info[VERSION_I]  # this is object version, not type version

            except Exception as e:
                raise ValueError('Unable to get read library object from workspace: (' + str(reads_ref) +')' + str(e))

            acceptable_types = ["KBaseSets.ReadsSet", "KBaseFile.PairedEndLibrary", "KBaseFile.SingleEndLibrary"]
            if input_reads_obj_type not in acceptable_types:
                raise ValueError ("Input reads of type: '"+input_reads_obj_type+"'.  Must be one of "+", ".join(acceptable_types))

            if input_reads_obj_type != "KBaseSets.ReadsSet":  # readsLib
                readsSet_ref_list.append(reads_ref)

            else:  # readsSet
                try:
                    setAPI_Client = SetAPI (url=self.serviceWizardURL, token=ctx['token'])  # for dynamic service
                    input_readsSet_obj = setAPI_Client.get_reads_set_v1 ({'ref':reads_ref,'include_item_info':1})
                except Exception as e:
                    raise ValueError('SetAPI FAILURE: Unable to get read library set object from workspace: (' + str(reads_ref)+")\n" + str(e))
                
                for readsLibrary_obj in input_readsSet_obj['data']['items']:
                    readsSet_ref_list.append(readsLibrary_obj['ref'])
#                    NAME_I = 1
#                    readsSet_names_list.append(readsLibrary_obj['info'][NAME_I])


        # check type of readsLibrary memebers of set
        #
        report = ''
        read_library_type = None
        for input_reads_library_ref in readsSet_ref_list:

            # make sure library types are consistent
            try:
                # object_info tuple
                [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)

                input_reads_obj_info = wsClient.get_object_info_new ({'objects':[{'ref':input_reads_library_ref}]})[0]
                input_reads_obj_type = input_reads_obj_info[TYPE_I]
                input_reads_obj_type = re.sub ('-[0-9]+\.[0-9]+$', "", input_reads_obj_type)  # remove trailing version
                readsSet_names_list.append(input_reads_obj_info[NAME_I])

            except Exception as e:
                raise ValueError('Unable to get read library object from workspace: (' + input_reads_library_ref +')' + str(e))

            acceptable_types = ["KBaseFile.PairedEndLibrary", "KBaseFile.SingleEndLibrary"]
            if input_reads_obj_type not in acceptable_types:
                raise ValueError ("Input reads of type: '"+input_reads_obj_type+"'.  Must be one of "+", ".join(acceptable_types))

            if read_library_type == None:
                read_library_type = input_reads_obj_type
            elif input_reads_obj_type != read_library_type:
                raise ValueError ("incompatible read library types in ReadsSet "+input_reads_library_ref)
            
        # combine read libraries
        #
        self.log (console, "CREATING COMBINED INPUT FASTQ FILES")

        # make dir
        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()*1000)
        input_dir = os.path.join(self.scratch,'input.'+str(timestamp))
        if not os.path.exists(input_dir):
            os.makedirs(input_dir)

        # connect to ReadsUtils Client
        try:
            readsUtils_Client = ReadsUtils (url=self.callbackURL, token=ctx['token'])  # SDK local
        except:
            raise ValueError("Unable to get readsUtils_Client\n" + str(e))

        # start combined file
        read_buf_size  = 65536
        write_buf_size = 65536
        combined_input_fwd_path = os.path.join (input_dir, 'input_reads_fwd.fastq')
        combined_input_rev_path = os.path.join (input_dir, 'input_reads_rev.fastq')
        combined_input_fwd_handle = open (combined_input_fwd_path, 'w', write_buf_size)
        combined_input_rev_handle = open (combined_input_rev_path, 'w', write_buf_size)


        # add libraries, one at a time
        sequencing_tech = None
        for this_input_reads_ref in readsSet_ref_list:
            clean_ref = re.sub("\/", "_", this_input_reads_ref)
            
            self.log (console, "DOWNLOADING FASTQ FILES FOR ReadsSet member: "+str(this_input_reads_ref))
            try:
                readsLibrary = readsUtils_Client.download_reads ({'read_libraries': [this_input_reads_ref],
                                                                  'interleaved': 'false'
                                                                  })
            except Exception as e:
                raise ValueError('Unable to get reads object from workspace: (' + this_input_reads_ref +")\n" + str(e))

            this_input_fwd_path = readsLibrary['files'][this_input_reads_ref]['files']['fwd']

            if input_reads_obj_type == "KBaseFile.PairedEndLibrary":
                this_input_rev_path = readsLibrary['files'][this_input_reads_ref]['files']['rev']

            this_sequencing_tech = readsLibrary['files'][this_input_reads_ref]['sequencing_tech']
            if sequencing_tech == None:
                sequencing_tech = this_sequencing_tech
            elif this_sequencing_tech != sequencing_tech:
                sequencing_tech = 'N/A'


            # append fwd
            self.log (console, "APPENDING FASTQ FILES FOR ReadsSet member: "+str(this_input_reads_ref))
            this_input_path = this_input_fwd_path
            cat_file_handle = combined_input_fwd_handle
            with open (this_input_path, 'r', read_buf_size) as this_input_handle:
                #while True:
                #    read_data = this_input_handle.read(read_buf_size)
                #    if read_data:
                #        cat_file_handle.write(read_data)
                #    else:
                #        break
                rec_line_i = -1
                for line in this_input_handle:
                    rec_line_i += 1
                    if rec_line_i == 3:
                        rec_line_i = -1
                    elif rec_line_i == 0:
                        if not line.startswith('@'):
                            raise ValueError ("badly formatted rec line: '"+line+"'")
                        
                        line = '@'+clean_ref+':'+line[1:]
                    cat_file_handle.write(line)
            os.remove (this_input_path)  # create space since we no longer need the piece file

            # append rev
            if input_reads_obj_type == "KBaseFile.PairedEndLibrary":
                this_input_path = this_input_rev_path
                cat_file_handle = combined_input_rev_handle
                with open (this_input_path, 'r', read_buf_size) as this_input_handle:
                    #while True:
                    #    read_data = this_input_handle.read(read_buf_size)
                    #    if read_data:
                    #        cat_file_handle.write(read_data)
                    #    else:
                    #        break
                    rec_line_i = -1
                    for line in this_input_handle:
                        rec_line_i += 1
                        if rec_line_i == 3:
                            rec_line_i = -1
                        elif rec_line_i == 0:
                            if not line.startswith('@'):
                                raise ValueError ("badly formatted rec line: '"+line+"'")
                        
                            line = '@'+clean_ref+':'+line[1:]
                        cat_file_handle.write(line)
                os.remove (this_input_path)  # create space since we no longer need the piece file

        combined_input_fwd_handle.close()
        combined_input_rev_handle.close()


        # upload reads
        #
        self.log (console, "UPLOADING MERGED READS LIB")  # DEBUG
        if not os.path.isfile (combined_input_fwd_path) \
                or os.path.getsize (combined_input_fwd_path) == 0:
            raise ValueError ("failed to create fwd read library output")
        if input_reads_obj_type == "KBaseFile.PairedEndLibrary":
            if not os.path.isfile (combined_input_rev_path) \
                or os.path.getsize (combined_input_rev_path) == 0:
                    
                raise ValueError ("failed to create rev read library output")

        output_obj_name = params['output_name']
        self.log(console, 'Uploading reads library: '+output_obj_name)

        if input_reads_obj_type == "KBaseFile.PairedEndLibrary":
            reads_library_ref = readsUtils_Client.upload_reads ({ 'wsname': str(params['workspace_name']),
                                                                  'name': output_obj_name,
                                                                  # remove sequencing_tech when source_reads_ref is working
                                                                  #'sequencing_tech': sequencing_tech,
                                                                  'source_reads_ref': readsSet_ref_list[0],
                                                                  'single_genome': 0,
                                                                  'fwd_file': combined_input_fwd_path,
                                                                  'rev_file': combined_input_rev_path
                                                                  })['obj_ref']
        else:
            reads_library_ref = readsUtils_Client.upload_reads ({ 'wsname': str(params['workspace_name']),
                                                                  'name': output_obj_name,
                                                                  # remove sequencing_tech when source_reads_ref is working
                                                                  #'sequencing_tech': sequencing_tech,
                                                                  'source_reads_ref': readsSet_ref_list[0],
                                                                  'single_genome': 0,
                                                                  'fwd_file': combined_input_fwd_path,
                                                                  })['obj_ref']
            

        # build report message
        report += "NUM READS LIBRARIES COMBINED INTO ONE READS LIBRARY: " + str(len(readsSet_ref_list))+"\n"

        # build report
        #
        self.log (console, "SAVING REPORT")  # DEBUG        
        reportObj = {'objects_created':[], 
                     'text_message': report}

        reportObj['objects_created'].append({'ref':reads_library_ref,
                                             'description':params['desc']})


        # save report object
        #
        report = KBaseReport(self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
        report_info = report.create({'report':reportObj, 'workspace_name':params['workspace_name']})

        returnVal = { 'report_name': report_info['name'], 'report_ref': report_info['ref'] }
        #END KButil_Merge_MultipleReadsLibs_to_OneLibrary

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method KButil_Merge_MultipleReadsLibs_to_OneLibrary return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def KButil_Extract_Unpaired_Reads_and_Synchronize_Pairs(self, ctx, params):
        """
        :param params: instance of type
           "KButil_Extract_Unpaired_Reads_Params"
           (KButil_Extract_Unpaired_Reads_and_Synchronize_Pairs() ** ** 
           Method for removing unpaired reads from a paired end library or
           set and matching the order of reads) -> structure: parameter
           "workspace_name" of type "workspace_name" (** The workspace object
           refs are of form: ** **    objects = ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_ref" of type "data_obj_ref", parameter
           "output_name" of type "data_obj_name", parameter "desc" of String
        :returns: instance of type "KButil_Extract_Unpaired_Reads_Output" ->
           structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN KButil_Extract_Unpaired_Reads_and_Synchronize_Pairs
        console = []
        report = ''
        self.log(console, 'Running KButil_Extract_Unpaired_Reads_and_Synchronize_Pairs with parameters: ')
        self.log(console, "\n"+pformat(params))

        token = ctx['token']
        wsClient = workspaceService(self.workspaceURL, token=token)
        headers = {'Authorization': 'OAuth '+token}
        env = os.environ.copy()
        env['KB_AUTH_TOKEN'] = token

        #SERVICE_VER = 'dev'  # DEBUG
        SERVICE_VER = 'release'

        # param checks
        required_params = ['workspace_name',
                           'input_ref', 
                           'output_name'
                           ]
        for required_param in required_params:
            if required_param not in params or params[required_param] == None:
                raise ValueError ("Must define required param: '"+required_param+"'")
            
        # load provenance
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects']=params['input_ref']


        # Determine whether read library or read set is input object
        #
        try:
            # object_info tuple
            [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)

            input_reads_obj_info = wsClient.get_object_info_new ({'objects':[{'ref':params['input_ref']}]})[0]
            input_reads_obj_type = input_reads_obj_info[TYPE_I]
            input_reads_obj_type = re.sub ('-[0-9]+\.[0-9]+$', "", input_reads_obj_type)  # remove trailing version
            #input_reads_obj_version = input_reads_obj_info[VERSION_I]  # this is object version, not type version

        except Exception as e:
            raise ValueError('Unable to get read library object from workspace: (' + str(params['input_ref']) +')' + str(e))

        acceptable_types = ["KBaseSets.ReadsSet", "KBaseFile.PairedEndLibrary", "KBaseAssembly.PairedEndLibrary"]
        if input_reads_obj_type not in acceptable_types:
            raise ValueError ("Input reads of type: '"+input_reads_obj_type+"'.  Must be one of "+", ".join(acceptable_types))


        # get set
        #
        readsSet_ref_list = []
        readsSet_names_list = []
        if input_reads_obj_type != "KBaseSets.ReadsSet":
            readsSet_ref_list = [params['input_ref']]
        else:
            try:
                setAPI_Client = SetAPI (url=self.serviceWizardURL, token=ctx['token'])  # for dynamic service
                input_readsSet_obj = setAPI_Client.get_reads_set_v1 ({'ref':params['input_ref'],'include_item_info':1})

            except Exception as e:
                raise ValueError('SetAPI FAILURE: Unable to get read library set object from workspace: (' + str(params['input_ref'])+")\n" + str(e))
            for readsLibrary_obj in input_readsSet_obj['data']['items']:
                readsSet_ref_list.append(readsLibrary_obj['ref'])
                NAME_I = 1
                readsSet_names_list.append(readsLibrary_obj['info'][NAME_I])

        # Make sure all libraries are PairedEnd
        #
        if input_reads_obj_type == "KBaseSets.ReadsSet":
            for lib_i,input_reads_ref in enumerate(readsSet_ref_list):
                try:
                    # object_info tuple
                    [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)

                    this_input_reads_obj_info = wsClient.get_object_info_new ({'objects':[{'ref':input_reads_ref}]})[0]
                    this_input_reads_obj_type = this_input_reads_obj_info[TYPE_I]
                    this_input_reads_obj_type = re.sub ('-[0-9]+\.[0-9]+$', "", this_input_reads_obj_type)  # remove trailing version

                except Exception as e:
                    raise ValueError('Unable to get read library object from workspace: (' + input_reads_ref +')' + str(e))

                acceptable_types = ["KBaseFile.PairedEndLibrary", "KBaseAssembly.PairedEndLibrary"]
                if this_input_reads_obj_type not in acceptable_types:
                    raise ValueError ("Input reads in set at index "+str(lib_i)+" and name "+readSet_names_list[lib_i]+" of type: '"+this_input_reads_obj_type+"'.  Must be one of "+", ".join(acceptable_types))


        # Iterate through readsLibrary members of set
        #
        report = ''
        paired_readsSet_ref        = None
        unpaired_fwd_readsSet_ref  = None
        unpaired_rev_readsSet_ref  = None
        paired_readsSet_refs       = []
        unpaired_fwd_readsSet_refs = []
        unpaired_rev_readsSet_refs = []
        paired_obj_refs            = []
        unpaired_fwd_obj_refs      = []
        unpaired_rev_obj_refs      = []

        for lib_i,input_reads_ref in enumerate(readsSet_ref_list):

            # Download Reads
            #
            self.log (console, "DOWNLOADING READS FOR READ LIBRARY: "+input_reads_ref)  # DEBUG
            try:
                readsUtils_Client = ReadsUtils (url=self.callbackURL, token=ctx['token'])  # SDK local
            except Exception as e:
                raise ValueError('Unable to get ReadsUtils Client' +"\n" + str(e))
            try:
                readsLibrary = readsUtils_Client.download_reads ({'read_libraries': [input_reads_ref],
                                                                  'interleaved': 'false'
                                                                  })
            except Exception as e:
                raise ValueError('Unable to download read library sequences from workspace: (' + str(input_reads_ref) +")\n" + str(e))


            # Download reads Libs to FASTQ files
            input_fwd_file_path = readsLibrary['files'][input_reads_ref]['files']['fwd']
            input_rev_file_path = readsLibrary['files'][input_reads_ref]['files']['rev']
            sequencing_tech     = readsLibrary['files'][input_reads_ref]['sequencing_tech']

            # file paths
            input_fwd_path = re.sub ("\.fastq$", "", input_fwd_file_path)
            input_fwd_path = re.sub ("\.FASTQ$", "", input_fwd_path)
            input_rev_path = re.sub ("\.fastq$", "", input_rev_file_path)
            input_rev_path = re.sub ("\.FASTQ$", "", input_rev_path)
            output_fwd_paired_file_path_base   = input_fwd_path+"_fwd_paired"
            output_rev_paired_file_path_base   = input_rev_path+"_rev_paired"
            output_fwd_unpaired_file_path_base = input_fwd_path+"_fwd_unpaired"
            output_rev_unpaired_file_path_base = input_rev_path+"_rev_unpaired"


            # set up for file io
            paired_read_cnt = 0
            unpaired_fwd_read_cnt = 0
            unpaired_rev_read_cnt = 0
            total_fwd_recs = 0
            total_rev_recs = 0
            pair_ids = dict()
            rev_ids = dict()
            fwd_id_pos = dict()
            rev_id_pos = dict()
            pair_ids_order = []
            rev_ids_order = []
            unpaired_buf_size = 100000
            paired_buf_size = 100000
            recs_beep_n = 1000000

            # read rev file to get rev ids and order
            rec_cnt = 0
            self.log (console, "GETTING REV IDS")  # DEBUG
            with open (input_rev_file_path, 'r', 0) as input_reads_file_handle:
                rec_line_i = -1
                for line in input_reads_file_handle:
                    rec_line_i += 1
                    if rec_line_i == 3:
                        rec_line_i = -1
                    elif rec_line_i == 0:
                        if not line.startswith('@'):
                            raise ValueError ("badly formatted rec line: '"+line+"'")
                        rec_cnt += 1 
                        read_id = line.rstrip('\n')
                        read_id = re.sub ("[ \t]+.*$", "", read_id)
                        read_id = re.sub ("[\/\.\_\-\:\;][012lrLRfrFR53]\'*$", "", read_id)
                        rev_ids[read_id] = True
                        rev_ids_order.append(read_id)
                        # DEBUG
#                        if rec_cnt % 100 == 0:
#                            self.log(console,"read_id: '"+str(read_id)+"'")
            total_rev_recs = rec_cnt

            # read fwd file to get pair ids and order based on fwd file
            rec_cnt = 0
            pair_pos = 0
            self.log (console, "GETTING FWD IDS")  # DEBUG
            with open (input_fwd_file_path, 'r', 0) as input_reads_file_handle:
                rec_line_i = -1
                for line in input_reads_file_handle:
                    rec_line_i += 1
                    if rec_line_i == 3:
                        rec_line_i = -1
                    elif rec_line_i == 0:
                        if not line.startswith('@'):
                            raise ValueError ("badly formatted rec line: '"+line+"'")
                        rec_cnt += 1
                        read_id = line.rstrip('\n')
                        read_id = re.sub ("[ \t]+.*$", "", read_id)
                        read_id = re.sub ("[\/\.\_\-\:\;][012lrLRfrFR53]\'*$", "", read_id)
                        try:
                            pair_ids[read_id] = rev_ids[read_id]
                            pair_pos += 1
                            fwd_id_pos[read_id] = pair_pos
                            pair_ids_order.append(read_id)
                            paired_read_cnt += 1
                        except:
                            unpaired_fwd_read_cnt += 1
                        # DEBUG
#                        if rec_cnt % 100 == 0:
#                            self.log(console,"read_id: '"+str(read_id)+"'")
            total_fwd_recs = rec_cnt
            unpaired_rev_read_cnt = total_rev_recs - paired_read_cnt

            # get rev pos
            for read_id in rev_ids_order:
                try:
                    pair = pair_ids[read_id]
                    pair_pos += 1
                    rev_id_pos[read_id] = pair_pos
                except:
                    pass


            # don't bother if there are no pairs
            if paired_read_cnt == 0:
                raise ValueError ("No pairs found in read library "+readsSet_name_list[lib_i]+" ("+readsSet_ref_list[lib_i]+")")

            # determine if pairs are already in order, or if they're too shuffled to fit in memory
            ordering_offset_upper_bound = 1000000   # only allow a million recs in buf
            ordering_offset_cnt = 0
            last_rev_pos = None
            last_fwd_pos = None
            # THIS LOGIC IS BAD
            for i,read_id in enumerate(pair_ids_order):
                fwd_pos = i+1
                rev_pos = rev_id_pos[read_id]
                if rev_pos > fwd_pos:
                    if last_fwd_pos != None:
                        ordering_offset_cnt -= last_rev_pos - last_fwd_pos

                    ordering_offset_cnt += rev_pos - fwd_pos

                    last_fwd_pos = fwd_pos
                    last_rev_pos = rev_pos

                    if i % 1000 == 0:
                        print (str(read_id)+"\t"+str(fwd_pos)+"\t"+str(rev_pos)+"\t"+str(rev_pos-fwd_pos)+"\t"+str(ordering_offset_cnt))
            # RESTORE when corrected
#            if ordering_offset_cnt > ordering_offset_upper_bound:
#                raise ValueError ("Too many shuffled pairs with too great a distance to fit in memory.  Ordering_offset_cnt="+str(ordering_offset_cnt)+" > Ordering_offset_upper_bound="+str(ordering_offset_upper_bound)+"\nPerhaps do a sort by record ID first")
                
            # determine if there's nothing to do
            if ordering_offset_cnt == 0 and unpaired_fwd_read_cnt == 0 and unpaired_rev_read_cnt == 0:
                self.log (console,"Read Libraries are already Paired and Synchronous")
                continue


            # write fwd paired and fwd unpaired
            #
            self.log (console, "WRITING FWD PAIRED and FWD UNPAIRED")  # DEBUG
            paired_output_reads_file_handle = open (output_fwd_paired_file_path_base+"-"+str(lib_i)+".fastq", 'w', paired_buf_size)
            unpaired_output_reads_file_handle = open (output_fwd_unpaired_file_path_base+"-"+str(lib_i)+".fastq", 'w', unpaired_buf_size)

            rec_buf = []
            unpaired_fwd_buf = []
            last_read_id = None
            paired_cnt = 0
            unpaired_cnt = 0
            capture_type_paired = False

            with open (input_fwd_file_path, 'r', 0) as input_reads_file_handle:
                rec_line_i = -1
                for line in input_reads_file_handle:
                    rec_line_i += 1
                    if rec_line_i == 3:
                        rec_line_i = -1
                    elif rec_line_i == 0:
                        if not line.startswith('@'):
                            raise ValueError ("badly formatted rec line: '"+line+"'")
                        if last_read_id != None:
                            if capture_type_paired:
                                paired_output_reads_file_handle.writelines(rec_buf)
                                paired_cnt += 1
                                #if paired_cnt != 0 and paired_cnt % recs_beep_n == 0:
                                #    self.log(console,"\t"+str(paired_cnt)+" recs processed")
                            else:
                                unpaired_output_reads_file_handle.writelines(rec_buf)
                                unpaired_cnt += 1
                            rec_buf = []
                        read_id = line.rstrip('\n')
                        read_id = re.sub ("[ \t]+.*$", "", read_id)
                        read_id = re.sub ("[\/\.\_\-\:\;][012lrLRfrFR53]\'*$", "", read_id)
                        last_read_id = read_id
                        try:
                            found = pair_ids[read_id]
                            capture_type_paired = True
                        except:
                            capture_type_paired = False
                    rec_buf.append(line)
                # last rec
                if len(rec_buf) > 0:
                    if capture_type_paired:
                        paired_output_reads_file_handle.writelines(rec_buf)
                        paired_cnt += 1
                        #if paired_cnt != 0 and paired_cnt % recs_beep_n == 0:
                        #    self.log(console,"\t"+str(paired_cnt)+" recs processed")
                    else:
                        unpaired_output_reads_file_handle.writelines(rec_buf)
                        unpaired_cnt += 1
                    rec_buf = []

            paired_output_reads_file_handle.close()
            unpaired_output_reads_file_handle.close()
            self.log(console,"\t"+str(paired_cnt)+" PAIRED READS processed")
            self.log(console,"\t"+str(unpaired_cnt)+" UNPAIRED FWD READS processed")
            os.remove (input_fwd_file_path)  # create space since we no longer need the input file
            if paired_cnt != paired_read_cnt:
                raise ValueError ("FAILURE: didn't find expected paired reads in fwd file for lib_i: "+str(lib_i))
            if unpaired_cnt != unpaired_fwd_read_cnt:
                raise ValueError ("FAILURE: didn't find expected unpaired reads in fwd file for lib_i: "+str(lib_i))


            # write rev paired (in order of fwd paired) and rev unpaired.  Store offset reads in memory buffer until correct turn
            #
            self.log (console, "WRITING REV PAIRED and REV UNPAIRED")  # DEBUG
            self.log (console, "USING ORDER FROM FWD PAIRED")  # DEBUG
            paired_output_reads_file_handle = open (output_rev_paired_file_path_base+"-"+str(lib_i)+".fastq", 'w', paired_buf_size)
            unpaired_output_reads_file_handle = open (output_rev_unpaired_file_path_base+"-"+str(lib_i)+".fastq", 'w', unpaired_buf_size)

            rec_buf = []
            paired_unsynch_bufs = dict()
            unpaired_rev_buf = []
            last_read_id = None
            pair_i = 0
            paired_cnt = 0
            unpaired_cnt = 0
            capture_type_paired = False

            with open (input_rev_file_path, 'r', 0) as input_reads_file_handle:
                rec_line_i = -1
                for line in input_reads_file_handle:
                    rec_line_i += 1
                    if rec_line_i == 3:
                        rec_line_i = -1
                    elif rec_line_i == 0:
                        if not line.startswith('@'):
                            raise ValueError ("badly formatted rec line: '"+line+"'")
                        if last_read_id != None:
                            if capture_type_paired:
                                paired_cnt += 1
                                if last_read_id == pair_ids_order[pair_i]:
                                    paired_output_reads_file_handle.writelines(rec_buf)
                                    pair_i += 1
                                    # clear what's available in unsynch buf
                                    while True:
                                        try:
                                            pair_id = pair_ids_order[pair_i]
                                            rec_buf = paired_unsynch_buf[pair_id]
                                            if rec_buf != None:
                                                paired_output_reads_file_handle.writelines(rec_buf)
                                                paired_unsynch_buf[pair_id] = None
                                                pair_i += 1
                                            else:
                                                break
                                        except:
                                            break
                                else:
                                    paired_unsynch_buf[last_read_id] = rec_buf

                                #if paired_cnt != 0 and paired_cnt % recs_beep_n == 0:
                                #    self.log(console,"\t"+str(paired_cnt)+" recs processed")
                            else:
                                unpaired_cnt += 1
                                unpaired_output_reads_file_handle.writelines(rec_buf)
                            rec_buf = []
                        read_id = line.rstrip('\n')
                        read_id = re.sub ("[ \t]+.*$", "", read_id)
                        read_id = re.sub ("[\/\.\_\-\:\;][012lrLRfrFR53]\'*$", "", read_id)
                        last_read_id = read_id
                        try:
                            found = pair_ids[read_id]
                            capture_type_paired = True
                        except:
                            capture_type_paired = False
                    rec_buf.append(line)
                # last rec
                if len(rec_buf) > 0:
                    if last_read_id != None:
                        if capture_type_paired:
                            paired_cnt += 1
                            if last_read_id == pair_ids_order[pair_i]:
                                paired_output_reads_file_handle.writelines(rec_buf)
                                pair_i += 1
                                # clear what's available in unsynch buf
                                while True:
                                    try:
                                        pair_id = pair_ids_order[pair_i]
                                        rec_buf = paired_unsynch_buf[pair_id]
                                        if rec_buf != None:
                                            paired_output_reads_file_handle.writelines(rec_buf)
                                            paired_unsynch_buf[pair_id] = None
                                            pair_i += 1
                                        else:
                                            break
                                    except:
                                        break
                            else:
                                paired_unsynch_buf[last_read_id] = rec_buf

                            #if paired_cnt != 0 and paired_cnt % recs_beep_n == 0:
                            #    self.log(console,"\t"+str(paired_cnt)+" recs processed")
                        else:
                            unpaired_cnt += 1
                            unpaired_output_reads_file_handle.writelines(rec_buf)
                        rec_buf = []

            paired_output_reads_file_handle.close()
            unpaired_output_reads_file_handle.close()
            self.log(console,"\t"+str(paired_cnt)+" PAIRED READS processed")
            self.log(console,"\t"+str(unpaired_rev_read_cnt)+" UNPAIRED REV READS processed")
            os.remove (input_rev_file_path)  # create space since we no longer need the piece file
            if paired_cnt != paired_read_cnt:
                raise ValueError ("FAILURE: didn't find expected paired reads in rev file for lib_i: "+str(lib_i))
            if unpaired_cnt != unpaired_rev_read_cnt:
                raise ValueError ("FAILURE: didn't find expected unpaired reads in rev file for lib_i: "+str(lib_i))


            # add to report
            #
            report += "PAIRED READS: "+str(paired_read_cnt)+"\n"
            report += "UNPAIRED FWD READS: "+str(unpaired_fwd_read_cnt)+"\n"
            report += "UNPAIRED REV READS: "+str(unpaired_rev_read_cnt)+"\n"
            report += "\n"

            if paired_read_cnt == 0:
                raise ValueError ("There were no paired reads")


            # upload paired reads
            #
            if paired_read_cnt > 0:
                self.log (console, "UPLOAD PAIRED READS LIBS")  # DEBUG
                #paired_obj_refs = []
                output_fwd_paired_file_path = output_fwd_paired_file_path_base+"-"+str(lib_i)+".fastq"
                output_rev_paired_file_path = output_rev_paired_file_path_base+"-"+str(lib_i)+".fastq"
                if not os.path.isfile (output_fwd_paired_file_path) \
                        or os.path.getsize (output_fwd_paired_file_path) == 0 \
                      or not os.path.isfile (output_rev_paired_file_path) \
                        or os.path.getsize (output_rev_paired_file_path) == 0:
                    
                    raise ValueError ("failed to create paired output")
                else:
                    output_obj_name = params['output_name']+'_paired-'+str(lib_i)
                    self.log(console, 'Uploading paired reads: '+output_obj_name)
                    paired_obj_refs.append (readsUtils_Client.upload_reads ({ 'wsname': str(params['workspace_name']),
                                                                              'name': output_obj_name,
                                                                              # remove sequencing_tech when source_reads_ref is working
                                                                              #'sequencing_tech': sequencing_tech,
                                                                              'source_reads_ref': input_reads_ref,
                                                                              'fwd_file': output_fwd_paired_file_path,
                                                                              'rev_file': output_rev_paired_file_path
                                                                              })['obj_ref'])
                    

            # upload reads forward unpaired
            #unpaired_fwd_obj_refs = []
            if unpaired_fwd_read_cnt > 0:
                self.log (console, "UPLOAD UNPAIRED FWD READS LIB")  # DEBUG
                unpaired_fwd_ref = None
                output_fwd_unpaired_file_path = output_fwd_unpaired_file_path_base+"-"+str(lib_i)+".fastq"
                if os.path.isfile (output_fwd_unpaired_file_path) \
                        and os.path.getsize (output_fwd_unpaired_file_path) != 0:
                    
                    output_obj_name = params['output_name']+'_unpaired-fwd'
                    self.log(console, '\nUploading trimmed unpaired forward reads: '+output_obj_name)
                    unpaired_fwd_obj_refs.append (readsUtils_Client.upload_reads ({ 'wsname': str(params['workspace_name']),
                                                                                    'name': output_obj_name,
                                                                                    # remove sequencing_tech when source_reads_ref is working
                                                                                    #'sequencing_tech': sequencing_tech,
                                                                                    'source_reads_ref': input_reads_ref,
                                                                                    'fwd_file': output_fwd_unpaired_file_path
                                                                                    })['obj_ref'])
                else:
                    unpaired_fwd_obj_refs.append (None)


            # upload reads reverse unpaired
            #unpaired_rev_obj_refs = []
            if unpaired_rev_read_cnt > 0:
                self.log (console, "UPLOAD UNPAIRED REV READS LIB")  # DEBUG
                unpaired_rev_ref = None
                output_rev_unpaired_file_path = output_rev_unpaired_file_path_base+"-"+str(lib_i)+".fastq"
                if os.path.isfile (output_rev_unpaired_file_path) \
                        and os.path.getsize (output_rev_unpaired_file_path) != 0:

                    output_obj_name = params['output_name']+'_unpaired-rev'
                    self.log(console, '\nUploading trimmed unpaired reverse reads: '+output_obj_name)
                    unpaired_rev_obj_refs.append (readsUtils_Client.upload_reads ({ 'wsname': str(params['workspace_name']),
                                                                                    'name': output_obj_name,
                                                                                    # remove sequencing_tech when source_reads_ref is working
                                                                                    #'sequencing_tech': sequencing_tech,
                                                                                    'source_reads_ref': input_reads_ref,
                                                                                    'fwd_file': output_rev_unpaired_file_path
                                                                                    })['obj_ref'])
                else:
                    unpaired_rev_obj_refs.append (None)

        
        # Create ReadsSets for paired libs and unpaired fwd and rev libs if input was ReadsSet
        #
        if input_reads_obj_type == "KBaseSets.ReadsSet":

            paired_readsSet_ref = None
            unpaired_fwd_readsSet_ref = None
            unpaired_rev_readsSet_ref = None

            # save paired readsSet
            some_paired_output_created = False
            items = []
            for i,lib_ref in enumerate(paired_obj_refs):   # FIX: assumes order maintained
                if lib_ref == None:
                    #items.append(None)  # can't have 'None' items in ReadsSet
                    continue
                else:
                    some_paired_output_created = True
                    try:
                        label = input_readsSet_obj['data']['items'][i]['label']
                    except:
                        NAME_I = 1
                        label = wsClient.get_object_info_new ({'objects':[{'ref':lib_ref}]})[0][NAME_I]
                    label = label + "_paired_synched"

                    items.append({'ref': lib_ref,
                                  'label': label
                                  #'data_attachment': ,
                                  #'info':
                                      })
            if some_paired_output_created:
                reads_desc_ext = " Synched paired reads"
                reads_name_ext = "_paired_synched"
                output_readsSet_obj = { 'description': input_readsSet_obj['data']['description']+reads_desc_ext,
                                        'items': items
                                        }
                output_readsSet_name = str(params['output_name'])+reads_name_ext
                paired_readsSet_ref = setAPI_Client.save_reads_set_v1 ({'workspace_name': params['workspace_name'],
                                                                        'output_object_name': output_readsSet_name,
                                                                        'data': output_readsSet_obj
                                                                        })['set_ref']
            else:
                raise ValueError ("No paired output created")

                              
            # save unpaired forward readsSet
            some_unpaired_fwd_output_created = False
            if len(unpaired_fwd_obj_refs) > 0:
                items = []
                for i,lib_ref in enumerate(unpaired_fwd_obj_refs):  # FIX: assumes order maintained
                    if lib_ref == None:
                        #items.append(None)  # can't have 'None' items in ReadsSet
                        continue
                    else:
                        some_unpaired_fwd_output_created = True
                        try:
                            if len(unpaired_fwd_readsSet_refs) == len(input_readsSet_obj['data']['items']):
                                label = input_readsSet_obj['data']['items'][i]['label']
                            else:
                                NAME_I = 1
                                label = wsClient.get_object_info_new ({'objects':[{'ref':lib_ref}]})[0][NAME_I]
                        except:
                            NAME_I = 1
                            label = wsClient.get_object_info_new ({'objects':[{'ref':lib_ref}]})[0][NAME_I]
                        label = label + "_unpaired_fwd"

                        items.append({'ref': lib_ref,
                                      'label': label
                                      #'data_attachment': ,
                                      #'info':
                                          })
                if some_unpaired_fwd_output_created:
                    reads_desc_ext = " Unpaired FWD reads"
                    reads_name_ext = "_unpaired_fwd"
                    output_readsSet_obj = { 'description': input_readsSet_obj['data']['description']+reads_desc_ext,
                                            'items': items
                                            }
                    output_readsSet_name = str(params['output_name'])+reads_name_ext
                    unpaired_fwd_readsSet_ref = setAPI_Client.save_reads_set_v1 ({'workspace_name': params['workspace_name'],
                                                                                  'output_object_name': output_readsSet_name,
                                                                                  'data': output_readsSet_obj
                                                                                  })['set_ref']
                else:
                    self.log (console, "no unpaired_fwd readsLibraries created")
                    unpaired_fwd_readsSet_ref = None
                              
            # save unpaired reverse readsSet
            some_unpaired_rev_output_created = False
            if len(unpaired_rev_obj_refs) > 0:
                items = []
                for i,lib_ref in enumerate(unpaired_fwd_obj_refs):  # FIX: assumes order maintained
                    if lib_ref == None:
                        #items.append(None)  # can't have 'None' items in ReadsSet
                        continue
                    else:
                        some_unpaired_rev_output_created = True
                        try:
                            if len(unpaired_rev_readsSet_refs) == len(input_readsSet_obj['data']['items']):
                                label = input_readsSet_obj['data']['items'][i]['label']
                            else:
                                NAME_I = 1
                                label = wsClient.get_object_info_new ({'objects':[{'ref':lib_ref}]})[0][NAME_I]

                        except:
                            NAME_I = 1
                            label = wsClient.get_object_info_new ({'objects':[{'ref':lib_ref}]})[0][NAME_I]
                        label = label + "_unpaired_rev"

                        items.append({'ref': lib_ref,
                                      'label': label
                                      #'data_attachment': ,
                                      #'info':
                                          })
                if some_unpaired_rev_output_created:
                    reads_desc_ext = " Unpaired REV reads"
                    reads_name_ext = "_unpaired_rev"
                    output_readsSet_obj = { 'description': input_readsSet_obj['data']['description']+reads_desc_ext,
                                            'items': items
                                            }
                    output_readsSet_name = str(params['output_name'])+reads_name_ext
                    unpaired_rev_readsSet_ref = setAPI_Client.save_reads_set_v1 ({'workspace_name': params['workspace_name'],
                                                                                  'output_object_name': output_readsSet_name,
                                                                                  'data': output_readsSet_obj
                                                                                  })['set_ref']
                else:
                    self.log (console, "no unpaired_rev readsLibraries created")
                    unpaired_rev_readsSet_ref = None



        # build report
        #
        self.log (console, "SAVING REPORT")  # DEBUG        
        reportObj = {'objects_created':[], 
                     'text_message': report}

        if input_reads_obj_type == "KBaseSets.ReadsSet":
            if paired_readsSet_ref != None:
                reportObj['objects_created'].append({'ref':paired_readsSet_ref,
                                                     'description':params['desc']+": PAIRED"})
            if unpaired_fwd_readsSet_ref != None:
                reportObj['objects_created'].append({'ref':unpaired_fwd_readsSet_ref,
                                                     'description':params['desc']+": UNPAIRED FWD"})
            if unpaired_rev_readsSet_ref != None:
                reportObj['objects_created'].append({'ref':unpaired_rev_readsSet_ref,
                                                     'description':params['desc']+": UNPAIRED REV"})
        else:  # Single Library
            if len(paired_obj_refs) > 0:
                reportObj['objects_created'].append({'ref':paired_obj_refs[0],
                                                     'description':params['desc']+": PAIRED"})
            if len(unpaired_fwd_obj_refs) > 0:
                reportObj['objects_created'].append({'ref':unpaired_fwd_obj_refs[0],
                                                     'description':params['desc']+": UNPAIRED FWD"})
            if len(unpaired_rev_obj_refs) > 0:
                reportObj['objects_created'].append({'ref':unpaired_rev_obj_refs[0],
                                                     'description':params['desc']+": UNPAIRED REV"})


        # save report object
        #
        report = KBaseReport(self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
        report_info = report.create({'report':reportObj, 'workspace_name':params['workspace_name']})

        returnVal = { 'report_name': report_info['name'], 'report_ref': report_info['ref'] }
        #END KButil_Extract_Unpaired_Reads_and_Synchronize_Pairs

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method KButil_Extract_Unpaired_Reads_and_Synchronize_Pairs return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def KButil_Translate_ReadsLibs_QualScores(self, ctx, params):
        """
        :param params: instance of type
           "KButil_Translate_ReadsLibs_QualScores_Params"
           (KButil_Translate_ReadsLibs_QualScores() ** **  Method for
           Translating ReadsLibs Qual scores) -> structure: parameter
           "workspace_name" of type "workspace_name" (** The workspace object
           refs are of form: ** **    objects = ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_refs" of type "data_obj_ref"
        :returns: instance of type
           "KButil_Translate_ReadsLibs_QualScores_Output" -> structure:
           parameter "report_name" of type "data_obj_name", parameter
           "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN KButil_Translate_ReadsLibs_QualScores
        console = []
        invalid_msgs = []
        report = ''
        self.log(console, 'Running KButil_Translate_ReadsLibs_QualScores with parameters: ')
        self.log(console, "\n"+pformat(params))

        token = ctx['token']
        wsClient = workspaceService(self.workspaceURL, token=token)
        headers = {'Authorization': 'OAuth '+token}
        env = os.environ.copy()
        env['KB_AUTH_TOKEN'] = token

        # internal Methods
        def qual33(qual64): return chr(ord(qual64)-31)

        #SERVICE_VER = 'dev'  # DEBUG
        SERVICE_VER = 'release'

        # param checks
        required_params = ['workspace_name',
                           'input_refs', 
                           ]
        for required_param in required_params:
            if required_param not in params or params[required_param] == None:
                raise ValueError ("Must define required param: '"+required_param+"'")
            
        # clean input_refs
        clean_input_refs = []
        for ref in params['input_refs']:
            if ref != None and ref != '':
                clean_input_refs.append(ref)
        params['input_refs'] = clean_input_refs


        # load provenance
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
        for input_ref in params['input_refs']:
            provenance[0]['input_ws_objects'].append(input_ref)

        # Determine whether read library or read set is input object
        #
        first_input_ref = params['input_refs']


        # get set
        #
        readsSet_ref_list = []
        readsSet_names_list = []
        for reads_ref in params['input_refs']:
            try:
                # object_info tuple
                [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)
                
                input_reads_obj_info = wsClient.get_object_info_new ({'objects':[{'ref':reads_ref}]})[0]
                input_reads_obj_type = input_reads_obj_info[TYPE_I]
                input_reads_obj_type = re.sub ('-[0-9]+\.[0-9]+$', "", input_reads_obj_type)  # remove trailing version
                #input_reads_obj_version = input_reads_obj_info[VERSION_I]  # this is object version, not type version

            except Exception as e:
                raise ValueError('Unable to get read library object from workspace: (' + str(reads_ref) +')' + str(e))
            
            acceptable_types = ["KBaseSets.ReadsSet", "KBaseFile.PairedEndLibrary"]
            if input_reads_obj_type not in acceptable_types:
                raise ValueError ("Input reads of type: '"+input_reads_obj_type+"'.  Must be one of "+", ".join(acceptable_types))
            
            if input_reads_obj_type != "KBaseSets.ReadsSet":  # readsLib
                readsSet_ref_list.append(reads_ref)

            else:  # readsSet
                try:
                    setAPI_Client = SetAPI (url=self.serviceWizardURL, token=ctx['token'])  # for dynamic service
                    input_readsSet_obj = setAPI_Client.get_reads_set_v1 ({'ref':reads_ref,'include_item_info':1})
                except Exception as e:
                    raise ValueError('SetAPI FAILURE: Unable to get read library set object from workspace: (' + str(reads_ref)+")\n" + str(e))
                
                for readsLibrary_obj in input_readsSet_obj['data']['items']:
                    readsSet_ref_list.append(readsLibrary_obj['ref'])
#                    NAME_I = 1
#                    readsSet_names_list.append(readsLibrary_obj['info'][NAME_I])
        # add names and types
        reads_obj_types_list = []
        for reads_ref in readsSet_ref_list:
            try:
                # object_info tuple
                [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)
                
                input_reads_obj_info = wsClient.get_object_info_new ({'objects':[{'ref':reads_ref}]})[0]
                input_reads_obj_name = input_reads_obj_info[NAME_I]
                input_readsLib_obj_type = input_reads_obj_info[TYPE_I]
                input_readsLib_obj_type = re.sub ('-[0-9]+\.[0-9]+$', "", input_reads_obj_type)  # remove trailing version

            except Exception as e:
                raise ValueError('Unable to get read library object from workspace: (' + str(reads_ref) +')' + str(e))

            readsSet_names_list.append (input_reads_obj_name)
            reads_obj_types_list.append (input_readsLib_obj_type)


        # translate qual scores for each read library
        #
        self.log (console, "CREATING Translated FASTQ FILES")

        # make dir
        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()*1000)
        input_dir = os.path.join(self.scratch,'input.'+str(timestamp))
        if not os.path.exists(input_dir):
            os.makedirs(input_dir)

        # connect to ReadsUtils Client
        try:
            readsUtils_Client = ReadsUtils (url=self.callbackURL, token=ctx['token'])  # SDK local
        except:
            raise ValueError("Unable to get readsUtils_Client\n" + str(e))


        # add libraries, one at a time
        #
        new_objects = []
        translated_cnt = 0
        for reads_i,this_input_reads_ref in enumerate(readsSet_ref_list):
            self.log (console, "DOWNLOADING FASTQ FILES FOR ReadsSet member: "+readsSet_names_list[reads_i]+" ("+str(this_input_reads_ref)+")")
            try:
                readsLibrary = readsUtils_Client.download_reads ({'read_libraries': [this_input_reads_ref],
                                                                  'interleaved': 'false'
                                                                  })
            except Exception as e:
                raise ValueError('Unable to get reads object from workspace: (' + this_input_reads_ref +")\n" + str(e))

            this_input_fwd_path = readsLibrary['files'][this_input_reads_ref]['files']['fwd']
            if reads_obj_types_list[reads_i] == "KBaseFile.PairedEndLibrary":
                this_input_rev_path = readsLibrary['files'][this_input_reads_ref]['files']['rev']

            # read through and translate qual scores
            self.log (console, "TRANSLATING FWD FASTQ FILE FOR ReadsSet member: "+str(this_input_reads_ref))
            read_buf_size  = 65536
            write_buf_size = 65536

            qual33_fwd_path = this_input_fwd_path + '.qual33'
            qual33_fwd_handle = open (qual33_fwd_path, 'w', write_buf_size)

            input_is_already_phred33 = False
            this_input_path = this_input_fwd_path
            with open (this_input_path, 'r', read_buf_size) as this_input_handle:
                while True:
                    buf = []
                    line = this_input_handle.readline()
                    if not line:
                        break
                    if input_is_already_phred33:
                        break
                    if line.startswith('@'):
                        buf.append(line)  # header
                        buf.append(this_input_handle.readline())  # seq
                        buf.append(this_input_handle.readline())  # '+'

                        qual_line = this_input_handle.readline().rstrip()
                        q33_line = ''
                        #def qual33(qual64): return chr(ord(qual64)-31)
                        #trans_report = ''  # DEBUG
                        #self.log (console, "ORIG_LINE: "+qual_line)  # DEBUG
                        for q64 in qual_line:
                            q64_ascii = ord(q64)
                            #trans_report += q64+'('+str(q64_ascii)+')'
                            if q64_ascii < 64:
                                input_is_already_phred33 = True
                                break
                            q33_line += chr(q64_ascii - 31)
                        buf.append(q33_line+"\n")
                        #self.log (console, "TRNS_LINE: "+trans_report)  # DEBUG
                        #self.log (console, "TRNS_LINE: "+q33_line)  # DEBUG
                        qual33_fwd_handle.write(''.join(buf))

            qual33_fwd_handle.close()
            os.remove (this_input_path)  # create space since we no longer need the piece file

            # append rev
            if reads_obj_types_list[reads_i] == "KBaseFile.PairedEndLibrary" and \
                    not input_is_already_phred33:
                
                self.log (console, "TRANSLATING REV FASTQ FILE FOR ReadsSet member: "+str(this_input_reads_ref))

                qual33_rev_path = this_input_rev_path + '.qual33'
                qual33_rev_handle = open (qual33_rev_path, 'w', write_buf_size)
                this_input_path = this_input_rev_path

                with open (this_input_path, 'r', read_buf_size) as this_input_handle:
                    while True:
                        buf = []
                        line = this_input_handle.readline()
                        if not line:
                            break
                        if input_is_already_phred33:
                            break
                        if line.startswith('@'):
                            buf.append(line)  # header
                            buf.append(this_input_handle.readline())  # seq
                            buf.append(this_input_handle.readline())  # '+'
                            
                            qual_line = this_input_handle.readline().rstrip()
                            q33_line = ''
                            #self.log (console, "ORIG_LINE: "+qual_line)  # DEBUG
                            for q64 in qual_line:
                                q64_ascii = ord(q64)
                                if q64_ascii < 64:
                                    input_is_already_phred33 = True
                                    break
                                q33_line += chr(q64_ascii - 31)
                            buf.append(q33_line+"\n")
                            #self.log (console, "TRNS_LINE: "+q33_line)  # DEBUG
                            qual33_rev_handle.write(''.join(buf))

                qual33_rev_handle.close()
                os.remove (this_input_path)  # create space since we no longer need the piece file

            # upload reads
            #
            if input_is_already_phred33:
                self.log(console, "WARNING: "+readsSet_names_list[reads_i]+" ("+str(this_input_reads_ref)+") is already phred33.  Skipping.")
                continue

            translated_cnt += 1
            self.log (console, "UPLOADING Translated 64->33 READS LIB")  # DEBUG
            if not os.path.isfile (qual33_fwd_path) \
                    or os.path.getsize (qual33_fwd_path) == 0:
                raise ValueError ("failed to create fwd read library output")
            if reads_obj_types_list[reads_i] == "KBaseFile.PairedEndLibrary":
                if not os.path.isfile (qual33_rev_path) \
                        or os.path.getsize (qual33_rev_path) == 0:
                    
                    raise ValueError ("failed to create rev read library output")

            output_obj_name = readsSet_names_list[reads_i]+".phred33"
            self.log(console, 'Uploading reads library: '+output_obj_name)

            if reads_obj_types_list[reads_i] == "KBaseFile.PairedEndLibrary":
                reads_library_ref = readsUtils_Client.upload_reads ({ 'wsname': str(params['workspace_name']),
                                                                      'name': output_obj_name,
                                                                      # remove sequencing_tech when source_reads_ref is working
                                                                      #'sequencing_tech': sequencing_tech,
                                                                      'source_reads_ref': readsSet_ref_list[0],
                                                                      'fwd_file': qual33_fwd_path,
                                                                      'rev_file': qual33_rev_path
                                                                      })['obj_ref']
            else:
                reads_library_ref = readsUtils_Client.upload_reads ({ 'wsname': str(params['workspace_name']),
                                                                      'name': output_obj_name,
                                                                      # remove sequencing_tech when source_reads_ref is working
                                                                      #'sequencing_tech': sequencing_tech,
                                                                      'source_reads_ref': readsSet_ref_list[0],
                                                                      'fwd_file': qual33_fwd_path,
                                                                      })['obj_ref']
            

            # add object to list
            desc = readsSet_names_list[reads_i]+' translated to phred33'
            new_objects.append({'ref':reads_library_ref,
                                'description':desc})


        # build report message
        report += "NUM READS LIBRARIES INPUT: " + str(len(readsSet_ref_list))+"\n"
        report += "NUM READS LIBRARIES TRANSLATED: " + str(translated_cnt)+"\n"


        # build report
        #
        self.log (console, "SAVING REPORT")  # DEBUG        
        reportObj = {'objects_created':new_objects, 
                     'text_message': report}


        # save report object
        #
        report = KBaseReport(self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
        report_info = report.create({'report':reportObj, 'workspace_name':params['workspace_name']})

        returnVal = { 'report_name': report_info['name'], 'report_ref': report_info['ref'] }
        #END KButil_Translate_ReadsLibs_QualScores

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method KButil_Translate_ReadsLibs_QualScores return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def KButil_AddInsertLen_to_ReadsLibs(self, ctx, params):
        """
        :param params: instance of type
           "KButil_AddInsertLen_to_ReadsLibs_Params"
           (KButil_AddInsertLen_to_ReadsLibs() ** **  Method for Adding
           Insert Len to PairedEnd ReadsLibs) -> structure: parameter
           "workspace_name" of type "workspace_name" (** The workspace object
           refs are of form: ** **    objects = ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_refs" of type "data_obj_ref", parameter
           "insert_len" of Double, parameter "insert_stddev" of Double
        :returns: instance of type "KButil_AddInsertLen_to_ReadsLibs_Output"
           -> structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN KButil_AddInsertLen_to_ReadsLibs
        console = []
        invalid_msgs = []
        report = ''
        self.log(console, 'Running KButil_AddInsertLen_to_ReadsLibs with parameters: ')
        self.log(console, "\n"+pformat(params))

        token = ctx['token']
        wsClient = workspaceService(self.workspaceURL, token=token)
        headers = {'Authorization': 'OAuth '+token}
        env = os.environ.copy()
        env['KB_AUTH_TOKEN'] = token

        #SERVICE_VER = 'dev'  # DEBUG
        SERVICE_VER = 'release'

        # param checks
        required_params = ['workspace_name',
                           'input_refs', 
                           'insert_len',
                           'insert_stddev'
                           ]
        for required_param in required_params:
            if required_param not in params or params[required_param] == None:
                raise ValueError ("Must define required param: '"+required_param+"'")
            
        # clean input_refs
        clean_input_refs = []
        for ref in params['input_refs']:
            if ref != None and ref != '':
                clean_input_refs.append(ref)
        params['input_refs'] = clean_input_refs


        # Determine whether read library or read set is input object
        #
        first_input_ref = params['input_refs']


        # get set
        #
        readsSet_ref_list = []
        readsSet_names_list = []
        for reads_ref in params['input_refs']:
            try:
                # object_info tuple
                [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)
                
                input_reads_obj_info = wsClient.get_object_info_new ({'objects':[{'ref':reads_ref}]})[0]
                input_reads_obj_type = input_reads_obj_info[TYPE_I]
                input_reads_obj_type = re.sub ('-[0-9]+\.[0-9]+$', "", input_reads_obj_type)  # remove trailing version
                #input_reads_obj_version = input_reads_obj_info[VERSION_I]  # this is object version, not type version

            except Exception as e:
                raise ValueError('Unable to get read library object from workspace: (' + str(reads_ref) +')' + str(e))
            
            acceptable_types = ["KBaseSets.ReadsSet", "KBaseFile.PairedEndLibrary"]
            if input_reads_obj_type not in acceptable_types:
                raise ValueError ("Input reads of type: '"+input_reads_obj_type+"'.  Must be one of "+", ".join(acceptable_types))
            
            if input_reads_obj_type != "KBaseSets.ReadsSet":  # readsLib
                readsSet_ref_list.append(reads_ref)

            else:  # readsSet
                try:
                    setAPI_Client = SetAPI (url=self.serviceWizardURL, token=ctx['token'])  # for dynamic service
                    input_readsSet_obj = setAPI_Client.get_reads_set_v1 ({'ref':reads_ref,'include_item_info':1})
                except Exception as e:
                    raise ValueError('SetAPI FAILURE: Unable to get read library set object from workspace: (' + str(reads_ref)+")\n" + str(e))
                
                for readsLibrary_obj in input_readsSet_obj['data']['items']:
                    readsSet_ref_list.append(readsLibrary_obj['ref'])
#                    NAME_I = 1
#                    readsSet_names_list.append(readsLibrary_obj['info'][NAME_I])
        # add names and types
        reads_obj_types_list = []
        for reads_ref in readsSet_ref_list:
            try:
                # object_info tuple
                [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)
                
                input_reads_obj_info = wsClient.get_object_info_new ({'objects':[{'ref':reads_ref}]})[0]
                input_reads_obj_name = input_reads_obj_info[NAME_I]
                input_readsLib_obj_type = input_reads_obj_info[TYPE_I]
                input_readsLib_obj_type = re.sub ('-[0-9]+\.[0-9]+$', "", input_reads_obj_type)  # remove trailing version

            except Exception as e:
                raise ValueError('Unable to get read library object from workspace: (' + str(reads_ref) +')' + str(e))

            readsSet_names_list.append (input_reads_obj_name)
            reads_obj_types_list.append (input_readsLib_obj_type)


        # Add insert len vals to read library objects
        #
        self.log (console, "CREATING UPDATED READ LIBRARY OBJECTS")

        new_objects = []
        updated_cnt = 0
        for reads_i,this_input_reads_ref in enumerate(readsSet_ref_list):
            self.log (console, "UPDATING ReadsSet member: "+readsSet_names_list[reads_i]+" ("+str(this_input_reads_ref)+")")

            try:
                objects = wsClient.get_objects2({'objects':[{'ref': this_input_reads_ref}]})['data']
                new_reads_lib_obj = objects[0]['data']
                #info = objects[0]['info']
            except Exception as e:
                raise ValueError('Unable to get reads object from workspace: '+str(readsSet_names_list[reads_i])+'(' + this_input_reads_ref +")\n" + str(e))

            new_reads_lib_obj['insert_size_mean'] = float(params['insert_len'])
            new_reads_lib_obj['insert_size_std_dev'] = float(params['insert_stddev'])

            updated_cnt += 1

            # Save updated object
            self.log (console, "SAVING READS LIB")  # DEBUG
            output_obj_name = readsSet_names_list[reads_i]

            # load provenance
            provenance = [{}]
            if 'provenance' in ctx:
                provenance = ctx['provenance']
            provenance[0]['input_ws_objects'] = [this_input_reads_ref]

            try:
                new_obj_info = wsClient.save_objects({
                        'workspace':params['workspace_name'],
                        'objects':[
                            {
                                'type':'KBaseFile.PairedEndLibrary',
                                'data':new_reads_lib_obj,
                                'name':output_obj_name,
                                'meta':{},
                                'provenance':provenance
                            }]
                        })[0]
            except Exception as e:
                raise ValueError('Unable to save reads object to workspace: '+str(readsSet_names_list[reads_i])+'(' + this_input_reads_ref +")\n" + str(e))

            # object_info tuple
            [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)
            # new ref
            new_reads_library_ref = str(new_obj_info[WSID_I]) +'/'+str(new_obj_info[OBJID_I]) +'/'+ str(new_obj_info[VERSION_I])

            # add object to list
            desc = readsSet_names_list[reads_i]+' (with Insert Length and Stddev Added)'
            new_objects.append({'ref':new_reads_library_ref,
                                'description':desc})

        # Create new ReadsSet with updated libs if more than one input lib
        if len(params['input_refs']) > 1:
            readsSet_desc = "ReadsLibs "+", ".join(readsSet_names_list)+" with InsertLen:"+str(params['insert_len'])+" Insert_STDDEV:"+str(params['insert_stddev'])
            items = []
            for reads_i,this_input_reads_ref in enumerate(readsSet_ref_list):
                items.append({ 'ref': new_objects[reads_i]['ref'],
                               'label': readsSet_names_list[reads_i]
                           })

            try:
                setAPI_Client = SetAPI (url=self.serviceWizardURL, token=ctx['token'])  # for dynamic service
            except Exception as e:
                raise ValueError('ERROR: unable to instantiate SetAPI' + str(e))

            output_readsSet_obj = { 'description': readsSet_desc,
                                    'items': items
                                    }
            output_readsSet_name = 'Updated_ReadsLibs_with_Insert_Len.ReadsSet'
            try:
                output_readsSet_ref = setAPI_Client.save_reads_set_v1 ({'workspace_name': params['workspace_name'],
                                                                        'output_object_name': output_readsSet_name,
                                                                        'data': output_readsSet_obj
                                                                        })['set_ref']
            except Exception as e:
                raise ValueError('SetAPI FAILURE: Unable to save read library set object to workspace: (' + params['workspace_name']+")\n" + str(e))
            new_objects = [{'ref':output_readsSet_ref,'description':readsSet_desc}] + new_objects


        # build report message
        report += "NUM READS LIBRARIES INPUT: " + str(len(readsSet_ref_list))+"\n"
        report += "NUM READS LIBRARIES TRANSLATED: " + str(updated_cnt)+"\n"


        # build report
        #
        self.log (console, "SAVING REPORT")  # DEBUG        
        reportObj = {'objects_created':new_objects, 
                     'text_message': report}


        # save report object
        #
        report = KBaseReport(self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
        report_info = report.create({'report':reportObj, 'workspace_name':params['workspace_name']})

        returnVal = { 'report_name': report_info['name'], 'report_ref': report_info['ref'] }
        #END KButil_AddInsertLen_to_ReadsLibs

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method KButil_AddInsertLen_to_ReadsLibs return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
