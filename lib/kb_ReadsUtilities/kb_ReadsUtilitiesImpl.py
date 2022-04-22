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
import math
import gzip
import random

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

from installed_clients.WorkspaceClient import Workspace as workspaceService

# SDK Utils
from installed_clients.ReadsUtilsClient import ReadsUtils
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.SetAPIServiceClient import SetAPI
from installed_clients.KBaseReportClient import KBaseReport

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
    VERSION = "1.2.2"
    GIT_URL = "https://github.com/kbaseapps/kb_ReadsUtilities"
    GIT_COMMIT_HASH = "e3e6239b4d6cf1f5b5cc7178c5447c473ae4198f"

    #BEGIN_CLASS_HEADER
    workspaceURL = None
    shockURL = None
    handleURL = None
    serviceWizardsURL = None
    callbackURL = None
    scratch = None


    def now_ISO(self):
        now_timestamp = datetime.now()
        now_secs_from_epoch = (now_timestamp - datetime(1970,1,1)).total_seconds()
        now_timestamp_in_iso = datetime.fromtimestamp(int(now_secs_from_epoch)).strftime('%Y-%m-%d_%T')
        return now_timestamp_in_iso

    def log(self, target, message):
        message = '['+self.now_ISO()+'] '+message
        if target is not None:
            target.append(message)
        print(message)
        sys.stdout.flush()

    def add_id_to_plus_line(self,readsfile):
        gzip_flag = False
        if readsfile.lower().endswith('.gz'):
            gzip_flag = True

        if gzip_flag:
            f = gzip.open(readsfile, 'rt')
        else:
            f = open(readsfile, 'r')        

        outputfile = readsfile+".with_IDs"
        if gzip_flag:
            outputfile += ".gz"
            out = gzip.open(outputfile, 'wt')
        else:
            out = open(outputfile, 'w')
    
        this_id = None
        counter = 0
        read_cnt = 0
        for line in f:
            line = line.rstrip()
            if counter == 4:
                counter = 0
            if line.startswith('@') and counter == 0:
                this_id = line.lstrip('@')
                read_cnt += 1
                if read_cnt % 1000000 == 0:
                    print ("reads processed {}".format(read_cnt))
            if line.startswith('+') and counter == 2 and len(line) < 3:
                line = '+'+this_id
            out.write(line+"\n")
            counter += 1

        print ("READS processed {}".format(read_cnt))
        f.close()
        out.close()
            
        os.rename (outputfile, readsfile)
        return

    def reverse_complement (self,seq):
        rev_seq = ''
        complement = {
            'A': 'T',
            'T': 'A',
            'G': 'C',
            'C': 'G',
            'N': 'N'
        }
        for c in seq[::-1]:
            if not complement.get(c):
                raise ValueError ("unknown character '"+c+"' in seq '"+seq+"'")
            rev_seq += complement[c]

        return rev_seq

    def seq_call_with_error(self,c=None, qual=None, qual_type='phred33'):
        if c == 'N': return c

        if not c or not qual:
            raise ValueError ("missing param for seq_call_with_error()")

        def qual33(qual64): return chr(ord(qual64)-31)

        def error_prob(qual33):
            # P = 10^-Q/10
            # Q = -10*log10(P)
            Q = ord(qual33)-33
            P = math.pow(10,-Q/10.0)
            return P

        base = ['A', 'T', 'G', 'C']

        if qual_type == 'phred64':
            qual = qual33(qual)
        elif qual_type != 'phred33':
            raise ValueError ("seq_call_with_error() can only handle phred33 and phred64")

        if random.uniform(0.0,1.0) <= error_prob(qual):
            #print("adding error to call.  QUAL: "+str(qual)+" prob: "+str(error_prob(qual))+" orig_base: "+c)  # DEBUG
            c = base[random.randint(0,3)]
            #print("new base: "+c)  # DEBUG
        return c


    # replace read with source genome sequence
    def overlay_source_genome_seq (self,
                                   read_rec=None,
                                   source_genome_buf=None,
                                   contig_mapping=None,
                                   lib_type='PE',  # 'PE/SE'
                                   read_dir='fwd',  # 'fwd/rev'
                                   fwd_insilico_pos=None,  # [contig_i,strand(+/-),beg_pos]
                                   pe_orientation='IN-IN',  # 'IN-IN' only at this time
                                   pe_insert_len=None,  # typically 450: fwd+gap+rev = 150+150+150
                                   add_errors_by_qual_freq=True):

        [POS_CONTIG_I, POS_STRAND_I, POS_BEG_I] = range(3)
        [READ_HEADER_LINE_I, READ_SEQ_LINE_I, READ_SPACER_LINE_I, READ_QUAL_LINE_I] = range(4)
        insilico_pos = []
        insilico_read_rec_buf = []

        # call check
        if not read_rec or not source_genome_buf or not contig_mapping:
            raise ValueError ("missing required args to overlay_source_genome_seq()")
        if lib_type == 'PE':
            if not pe_orientation or not pe_insert_len:
                raise ValueError ("missing pe_orientation or pe_insert_len for overlay_source_genome_seq()")
            if read_dir == 'rev':
                if not fwd_insilico_pos:
                    raise ValueError ("missing fwd_insilico_pos for reverse read in overlay_source_genome_seq()")
        if pe_insert_len:
            pe_insert_len = int(pe_insert_len)

        # break up rec
        header_line = read_rec[READ_HEADER_LINE_I]
        seq_line    = read_rec[READ_SEQ_LINE_I].rstrip()
        qual_line   = read_rec[READ_QUAL_LINE_I].rstrip()
        if len(seq_line) != len(qual_line):
            raise ValueError ("inconsistent lengths for seq and qual in read rec: "+header_line)
        read_len = len(seq_line)

        # create new rec
        insilico_read_rec_buf.append(header_line)
        insilico_read_rec_buf.append('')  # add sequence later
        insilico_read_rec_buf.append(read_rec[READ_SPACER_LINE_I])  # should just be the '+' symbol
        insilico_read_rec_buf.append(qual_line+"\n")

        # if reverse read, assign contig, strand, and pos from fwd_insilico_pos
        if lib_type == 'PE' and read_dir == 'rev':
            contig_i = fwd_insilico_pos[POS_CONTIG_I]
            if pe_orientation == 'IN-IN':
                if fwd_insilico_pos[POS_STRAND_I] == '+':
                    strand = '-'
                    beg_pos = fwd_insilico_pos[POS_BEG_I] + (pe_insert_len - 1)
                else:
                    strand = '+'
                    beg_pos = fwd_insilico_pos[POS_BEG_I] - (pe_insert_len - 1)
            else:
                raise ValueError ("do not yet support Paired-End orientations other than IN-IN in overlay_source_genome_seq()")

        # pick contig
        else:  # read_dir == 'fwd'
            acceptable_contig = False
            max_tries = 1000
            contig_i = -1
            contig_len = -1
            for try_i in range(max_tries):
                contig_i = contig_mapping[random.randint (0,len(contig_mapping)-1)]
                contig_len = len(source_genome_buf[contig_i])
                if lib_type == 'SE':
                    if read_len < contig_len:
                        acceptable_contig = True
                        break
                else:
                    if 2*read_len < contig_len and pe_insert_len < contig_len:
                        acceptable_contig = True
                        break
            if not acceptable_contig:
                raise ValueError ("unable to find long enough contig for "+header_line)

            # pick strand
            if random.randint(0,1) == 0:
                strand = '+'
            else:
                strand = '-'

            # pick beg_pos
            max_tries = 1000
            beg_pos = -1
            acceptable_pos = False
            for try_i in range(max_tries):
                if strand == '+':
                    if lib_type == 'SE':
                        beg_pos = random.randint(0,contig_len - read_len - 1)
                        acceptable_pos = True
                        break
                    else:
                        beg_pos = random.randint(0,contig_len - pe_insert_len - 1)
                        if contig_len - beg_pos > 2*read_len:
                            acceptable_pos = True
                            break
                else:  # strand == '-'
                    if lib_type == 'SE':
                        beg_pos = random.randint(read_len-1, contig_len-1)
                        acceptable_pos = True
                        break
                    else:
                        beg_pos = random.randint(pe_insert_len-1, contig_len-1)
                        if beg_pos > 2*read_len:
                            acceptable_pos = True
                            break
            if not acceptable_pos:
                raise ValueError ("unable to find acceptable position for "+header_line)

        # set insilico_pos to return (only matters for fwd reads of PE lib)
        insilico_pos = [contig_i, strand, beg_pos]

        # overlay genome sequence onto read
        if strand == '+':
            genome_slice = source_genome_buf[contig_i][beg_pos:beg_pos+read_len]
            #print("POSSEQ: "+source_genome_buf[contig_i][beg_pos:beg_pos+read_len])  # DEBUG
            if len(genome_slice) != read_len:
                raise ValueError ("unequal read length.  target read len: "+str(read_len)+" insilico read len: "+str(len(genome_slice)))
        else:
            #print("SOURCE: "+source_genome_buf[contig_i][beg_pos+1-read_len:beg_pos+1])  # DEBUG

            genome_slice = self.reverse_complement(source_genome_buf[contig_i][beg_pos+1-read_len:beg_pos+1])
            #print("COMPLT: "+genome_slice)  # DEBUG
            if len(genome_slice) != read_len:
                raise ValueError ("unequal read length.  target read len: "+str(read_len)+" insilico read len: "+str(len(genome_slice)))
        insilico_read_rec_buf[READ_SEQ_LINE_I] = genome_slice+"\n"

        # add in errors
        if add_errors_by_qual_freq:
            new_seq = ''
            for c_i,c in enumerate(genome_slice):
                qual = qual_line[c_i]
                new_seq += self.seq_call_with_error(c, qual, 'phred33')

            insilico_read_rec_buf[READ_SEQ_LINE_I] = new_seq+"\n"
                    
        # return
        return (insilico_pos, insilico_read_rec_buf)

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.handleURL = config['handle-service-url']
        self.serviceWizardURL = config['service-wizard-url']

        self.callbackURL = os.environ.get('SDK_CALLBACK_URL')
        if self.callbackURL == None:
            raise ValueError ("SDK_CALLBACK_URL not set in environment")

        self.scratch = os.path.abspath(config['scratch'])
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
        #sequencing_tech = 'N/A'  # no longer needed
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
        #sequencing_tech     = readsLibrary['files'][input_ref]['sequencing_tech'] # no longer needed


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
            try:
                readsUtils_Client = ReadsUtils (url=self.callbackURL, token=ctx['token'])  # SDK local
            except Exception as e:
                raise ValueError('Unable to get ReadsUtils Client' +"\n" + str(e))

            #self.add_id_to_plus_line(output_file_path)
            readsUtils_Client.upload_reads ({ 'wsname': str(params['workspace_name']),
                                              'name': params['output_name'],
                                              # remove sequencing_tech when source_reads_ref working
                                              #'sequencing_tech': sequencing_tech,
                                              'source_reads_ref': input_ref,
                                              'fwd_file': output_file_path
            })

                

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
            with open (input_fwd_file_path, 'r') as input_reads_file_handle:
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

            with open (input_rev_file_path, 'r') as input_reads_file_handle:
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
            output_reads_file_handle = open (output_rev_unpaired_file_path, 'w')
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

            with open (input_fwd_file_path, 'r') as input_reads_file_handle:
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
            output_reads_file_handle = open (output_fwd_unpaired_file_path, 'w')
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
                    #self.add_id_to_plus_line(output_fwd_paired_file_path)
                    #self.add_id_to_plus_line(output_rev_paired_file_path)
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
                #self.add_id_to_plus_line(output_fwd_unpaired_file_path)
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
                #self.add_id_to_plus_line(output_rev_unpaired_file_path)
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
            with open (input_fwd_file_path, 'r') as input_reads_file_handle:
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
                    #self.add_id_to_plus_line(output_fwd_paired_file_path)
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
            orig_reads_num = params['subsample_fraction']['reads_num']
            params['subsample_fraction']['reads_num'] = int (orig_reads_num/2.0 + 0.5)
            self.log (console, "Adjusting reads num to number of pairs.  Input reads num: "+str(orig_reads_num)+" new pairs num: "+str(params['subsample_fraction']['reads_num']))


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
            with open (input_fwd_file_path, 'r') as input_reads_file_handle:
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
                        # manage read_id edge case: e.g. @SRR5891520.1.1 (forward) & @SRR5891520.1.2 (reverse)                                                               
                        read_id = ''.join(read_id.rsplit('.',1)) # replace last '.' with ''   
                        read_id = re.sub ("[\/\.\_\-\:\;][012lrLRfrFR53]\'*$", "", read_id)
                        fwd_ids[read_id] = True

                        # DEBUG
#                        if rec_cnt % 100 == 0:
#                            self.log(console,"read_id: '"+str(read_id)+"'")
#                        rec_cnt += 1 


            # read reverse to determine paired
            self.log (console, "DETERMINING PAIRED IDS")  # DEBUG
            with open (input_rev_file_path, 'r') as input_reads_file_handle:
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
                        # manage read_id edge case: e.g. @SRR5891520.1.1 (forward) & @SRR5891520.1.2 (reverse)                                                               
                        read_id = ''.join(read_id.rsplit('.',1)) # replace last '.' with ''   
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

            #self.log (console, "len PAIRED_IDS_LIST: "+str(len(paired_ids_list)))
            #self.log (console, "reads_per_lib: "+str(reads_per_lib))
            #self.log (console, "split_num: "+str(params['subsample_fraction']['split_num']))

            for i,read_id in enumerate(random.sample (paired_ids_list, reads_per_lib * params['subsample_fraction']['split_num'])):
                if read_id not in paired_lib_i:
                    lib_i = i % params['subsample_fraction']['split_num']
                    paired_lib_i[read_id] = lib_i
                else:
                    raise ValueError ("repeated random sample for read id "+read_id)
            # not all reads are necessarily included.  Shouldn't check
            #for read_id in paired_ids_list:
            #    if read_id not in paired_lib_i:
            #        raise ValueError ("failed to assign output lib for read id "+read_id)

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

            with open (input_fwd_file_path, 'r') as input_reads_file_handle:
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
                        # manage read_id edge case: e.g. @SRR5891520.1.1 (forward) & @SRR5891520.1.2 (reverse)                                                               
                        read_id = ''.join(read_id.rsplit('.',1)) # replace last '.' with ''   
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

            with open (input_rev_file_path, 'r') as input_reads_file_handle:
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
                        # manage read_id edge case: e.g. @SRR5891520.1.1 (forward) & @SRR5891520.1.2 (reverse)                                                               
                        read_id = ''.join(read_id.rsplit('.',1)) # replace last '.' with ''   
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
                    #self.add_id_to_plus_line(output_fwd_paired_file_path)
                    #self.add_id_to_plus_line(output_rev_paired_file_path)
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

            with open (input_fwd_file_path, 'r') as input_reads_file_handle:
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
                        # manage read_id edge case: e.g. @SRR5891520.1.1 (forward) & @SRR5891520.1.2 (reverse)                                                               
                        read_id = ''.join(read_id.rsplit('.',1)) # replace last '.' with ''
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

            #self.log (console, "len PAIRED_IDS_LIST: "+str(len(paired_ids_list)))
            #self.log (console, "reads_per_lib: "+str(reads_per_lib))
            #self.log (console, "split_num: "+str(params['subsample_fraction']['split_num']))
            
            for i,read_id in enumerate(random.sample (paired_ids_list, reads_per_lib * params['subsample_fraction']['split_num'])):
                if read_id not in paired_lib_i:
                    lib_i = i % params['subsample_fraction']['split_num']
                    paired_lib_i[read_id] = lib_i
                else:
                    raise ValueError ("repeated random sample for read id "+read_id)
            # not all reads are necessarily included.  Shouldn't check
            #for read_id in paired_ids_list:
            #    if read_id not in paired_lib_i:
            #        raise ValueError ("failed to assign output lib for read id "+read_id)


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
            with open (input_fwd_file_path, 'r') as input_reads_file_handle:
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
                        # manage read_id edge case: e.g. @SRR5891520.1.1 (forward) & @SRR5891520.1.2 (reverse)                                                               
                        read_id = ''.join(read_id.rsplit('.',1)) # replace last '.' with ''   
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
                    #self.add_id_to_plus_line(output_fwd_paired_file_path)
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
            #self.add_id_to_plus_line(combined_input_fwd_path)
            #self.add_id_to_plus_line(combined_input_rev_path)
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
            #self.add_id_to_plus_line(combined_input_fwd_path)
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
        #combined_input_rev_path = os.path.join (input_dir, 'input_reads_rev.fastq')
        combined_input_fwd_handle = open (combined_input_fwd_path, 'w', write_buf_size)
        #combined_input_rev_handle = open (combined_input_rev_path, 'w', write_buf_size)


        # add libraries, one at a time
        sequencing_tech = None
        for this_input_reads_ref in readsSet_ref_list:
            clean_ref = re.sub("\/", "_", this_input_reads_ref)
            
            self.log (console, "DOWNLOADING FASTQ FILES FOR ReadsSet member: "+str(this_input_reads_ref))
            try:
                readsLibrary = readsUtils_Client.download_reads ({'read_libraries': [this_input_reads_ref],
                                                                  #'interleaved': 'false'
                                                                  'interleaved': 'true'
                                                                  })
            except Exception as e:
                raise ValueError('Unable to get reads object from workspace: (' + this_input_reads_ref +")\n" + str(e))

            this_input_fwd_path = readsLibrary['files'][this_input_reads_ref]['files']['fwd']

            # doing interleaved now
            #if input_reads_obj_type == "KBaseFile.PairedEndLibrary":
            #    this_input_rev_path = readsLibrary['files'][this_input_reads_ref]['files']['rev']

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
                this_header = None
                for line in this_input_handle:
                    rec_line_i += 1
                    if rec_line_i == 3:
                        rec_line_i = -1
                    elif rec_line_i == 0:
                        if not line.startswith('@'):
                            raise ValueError ("badly formatted rec line: '"+line+"'")
                        this_header = clean_ref+':'+line[1:]
                        line = '@'+this_header
                    elif rec_line_i == 2:
                        if not line.startswith('+'):
                            raise ValueError ("badly formatted rec line: '"+line+"'")
                        line = '+'+this_header
                    cat_file_handle.write(line)
            os.remove (this_input_path)  # create space since we no longer need the piece file

            # doing interleaved now
            """
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
            """

        combined_input_fwd_handle.close()
        #combined_input_rev_handle.close()  # doing interleaved now


        # upload reads
        #
        self.log (console, "UPLOADING MERGED READS LIB")  # DEBUG
        if not os.path.isfile (combined_input_fwd_path) \
                or os.path.getsize (combined_input_fwd_path) == 0:
            raise ValueError ("failed to create fwd read library output")
        #if input_reads_obj_type == "KBaseFile.PairedEndLibrary":
        #    if not os.path.isfile (combined_input_rev_path) \
        #        or os.path.getsize (combined_input_rev_path) == 0:
        #            
        #        raise ValueError ("failed to create rev read library output")

        output_obj_name = params['output_name']
        self.log(console, 'Uploading reads library: '+output_obj_name)

        if input_reads_obj_type == "KBaseFile.PairedEndLibrary":
            #self.add_id_to_plus_line(combined_input_fwd_path)
            reads_library_ref = readsUtils_Client.upload_reads ({ 'wsname': str(params['workspace_name']),
                                                                  'name': output_obj_name,
                                                                  # remove sequencing_tech when source_reads_ref is working
                                                                  #'sequencing_tech': sequencing_tech,
                                                                  'source_reads_ref': readsSet_ref_list[0],
                                                                  'single_genome': 0,
                                                                  'fwd_file': combined_input_fwd_path,
                                                                  #'rev_file': combined_input_rev_path
                                                                  'interleaved': 'true'
                                                                  })['obj_ref']
        else:
            #self.add_id_to_plus_line(combined_input_fwd_path)
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
            with open (input_rev_file_path, 'r') as input_reads_file_handle:
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
            with open (input_fwd_file_path, 'r') as input_reads_file_handle:
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

            with open (input_fwd_file_path, 'r') as input_reads_file_handle:
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

            with open (input_rev_file_path, 'r') as input_reads_file_handle:
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
                    #self.add_id_to_plus_line(output_fwd_paired_file_path)
                    #self.add_id_to_plus_line(output_rev_paired_file_path)
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
                    #self.add_id_to_plus_line(output_fwd_unpaired_file_path)
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
                    #self.add_id_to_plus_line(output_rev_unpaired_file_path)
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

            #qual33_fwd_path = this_input_fwd_path + '.qual33'
            qual33_fwd_path = re.sub ('.fastq$', '-q33.fastq', this_input_fwd_path)
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

                #qual33_rev_path = this_input_rev_path + '.qual33'
                qual33_rev_path = re.sub ('.fastq$', '-q33.fastq', this_input_rev_path)
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
                #self.add_id_to_plus_line(qual33_fwd_path)
                #self.add_id_to_plus_line(qual33_rev_path)
                reads_library_ref = readsUtils_Client.upload_reads ({ 'wsname': str(params['workspace_name']),
                                                                      'name': output_obj_name,
                                                                      # remove sequencing_tech when source_reads_ref is working
                                                                      #'sequencing_tech': sequencing_tech,
                                                                      'source_reads_ref': readsSet_ref_list[0],
                                                                      'fwd_file': qual33_fwd_path,
                                                                      'rev_file': qual33_rev_path
                                                                      })['obj_ref']
            else:
                #self.add_id_to_plus_line(qual33_fwd_path)
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

    def KButil_Generate_Microbiome_InSilico_Reads_From_Real_Reads(self, ctx, params):
        """
        :param params: instance of type
           "KButil_Generate_Microbiome_InSilico_Reads_From_Real_Reads_Params"
           (KButil_Generate_Microbiome_InSilico_Reads_From_Real_Reads() ** < 
           **  Method for random subsampling of reads library combined with
           overlay of configured genomes) -> structure: parameter
           "workspace_name" of type "workspace_name" (** The workspace object
           refs are of form: ** **    objects = ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_genomeSet_ref" of type "data_obj_ref", parameter
           "genome_abundances" of String, parameter "input_reads_ref" of type
           "data_obj_ref", parameter "output_name" of type "data_obj_name",
           parameter "subsample_fraction" of type "Fractionate_Options"
           (KButil_Random_Subsample_Reads() ** **  Method for random
           subsampling of reads library) -> structure: parameter "split_num"
           of Long, parameter "reads_num" of Long, parameter "reads_perc" of
           Double, parameter "genome_length_bias" of type "bool", parameter
           "desc" of String, parameter "pe_insert_len" of Long, parameter
           "pe_orientation" of String, parameter "seed" of Long
        :returns: instance of type
           "KButil_Generate_Microbiome_InSilico_Reads_From_Real_Reads_Output"
           -> structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN KButil_Generate_Microbiome_InSilico_Reads_From_Real_Reads
        console = []
        invalid_msgs = []
        self.log(console, 'Running KButil_Generate_Microbiome_InSilico_Reads_From_Real_Reads() with parameters: ')
        self.log(console, "\n"+pformat(params))
        report = ''
        
        token = ctx['token']
        headers = {'Authorization': 'OAuth '+token}
        env = os.environ.copy()
        env['KB_AUTH_TOKEN'] = token

        # dirs
        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()*1000)
        input_dir = os.path.join(self.scratch,'input.'+str(timestamp))
        if not os.path.exists(input_dir):
            os.makedirs(input_dir)
        output_dir = os.path.join(self.scratch,'output.'+str(timestamp))
        if not os.path.exists(input_dir):
            os.makedirs(input_dir)

        
        # Init clients
        #
        #SERVICE_VER = 'dev'  # DEBUG
        SERVICE_VER = 'release'
        try:
            wsClient = workspaceService(self.workspaceURL, token=token)            
            [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple
        except Exception as e:
            raise ValueError('Unable to get Workspace Client' +"\n" + str(e))
        try:
            readsUtils_Client = ReadsUtils (url=self.callbackURL, token=token, service_ver=SERVICE_VER)
        except Exception as e:
            raise ValueError('Unable to get ReadsUtils Client' +"\n" + str(e))
        try:
            dfuClient = DataFileUtil (url=self.callbackURL, token=token, service_ver=SERVICE_VER)
        except Exception as e:
            raise ValueError('Unable to get AssemblyUtil Client' +"\n" + str(e))
        try:
            AssemblyUtilClient = AssemblyUtil (url=self.callbackURL, token=token, service_ver=SERVICE_VER)
        except Exception as e:
            raise ValueError('Unable to get AssemblyUtil Client' +"\n" + str(e))
        try:
            setAPI_Client = SetAPI (url=self.serviceWizardURL, token=token, service_ver=SERVICE_VER)
        except Exception as e:
            raise ValueError('Unable to get SetAPI Client' +"\n" + str(e))
        try:
            reportClient = KBaseReport (self.callbackURL, token=token, service_ver=SERVICE_VER)
        except Exception as e:
            raise ValueError('Unable to get report Client' +"\n" + str(e))
        
        # init randomizer
        if 'seed' in params and params['seed'] != None:
            random.seed(params['seed'])
        else:
            random.seed()

        # param checks
        required_params = ['workspace_name',
                           'input_genomeSet_ref', 
                           'genome_abundances',
                           'input_reads_ref', 
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

        if not params.get('subsample_fraction'):
            raise ValueError ("Missing subsample_fraction params")
        if not int(params['subsample_fraction'].get('split_num',0)):
            raise ValueError ("Missing split_num")

        # use split_num to create reads_perc if neither reads_num or reads_perc defined
        use_reads_num  = False
        use_reads_perc = False
        if int(params['subsample_fraction'].get('reads_num',0)):
            self.log (console, "Ignoring reads_perc and just using reads_num: "+str(params['subsample_fraction']['reads_num']))
            use_reads_num  = True
            
        elif int(params['subsample_fraction'].get('reads_perc',0)) and float(params['subsample_fraction']['reads_perc']) <= 100:
            self.log (console, "Ignoring reads_num and just using reads_perc: "+str(params['subsample_fraction']['reads_perc']))
            use_reads_perc = True

        elif not int(params['subsample_fraction'].get('reads_num',0)) and not int(params['subsample_fraction'].get('reads_perc',0)):
            params['subsample_fraction']['reads_perc'] = int(100.0 * 1.0/float(params['subsample_fraction']['split_num']))
            use_reads_perc = True

        else:
            raise ValueError ("Badly configured subsample_fraction params")
            

        # load provenance
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        provenance[0]['input_ws_objects']=[]
        provenance[0]['input_ws_objects'].append(str(params['input_genomeSet_ref']))
        provenance[0]['input_ws_objects'].append(str(params['input_reads_ref']))


        # Determine whether read library is of correct type
        #
        try:
            input_reads_ref = params['input_reads_ref']
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


        # Determine whether paired end lib input params are defined
        #
        if input_reads_obj_type in PE_types:
            if not params.get('pe_orientation') \
               or not int(params.get('pe_insert_len',0)):
                raise ValueError ("require Paired-End Orientation and Paired-End insert length for Paired-End library input")


        # divide reads_num by 2 if paired end library because only counting pairs
        #
        if input_reads_obj_type in PE_types and int(params['subsample_fraction'].get('reads_num',0)):
            orig_reads_num = params['subsample_fraction']['reads_num']
            params['subsample_fraction']['reads_num'] = int (orig_reads_num/2.0 + 0.5)
            self.log (console, "Adjusting reads num to number of pairs.  Input reads num: "+str(orig_reads_num)+" new pairs num: "+str(params['subsample_fraction']['reads_num']))


        # Get Genome Set and object name <-> ref mapping
        #
        genome_name_by_ref = dict()
        genome_ref_by_name = dict()
        input_genomeSet_name = None
        try:
            #objects = wsClient.get_objects([{'ref': params['input_genomeset_ref']}])
            objects = wsClient.get_objects2({'objects': [{'ref': params['input_genomeSet_ref']}]})['data']
            genomeSet = objects[0]['data']
            info = objects[0]['info']

            input_genomeSet_name = info[NAME_I]
            type_name = info[TYPE_I].split('.')[1].split('-')[0]
            if type_name != 'GenomeSet':
                raise ValueError("Bad Type: Should be GenomeSet instead of '" + type_name + "'")
        except Exception as e:
            raise ValueError('Unable to fetch input_genomeSet_ref object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()

        for gId in genomeSet['elements'].keys():
            genomeRef = genomeSet['elements'][gId]['ref']
            try:
                already_included = genome_name_by_ref[genomeRef]
            except:
                try:
                    genome_obj_info = wsClient.get_object_info_new ({'objects':[{'ref':genomeRef}]})[0]
                    genomeName = genome_obj_info[NAME_I]
                except Exception as e:
                    raise ValueError('Unable to fetch genome object from workspace: ' + str(genomeRef) + str(e))
                    #to get the full stack trace: traceback.format_exc()
                genome_name_by_ref[genomeRef] = genomeName
                genome_ref_by_name[genomeName] = genomeRef

            
        # Determine relative unbiased percentages of genomes
        #
        genome_abundances = dict()
        genome_ref_order = []
        summed_abundance = 0.0
        if not params.get('genome_abundances'):
            for genome_ref in sorted(genome_name_by_ref.keys()):
                genome_ref_order.append(genome_ref)
            n_genomes = len(genome_ref_order)
            even_perc = 100.0 / n_genomes
            for genome_ref in genome_ref_order:
                genome_abundances[genome_ref] = even_perc
        else:
            ws_ref_pattern = re.compile('^\d+/\d+/\d+$')
            for genome_abundance_row in params['genome_abundances'].split("\n"):
                genome_abundance_row = genome_abundance_row.strip()
                [genome_id, abundance] = genome_abundance_row.split()

                # genome_ref
                if ws_ref_pattern.match(genome_id):
                    genome_ref = genome_id
                else:
                    try:
                        genome_ref = genome_ref_by_name[genome_id]
                    except:
                        raise ValueError ("Genome object: "+genome_id+" in genome abundances configuration not found in Genome Set: "+input_genomeSet_name)
                genome_ref_order.append(genome_ref)
                
                # abundance
                abundance = float(abundance)
                if abundance < 0  or  abundance > 100:
                    raise ValueError ("Genome object: "+genome_id+" ("+str(genome_ref)+") in genome abundances configuration has a non-percentile abundance: "+str(abundance))
                genome_abundances[genome_ref] = abundance
                summed_abundance += abundance
                
            # add unspecified genomes
            if summed_abundance > 100.1:
                raise ValueError ("Genome abundances sum to: "+str(summed_abundance)+".  Must be <= 100.0")
            elif summed_abundance < 99.9:
                abundance_remainder = 100.0 - summed_abundance
                genome_ref_remainder = []
                for genome_ref in sorted(genome_name_by_ref.keys()):
                    if genome_ref in genome_abundances:
                        continue
                    genome_ref_remainder.append(genome_ref)
                if len(genome_ref_remainder) == 0:
                    raise ValueError ("There were no unspecified genomes in Genome Set: "+input_genomeSet_name+" ("+str(input_genomeSet_ref)+") to assign the remaining abundance: "+str(abundance_remainder))

                # split remaining abundance evenly
                n_remaining_genomes = len(genome_ref_remainder)
                even_perc = abundance_remainder / n_remaining_genomes
                for genome_ref in genome_ref_remainder:
                    genome_ref_order.append(genome_ref)
                    genome_abundances[genome_ref] = even_perc


        # Check on progress
        #
        self.log(console, "TOTAL ABUNDANCE: "+str(summed_abundance))
        for genome_ref in genome_ref_order:
            self.log(console, "Genome Abundance: "+str(genome_abundances[genome_ref])+" for "+genome_name_by_ref[genome_ref])


        # Determine assembly lengths and Download source Assemblies
        #
        assembly_file_by_genome_ref = dict()
        assembly_len_by_genome_ref = dict()
        total_assembly_len = 0
        for genome_ref in genome_ref_order:
            assembly_len_by_genome_ref[genome_ref] = 0
            try:
                genome_obj = wsClient.get_objects([{'ref': genome_ref}])[0]['data']
            except:
                raise ValueError("unable to fetch genome: " + genome_ref)

            # contig lengths may already be in Genome obj
            if genome_obj.get('contig_ids') and genome_obj.get('contig_lengths'):
                for contig_i, contig_id in enumerate(genome_obj['contig_ids']):
                    contig_len = genome_obj['contig_lengths'][contig_i]
                    assembly_len_by_genome_ref[genome_ref] += contig_len
                    total_assembly_len += contig_len

            # Get genome_assembly_refs
            genome_assembly_ref = None
            genome_assembly_type = None
            if not genome_obj.get('contigset_ref') and not genome_obj.get('assembly_ref'):
                msg = "Genome " + genome_name_by_ref[genome_ref] + " ("+genome_ref+") "+ \
                    " MISSING BOTH contigset_ref AND assembly_ref.  Cannot process.  Exiting."
                raise ValueError(msg)
            elif genome_obj.get('assembly_ref'):
                assembly_ref = genome_obj['assembly_ref']
                genome_assembly_type = 'assembly'
                msg = "Genome " + genome_name_by_ref[genome_ref] + " ("+genome_ref+") "+ \
                    " USING assembly_ref: "+assembly_ref
                self.log(console, msg)
            elif genome_obj.get('contigset_ref'):
                assembly_ref = genome_obj['contigset_ref']
                genome_assembly_type = 'contigset'
                msg = "Genome " + genome_name_by_ref[genome_ref] + " ("+genome_ref+") "+ \
                    " USING contigset_ref: "+assembly_ref
                self.log(console, msg)
            else:
                raise ValueError ("Bad logic in finding assembly or contigset in genome object: "+genome_name_by_ref[genome_ref]+" ("+genome_ref+")")

            # Get and save assemblies
            contig_file = AssemblyUtilClient.get_assembly_as_fasta({'ref':assembly_ref}).get('path')
            sys.stdout.flush()
            contig_file_path = dfuClient.unpack_file({'file_path': contig_file})['file_path']
            assembly_file_by_genome_ref[genome_ref] = contig_file_path


        # Read assemblies into buffers
        #
        assembly_bufs_by_genome_ref = dict()
        contig_i_weighted_mapping_by_genome_ref = dict()
        for genome_ref in genome_ref_order:
            assembly_bufs_by_genome_ref[genome_ref] = []
            contig_i_weighted_mapping_by_genome_ref[genome_ref] = []

            # score fasta lens in contig files
            read_buf_size  = 65536
            write_buf_size = 65536

            with open (assembly_file_by_genome_ref[genome_ref], 'r', read_buf_size) as ass_handle:
                seq_buf = ''
                last_header = ''
                contig_i = 0
                for fasta_line in ass_handle:
                    fasta_line = fasta_line.strip()
                    if fasta_line.startswith('>'):
                        if seq_buf != '':
                            assembly_bufs_by_genome_ref[genome_ref].append(seq_buf)
                            for pos_i in range(len(seq_buf)):
                                contig_i_weighted_mapping_by_genome_ref[genome_ref].append(contig_i)
                            contig_i += 1
                            seq_buf = ''
                            last_header = fasta_line
                    else:
                        seq_buf += ''.join(fasta_line.split())
                if seq_buf != '':
                    assembly_bufs_by_genome_ref[genome_ref].append(seq_buf)
                    for pos_i in range(len(seq_buf)):
                        contig_i_weighted_mapping_by_genome_ref[genome_ref].append(contig_i)
                    seq_buf = ''


        # Add assembly lengths if not yet captured
        #
        for genome_ref in genome_ref_order:
            if assembly_len_by_genome_ref[genome_ref] == 0:
                for ass_buf in assembly_bufs_by_genome_ref[genome_ref]:
                    contig_len = len(ass_buf)
                    assembly_len_by_genome_ref[genome_ref] += contig_len
                    total_assembly_len += contig_len


        # Adjust relative abundance to include genome length bias
        #
        if int(params.get('genome_length_bias',0)):
            scaled_total = 0
            for genome_ref in genome_ref_order:
                length_bias = assembly_len_by_genome_ref[genome_ref] / float(total_assembly_len)
                self.log(console, "Genome Bias: "+str(length_bias)+" for "+genome_name_by_ref[genome_ref])  # DEBUG
                genome_abundances[genome_ref] *= length_bias
                self.log(console, "Biased Abundance: "+str(genome_abundances[genome_ref])+" for "+genome_name_by_ref[genome_ref])  # DEBUG


        # Rescale relative abundances to fix biased abundances (and fix near but off 100% totals)
        #
        adjusted_summed_abundance = 0
        corrected_summed_abundance = 0
        for genome_ref in genome_ref_order:
            adjusted_summed_abundance += genome_abundances[genome_ref]
        for genome_ref in genome_ref_order:
            genome_abundances[genome_ref] = 100 * genome_abundances[genome_ref] / adjusted_summed_abundance
            corrected_summed_abundance += genome_abundances[genome_ref]


        # Check on progress
        #
        self.log(console, "REVISED TOTAL ABUNDANCE: "+str(adjusted_summed_abundance))
        self.log(console, "CORRECTED TOTAL ABUNDANCE: "+str(corrected_summed_abundance))
        for genome_ref in genome_ref_order:
            self.log(console, "Corrected Genome Abundance: "+str(genome_abundances[genome_ref]))


        # Build genome lookup for random sampler
        #
        genome_lookup = []
        genome_lookup_resolution = 1000000
        for genome_i,genome_ref in enumerate(genome_ref_order):
            for i in range(int(genome_lookup_resolution * genome_abundances[genome_ref]/100.0)):
                genome_lookup.append(genome_i)
        filled = len(genome_lookup)
        if filled < genome_lookup_resolution:
            for i in range(genome_lookup_resolution-filled):
                genome_i = i % len(genome_ref_order)
                genome_lookup.append(genome_i)


        # Download Reads
        #
        self.log (console, "DOWNLOADING READS")  # DEBUG
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
            source_genome_i = dict()
            read_ids_by_lib = []
            paired_buf_size = 100000
            recs_beep_n = 1000000

            # read fwd file to get fwd ids
#            rec_cnt = 0  # DEBUG
            self.log (console, "GETTING IDS")  # DEBUG
            with open (input_fwd_file_path, 'r') as input_reads_file_handle:
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
            with open (input_rev_file_path, 'r') as input_reads_file_handle:
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
            self.log (console, "GETTING RANDOM READ SUBSAMPLES")  # DEBUG
            for lib_i in range(int(params['subsample_fraction']['split_num'])):
                read_ids_by_lib.append([])
            for i,read_id in enumerate(random.sample (paired_ids_list, reads_per_lib * int(params['subsample_fraction']['split_num']))):
                if read_id not in paired_lib_i:
                    lib_i = i % int(params['subsample_fraction']['split_num'])
                    paired_lib_i[read_id] = lib_i
                    read_ids_by_lib[lib_i].append(read_id)
                else:
                    raise ValueError ("repeated random sample for read id "+read_id)
            for read_id in paired_ids_list:
                if read_id not in paired_lib_i:
                    raise ValueError ("failed to assign output lib for read id "+read_id)

            # Determine random source genome for each read
            self.log (console, "GETTING RANDOM GENOME SOURCE ASSIGNMENT FOR READS")  # DEBUG
            for lib_i in range(len(read_ids_by_lib)):
                for read_i,genome_lookup_i_map in enumerate(random.sample (genome_lookup, reads_per_lib)):
                    source_genome_i[read_ids_by_lib[lib_i][read_i]] = genome_lookup_i_map

            # split fwd paired
            self.log (console, "WRITING FWD SPLIT PAIRED")  # DEBUG
            paired_output_reads_file_handles = []
            for lib_i in range(int(params['subsample_fraction']['split_num'])):
                paired_output_reads_file_handles.append(open (output_fwd_paired_file_path_base+"-"+str(lib_i)+".fastq", 'w', paired_buf_size))
                total_paired_reads_by_set.append(0)

            rec_buf = []
            last_read_id = None
            paired_cnt = 0
            capture_type_paired = False
            fwd_insilico_position = dict()

            with open (input_fwd_file_path, 'r') as input_reads_file_handle:
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
                                try:
                                    lib_i = paired_lib_i[last_read_id]
                                except:
                                    raise ValueError ("Failed to find lib for read rec "+last_read_id)

                                #if paired_cnt % recs_beep_n == 0:
                                #    self.log(console,"OLD FWD READ:")  # DEBUG
                                #    self.log(console,"\n".join(rec_buf)) # DEBUG

                                try:
                                    # replace read with source genome sequence
                                    source_genome_ref = genome_ref_order[source_genome_i[last_read_id]]
                                except:
                                    raise ValueError ("Failed to find source genome for read rec "+last_read_id)

                                try:
                                    (fwd_insilico_position[last_read_id], insilico_read_rec_buf) = self.overlay_source_genome_seq (
                                        read_rec=rec_buf, 
                                        source_genome_buf=assembly_bufs_by_genome_ref[source_genome_ref],
                                        contig_mapping=contig_i_weighted_mapping_by_genome_ref[source_genome_ref],
                                    
                                        lib_type='PE',
                                        read_dir='fwd',
                                        fwd_insilico_pos=None,
                                        pe_orientation=params['pe_orientation'],
                                        pe_insert_len=int(params['pe_insert_len']),
                                        add_errors_by_qual_freq=True
                                    )
                                except:
                                    raise ValueError ("Failed to transform read rec "+last_read_id+"\n"+"".join(rec_buf))

                                #if paired_cnt % recs_beep_n == 0:
                                #    self.log(console,"NEW FWD READ:")  # DEBUG
                                #    self.log(console,"\n".join(insilico_read_rec_buf)) # DEBUG

                                try:
                                    # write rec to file
                                    paired_output_reads_file_handles[lib_i].writelines(insilico_read_rec_buf)
                                except:
                                    raise ValueError ("Failed to write read rec "+last_read_id+"\n"+"".join(insilico_read_rec_buf))

                                    
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
                        try:
                            lib_i = paired_lib_i[last_read_id]
                        except:
                            raise ValueError ("Failed to find lib for read rec "+last_read_id)

                        try:
                            # replace read with source genome sequence
                            source_genome_ref = genome_ref_order[source_genome_i[last_read_id]]
                        except:
                            raise ValueError ("Failed to find source genome for read rec "+last_read_id)
                        try:
                            (fwd_insilico_position[last_read_id], insilico_read_rec_buf) = self.overlay_source_genome_seq (read_rec=rec_buf, 
                                source_genome_buf=assembly_bufs_by_genome_ref[source_genome_ref],
                                contig_mapping=contig_i_weighted_mapping_by_genome_ref[source_genome_ref],
                                lib_type='PE',
                                read_dir='fwd',
                                fwd_insilico_pos=None,
                                pe_orientation=params['pe_orientation'],
                                pe_insert_len=int(params['pe_insert_len']),
                                add_errors_by_qual_freq=True
                            )
                        except:
                            raise ValueError ("Failed to transform read rec "+last_read_id+"\n"+"".join(rec_buf))
                        try:
                            # write rec to file
                            paired_output_reads_file_handles[lib_i].writelines(insilico_read_rec_buf)
                        except:
                            raise ValueError ("Failed to write read rec "+last_read_id+"\n"+"".join(insilico_read_rec_buf))

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

            with open (input_rev_file_path, 'r') as input_reads_file_handle:
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
                                try:
                                    lib_i = paired_lib_i[last_read_id]
                                except:
                                    raise ValueError ("Failed to find lib for read rec "+last_read_id)

                                #if paired_cnt % recs_beep_n == 0:
                                #    self.log(console,"OLD REV READ:")  # DEBUG
                                #    self.log(console,"\n".join(rec_buf)) # DEBUG

                                try:
                                    # replace read with source genome sequence
                                    source_genome_ref = genome_ref_order[source_genome_i[last_read_id]]
                                except:
                                    raise ValueError ("Failed to find source genome for read rec "+last_read_id)

                                try:
                                    (unused, insilico_read_rec_buf) = self.overlay_source_genome_seq (read_rec=rec_buf, 
                                        source_genome_buf=assembly_bufs_by_genome_ref[source_genome_ref],
                                        contig_mapping=contig_i_weighted_mapping_by_genome_ref[source_genome_ref],
                                        lib_type='PE',
                                        read_dir='rev',
                                        fwd_insilico_pos=fwd_insilico_position[last_read_id],
                                        pe_orientation=params['pe_orientation'],
                                        pe_insert_len=int(params['pe_insert_len']),
                                        add_errors_by_qual_freq=True
                                    )
                                except:
                                    raise ValueError ("Failed to transform read rec "+last_read_id+"\n"+"".join(rec_buf))

                                #if paired_cnt % recs_beep_n == 0:
                                #    self.log(console,"NEW REV READ:")  # DEBUG
                                #    self.log(console,"\n".join(insilico_read_rec_buf)) # DEBUG

                                try:
                                    # write rec to file
                                    paired_output_reads_file_handles[lib_i].writelines(insilico_read_rec_buf)
                                except:
                                    raise ValueError ("Failed to write read rec "+last_read_id+"\n"+"".join(insilico_read_rec_buf))

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
                        try:
                            lib_i = paired_lib_i[last_read_id]
                        except:
                            raise ValueError ("Failed to find lib for read rec "+last_read_id)

                        try:
                            # replace read with source genome sequence
                            source_genome_ref = genome_ref_order[source_genome_i[last_read_id]]
                        except:
                            raise ValueError ("Failed to find source genome for read rec "+last_read_id)

                        try:
                            (unused, insilico_read_rec_buf) = self.overlay_source_genome_seq (read_rec=rec_buf, 
                                source_genome_buf=assembly_bufs_by_genome_ref[source_genome_ref],
                                contig_mapping=contig_i_weighted_mapping_by_genome_ref[source_genome_ref],
                                lib_type='PE',
                                read_dir='rev',
                                fwd_insilico_pos=fwd_insilico_position[last_read_id],
                                pe_orientation=params['pe_orientation'],
                                pe_insert_len=int(params['pe_insert_len']),
                                add_errors_by_qual_freq=True
                            )
                        except:
                            raise ValueError ("Failed to transform read rec "+last_read_id+"\n"+"".join(rec_buf))

                        try:
                            # write rec to file
                            paired_output_reads_file_handles[lib_i].writelines(insilico_read_rec_buf)
                        except:
                            raise ValueError ("Failed to write read rec "+last_read_id+"\n"+"".join(insilico_read_rec_buf))

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
            self.log (console, "UPLOAD In Silico PAIRED READS LIBS")  # DEBUG
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
                    #self.add_id_to_plus_line(output_fwd_paired_file_path)
                    #self.add_id_to_plus_line(output_rev_paired_file_path)
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
            total_paired_reads = 0
            total_paired_reads_by_set = []
            paired_ids = dict()
            paired_ids_list = []
            paired_lib_i = dict()
            source_genome_i = dict()
            read_ids_by_lib = []
            paired_buf_size = 100000
            recs_beep_n = 1000000

            with open (input_fwd_file_path, 'r') as input_reads_file_handle:
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
            self.log (console, "GETTING RANDOM READ SUBSAMPLES")  # DEBUG
            for lib_i in range(int(params['subsample_fraction']['split_num'])):
                read_ids_by_lib.append([])
            for i,read_id in enumerate(random.sample (paired_ids_list, reads_per_lib * params['subsample_fraction']['split_num'])):
                if read_id not in paired_lib_i:
                    lib_i = i % int(params['subsample_fraction']['split_num'])
                    paired_lib_i[read_id] = lib_i
                    read_ids_by_lib[lib_i].append(read_id)
                else:
                    raise ValueError ("repeated random sample for read id "+read_id)
            for read_id in paired_ids_list:
                if read_id not in paired_lib_i:
                    raise ValueError ("failed to assign output lib for read id "+read_id)

            # Determine random source genome for each read
            self.log (console, "GETTING RANDOM GENOME SOURCE ASSIGNMENT FOR READS")  # DEBUG
            for lib_i in range(len(read_ids_by_lib)):
                for read_i,genome_lookup_i_map in enumerate(random.sample (genome_lookup, reads_per_lib)):
                    source_genome_i[read_ids_by_lib[lib_i][read_i]] = genome_lookup_i_map


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
            with open (input_fwd_file_path, 'r') as input_reads_file_handle:
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
                            except:
                                raise ValueError ("Failed to find lib for read rec "+last_read_id)

                            # DEBUG
                            print ("LIB_I: "+str(lib_i)+" LAST_READ_ID: "+last_read_id)

                            #if paired_cnt % recs_beep_n == 0:
                            #    self.log(console,"OLD FWD READ:")  # DEBUG
                            #    self.log(console,"\n".join(rec_buf)) # DEBUG

                            try:
                                # replace read with source genome sequence
                                source_genome_ref = genome_ref_order[source_genome_i[last_read_id]]
                            except:
                                raise ValueError ("Failed to find source genome for read rec "+last_read_id)

                            try:
                                (unused, insilico_read_rec_buf) = self.overlay_source_genome_seq (
                                    read_rec=rec_buf, 
                                    source_genome_buf=assembly_bufs_by_genome_ref[source_genome_ref],
                                    contig_mapping=contig_i_weighted_mapping_by_genome_ref[source_genome_ref],
                                    lib_type='SE',
                                    read_dir='fwd',
                                    fwd_insilico_pos=None,
                                    pe_orientation=None,
                                    pe_insert_len=None,
                                    add_errors_by_qual_freq=True
                                )
                            except:
                                raise ValueError ("Failed to transform read rec "+last_read_id+"\n"+"".join(rec_buf))

                            #if paired_cnt % recs_beep_n == 0:
                            #    self.log(console,"NEW FWD READ:")  # DEBUG
                            #    self.log(console,"\n".join(insilico_read_rec_buf)) # DEBUG

                            try:
                                # write rec to file
                                paired_output_reads_file_handles[lib_i].writelines(insilico_read_rec_buf)                                
                            except:
                                raise ValueError ("Failed to write read rec "+last_read_id+"\n"+"".join(insilico_read_rec_buf))
                                
                            paired_cnt += 1
                            total_paired_reads_by_set[lib_i] += 1
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
                        except:
                            raise ValueError ("Failed to find lib for read rec "+last_read_id)

                        # DEBUG
                        print ("LIB_I: "+str(lib_i)+" LAST_READ_ID: "+last_read_id)

                        try:
                            # replace read with source genome sequence
                            source_genome_ref = genome_ref_order[source_genome_i[last_read_id]]
                        except:
                            raise ValueError ("Failed to find source genome for read rec "+last_read_id)

                        try:
                            (unused, insilico_read_rec_buf) = self.overlay_source_genome_seq (
                                read_rec=rec_buf, 
                                source_genome_buf=assembly_bufs_by_genome_ref[source_genome_ref],
                                contig_mapping=contig_i_weighted_mapping_by_genome_ref[source_genome_ref],
                                lib_type='SE',
                                read_dir='fwd',
                                fwd_insilico_pos=None,
                                pe_orientation=None,
                                pe_insert_len=None,
                                add_errors_by_qual_freq=True
                            )
                        except:
                            raise ValueError ("Failed to transform read rec "+last_read_id+"\n"+"".join(rec_buf))

                        try:
                            # write rec to file
                            paired_output_reads_file_handles[lib_i].writelines(insilico_read_rec_buf)                                
                        except:
                            raise ValueError ("Failed to write read rec "+last_read_id+"\n"+"".join(insilico_read_rec_buf))
                            
                        paired_cnt += 1
                        total_paired_reads_by_set[lib_i] += 1
                        if paired_cnt != 0 and paired_cnt % recs_beep_n == 0:
                            self.log(console,"\t"+str(paired_cnt)+" recs processed")
                    rec_buf = []
                    last_read_id = None

            for output_handle in paired_output_reads_file_handles:
                output_handle.close()


            # store report
            #
            report += "TOTAL READS: "+str(total_paired_reads)+"\n"
            for lib_i in range(params['subsample_fraction']['split_num']):
                report += "SINGLE END READS IN SET "+str(lib_i)+": "+str(total_paired_reads_by_set[lib_i])+"\n"


            # upload reads
            #
            self.log (console, "UPLOAD In Silico SINGLE END READS")  # DEBUG
            paired_obj_refs = []
            for lib_i in range(params['subsample_fraction']['split_num']):
                output_fwd_paired_file_path = output_fwd_paired_file_path_base+"-"+str(lib_i)+".fastq"
                if not os.path.isfile (output_fwd_paired_file_path) \
                        or os.path.getsize (output_fwd_paired_file_path) == 0:
                    
                    raise ValueError ("failed to create single end library output")
                else:
                    output_obj_name = params['output_name']+'-'+str(lib_i)
                    self.log(console, 'Uploading single end reads: '+output_obj_name)
                    #self.add_id_to_plus_line(output_fwd_paired_file_path)
                    paired_obj_refs.append (readsUtils_Client.upload_reads ({ 'wsname': str(params['workspace_name']),
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
        report_info = reportClient.create({'report':reportObj, 'workspace_name':params['workspace_name']})

        returnVal = { 'report_name': report_info['name'], 'report_ref': report_info['ref'] }
        #END KButil_Generate_Microbiome_InSilico_Reads_From_Real_Reads

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method KButil_Generate_Microbiome_InSilico_Reads_From_Real_Reads return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def KButil_Fractionate_Reads_by_Contigs(self, ctx, params):
        """
        :param params: instance of type
           "KButil_Fractionate_Reads_by_Contigs_Params"
           (KButil_Fractionate_Reads_by_Contigs() ** **  Split reads library
           into two based on whether they match contigs) -> structure:
           parameter "workspace_name" of type "workspace_name" (** The
           workspace object refs are of form: ** **    objects =
           ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_reads_ref" of type "data_obj_ref", parameter
           "input_assembly_ref" of type "data_obj_ref", parameter
           "output_name" of type "data_obj_name", parameter
           "fractionate_mode" of String
        :returns: instance of type
           "KButil_Fractionate_Reads_by_Contigs_Output" -> structure:
           parameter "report_name" of type "data_obj_name", parameter
           "report_ref" of type "data_obj_ref", parameter
           "source_reads_count" of Long, parameter "positive_reads_count" of
           Long, parameter "negative_reads_count" of Long, parameter
           "source_reads_sum_length" of Long, parameter
           "positive_reads_sum_length" of Long, parameter
           "negative_reads_sum_length" of Long
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN KButil_Fractionate_Reads_by_Contigs

        # very strange, re import from above isn't being retained in this scope
        import re

        #### Step 0: basic init
        ##
        console = []
        invalid_msgs = []
        report_text = ''
        self.log(console, 'Running run_fractionate_contigs(): ')
        self.log(console, "\n"+pformat(params))

        # Auth
        token = ctx['token']
        headers = {'Authorization': 'OAuth '+token}
        env = os.environ.copy()
        env['KB_AUTH_TOKEN'] = token

        # API Clients
        #SERVICE_VER = 'dev'  # DEBUG
        SERVICE_VER = 'release'
        # wsClient
        try:
            wsClient = Workspace(self.workspaceURL, token=token)
        except Exception as e:
            raise ValueError('Unable to instantiate wsClient with workspaceURL: '+ self.workspaceURL +' ERROR: ' + str(e))
        # dfuClient
        try:
            dfuClient = DataFileUtil(self.callbackURL)
        except Exception as e:
            raise ValueError('Unable to instantiate dfuClient with callbackURL: '+ self.callbackURL +' ERROR: ' + str(e))
        # gfuClient
        try:
            gfuClient = GenomeFileUtil(self.callbackURL)
        except Exception as e:
            raise ValueError('Unable to instantiate gfuClient with callbackURL: '+ self.callbackURL +' ERROR: ' + str(e))
        # setAPI_Client
        try:
            #setAPI_Client = SetAPI (url=self.callbackURL, token=ctx['token'])  # for SDK local.  local doesn't work for SetAPI
            setAPI_Client = SetAPI (url=self.serviceWizardURL, token=ctx['token'])  # for dynamic service
        except Exception as e:
            raise ValueError('Unable to instantiate setAPI_Client with serviceWizardURL: '+ self.serviceWizardURL +' ERROR: ' + str(e))
        # auClient
        try:
            auClient = AssemblyUtil(self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
        except Exception as e:
            raise ValueError('Unable to instantiate auClient with callbackURL: '+ self.callbackURL +' ERROR: ' + str(e))

        # param checks
        required_params = ['workspace_name',
                           'input_reads_ref',
                           'input_assembly_ref',
                           'fractionate_mode',
                           'output_name'
                          ]
        for arg in required_params:
            if arg not in params or params[arg] == None or params[arg] == '':
                raise ValueError ("Must define required param: '"+arg+"'")

        # load provenance
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        provenance[0]['input_ws_objects']=[]
        provenance[0]['input_ws_objects'].append(params['input_reads_ref'])
        provenance[0]['input_ws_objects'].append(params['input_assembly_ref'])

        # set the output paths
        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()*1000)
        output_dir = os.path.join(self.scratch,'output.'+str(timestamp))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        html_output_dir = os.path.join(output_dir,'html')
        if not os.path.exists(html_output_dir):
            os.makedirs(html_output_dir)

        # obj info fields
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple
            
        # configure data types
        assembly_obj_type = "KBaseGenomeAnnotations.Assembly"
        #assembly_set_obj_type = "KBaseSets.AssemblySet"
        genome_obj_type = "KBaseGenomes.Genome"
        #genome_set_obj_type = "KBaseSearch.GenomeSet"
        #binned_contigs_obj_type = "KBaseMetagenomes.BinnedContigs"
        ama_assembly_obj_type = "KBaseMetagenomes.AnnotatedMetagenomeAssembly"

            
        #### STEP 1: Get assembly to fractionate as file
        ##
        input_assembly_name = None
        input_assembly_type = None
        contig_file_path = None
        gff_file_path = None
        src_contig_ids = []

        if len(invalid_msgs) == 0:
            self.log (console, "Retrieving Assembly for pos filter")  # DEBUG

            # assembly obj info
            input_ref = params['input_assembly_ref']
            accepted_input_types = [assembly_obj_type, ama_assembly_obj_type]
            try:
                #input_obj_info = wsClient.get_object_info_new ({'objects':[{'ref':input_ref}]})[0]
                assembly_obj = wsClient.get_objects2({'objects':[{'ref':input_ref}]})['data'][0]
                input_obj_data = assembly_obj['data']
                input_obj_info = assembly_obj['info']

                input_obj_type = re.sub ('-[0-9]+\.[0-9]+$', "", input_obj_info[TYPE_I])  # remove trailing version
                input_obj_name = input_obj_info[NAME_I]
            except Exception as e:
                raise ValueError('Unable to get object from workspace: (' + input_ref +'): ' + str(e))

            #self.log(console, "\nASSEMBLY_OBJ:\n"+pformat(input_obj_data))  # DEBUG
            input_assembly_type = input_obj_type
            input_assembly_name = input_obj_name

            if input_obj_type not in accepted_input_types:
                raise ValueError ("Input object of type '"+input_obj_type+"' not accepted.  Must be one of "+", ".join(accepted_input_types))

            # Assembly type
            if input_obj_type == assembly_obj_type:
                self.log(console, "getting contigs for Assembly "+input_obj_name)
                for contig_key in sorted (input_obj_data['contigs'].keys()):
                    src_contig_ids.append(input_obj_data['contigs'][contig_key]['contig_id'])

                contig_file = auClient.get_assembly_as_fasta({'ref':input_ref}).get('path')
                sys.stdout.flush()
                contig_file_path = dfuClient.unpack_file({'file_path': contig_file})['file_path']

                
            # AnnotatedMetagenomeAssembly type
            else:
                self.log(console, "getting contigs for AnnotatedMetagenomeAssembly "+input_obj_name)
                src_contig_ids.extend(input_obj_data['contig_ids'])

                contig_file = auClient.get_fastas({'ref_lst':[input_ref]})[input_ref]['paths'][0]
                sys.stdout.flush()
                contig_file_path = dfuClient.unpack_file({'file_path': contig_file})['file_path']
                gff_file = gfuClient.metagenome_to_gff({'metagenome_ref':input_ref}).get('file_path')
                sys.stdout.flush()
                gff_file_path = dfuClient.unpack_file({'file_path': gff_file})['file_path']


        #### STEP 2: get positive filter object contig IDs as contig files
        ##
        if len(invalid_msgs) == 0:

            accepted_input_types = [assembly_obj_type,
                                    assembly_set_obj_type,
                                    genome_obj_type,
                                    genome_set_obj_type,
                                    binned_contigs_obj_type,
                                    ama_assembly_obj_type]
            pos_filter_contig_ids = dict()
        
            for i,input_ref in enumerate(params['input_pos_filter_obj_refs']):
                # get object
                try:
                    #input_obj_info = wsClient.get_object_info_new ({'objects':[{'ref':input_ref}]})[0]
                    input_obj = wsClient.get_objects2({'objects':[{'ref':input_ref}]})['data'][0]
                    input_obj_data = input_obj['data']
                    input_obj_info = input_obj['info']
                    input_obj_type = re.sub ('-[0-9]+\.[0-9]+$', "", input_obj_info[TYPE_I])  # remove trailing version
                    input_obj_name = input_obj_info[NAME_I]
                except Exception as e:
                    raise ValueError('Unable to get object from workspace: (' + input_ref +'): ' + str(e))
                if input_obj_type not in accepted_input_types:
                    raise ValueError ("Input object "+input_obj_name+" of type '"+input_obj_type+"' not accepted.  Must be one of "+", ".join(accepted_input_types))

                # Assembly
                if input_obj_type == assembly_obj_type:
                    for contig_key in input_obj_data['contigs'].keys():
                        pos_filter_contig_ids[input_obj_data['contigs'][contig_key]['contig_id']] = True

                # Genome or AnnotatedMetagenomeAssembly
                elif input_obj_type == genome_obj_type \
                     or input_obj_type == ama_assembly_obj_type:

                    for contig_id in input_obj_data['contig_ids']:
                        pos_filter_contig_ids[contig_id] = True

                # AssemblySet
                elif input_obj_type == assembly_set_obj_type:
                    try:
                        set_obj = setAPI_Client.get_assembly_set_v1 ({'ref':input_ref, 'include_item_info':1})
                    except Exception as e:
                        raise ValueError('Unable to get object '+input_obj_name+' from workspace: (' + input_ref +')' + str(e))

                    for item_obj in set_obj['data']['items']:
                        this_input_ref = item_obj['ref']
                        try:
                            this_input_obj = wsClient.get_objects2({'objects':[{'ref':this_input_ref}]})['data'][0]
                            this_input_obj_data = this_input_obj['data']
                        except Exception as e:
                            raise ValueError('Unable to get object from workspace: (' + this_input_ref +'): ' + str(e))
                        for contig_key in this_input_obj_data['contigs'].keys():
                            pos_filter_contig_ids[this_input_obj_data['contigs'][contig_key]['contig_id']] = True

                # GenomeSet
                elif input_obj_type == genome_set_obj_type:
                    # use SetAPI when it sends back 'items' for KBaseSearch.GenomeSet
                    #try:
                    #    set_obj = setAPI_Client.get_genome_set_v1 ({'ref':input_ref, 'include_item_info':1})
                    #except Exception as e:
                    #    raise ValueError('Unable to get object '+input_obj_name+' from workspace: (' + input_ref +')' + str(e))
                    # for now use this mimic to get the same 'items' structure
                    set_obj = wsClient.get_objects2({'objects':[{'ref':input_ref}]})['data'][0]
                    set_obj['data']['items'] = []
                    for element_key in set_obj['data']['elements'].keys():
                        set_obj['data']['items'].append(set_obj['data']['elements'][element_key])

                    for item_obj in set_obj['data']['items']:
                        this_input_ref = item_obj['ref']
                        try:
                            this_input_obj = wsClient.get_objects2({'objects':[{'ref':this_input_ref}]})['data'][0]
                            this_input_obj_data = this_input_obj['data']
                        except Exception as e:
                            raise ValueError('Unable to get object from workspace: (' + this_input_ref +'): ' + str(e))
                        for contig_id in this_input_obj_data['contig_ids']:
                            pos_filter_contig_ids[contig_id] = True

                # BinnedContigs
                elif input_obj_type == binned_contigs_obj_type:
                    for bin in input_obj_data['bins']:
                        for contig_id in bin['contigs'].keys():
                            pos_filter_contig_ids[contig_id] = True
                
                else:
                    raise ValueError ("data type "+input_obj_type+" for object "+input_obj_name+" not handled")
                    

        #### STEP 3: Determine pos and neg contig ids.
        ## Note: pos_contig_ids may be smaller than pos_filter_contig_ids
        ##
        if len(invalid_msgs) == 0:
            pos_contig_ids = dict()
            neg_contig_ids = dict()
            if len(invalid_msgs) == 0:
                self.log (console, "Filtering Contig IDs")
                for contig_id in src_contig_ids:
                    if pos_filter_contig_ids.get(contig_id):
                        pos_contig_ids[contig_id] = True
                    else:
                        neg_contig_ids[contig_id] = True


        #### STEP 4: Create pos and neg output files
        ##
        if len(invalid_msgs) == 0:
            pos_contig_file_path = None
            neg_contig_file_path = None
            orig_contig_count = 0
            pos_contig_count = 0
            neg_contig_count = 0
            orig_contig_len = 0
            pos_contig_len = 0
            neg_contig_len = 0

            self.log (console, "Fractionating contigs in assembly: "+input_assembly_name)
            pos_contig_file_path = contig_file_path+".positive_fraction.fa"
            neg_contig_file_path = contig_file_path+".negative_fraction.fa"
            with open (contig_file_path, 'r') as src_handle, \
                 open (pos_contig_file_path, 'w') as pos_handle, \
                 open (neg_contig_file_path, 'w') as neg_handle:
                seq_buf = ''
                last_header = ''
                for fasta_line in src_handle:
                    if fasta_line.startswith('>'):
                        if seq_buf != '':
                            seq_len = len(seq_buf)
                            orig_contig_count += 1
                            orig_contig_len += seq_len
                            contig_id = last_header.split()[0].replace('>','')
                            if pos_contig_ids.get(contig_id):
                                pos_contig_count += 1
                                pos_contig_len += seq_len
                                pos_handle.write(last_header)  # last_header already has newline
                                pos_handle.write(seq_buf+"\n")
                            else:
                                neg_contig_count += 1
                                neg_contig_len += seq_len
                                neg_handle.write(last_header)  # last_header already has newline
                                neg_handle.write(seq_buf+"\n")
                        seq_buf = ''
                        last_header = fasta_line
                    else:
                        seq_buf += ''.join(fasta_line.split())
                # last rec
                if seq_buf != '':
                    seq_len = len(seq_buf)
                    orig_contig_count += 1
                    orig_contig_len += seq_len
                    contig_id = last_header.split()[0].replace('>','')
                    if pos_contig_ids.get(contig_id):
                        pos_contig_count += 1
                        pos_contig_len += seq_len
                        pos_handle.write(last_header)  # last_header already has newline
                        pos_handle.write(seq_buf+"\n")
                    else:
                        neg_contig_count += 1
                        neg_contig_len += seq_len
                        neg_handle.write(last_header)  # last_header already has newline
                        neg_handle.write(seq_buf+"\n")
                    seq_buf = ''
                    last_header = ''


        #### STEP 5: for AnnotatedMetagenomeAssembly input, Create pos and neg gff output files
        ##
        if len(invalid_msgs) == 0 and input_assembly_type == ama_assembly_obj_type:
            pos_gff_file_path = None
            neg_gff_file_path = None
            orig_feature_count = 0
            pos_feature_count = 0
            neg_feature_count = 0

            self.log (console, "Fractionating genes in AnnotatedMetagenomeAssembly: "+input_assembly_name)
            pos_gff_file_path = contig_file_path+".positive_fraction.gff"
            neg_gff_file_path = contig_file_path+".negative_fraction.gff"
            with open (gff_file_path, 'r') as src_handle, \
                 open (pos_gff_file_path, 'w') as pos_handle, \
                 open (neg_gff_file_path, 'w') as neg_handle:

                for gff_line in src_handle:
                    if gff_line.startswith('#'):
                        pos_handle.write(gff_line)  # gff_line already has newline
                        neg_handle.write(gff_line)
                    else:
                        orig_feature_count += 1
                        contig_id = gff_line.split()[0]
                        if pos_contig_ids.get(contig_id):
                            pos_handle.write(gff_line)  # gff_line already has newline
                            pos_feature_count += 1
                        else:
                            neg_handle.write(gff_line)  # gff_line already has newline
                            neg_feature_count += 1
            
                            
        #### STEP 6: save the fractionated assemblies
        ##
        if pos_contig_count == 0 and neg_contig_count == 0:
            self.log (invalid_msgs, "Contig IDs don't match between "+input_assembly_name+" and positive filter objects.")
        if pos_contig_count == 0:
            self.log (invalid_msgs, "No positive matches for Contig IDs.  Fractionation is unnecessary.")
        if neg_contig_count == 0:
            self.log (invalid_msgs, "No negative matches for Contig IDs.  Fractionation is unnecessary.")

        if len(invalid_msgs) == 0:
            fractionated_obj_refs = []
            fractionated_obj_names = []

            # Pos fraction
            if not params.get('fractionate_mode') \
               or params['fractionate_mode'] == 'both' \
               or params['fractionate_mode'] == 'pos':

                if params.get('fractionate_mode') == 'pos':
                    output_obj_name = params['output_name']
                else:
                    output_obj_name = params['output_name']+'-positive_fraction'

                if input_assembly_type == assembly_obj_type:
                    output_data_ref = auClient.save_assembly_from_fasta({
                        'file': {'path': pos_contig_file_path},
                        'workspace_name': params['workspace_name'],
                        'assembly_name': output_obj_name
                    })
                elif input_assembly_type == ama_assembly_obj_type:
                    output_data_ref = gfuClient.fasta_gff_to_metagenome({
                        'fasta_file': {'path': pos_contig_file_path},
                        'gff_file': {'path': pos_gff_file_path},
                        'source': 'GFF',
                        'generate_missing_genes': 1,
                        'scientific_name': output_obj_name,  # not used
                        'workspace_name': params['workspace_name'],
                        'genome_name': output_obj_name
                    })['metagenome_ref']
                else:
                    raise ValueError("unknown type "+input_assembly_type)

                fractionated_obj_refs.append(output_data_ref)
                fractionated_obj_names.append(output_obj_name)

            # Neg fraction
            if not params.get('fractionate_mode') \
               or params['fractionate_mode'] == 'both' \
               or params['fractionate_mode'] == 'neg':

                if params.get('fractionate_mode') == 'neg':
                    output_obj_name = params['output_name']
                else:
                    output_obj_name = params['output_name']+'-negative_fraction'

                if input_assembly_type == assembly_obj_type:
                    output_data_ref = auClient.save_assembly_from_fasta({
                        'file': {'path': neg_contig_file_path},
                        'workspace_name': params['workspace_name'],
                        'assembly_name': output_obj_name
                    })
                elif input_assembly_type == ama_assembly_obj_type:
                    output_data_ref = gfuClient.fasta_gff_to_metagenome({
                        'fasta_file': {'path': neg_contig_file_path},
                        'gff_file': {'path': neg_gff_file_path},
                        'source': 'GFF',
                        'generate_missing_genes': 1,
                        'scientific_name': output_obj_name,  # not used
                        'workspace_name': params['workspace_name'],
                        'genome_name': output_obj_name
                    })['metagenome_ref']
                else:
                    raise ValueError("unknown type "+input_assembly_type)

                fractionated_obj_refs.append(output_data_ref)
                fractionated_obj_names.append(output_obj_name)


        #### STEP 7: generate and save the report
        ##
        if len(invalid_msgs) > 0:
            report_text += "\n".join(invalid_msgs)
            objects_created = []
        else:
            # report text
            report_text += 'Assembly '+input_assembly_name+"\n"
            report_text += "\t"+'ORIGINAL contig count: '+str(orig_contig_count)+"\n"
            report_text += "\t"+'ORIGINAL contig length sum: '+str(orig_contig_len)+"\n"
            if input_assembly_type == ama_assembly_obj_type:
                report_text += "\t"+'ORIGINAL feature count: '+str(orig_feature_count)+"\n"
            report_text += "\n"
            report_text += "\t"+'POSITIVE FRACTION contig count: '+str(pos_contig_count)+"\n"
            report_text += "\t"+'POSITIVE FRACTION contig length sum: '+str(pos_contig_len)+"\n"
            if input_assembly_type == ama_assembly_obj_type:
                report_text += "\t"+'POSITIVE FRACTION feature count: '+str(pos_feature_count)+"\n"
            report_text += "\n"
            report_text += "\t"+'NEGATIVE FRACTION contig count: '+str(neg_contig_count)+"\n"
            report_text += "\t"+'NEGATIVE FRACTION contig length sum: '+str(neg_contig_len)+"\n"
            if input_assembly_type == ama_assembly_obj_type:
                report_text += "\t"+'NEGATIVE FRACTION feature count: '+str(neg_feature_count)+"\n"

            # created objects
            objects_created = []
            for ass_i,fractionated_obj_ref in enumerate(fractionated_obj_refs):
                objects_created.append({'ref': fractionated_obj_ref, 'description': fractionated_obj_names[ass_i]+" fractionated contigs"})

        # Save report
        print('Saving report')
        kbr = KBaseReport(self.callbackURL)
        try:
            report_info = kbr.create_extended_report(
                {'message': report_text,
                 'objects_created': objects_created,
                 'direct_html_link_index': None,
                 'html_links': [],
                 'file_links': [],
                 'report_object_name': 'kb_fractionate_contigs_report_' + str(uuid.uuid4()),
                 'workspace_name': params['workspace_name']
                 })
        #except _RepError as re:
        except Exception as re:
            # not really any way to test this, all inputs have been checked earlier and should be
            # ok 
            print('Logging exception from creating report object')
            print(str(re))
            # TODO delete shock node
            raise

        # STEP 6: contruct the output to send back
        output_info = {'report_name': report_info['name'],
                       'report_ref': report_info['ref'],
                       'source_contigs_count': orig_contig_count,
                       'positive_contigs_count': pos_contig_count,
                       'negative_contigs_count': neg_contig_count,
                       'source_contigs_sum_length': orig_contig_len,
                       'positive_contigs_sum_length': pos_contig_len,
                       'negative_contigs_sum_length': neg_contig_len
        }
        if input_assembly_type == ama_assembly_obj_type:
            output_info['source_contigs_feature_count'] = orig_feature_count
            output_info['positive_contigs_feature_count'] = pos_feature_count
            output_info['negative_contigs_feature_count'] = neg_feature_count
            
        returnVal = output_info

        #END KButil_Fractionate_Reads_by_Contigs

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method KButil_Fractionate_Reads_by_Contigs return value ' +
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
