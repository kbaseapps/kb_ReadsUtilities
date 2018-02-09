# -*- coding: utf-8 -*-
#BEGIN_HEADER
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
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/dcchivian/kb_ReadsUtilities"
    GIT_COMMIT_HASH = "b59b0faeb860c56b556d8ee7417cdb3930623837"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
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
