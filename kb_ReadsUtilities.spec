/*
** A KBase module: kb_ReadsUtilities
**
** This module contains basic utility Apps for manipulating Reads Libraries
*/

module kb_ReadsUtilities {

    /* 
    ** The workspace object refs are of form:
    **
    **    objects = ws.get_objects([{'ref': params['workspace_id']+'/'+params['obj_name']}])
    **
    ** "ref" means the entire name combining the workspace id and the object name
    ** "id" is a numerical identifier of the workspace or object, and should just be used for workspace
    ** "name" is a string identifier of a workspace or object.  This is received from Narrative.
    */
    typedef string workspace_name;
    typedef string sequence;
    typedef string data_obj_name;
    typedef string data_obj_ref;
    typedef int    bool;


    /* KButil_FASTQ_to_FASTA()
    **
    ** Method for Converting a FASTQ SingleEndLibrary to a FASTA SingleEndLibrary
    */
    typedef structure {
        workspace_name workspace_name;
	data_obj_ref   input_ref;
        data_obj_name  output_name;
    } KButil_FASTQ_to_FASTA_Params;

    typedef structure {
	data_obj_name report_name;
	data_obj_ref  report_ref;
/*       data_obj_ref  output_filtered_ref;
*
*        int n_initial_seqs;
*        int n_seqs_matched;
*        int n_seqs_notmatched;
*/
    } KButil_FASTQ_to_FASTA_Output;
	
    funcdef KButil_FASTQ_to_FASTA (KButil_FASTQ_to_FASTA_Params params)  returns (KButil_FASTQ_to_FASTA_Output) authentication required;


    /* KButil_Split_Reads()
    **
    **  Method for spliting a ReadsLibrary into evenly sized ReadsLibraries
    */
    typedef structure {
        workspace_name workspace_name;
	data_obj_ref   input_ref;    /* ReadsLibrary */
        data_obj_name  output_name;  /* ReadsSet */
	int            split_num;
	string         desc;
    } KButil_Split_Reads_Params;

    typedef structure {
	data_obj_name report_name;
	data_obj_ref  report_ref;
    } KButil_Split_Reads_Output;

    funcdef KButil_Split_Reads (KButil_Split_Reads_Params params)  returns (KButil_Split_Reads_Output) authentication required;


    /* KButil_Random_Subsample_Reads()
    **
    **  Method for random subsampling of reads library
    */
    typedef structure {
	int            split_num;
	int            reads_num;
	float          reads_perc;
    } Fractionate_Options;
    typedef structure {
        workspace_name workspace_name;
	data_obj_ref   input_ref;    /* ReadsLibrary */
        data_obj_name  output_name;  /* ReadsSet */
	Fractionate_Options subsample_fraction;
	/*bool           reads_uniq;*/  /* sampling without replacement */
	string         desc;
	int            seed;
    } KButil_Random_Subsample_Reads_Params;

    typedef structure {
	data_obj_name report_name;
	data_obj_ref  report_ref;
    } KButil_Random_Subsample_Reads_Output;

    funcdef KButil_Random_Subsample_Reads (KButil_Random_Subsample_Reads_Params params)  returns (KButil_Random_Subsample_Reads_Output) authentication required;


    /* KButil_Merge_ReadsSet_to_OneLibrary()
    **
    **  Method for merging a ReadsSet into one library
    */
    typedef structure {
        workspace_name workspace_name;
	data_obj_ref   input_ref;    /* ReadsSet */
        data_obj_name  output_name;  /* ReadsLibrary */
	string         desc;
    } KButil_Merge_ReadsSet_to_OneLibrary_Params;

    typedef structure {
	data_obj_name report_name;
	data_obj_ref  report_ref;
    } KButil_Merge_ReadsSet_to_OneLibrary_Output;

    funcdef KButil_Merge_ReadsSet_to_OneLibrary (KButil_Merge_ReadsSet_to_OneLibrary_Params params)  returns (KButil_Merge_ReadsSet_to_OneLibrary_Output) authentication required;


    /* KButil_Merge_MultipleReadsLibs_to_OneLibrary()
    **
    **  Method for merging ReadsLibs into one library
    */
    typedef structure {
        workspace_name workspace_name;
	data_obj_ref   input_refs;    /* ReadsLibraries */
        data_obj_name  output_name;  /* ReadsLibrary */
	string         desc;
    } KButil_Merge_MultipleReadsLibs_to_OneLibrary_Params;

    typedef structure {
	data_obj_name report_name;
	data_obj_ref  report_ref;
    } KButil_Merge_MultipleReadsLibs_to_OneLibrary_Output;

    funcdef KButil_Merge_MultipleReadsLibs_to_OneLibrary (KButil_Merge_MultipleReadsLibs_to_OneLibrary_Params params)  returns (KButil_Merge_MultipleReadsLibs_to_OneLibrary_Output) authentication required;


    /* KButil_Extract_Unpaired_Reads_and_Synchronize_Pairs()
    **
    **  Method for removing unpaired reads from a paired end library or set and matching the order of reads
    */
    typedef structure {
        workspace_name workspace_name;
	data_obj_ref   input_ref;    /* ReadsSet or ReadLibrary */
        data_obj_name  output_name;  /* ReadsSet or ReadLibrary */
	string         desc;
    } KButil_Extract_Unpaired_Reads_Params;

    typedef structure {
	data_obj_name report_name;
	data_obj_ref  report_ref;
    } KButil_Extract_Unpaired_Reads_Output;

    funcdef KButil_Extract_Unpaired_Reads_and_Synchronize_Pairs (KButil_Extract_Unpaired_Reads_Params params)  returns (KButil_Extract_Unpaired_Reads_Output) authentication required;


    /* KButil_Translate_ReadsLibs_QualScores()
    **
    **  Method for Translating ReadsLibs Qual scores
    */
    typedef structure {
        workspace_name workspace_name;
	data_obj_ref   input_refs;    /* ReadsLibraries */
    } KButil_Translate_ReadsLibs_QualScores_Params;

    typedef structure {
	data_obj_name report_name;
	data_obj_ref  report_ref;
    } KButil_Translate_ReadsLibs_QualScores_Output;

    funcdef KButil_Translate_ReadsLibs_QualScores (KButil_Translate_ReadsLibs_QualScores_Params params)  returns (KButil_Translate_ReadsLibs_QualScores_Output) authentication required;


    /* KButil_AddInsertLen_to_ReadsLibs()
    **
    **  Method for Adding Insert Len to PairedEnd ReadsLibs
    */
    typedef structure {
        workspace_name workspace_name;
	data_obj_ref   input_refs;    /* ReadsLibraries */
	float          insert_len;
	float          insert_stddev;
    } KButil_AddInsertLen_to_ReadsLibs_Params;

    typedef structure {
	data_obj_name report_name;
	data_obj_ref  report_ref;
    } KButil_AddInsertLen_to_ReadsLibs_Output;

    funcdef KButil_AddInsertLen_to_ReadsLibs (KButil_AddInsertLen_to_ReadsLibs_Params params)  returns (KButil_AddInsertLen_to_ReadsLibs_Output) authentication required;


    /* KButil_Generate_Microbiome_InSilico_Reads_From_Real_Reads()
    **
<    **  Method for random subsampling of reads library combined with overlay of configured genomes
    */
    typedef structure {
        workspace_name workspace_name;
	data_obj_ref   input_genomeSet_ref;  /* GenomeSet */
	string         genome_abundances;    /* Genome abundances */
	data_obj_ref   input_reads_ref;      /* ReadsLibrary */
        data_obj_name  output_name;          /* ReadsSet */
	Fractionate_Options subsample_fraction;
	/*bool           reads_uniq;*/  /* sampling without replacement */
	bool           genome_length_bias;
	string         desc;
	int            pe_insert_len;
        string         pe_orientation;
	int            seed;
    } KButil_Generate_Microbiome_InSilico_Reads_From_Real_Reads_Params;

    typedef structure {
	data_obj_name report_name;
	data_obj_ref  report_ref;
    } KButil_Generate_Microbiome_InSilico_Reads_From_Real_Reads_Output;

    funcdef KButil_Generate_Microbiome_InSilico_Reads_From_Real_Reads (KButil_Generate_Microbiome_InSilico_Reads_From_Real_Reads_Params params)  returns (KButil_Generate_Microbiome_InSilico_Reads_From_Real_Reads_Output) authentication required;

    
    /* KButil_Fractionate_Reads_by_Contigs()
    **
    **  Split reads library into two based on whether they match contigs
    */
    typedef structure {
        workspace_name workspace_name;
	data_obj_ref   input_reads_ref;     /* ReadsLibrary */
	data_obj_ref   input_assembly_ref;  /* Assembly, AMA, or Genome */
        data_obj_name  output_name;         /* BaseName for Reads */
	string         fractionate_mode;
	/*int            seed;*/
    } KButil_Fractionate_Reads_by_Contigs_Params;

    typedef structure {
	data_obj_name report_name;
	data_obj_ref  report_ref;
        int    source_reads_count;
        int    positive_reads_count;
        int    negative_reads_count;
        int    source_reads_sum_length;
        int    positive_reads_sum_length;
        int    negative_reads_sum_length;
    } KButil_Fractionate_Reads_by_Contigs_Output;

    funcdef KButil_Fractionate_Reads_by_Contigs (KButil_Fractionate_Reads_by_Contigs_Params params)  returns (KButil_Fractionate_Reads_by_Contigs_Output) authentication required;


};
