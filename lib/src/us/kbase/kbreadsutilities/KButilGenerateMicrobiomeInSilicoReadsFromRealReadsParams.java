
package us.kbase.kbreadsutilities;

import java.util.HashMap;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: KButil_Generate_Microbiome_InSilico_Reads_From_Real_Reads_Params</p>
 * <pre>
 * KButil_Generate_Microbiome_InSilico_Reads_From_Real_Reads()
 * **
 * <    **  Method for random subsampling of reads library combined with overlay of configured genomes
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "workspace_name",
    "input_genomeSet_ref",
    "genome_abundances",
    "input_reads_ref",
    "output_name",
    "subsample_fraction",
    "genome_length_bias",
    "desc",
    "pe_insert_len",
    "pe_orientation",
    "seed"
})
public class KButilGenerateMicrobiomeInSilicoReadsFromRealReadsParams {

    @JsonProperty("workspace_name")
    private String workspaceName;
    @JsonProperty("input_genomeSet_ref")
    private String inputGenomeSetRef;
    @JsonProperty("genome_abundances")
    private String genomeAbundances;
    @JsonProperty("input_reads_ref")
    private String inputReadsRef;
    @JsonProperty("output_name")
    private String outputName;
    /**
     * <p>Original spec-file type: Fractionate_Options</p>
     * <pre>
     * KButil_Random_Subsample_Reads()
     * **
     * **  Method for random subsampling of reads library
     * </pre>
     * 
     */
    @JsonProperty("subsample_fraction")
    private FractionateOptions subsampleFraction;
    @JsonProperty("genome_length_bias")
    private Long genomeLengthBias;
    @JsonProperty("desc")
    private String desc;
    @JsonProperty("pe_insert_len")
    private Long peInsertLen;
    @JsonProperty("pe_orientation")
    private String peOrientation;
    @JsonProperty("seed")
    private Long seed;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("workspace_name")
    public String getWorkspaceName() {
        return workspaceName;
    }

    @JsonProperty("workspace_name")
    public void setWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
    }

    public KButilGenerateMicrobiomeInSilicoReadsFromRealReadsParams withWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
        return this;
    }

    @JsonProperty("input_genomeSet_ref")
    public String getInputGenomeSetRef() {
        return inputGenomeSetRef;
    }

    @JsonProperty("input_genomeSet_ref")
    public void setInputGenomeSetRef(String inputGenomeSetRef) {
        this.inputGenomeSetRef = inputGenomeSetRef;
    }

    public KButilGenerateMicrobiomeInSilicoReadsFromRealReadsParams withInputGenomeSetRef(String inputGenomeSetRef) {
        this.inputGenomeSetRef = inputGenomeSetRef;
        return this;
    }

    @JsonProperty("genome_abundances")
    public String getGenomeAbundances() {
        return genomeAbundances;
    }

    @JsonProperty("genome_abundances")
    public void setGenomeAbundances(String genomeAbundances) {
        this.genomeAbundances = genomeAbundances;
    }

    public KButilGenerateMicrobiomeInSilicoReadsFromRealReadsParams withGenomeAbundances(String genomeAbundances) {
        this.genomeAbundances = genomeAbundances;
        return this;
    }

    @JsonProperty("input_reads_ref")
    public String getInputReadsRef() {
        return inputReadsRef;
    }

    @JsonProperty("input_reads_ref")
    public void setInputReadsRef(String inputReadsRef) {
        this.inputReadsRef = inputReadsRef;
    }

    public KButilGenerateMicrobiomeInSilicoReadsFromRealReadsParams withInputReadsRef(String inputReadsRef) {
        this.inputReadsRef = inputReadsRef;
        return this;
    }

    @JsonProperty("output_name")
    public String getOutputName() {
        return outputName;
    }

    @JsonProperty("output_name")
    public void setOutputName(String outputName) {
        this.outputName = outputName;
    }

    public KButilGenerateMicrobiomeInSilicoReadsFromRealReadsParams withOutputName(String outputName) {
        this.outputName = outputName;
        return this;
    }

    /**
     * <p>Original spec-file type: Fractionate_Options</p>
     * <pre>
     * KButil_Random_Subsample_Reads()
     * **
     * **  Method for random subsampling of reads library
     * </pre>
     * 
     */
    @JsonProperty("subsample_fraction")
    public FractionateOptions getSubsampleFraction() {
        return subsampleFraction;
    }

    /**
     * <p>Original spec-file type: Fractionate_Options</p>
     * <pre>
     * KButil_Random_Subsample_Reads()
     * **
     * **  Method for random subsampling of reads library
     * </pre>
     * 
     */
    @JsonProperty("subsample_fraction")
    public void setSubsampleFraction(FractionateOptions subsampleFraction) {
        this.subsampleFraction = subsampleFraction;
    }

    public KButilGenerateMicrobiomeInSilicoReadsFromRealReadsParams withSubsampleFraction(FractionateOptions subsampleFraction) {
        this.subsampleFraction = subsampleFraction;
        return this;
    }

    @JsonProperty("genome_length_bias")
    public Long getGenomeLengthBias() {
        return genomeLengthBias;
    }

    @JsonProperty("genome_length_bias")
    public void setGenomeLengthBias(Long genomeLengthBias) {
        this.genomeLengthBias = genomeLengthBias;
    }

    public KButilGenerateMicrobiomeInSilicoReadsFromRealReadsParams withGenomeLengthBias(Long genomeLengthBias) {
        this.genomeLengthBias = genomeLengthBias;
        return this;
    }

    @JsonProperty("desc")
    public String getDesc() {
        return desc;
    }

    @JsonProperty("desc")
    public void setDesc(String desc) {
        this.desc = desc;
    }

    public KButilGenerateMicrobiomeInSilicoReadsFromRealReadsParams withDesc(String desc) {
        this.desc = desc;
        return this;
    }

    @JsonProperty("pe_insert_len")
    public Long getPeInsertLen() {
        return peInsertLen;
    }

    @JsonProperty("pe_insert_len")
    public void setPeInsertLen(Long peInsertLen) {
        this.peInsertLen = peInsertLen;
    }

    public KButilGenerateMicrobiomeInSilicoReadsFromRealReadsParams withPeInsertLen(Long peInsertLen) {
        this.peInsertLen = peInsertLen;
        return this;
    }

    @JsonProperty("pe_orientation")
    public String getPeOrientation() {
        return peOrientation;
    }

    @JsonProperty("pe_orientation")
    public void setPeOrientation(String peOrientation) {
        this.peOrientation = peOrientation;
    }

    public KButilGenerateMicrobiomeInSilicoReadsFromRealReadsParams withPeOrientation(String peOrientation) {
        this.peOrientation = peOrientation;
        return this;
    }

    @JsonProperty("seed")
    public Long getSeed() {
        return seed;
    }

    @JsonProperty("seed")
    public void setSeed(Long seed) {
        this.seed = seed;
    }

    public KButilGenerateMicrobiomeInSilicoReadsFromRealReadsParams withSeed(Long seed) {
        this.seed = seed;
        return this;
    }

    @JsonAnyGetter
    public Map<String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public String toString() {
        return ((((((((((((((((((((((((("KButilGenerateMicrobiomeInSilicoReadsFromRealReadsParams"+" [workspaceName=")+ workspaceName)+", inputGenomeSetRef=")+ inputGenomeSetRef)+", genomeAbundances=")+ genomeAbundances)+", inputReadsRef=")+ inputReadsRef)+", outputName=")+ outputName)+", subsampleFraction=")+ subsampleFraction)+", genomeLengthBias=")+ genomeLengthBias)+", desc=")+ desc)+", peInsertLen=")+ peInsertLen)+", peOrientation=")+ peOrientation)+", seed=")+ seed)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
