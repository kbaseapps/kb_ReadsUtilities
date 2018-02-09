
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
 * <p>Original spec-file type: KButil_Random_Subsample_Reads_Params</p>
 * 
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "workspace_name",
    "input_ref",
    "output_name",
    "subsample_fraction",
    "desc",
    "seed"
})
public class KButilRandomSubsampleReadsParams {

    @JsonProperty("workspace_name")
    private String workspaceName;
    @JsonProperty("input_ref")
    private String inputRef;
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
    @JsonProperty("desc")
    private String desc;
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

    public KButilRandomSubsampleReadsParams withWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
        return this;
    }

    @JsonProperty("input_ref")
    public String getInputRef() {
        return inputRef;
    }

    @JsonProperty("input_ref")
    public void setInputRef(String inputRef) {
        this.inputRef = inputRef;
    }

    public KButilRandomSubsampleReadsParams withInputRef(String inputRef) {
        this.inputRef = inputRef;
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

    public KButilRandomSubsampleReadsParams withOutputName(String outputName) {
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

    public KButilRandomSubsampleReadsParams withSubsampleFraction(FractionateOptions subsampleFraction) {
        this.subsampleFraction = subsampleFraction;
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

    public KButilRandomSubsampleReadsParams withDesc(String desc) {
        this.desc = desc;
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

    public KButilRandomSubsampleReadsParams withSeed(Long seed) {
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
        return ((((((((((((((("KButilRandomSubsampleReadsParams"+" [workspaceName=")+ workspaceName)+", inputRef=")+ inputRef)+", outputName=")+ outputName)+", subsampleFraction=")+ subsampleFraction)+", desc=")+ desc)+", seed=")+ seed)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
