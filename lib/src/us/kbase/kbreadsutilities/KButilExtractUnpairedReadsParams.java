
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
 * <p>Original spec-file type: KButil_Extract_Unpaired_Reads_Params</p>
 * <pre>
 * KButil_Extract_Unpaired_Reads_and_Synchronize_Pairs()
 * **
 * **  Method for removing unpaired reads from a paired end library or set and matching the order of reads
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "workspace_name",
    "input_ref",
    "output_name",
    "desc"
})
public class KButilExtractUnpairedReadsParams {

    @JsonProperty("workspace_name")
    private String workspaceName;
    @JsonProperty("input_ref")
    private String inputRef;
    @JsonProperty("output_name")
    private String outputName;
    @JsonProperty("desc")
    private String desc;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("workspace_name")
    public String getWorkspaceName() {
        return workspaceName;
    }

    @JsonProperty("workspace_name")
    public void setWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
    }

    public KButilExtractUnpairedReadsParams withWorkspaceName(String workspaceName) {
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

    public KButilExtractUnpairedReadsParams withInputRef(String inputRef) {
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

    public KButilExtractUnpairedReadsParams withOutputName(String outputName) {
        this.outputName = outputName;
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

    public KButilExtractUnpairedReadsParams withDesc(String desc) {
        this.desc = desc;
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
        return ((((((((((("KButilExtractUnpairedReadsParams"+" [workspaceName=")+ workspaceName)+", inputRef=")+ inputRef)+", outputName=")+ outputName)+", desc=")+ desc)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
