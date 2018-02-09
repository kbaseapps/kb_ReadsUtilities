
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
 * <p>Original spec-file type: KButil_AddInsertLen_to_ReadsLibs_Params</p>
 * <pre>
 * KButil_AddInsertLen_to_ReadsLibs()
 * **
 * **  Method for Adding Insert Len to PairedEnd ReadsLibs
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "workspace_name",
    "input_refs",
    "insert_len",
    "insert_stddev"
})
public class KButilAddInsertLenToReadsLibsParams {

    @JsonProperty("workspace_name")
    private String workspaceName;
    @JsonProperty("input_refs")
    private String inputRefs;
    @JsonProperty("insert_len")
    private Double insertLen;
    @JsonProperty("insert_stddev")
    private Double insertStddev;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("workspace_name")
    public String getWorkspaceName() {
        return workspaceName;
    }

    @JsonProperty("workspace_name")
    public void setWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
    }

    public KButilAddInsertLenToReadsLibsParams withWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
        return this;
    }

    @JsonProperty("input_refs")
    public String getInputRefs() {
        return inputRefs;
    }

    @JsonProperty("input_refs")
    public void setInputRefs(String inputRefs) {
        this.inputRefs = inputRefs;
    }

    public KButilAddInsertLenToReadsLibsParams withInputRefs(String inputRefs) {
        this.inputRefs = inputRefs;
        return this;
    }

    @JsonProperty("insert_len")
    public Double getInsertLen() {
        return insertLen;
    }

    @JsonProperty("insert_len")
    public void setInsertLen(Double insertLen) {
        this.insertLen = insertLen;
    }

    public KButilAddInsertLenToReadsLibsParams withInsertLen(Double insertLen) {
        this.insertLen = insertLen;
        return this;
    }

    @JsonProperty("insert_stddev")
    public Double getInsertStddev() {
        return insertStddev;
    }

    @JsonProperty("insert_stddev")
    public void setInsertStddev(Double insertStddev) {
        this.insertStddev = insertStddev;
    }

    public KButilAddInsertLenToReadsLibsParams withInsertStddev(Double insertStddev) {
        this.insertStddev = insertStddev;
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
        return ((((((((((("KButilAddInsertLenToReadsLibsParams"+" [workspaceName=")+ workspaceName)+", inputRefs=")+ inputRefs)+", insertLen=")+ insertLen)+", insertStddev=")+ insertStddev)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
