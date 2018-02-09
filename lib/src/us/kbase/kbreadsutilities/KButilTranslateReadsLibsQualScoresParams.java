
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
 * <p>Original spec-file type: KButil_Translate_ReadsLibs_QualScores_Params</p>
 * <pre>
 * KButil_Translate_ReadsLibs_QualScores()
 * **
 * **  Method for Translating ReadsLibs Qual scores
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "workspace_name",
    "input_refs"
})
public class KButilTranslateReadsLibsQualScoresParams {

    @JsonProperty("workspace_name")
    private String workspaceName;
    @JsonProperty("input_refs")
    private String inputRefs;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("workspace_name")
    public String getWorkspaceName() {
        return workspaceName;
    }

    @JsonProperty("workspace_name")
    public void setWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
    }

    public KButilTranslateReadsLibsQualScoresParams withWorkspaceName(String workspaceName) {
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

    public KButilTranslateReadsLibsQualScoresParams withInputRefs(String inputRefs) {
        this.inputRefs = inputRefs;
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
        return ((((((("KButilTranslateReadsLibsQualScoresParams"+" [workspaceName=")+ workspaceName)+", inputRefs=")+ inputRefs)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
