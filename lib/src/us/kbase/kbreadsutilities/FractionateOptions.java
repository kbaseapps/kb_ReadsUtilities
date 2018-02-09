
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
 * <p>Original spec-file type: Fractionate_Options</p>
 * <pre>
 * KButil_Random_Subsample_Reads()
 * **
 * **  Method for random subsampling of reads library
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "split_num",
    "reads_num",
    "reads_perc"
})
public class FractionateOptions {

    @JsonProperty("split_num")
    private Long splitNum;
    @JsonProperty("reads_num")
    private Long readsNum;
    @JsonProperty("reads_perc")
    private Double readsPerc;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("split_num")
    public Long getSplitNum() {
        return splitNum;
    }

    @JsonProperty("split_num")
    public void setSplitNum(Long splitNum) {
        this.splitNum = splitNum;
    }

    public FractionateOptions withSplitNum(Long splitNum) {
        this.splitNum = splitNum;
        return this;
    }

    @JsonProperty("reads_num")
    public Long getReadsNum() {
        return readsNum;
    }

    @JsonProperty("reads_num")
    public void setReadsNum(Long readsNum) {
        this.readsNum = readsNum;
    }

    public FractionateOptions withReadsNum(Long readsNum) {
        this.readsNum = readsNum;
        return this;
    }

    @JsonProperty("reads_perc")
    public Double getReadsPerc() {
        return readsPerc;
    }

    @JsonProperty("reads_perc")
    public void setReadsPerc(Double readsPerc) {
        this.readsPerc = readsPerc;
    }

    public FractionateOptions withReadsPerc(Double readsPerc) {
        this.readsPerc = readsPerc;
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
        return ((((((((("FractionateOptions"+" [splitNum=")+ splitNum)+", readsNum=")+ readsNum)+", readsPerc=")+ readsPerc)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
