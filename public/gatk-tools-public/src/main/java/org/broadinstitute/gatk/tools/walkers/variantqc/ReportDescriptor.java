package org.broadinstitute.gatk.tools.walkers.variantqc;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import com.google.gson.JsonPrimitive;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.gatk.utils.exceptions.GATKException;
import org.broadinstitute.gatk.utils.report.GATKReportColumn;
import org.broadinstitute.gatk.utils.report.GATKReportTable;

import java.util.*;

/**
 * Created by bimber on 5/22/2017.
 */
abstract class ReportDescriptor {
    protected final String label;
    protected final SectionJsonDescriptor.PlotType plotType;
    protected final String evaluatorModuleName;
    protected Map<String, JsonObject> columnInfoMap;
    protected SectionJsonDescriptor sectionConfig;
    protected GATKReportTable table;
    protected Set<String> skippedSamples = new HashSet<>();
    protected Map<String, String> descriptionMap;

    protected ReportDescriptor(String label, SectionJsonDescriptor.PlotType plotType, String evaluatorModuleName) {
        this.label = label;
        this.plotType = plotType;
        this.evaluatorModuleName = evaluatorModuleName;
        this.columnInfoMap = new HashMap<>();
    }

    public void addColumnInfo(String colName, JsonObject columnInfo) {
        this.columnInfoMap.put(colName, columnInfo);
    }

    public void bindSection(SectionJsonDescriptor sectionConfig, GATKReportTable table, Map<String, String> descriptionMap){
        this.sectionConfig = sectionConfig;
        this.table = table;
        this.descriptionMap = descriptionMap;
    }

    abstract JsonObject getReportJson(String sectionTitle);

    protected List<String> columnsInSampleName = null;
    /**
     * This provides the opportunity for subclasses to supply custom logic to parse the sampleName from rows.
     * This can return null, in which case that row will be ignored, for example if specific states should not be included.
     * @param rowId
     * @return
     */
    protected String getSampleNameForRow(Object rowId){
        if (columnsInSampleName == null){
            List<String> ret = new ArrayList<>();
            for (GATKReportColumn col : table.getColumnInfo()){
                if (sectionConfig.stratifications.contains(col.getColumnName())){
                    ret.add(col.getColumnName());
                }
            }

            columnsInSampleName = ret;
        }

        List<String> tokens = new ArrayList<>();
        for (String colName : columnsInSampleName){
            tokens.add(table.get(rowId, colName).toString());
        }

        return StringUtil.join(" / ", tokens);
    }

    protected JsonArray getSampleNames(){
        Set<String> sampleNames = new LinkedHashSet<>();
        for (Object rowId : table.getRowIDs()){
            String sn = getSampleNameForRow(rowId);
            if (sn != null && !skippedSamples.contains(sn)){
                sampleNames.add(sn);
            }
        }

        JsonArray ret = new JsonArray();
        for (String sn : sampleNames){
            ret.add(new JsonPrimitive(sn));
        }

        return ret;
    }
}
