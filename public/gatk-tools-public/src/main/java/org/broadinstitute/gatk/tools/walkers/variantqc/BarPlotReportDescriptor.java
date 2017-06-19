package org.broadinstitute.gatk.tools.walkers.variantqc;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import com.google.gson.JsonPrimitive;
import org.apache.commons.lang.math.NumberUtils;
import org.broadinstitute.gatk.utils.report.GATKReportColumn;
import org.broadinstitute.gatk.utils.report.GATKReportTable;

import java.util.*;


/**
 * Created by bimber on 5/22/2017.
 */
public class BarPlotReportDescriptor extends ReportDescriptor {
    private String[] columnsToPlot;
    private String yLabel;

    public BarPlotReportDescriptor(String plotTitle, SectionJsonDescriptor.PlotType plotType, String evaluatorModuleName, String[] columnsToPlot, String yLabel, Collection<String> skippedSamples) {
        super(plotTitle, plotType, evaluatorModuleName);
        this.columnsToPlot = columnsToPlot;
        this.yLabel = yLabel;
        if (skippedSamples != null){
            this.skippedSamples.addAll(skippedSamples);
        }
    }

    public List<String> getColumnsToPlot(GATKReportTable table){
        return Arrays.asList(columnsToPlot);
    }

    public static BarPlotReportDescriptor getVariantTypeBarPlot() {
        return new BarPlotReportDescriptor("Variant Type", SectionJsonDescriptor.PlotType.bar_graph, "CountVariants", new String[]{"nSNPs", "nMNPs", "nInsertions", "nDeletions", "nComplex", "nSymbolic", "nMixed"}, "# Variants", Arrays.asList("all"));
    }

    public static BarPlotReportDescriptor getSiteFilterTypeBarPlot() {
        return new BarPlotReportDescriptor("Filter Type", SectionJsonDescriptor.PlotType.bar_graph, "CountVariants", null, "# Variants", Arrays.asList("all")){
            @Override
            public List<String> getColumnsToPlot(GATKReportTable table){
                List<String> ret = new ArrayList<>();
                for (GATKReportColumn col : table.getColumnInfo()){
                    if (!sectionConfig.stratifications.contains(col.getColumnName())){
                        ret.add(col.getColumnName());
                    }
                }

                return ret;
            }
        };
    }

    @Override
    public JsonObject getReportJson(String sectionTitle) {
        JsonObject ret = new JsonObject();
        ret.addProperty("label", label);

        JsonObject dataObj = new JsonObject();
        ret.add("data", dataObj);

        dataObj.addProperty("plot_type", plotType.name());

        dataObj.add("samples", new JsonArray());
        dataObj.getAsJsonArray("samples").add(getSampleNames());

        JsonArray datasetsJson = new JsonArray();

        for (String colName : getColumnsToPlot(table)) {
            JsonObject datasetJson = new JsonObject();
            datasetJson.addProperty("name", colName);

            JsonArray data = new JsonArray();
            for (Object rowId : table.getRowIDs()) {
                String sampleName = getSampleNameForRow(rowId);
                if (skippedSamples.contains(sampleName)){
                    continue;
                }

                data.add(new JsonPrimitive(NumberUtils.createNumber(table.get(rowId, colName).toString())));
            }
            datasetJson.add("data", data);

            datasetsJson.add(datasetJson);
        }
        dataObj.add("datasets", new JsonArray());
        dataObj.getAsJsonArray("datasets").add(datasetsJson);

        JsonObject configJson = new JsonObject();
        configJson.addProperty("ylab", this.yLabel);
        configJson.addProperty("title", label);

        dataObj.add("config", configJson);

        return ret;
    }
}
