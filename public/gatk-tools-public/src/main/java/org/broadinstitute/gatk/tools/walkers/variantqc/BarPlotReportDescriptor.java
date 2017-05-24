package org.broadinstitute.gatk.tools.walkers.variantqc;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import com.google.gson.JsonPrimitive;
import org.apache.commons.lang.math.NumberUtils;

/**
 * Created by bimber on 5/22/2017.
 */
public class BarPlotReportDescriptor extends ReportDescriptor {
    private String[] columnsToPlot;
    private String yLabel;

    public BarPlotReportDescriptor(String plotTitle, SectionJsonDescriptor.PlotType plotType, String evaluatorModuleName, String[] columnsToPlot, String yLabel) {
        super(plotTitle, plotType, evaluatorModuleName);
        this.columnsToPlot = columnsToPlot;
        this.yLabel = yLabel;
    }

    public static BarPlotReportDescriptor getVariantTypeBarPlot() {
        return new BarPlotReportDescriptor("Variant Type", SectionJsonDescriptor.PlotType.bar_graph, "CountVariants", new String[]{"nSNPs", "nMNPs", "nInsertions", "nDeletions", "nComplex", "nSymbolic", "nMixed"}, "# Variants");
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
        for (String colName : columnsToPlot) {
            int colIdx = getColumnByName(colName);

            JsonObject datasetJson = new JsonObject();
            datasetJson.addProperty("name", colName);

            JsonArray data = new JsonArray();
            for (int i = 0; i < table.getNumRows(); i++) {
                data.add(new JsonPrimitive(NumberUtils.createNumber(table.get(i, colIdx).toString())));
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
