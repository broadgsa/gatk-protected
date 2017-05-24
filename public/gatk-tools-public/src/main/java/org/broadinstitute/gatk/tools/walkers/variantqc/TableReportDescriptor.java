package org.broadinstitute.gatk.tools.walkers.variantqc;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import org.broadinstitute.gatk.utils.report.GATKReportColumn;

import java.util.Map;

/**
 * Created by bimber on 5/22/2017.
 */
public class TableReportDescriptor extends ReportDescriptor {
    public TableReportDescriptor(String label, SectionJsonDescriptor.PlotType plotType, String evaluatorModuleName) {
        super(label, plotType, evaluatorModuleName);
    }

    public static TableReportDescriptor getCountVariantsTable() {
        return new TableReportDescriptor("Variant Summary", SectionJsonDescriptor.PlotType.data_table, "CountVariants");
    }

    @Override
    public JsonObject getReportJson(String sectionTitle) {
        JsonObject ret = new JsonObject();
        ret.addProperty("label", label);

        JsonObject data = new JsonObject();
        ret.add("data", data);

        data.addProperty("plot_type", plotType.name());
        data.add("samples", getSampleNames());
        data.add("datasets", new JsonArray());

        data.add("columns", new JsonArray());
        for (GATKReportColumn col : table.getColumnInfo()) {
            JsonObject colJson = new JsonObject();
            colJson.addProperty("name", col.getColumnName());
            colJson.addProperty("label", col.getColumnName());

            //we will probably need a way to provide information beyond just what is in the GATKReportTable itself
            if (columnInfoMap.containsKey(col.getColumnName())) {
                for (Map.Entry<String, JsonElement> entry : columnInfoMap.get(col.getColumnName()).entrySet()) {
                    colJson.add(entry.getKey(), entry.getValue());
                }
            }


            if (!colJson.has("dmin")) {
                //TODO: infer based on data.  perhaps find the min/max and take +/- 10%?
            }

            data.getAsJsonArray("columns").add(colJson);
        }

        return ret;
    }
}
