package org.broadinstitute.gatk.tools.walkers.variantqc;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import org.broadinstitute.gatk.utils.report.GATKReportTable;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * Created by bimber on 5/12/2017.
 */
public class SectionJsonDescriptor {
    protected final String sectionLabel;
    protected final List<ReportDescriptor> rds;
    protected final List<String> stratifications;

    protected enum PlotType {
        data_table(),
        bar_graph(),
        xy_line();
    };

    public SectionJsonDescriptor(String sectionLabel, String[] stratifications){
        this.sectionLabel = sectionLabel;
        this.rds = new ArrayList<>();
        this.stratifications = Arrays.asList(stratifications);
    }

    public void addReportDescriptor(ReportDescriptor rd, GATKReportTable table, Map<String, String> descriptionMap){
        rd.bindSection(this, table, descriptionMap);
        rds.add(rd);
    }

    public JsonObject getConfig(){
        JsonObject ret = new JsonObject();
        ret.addProperty("label", sectionLabel);

        JsonArray reports = new JsonArray();
        for (ReportDescriptor rd : rds){
            reports.add(rd.getReportJson(sectionLabel));
        }
        ret.add("reports", reports);
        
        return ret;
    }

    //    public static class XYLineReportDescriptor extends ReportDescriptor
//    {
//        private String[] columnsToPlot;
//        private String yLabel;
//
//        public XYLineReportDescriptor(String plotTitle, SectionJsonDescriptor.PlotType plotType, String evaluatorModuleName){
//            super(plotTitle, plotType, evaluatorModuleName);
//        }
//
//        @Override
//        void JsonObject getReportJson(String sectionTitle, GATKReportTable table) {
//            JsonArray samples = inferSampleNames(table);
//            ret.add("samples", samples);
//
//            JsonArray datasetsJson = new JsonArray();
//            JsonObject colorJson = new JsonObject();
//            for (String colName : columnsToPlot){
//                int colIdx = getColumnByName(table, colName);
//
//                JsonObject datasetJson = new JsonObject();
//                //datasetJson.addProperty("color", "");
//                datasetJson.addProperty("name", colName);
//
//                JsonArray data = new JsonArray();
//                for (int i=0;i<table.getNumRows();i++){
//                    data.add(new JsonPrimitive(NumberUtils.createNumber(table.get(i, colIdx).toString())));
//                }
//                datasetJson.add("data", data);
//            }
//
//
//
//            //colorJson.addProperty(sample, null);
//
//            //"data": [[1, 32.385829384477994], [2, 32.51257319784182], [3, 32.53051509589336], [4, 35.930739475597946], [5, 35.91741590046593], [6, 35.928392537730666], [7, 35.91026610496744], [8, 35.897712428885164], [9, 37.57662158757583], [10, 37.585727700319836], [12, 37.53505206612946], [14, 38.92510002043037], [16, 38.86395021022102], [18, 38.8059641509619], [20, 38.73872162799532], [22, 38.63794941362771], [24, 38.520917890439165], [26, 38.36361591192953], [28, 38.230401956153074], [30, 38.06566984593516], [32, 37.91252182451145], [34, 37.841280868131435], [36, 37.878115239254726], [38, 37.8159351759587], [40, 37.661726689574806], [42, 37.51950302387578], [44, 37.33677994331231], [46, 37.13236408525927], [48, 36.92583309222985], [50, 36.258499335319016], [52, 36.14480792257834], [54, 36.53248978790511], [56, 36.499392821520985], [58, 36.315043984495425], [60, 36.074915082678686], [62, 35.79003837432995], [64, 35.528671738335916], [66, 35.24913007505768], [68, 34.97733730062913], [70, 34.72833874478797], [72, 34.4612219301588], [74, 34.21380766485854], [76, 33.99092253379628], [78, 33.77986131866412], [80, 33.591852978159906], [82, 33.43427846620381], [84, 33.291657585804636], [86, 33.19623853098773], [88, 33.10817420745926], [90, 33.005813719535894], [92, 32.89738799884225], [94, 32.93315393977193], [96, 33.019364145676036], [98, 33.18226949894594], [100, 32.70148968981316]],
//
//            JsonObject configJson = new JsonObject();
//            configJson.addProperty("xlab", this.yLabel);
//
//            configJson.add("colors", colorJson);
//            configJson.addProperty("tt_label", "<b>Base {point.x}</b>: {point.y:.2f}");
//            configJson.addProperty("xDecimals", false);
//            configJson.addProperty("title", label);
//            configJson.addProperty("ylab", this.yLabel);
//
//            ret.add("config", configJson);
//        }
//    }
}
