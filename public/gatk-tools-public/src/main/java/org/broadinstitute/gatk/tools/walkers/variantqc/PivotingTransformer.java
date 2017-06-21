package org.broadinstitute.gatk.tools.walkers.variantqc;

import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.gatk.engine.samples.Sample;
import org.broadinstitute.gatk.engine.samples.SampleDB;
import org.broadinstitute.gatk.utils.exceptions.GATKException;
import org.broadinstitute.gatk.utils.report.GATKReportColumn;
import org.broadinstitute.gatk.utils.report.GATKReportDataType;
import org.broadinstitute.gatk.utils.report.GATKReportTable;

import java.util.*;

/**
 * Created by bimber on 5/31/2017.
 */
public class PivotingTransformer implements GATKReportTableTransformer {
    private final List<String> groupBy;
    private final List<Pivot> colsToPivot;
    private final String evalModuleName;
    private final SampleDB sampleDB;

    public PivotingTransformer(String evalModuleName, List<String> groupBy, List<Pivot> colsToPivot){
        this(evalModuleName, groupBy, colsToPivot, null);
    }

    public PivotingTransformer(String evalModuleName, List<String> groupBy, List<Pivot> colsToPivot, SampleDB sampleDB){
        this.evalModuleName = evalModuleName;
        this.groupBy = groupBy;
        this.colsToPivot = colsToPivot;
        this.sampleDB = sampleDB;
    }

    @Override
    public GATKReportTable transform(GATKReportTable table) {
        Map<String, Map<String, Object>> rowMap = new LinkedHashMap<>();
        Set<String> distinctColNames = new LinkedHashSet<>();
        Map<String, String> colFormatMap = new HashMap<>();

        distinctColNames.addAll(groupBy);
        for (String colName : groupBy){
            colFormatMap.put(colName, GATKReportDataType.String.getDefaultFormatString());
        }

        if (sampleDB != null){
            distinctColNames.add("Gender");
            colFormatMap.put("Gender", GATKReportDataType.String.getDefaultFormatString());
        }

        Map<String, GATKReportColumn> colMap = new HashMap<>();
        for (GATKReportColumn col : table.getColumnInfo()){
            colMap.put(col.getColumnName(), col);
        }

        for (int i = 0;i<table.getNumRows();i++){
            List<String> keys = new ArrayList<>();
            for (String colName : groupBy){
                keys.add(String.valueOf(table.get(i, colName)));
            }

            String key = StringUtil.join(";", keys);
            if (!rowMap.containsKey(key)){
                rowMap.put(key, new HashMap<>());
                for (String colName : groupBy){
                    rowMap.get(key).put(colName, table.get(i, colName));
                }
            }

            for (Pivot p : colsToPivot){
                if (table.get(i, p.colNameSource) == null){
                    throw new GATKException("row lacks a value for: " + p.colNameSource);
                }

                String colName = p.getTargetColName(String.valueOf(table.get(i, p.colNameSource)));
                if (!distinctColNames.contains(colName)){
                    distinctColNames.add(colName);
                    colFormatMap.put(colName, colMap.get(p.colValueSource).getDataType().getDefaultFormatString());
                }

                rowMap.get(key).put(colName, table.get(i, p.colValueSource));
            }
        }

        GATKReportTable ret = new GATKReportTable(table.getTableName(), table.getTableDescription(), rowMap.keySet().size());
        for (String colName : distinctColNames){
            ret.addColumn(colName, colFormatMap.get(colName));
        }

        int rowIdx = 0;
        for (String key : rowMap.keySet()){
            for (GATKReportColumn col : ret.getColumnInfo()){
                Object val = null;
                if ("Gender".equals(col.getColumnName())){
                    Object sample = rowMap.get(key).get("Subject");
                    if (sample != null){
                        Sample s = sampleDB.getSample(String.valueOf(sample));
                        if (s != null){
                            val = s.getGender().name();
                        }
                    }
                }
                else {
                    val = rowMap.get(key).get(col.getColumnName());
                }

                if (val != null) {
                    ret.set(rowIdx, col.getColumnName(), val);
                }
            }

            rowIdx++;
        }

        return ret;
    }

    @Override
    public String getEvalModuleName() {
        return evalModuleName;
    }

    public static class Pivot {
        String colNameSource;
        String colValueSource;
        String suffix = null;

        public Pivot(String colNameSource, String colValueSource){
            this(colNameSource, colValueSource, null);
        }

        public Pivot(String colNameSource, String colValueSource, String suffix){
            this.colNameSource = colNameSource;
            this.colValueSource = colValueSource;
            this.suffix = suffix;
        }

        public String getTargetColName(String sourceCol){
            return sourceCol + (this.suffix == null ? "" : suffix);
        }
    }
}
