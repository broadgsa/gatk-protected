package org.broadinstitute.gatk.tools.walkers.variantqc;

import com.google.gson.JsonArray;
import org.apache.commons.io.IOUtils;
import org.broadinstitute.gatk.utils.io.Resource;

import java.io.IOException;
import java.io.PrintStream;
import java.text.SimpleDateFormat;
import java.util.Collection;
import java.util.Date;

/**
 * Created by bimber on 5/18/2017.
 */
public class HtmlGenerator {
    public static final String[] JS_SCRIPTS = new String[]{
            "templates/assets/js/packages/jquery-3.1.1.min.js",
            "templates/assets/js/packages/jquery-ui.min.js",
            "templates/assets/js/packages/numeral.min.js",
            //"https://cdnjs.cloudflare.com/ajax/libs/numeral.js/2.0.6/numeral.min.js",
            "templates/assets/js/packages/bootstrap.min.js",
            "templates/assets/js/packages/highcharts.js",
            "templates/assets/js/packages/highcharts.exporting.js",
            "templates/assets/js/packages/highcharts.heatmap.js",
            "templates/assets/js/packages/highcharts.offline-exporting.js",
            "templates/assets/js/packages/jquery.tablesorter.min.js",
            "templates/assets/js/packages/clipboard.min.js",
            "templates/assets/js/packages/FileSaver.min.js",
            //"https://cdnjs.cloudflare.com/ajax/libs/chroma-js/1.3.3/chroma.min.js",
            "templates/assets/js/packages/chroma.min.js",
            "templates/assets/js/packages/lz-string.min.js",
            "templates/assets/js/summaryTable.js",
            "templates/assets/js/shim.js"
    };

    public static final String[] JS_SCRIPTS2 = new String[]{
            "templates/assets/js/multiqc_tables.js",
            "templates/assets/js/multiqc_toolbox.js",
            "templates/assets/js/multiqc.js",
            "templates/assets/js/multiqc_plotting.js"
    };

    public static final String[] CSS_FILES = new String[]{
            "templates/assets/css/bootstrap.min.css",
            "templates/assets/css/default_multiqc.css",
            "templates/assets/css/font.css"
    };

    public HtmlGenerator() {

    }

    public void generateHtml(Collection<SectionJsonDescriptor> translatorList, PrintStream out) throws IOException {

        //append header
        Resource header = new Resource("templates/template1.html", VariantQC.class);
        IOUtils.copy(header.getResourceContentsAsStream(), out);

        //scripts:
        for (String script : CSS_FILES){
            Resource r = new Resource(script, VariantQC.class);
            out.println("<style>");
            IOUtils.copy(r.getResourceContentsAsStream(), out);
            out.println("</style>");
        }

        for (String script : JS_SCRIPTS){
            appendScript(script, out);
        }

        //config:
        out.println("<script type=\"text/javascript\">");
        out.println("mqc_plots = {};");
        out.println("num_datasets_plot_limit = 50;");
        out.println("$(function() {");
        String dateStr = new SimpleDateFormat("yyyy-MM-dd HH:mm").format(new Date());
        out.println("$('#dateTime').html('Generated on: " + dateStr + "');");
        out.println("processPlots({");
        out.println("sections:");

        JsonArray arr = new JsonArray();
        for (SectionJsonDescriptor t : translatorList){
            arr.add(t.getConfig());
        }
        out.println(arr.toString());
        out.println("});");
        out.println("});");
        out.println("</script>");

        for (String script : JS_SCRIPTS2){
            appendScript(script, out);
        }

        //append header
        Resource header2 = new Resource("templates/template2.html", VariantQC.class);
        IOUtils.copy(header2.getResourceContentsAsStream(), out);
    }

    private void appendScript(String script, PrintStream out) throws IOException{
        Resource r = new Resource(script, VariantQC.class);
        out.println("<script type=\"text/javascript\">");
        IOUtils.copy(r.getResourceContentsAsStream(), out);
        out.println("</script>");
    }
}
