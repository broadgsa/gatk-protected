function processPlots(config){
    var toAppend = [];
    var navConfig = {};
    var plotData = {};
    var sectionMap = {};
    $.each(config.sections, function(idx, section){
        section.id = uniqID.get(section.label.toLowerCase());
        section.id = section.id.replace(/ /g, '_');

        //add to nav panel
        navConfig[section.label] = navConfig[section.label] || {label: section.label, id: section.id, children: []};

        //create div for section
        if (!sectionMap[section.label]){
            sectionMap[section.label] = section;
        }

        $.each(section.reports, function(idx, r){
            r.id = uniqID.get(r.label.toLowerCase());
            r.id = r.id.replace(/ /g, '_');
            r.plotDivId = r.id + '_plot';

            //add to NavPanel
            navConfig[section.label].children.push({
                id: r.id,
                label: r.label
            });

            //merge for nav panel
            r.data.config = r.data.config || {};
            r.data.config.id = r.plotDivId;
            plotData[r.plotDivId] = r.data;
        });
    });

    $('div.mainpage').append(toAppend);
    updateNavPanel(navConfig, sectionMap);

    window.mqc_compressed_plotdata = mqc_compressed_plotdata = plotData;
}

function updateNavPanel(navConfig, sectionMap){
    var toAppend = [];
    $.each(navConfig, function(key, val){
        var newEl = $('<li></li>', {
            html: '<a href="javascript:void(0)" class="nav-l1">' + val.label + '</a><ul></ul>'
        });

        newEl.click(function(){
            var sectionDiv = buildSectionDiv(sectionMap[val.label]);
            $.each(sectionMap[val.label].reports, function (idx, r) {
                //make plot div
                buildReportDiv(sectionDiv.get(1), r);
            });
            $('div #section_wrapper').empty();
            $('div #section_wrapper').append(sectionDiv);

            render_plots();
            render_tables();
            configure_toolbox();
        });

        if (val.children.length > 1) {
            $.each(val.children, function (idx, child) {
                var el = $('<li/>', {
                    html: '<a href="#' + child.id + '_div" class="nav-l2">' + child.label + '</a>'
                });
                newEl.find('ul').append(el);
            });
        }

        toAppend.push(newEl);
    });

    $('ul.mqc-nav').append(toAppend);
}

function buildSectionDiv(section){
    var sectionId = section.id || '';
    var sectionLabel = section.label || '';
    var descriptionHTML = section.descriptionHTML || '';

    var html = '<hr>' +
        '<div id="mqc-module-section-' + sectionId + '" class="mqc-module-section">' +
        '<h2 id="' + sectionId + '-header">' + sectionLabel + '</h2>' +
        '<p>' + descriptionHTML + '</p>' +
        '</div>';

    return $(html);
}

function buildReportDiv(parentDiv, config) {
    var buttonCfg = '';
    switch (config.data.plot_type) {
        case 'xy_line':
            buttonCfg = '<button class="btn btn-default btn-sm active" data-action="set_data" data-ylab="Count"  data-newdata="0" data-target="' + config.plotDivId + '">Counts</button>' +
                '<button class="btn btn-default btn-sm " data-action="set_data" data-ylab="Observed / Expected"  data-newdata="1" data-target="' + config.plotDivId + '">Obs/Exp</button>';
            return buildPlotDiv(parentDiv, config, 'hc-line-plot', buttonCfg);
        case 'bar_graph':
            buttonCfg = '<button class="btn btn-default btn-sm active" data-action="set_numbers" data-target="' + config.plotDivId + '" data-ylab="Number of Reads">Number of Reads</button>' +
                '<button class="btn btn-default btn-sm " data-action="set_percent" data-target="' + config.plotDivId + '" data-ylab="Percentages">Percentages</button>';
            return buildPlotDiv(parentDiv, config, 'hc-bar-plot', buttonCfg);
        case 'data_table':
            return buildTableDiv(parentDiv, config);

        //default:
        //    throw 'Unknown type: ' + config.data.plot_type;
    }
}

function buildPlotDiv(parentDiv, config, plotType, buttonCfg){
    var html = '<hr><h3 id="' + config.id + '_div">' + config.label + '</h3>' +
        '<div class="mqc_hcplot_plotgroup">' +
        '<div class="btn-group hc_switch_group">' + buttonCfg + '</div>' +
        '<div class="hc-plot-wrapper">' +
        '<div id="' + config.plotDivId + '" class="hc-plot not_rendered ' + plotType + '"><small>loading..</small></div>' +
        '</div>' +
        '</div>';

    $(parentDiv).append($(html));
}

function buildTableDiv(parentDiv, config){
    var html = generateTableHtml(config);
    $(parentDiv).append($(html));
}

var uniqID = {
    counter:0,
    get:function(prefix) {
        if(!prefix) {
            prefix = "uniqid";
        }
        var id =  prefix+""+uniqID.counter++;
        id = id.replace(/ /g, '_');
        id = id.replace(/\//g, '_');

        if($("#"+id).length == 0)
            return id;
        else
            return uniqID.get()
    }
};