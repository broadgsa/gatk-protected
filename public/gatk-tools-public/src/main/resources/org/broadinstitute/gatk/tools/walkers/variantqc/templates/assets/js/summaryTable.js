function generateTableHtml(config){

    var rowCount = config.data.samples.length;
    var colCount = config.data.columns.length + 1;

    var html = '<h3 id="' + config.id + '_div">' + config.label + '</h3>' +
        '<button type="button" class="mqc_table_copy_btn btn btn-default btn-sm" data-clipboard-target="#' + config.id + '">' +
        '<span class="glyphicon glyphicon-copy"></span> Copy table' +
        '</button>' +
        '<button type="button" class="mqc_table_configModal_btn btn btn-default btn-sm" data-toggle="modal" data-target="#' + config.id + '_configModal">' +
        '<span class="glyphicon glyphicon-th"></span> Configure Columns' +
        '</button>' +
        '<button type="button" class="mqc_table_sortHighlight btn btn-default btn-sm" data-target="#' + config.id + '" data-direction="desc" style="display:none;">' +
        '<span class="glyphicon glyphicon-sort-by-attributes-alt"></span> Sort by highlight' +
        '</button>' +
        '<button type="button" class="mqc_table_makeScatter btn btn-default btn-sm" data-toggle="modal" data-target="#tableScatterModal" data-table="#' + config.id + '">' +
        '<span class="glyphicon glyphicon glyphicon-stats"></span> Plot' +
        '</button>' +
        '<small id="' + config.id + '_numrows_text" class="mqc_table_numrows_text">Showing <sup id="' + config.id + '_numrows" class="mqc_table_numrows">' + rowCount + '</sup>/<sub>' + rowCount + '</sub> rows and <sup id="' + config.id + '_numcols" class="mqc_table_numcols">' + colCount + '</sup>/<sub>' + colCount + '</sub> columns.</small>' +
        '<div id="' + config.id + '_container" class="mqc_table_container">' +
        '<div class="table-responsive">' +
        '<table id="' + config.id + '" class="table table-condensed mqc_table">';

    //<th id="header_spor_Assigned" class="chroma-col spor_Assigned " data-chroma-scale="PuBu" data-chroma-max="104.413184" data-chroma-min="0.0" data-namespace="featureCounts"  data-shared-key=read_count><span data-toggle="tooltip" title="featureCounts: Assigned reads (millions)">M Assigned</span></th><th id="header_pqcg_uniquely_mapped_percent" class="chroma-col pqcg_uniquely_mapped_percent " data-chroma-scale="YlGn" data-chroma-max="100.0" data-chroma-min="0.0" data-namespace="STAR" ><span data-toggle="tooltip" title="STAR: % Uniquely mapped reads">% Aligned</span></th><th id="header_liyq_uniquely_mapped" class="chroma-col liyq_uniquely_mapped " data-chroma-scale="PuRd" data-chroma-max="104.413184" data-chroma-min="0.0" data-namespace="STAR"  data-shared-key=read_count><span data-toggle="tooltip" title="STAR: Uniquely mapped reads (millions)">M Aligned</span></th><th id="header_ejtv_percent_trimmed" class="chroma-col ejtv_percent_trimmed " data-chroma-scale="RdYlBu-rev" data-chroma-max="100.0" data-chroma-min="0.0" data-namespace="Cutadapt" ><span data-toggle="tooltip" title="Cutadapt: % Total Base Pairs trimmed">% Trimmed</span></th><th id="header_tjid_percent_duplicates" class="chroma-col tjid_percent_duplicates " data-chroma-scale="RdYlGn-rev" data-chroma-max="100.0" data-chroma-min="0.0" data-namespace="FastQC" ><span data-toggle="tooltip" title="FastQC: % Duplicate Reads">% Dups</span></th><th id="header_spka_percent_gc" class="chroma-col spka_percent_gc " data-chroma-scale="Set1" data-chroma-max="100.0" data-chroma-min="0.0" data-namespace="FastQC" ><span data-toggle="tooltip" title="FastQC: Average % GC Content">% GC</span></th><th id="header_bvyk_avg_sequence_length" class="chroma-col bvyk_avg_sequence_length hidden" data-chroma-scale="RdYlGn" data-chroma-max="97.5523706167" data-chroma-min="0.0" data-namespace="FastQC" ><span data-toggle="tooltip" title="FastQC: Average Sequence Length (bp)">Length</span></th><th id="header_enxo_percent_fails" class="chroma-col enxo_percent_fails hidden" data-chroma-scale="Reds" data-chroma-max="100.0" data-chroma-min="0.0" data-namespace="FastQC" ><span data-toggle="tooltip" title="FastQC: Percentage of modules failed in FastQC report (includes those not plotted here)">% Failed</span></th><th id="header_jwkr_total_sequences" class="chroma-col jwkr_total_sequences " data-chroma-scale="Blues" data-chroma-max="104.413184" data-chroma-min="0.0" data-namespace="FastQC"  data-shared-key=read_count><span data-toggle="tooltip" title="FastQC: Total Sequences (millions)">M Seqs</span></th></tr></thead>';
    var data = config.data;
    html += '<thead><tr><th class="rowheader">Sample Name</th>';
    $.each(data.columns, function(colIdx, column){
        column.name = column.name || column.label.toLowerCase().replace(/ /g, '');
        column.colId = 'c' + colIdx;
        column.chroma = column.chroma || '';
        var showBar = !!column.dmax;
        var chromaScale = brewer_scales && brewer_scales.length > colIdx ? brewer_scales[colIdx] : 'RdYlGn';

        html += '<th id="header_' + column.colId + '_' + column.name + '" class="' + (showBar ? 'chroma-col ' : '') + column.colId + '_' + column.name + ' " ' + (showBar ? 'data-chroma-scale="' + chromaScale + '" data-chroma-max="' + column.dmin + '" data-chroma-min="' + column.dmax + '"' : '') + ' data-namespace="' + (column.category || column.colId) + '" ><span data-toggle="tooltip" title="' + (column.category ? column.category + ': ' : '') + (column.description || column.label) + '">' + column.label + '</span></th>';
    });

    html += '</thead><tbody>';

    //rows
    $.each(data.samples, function(sampleIdx, sample){
        html += '<tr><th class="rowheader" data-original-sn="' + sample + '">' + sample + '</th>';

        var rowData = data.datasets[sampleIdx];

        $.each(data.columns, function(colIdx, column){
            if (rowData && rowData.length > colIdx){
                var val = rowData[colIdx];
                var formattedVal = column.formatString ? numeral(val).format(column.formatString) : val;
                var percentage = 100.0;
                var showBar = !!column.dmax;
                if (column.dmax){
                    percentage = ((parseFloat(val) - (column.dmin || 0)) / (column.dmax - (column.dmin || 0))) * 100;
                }

                if (val < column.flagBelow || val > column.flagAbove){
                    html += '<td ' + (showBar ? 'class="data-coloured-flagged ' + column.colId + '_' + column.name + ' "' : '') + '><div class="wrapper">' + (showBar ? '<span class="bar" style="width:' + percentage + '%;"></span>' : '') + '<span class="val">' + formattedVal + '</span></div></td>';
                } else {
                    html += '<td ' + (showBar ? 'class="data-coloured ' + column.colId + '_' + column.name + ' "' : '') + '><div class="wrapper">' + (showBar ? '<span class="bar" style="width:' + percentage + '%;"></span>' : '') + '<span class="val">' + formattedVal + '</span></div></td>';
                }

            }
            else {
                html += '<td>ND</td>';
            }
        });
    });

    html += '</tbody>';
    html += '</table>';

    var modalHtml = '' +
        '<!-- MultiQC Table Columns Modal -->' +
        '<div class="modal fade" id="' + config.id + '_configModal" tabindex="-1">' +
            '<div class="modal-dialog modal-lg">' +
                '<div class="modal-content">' +
                '<div class="modal-header">' +
                    '<button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
                    '<h4 class="modal-title">General Statistics: Columns</h4>' +
                '</div>' +
                '<div class="modal-body">' +
                    '<p>Uncheck the tick box to hide columns. Click and drag the handle on the left to change order.</p>' +
                    '<p>' +
                        '<button class="btn btn-default btn-sm mqc_configModal_bulkVisible" data-target="#' + config.id + '" data-action="showAll">Show All</button>' +
                        '<button class="btn btn-default btn-sm mqc_configModal_bulkVisible" data-target="#' + config.id + '" data-action="showNone">Show None</button>' +
                    '</p>' +
                    '<table class="table mqc_table mqc_sortable mqc_configModal_table" id="' + config.id + '_configModal_table">' +
                        '<thead>' +
                            '<tr>' +
                                '<th class="sorthandle" style="text-align:center;">Sort</th>' +
                                '<th style="text-align:center;">Visible</th>' +
                                '<th>Group</th>' +
                                '<th>Column</th>' +
                                '<th>Description</th>' +
                                '<th>ID</th>' +
                                '<th>Scale</th>' +
                            '</tr>' +
                        '</thead>' +
                        '<tbody>';

    $.each(data.columns, function(colIdx, column) {
        modalHtml += '<tr class="enxo_percent_fails" style="background-color: rgba(255,127,0, 0.15);">' +
            '<td class="sorthandle ui-sortable-handle">||</span></td>' +
            '<td style="text-align:center;">' +
            '<input class="mqc_table_col_visible" type="checkbox"  checked="checked" value="' + (column.colId + '_' + column.name) + '" data-target="#' + config.id + '">' +
            '</td>' +
            '<td>' + (column.category || 'None') + '</td>' +
            '<td>' + column.label + '</td>' +
            '<td>' + (column.description || '') + '</td>' +
            '<td><code>' + (column.name || '') + '</code></td>' +
            '<td>None</td>' +
            '</tr>';
    });

    modalHtml += '</tbody>' +
                    '</table>' +
                '</div>' +
                '<div class="modal-footer"> <button type="button" class="btn btn-default" data-dismiss="modal">Close</button> </div>' +
            '</div>' +
        '</div>' +
    '</div>';

    return html + modalHtml;
}