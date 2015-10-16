var textAlign = ["center", "center", "center", "start", "start", "center", "center", "start", "start"];
var dynatable = null,
    dt_dirty = false;
var alignmentViewer = null;

function myAttributeWriter(record) {
    // `this` is the column object in settings.columns
    var html = "";
    if (this.id == "export") {
        html = /*JSON.stringify(record)+*/'<input type="checkbox" id="' + record.index + '" class="selector" ';
        if (record['orthologId'] == "413057") { console.log(record); }
        if (record[this.id]=="1") {
            console.log("on!");
            html += "checked ";
        }
        html += '/>';
    }
    else {
        html = record[this.id];
    }
    return html;
};
  
function simpleCellWriter(column, record) {
    var html = column.attributeWriter(record)
        td = '<td';
    
    // add css style
    td += ' style="text-align: ' + textAlign[column.index] + ';"';
    
    console.log(record['export']);
    return td +' id="' +JSON.stringify(record)+ '">' + html + '</td>';
};

function dataRowWriter(rowIndex, record, columns, cellWriter) {
    var tr = '';

    // grab the record's attribute for each column
    for (var i = 0, len = columns.length; i < len; i++) {
      tr += cellWriter(columns[i], record);
    }

    return '<tr data-id="' + record['orthologId'] + '">' + tr + '</tr>';
  };

function selectableRowWriter(rowIndex, record, columns, cellWriter) {
    var tr = '<td>Möp!</td>';

    // grab the record's attribute for each column
    for (var i = 0, len = columns.length; i < len; i++) {
      tr += cellWriter(columns[i], record);
    }

    return '<tr>' + tr + '</tr>';
}

function updateDynatable(newRecords) {
    var dt = dynatable.data('dynatable');
    dt.records.updateFromJson({records: myRecords});
    dt.records.init();
    dt.process();
}

function updateAlignmentViewer(orthoId, showSeq) {
    console.log(orthoId);
    $('#ortholog-id').html(orthoId);
    alignmentViewer.records = alignments[orthoId];
    alignmentViewer.drawAlignment(showSeq); 
}

function record2Fasta(rec) {
    var fastaStr = "";
    fastaStr += ">" + rec.orthologId + "_fw\n";
    fastaStr += rec.fwSequence + "\n";
    fastaStr += ">" + rec.orthologId + "_rv\n";
    fastaStr += rec.rvSequence + "\n";
    
    return fastaStr;
}

function onDownloadFasta() {
    // get selected records
    var downloadStr = "";
    for (var i=0; i<myRecords.length; i++) {
        if (myRecords[i].export == '1') {
            downloadStr += record2Fasta(myRecords[i]);
        }
    }
    // provide export records for download
    document.location = 'data:Application/octet-stream,' +
                         encodeURIComponent(downloadStr);
}

function dtProcessingComplete() {
    $('#primer-table tr').click(function() {
        var orthoId = $(this).data("id");
        updateAlignmentViewer(orthoId, false);
    });
    $('input.selector').change(function() {
        var idx = parseInt($(this).attr('id'));
        if ($(this).is(':checked')) {
            myRecords[idx].export = "1";
        } else {
            myRecords[idx].export = "0";
        }
        dt_dirty = true; // remember to update the records later
    });
};

function finalizeSummary() {
    // generate summary text from template
    var template = $.templates("#tmplSummary");
    template.link("#resSummary", summary);
    
    // whip up pie charts
    var plotCat = $.jqplot('chartCategories', [categories], {
        title: 'Functional categories for discovered orthologs',
        seriesDefaults: {
          // make this a donut chart.
          renderer:$.jqplot.DonutRenderer,
          rendererOptions:{
            // Donut's can be cut into slices like pies.
            sliceMargin: 3,
            // Pies and donuts can start at any arbitrary angle.
            startAngle: -90,
            showDataLabels: true,
            // By default, data labels show the percentage of the donut/pie.
            // You can show the data 'value' or data 'label' instead.
            dataLabels: 'value'
          }
        }, 
        legend: { show:true, location: 'e' }
    });
/*
    var plotCat = $.jqplot('chartSubCategories', [subcats], {
        seriesDefaults: {
          // make this a donut chart.
          renderer:$.jqplot.DonutRenderer,
          rendererOptions:{
            // Donut's can be cut into slices like pies.
            sliceMargin: 2,
            // Pies and donuts can start at any arbitrary angle.
            startAngle: -90,
            showDataLabels: true,
            // By default, data labels show the percentage of the donut/pie.
            // You can show the data 'value' or data 'label' instead.
            dataLabels: 'value'
          }
        }, 
        legend: { 
            renderer: $.jqplot.EnhancedLegendRenderer,
            show:true, 
            location: 'e',
            rendererOptions: { numberRows : 12 } }
            
    });
*/
}


$( document ).ready(function() {
    $('#tabs').tabs();
    finalizeSummary();
//    alert(JSON.stringify(myRecords));    
    dynatable = $('#primer-table').dynatable({
        features: {
            //paginate: false,
            search: false
        },
        dataset: {
            perPageDefault: 10,
            records: myRecords
        },
        writers: {
            _attributeWriter: myAttributeWriter,
            _cellWriter: simpleCellWriter,
            _rowWriter: dataRowWriter
        }
    });
    dynatable.bind('dynatable:afterProcess', dtProcessingComplete);
    dynatable.bind('dynatable:afterUpdate', function() {
        if (dt_dirty) {
            dt_dirty = false;
            updateDynatable(myRecords);
        }
    });
    dtProcessingComplete(); // needs to be called manually once

    // get ortholog id for first primer pair
    var orthologId = myRecords[0]['orthologId'];

    //$('#alignment').html( alignments['413291'] );
    alignmentViewer = new CanvasState(document.getElementById('alignmentCanvas'));
    updateAlignmentViewer(orthologId, false);
    
    $('td.export-toggle').click(function() {
        this.toggleClass('selected');
    });
    $('#cb_showSeq').change(function() {
        var orthoId = $('#ortholog-id').text();
        if($(this).is(":checked")) {
            updateAlignmentViewer(orthoId, true);
        } else {
            updateAlignmentViewer(orthoId, false);
        }
    });
});



$('#chckHead').click(function () {
    if (this.checked == false) {
        $('.chcktbl:checked').attr('checked', false);
    } 
    else {
        $('.chcktbl:not(:checked)').attr('checked', true);
    }
});
