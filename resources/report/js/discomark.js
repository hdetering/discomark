var textAlign = ["center", "center", "center", "center", "start", "start", "center", "center", "start", "start"];
var dynatable = null,
    dt_dirty = false;
var alignmentViewer = null;

function myAttributeWriter(record) {
    // `this` is the column object in settings.columns
    var html = "";
    if (this.id == "export") {
        html = /*JSON.stringify(record)+*/'<input type="checkbox" id="' + record.index + '" class="selector" ';
        if (record['markerId'] == "413057") { console.log(record); }
        if (record[this.id]=="1") {
            console.log("on!");
            html += "checked ";
        }
        html += '/>';
    }
    // insert a link for NCBI records
    else if (this.id == "fwBlastHit" || this.id == "rvBlastHit") {
      if (record[this.id] != "None") {
        html = "<a href='http://www.ncbi.nlm.nih.gov/nuccore/" + record[this.id] + "' target='_blank'>" + record[this.id] + "</a>";
      }
      else {
        html = record[this.id];
      }
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

    return '<tr data-id="' + record['markerId'] + '">' + tr + '</tr>';
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

function updateAlignmentViewer(markerId, showSeq) {
    console.log(markerId);
    $('#marker-id').html(markerId);
    alignmentViewer.records = alignments[markerId];
    alignmentViewer.drawAlignment(showSeq);
}

function record2Fasta(rec) {
    var fastaStr = "";
    fastaStr += ">" + rec.markerId + "_fw\n";
    fastaStr += rec.fwSequence + "\n";
    fastaStr += ">" + rec.markerId + "_rv\n";
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
        var markerId = $(this).data("id");
        updateAlignmentViewer(markerId, false);
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
        title: 'Functional categories for discovered markers',
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
            pushState: false,
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

    // get marker id for first primer pair
    var markerId = myRecords[0]['markerId'];

    //$('#alignment').html( alignments['413291'] );
    alignmentViewer = new CanvasState(document.getElementById('alignmentCanvas'));
    updateAlignmentViewer(markerId, false);

    $('td.export-toggle').click(function() {
        this.toggleClass('selected');
    });
    $('#cb_showSeq').change(function() {
        var mId = $('#marker-id').text();
        if($(this).is(":checked")) {
            updateAlignmentViewer(mId, true);
        } else {
            updateAlignmentViewer(mId, false);
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
