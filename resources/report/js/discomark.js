//var textAlign = ["center", "center", "center", "center", "center", "start", "start", "center", "center", "start", "start"];
//var dynatable = null,
//    dt_dirty = false;
var alignmentViewer = null;

// obsolete (DynaTable)
// function myAttributeWriter(record) {
//     // `this` is the column object in settings.columns
//     var html = "";
//     if (this.id == "export") {
//         html = /*JSON.stringify(record)+*/'<input type="checkbox" id="' + record.index + '" class="selector" ';
//         if (record['markerId'] == "412698") { console.log(record); }
//         if (record[this.id]=="1") {
//             console.log("on!");
//             html += "checked ";
//         }
//         html += '/>';
//     }
//     // include primer pair index in ortholog_id field
//     else if (this.id == "markerId") {
//         html = record[this.id] + "_" + record['ps_idx'];
//     }
//     // insert a link for NCBI records
//     else if ((this.id == "fwBlastHit" || this.id == "rvBlastHit") && (record[this.id] != "None")) {
//         html = "<a href='http://www.ncbi.nlm.nih.gov/nuccore/" + record[this.id] + "' target='_blank'>" + record[this.id] + "</a>";
//     }
//     else {
//         html = record[this.id];
//     }
//     return html;
// };

// obsolete (DynaTable)
// function simpleCellWriter(column, record) {
//     var html = column.attributeWriter(record)
//         td = '<td';
//
//     // add css style
//     td += ' style="text-align: ' + textAlign[column.index] + ';"';
//
//     console.log(record['export']);
//     return td +' id="' +JSON.stringify(record)+ '">' + html + '</td>';
// };

// obsolete (DynaTable)
// function dataRowWriter(rowIndex, record, columns, cellWriter) {
//     var tr = '';
//
//     // grab the record's attribute for each column
//     for (var i = 0, len = columns.length; i < len; i++) {
//       tr += cellWriter(columns[i], record);
//     }
//
//     return '<tr data-id="' + record['markerId'] + '">' + tr + '</tr>';
//   };

// obsolete (DynaTable)
// function selectableRowWriter(rowIndex, record, columns, cellWriter) {
//     var tr = '<td>Möp!</td>';
//
//     // grab the record's attribute for each column
//     for (var i = 0, len = columns.length; i < len; i++) {
//       tr += cellWriter(columns[i], record);
//     }
//
//     return '<tr>' + tr + '</tr>';
// }

// obsolete (DynaTable)
// function updateDynatable(newRecords) {
//     var dt = dynatable.data('dynatable');
//     dt.records.updateFromJson({records: myRecords});
//     dt.records.init();
//     dt.process();
// }

function updateAlignmentViewer(markerId, showSeq) {
    console.log(markerId);
    $('#marker-id').html(markerId);
    alignmentViewer.records = alignments[markerId];
    alignmentViewer.drawAlignment(showSeq);
}

function record2Fasta(rec) {
    var fastaStr = "";
    fastaStr += ">" + rec.markerId + "_" + rec.ps_idx + "_fw\n";
    fastaStr += rec["fwSequence-(5'-3')"] + "\n";
    fastaStr += ">" + rec.markerId + "_" + rec.ps_idx + "_rv\n";
    fastaStr += rec["rvSequence-(5'-3')"] + "\n";

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

function setupPrimerTable(tableId) {
  $(tableId).DataTable( {
      data: myRecords,
      columns: [
          { title: "idx" },
          { title: "export" },
          { title: "markerId" },
          { title: "primerSet"},
          { title: "species" },
          { title: "product" },
          { title: "fw sequence" },
          { title: "rv sequence" },
          { title: "Tm" },
          { title: "primer length" },
          { title: "fw BLAST hit" },
          { title: "rv BLAST hit" },
          { title: "annotation"}
      ],
      columnDefs: [
        { visible: false, targets: [0,2] },
        { orderData: [4,2], targets: [4] },
        { targets: 1,
          searchable: false,
          orderable: false,
          className: 'dt-body-center',
          render: function (data, type, row, meta){
            return '<input type="checkbox" name="id[]" value="' + $('<div/>').text(data).html() + '">';
          }
        },
        { targets: 3,
          render: function (data, type, row, meta){
            return row[2]+ '_' + data;
          }
        },
        { targets: [10,11],
          render: function (data, type, row, meta){
            if (data != 'None') {
              return "<a href='http://www.ncbi.nlm.nih.gov/nucleotide/" + data + "' target='_blank'>" + data +  "</a>";
            }
            else {
              return data;
            }
          }
        },
        { targets: 12,
          render: function (data, type, row, meta){
            var html = '';
            var terms = data.split(',');
            for (i=0; i<terms.length; ++i) {
              if (terms[i].substring(0, 3) == 'GO:') {
                html += "<a href='http://www.ebi.ac.uk/QuickGO/GTerm?id=" + terms[i] + "' target='_blank'>" + terms[i] +  "</a>";
              } else {
                html += data;
              }
            }
            return html;
          }
        }
      ],
      order: [[4, 'desc'], [2, 'asc']],
      displayLength: 25,
      // add checkbox to mark rows for export
      drawCallback: function ( settings ) {
          var api = this.api();
          var rows = api.rows( {page:'current'} ).nodes();
          var last = null;

          api.column(2, {page:'current'} ).data().each( function ( group, i ) {
              if ( last !== group ) {
                  $(rows).eq( i ).before(
                      '<tr class="group"><td colspan="11">'+group+'</td></tr>'
                  );

                  last = group;
              }
          } );
        }
  } );

  // make rows selectable
  var table = $(tableId).DataTable();
  $(tableId + ' tbody').on( 'click', 'tr', function () {
    if ( $(this).hasClass('selected') ) {
      $(this).removeClass('selected');
    }
    else {
      if (Object.keys(this).length > 0) {
        table.$('tr.selected').removeClass('selected');
        $(this).addClass('selected');
        var markerId = myRecords[this._DT_RowIndex][2];
        updateAlignmentViewer(markerId, false);
      }
    }
  } );
}

function setupVenn(div_id, dataset) {
  var chart = venn.VennDiagram()
                  .width(400)
                  .height(400);

  var div = d3.select(div_id);
  div.datum(dataset).call(chart);

  var tooltip = d3.select("body").append("div")
      .attr("class", "venntooltip")

  div.selectAll("path")
      .style("stroke-opacity", 0)
      .style("stroke", "#fff")
      .style("stroke-width", 0)

  div.selectAll("g")
  .on("mouseover", function(d, i) {
      // sort all the areas relative to the current item
      venn.sortAreas(div, d);

      // Display a tooltip with the current size
      tooltip.transition().duration(400).style("opacity", .9);
      tooltip.text("{" + d.sets + "}: " + d.size);

      // highlight the current path
      var selection = d3.select(this).transition("tooltip").duration(400);
      selection.select("path")
          .style("stroke-width", 3)
          .style("fill-opacity", d.sets.length == 1 ? .4 : .1)
          .style("stroke-opacity", 1);
  })

  .on("mousemove", function() {
      tooltip.style("left", (d3.event.pageX) + "px")
             .style("top", (d3.event.pageY - 28) + "px");
  })

  .on("mouseout", function(d, i) {
      tooltip.transition().duration(400).style("opacity", 0);
      var selection = d3.select(this).transition("tooltip").duration(400);
      selection.select("path")
          .style("stroke-width", 0)
          .style("fill-opacity", d.sets.length == 1 ? .25 : .0)
          .style("stroke-opacity", 0);
  });
}

function finalizeSummary() {
    // generate summary text from template
    var template = $.templates("#tmplSummary");
    template.link("#resSummary", summary);

    // whip up venn diagrams (data included in counts.js)
    setupVenn("#chart-venn-markers-input", marker_sets_input);
    $('#tab-venn-markers-input').DataTable( {
      data: species_key,
      paging: false,
      searching: false,
      info: false,
      columns: [
        { title: "ID" },
        { title: "Species" },
        { title: "#Input markers" }
      ]
    } );

    // whip up pie charts
    /*var plotCat = $.jqplot('chartCategories', [categories], {
        title: 'Functional categories of markers',
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
    });*/
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
    // populate species vs. primers table
    $('#tabSumSpecies').DataTable( {
      data: species,
      paging: false,
      searching: false,
      info: false,
      columns: [
        { title: "#Species" },
        { title: "#Orthologs" },
        { title: "#PrimerPairs" }
      ]
    } );
}

$( document ).ready(function() {
    setupPrimerTable('#primer-t');
    finalizeSummary();
    $('#tabs').tabs();
//    alert(JSON.stringify(myRecords));
    /*dynatable = $('#primer-table').dynatable({
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
    */

    // get marker id for first primer pair
    //var markerId = myRecords[0]['markerId'];
    var markerId = myRecords[0][2];

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
