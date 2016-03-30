//var textAlign = ["center", "center", "center", "center", "center", "start", "start", "center", "center", "start", "start"];
//var dynatable = null,
//    dt_dirty = false;
var alignmentViewer = null;



function updateAlignmentViewer(markerId, showSeq) {
    console.log(markerId);
    $('#marker-id').html(markerId);
    alignmentViewer.records = alignments[markerId];
    alignmentViewer.drawAlignment(showSeq);
}

function record2Fasta(rec) {
    var fastaStr = "";
    fastaStr += ">" + rec[2] + "_" + rec[3] + "_fw\n";
    fastaStr += rec[7] + "\n";
    fastaStr += ">" + rec[2] + "_" + rec[3] + "_rv\n";
    fastaStr += rec[8] + "\n";

    return fastaStr;
}

function onDownloadFasta() {
    // get selected records
    var downloadStr = "";
    for (var i=0; i<myRecords.length; i++) {
        if (myRecords[i][1] == 1) {
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
          { title: "snps" },
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
            return '<input type="checkbox" value="' + $('<div/>').text(data).html() + '">';
          }
        },
        { targets: 3,
          render: function (data, type, row, meta){
            return row[2]+ '_' + data;
          }
        },
        { targets: [11,12],
          render: function (data, type, row, meta){
            if (data != 'None') {
              return "<a href='http://www.ncbi.nlm.nih.gov/nucleotide/" + data + "' target='_blank'>" + data +  "</a>";
            }
            else {
              return data;
            }
          }
        },
        { targets: 13,
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
                      '<tr class="group"><td colspan="12">'+group+'</td></tr>'
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

  // Handle click on checkbox
  $(tableId + ' tbody').on('click', 'input[type="checkbox"]', function(e) {
    var $row = $(this).closest('tr');
    // Get row data
    var data = table.row($row).data();
    // Get row ID
    var rowId = data[0];
    myRecords[rowId][1] = (myRecords[rowId][1]==0 ? 1 : 0);
  } );
}


function setupScatterSnps(data) {
  console.log("initializing scatter plot...");
  var margin = {top: 20, right: 20, bottom: 30, left: 40},
      width = 960 - margin.left - margin.right,
      height = 500 - margin.top - margin.bottom;

  var xValue = function(d) { return d.prod_len; },
      xScale = d3.scale.linear().range([0, width]);

  var yValue = function(d) { return d.n_snps; },
      yScale = d3.scale.linear().range([height, 0]);

  // fill color
  var cValue = function(d) { return +d.n_species; },
      color = d3.scale.category10();

  var xAxis = d3.svg.axis()
      .scale(xScale)
      .orient("bottom");

  var yAxis = d3.svg.axis()
      .scale(yScale)
      .orient("left");

  // add graph canvas to DOM
  var svg = d3.select("#chart-scatter-snps").append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom)
      .append("g")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  // add tooltip
  var tooltip = d3.select('#chart-scatter-snps').append('div')
    .attr('class', 'tooltip');
  tooltip.append('div')
    .attr('class', 'label'); // data field to display
  tooltip.append('div')
    .attr('class', 'x-value'); // data field to display
  tooltip.append('div')
    .attr('class', 'y-value'); // data field to display

  // load data
  //d3.tsv(filename, function(error, data) {
  //  if (error) throw error;

    // convert string values into numbers
    data.forEach(function(d) {
      d.prod_len = +d.prod_len;
      d.n_snps = +d.n_snps;
    });

    // don't want dots overlapping axis, so add in buffer to data domain
    xScale.domain([d3.min(data, xValue)-1, d3.max(data, xValue)+1]);
    yScale.domain([d3.min(data, yValue)-1, d3.max(data, yValue)+1]);
    xScale.domain(d3.extent(data, function(d) { return d.prod_len; })).nice();
    yScale.domain(d3.extent(data, function(d) { return d.n_snps; })).nice();

    svg.append("g")
      .attr("class", "x axis")
      .attr("transform", "translate(0," + height + ")")
      .call(xAxis)
    .append("text")
      .attr("class", "label")
      .attr("x", width)
      .attr("y", -6)
      .style("text-anchor", "end")
      .text("Product length (bp)");

    svg.append("g")
        .attr("class", "y axis")
        .call(yAxis)
      .append("text")
        .attr("class", "label")
        .attr("transform", "rotate(-90)")
        .attr("y", 6)
        .attr("dy", ".71em")
        .style("text-anchor", "end")
        .text("#SNPs")

    var dot = svg.selectAll(".dot")
        .data(data)
      .enter().append("circle")
        .attr("class", "dot")
        .attr("r", 3.5)
        .attr("cx", function(d) { return xScale(d.prod_len); })
        .attr("cy", function(d) { return yScale(d.n_snps); })
        .style("fill", function(d) { return color(d.n_species); });

    // show tooltip on mouseover
    dot.on('mouseover', function(d) {
      tooltip.select('.label').html('primer set <b>' + d.id + '</b>');
      tooltip.select('.x-value').html('product <b>' + d.prod_len + 'bp</b>');
      tooltip.select('.y-value').html('<b>' + d.n_snps + '</b> SNPs');
      tooltip.style('display', 'block');
    });
    dot.on('mouseout', function() {
      tooltip.style('display', 'none');
    });
    dot.on('mousemove', function(d) {
      tooltip.style('top', (d3.event.layerY + 10) + 'px')
        .style('left', (d3.event.layerX + 10) + 'px');
    });

    var legend = svg.selectAll(".legend")
        .data(color.domain())
        .enter().append("g")
          .attr("class", "legend")
          .attr("transform", function(d, i) { return "translate(0," + d * 20 + ")"; });

    legend.append("rect")
          .attr("x", width - 18)
          .attr("width", 18)
          .attr("height", 18)
          .style("fill", color)
          .style("stroke", color)
          .on("click", function(label) {
              var rect = d3.select(this);
              var enabled = true;
              if (rect.attr("class") === "disabled") {
                rect.attr("class", "");
              } else {
                rect.attr("class", "disabled");
                enabled = false;
              }

              display = enabled ? "inline" : "none";

              svg.selectAll(".dot")
                .filter(function(d) {return label == d.n_species;})
                .attr("display", display);
           });

    legend.append("text")
          .attr("x", width - 24)
          .attr("y", 9)
          .attr("dy", ".35em")
          .style("text-anchor", "end")
          .text(function(d) { return d; });
  //});
}


function setupBarMarkers(data) {

  var margin = {top: 20, right: 30, bottom: 30, left: 80},
      width = 420 - margin.left - margin.right,
      height = 100 - margin.top - margin.bottom,
      barHeight = height / data.length;

  var x = d3.scale.linear()
      .domain([0, d3.max(data, function(d) { return d.value; })])
      .range([0, width]);

  var y = d3.scale.ordinal()
      .domain(data.map(function(d) { return d.name; }))
      .rangeRoundBands([height, 0], .1);

  var xAxis = d3.svg.axis()
      .scale(x)
      .orient("bottom");

  var yAxis = d3.svg.axis()
      .scale(y)
      .orient("left");

  var chart = d3.select(".barchart")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom)
      .append("g")
      .attr("transform", "translate(" + margin.left + ", " + margin.top + ")");

  chart.append("g")
      .attr("class", "x axis")
      .attr("transform", "translate(0," + height + ")")
      .call(xAxis);

  chart.append("g")
      .attr("class", "y axis")
      .call(yAxis);

  var bar = chart.selectAll(".bar")
      .data(data)
    .enter().append("g")
      .attr("transform", function(d, i) { return "translate(0," + i * barHeight + ")" });

  bar.append("rect")
      .attr("class", "bar")
      .attr("height", y.rangeBand())
      .attr("width", function(d) { return x(d.value); });
      //.attr("transform", function(d) { return "translate(0," + y(d.name) + ")"; });



  /*bar.append("text")
      .attr("x", function(d) { return x(d.value) - 3; })
      .attr("y", y.rangeBand() / 2)
      .attr("dy", ".35em")
      .text(function(d) { return d.value; });*/
}

function type(d) {
  d.value = +d.value; // coerce to number
  return d;
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
        { title: "#Input files" }
      ]
    } );

    setupScatterSnps(primers);
    setupBarMarkers(species_markers_output);

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

    // initialize alignment viewer
    alignmentViewer = new CanvasState(document.getElementById('alignmentCanvas'));
    // select first primer pair displayed in table
    $('#primer-t tbody tr:eq(1)').click();

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
