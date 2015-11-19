function getAlignmentStats(aln) {
    var l = 0,
        n = 0,
        snps = [];

    // get maximum seq length as alignment length
    for (var seq_id in aln) {
        l = Math.max(l, aln[seq_id].length);
        n++;
    }
    // find variants in alignment
    for (var i=0; i<l; i++) {
        var chars = Object.create(null);
        for (var seq_id in aln) {
            if (aln[seq_id].charAt(i) != '-') {
                chars[aln[seq_id][i]] = true;
            }
        }
        if (Object.keys(chars).length > 1) {
            snps.push(i);
        }
    }

    return {'numSeqs': n, 'maxLen': l, 'snps': snps};
}

function getSeqCoords(seq) {
    var l=seq.length, pos = [], a=-1, gap=true;
    for (var i=0; i<l; i++) {
        // gap -> sequence change?
        if (gap && seq.charAt(i)!='-') {
            gap = false;
            a = i;
        // sequence -> gap change?
        } else if (!gap && (seq.charAt(i)=='-' || i==l-1)) {
            gap = true;
            pos.push([a,i])
        }
    }

    return pos;
}

function cartoonSeqWriter(id, seq, snps, context, x, y) {
    var l = seq.length;
    // draw guide line
    context.beginPath();
    context.moveTo(x, y);
    context.lineTo(x+l, y);
    context.strokeStyle = '#000000';
    context.lineWidth = 1;
    context.stroke();
    // identify sequence segments
    pos = getSeqCoords(seq);
    // primer or ref?
    if (id.match(/\d{6}_\d+/)) {
        context.strokeStyle = '#549499';
    } else {
        context.strokeStyle = '#222222';
    }
    // draw sequence segments
    for (var i=0; i<pos.length; i++) {
        var x1 = pos[i][0],
            x2 = pos[i][1];
//            alert(seq_id + ": " + x1 + "-" + x2);
        context.beginPath();
        context.moveTo(x+x1, y);
        context.lineTo(x+x2, y);
        context.lineWidth = 8;
        context.stroke();
    }

    // print variants in different color
    for (var i=0; i<snps.length; i++) {
        if (seq.charAt(snps[i]) != '-') {
            context.strokeStyle = '#dd0000';
            context.beginPath();
            context.moveTo(x+snps[i], y);
            context.lineTo(x+snps[i]+2, y);
            context.lineWidth = 8;
            context.stroke();
        }
    }
}

function simpleSeqWriter(id, seq, snps, context, x, y) {
    context.font = '6pt Monospace';
    console.log(id);
    if (id.match(/\d{6}_\d+/)) {
        context.fillStyle = '#22dd00';
    } else {
        context.fillStyle = '#222222';
    }
    context.fillText(seq, x, y);

    // print variants in different color
    for (var i=0; i< snps.length; i++) {
        if (seq.charAt(snps[i]) != '-') {
            var x_snp = x + context.measureText(seq.substr(0,snps[i])).width;
            context.fillStyle = '#dd0000';
            context.fillText(seq.charAt(snps[i]), x_snp, y);
        }
    }
}

function CanvasState(canvas) {
  // **** First some setup! ****

  this.records = null;
  this.container = canvas.parentNode;
  this.canvas = canvas;
  this.width = canvas.width;
  this.height = canvas.height;
  this.ctx = canvas.getContext('2d');
  // This complicates things a little but but fixes mouse co-ordinate problems
  // when there's a border or padding. See getMouse for more detail
  var stylePaddingLeft, stylePaddingTop, styleBorderLeft, styleBorderTop;
  if (document.defaultView && document.defaultView.getComputedStyle) {
    this.stylePaddingLeft = parseInt(document.defaultView.getComputedStyle(canvas, null)['paddingLeft'], 10)      || 0;
    this.stylePaddingTop  = parseInt(document.defaultView.getComputedStyle(canvas, null)['paddingTop'], 10)       || 0;
    this.styleBorderLeft  = parseInt(document.defaultView.getComputedStyle(canvas, null)['borderLeftWidth'], 10)  || 0;
    this.styleBorderTop   = parseInt(document.defaultView.getComputedStyle(canvas, null)['borderTopWidth'], 10)   || 0;
  }
  // Some pages have fixed-position bars (like the stumbleupon bar) at the top or left of the page
  // They will mess up mouse coordinates and this fixes that
  var html = document.body.parentNode;
  this.htmlTop = html.offsetTop;
  this.htmlLeft = html.offsetLeft;

  // **** Keep track of state! ****

  this.zoom = false;
  this.dragging = false; // Keep track of when we are dragging
  this.dragoffx = 0; // See mousedown and mousemove events for explanation
  this.dragoffy = 0;

  // **** The events! ****

  // This is an example of a closure!
  // Right here "this" means the CanvasState. But we are making events on the Canvas itself,
  // and when the events are fired on the canvas the variable "this" is going to mean the canvas!
  // Since we still want to use this particular CanvasState in the events we have to save a reference to it.
  // This is our reference!
  var myState = this;

  // fixes a problem where double clicking causes text to get selected on the canvas
  canvas.addEventListener('selectstart', function(e) { e.preventDefault(); return false; }, false);
  // Up, down, and move are for dragging
  canvas.addEventListener('mousedown', function(e) {
    var mouse = myState.getMouse(e);
    myState.dragoffx = mouse.x;
    myState.dragoffy = mouse.y;
    myState.dragging = true;
    console.log('mousedown ('+mouse.x+','+mouse.y+')');
  }, true);
  canvas.addEventListener('mousemove', function(e) {
    if (myState.dragging){
      var mouse = myState.getMouse(e);
      // scroll in the direction the mouse moved
      var posx = myState.container.scrollLeft;
      myState.container.scrollLeft = posx - (mouse.x - myState.dragoffx);
      myState.dragoffx = mouse.x;
    }
  }, true);
  canvas.addEventListener('mouseup', function(e) {
    myState.dragging = false;
  }, true);
  // double click to toggle zoom
  canvas.addEventListener('dblclick', function(e) {
    var mouse = myState.getMouse(e);
    var scrolled = myState.container.scrollLeft / myState.width;
    if (myState.zoom) {
        myState.drawAlignment(false);
        //myState.container.scrollLeft = scrolled * myState.width;
    } else {
        myState.drawAlignment(true);
    }
    myState.zoom = !myState.zoom;
    //myState.addShape(new Shape(mouse.x - 10, mouse.y - 10, 20, 20, 'rgba(0,255,0,.6)'));
  }, true);

  // **** Options! ****

  this.selectionColor = '#CC0000';
  this.selectionWidth = 2;
  this.interval = 30;
  //setInterval(function() { myState.draw(); }, myState.interval);
}

CanvasState.prototype.drawAlignment = function(showSeq) {
    showSeq = typeof showSeq !== 'undefined' ? showSeq : false;

    var aln = this.records;
    // get alignment length and depth
    var aln_stats = getAlignmentStats(aln);
//    alert(JSON.stringify(aln_stats));
    var l = aln_stats['maxLen'],
        n = aln_stats['numSeqs'],
        snps = aln_stats['snps'];

    // initialize canvas
    var x = 200,
        y = 50,
        ySpacing = 20
        xSpacing = showSeq ? 6 : 1;
    var canvas = document.getElementById('alignmentCanvas');
    var context = canvas.getContext('2d');
    var width
    canvas.width  = xSpacing*l+x;
    canvas.height = y+((n+1)*ySpacing);
    //var styles = ['#000000', ''];

    // clear canvas
    canvas.width = canvas.width;

    for (var seq_id in aln) {
        // write seq name
        context.fillStyle = '#000000';
        //var new_seq_id = seq_id.replace(/(\d{6})_\d-/, '$1_primers_');
        context.fillText(seq_id, 0, y);
        if (showSeq) {
            simpleSeqWriter(seq_id, aln[seq_id], snps, context, x, y);
        } else {
            cartoonSeqWriter(seq_id, aln[seq_id], snps, context, x, y);
        }
        y += ySpacing;
    }
}

// Creates an object with x and y defined, set to the mouse position relative to the state's canvas
// If you wanna be super-correct this can be tricky, we have to worry about padding and borders
CanvasState.prototype.getMouse = function(e) {
  var element = this.canvas, offsetX = 0, offsetY = 0, mx, my;

  // Compute the total offset
  if (element.offsetParent !== undefined) {
    do {
      offsetX += element.offsetLeft;
      offsetY += element.offsetTop;
    } while ((element = element.offsetParent));
  }

  // Add padding and border style widths to offset
  // Also add the <html> offsets in case there's a position:fixed bar
  offsetX += this.stylePaddingLeft + this.styleBorderLeft + this.htmlLeft;
  offsetY += this.stylePaddingTop + this.styleBorderTop + this.htmlTop;

  mx = e.pageX - offsetX;
  my = e.pageY - offsetY;

  // We return a simple javascript object (a hash) with x and y defined
  return {x: mx, y: my};
}
