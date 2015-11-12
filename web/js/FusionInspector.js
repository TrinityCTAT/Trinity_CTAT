"use strict;"

/**
 * State of the page.
 */
var fusionInspectorState = {
  cache : { tabs : [] },
  browserMade : false
}

/**
 * Make the IGVBrowser and dispay at a certain loci.
 */
function makeIGVBrowser( curFusion )
{
    // Load browser
    var divBrowser = $("#igvBrowser")[0],
    options = {
        showKaryo: false,
        showNavigation: true,
        fastaURL: fusionInspectorState.cache.json.reference,
        cytobandURL: fusionInspectorState.cache.json.cytoband,
        locus: curFusion,
        tracks: [{
                     type: "sequence",
                     order: 1
                 },
                 {
                     label: "Reference_Annotations",
                     url: fusionInspectorState.cache.json.referenceBed,
                     displayMode: "SQUISHED",
                     indexed: false,
                     order: 10
                 },
                 {
                     label: "Junction_Spanning_View",
                     type: "FusionJuncSpan",
                     url: fusionInspectorState.cache.json.junctionSpanning,
                     displayMode: "SQUISHED",
                     indexed: false,
                     order: 20
                 },
                 {
                     url: fusionInspectorState.cache.json.junctionReads,
                     label: "Junction_Reads",
                     displayMode: "SQUISHED",
                     minHeight: 30,
                     maxHeight: 100,
                     color: "cyan",
                     indexed: false,
                     order: 40
                 },
                 {
                     url: fusionInspectorState.cache.json.junctionReadsBam,
                     indexURL: fusionInspectorState.cache.json.junctionReadsBai,
                     label: "Junction_Reads_Alignments",
                     type: "bam",
                     height: 100,
                     indexed: true,
                     visibilityWindow: 2000000,
                     order: 45
                 },
                 {
                     url: fusionInspectorState.cache.json.spanningReads,
                     label: "Spanning_Reads",
                     displayMode: "SQUISHED",
                     minHeight: 30,
                     maxHeight: 100,
                     color: "purple",
                     indexed: false,
                     order: 50
                 },
                 {
                     url: fusionInspectorState.cache.json.spanningReadsBam,
                     indexURL: fusionInspectorState.cache.json.spanningReadsBai,
                     type: "bam",
                     label: "Spanning_Reads_Alignments",
                     height: 100,
                     indexed: true,
                     visibilityWindow: 2000000,
                     order: 55
                 }]};
    if( ! ( fusionInspectorState.cache.json.trinityBed === "NA" ) ){
        options.tracks.push({
            label: "Trinity_Fusion",
            url: fusionInspectorState.cache.json.trinityBed,
            displayMode: "EXPANDED",
            color: "green",
            indexed: false,
            order: 30
        });
    }
    fusionInspectorState.cache[ "curBrowser" ] = igv.createBrowser(divBrowser, options);
    fusionInspectorState.browserMade = true;
}

/**
 * If the IGV browser is not made then add a tab and load the igv browser.
 * Then update the position of the browser if needed and update the header information.
 */
function loadIGVBrowser( curFusion, curJunctionReads, curSpanningReads, curLeftPos, curRightPos )
{
    // Update header no matter if the browser is being made.
    setFusionDetails( curFusion, curJunctionReads, curSpanningReads );
    // Do not remake browser.
    if( ! fusionInspectorState.browserMade ){

        // Add tab for IGV and click on tab
        $('#tabDescription').append('<li role="presentation" id="igvTab"><a href="#igvBrowser" data-toggle="tab">IGV Detail</a></li>')
        $('#tabContent').append('<div role="tabpanel" id="igvBrowser" class="tab-pane fade"></div>')
        $('.nav-tabs a[href="#igvBrowser"]').tab('show');
 
        // Need a wait to make sure the tabs have switched.
        setTimeout(function(){
            makeIGVBrowser( curFusion );
        }, 1000);

    } else {
        goToFusion( curFusion, curRightPos, curLeftPos );
    }
}

/**
* Move IGV browser to a fusion of choice by genomic location.
*/
function goToFusion( fusionChr ){ //, fusionBreakRight, fusionBreakLeft ){
    //var numBreakRight = parseInt( fusionBreakRight )
    //var numBreakLeft = parseInt( fusionBreakLeft )
    //var padLength = ( numBreakRight - numBreakLeft ) / 2
    //fusionInspectorState.cache.curBrowser.search( fusionChr + ":" + Math.max( 0, numBreakLeft - padLength ) + "-" + ( numBreakRight + padLength ) );
    var location = fusionChr
    $('.nav-tabs a[href="#igvBrowser"]').tab('show');
    setTimeout(function(){
        fusionInspectorState.cache.curBrowser.search( location );
    }, 1000);
}

/**
 * Update the header information on the page.
 */
function setFusionDetails( fusionName, fusionJunctReads, fusionSpanFrags ){
    $( "#FusionNameDetail" ).html( "<p class='navbar-text'><b>Fusion Name:</b> " + fusionName + "</p>" );
    $( "#FusionJunctionDetail" ).html( "<p class='navbar-text'><b>Junction Read Count:</b> " + fusionJunctReads + "</p>" );
    $( "#FusionSpanningDetail" ).html( "<p class='navbar-text'><b>Spanning Read Count:</b> " + fusionSpanFrags + "</p>" );
}

/**
 * Order the mutation table keys so certain element are in front in a specific order
 * @param {array} ArrayToOrder - Array of elements to reorder.
 * @param {array} forcedOrder - Array of elements that should be at the beginning and in this order.
 */
function orderTableKeysBeginning( arrayToOrder, forcedOrder ){
  newArray = [];
  for( arrayElement = 0; arrayElement < arrayToOrder.length; arrayElement++ ){
    if( !( forcedOrder.indexOf( arrayToOrder[ arrayElement ] ) >= 0 )){
      newArray.push( arrayToOrder[ arrayElement ] );
    }
  }
  return( forcedOrder.concat( newArray ));
}

/**
 * Make a table header row containing a given value.
 * Helper function for the map call.
 * @param {string} tableRowValue - Value to put in the header row.
 */
function toTableRowHeaderElement( tableRowValue ){
    return '<th>' + tableRowValue + '</th>';
}

/**
 * Makes a table row from fusion information, making usre to order the data as given.
 */
function toTableBodyElement( fusionEntry, orderedHeaderKeys ){
  var bodyRow = '<tr>';
  for( var headerKeyIndex = 0; headerKeyIndex < orderedHeaderKeys.length; headerKeyIndex++ ){
    bodyRow = bodyRow + '<td>' + fusionEntry[ orderedHeaderKeys[ headerKeyIndex ] ] + '</td>';
  }
  return( bodyRow + '</tr>');
}

/**
 * Get the annotation info form the data table row.
 */
function getFusionAnnotationFromRow( infoHeader, dataRow ){
  var index = fusionInspectorState.cache.fusionKeys.indexOf( infoHeader );
  if( index == -1 ){
    return( None );
  }
  return( dataRow[ index ] );
}

/**
* Load the fusion list from the json array to a html list
*/
function loadFusionDataTable( ){

    // Forced order of the mutation table elements.
    // Any element in the table and not in this array
    // will be after these elements in the table and will
    // be in no specific order
    var forcedHeaderKeyOrder = ['Fusion', 'Junction Reads', 'Spanning Fragments', 'Splice Type', 'Left Gene', 'Left Chr', 'Left Pos', 'Left Strand', 'Right Gene', 'Right Chr', 'Right Pos', "Right Strand"];

    var fusionKeys = [];
    for( var fusionKey in fusionInspectorState.cache.json.fusions[ 0 ] ){
      if( fusionInspectorState.cache.json.fusions[ 0 ].hasOwnProperty( fusionKey ) ){
        fusionKeys.push( fusionKey );
      }
    }
    fusionInspectorState.cache[ "fusionKeys" ] = orderTableKeysBeginning( forcedHeaderKeyOrder, fusionKeys );
    // Make data table header and footer
    var fusionTable = $('#fusionTable');
    var fusionHeader = fusionInspectorState.cache.fusionKeys.map( toTableRowHeaderElement );
    fusionTable.append( '<thead><tr>' + fusionHeader.join('') + '</tr></thead>' );
    fusionTable.append( '<tfoot><tr>' + fusionHeader.join('') + '</tr></tfoot>' );

    // Add data table body
    fusionTable.append( '<tbody>' );
    for( var fusionIndex = 0; fusionIndex < fusionInspectorState.cache.json.fusions.length; fusionIndex++ ){
      var fusionEntry = fusionInspectorState.cache.json.fusions[ fusionIndex ];
      fusionTable.append( toTableBodyElement( fusionEntry, fusionInspectorState.cache.fusionKeys ) );
    }
    fusionTable.append( '</tbody>' );
}

/**
 * Set sample name in header.
 */
//function setSampleName( sampleInfo ){
//    $( "#sampleId" ).html( "<p class='navbar-text'><b>Sample:</b> "+sampleInfo.sampleName+"</p>" )
//}
