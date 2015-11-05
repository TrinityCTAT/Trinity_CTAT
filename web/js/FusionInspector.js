"use strict;"

var fusionInspectorState = {
  cache : { tabs : [] }
}

function loadIGVBrowser()
{
    var divBrowser = $("#igvBrowser")[0],
    options = {
        showKaryo: false,
        showNavigation: true,
        fastaURL: fusionInspectorState.cache.json.reference,
        cytobandURL: fusionInspectorState.cache.json.cytoband,
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
//                 {
//                     label: "Trinity_Fusion",
//                     url: fusionInspectorState.cache.json.trinityBed,
//                     displayMode: "EXPANDED",
//                     color: "green",
//                     indexed: false,
//                     order: 30
//                 },
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
                     label: "Spanning_Reads_Alignments",
                     height: 100,
                     indexed: true,
                     visibilityWindow: 2000000,
                     order: 55
                 }]};
    fusionInspectorState.cache[ "curBrowser" ] = igv.createBrowser(divBrowser, options);
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
    }, 2000);
}

function setFusionDetails( fusionName, fusionJunctReads, fusionSpanFrags ){
    $( "#FusionNameDetail" ).html( "<p class='navbar-text'><b>Fusion Name:</b> " + fusionName + "</p>" );
    $( "#FusionJunctionDetail" ).html( "<p class='navbar-text'><b>Junction Read Count:</b> " + fusionJunctReads + "</p>" );
    $( "#FusionSpanningDetail" ).html( "<p class='navbar-text'><b>Spanning Read Count:</b> " + fusionSpanFrags + "</p>" );
}

/** 
* Change array of fusion names to a html list.
*/
//function toFusionList( fusionName ){
//    return '<li id='+fusionName+'><a href="#">' + fusionName + '</a></li>';
//} 

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

function toTableBodyElement( fusionEntry, orderedHeaderKeys ){
  var bodyRow = '<tr>';
  for( var headerKeyIndex = 0; headerKeyIndex < orderedHeaderKeys.length; headerKeyIndex++ ){
    bodyRow = bodyRow + '<td>' + fusionEntry[ orderedHeaderKeys[ headerKeyIndex ] ] + '</td>';
  }
  return( bodyRow + '</tr>');
}

function getFusionAnnotationFromRow( infoHeader, dataRow ){
  var index = fusionInspectorState.cache.fusionKeys.indexOf( infoHeader );
  if( index == -1 ){
    return( None );
  }
  return( dataRow[ index ] );
}

//TODO
//function addFusionTab( leftGene, rightGene, leftChr, rightChr, leftPos, rightPos ){
//  var tabName = leftGene + "_" + rightGene + "_" + leftChr + "_" + rightChr + "_" + leftPos + "_" + rightPos;
//  if( fusionInspectorState.cache.tabs.indexOf( tabName ) < 0 ){
//    // Add tab name to list
//    fusionInspectorState.cache.tabs.push( tabName );
//    // Make tab
//    var tabTitle = leftGene + ":" + rightGene;
//    var tabDescription = $( '#tabDescription' );
//    var tabContent = $( 'tabContent' );
//    tabDescription.append( '<li id="'+tabName+'_tab"><a hreaf="#'+tabName+'# data-toggle="tab"><button id="'+tabName+'_close" class="close closeTab" type="button">x</button>'+tabTitle+'</a></li>' );
//    tabContent.append( '<div id="'+tabName+'" class="tab-pane fade"></div>' );
//    addIGVBrowser( tabName, leftChr, leftPos );
//    registerCloseEvent( tabName+'_close', tabName+'_tab', tabName );
//  }
//  //Focus on Tab
//  $('.nav-tabs a[href="#' + tabName + '"]' ).tab( 'show' );
//}

//function registerCloseEvent( closeButtonId, closeTabId, closeBodyId ){
//  // Add close button solution from
//  // Hardcoded and not dynamic but works for now.
//  $( "#"+closeButtonId ).click( function() {
//    $( '#' + closeTabId ).remove();
//    $( '#' + closeBodyId ).remove();
//    $( '#tabDescription a[href="#tabBrowser"]' ).tab('show'); // Show the default tab body
//    $( "#tabBrowser_tab" ).click();
//    var tabIndex = fusionInspectorState.cache.tabs.indexOf( closeBodyId );
//    if( tabIndex > -1 ){
//      fusionInspectorState.cache.tabs.splice( tabIndex, 1 );
//    }
//  });
//}

//function isExistingTab( leftGene, rightGene, leftChr, rightChr, leftPos, rightPos ){
//  var tabName = leftGene + "_" + rightGene + "_" + leftChr + "_" + rightChr + "_" + leftPos + "_" + rightPos;
//  for( var tabIndex = 0; tabIndex < fusionInspectorState.cache.tabs.length; tabIndex++ ){
//    if( fusionInspectorState.cache.tabs[ tabIndex ] === tabName ){
//      return true;
//    }   
//  }
//  return false; 
//}

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
function setSampleName( sampleInfo ){
    $( "#sampleId" ).html( "<p class='navbar-text'><b>Sample:</b> "+sampleInfo.sampleName+"</p>" )
}

//function registerClickFusionMenu( fusionId ){
//    $( "#" + fusionId ).click( function() {
//       goToFusion( fusionId, fusionInspectorState.cache.json.fusions[ fusionId ] );
//    })
//}
