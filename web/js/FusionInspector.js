var fusionInspectorState = {
  cache : {}
}

function loadFusionIGV( sampleInfo )
{
    var div = $("#igvBrowser")[0],
    options = {
        showKaryo: false,
        showNavigation: true,
        fastaURL: sampleInfo.reference,
        cytobandURL: sampleInfo.cytoband,
        tracks: [{
                     type: "sequence",
                     order: 1
                 },
                 {
                     label: "Reference_Annotations",
                     url: sampleInfo.referenceBed,
                     displayMode: "SQUISHED",
                     indexed: false,
                     order: 10
                 },
                 {
                     label: "Junction_Spanning_View",
                     type: "FusionJuncSpan",
                     url: sampleInfo.junctionSpanning,
                     displayMode: "SQUISHED",
                     indexed: false,
                     order: 20
                 },
//                 {
//                     label: "Trinity_Fusion",
//                     url: sampleInfo.trinityBed,
//                     displayMode: "EXPANDED",
//                     color: "green",
//                     indexed: false,
//                     order: 30
//                 },
                 {
                     url: sampleInfo.junctionReads,
                     label: "Junction_Reads",
                     displayMode: "SQUISHED",
                     minHeight: 30,
                     maxHeight: 100,
                     color: "cyan",
                     indexed: false,
                     order: 40
                 },
                 {
                     url: sampleInfo.junctionReadsBam,
                     indexURL: sampleInfo.junctionReadsBai,
                     label: "Junction_Reads_Alignments",
                     height: 100,
                     indexed: true,
                     visibilityWindow: 2000000,
                     order: 45
                 },
                 {
                     url: sampleInfo.spanningReads,
                     label: "Spanning_Reads",
                     displayMode: "SQUISHED",
                     minHeight: 30,
                     maxHeight: 100,
                     color: "purple",
                     indexed: false,
                     order: 50
                 },
                 {
                     url: sampleInfo.spanningReadsBam,
                     indexURL: sampleInfo.spanningReadsBai,
                     label: "Spanning_Reads_Alignments",
                     height: 100,
                     indexed: true,
                     visibilityWindow: 2000000,
                     order: 55
                 }]};
    igv.createBrowser(div, options);
}

/**
* Move IGV browser to a fusion of choice by genomic location.
*/
function goToFusion( fusionChr, fusionInfo ){
    var numBreakRight = parseInt( fusionInfo.breakRight )
    var numBreakLeft = parseInt( fusionInfo.breakLeft )
    var padLength = ( numBreakRight - numBreakLeft ) / 2
    
    igv.browser.search( fusionChr + ":" + Math.max( 0, numBreakLeft - padLength ) + "-" + ( numBreakRight + padLength ) );
    setFusionDetails( fusionChr, fusionInfo );
}

function setFusionDetails( fusionName, fusionInfo ){
    $( "#FusionNameDetail" ).html( "<p class='navbar-text'><b>Fusion Name:</b> " + fusionName + "</p>" );
    $( "#FusionBreakLeftDetail" ).html( "<p class='navbar-text'><b>Left Break Position:</b> " + fusionInfo.breakLeft  + "</p>" );
    $( "#FusionBreakRightDetail" ).html( "<p class='navbar-text'><b>Right Break Position:</b> " + fusionInfo.breakRight + "</p>" );
    $( "#FusionJunctionDetail" ).html( "<p class='navbar-text'><b>Junction Read Count:</b> " + fusionInfo.numJunctionReads + "</p>" );
    $( "#FusionSpanningDetail" ).html( "<p class='navbar-text'><b>Spanning Read Count:</b> " + fusionInfo.numSpanningFrags + "</p>" );
}

/** 
* Change array of fusion names to a html list.
*/
function toFusionList( fusionName ){
    return '<li id='+fusionName+'><a href="#">' + fusionName + '</a></li>';
} 

/**
* Load the fusion list from the json array to a html list
*/
function loadFusionList( sampleInfo ){
    var fusionList = []
    for( var fusionName in sampleInfo.fusions ){
        if( sampleInfo.fusions.hasOwnProperty( fusionName ) ){
            fusionList.push( fusionName );
        }
    }
    fusionInspectorState.cache[ "fusionList" ] = fusionList
    $( "#fusionList" ).html( fusionInspectorState.cache.fusionList.map( toFusionList ) );
    for( var fusionNameItr=0; fusionNameItr < fusionInspectorState.cache.fusionList.length; fusionNameItr ++ ){
        registerClickFusionMenu( fusionInspectorState.cache.fusionList[ fusionNameItr ] );
    }
}

/**
 * Set sample name in header.
 */
function setSampleName( sampleInfo ){
    $( "#sampleId" ).html( "<p class='navbar-text'><b>Sample:</b> "+sampleInfo.sampleName+"</p>" )
}

function registerClickFusionMenu( fusionId ){
    $( "#" + fusionId ).click( function() {
       goToFusion( fusionId, fusionInspectorState.cache.json.fusions[ fusionId ] );
    })
}
