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
                 {
                     label: "Trinity_Fusion",
                     url: sampleInfo.trinityBed,
                     displayMode: "EXPANDED",
                     color: "green",
                     indexed: false,
                     order: 30
                 },
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
fuction goToFusion( fusionChr, fusionPos ){
    igv.browser.search( fusionChr + ":" + Math.max( 0, parseInt( fusionPos ) - 50 ) + "-" + ( parseInt( fusionPos ) + 50 ) );
}

/** 
* Change array of fusion names to a html list.
*/
function toFusionList( fusionName ){
    return '<il>' + fusionName + '</il>;
} 

/**
* Load the fusion list from the json array to a html list
*/
function loadFusionList( sampleInfo ){
    $( "#fusionList" ).html( sampleInfo.fusionList.map( toFusionList ) );
}
