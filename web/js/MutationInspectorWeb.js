"use strict;"

///////////////////////
// Data
//////////////////////

/**
 * A small cach of SNP related info
 */
var mutationInspectorState = {
  cache : {}
};


///////////////////////
// Mutation Table Associated
//
///////////////////////

/**
 * Read in the JSON object that points to the bams and variants of interest.
 * @param {string} mutationTabFile - The path / URL to the JSON object
 */
function loadMutationTable( mutationTabFile ) {

  // Holds all information about the mutation view
  mutationInspectorView = {}

  // Read in the JSON file
  mutationInspectorView.json = dummy_data; // readMutationJSON( mutationTabFile )

  // Create table from json file
  // Get an array of the keys
  mutationHeaderKeys = []
  for( mutationHeader in mutationInspectorView.json.SNV[0] ) {
    if( mutationInspectorView.json.SNV[0].hasOwnProperty( mutationHeader ) ) {
      mutationHeaderKeys.push( mutationHeader );
    }
  }

  // Store locations of certain key row elements used later
  // These elements can be considered REQUIRED in the data
  mutationInspectorView.headerKeys = mutationHeaderKeys
  mutationInspectorView.chrKey = mutationHeaderKeys.indexOf( "CHROM" );
  mutationInspectorView.posKey = mutationHeaderKeys.indexOf( "POS" );
  mutationInspectorView.refKey = mutationHeaderKeys.indexOf( "REF" );
  mutationInspectorView.altKey = mutationHeaderKeys.indexOf( "ALT" );

  // Add header and footer elements to the table
  var mutationTable = $('#mutationTable');
  var mutationHeader = mutationHeaderKeys.map( toTableRowHeaderElement );
  mutationTable.append( '<thead><tr>' + mutationHeader.join('') + '</tr></thead>' );
  mutationTable.append( '<tfoot><tr>' + mutationHeader.join('') + '</tr></tfoot>' );

  // Add table body
  mutationTable.append( '<tbody>' );
  for( snvIndex = 0; snvIndex < mutationInspectorView.json.SNV.length; snvIndex++ ){
      snvEntry = mutationInspectorView.json.SNV[ snvIndex ];
      mutationEntryValues = mutationHeaderKeys.map( function( key ){
        return snvEntry[ key ]; } )
      mutationTable.append( '<tr>' + mutationEntryValues.map( toTableRowBodyElement ) + '</tr>' );
  }
  mutationTable.append( '</tbody>' );

  // Return the inspector
  return mutationInspectorView;
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
 * Make a table row containing a given value.
 * Helper function for the map call.
 * @param {string} tableRowValue - Value to put in the row.
 */
function toTableRowBodyElement( tableRowValue ){
    return '<td>' + tableRowValue + '</td>';
}


//////////////////////
// Associated with the specific view tab.
//
//////////////////////

/**
 * Add a tab for a specific genomic location / SNP.
 * Update the full UI to be consistent with selection.
 * Populate the body of the table with specific information.
 * @param {string} curRowChr - Chromsome location
 * @param {string} curRowPos - Genomic position on chromosome
 * @param {string} curRowRef - Reference base
 * @param {string} curRowAlt - Alternative base
 */
function addSpecificTab( curRowChr, curRowPos, curRowRef, curRowAlt ){
  var newTabName = curRowChr+"_"+curRowPos
  var tabArea = $( '#tabContent' );
  // If the the tab already exists, go to tab, do not make a new one.
  if( isExistingSpecificTab( curRowChr, curRowPos )){
    clickSpecificViewTab( newTabName )
    retrieveCRAVATInfo( curRowChr, curRowPos, curRowRef, curRowAlt )
    return;
  }
  // Make a new tab.
  var tabDescArea = $( '#tabDescription' );
  var chromLocation = curRowChr + ":" + curRowPos
  var closeButton = newTabName + "_close"
  var tabHeader = newTabName+"_tab"
  tabDescArea.append( '<li id="'+tabHeader+'"><a href="#'+newTabName+'" data-toggle="tab"><button id="'+closeButton+'" class="close closeTab" type="button">x</button>'+chromLocation+'</a></li>' );
  // Add in the tab area ( by default add in area for browser and cravat info )
  tabArea.append( '<div id="'+newTabName+'" class="tab-pane fade"></div>' );
  clickSpecificViewTab( newTabName );
  // Update cache of SNP info per genomic location.
  mutationInspectorState.cache[ chromLocation ] = { 
    'Chromosome' : curRowChr,
    'Position' : curRowPos,
    'alt' : curRowRef,
    'ref' : curRowAlt,
  }
  // MuPIT link will be added by an asynchronous call
  mutationInspectorState.cache[ chromLocation ][ "MuPIT Link" ] = null
  var currentCravatData = retrieveCRAVATInfo( curRowChr, curRowPos, curRowRef, curRowAlt );
  // Add click event for close button and tab.
  registerCloseEvent( closeButton, tabHeader, newTabName );
  registerOnClickEvent( tabHeader, chromLocation );
}

/**
 * Indicates if the tab already exists.
 * @param {string} curCheckChr - Chromsome to check
 * @param {string} curCheckPos - Postion on chromosome to check.
 */
function isExistingSpecificTab( curCheckChr, curCheckPos ){
  var newTabName = curCheckChr+"_"+curCheckPos
  var mutationTabs = $( '.tab-pane' );
  for( var tabIndex = 0; tabIndex < mutationTabs.length; tabIndex++ ){
    if( mutationTabs[ tabIndex ].id === newTabName ){
      return true;
    }
  }
  return false; 
}

/**
 * clicks on a specific tab to make it active.
 * @param {string} curSpecificViewTabId - The id of the tab to click and make active.
 */
function clickSpecificViewTab( curSpecificViewTabId ){
  $('.nav-tabs a[href="#' + curSpecificViewTabId + '"]').tab('show');
}

/**
 * Create custom close button event.
 * Closes associated tab and changes the active tab to the browsing tab.
 * @param {string } closeButtonId - Id of close button to which to add the event.
 * @param {string } closeTabId - Id of tab header to remove.
 * @param {string } closeBodyId - Id of tab content to remove.
 */
function registerCloseEvent( closeButtonId, closeTabId, closeBodyId ){
  // Add close button solution from
  // Hardcoded and not dynamic but works for now.
  $( "#"+closeButtonId ).click( function() {
    $( '#' + closeTabId ).remove();
    $( '#' + closeBodyId ).remove();
    $( '#tabDescription a[href="#tabBrowser"]' ).tab('show'); // Show the default tab body
    $( "#tabBrowser_tab" ).click();
  });
}

/**
 * Create a custom click event for the tabs.
 * Updates the page to be consistent with the active tab.
 * @param {string} tabHeader - The id of the tab to which to add the click event.
 * @param {string} registerChrLoc - The genomic location of the SNP being viewed (format= Chr:Pos)
 */
function registerOnClickEvent( tabHeader, registerChrLoc ){
  $( "#" + tabHeader ).click( function() {
    var currentState = mutationInspectorState.cache[ registerChrLoc ];
    goToSNP( currentState.Chromosome, currentState.Position );
    updateSNPInfo( currentState.Chromosome, currentState.Position, currentState.ref, currentState.alt );
    updateMupitLink( currentState );
  });
}

/**
 * Create a custom click event for the default tab.
 * This tab does not represent a specific location so some
 * of the page is cleared of info. The browser is not removed but stays
 * in it's last state.
 * @param {string} tabHeader - The id of the tab to which to add the click event.
 */
function registerDefaultTabClick( tabHeader ){
  $( "#"+tabHeader ).click( function() {
    updateSNPInfo( "NA", "NA", "NA", "NA" );
    updateMupitLink( { 'MuPIT Link' : null,
                       'Chromosome' : null 
    });
  });
}

/**
 * Update the top of the page with a summary of the location being viewed.
 * Also set the MuPIT link to a spinner as it will be loading.
 * @params {string} curSNPChr - Current view's chromosome.
 * @params {string} curSNPChr - Current view's position.
 * @params {string} curSNPChr - Current view's reference base.
 * @params {string} curSNPChr - Current view's alternative base.
 */
function updateSNPInfo( curSNPChr, curSNPPos, curSNPRef, curSNPAlt ){
  $( '#currentChr' ).text( curSNPChr );
  $( '#currentPosition' ).text( curSNPPos );
  $( '#currentRef' ).text( curSNPRef );
  $( '#currentAlt' ).text( curSNPAlt );
  $( '#currentMupit' ).text( '' );
  $( '#currentMupit' ).append( '<div class=\"spinner-loader\">Loading...</div>' );
}


///////////////////////
// IGV browser Associated
//
//////////////////////

/**
 * Initializes a IGV browser instance.
 * @params {string} sampleInfo - Object holding the bam url/path, bam index url/path, and sample name.
 */
function createIGVBrowser( sampleInfo ){
  // Create a browser
  var div = $("#igvBrowser")[0],
  options = {
    showNavigation: true,
    genome: "hg19",
    tracks: [{ url: sampleInfo.BAM,
               indexURL: sampleInfo.BAM_INDEX,
               type: "bam",
               label: sampleInfo.SAMPLE,
               height: 150 }]
  };
  igv.createBrowser( div, options );
}

/**
 * Go to SNP location.
 * @params {string} dataTableRowChr - Chromosomal location to which to move.
 * @params {string} dataTableRowPos - Position of interest
 */
function goToSNP( dataTableRowChr, dataTableRowPos ){
  // Move the igv browser to a specific location
  // Example "chr1:181,413,875-181,413,925"
  // The position needs to be a span so we are adding a window around the given position.
  igv.browser.search( dataTableRowChr + ":" + Math.max( 0, parseInt( dataTableRowPos ) - 50 ) + "-" + ( parseInt( dataTableRowPos ) + 50 ) );
}


//////////////////////
// Data IO
//
//////////////////////

/**
 * Reads in the mutation JSON file.
 * @params {string} readInFile - Path or URL to file.
 */
function readMutationJSON( readInFile ){
  $.getJSON( readInFile , function( jsonInfo ){
    mutationInspectorView.json = jsonInfo
  })
  .done( function(){ console.log( 'Completed reading file:' + readInFile ); } )
  .fail( function(){ console.log( 'Failed to read file:' + readInFile ); } )
}


//////////////////////
// CRAVAT Associated
//
/////////////////////

/**
 * Sets the area to contain the detailed (CRAVAT) info to a spinner
 * given we will wait for the associated asynchronous call.
 * @params {string} retrieveAnnotationTabName - The id of the tab content div to set to a spinner as we wait.
 */
function setAnnotationTabToLoad( retrieveAnnotationTabName ){
  $( '#' + retrieveAnnotationTabName ).html( "" );
  $( '#' + retrieveAnnotationTabName ).append( "<div class=\"spinner-loader\">Loading...</div>" );
}

/**
 * Query the CRAVAT web service for information about the genomic location of interest.
 * Update the page when the data is recevied.
 * Asyncronous call.
 * @params {string} retrieveChr - Chromosome of interest, used in the cravat call.
 * @params {string} retrievePos - Position of interest, used in the cravat call.
 * @params {string} retrieveRef - Reference base of interest, used in the cravat call.
 * @params {string} retrieveAlt - Alternative base of interest, used in the cravat call.
 */
function retrieveCRAVATInfo( retrieveChr, retrievePos, retrieveRef, retrieveAlt ){
  // Performs an asynchronous call to the CRAVAT web service
  // Updates both the CRAVAT info header and the info tab
  // Puts a loading logo up while waiting
  var positionKey = retrieveChr + "_" + retrievePos
  setAnnotationTabToLoad( positionKey );
  $.ajax({ type: 'GET',
           dataType: 'json',
           success: function( cravatData ){
    if( cravatData ){
      updateMupitLink( cravatData );
      updateCravatTab( retrieveChr + "_" + retrievePos, cravatData );
      mutationInspectorState.cache[ retrieveChr+':'+retrievePos ]["MuPIT Link"] = cravatData[ "MuPIT Link" ];
      }
    },
           url: "http://staging.cravat.us/rest/service/query?mutation="+retrieveChr+"_"+retrievePos+"_+_"+retrieveRef+"_"+retrieveAlt
  });
  return null;
}

/**
 * Write all information in a CRAVAT object received from the CRAVAT web service to a content tab /table.
 * @params {string} updateTab - Tab to add content to from CRAVAT.
 * @params {object} cravatItem - Obect of annotation, all members and value of the object are written to the table.
 */
function updateCravatTab( updateTab, cravatItem ){
  // Make CRAVAT annotation table for data
  var newTable = "<div class=\"table-responsive\"><table class=\"table-hover\">"
  for( var cravatElement in cravatItem ){
    if( cravatItem.hasOwnProperty( cravatElement ) ){
      var curValue = cravatItem[ cravatElement ]
      newTable = newTable + "<tr><td><b>" + cravatElement + ":</b> " + ( curValue ? curValue : "Not Specified" ) + "</td></tr>";
    }
  }
  newTable = newTable + "</table></div>"
  $( '#'+updateTab ).html( "" );
  $( '#'+updateTab ).append( newTable );
}


/////////////////////
// MuPIT Link / Button
//
/////////////////////

/**
 * Update the MuPIT link, handling cases where there was a link, there was no link, or a bad call occured.
 * @param {object} cravatItem - Object from the CRAVAT web service.
 */
function updateMupitLink( cravatItem ){
  var mupit = cravatItem[ "MuPIT Link" ];
  // If there is a MuPIT link add as a button (update the label to a label)
  if( mupit ){
    $( '#currentMupit' ).html( "" );
    $( '#currentMupit' ).append( '<button id=\"mupitButton\" class=\"btn btn-info\"> View in MuPIT</button>');
    $( '#mupitButton' ).click( function() {
      window.open( mupit );
    });
  } else if( cravatItem[ 'Chromosome' ] ){
    $( '#currentMupit' ).html( "" );
    $( '#currentMupit' ).html( 'No Link for '+cravatItem[ 'Chromosome' ]+':'+cravatItem[ 'Position' ] );
  } else {
    $( '#currentMupit' ).html( "" );
    $( '#currentMupit' ).html( 'Please select a variant' );
  }
}
