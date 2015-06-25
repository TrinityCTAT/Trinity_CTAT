function createBrowser( inspectorView ){
  // Create a browser with default (first listed sample in the sample menu)
  var defaultSample = $("#sampleMenu li")[0].innerText;
  var div = $("#igvDiv")[0],
  options = {
    showNavigation: true,
    genome: "hg19",
    tracks: [{ url: inspectorView[ defaultSample ][ "DNA" ],
               label: "Exome Sequencing (" + defaultSample + ")",
               height: 150 },
             { url: inspectorView[ defaultSample ][ "RNA" ],
               label: "RNA-Seq Sequencing (" + defaultSample + ")",
               color: "rgb( 102, 153, 255 )",
               height: 150 }]
  };
  igv.createBrowser( div, options );
}

function createMenus( InspectorViewCurrent ){
  // Make sample menu
  var defaultSample = Object.keys( inspectorView )[ 0 ];
  var sampleMenuDropdown =  $("#sampleMenu") 
  for( var sample in inspectorView ){
    if( inspectorView.hasOwnProperty( sample )){
      sampleMenuDropdown.append( "<li><a href=\"#\">" + sample + "</a></li>" );
    }
  }
  // Add sample menu event
  // Trigger the event on the defaul sample menu
  $("#sampleMenu li").click( function(){
    switchSample( $(this).text() );
  });

  // Update error class menus
  updateErrorMenus( defaultSample );
}

function updateErrorMenus( sample ){
  // Update the menus for the error classes based on the given sample.
  var TPMenuDropdown =  $("#truePositive")
  var FPMenuDropdown =  $("#falsePositive") 
  var FNMenuDropdown =  $("#falseNegative") 

  var TPObj = inspectorView[ sample ].TP;
  var FPObj = inspectorView[ sample ].FP;
  var FNObj = inspectorView[ sample ].FN;

  // Clear the previous list entries
  // This clears the event as well so you have to add them back
  TPMenuDropdown.empty();
  FPMenuDropdown.empty();
  FNMenuDropdown.empty();
  
  // Make TP Menu
  for( var TPLabel in TPObj ){
    if( TPObj.hasOwnProperty( TPLabel )){
      TPMenuDropdown.append("<li><a href=\"#\">"+TPLabel+"</a></li>" );
    }
  }
  // Make FP Menu
  for( var FPLabel in FPObj ){
    if( FPObj.hasOwnProperty( FPLabel )){
      FPMenuDropdown.append("<li><a href=\"#\">"+FPLabel+"</a></li>" );
    }
  }
  // Make FN Menu
  for( var FNLabel in FNObj ){
    if( FNObj.hasOwnProperty( FNLabel )){
      FNMenuDropdown.append("<li><a href=\"#\">"+FNLabel+"</a></li>" );
    }
  }

  // Add click event to the error class locations
  $("#truePositive li").click( function(){
    moveToSNP( $(this).text(), TPObj[ $(this).text() ]);
    snpClass.innerHTML = "<b>True Positive</b>"
  });
  $("#falsePositive li").click( function(){
    moveToSNP( $(this).text(), FPObj[ $(this).text() ]);
    snpClass.innerHTML = "<b>False Positive</b>";
  });
  $("#falseNegative li").click( function(){
    moveToSNP( $(this).text(), FNObj[ $(this).text() ]);
    snpClass.innerHTML = "<b>False Negative</b>";
  });
}

function formatSNPLocationForBrowser( locationToFormat ){
  // Changes the format of a SNP location from what is viewed to what is needed for the browser.
  var SNPLocation = parseInt( getSNPLocation( locationToFormat ) );
  return getChr( locationToFormat ) + ":" + Math.max( 0, SNPLocation - 30 ) + "-" + ( SNPLocation + 30 )
}

function getChrCoverage( locationToFormat ){
  // Get Chr coverage from UI formated location
  return locationToFormat.split(" ")[1].replace("(","").replace(")","");
}

function getChr( locationToFormat ){
  // Get just the Chr location from UI formated location
  chrTemp = locationToFormat.split(" ")[0].split("-")[0].toLowerCase();
  if ( chrTemp[ 3 ] === "x" ){
    chrTemp = chrTemp.substring( 0, 3 ) + "X" + chrTemp.substring( 4, chrTemp.length );
  } else if ( chrTemp[ 3 ] === "y" ){
    chrTemp = chrTemp.substring( 0, 3 ) + "Y" + chrTemp.substring( 4, chrTemp.length );
  } else if ( chrTemp[ 3 ] === "m" ){
    chrTemp = chrTemp.substring( 0, 3 ) + "M" + chrTemp.substring( 4, chrTemp.length );
  }
  return chrTemp;
}

function getChrLocation( locationToFormat ){
  // Get the Chr location from UI formated location
  return locationToFormat.split(" ")[0].replace("-",":");
}

function updateCRAVATAnnotationTable( cravatItem ){
  // Make CRAVAT annotation table for data

  var annotationArea = $("#annotationTableDiv")
  annotationArea.empty();

  // Add table
  var hugoSymbol = cravatItem["HUGO symbol"];
  var thousandFreq = cravatItem["1000 Genomes allele frequency"];
  var cosmic = cravatItem["Occurences in COSMIC [exact nucleotide change]"];
  var mupit = cravatItem["MuPIT Link"];
  var mupitButton = "";
  var geneCards = cravatItem["GeneCards summary"];

  //If there is a MuPIT link add as button and add to annotation header html
  if( mupit ){
    mupitButton = "<button id=\"mupitButton\" class=\"btn btn-info\"><span class=\"glyphicon glyphicon-eye-open\"></span> View in MuPIT</button>";
  }
  annotationArea.append(
    "<div class=\"col-xs-3\"><div class=\"table-responsive\"><table class=\"table-hover\">" +
      "<tr><td><b>HUGO Symbol:</b> " + ( hugoSymbol ? hugoSymbol : "Not Specified" ) + "</td></tr>" +
      "<tr><td><b>1000 Genomes Freq:</b> " + parseFloat( thousandFreq ? thousandFreq : "Not Specified" ).toFixed(4) + "</td></tr>" +
      "<tr><td><b>COSMIC Occurences:</b> " + ( cosmic ? cosmic : "Not Specified" ) + "</td></tr>" +
    "</table></div></div>" +
    "<div class=\"col-xs-6\"><div class=\"table-responsive\"><table class=\"table-hover\">" +
      "<tr><td><b>MuPIT Link:</b> "+ ( mupitButton ? mupitButton : "Not Specified" )+"</td></tr>" +
      "<tr><td><b>GeneCards Summary:</b> "+ ( geneCards ? geneCards : "Not Specified" ) +"</td></tr>" +
    "</table></div></div>" );
  //Add click event for MuPIT
  if( mupit ){
    $('#mupitButton').click( function(){
      window.open( mupit );
    });
  }
}

function addCravatTab( inputCravatData ){
  // Given an object form a CRAVAT reponse.
  // Drop the contents of the CRAVAT object into a table
  // and add it to the CRAVAT div. 
  var cravatTab = $("#cravatTab");
  cravatTab.empty()

  // Make table from object
  var newTable = "<div class=\"table-responsive\"><table class=\"table-hover\">"
  for( var cravatElement in inputCravatData ){
    if( inputCravatData.hasOwnProperty( cravatElement ) ){
      var curValue = inputCravatData[ cravatElement ]
      newTable = newTable + "<tr><td><b>" + cravatElement + ":</b> " + ( curValue ? curValue : "Not Specified" ) + "</td></tr>";
    }
  }
  newTable = newTable + "</table></div>"

  // Add to div
  cravatTab.append( newTable ); 
}

function resetCRAVATArea(){
  // Reset the content of the CRAVAT Tab and header
  var cravatTab = $("#cravatTab");
  var cravatHeader = $("#annotationTableDiv");
  cravatTab.empty();
  cravatHeader.empty();

  cravatTab.append("<p>No Variant Annotation Loaded.</p>");
  cravatHeader.append("<p>No Variant Annotation Loaded.</p>");
}

function setCRAVATAreaToLoading(){
  // Set the CRAVAT area to loading

  // Clear annotation area
  var cravatTab = $("#cravatTab");
  var cravatHeader = $("#annotationTableDiv");
  cravatTab.empty();
  cravatHeader.empty();
 
  cravatTab.append("<div class=\"spinner-loader\">Loading...</div>");
  cravatHeader.append("<div class=\"spinner-loader\">Loading...</div>");
}

function getSNPLocation( locationToFormat ){
  // Get the SNP location from UI formated location
  return locationToFormat.split(" ")[0].split("-")[1];
}

function moveToSNP( location, SNPInfo ){
  // Change the view on the genomic tracks to the selected location

  // Update the location information
  snpLocation.innerHTML = "<b>SNP Location:</b> " + getChrLocation( location );
  snpCoverage.innerHTML = "<b>SNP Coverage:</b> " + SNPInfo[ "Cov_dna" ] + " (DNA) " + SNPInfo[ "Cov" ] + " (RNA)";
  snpRef.innerHTML = "<b>Ref:</b> " + SNPInfo[ "Ref" ];
  snpAlt.innerHTML = "<b>Alt:</b> " + SNPInfo[ "Alt" ];

  // Move to the SNP in the track
  // Incoming formate Chr#-loc# (coverage#)
  // Chr1-34 (43)
  // Needed format
  // chr1:43
  igv.browser.search( formatSNPLocationForBrowser( location ) );

  // Update CRAVAT info
  setCRAVATAreaToLoading();
  // Retrieve and update CRAVAT areas 
  retrieveCRAVATInfo( SNPInfo );
}

function readData( pipelineFile ){

  $.getJSON( pipelineFile, function(jsonFile){

    console.log( "Successfully read file.");
    console.log( jsonFile );
  
    // Set the global inspector view 
    inspectorView = jsonFile;
 
    // Initialize View
    createMenus( jsonFile );
    // Trigger the even on the default sample menu
    $("#sampleMenu li")[0].click();

    // Update the CRAVAT area
    resetCRAVATArea();
  })
  .done( function(){ console.log( "Completed reading file:"+pipelineFile );})
  .fail( function(){ console.log( "Failed to read file." );});
}

function retrieveCRAVATInfo( SNPinfo ){
  // Performs an asynchronous call to the CRAVAT web service
  // Updates both the CRAVAT info header and the info tab
  // Puts a loading logo up while waiting.
  // SNPInfo needs to be { chr: '22', Loc: '30421786', Ref: 'A', Alt: 'T', Strand: '+' }

  $.ajax({ type: 'GET', 
           dataType: 'json', 
           success: function( cravatData ) {
                                     if( cravatData ){
                                       updateCRAVATAnnotationTable( cravatData );
                                       addCravatTab( cravatData );
                                                     }
                                  }, 
           url: "http://staging.cravat.us/rest/service/query?mutation=chr"+SNPinfo.Chr+"_"+SNPinfo.Loc+"_"+SNPinfo.Strand+"_"+SNPinfo.Ref+"_"+SNPinfo.Alt 
        });
}

function switchSample( sampleName ){
  // Switch the sample label to the given sample.
  // Switch the samples in the browser.
  activeSample.innerText = "Sample: " + sampleName

  // Create browser
  if(! igv.browser ){
    createBrowser( inspectorView );
  } else {
    // Switch sample
    // Remove old tracks
    igv.browser.removeTrack( igv.browser.trackViews[2].track)
    igv.browser.removeTrack( igv.browser.trackViews[2].track)
    // Add new tracks
    sampleTrackOptionsRNA = { url: inspectorView[ sampleName ].RNA,
                            label: "RNA-Seq Sequencing (" + sampleName + ")",
                            height: 150,
                            color: "rgb( 102, 153, 255 )" }
    sampleTrackOptionsDNA = { url: inspectorView[ sampleName ].DNA,
                            label: "exome Sequencing (" +sampleName + ")",
                            height: 150 }
    // Update the Genome
    igv.browser.loadTrack( sampleTrackOptionsDNA );
    igv.browser.loadTrack( sampleTrackOptionsRNA );
  }
  // Update True Positive 
  // Update False Positive
  // Update False Negative
  updateErrorMenus( sampleName );
}
