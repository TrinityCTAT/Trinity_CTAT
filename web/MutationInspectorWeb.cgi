#!/usr/bin/python

print( "Content-Type: text/html\n\n" )
print( """<!DOCTYPE html>

<html lang="en">
<head>
  <meta charset="utf-8">
  <meta http-equiv="X-UA-Compatible" content="IE=edge">
  <meta name="viewport" content="width=device-width, initial-scale=1, maximum-scale=1, user-scalable=no">
  <meta name="description" content="Trinity CTAT visualization for SNVs in discovered RNA-Seq.">
  <meta name="author" content="Timothy Tickle">
  <title>Trinity CTAT Mutation Inspector</title>

  <!-- CSS -->
  <!-- Bootstrap -->
  <link rel="stylesheet" type="text/css" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/css/bootstrap.min.css">
  <!-- jQuery UI CSS -->
  <link rel="stylesheet" type="text/css" href="http://ajax.googleapis.com/ajax/libs/jqueryui/1.11.3/themes/smoothness/jquery-ui.css"/>
  <!-- Font Awesome CSS -->
  <link rel="stylesheet" type="text/css" href="http://maxcdn.bootstrapcdn.com/font-awesome/4.2.0/css/font-awesome.min.css">
  <!-- IGV CSS -->
  <link rel="stylesheet" type="text/css" href="http://www.broadinstitute.org/igv/projects/igv-web/css/igv.css">
  <!-- Spinner from http://www.css-spinners.com/spinner/spinner -->
  <link rel="stylesheet" href="css/spinner.css">
  <!-- Associated with the Data Table -->
  <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.8/css/dataTables.bootstrap.min.css">
  <!-- Specific CSS to inspectors -->
  <link rel="stylesheet" type="text/css" href="css/inspector.css">

  <!-- Scripts -->
  <!-- jQuery JS -->
  <script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.11.3/jquery.min.js"></script>
  <script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jqueryui/1.11.3/jquery-ui.min.js"></script>
  <!-- Bootstrap -->
  <script type="text/javascript" src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.4/js/bootstrap.min.js"></script>
  <!-- IGV JS -->
  <script type="text/javascript" src="http://www.broadinstitute.org/igv/projects/igv-web/dist/igv-all.js"></script>
  <!-- Data Table -->
  <script type="text/javascript" src="https://cdn.datatables.net/1.10.8/js/jquery.dataTables.min.js"></script>
  <script type="text/javascript" src="https://cdn.datatables.net/1.10.8/js/dataTables.bootstrap.min.js"></script>

</head>
<body>
  <!-- Header (info) -->
  <div class="container-fluid" id="sampleHeader">
    <div class="row" style="background-color: #E7E7EF">
      <div class="col-xs-6 col-md-3" id="activeSample"><p><b>Sample: </b>None</p></div>
      <div class="col-xs-6 col-md-2"><p><b>Chr: </b><span id="currentChr">NA</span></p></div>
      <div class="col-xs-6 col-md-2"><p><b>Position: </b><span id="currentPosition">NA</span></p></div>
      <div class="col-xs-6 col-md-1"><p><b>Ref: </b><span id="currentRef">NA</span></p></td></div>
      <div class="col-xs-6 col-md-1"><p><b>Alt: </b><span id="currentAlt">NA</span></p></td></div>
      <div class="col-xs-6 col-md-1"><b>MuPIT: </b></div>
      <div class="col-xs-6 col-md-2" id="currentMupit">Not Available</div>
    </div>
  </div>
  <hr>
  <!-- End Header (info) -->

  <!-- IGV browser area -->
  <div id="igvBrowser">
  </div>
  <hr>
  <!-- End IGV browser area -->

  <!-- Start tabs -->
    <!-- Start tabs header -->
    <ul id="tabDescription" class="nav nav-tabs">
      <li id="tabBrowser_tab" class="active"><a href="#tabBrowser" data-toggle="tab">Browse All</a></li>
    </ul>
    <!-- End tabs header -->

    <!-- Start tabs content -->
    <div class="tab-content" id="tabContent">
      <!-- Start Browser -->
      <div role="tabpanel" id="tabBrowser" class="tab-pane fade in active">
        <div class="container-fluid" id="igvDiv" style="padding:5px; border:1px solid lightgray">
          <!-- Start data table -->
            <table id="mutationTable" class="table table-striped table-bordered table-hover" cell spacing="0" width="100%">
            </table>
          <!-- End data table -->
        </div>
      </div>
      <!-- End Browser -->
    </div>
    <!-- End tab content -->

  <!-- Scripts -->
  <script src='js/MutationInspectorWeb.js'></script>
  <script src='cancer.json'></script>
  <script>

    // Add entries to the table
    // Load Mutation table
    $(document).ready(function() {
        // Load mutation table
        var mutationInspector = loadMutationTable( "mutations.json" );
        var mutationTable = $('#mutationTable').DataTable({
          'order': [[ 0, 'asc' ]]
        });
        // Add click events to the table rows
        $('#mutationTable tbody').on('click', 'tr', function() {
          curMutationRow = mutationTable.row( this ).data()
          goToSNP( curMutationRow[ mutationInspectorView.chrKey ], curMutationRow[ mutationInspectorView.posKey ] );
          updateSNPInfo( curMutationRow[ mutationInspectorView.chrKey ],
                          curMutationRow[ mutationInspectorView.posKey ],
                          curMutationRow[ mutationInspectorView.refKey ],
                          curMutationRow[ mutationInspectorView.altKey ] );
          addSpecificTab( curMutationRow[ mutationInspectorView.chrKey ],
                          curMutationRow[ mutationInspectorView.posKey ],
                          curMutationRow[ mutationInspectorView.refKey ],
                          curMutationRow[ mutationInspectorView.altKey ] );
        })

        // Update Sample name
        $('#activeSample').html( '<p><b>Sample: </b>'+ mutationInspector.json.SAMPLE+'</p>' );
        // Load browser
        var igvBrowser = createIGVBrowser( mutationInspector.json );

        // Add state change event for default tab
        registerDefaultTabClick( "tabBrowser_tab" );
    } )

  </script>
</body>
</html>""" )
