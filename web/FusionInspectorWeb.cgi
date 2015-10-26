#!/usr/bin/python

print( "Content-Type: text/html\n\n" )
print( """<!DOCTYPE html>

<html lang="en">
<head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1, maximum-scale=1, user-scalable=no">
    <meta name="description" content="Trintiy CRTAT visualization for SNVs in discovered RNA-Seq.">
    <meta name="author" content="Brian Haas,Timothy Tickle">
    <title>Trinity CTAT Fusion Inspector</title>

    <!-- CSS -->
    <!-- jQuery UI CSS -->
    <link rel="stylesheet" type="text/css" href="//ajax.googleapis.com/ajax/libs/jqueryui/1.11.3/themes/smoothness/jquery-ui.css"/>
    <!-- Bootstrap CSS -->
    <link rel="stylesheet" type="text/css" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/css/bootstrap.min.css">
    <!-- Font Awesome CSS -->
    <link rel="stylesheet" type="text/css" href="//maxcdn.bootstrapcdn.com/font-awesome/4.2.0/css/font-awesome.min.css">
    <!-- IGV CSS -->
    <link rel="stylesheet" type="text/css" href="//igv.org/web/beta/igv.css">
    <!-- inspector css -->
    <link rel="stylesheet" type="text/css" href="css/inspector.css">
    <!-- Spinner from http://www.css-spinners.com/spinner/spinner -->
    <link rel="stylesheet" href="css/spinner.css">
    <!-- Associated with the Data Table -->
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.8/css/dataTables.bootstrap.min.css">

    <!-- Scripts -->
    <!-- jQuery JS -->
    <script type="text/javascript" src="//ajax.googleapis.com/ajax/libs/jquery/1.11.3/jquery.min.js"></script>
    <script type="text/javascript" src="//ajax.googleapis.com/ajax/libs/jqueryui/1.11.3/jquery-ui.min.js"></script>
    <!-- IGV JS-->
    <script type="text/javascript" src="//igv.org/web/beta/igv-all.min.js"></script>
    <!-- Data Table -->
    <script type="text/javascript" src="https://cdn.datatables.net/1.10.8/js/jquery.dataTables.min.js"></script>
    <script type="text/javascript" src="https://cdn.datatables.net/1.10.8/js/dataTables.bootstrap.min.js"></script>
    <!-- Bootstrap JS -->
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.4/js/bootstrap.min.js"></script>
</head>    
<body>
    <!-- Header (info) -->
    <div class="constainer-fluid" id="sampleHeader" style="background-color: #E7E7EF">
        <nav class="navbar navbar-default" style="background-color: #E7E7EF">
            <div class="col-xs-3" id="sampleId"><p class="navbar-text"><b>Sample:</b></p></div>
            <div class="col-xs-3" id="FusionNameDetail"><p class="navbar-text"><b>Fusion Name:</b> Not Selected</p></div>
            <div class="col-xs-3" id="FusionJunctionDetail"><p class="navbar-text"><b>Junction Read Count:</b> Not Selected</p></div>
            <div class="col-xs-3" id="FusionSpanningDetail"><p class="navbar-text"><b>Spanning Read Count:</b> Not Selected</p></div>
        </nav>
    </div>
    <!-- End fusion details -->
    <hr>

    <!-- Start tabs -->
    <!-- Start tabs header -->
    <ul id="tabDescription" class="nav nav-tabs">
      <li role="presentation" id="tabBrowser_tab" class="active"><a href="#tabBrowser" data-toggle="tab">Browse All Fusions</a></li>
      <li role="presentation" id="igvTab"><a href="#igvBrowser" data-toggle="tab">IGV Detail</a></li>
    </ul>
    <!-- End tabs header -->

    <!-- Start tabs content -->
    <div class="tab-content" id="tabContent">
      <!-- Start Data Table Tab -->
      <div role="tabpanel" id="tabBrowser" class="tab-pane fade in">
        <!-- Start data table -->
          <div class="table-responsive">
            <table id="fusionTable" class="table table-striped table-bordered table-hover active" cell spacing="0" width="100%"></table>
          </div>
        <!-- End data table -->
      </div>
      <div role="tabpanel" id="igvBrowser" class="tab-pane fade"></div>
      <!-- End Browser -->
    </div>
    <!-- End tab content -->

    <!-- Scripts -->
    <script src='fusion_inspector.json'></script>
    <script src='js/FusionInspector.js'></script>
    <script>
    $(document).ready(function() {
      // Load data (MOCKED)
      fusionInspectorState.cache[ "json" ] = fusionInspector;
      // Set sample name in header
      setSampleName( fusionInspectorState.cache.json );
      // Load the data table
      loadFusionDataTable();
      fusionInspectorState.cache.fusionTable = $('#fusionTable').DataTable({
          'order': [[ 0, 'asc' ]],
          'scrollX': true
      });
      $('#fusionTable tbody').on('click', 'tr', function() {
        curFusionRow = fusionInspectorState.cache.fusionTable.row( this ).data();
        setFusionDetails( getFusionAnnotationFromRow( 'Fusion', curFusionRow ),
                          getFusionAnnotationFromRow( 'Junction Reads', curFusionRow ),
                          getFusionAnnotationFromRow( 'Spanning Fragments', curFusionRow ) );
        goToFusion( getFusionAnnotationFromRow( 'Fusion', curFusionRow ),
                    getFusionAnnotationFromRow( 'Right Pos', curFusionRow ),
                    getFusionAnnotationFromRow( 'Left Pos', curFusionRow ));
      });
      loadIGVBrowser();
      // IGV browser has to be visible when the files are loaded.
      // If it is hidden the files load as 200 (full file) as opposed
      // to 206, whichis needed for indexed reading as igv web needs it.
      $('.nav-tabs a[href="#igvBrowser"]').tab('show');
      // This hooks into the event fired off by tabs being selected.
      // It forces a redraw of the tab. Because the data table is originally
      // draw in a hidden (height = 0) div, the table is misdrawn. You have to
      // Trigger a redraw when the tab is visible so the height of the data table
      // can be correctly calculated.
      // Thanks to http://stackoverflow.com/questions/20705905/bootstrap-3-jquery-event-for-active-tab-change
      $('a[data-toggle="tab"]').on('shown.bs.tab', function (e) {
        $('#fusionTable').DataTable().columns.adjust().draw();
      })
    });
    </script>
</body>
</html> """ )
