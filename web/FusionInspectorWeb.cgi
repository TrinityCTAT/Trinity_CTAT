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
    <!-- Bootstrap CSS -->
    <link rel="stylesheet" type="text/css" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/css/bootstrap.min.css">
    <!-- jQuery UI CSS -->
    <link rel="stylesheet" type="text/css" href="//ajax.googleapis.com/ajax/libs/jqueryui/1.11.3/themes/smoothness/jquery-ui.css"/>
    <!-- Font Awesome CSS -->
    <link rel="stylesheet" type="text/css" href="//maxcdn.bootstrapcdn.com/font-awesome/4.2.0/css/font-awesome.min.css">
    <!-- IGV CSS -->
    <link rel="stylesheet" type="text/css" href="//igv.org/web/beta/igv.css">
    <!-- inspector css -->
    <link rel="stylesheet" type="text/css" href="css/inspector.css">

    <!-- Scripts -->
    <!-- jQuery JS -->
    <script type="text/javascript" src="//ajax.googleapis.com/ajax/libs/jquery/1.11.3/jquery.min.js"></script>
    <script type="text/javascript" src="//ajax.googleapis.com/ajax/libs/jqueryui/1.11.3/jquery-ui.min.js"></script>
    <!-- Bootstrap JS -->
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.4/js/bootstrap.min.js"></script>
    <!-- IGV JS-->
    <script type="text/javascript" src="//igv.org/web/beta/igv-all.min.js"></script>
</head>    
<body>
    <!-- Header (info) -->
    <div class-"constainer-fluid" id="sampleHeader">
        <div class="row" style="background-color: #E7E7EF">
            <div class="col-xs-6" id="sampleId"><p>Sample:</p></div>
            <div class="col=xs-6">
                <ul id="fusionList">
                    <li>No fusions loaded</li>
                </ul>
            </div>
        </div>
    </div>
    <hr>
    <!--- End header )info) -->


    <!-- Start tabs header -->
    
    <!-- End tabs header -->
    <!-- Start tabs content -->

    <!-- End tabs content -->


    <!-- IGV browser area -->
    <div id="igvBrowser"></div>

    <!-- Scripts -->
    <script src='fusion_inspector.json'></script>
    <script src='js/FusionInspector.js'></script>
    <script>
    var igvBrowser = loadFusionIGV( fusionInspector )
    </script>
</body>
</html> """ )
