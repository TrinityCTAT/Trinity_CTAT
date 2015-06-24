#!/usr/bin/python

print( "Content-Type: text/html\n\n" )
print( """<!DOCTYPE html>

<html lang="en">
<head>
  <meta charset="utf-8">
  <meta http-equiv="X-UA-Compatible" content="IE=edge">
  <meta name="viewport" content="width=device-width, initial-scale=1, maximum-scale=1, user-scalable=no">
  <meta name="description" content="Trinity CTAT visualization for error classes in exome to RNA-Seq comparisons">
  <meta name="author" content="Timothy Tickle">
  <title>Trinity CTAT Mutation Inspector</title>

  <!-- CSS -->
  <!-- Bootstrap -->
  <link rel="stylesheet" type="text/css" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.4/css/bootstrap.min.css">
  <!-- jQuery UI CSS -->
  <link rel="stylesheet" type="text/css" href="http://ajax.googleapis.com/ajax/libs/jqueryui/1.11.2/themes/smoothness/jquery-ui.css"/>
  <!-- Font Awesome CSS -->
  <link rel="stylesheet" type="text/css" href="http://maxcdn.bootstrapcdn.com/font-awesome/4.2.0/css/font-awesome.min.css">
  <!-- IGV CSS -->
  <link rel="stylesheet" type="text/css" href="http://www.broadinstitute.org/igv/projects/igv-web/css/igv.css">
  <!-- Spinner from http://www.css-spinners.com/spinner/spinner -->
  <link rel="stylesheet" href="css/spinner.css">

  <!-- Scripts -->
  <!-- jQuery JS -->
  <script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js"></script>
  <script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jqueryui/1.11.2/jquery-ui.min.js"></script>
  <!-- Bootstrap -->
  <script type="text/javascript" src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.4/js/bootstrap.min.js"></script>
  <!-- IGV JS -->
  <script type="text/javascript" src="http://www.broadinstitute.org/igv/projects/igv-web/dist/igv-all.js"></script>

</head>
<body>
  <!-- nav bar -->
  <ul class="nav nav-pills" style="background-color: #E0E0E0">
    <li role="presentation" class="dropdown">
      <a class="dropdown-toggle" data-toggle="dropdown" href="#" role="button" aria-expanded="false">
        Samples <span class="caret"></span>
      </a>
      <ul id="sampleMenu" class="dropdown-menu" role="menu">
      </ul>
    </li>

    <li role="presentation" class="dropdown">
      <a class="dropdown-toggle" data-toggle="dropdown" href="#" role="button" aria-expanded="false">
        True Positive <span class="caret"></span>
      </a>
      <ul id="truePositive" class="dropdown-menu" role="menu">
      </ul>
    </li>

    <li role="presentation" class="dropdown">
      <a class="dropdown-toggle" data-toggle="dropdown" href="#" role="button" aria-expanded="false">
        False Positive <span class="caret"></span>
      </a>
      <ul id="falsePositive" class="dropdown-menu" role="menu">
      </ul>
    </li>

    <li role="presentation" class="dropdown">
      <a class="dropdown-toggle" data-toggle="dropdown" href="#" role="button" aria-expanded="false">
        False Negative <span class="caret"></span>
      </a>
      <ul id="falseNegative" class="dropdown-menu" role="menu">
      </ul>
    </li>

  </ul>
  <!-- End nav bar -->

  <!-- Header (info) -->
  <div class="container-fluid" id="sampleHeader">
    <div class="page-header">
      <h2>Visualize Error Classes for RNA-Seq Validation Runs</h2>
      <h3 id="activeSample">Sample: CW104</h3>
    </div>
  </div>
  <div class="container-fluid" id="annotation">
    <div class="row">
      <div class="col-xs-3" id="sampleInfo">
        <p id="snpLocation"><b>SNP Location:</b> None Selected</p>
        <p id="snpCoverage"><b>SNP Coverage:</b> None Selected</p>
        <p id="snpRef"><b>Ref:</b> None Selected</p>
        <p id="snpAlt"><b>Alt:</b> None Selected</p>
        <p id="snpClass"></p>
      </div>
      <div id="annotationTableDiv" class="col-xs-9">
      </div>
    </div>
    <!-- End Annotations -->
  </div>
  <!-- End Header (info) -->

  <!-- Start tabs header -->
  <ul class="nav nav-tabs">
    <li role="presentation" class="active"><a data-toggle="tab" href="#tabBrowser">Browse</a></li>
    <li role="presentation"><a data-toggle="tab" href="#tabCravat">Annotation</a></li>
  </ul>
  <!-- End tabs header -->

  <!-- Start tabs content -->
  <div class="tab-content">
    <!-- Start Bowser -->
    <div id="tabBrowser" class="tab-pane fade in active">
      <div class="container-fluid" id="igvDiv" style="padding:5px; border:1px solid lightgray"></div>
    </div>
    <!-- End Browser -->
    <!-- Annotations -->
    <div id="tabCravat" class="tab-pane fade in">
      <div class="constainer-fluid" id="cravatTab" style="padding:5px; border:1px solid lightgray">
        <p>Not Specified</p>
      </div>
    </div>
    <!-- End Annotations -->
  </div>
  <!-- End tab content -->

  <!-- Scripts -->
  <script src='js/MutationInspectorWeb.js'></script>
  <script>

    // Read in data from file
    // Build website from file
    inspectorView = null
    readData( "pipeline_inspector.json" );

  </script>
</body>
</html>""" )
