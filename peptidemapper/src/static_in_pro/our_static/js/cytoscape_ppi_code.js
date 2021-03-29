Promise.all([
])
  .then(function(dataArray) {

    var cy = window.cy = cytoscape({
      container: document.getElementById('cy'),
      style: cystylejson,
      elements: ppidata
    });
    var params = {
      name: 'cola',
      nodeSpacing: 5,
      edgeLengthVal: 3,
      edgeScore:0,
      animate: true,
      randomize: false,
      maxSimulationTime: 1500
    };
    var layout = makeLayout();
    var running = false;

    cy.on('layoutstart', function(){
      running = true;
    }).on('layoutstop', function(){
      running = false;
    });

    layout.run();

    var $config = $('#config');
    var $btnParam = $('<div class="param"></div>');
    $config.append( $btnParam );

    var sliders = [
      {
        label: 'Edge length',
        param: 'edgeLengthVal',
        min: 1,
        max: 200
      },

      {
        label: 'Node spacing',
        param: 'nodeSpacing',
        min: 1,
        max: 50
      },
      {
        label: 'Score',
        param: 'edgeScore',
        min: 0,
        max: 1
      }
    ];

    var buttons = [
      {
        label: '<i class="fa fa-random"></i>',
        layoutOpts: {
          'name': 'cola',
          randomize: true,
          flow: null
        }
      },
      {
        label: '<i class="fa fa-circle-o"></i>',
        layoutOpts: {
          'name': 'circle'
        }
      },
      {
        label: '<i class="fa fa-th"></i>',
        layoutOpts: {
          'name': 'grid'
        }
      },
      {
        label: '<i class="fa fa-dot-circle-o"></i>',
        layoutOpts: {
          'name': 'concentric'
        }
      }
    ];
    var intactType=new Set();
    cy.edges().forEach(function( edg ) {
      intactType.add(edg.data('group'));
    });
    intactType=[...intactType];
    var checkboxes=[];
    for(var i=0; i<intactType.length; i++) {
      checkboxes.push(
      {
        label: intactType[i],
        layoutOpts: {
          randomize: true,
          flow: null
        }
      }
      )
    }
    sliders.forEach( makeSlider );

    buttons.forEach( makeButton );

    $config.append( 'Evidence suggesting a functional link' );
    checkboxes.forEach( makeCheckbox );

    function capitalizeFirstLetter(string) {
      return string.charAt(0).toUpperCase() + string.slice(1);
    }

    function makeLayout( opts ){
      params.randomize = false;
      params.edgeLength = function(e){ return params.edgeLengthVal / e.data('weight'); };

      for( var i in opts ){
        params[i] = opts[i];
      }

      return cy.layout( params );
    }

    function makeSlider( opts ){
      var $input = $('<input></input>');
      var $param = $('<div class="param"></div>');

      $param.append('<span class="label label-default">'+ opts.label +'</span>');
      $param.append( $input );
      $config.append( $param );
      var step=1;
      if (opts.label == 'Score'){
        step=0.01
      }
      var p = $input.slider({
        min: opts.min,
        max: opts.max,
        value: params[ opts.param ],
        step: step
      }).on('slide', _.throttle( function(){
        params[ opts.param ] = p.getValue();
        if (opts.label == 'Score'){
          makeScoreSlider(p.getValue());
        } else {
          layout.stop();
          layout = makeLayout();
          layout.run();
        }
      }, 16 ) ).data('slider');
    }

    function makeButton( opts ){
      var $button = $('<button class="btn btn-default">'+ opts.label +'</button>');

      $btnParam.append( $button );

      $button.on('click', function(){
        layout.stop();

        if( opts.fn ){ opts.fn(); }

        layout = makeLayout( opts.layoutOpts );
        layout.run();
      });
    }

    function makeCheckbox( opts ){
      var $checkbox = $('<input class="my-check" type="checkbox" checked />');
      var $checkParam = $('<div class="param"></div>');
      $checkParam.append( $checkbox );
      $checkParam.append( '&nbsp;'+capitalizeFirstLetter(opts.label ));
      if(opts.label == "combined"){
        $checkParam.append('&nbsp;'+'<i class="fa fa-circle" style="color:#D0D0D0;"></i>');
      }
      if (opts.label == "textmining"){
        $checkParam.append('&nbsp;'+'<i class="fa fa-circle" style="color:#dad4a2;"></i>');
      }
      if (opts.label == "database"){
        $checkParam.append('&nbsp;'+'<i class="fa fa-circle" style="color:#f6c384;"></i>');
      }
      if (opts.label == "experimental"){
        $checkParam.append('&nbsp;'+'<i class="fa fa-circle" style="color:#eaa2a2;"></i>');
      }
      if (opts.label == "phylogenetic profile"){
        $checkParam.append('&nbsp;'+'<i class="fa fa-circle" style="color:#9bd8de;"></i>');
      }
      if (opts.label == "gene fusion"){
        $checkParam.append('&nbsp;'+'<i class="fa fa-circle" style="color:#90e190;"></i>');
      }
      if (opts.label == "gene neighborhood"){
        $checkParam.append('&nbsp;'+'<i class="fa fa-circle" style="color:#a0b3dc;"></i>');
      }
      if (opts.label == "coexpression"){
        $checkParam.append('&nbsp;'+'<i class="fa fa-circle" style="color:#d0b7d5;"></i>');
      }
      $config.append( $checkParam );

      var removedData;
      $checkbox.on('click', function(){
        var checkBoxStatus=$checkbox.is(':checked');
        var edgeGroup=opts.label;
        if (checkBoxStatus){
           removedData.restore();
           cy.add(removedData);
         } else {
           removedData=cy.remove('edge[group=\'' + edgeGroup + '\']');
         }
      });

    }

    var filteredEdges;
    function makeScoreSlider( sliderVal){
      if (typeof filteredEdges !== 'undefined'){
        filteredEdges.restore();
        cy.add(filteredEdges);
      }
      filteredEdges =cy.remove('edge[weight < ' + sliderVal + ']');
    }
    cy.nodes().forEach(function(n){
      var ncolor = n.data('nodecolor');
      if (ncolor !== 'skyblue') {
        var g = n.data('name');
        var orgid = n.data('organismid');
        var hrf = n.data('href');
        var uni = n.data('Uniprotkb');
        n.qtip({
          content: [
            {
              name: 'MRMAssayDB search',
              url: hrf+'UniProtKB Accession='+ uni
            },
            {
              name: 'UniProt search',
              url: 'https://www.uniprot.org/uniprot/'+ uni
            }
          ].map(function( link ){
            return '<a target="_blank" href="' + link.url + '">' + link.name + '</a>';
          }).join('<br />\n'),
          position: {
            my: 'top center',
            at: 'bottom center'
          },
          style: {
            classes: 'qtip-bootstrap',
            tip: {
              width: 16,
              height: 8
            }
          }
        });
      }
    });

   // Find its connected neighbours.
  cy.edges().forEach(function( edg ) {
    var node1;
    var node2;
    var x;
    var y;
    cy.nodes().forEach(function(nd){
       if (edg.data('source') == nd.data('id')){
         node1=nd.data('name')
       }

    });
    cy.nodes().forEach(function(nd){

       if (edg.data('target') == nd.data('id')){
         node2=nd.data('name')
       }
    });
    x= (edg._private.bodyBounds.x1+edg._private.bodyBounds.x2)/2;
    y= (edg._private.bodyBounds.y1+edg._private.bodyBounds.y2)/2;
    edg.qtip({
        content: [
          {
            name: 'Type of Interaction:'+edg.data('group'),
            info: 'Gene1='+ node1 +'&Gene2='+node2
          }
        ].map(function( d ){
          return  d.name + '<br />\n' + d.info;
        }).join('<br />\n'),
        position: {
          my: 'center center',
          at: 'top center',
          adjust:{
            y:150,
            x:20
          }
        },
        style: {
          classes: 'qtip-bootstrap',
          tip: {
            width: 16,
            height: 8
          }
        }
    });
  });
    $('#config-toggle').on('click', function(){
      $('body').toggleClass('config-closed');

      cy.resize();
    });

  })
  .then(function() {
    FastClick.attach( document.body );
  });
