define(["three"], function(THREE) {
    var TV3 = THREE.Vector3, TF3 = THREE.Face3, TCo = THREE.Color;

    function RMol ($el) {
        this.create($el);
        return true;
    }
    
/**
 * Constructor
 */
RMol.prototype.create = function($el) {
    this.$el = $el;
    this.aaScale = 1;
    // this.WIDTH = this.$el.width() * this.aaScale;
    // this.HEIGHT = this.$el.height() * this.aaScale;   
    this.WIDTH = 300;
    this.HEIGHT = 300;
    this.ASPECT = this.WIDTH / this.HEIGHT;
    this.NEAR = 1;
    this.FAR = 800;

    /* set up the camera */
    this.camera = new THREE.PerspectiveCamera(50, this.ASPECT, 1, 10000);
    this.CAMERA_Z = -150;
    this.camera.position = new TV3(0, 0, this.CAMERA_Z);
    this.camera.lookAt(new TV3(0, 0, 0));

    /* set up the renderer */
    // this.renderer = new THREE.WebGLRenderer({antialias: true});
    this.renderer = new THREE.WebGLRenderer();
    this.renderer.setClearColor(0xffffff, 1 );
    this.renderer.domElement.style.width = "100%";
    this.renderer.domElement.style.height = "100%";
    this.renderer.setSize(this.WIDTH, this.HEIGHT);

    /* */
    this.scene = new THREE.Scene();
    var geometry = new THREE.BoxGeometry(200, 200, 200);
    var material = new THREE.MeshNormalMaterial();
    var mesh = new THREE.Mesh(geometry, material);
    mesh.rotation.x += 0.5;
    mesh.rotation.y += 0.5;
    this.scene.add(mesh);
    
    // this.initializeScene();
    // // this.initializeLights();
    // this.initializeMesh();

    this.renderer.render(this.scene, this.camera);
    this.$el.append(this.renderer.domElement);
};


RMol.prototype.initializeScene = function() {
   // CHECK: Should I explicitly call scene.deallocateObject?
   this.scene = new THREE.Scene();
   // this.scene.fog = new THREE.Fog(this.bgColor, 100, 200);
   //
   // this.modelGroup = new THREE.Object3D();
   // this.rotationGroup = new THREE.Object3D();
   // this.rotationGroup.add(this.modelGroup);
   // this.scene.add(this.rotationGroup);
};

RMol.prototype.initializeLights = function() {
   var directionalLight =  new THREE.DirectionalLight(0xFFFFFF);
   directionalLight.position = new TV3(0.2, 0.2, -1).normalize();
   directionalLight.intensity = 1.2;
   this.scene.add(directionalLight);
   var ambientLight = new THREE.AmbientLight(0x202020);
   this.scene.add(ambientLight);
};

RMol.prototype.initializeMesh = function() {
    var geometry = new THREE.BoxGeometry(200, 200, 200);
    var material = new THREE.MeshNormalMaterial();
    var mesh = new THREE.Mesh(geometry, material);
    mesh.rotation.x += 0.5;
    mesh.rotation.y += 0.5;
    this.scene.add(mesh);
};

/* export the newly created object */
console.log("RMol loaded");
return RMol;
});


//
// require(["widgets/js/widget",
//          "RMol.js",
//          "three"],
// function(WidgetManager, RMol, THREE){
//     var ThreeWidgetView = IPython.DOMWidgetView.extend({
//         render : function(){
//             var $el = $("<div>").css({
//                 width: '302px',
//                 height: '302px', border:'1px solid #ccc',
//                 margin: 'auto'}).append(
//                     $('<div id="inner">').css({
//                         position: 'absolute',
//                         width: '300px',
//                         height: '300px',
//                     })
//             );
//             this.setElement($el);
//             this.rmol = new RMol($el);
//         }
//     });
//
//     // Register the DatePickerView with the widget manager.
//     WidgetManager.register_widget_view('ThreeWidgetView', ThreeWidgetView);
// });