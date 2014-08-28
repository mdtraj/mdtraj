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
    this.WIDTH = this.$el.width() * this.aaScale;
    this.HEIGHT = this.$el.height() * this.aaScale;
    this.ASPECT = this.WIDTH / this.HEIGHT;
    this.NEAR = 1;
    this.FAR = 800;
    this.CAMERA_Z = 500;

    /* setup the camera */
    this.camera = new THREE.PerspectiveCamera(50, this.ASPECT, 1, 10000);
    this.camera.position.z = this.CAMERA_Z;
    this.camera.lookAt(new TV3(0, 0, 0));
    
    /* set up the renderer */
    this.renderer = new THREE.WebGLRenderer({antialiased: true});
    this.renderer.setClearColor(0xffffff, 1 );
    this.renderer.domElement.style.width = "100%";
    this.renderer.domElement.style.height = "100%";
    this.renderer.setClearColor( 0xffffff, 1 );
    this.renderer.setSize(this.WIDTH, this.HEIGHT);
    
    this.initializeScene();
    this.initializeLights();
    this.initializeMesh();
    this.scene.add(this.camera);

    this.renderer.render(this.scene, this.camera);
    this.$el.append(this.renderer.domElement);
};

RMol.prototype.initializeScene = function() {
    // CHECK: Should I explicitly call scene.deallocateObject?
    this.scene = new THREE.Scene();
    this.scene.fog = new THREE.Fog(this.bgColor, 100, 200);
    this.modelGroup = new THREE.Object3D();
    this.rotationGroup = new THREE.Object3D();
    this.rotationGroup.add(this.modelGroup);
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

/* Return */
console.log("Loaded RMol");
return RMol;
});