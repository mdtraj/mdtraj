// define(["three", "three/trackball"], function(THREE) {
define([], function() {

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
    this.CAMERA_Z = -150;

    /* setup the camera */
    this.camera = new THREE.PerspectiveCamera(20, this.ASPECT, this.NEAR, this.FAR);
    this.camera.position.set(0, 0, this.CAMERA_Z);
    this.camera.lookAt(new TV3(0, 0, 0));
    
    /* set up the renderer */
    this.renderer = new THREE.WebGLRenderer({antialiased: true});
    this.renderer.setClearColor(0xffffff, 1 );
    // this.renderer.setClearColor(0x000000, 1 );
    this.renderer.domElement.style.width = "100%";
    this.renderer.domElement.style.height = "100%";
    this.renderer.setClearColor( 0xffffff, 1 );
    this.renderer.setSize(this.WIDTH, this.HEIGHT);
    
    this.initializeScene();
    this.initializeLights();
    this.initializeDefaultValues();
    this.scene.add(this.camera);
    
    this.$el.append(this.renderer.domElement);
    this.enableMouse();
};

RMol.prototype.render = function() {
    this.renderer.render(this.scene, this.camera);
}

/** 
 * Initialization of the core elements
 */
RMol.prototype.initializeScene = function() {

    // CHECK: Should I explicitly call scene.deallocateObject?
    this.scene = new THREE.Scene();
    // this.scene.fog = new THREE.Fog(this.bgColor, 100, 200);
    this.modelGroup = new THREE.Object3D();
    this.rotationGroup = new THREE.Object3D();
    this.rotationGroup.add(this.modelGroup);
    
    // var axes = new THREE.AxisHelper(10);
    // this.rotationGroup.add(axes);
    
    this.scene.add(this.rotationGroup);
};

RMol.prototype.initializeLights = function() {
    var directionalLight =  new THREE.DirectionalLight(0xFFFFFF);
    directionalLight.position.set(0.2, 0.2, -1);
    directionalLight.intensity = 1.2;

    var ambientLight = new THREE.AmbientLight(0x202020);

    this.scene.add(directionalLight);
    this.scene.add(ambientLight);
};


RMol.prototype.initializeDefaultValues = function() {
    this.sphereRadius = 1.5; 
    this.cylinderRadius = 0.4;
    this.lineWidth = 1.5 * this.aaScale;
    this.curveWidth = 10 * this.aaScale;
    this.defaultColor = 0xCCCCCC;
    this.sphereQuality = 16; //16;
    this.cylinderQuality = 16; //8;
    this.axisDIV = 5; // 3 still gives acceptable quality
    this.strandDIV = 6;
    this.nucleicAcidStrandDIV = 4;
    this.tubeDIV = 8;
    this.coilWidth = 0.3;
    this.helixSheetWidth = 1.3;
    this.nucleicAcidWidth = 0.8;
    this.thickness = 0.4;
    
    this.ElementColors = {
        "H": 0xCCCCCC, "C": 0xAAAAAA, "O": 0xCC0000,
        "N": 0x0000CC, "S": 0xCCCC00, "P": 0x6622CC,
        "F": 0x00CC00, "CL": 0x00CC00, "BR": 0x882200, "I": 0x6600AA,
        "FE": 0xCC6600, "CA": 0x8888AA};
    this.vdwRadii = {
        "H": 1.2, "Li": 1.82, "Na": 2.27, "K": 2.75, "C": 1.7, "N": 1.55,
        "O": 1.52, "F": 1.47, "P": 1.80, "S": 1.80, "CL": 1.75, "BR": 1.85,
        "SE": 1.90, "ZN": 1.39, "CU": 1.4, "NI": 1.63};
};

/**
 * Set the system topology
 */
RMol.prototype.setTopology = function(topology) {
    var chains = topology.chains;
    var bonds = topology.bonds;
    var atoms = [];
    
    for (var i = 0; i < chains.length; i++) {
        var chain = chains[i];
        for (var j = 0; j < chain.residues.length; j++) {
            var residue = chain.residues[j];
            for (var k = 0; k < residue.atoms.length; k++) {
                var atom = chains[i].residues[j].atoms[k];
                atoms[atom.index] = {
                    'resn': residue.name,
                    'elem': atom.element,
                    'hetflag': false,
                    'chain': chain.index,
                    'resi': residue.index,
                    'serial': atom.index,
                    'atom': atom.name,
                    'ss': 'c',
                    'color': 0xFFFFFF,
                };
            }
        }
    }
    this.atoms = atoms;
};

RMol.prototype.setXYZ = function(xyz) {
    var natoms = xyz.length;
    for (var i = 0; i < natoms; i++) {
        this.atoms[i].x = 10 * xyz[i][0];
        this.atoms[i].y = 10 * xyz[i][1];
        this.atoms[i].z = 10 * xyz[i][2];
    }
}


/** 
 * Create all of the geometry elements based on a particular set
 * of representationOptions
 */ 
RMol.prototype.setRepresentation = function(representationOptions) {
    var all = this.getAllAtoms();
    // var hetatm = this.removeSolvents(this.getHetatms(all));
    var doNotSmoothen = false;

    this.colorByAtom(all, {});

    var asu = this.modelGroup;
    var colorMode = representationOptions.color;
    var mainchainMode = representationOptions.mainChain;
    var hetatmMode = representationOptions.heteroAtomsOptions;
    var projectionMode = representationOptions.projection;
    var sideChains = representationOptions.sideChains;

    if (colorMode == 'ss') {
        this.colorByStructure(all, 0xcc00cc, 0x00cccc);
    } else if (colorMode == 'chain') {
        this.colorByChain(all);
    } else if (colorMode == 'chainbow') {
        this.colorChainbow(all);
    } else if (colorMode == 'b') {
        this.colorByBFactor(all);
    } else if (colorMode == 'polarity') {
        this.colorByPolarity(all, 0xcc0000, 0xcccccc);
    }

   //  if (mainchainMode == 'ribbon') {
   //      this.drawCartoon(asu, all, doNotSmoothen);
   //      this.drawCartoonNucleicAcid(asu, all);
   //  } else if (mainchainMode == 'thickRibbon') {
   //      this.drawCartoon(asu, all, doNotSmoothen, this.thickness);
   //      // this.drawCartoonNucleicAcid(asu, all, null, this.thickness);
   //  } else if (mainchainMode == 'strand') {
   //      this.drawStrand(asu, all, null, null, null, null, null, doNotSmoothen);
   //      this.drawStrandNucleicAcid(asu, all);
   //  } else if (mainchainMode == 'chain') {
   //      this.drawMainchainCurve(asu, all, this.curveWidth, 'CA', 1);
   //      this.drawMainchainCurve(asu, all, this.curveWidth, 'O3\'', 1);
   //  } else if (mainchainMode == 'cylinderHelix') {
   //      this.drawHelixAsCylinder(asu, all, 1.6);
   //      this.drawCartoonNucleicAcid(asu, all);
   //  } else if (mainchainMode == 'tube') {
   //      this.drawMainchainTube(asu, all, 'CA');
   //      this.drawMainchainTube(asu, all, 'O3\''); // FIXME: 5' end problem!
   //  } else if (mainchainMode == 'bonds') {
   //      this.drawBondsAsLine(asu, all, this.lineWidth);
   //  }
   //
   // if (hetatmMode == 'stick') {
   //    this.drawBondsAsStick(target, hetatm, this.cylinderRadius, this.cylinderRadius, true);
   // } else if (hetatmMode == 'sphere') {
   //     this.drawAtomsAsSphere(target, hetatm, this.sphereRadius);
   // } else if (hetatmMode == 'line') {
   //     this.drawBondsAsLine(target, hetatm, this.curveWidth);
   // } else if (hetatmMode == 'icosahedron') {
   //     this.drawAtomsAsIcosahedron(target, hetatm, this.sphereRadius);
   // } else if (hetatmMode == 'ballAndStick') {
   //     this.drawBondsAsStick(target, hetatm, this.cylinderRadius / 2.0, this.cylinderRadius, true, false, 0.3);
   // } else if (hetatmMode == 'ballAndStick2') {
   //     this.drawBondsAsStick(target, hetatm, this.cylinderRadius / 2.0, this.cylinderRadius, true, true, 0.3);
   // }
   //
   // if (sideChains == 'line') {
   //     this.drawBondsAsLine(this.modelGroup, this.getSidechains(all), this.lineWidth);
   // }
   
   this.drawAtomsAsSphere(this.modelGroup, this.getAllAtoms(), this.sphereRadius);
   

   // if (projectionMode == 'perspective') {
   //     this.camera = this.perspectiveCamera;
   // } else if (projectionMode == 'orthoscopic') {
   //     this.camera = this.orthoscopicCamera;
   // }

   // this.setBackground(representationOptions.backgroundColor);
};


RMol.prototype.zoomInto = function(atomlist, keepSlab) {
    var tmp = this.getExtent(atomlist);
    var center = new TV3(tmp[2][0], tmp[2][1], tmp[2][2]);
    this.modelGroup.position = center.multiplyScalar(-1);

    var x = tmp[1][0] - tmp[0][0];
    var y = tmp[1][1] - tmp[0][1];
    var z = tmp[1][2] - tmp[0][2];

    var maxD = Math.sqrt(x * x + y * y + z * z);
    if (maxD < 25) maxD = 25;

    if (!keepSlab) {
      this.slabNear = -maxD / 1.9;
      this.slabFar = maxD / 3;
    }

    this.rotationGroup.position.z = maxD * 0.35 / Math.tan(Math.PI / 180.0 * this.camera.fov / 2) - 150;
    this.rotationGroup.quaternion = new THREE.Quaternion(1, 0, 0, 0);

};

/* ----------------------------------------------------------- */


RMol.prototype.getAllAtoms = function() {
   var ret = [];
   for (var i in this.atoms) {
      ret.push(this.atoms[i].serial);
   }
   return ret;
};

RMol.prototype.getExtent = function(atomlist) {
    var xmin = 9999, ymin = 9999, zmin = 9999;
    var xmax = -9999, ymax = -9999, zmax = -9999;
    var xsum = 0, ysum = 0, zsum = 0, cnt = 0;

    for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;
      cnt++;
      xsum += atom.x; ysum += atom.y; zsum += atom.z;

      xmin = (xmin < atom.x) ? xmin : atom.x;
      ymin = (ymin < atom.y) ? ymin : atom.y;
      zmin = (zmin < atom.z) ? zmin : atom.z;
      xmax = (xmax > atom.x) ? xmax : atom.x;
      ymax = (ymax > atom.y) ? ymax : atom.y;
      zmax = (zmax > atom.z) ? zmax : atom.z;
    }
    return [[xmin, ymin, zmin], [xmax, ymax, zmax], [xsum / cnt, ysum / cnt, zsum / cnt]];
};


/** 
 * Color the atoms. These routines modifty the color attribute of each
 * atom in the topology.
 */ 
RMol.prototype.colorByAtom = function(atomlist, colors) {
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      var c = colors[atom.elem];
      if (c == undefined) c = this.ElementColors[atom.elem];
      if (c == undefined) c = this.defaultColor;
      atom.color = c;
   }
};

RMol.prototype.colorChainbow = function(atomlist, colorSidechains) {
   var cnt = 0;
   var atom, i;
   for (i in atomlist) {
      atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if ((colorSidechains || atom.atom != 'CA' || atom.atom != 'O3\'') && !atom.hetflag)
         cnt++;
   }

   var total = cnt;
   cnt = 0;
   for (i in atomlist) {
      atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if ((colorSidechains || atom.atom != 'CA' || atom.atom != 'O3\'') && !atom.hetflag) {
         var color = new TCo(0);
         color.setHSL(240.0 / 360 * (1 - cnt / total), 1, 0.9);
         atom.color = color.getHex();
         cnt++;
      }
   }
};


/**
 * Drawing methods
 */

RMol.prototype.drawCartoon = function(group, atomlist, doNotSmoothen, thickness) {
   this.drawStrand(group, atomlist, 2, undefined, true, undefined, undefined, doNotSmoothen, thickness);
};

RMol.prototype.drawStrand = function(group, atomlist, num, div, fill, coilWidth, helixSheetWidth, doNotSmoothen, thickness) {
   num = num || this.strandDIV;
   div = div || this.axisDIV;
   coilWidth = coilWidth || this.coilWidth;
   doNotSmoothen == (doNotSmoothen == undefined) ? false : doNotSmoothen;
   helixSheetWidth = helixSheetWidth || this.helixSheetWidth;
   var points = []; for (var k = 0; k < num; k++) points[k] = [];
   var colors = [];
   var currentChain, currentResi, currentCA;
   var prevCO = null, ss=null, ssborder = false;

   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined) continue;

      if ((atom.atom == 'O' || atom.atom == 'CA') && !atom.hetflag) {
         if (atom.atom == 'CA') {
            if (currentChain != atom.chain || currentResi + 1 != atom.resi) {
               for (var j = 0; !thickness && j < num; j++)
                  this.drawSmoothCurve(group, points[j], 1 ,colors, div);
               if (fill) this.drawStrip(group, points[0], points[num - 1], colors, div, thickness);
               var points = []; for (var k = 0; k < num; k++) points[k] = [];
               colors = [];
               prevCO = null; ss = null; ssborder = false;
            }
            currentCA = new TV3(atom.x, atom.y, atom.z);
            currentChain = atom.chain;
            currentResi = atom.resi;
            ss = atom.ss; ssborder = atom.ssstart || atom.ssend;
            colors.push(atom.color);
         } else { // O
            var O = new TV3(atom.x, atom.y, atom.z);
            O.sub(currentCA);
            O.normalize(); // can be omitted for performance
            O.multiplyScalar((ss == 'c') ? coilWidth : helixSheetWidth); 
            if (prevCO != undefined && O.dot(prevCO) < 0) O.negate();
            prevCO = O;
            for (var j = 0; j < num; j++) {
               var delta = -1 + 2 / (num - 1) * j;
               var v = new TV3(currentCA.x + prevCO.x * delta, 
                               currentCA.y + prevCO.y * delta, currentCA.z + prevCO.z * delta);
               if (!doNotSmoothen && ss == 's') v.smoothen = true;
               points[j].push(v);
            }                         
         }
      }
   }
   for (var j = 0; !thickness && j < num; j++)
      this.drawSmoothCurve(group, points[j], 1 ,colors, div);
   if (fill) this.drawStrip(group, points[0], points[num - 1], colors, div, thickness);
};


RMol.prototype.drawStrip = function(group, p1, p2, colors, div, thickness) {
   if ((p1.length) < 2) return;
   div = div || this.axisDIV;
   p1 = this.subdivide(p1, div);
   p2 = this.subdivide(p2, div);

   if (!thickness) return this.drawThinStrip(group, p1, p2, colors, div);

   var geo = new THREE.Geometry();
   var vs = geo.vertices, fs = geo.faces;
   var axis, p1v, p2v, a1v, a2v;
   for (var i = 0, lim = p1.length; i < lim; i++) {
      vs.push(p1v = p1[i]); // 0
      vs.push(p1v); // 1
      vs.push(p2v = p2[i]); // 2
      vs.push(p2v); // 3
      if (i < lim - 1) {
         var toNext = p1[i + 1].clone().sub(p1[i]);
         var toSide = p2[i].clone().sub(p1[i]);
         axis = toSide.cross(toNext).normalize().multiplyScalar(thickness);
      }
      vs.push(a1v = p1[i].clone().add(axis)); // 4
      vs.push(a1v); // 5
      vs.push(a2v = p2[i].clone().add(axis)); // 6
      vs.push(a2v); // 7
   }
   var faces = [[0, 2, -6, -8], [-4, -2, 6, 4], [7, 3, -5, -1], [-3, -7, 1, 5]];
   for (var i = 1, lim = p1.length; i < lim; i++) {
      var offset = 8 * i, color = new TCo(colors[Math.round((i - 1)/ div)]);
      for (var j = 0; j < 4; j++) {
         var f = new THREE.Face4(offset + faces[j][0], offset + faces[j][1], offset + faces[j][2], offset + faces[j][3], undefined, color);
         fs.push(f);
      }
   }
   var vsize = vs.length - 8; // Cap
   for (var i = 0; i < 4; i++) {vs.push(vs[i * 2]); vs.push(vs[vsize + i * 2])};
   vsize += 8;
   fs.push(new THREE.Face4(vsize, vsize + 2, vsize + 6, vsize + 4, undefined, fs[0].color));
   fs.push(new THREE.Face4(vsize + 1, vsize + 5, vsize + 7, vsize + 3, undefined, fs[fs.length - 3].color));
   geo.computeFaceNormals();
   geo.computeVertexNormals(false);
   var material =  new THREE.MeshLambertMaterial();
   material.vertexColors = THREE.FaceColors;
   var mesh = new THREE.Mesh(geo, material);
   mesh.doubleSided = true;
   group.add(mesh);
};

RMol.prototype.drawAtomsAsSphere = function(group, atomlist, defaultRadius, forceDefault, scale) {

   var sphereGeometry = new THREE.SphereGeometry(1, this.sphereQuality, this.sphereQuality); // r, seg, ring

   for (var i = 0; i < atomlist.length; i++) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined) continue;

      var sphereMaterial = new THREE.MeshLambertMaterial({color: atom.color});
      var sphere = new THREE.Mesh(sphereGeometry, sphereMaterial);
      group.add(sphere);
      var r = (!forceDefault && this.vdwRadii[atom.elem] != undefined) ? this.vdwRadii[atom.elem] : defaultRadius;
      if (!forceDefault && scale) r *= scale;
      sphere.scale.x = sphere.scale.y = sphere.scale.z = r;
      sphere.position.x = atom.x;
      sphere.position.y = atom.y;
      sphere.position.z = atom.z;
      // console.log(sphere.position);
      
   }
};



// Catmull-Rom subdivision
RMol.prototype.subdivide = function(_points, DIV) { // points as Vector3
   var ret = [];
   var points = _points;
   points = new Array(); // Smoothing test
   points.push(_points[0]);
   for (var i = 1, lim = _points.length - 1; i < lim; i++) {
      var p1 = _points[i], p2 = _points[i + 1];
      if (p1.smoothen) points.push(new TV3((p1.x + p2.x) / 2, (p1.y + p2.y) / 2, (p1.z + p2.z) / 2));
      else points.push(p1);
   }
   points.push(_points[_points.length - 1]);

   for (var i = -1, size = points.length; i <= size - 3; i++) {
      var p0 = points[(i == -1) ? 0 : i];
      var p1 = points[i + 1], p2 = points[i + 2];
      var p3 = points[(i == size - 3) ? size - 1 : i + 3];
      var v0 = new TV3().subVectors(p2, p0).multiplyScalar(0.5);
      var v1 = new TV3().subVectors(p3, p1).multiplyScalar(0.5);
      for (var j = 0; j < DIV; j++) {
         var t = 1.0 / DIV * j;
         var x = p1.x + t * v0.x 
                  + t * t * (-3 * p1.x + 3 * p2.x - 2 * v0.x - v1.x)
                  + t * t * t * (2 * p1.x - 2 * p2.x + v0.x + v1.x);
         var y = p1.y + t * v0.y 
                  + t * t * (-3 * p1.y + 3 * p2.y - 2 * v0.y - v1.y)
                  + t * t * t * (2 * p1.y - 2 * p2.y + v0.y + v1.y);
         var z = p1.z + t * v0.z 
                  + t * t * (-3 * p1.z + 3 * p2.z - 2 * v0.z - v1.z)
                  + t * t * t * (2 * p1.z - 2 * p2.z + v0.z + v1.z);
         ret.push(new TV3(x, y, z));
      }
   }
   ret.push(points[points.length - 1]);
   return ret;
};


RMol.prototype.enableMouse = function() {
	var controls = new THREE.TrackballControls(this.camera);

	controls.rotateSpeed = 1.0;
	controls.zoomSpeed = 1.2;
	controls.panSpeed = 0.8;

	controls.noZoom = false;
	controls.noPan = false;

	controls.staticMoving = true;
	controls.dynamicDampingFactor = 0.3;

	controls.keys = [ 65, 83, 68 ];
    this.controls = controls;
    this.controls.addEventListener('change', this.render);
    
};


/* Return */
console.log("Loaded RMol!");
return RMol;
});