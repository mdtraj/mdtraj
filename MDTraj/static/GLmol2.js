/*
 GLmol - Molecular Viewer on WebGL/Javascript (0.47)
  (C) Copyright 2011-2012, biochem_fan
      License: dual license of MIT or LGPL3

  Contributors:
    Robert Hanson for parseXYZ, deferred instantiation

  This program uses
      Three.js 
         https://github.com/mrdoob/three.js
         Copyright (c) 2010-2012 three.js Authors. All rights reserved.
      jQuery
         http://jquery.org/
         Copyright (c) 2011 John Resig
 */
define(["Three49custom"], function(THREE) {

// Workaround for Intel GMA series (gl_FrontFacing causes compilation error)
THREE.ShaderLib.lambert.fragmentShader = THREE.ShaderLib.lambert.fragmentShader.replace("gl_FrontFacing", "true");
THREE.ShaderLib.lambert.vertexShader = THREE.ShaderLib.lambert.vertexShader.replace(/\}$/, "#ifdef DOUBLE_SIDED\n if (transformedNormal.z < 0.0) vLightFront = vLightBack;\n #endif\n }");

var TV3 = THREE.Vector3, TF3 = THREE.Face3, TCo = THREE.Color;

THREE.Geometry.prototype.colorAll = function (color) {
   for (var i = 0; i < this.faces.length; i++) {
      this.faces[i].color = color;
   }
};

THREE.Matrix4.prototype.isIdentity = function() {
   for (var i = 0; i < 4; i++)
      for (var j = 0; j < 4; j++) 
         if (this.elements[i * 4 + j] != (i == j) ? 1 : 0) return false;
   return true;
};

function GLmol(el, suppressAutoload) {
   if (el) this.create(el, suppressAutoload);
   return true;
}


GLmol.prototype.create = function(el, suppressAutoload) {
   this.Nucleotides = ['  G', '  A', '  T', '  C', '  U', ' DG', ' DA', ' DT', ' DC', ' DU'];
   this.ElementColors = {"H": 0xCCCCCC, "C": 0xAAAAAA, "O": 0xCC0000, "N": 0x0000CC, "S": 0xCCCC00, "P": 0x6622CC,
                         "F": 0x00CC00, "CL": 0x00CC00, "BR": 0x882200, "I": 0x6600AA,
                         "FE": 0xCC6600, "CA": 0x8888AA};
// Reference: A. Bondi, J. Phys. Chem., 1964, 68, 441.
   this.vdwRadii = {"H": 1.2, "Li": 1.82, "Na": 2.27, "K": 2.75, "C": 1.7, "N": 1.55, "O": 1.52,
                   "F": 1.47, "P": 1.80, "S": 1.80, "CL": 1.75, "BR": 1.85, "SE": 1.90,
                   "ZN": 1.39, "CU": 1.4, "NI": 1.63};

   this.aaScale = 1; // or 2

   this.container = el;
   this.WIDTH = this.container.width() * this.aaScale, this.HEIGHT = this.container.height() * this.aaScale;
   this.ASPECT = this.WIDTH / this.HEIGHT;
   this.NEAR = 1, FAR = 800;
   this.CAMERA_Z = -150;
   this.renderer = new THREE.WebGLRenderer({antialias: true});
   this.renderer.sortObjects = false; // hopefully improve performance
   // 'antialias: true' now works in Firefox too!
   // setting this.aaScale = 2 will enable antialias in older Firefox but GPU load increases.
   this.renderer.domElement.style.width = "100%";
   this.renderer.domElement.style.height = "100%";
   this.container.append(this.renderer.domElement);
   this.renderer.setSize(this.WIDTH, this.HEIGHT);

   this.camera = new THREE.PerspectiveCamera(20, this.ASPECT, 1, 800); // will be updated anyway
   this.camera.position = new TV3(0, 0, this.CAMERA_Z);
   this.camera.lookAt(new TV3(0, 0, 0));
   this.perspectiveCamera = this.camera;
   this.orthoscopicCamera = new THREE.OrthographicCamera();
   this.orthoscopicCamera.position.z = this.CAMERA_Z;
   this.orthoscopicCamera.lookAt(new TV3(0, 0, 0));

   var self = this;
   $(window).resize(function() { // only window can capture resize event
      self.WIDTH = self.container.width() * self.aaScale;
      self.HEIGHT = self.container.height() * self.aaScale;
      self.ASPECT = self.WIDTH / self.HEIGHT;
      self.renderer.setSize(self.WIDTH, self.HEIGHT);
      self.camera.aspect = self.ASPECT;
      self.camera.updateProjectionMatrix();
      self.show();
   });

   this.scene = null;
   this.rotationGroup = null; // which contains modelGroup
   this.modelGroup = null;

   this.bgColor = 0x000000;
   this.fov = 20;
   this.fogStart = 0.4;
   this.slabNear = -50; // relative to the center of rotationGroup
   this.slabFar = +50;

   // Default values
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
 
   // UI variables
   this.cq = new THREE.Quaternion(1, 0, 0, 0);
   this.dq = new THREE.Quaternion(1, 0, 0, 0);
   this.isDragging = false;
   this.mouseStartX = 0;
   this.mouseStartY = 0;
   this.currentModelPos = 0;
   this.cz = 0;
   this.enableMouse();

   if (suppressAutoload) return;
   this.loadMolecule();
}


GLmol.prototype.setupLights = function(scene) {
   var directionalLight =  new THREE.DirectionalLight(0xFFFFFF);
   directionalLight.position = new TV3(0.2, 0.2, -1).normalize();
   directionalLight.intensity = 1.2;
   scene.add(directionalLight);
   var ambientLight = new THREE.AmbientLight(0x202020);
   scene.add(ambientLight);
};

GLmol.prototype.parsePDB2 = function(str) {
   var atoms = this.atoms;
   var protein = this.protein;
   var molID;
   
   var atoms_cnt = 0;
   lines = str.split("\n");
   for (var i = 0; i < lines.length; i++) {
       line = lines[i].replace(/^\s*/, ''); // remove indent
       var recordName = line.substr(0, 6);
       if (recordName == 'ATOM  ' || recordName == 'HETATM') {
           var atom, resn, chain, resi, x, y, z, hetflag, elem, serial, altLoc, b;
           altLoc = line.substr(16, 1);
           if (altLoc != ' ' && altLoc != 'A') continue; // FIXME: ad hoc
           serial = parseInt(line.substr(6, 5));
           atom = line.substr(12, 4).replace(/ /g, "");
           resn = line.substr(17, 3);
           chain = line.substr(21, 1);
           resi = parseInt(line.substr(22, 5));
           x = parseFloat(line.substr(30, 8));
           y = parseFloat(line.substr(38, 8));
           z = parseFloat(line.substr(46, 8));
           b = parseFloat(line.substr(60, 8));
           elem = line.substr(76, 2).replace(/ /g, "");
           if (elem == '') { // for some incorrect PDB files
               elem = line.substr(12, 4).replace(/ /g,"");
            }
            if (line[0] == 'H') hetflag = true;
            else hetflag = false;
            atoms[serial] = {'resn': resn, 'x': x, 'y': y, 'z': z, 'elem': elem,
                             'hetflag': hetflag, 'chain': chain, 'resi': resi,
                             'serial': serial, 'atom': atom, 'ss': 'c',
                             'color': 0xFFFFFF, 'bonds': [], 'bondOrder': [],
                             'b': b /*', altLoc': altLoc*/};
       }
   }

   return true;
};


// Catmull-Rom subdivision
GLmol.prototype.subdivide = function(_points, DIV) { // points as Vector3
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

GLmol.prototype.drawAtomsAsSphere = function(group, atomlist, defaultRadius, forceDefault, scale) {
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
   }
};

// about two times faster than sphere when div = 2
GLmol.prototype.drawAtomsAsIcosahedron = function(group, atomlist, defaultRadius, forceDefault) {
   var geo = this.IcosahedronGeometry();
   for (var i = 0; i < atomlist.length; i++) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined) continue;

      var mat = new THREE.MeshLambertMaterial({color: atom.color});
      var sphere = new THREE.Mesh(geo, mat);
      sphere.scale.x = sphere.scale.y = sphere.scale.z = (!forceDefault && this.vdwRadii[atom.elem] != undefined) ? this.vdwRadii[atom.elem] : defaultRadius;
      group.add(sphere);
      sphere.position.x = atom.x;
      sphere.position.y = atom.y;
      sphere.position.z = atom.z;
   }
};

GLmol.prototype.isConnected = function(atom1, atom2) {
   var s = atom1.bonds.indexOf(atom2.serial);
   if (s != -1) return atom1.bondOrder[s];

   if (this.protein.smallMolecule && (atom1.hetflag || atom2.hetflag)) return 0; // CHECK: or should I ?

   var distSquared = (atom1.x - atom2.x) * (atom1.x - atom2.x) + 
                     (atom1.y - atom2.y) * (atom1.y - atom2.y) + 
                     (atom1.z - atom2.z) * (atom1.z - atom2.z);

//   if (atom1.altLoc != atom2.altLoc) return false;
   if (isNaN(distSquared)) return 0;
   if (distSquared < 0.5) return 0; // maybe duplicate position.

   if (distSquared > 1.3 && (atom1.elem == 'H' || atom2.elem == 'H' || atom1.elem == 'D' || atom2.elem == 'D')) return 0;
   if (distSquared < 3.42 && (atom1.elem == 'S' || atom2.elem == 'S')) return 1;
   if (distSquared > 2.78) return 0;
   return 1;
};

GLmol.prototype.drawBondAsStickSub = function(group, atom1, atom2, bondR, order) {
   var delta, tmp;
   if (order > 1) delta = this.calcBondDelta(atom1, atom2, bondR * 2.3);
   var p1 = new TV3(atom1.x, atom1.y, atom1.z);
   var p2 = new TV3(atom2.x, atom2.y, atom2.z);
   var mp = p1.clone().add(p2).multiplyScalar(0.5);

   var c1 = new TCo(atom1.color), c2 = new TCo(atom2.color);
   if (order == 1 || order == 3) {
      this.drawCylinder(group, p1, mp, bondR, atom1.color);
      this.drawCylinder(group, p2, mp, bondR, atom2.color);
   }
   if (order > 1) {
      tmp = mp.clone().add(delta);
      this.drawCylinder(group, p1.clone().add(delta), tmp, bondR, atom1.color);
      this.drawCylinder(group, p2.clone().add(delta), tmp, bondR, atom2.color);
      tmp = mp.clone().sub(delta);
      this.drawCylinder(group, p1.clone().sub(delta), tmp, bondR, atom1.color);
      this.drawCylinder(group, p2.clone().sub(delta), tmp, bondR, atom2.color);
   }
};

GLmol.prototype.drawBondsAsStick = function(group, atomlist, bondR, atomR, ignoreNonbonded, multipleBonds, scale) {
   var sphereGeometry = new THREE.SphereGeometry(1, this.sphereQuality, this.sphereQuality);
   var nAtoms = atomlist.length, mp;
   var forSpheres = [];
   if (!!multipleBonds) bondR /= 2.5;
   for (var _i = 0; _i < nAtoms; _i++) {
      var i = atomlist[_i];
      var atom1 = this.atoms[i];
      if (atom1 == undefined) continue;
      for (var _j = _i + 1; _j < _i + 30 && _j < nAtoms; _j++) {
         var j = atomlist[_j];
         var atom2 = this.atoms[j];
         if (atom2 == undefined) continue;
         var order = this.isConnected(atom1, atom2);
         if (order == 0) continue;
         atom1.connected = atom2.connected = true;
         this.drawBondAsStickSub(group, atom1, atom2, bondR, (!!multipleBonds) ? order : 1);
      }
      for (var _j = 0; _j < atom1.bonds.length; _j++) {
         var j = atom1.bonds[_j];
         if (j < i + 30) continue; // be conservative!
         if (atomlist.indexOf(j) == -1) continue;
         var atom2 = this.atoms[j];
         if (atom2 == undefined) continue;
         atom1.connected = atom2.connected = true;
         this.drawBondAsStickSub(group, atom1, atom2, bondR, (!!multipleBonds) ? atom1.bondOrder[_j] : 1);
      }
      if (atom1.connected) forSpheres.push(i);
   }
   this.drawAtomsAsSphere(group, forSpheres, atomR, !scale, scale);
};

GLmol.prototype.defineCell = function() {
    var p = this.protein;
    if (p.a == undefined) return;

    p.ax = p.a;
    p.ay = 0;
    p.az = 0;
    p.bx = p.b * Math.cos(Math.PI / 180.0 * p.gamma);
    p.by = p.b * Math.sin(Math.PI / 180.0 * p.gamma);
    p.bz = 0;
    p.cx = p.c * Math.cos(Math.PI / 180.0 * p.beta);
    p.cy = p.c * (Math.cos(Math.PI / 180.0 * p.alpha) - 
               Math.cos(Math.PI / 180.0 * p.gamma) 
             * Math.cos(Math.PI / 180.0 * p.beta))
             / Math.sin(Math.PI / 180.0 * p.gamma);
    p.cz = Math.sqrt(p.c * p.c * Math.sin(Math.PI / 180.0 * p.beta)
               * Math.sin(Math.PI / 180.0 * p.beta) - p.cy * p.cy);
};

GLmol.prototype.drawUnitcell = function(group) {
    var p = this.protein;
    if (p.a == undefined) return;

    var vertices = [[0, 0, 0], [p.ax, p.ay, p.az], [p.bx, p.by, p.bz], [p.ax + p.bx, p.ay + p.by, p.az + p.bz],
          [p.cx, p.cy, p.cz], [p.cx + p.ax, p.cy + p.ay,  p.cz + p.az], [p.cx + p.bx, p.cy + p.by, p.cz + p.bz], [p.cx + p.ax + p.bx, p.cy + p.ay + p.by, p.cz + p.az + p.bz]];
    var edges = [0, 1, 0, 2, 1, 3, 2, 3, 4, 5, 4, 6, 5, 7, 6, 7, 0, 4, 1, 5, 2, 6, 3, 7];    

    var geo = new THREE.Geometry();
    for (var i = 0; i < edges.length; i++) {
       geo.vertices.push(new TV3(vertices[edges[i]][0], vertices[edges[i]][1], vertices[edges[i]][2]));
    }
   var lineMaterial = new THREE.LineBasicMaterial({linewidth: 1, color: 0xcccccc});
   var line = new THREE.Line(geo, lineMaterial);
   line.type = THREE.LinePieces;
   group.add(line);
};

// TODO: Find inner side of a ring
GLmol.prototype.calcBondDelta = function(atom1, atom2, sep) {
   var dot;
   var axis = new TV3(atom1.x - atom2.x, atom1.y - atom2.y, atom1.z - atom2.z).normalize();
   var found = null;
   for (var i = 0; i < atom1.bonds.length && !found; i++) {
      var atom = this.atoms[atom1.bonds[i]]; if (!atom) continue;
      if (atom.serial != atom2.serial && atom.elem != 'H') found = atom;
   }
   for (var i = 0; i < atom2.bonds.length && !found; i++) {
      var atom = this.atoms[atom2.bonds[i]]; if (!atom) continue;
      if (atom.serial != atom1.serial && atom.elem != 'H') found = atom;
   }
   if (found) {
      var tmp = new TV3(atom1.x - found.x, atom1.y - found.y, atom1.z - found.z).normalize();
      dot = tmp.dot(axis);
      delta = new TV3(tmp.x - axis.x * dot, tmp.y - axis.y * dot, tmp.z - axis.z * dot);
   }
   if (!found || Math.abs(dot - 1) < 0.001 || Math.abs(dot + 1) < 0.001) {
      if (axis.x < 0.01 && axis.y < 0.01) {
         delta = new TV3(0, -axis.z, axis.y);
      } else {
         delta = new TV3(-axis.y, axis.x, 0);
      }
   }
   delta.normalize().multiplyScalar(sep);
   return delta;
};

GLmol.prototype.drawBondsAsLineSub = function(geo, atom1, atom2, order) {
   var delta, tmp, vs = geo.vertices, cs = geo.colors;
   if (order > 1) delta = this.calcBondDelta(atom1, atom2, 0.15);
   var p1 = new TV3(atom1.x, atom1.y, atom1.z);
   var p2 = new TV3(atom2.x, atom2.y, atom2.z);
   var mp = p1.clone().add(p2).multiplyScalar(0.5);

   var c1 = new TCo(atom1.color), c2 = new TCo(atom2.color);
   if (order == 1 || order == 3) {
      vs.push(p1); cs.push(c1); vs.push(mp); cs.push(c1);
      vs.push(p2); cs.push(c2); vs.push(mp); cs.push(c2);
   }
   if (order > 1) {
      vs.push(p1.clone().add(delta)); cs.push(c1);
      vs.push(tmp = mp.clone().add(delta)); cs.push(c1);
      vs.push(p2.clone().add(delta)); cs.push(c2);
      vs.push(tmp); cs.push(c2);
      vs.push(p1.clone().sub(delta)); cs.push(c1);
      vs.push(tmp = mp.clone().sub(delta)); cs.push(c1);
      vs.push(p2.clone().sub(delta)); cs.push(c2);
      vs.push(tmp); cs.push(c2);
   }
};

GLmol.prototype.drawBondsAsLine = function(group, atomlist, lineWidth) {
   var geo = new THREE.Geometry();
   var nAtoms = atomlist.length;

   for (var _i = 0; _i < nAtoms; _i++) {
      var i = atomlist[_i];
      var  atom1 = this.atoms[i];
      if (atom1 == undefined) continue;
      for (var _j = _i + 1; _j < _i + 30 && _j < nAtoms; _j++) {
         var j = atomlist[_j];
         var atom2 = this.atoms[j];
         if (atom2 == undefined) continue;
         var order = this.isConnected(atom1, atom2);
         if (order == 0) continue;

         this.drawBondsAsLineSub(geo, atom1, atom2, order);
      }
      for (var _j = 0; _j < atom1.bonds.length; _j++) {
          var j = atom1.bonds[_j];
          if (j < i + 30) continue; // be conservative!
          if (atomlist.indexOf(j) == -1) continue;
          var atom2 = this.atoms[j];
          if (atom2 == undefined) continue;
          this.drawBondsAsLineSub(geo, atom1, atom2, atom1.bondOrder[_j]);
      }
    }
   var lineMaterial = new THREE.LineBasicMaterial({linewidth: lineWidth});
   lineMaterial.vertexColors = true;

   var line = new THREE.Line(geo, lineMaterial);
   line.type = THREE.LinePieces;
   group.add(line);
};

GLmol.prototype.drawSmoothCurve = function(group, _points, width, colors, div) {
   if (_points.length == 0) return;

   div = (div == undefined) ? 5 : div;

   var geo = new THREE.Geometry();
   var points = this.subdivide(_points, div);

   for (var i = 0; i < points.length; i++) {
      geo.vertices.push(points[i]);
      geo.colors.push(new TCo(colors[(i == 0) ? 0 : Math.round((i - 1) / div)]));
  }
  var lineMaterial = new THREE.LineBasicMaterial({linewidth: width});
  lineMaterial.vertexColors = true;
  var line = new THREE.Line(geo, lineMaterial);
  line.type = THREE.LineStrip;
  group.add(line);
};

GLmol.prototype.drawAsCross = function(group, atomlist, delta) {
   var geo = new THREE.Geometry();
   var points = [[delta, 0, 0], [-delta, 0, 0], [0, delta, 0], [0, -delta, 0], [0, 0, delta], [0, 0, -delta]];
 
   for (var i = 0, lim = atomlist.length; i < lim; i++) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      var c = new TCo(atom.color);
      for (var j = 0; j < 6; j++) {
         geo.vertices.push(new TV3(atom.x + points[j][0], atom.y + points[j][1], atom.z + points[j][2]));
         geo.colors.push(c);
      }
  }
  var lineMaterial = new THREE.LineBasicMaterial({linewidth: this.lineWidth});
  lineMaterial.vertexColors = true;
  var line = new THREE.Line(geo, lineMaterial, THREE.LinePieces);
  group.add(line);
};

// FIXME: Winkled...
GLmol.prototype.drawSmoothTube = function(group, _points, colors, radii) {
   if (_points.length < 2) return;

   var circleDiv = this.tubeDIV, axisDiv = this.axisDIV;
   var geo = new THREE.Geometry();
   var points = this.subdivide(_points, axisDiv);
   var prevAxis1 = new TV3(), prevAxis2;

   for (var i = 0, lim = points.length; i < lim; i++) {
      var r, idx = (i - 1) / axisDiv;
      if (i == 0) r = radii[0];
      else { 
         if (idx % 1 == 0) r = radii[idx];
         else {
            var floored = Math.floor(idx);
            var tmp = idx - floored;
            r = radii[floored] * tmp + radii[floored + 1] * (1 - tmp);
         }
      }
      var delta, axis1, axis2;

      if (i < lim - 1) {
         delta = new TV3().subVectors(points[i], points[i + 1]);
         axis1 = new TV3(0, - delta.z, delta.y).normalize().multiplyScalar(r);
         axis2 = new TV3().cross(delta, axis1).normalize().multiplyScalar(r);
//      var dir = 1, offset = 0;
         if (prevAxis1.dot(axis1) < 0) {
                 axis1.negate(); axis2.negate();  //dir = -1;//offset = 2 * Math.PI / axisDiv;
         }
         prevAxis1 = axis1; prevAxis2 = axis2;
      } else {
         axis1 = prevAxis1; axis2 = prevAxis2;
      }

      for (var j = 0; j < circleDiv; j++) {
         var angle = 2 * Math.PI / circleDiv * j; //* dir  + offset;
         var c = Math.cos(angle), s = Math.sin(angle);
         geo.vertices.push(new TV3(
         points[i].x + c * axis1.x + s * axis2.x,
         points[i].y + c * axis1.y + s * axis2.y, 
         points[i].z + c * axis1.z + s * axis2.z));
      }
   }

   var offset = 0;
   for (var i = 0, lim = points.length - 1; i < lim; i++) {
      var c =  new TCo(colors[Math.round((i - 1)/ axisDiv)]);

      var reg = 0;
      var r1 = new TV3().subVectors(geo.vertices[offset], geo.vertices[offset + circleDiv]).lengthSq();
      var r2 = new TV3().subVectors(geo.vertices[offset], geo.vertices[offset + circleDiv + 1]).lengthSq();
      if (r1 > r2) {r1 = r2; reg = 1;};
      for (var j = 0; j < circleDiv; j++) {
          geo.faces.push(new TF3(offset + j, offset + (j + reg) % circleDiv + circleDiv, offset + (j + 1) % circleDiv));
          geo.faces.push(new TF3(offset + (j + 1) % circleDiv, offset + (j + reg) % circleDiv + circleDiv, offset + (j + reg + 1) % circleDiv + circleDiv));
          geo.faces[geo.faces.length -2].color = c;
          geo.faces[geo.faces.length -1].color = c;
      }
      offset += circleDiv;
   }
   geo.computeFaceNormals();
   geo.computeVertexNormals(false);
   var mat = new THREE.MeshLambertMaterial();
   mat.vertexColors = THREE.FaceColors;
   var mesh = new THREE.Mesh(geo, mat);
   mesh.doubleSided = true;
   group.add(mesh);
};


GLmol.prototype.drawMainchainCurve = function(group, atomlist, curveWidth, atomName, div) {
   var points = [], colors = [];
   var currentChain, currentResi;
   if (div == undefined) div = 5;

   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined) continue;

      if ((atom.atom == atomName) && !atom.hetflag) {
         if (currentChain != atom.chain || currentResi + 1 != atom.resi) {
            this.drawSmoothCurve(group, points, curveWidth, colors, div);
            points = [];
            colors = [];
         }
         points.push(new TV3(atom.x, atom.y, atom.z));
         colors.push(atom.color);
         currentChain = atom.chain;
         currentResi = atom.resi;
      }
   }
    this.drawSmoothCurve(group, points, curveWidth, colors, div);
};

GLmol.prototype.drawMainchainTube = function(group, atomlist, atomName, radius) {
   var points = [], colors = [], radii = [];
   var currentChain, currentResi;
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined) continue;

      if ((atom.atom == atomName) && !atom.hetflag) {
         if (currentChain != atom.chain || currentResi + 1 != atom.resi) {
            this.drawSmoothTube(group, points, colors, radii);
            points = []; colors = []; radii = [];
         }
         points.push(new TV3(atom.x, atom.y, atom.z));
         if (radius == undefined) {
            radii.push((atom.b > 0) ? atom.b / 100 : 0.3);
         } else {
            radii.push(radius);
         }
         colors.push(atom.color);
         currentChain = atom.chain;
         currentResi = atom.resi;
      }
   }
   this.drawSmoothTube(group, points, colors, radii);
};

GLmol.prototype.drawStrip = function(group, p1, p2, colors, div, thickness) {
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


GLmol.prototype.drawThinStrip = function(group, p1, p2, colors, div) {
   var geo = new THREE.Geometry();
   for (var i = 0, lim = p1.length; i < lim; i++) {
      geo.vertices.push(p1[i]); // 2i
      geo.vertices.push(p2[i]); // 2i + 1
   }
   for (var i = 1, lim = p1.length; i < lim; i++) {
      var f = new THREE.Face4(2 * i, 2 * i + 1, 2 * i - 1, 2 * i - 2);
      f.color = new TCo(colors[Math.round((i - 1)/ div)]);
      geo.faces.push(f);
   }
   geo.computeFaceNormals();
   geo.computeVertexNormals(false);
   var material =  new THREE.MeshLambertMaterial();
   material.vertexColors = THREE.FaceColors;
   var mesh = new THREE.Mesh(geo, material);
   mesh.doubleSided = true;
   group.add(mesh);
};


GLmol.prototype.IcosahedronGeometry = function() {
   if (!this.icosahedron) this.icosahedron = new THREE.IcosahedronGeometry(1);
   return this.icosahedron;
};

GLmol.prototype.drawCylinder = function(group, from, to, radius, color, cap) {
   if (!from || !to) return;

   var midpoint = new TV3().add(from, to).multiplyScalar(0.5);
   var color = new TCo(color);

   if (!this.cylinderGeometry) {
      this.cylinderGeometry = new THREE.CylinderGeometry(1, 1, 1, this.cylinderQuality, 1, !cap);
      this.cylinderGeometry.faceUvs = [];
      this.faceVertexUvs = [];
   }
   var cylinderMaterial = new THREE.MeshLambertMaterial({color: color.getHex()});
   var cylinder = new THREE.Mesh(this.cylinderGeometry, cylinderMaterial);
   cylinder.position = midpoint;
   cylinder.lookAt(from);
   cylinder.updateMatrix();
   cylinder.matrixAutoUpdate = false;
   var m = new THREE.Matrix4().makeScale(radius, radius, from.distanceTo(to));
   m.rotateX(Math.PI / 2);
   cylinder.matrix.multiply(m);
   group.add(cylinder);
};

// FIXME: transition!
GLmol.prototype.drawHelixAsCylinder = function(group, atomlist, radius) {
   var start = null;
   var currentChain, currentResi;

   var others = [], beta = [];

   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined || atom.hetflag) continue;
      if ((atom.ss != 'h' && atom.ss != 's') || atom.ssend || atom.ssbegin) others.push(atom.serial);
      if (atom.ss == 's') beta.push(atom.serial);
      if (atom.atom != 'CA') continue;

      if (atom.ss == 'h' && atom.ssend) {
         if (start != null) this.drawCylinder(group, new TV3(start.x, start.y, start.z), new TV3(atom.x, atom.y, atom.z), radius, atom.color, true);
         start = null;
      }
      currentChain = atom.chain;
      currentResi = atom.resi;
      if (start == null && atom.ss == 'h' && atom.ssbegin) start = atom;
   }
   if (start != null) this.drawCylinder(group, new TV3(start.x, start.y, start.z), new TV3(atom.x, atom.y, atom.z), radius, atom.color);
   this.drawMainchainTube(group, others, "CA", 0.3);
   this.drawStrand(group, beta, undefined, undefined, true,  0, this.helixSheetWidth, false, this.thickness * 2);
};

GLmol.prototype.drawCartoon = function(group, atomlist, doNotSmoothen, thickness) {
   this.drawStrand(group, atomlist, 2, undefined, true, undefined, undefined, doNotSmoothen, thickness);
};

GLmol.prototype.drawStrand = function(group, atomlist, num, div, fill, coilWidth, helixSheetWidth, doNotSmoothen, thickness) {
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
            console.log(currentCA);
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

GLmol.prototype.drawNucleicAcidLadderSub = function(geo, lineGeo, atoms, color) {
//        color.r *= 0.9; color.g *= 0.9; color.b *= 0.9;
   if (atoms[0] != undefined && atoms[1] != undefined && atoms[2] != undefined &&
       atoms[3] != undefined && atoms[4] != undefined && atoms[5] != undefined) {
      var baseFaceId = geo.vertices.length;
      for (var i = 0; i <= 5; i++) geo.vertices.push(atoms[i]);
          geo.faces.push(new TF3(baseFaceId, baseFaceId + 1, baseFaceId + 2));
          geo.faces.push(new TF3(baseFaceId, baseFaceId + 2, baseFaceId + 3));
          geo.faces.push(new TF3(baseFaceId, baseFaceId + 3, baseFaceId + 4));
          geo.faces.push(new TF3(baseFaceId, baseFaceId + 4, baseFaceId + 5));
          for (var j = geo.faces.length - 4, lim = geo.faces.length; j < lim; j++) geo.faces[j].color = color;
    }
    if (atoms[4] != undefined && atoms[3] != undefined && atoms[6] != undefined &&
       atoms[7] != undefined && atoms[8] != undefined) {
       var baseFaceId = geo.vertices.length;
       geo.vertices.push(atoms[4]);
       geo.vertices.push(atoms[3]);
       geo.vertices.push(atoms[6]);
       geo.vertices.push(atoms[7]);
       geo.vertices.push(atoms[8]);
       for (var i = 0; i <= 4; i++) geo.colors.push(color);
       geo.faces.push(new TF3(baseFaceId, baseFaceId + 1, baseFaceId + 2));
       geo.faces.push(new TF3(baseFaceId, baseFaceId + 2, baseFaceId + 3));
       geo.faces.push(new TF3(baseFaceId, baseFaceId + 3, baseFaceId + 4));
       for (var j = geo.faces.length - 3, lim = geo.faces.length; j < lim; j++) geo.faces[j].color = color;
    }
};

GLmol.prototype.drawNucleicAcidLadder = function(group, atomlist) {
   var geo = new THREE.Geometry();
   var lineGeo = new THREE.Geometry();
   var baseAtoms = ["N1", "C2", "N3", "C4", "C5", "C6", "N9", "C8", "N7"];
   var currentChain, currentResi, currentComponent = new Array(baseAtoms.length);
   var color = new TCo(0xcc0000);
   
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined || atom.hetflag) continue;

      if (atom.resi != currentResi || atom.chain != currentChain) {
         this.drawNucleicAcidLadderSub(geo, lineGeo, currentComponent, color);
         currentComponent = new Array(baseAtoms.length);
      }
      var pos = baseAtoms.indexOf(atom.atom);
      if (pos != -1) currentComponent[pos] = new TV3(atom.x, atom.y, atom.z);
      if (atom.atom == 'O3\'') color = new TCo(atom.color);
      currentResi = atom.resi; currentChain = atom.chain;
   }
   this.drawNucleicAcidLadderSub(geo, lineGeo, currentComponent, color);
   geo.computeFaceNormals();
   var mat = new THREE.MeshLambertMaterial();
   mat.vertexColors = THREE.VertexColors;
   var mesh = new THREE.Mesh(geo, mat);
   mesh.doubleSided = true;
   group.add(mesh);
};

GLmol.prototype.drawNucleicAcidStick = function(group, atomlist) {
   var currentChain, currentResi, start = null, end = null;
   
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined || atom.hetflag) continue;

      if (atom.resi != currentResi || atom.chain != currentChain) {
         if (start != null && end != null)
            this.drawCylinder(group, new TV3(start.x, start.y, start.z), 
                              new TV3(end.x, end.y, end.z), 0.3, start.color, true);
         start = null; end = null;
      }
      if (atom.atom == 'O3\'') start = atom;
      if (atom.resn == '  A' || atom.resn == '  G' || atom.resn == ' DA' || atom.resn == ' DG') {
         if (atom.atom == 'N1')  end = atom; //  N1(AG), N3(CTU)
      } else if (atom.atom == 'N3') {
         end = atom;
      }
      currentResi = atom.resi; currentChain = atom.chain;
   }
   if (start != null && end != null)
      this.drawCylinder(group, new TV3(start.x, start.y, start.z), 
                        new TV3(end.x, end.y, end.z), 0.3, start.color, true);
};

GLmol.prototype.drawNucleicAcidLine = function(group, atomlist) {
   var currentChain, currentResi, start = null, end = null;
   var geo = new THREE.Geometry();

   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined || atom.hetflag) continue;

      if (atom.resi != currentResi || atom.chain != currentChain) {
         if (start != null && end != null) {
            geo.vertices.push(new TV3(start.x, start.y, start.z));
            geo.colors.push(new TCo(start.color));
            geo.vertices.push(new TV3(end.x, end.y, end.z));
            geo.colors.push(new TCo(start.color));
         }
         start = null; end = null;
      }
      if (atom.atom == 'O3\'') start = atom;
      if (atom.resn == '  A' || atom.resn == '  G' || atom.resn == ' DA' || atom.resn == ' DG') {
         if (atom.atom == 'N1')  end = atom; //  N1(AG), N3(CTU)
      } else if (atom.atom == 'N3') {
         end = atom;
      }
      currentResi = atom.resi; currentChain = atom.chain;
   }
   if (start != null && end != null) {
      geo.vertices.push(new TV3(start.x, start.y, start.z));
      geo.colors.push(new TCo(start.color));
      geo.vertices.push(new TV3(end.x, end.y, end.z));
      geo.colors.push(new TCo(start.color));
    }
   var mat =  new THREE.LineBasicMaterial({linewidth: 1, linejoin: false});
   mat.linewidth = 1.5; mat.vertexColors = true;
   var line = new THREE.Line(geo, mat, THREE.LinePieces);
   group.add(line);
};

GLmol.prototype.drawCartoonNucleicAcid = function(group, atomlist, div, thickness) {
        this.drawStrandNucleicAcid(group, atomlist, 2, div, true, undefined, thickness);
};

GLmol.prototype.drawStrandNucleicAcid = function(group, atomlist, num, div, fill, nucleicAcidWidth, thickness) {
   nucleicAcidWidth = nucleicAcidWidth || this.nucleicAcidWidth;
   div = div || this.axisDIV;
   num = num || this.nucleicAcidStrandDIV;
   var points = []; for (var k = 0; k < num; k++) points[k] = [];
   var colors = [];
   var currentChain, currentResi, currentO3;
   var prevOO = null;

   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]];
      if (atom == undefined) continue;

      if ((atom.atom == 'O3\'' || atom.atom == 'OP2') && !atom.hetflag) {
         if (atom.atom == 'O3\'') { // to connect 3' end. FIXME: better way to do?
            if (currentChain != atom.chain || currentResi + 1 != atom.resi) {               
               if (currentO3) {
                  for (var j = 0; j < num; j++) {
                     var delta = -1 + 2 / (num - 1) * j;
                     points[j].push(new TV3(currentO3.x + prevOO.x * delta, 
                      currentO3.y + prevOO.y * delta, currentO3.z + prevOO.z * delta));
                  }
               }
               if (fill) this.drawStrip(group, points[0], points[1], colors, div, thickness);
               for (var j = 0; !thickness && j < num; j++)
                  this.drawSmoothCurve(group, points[j], 1 ,colors, div);
               var points = []; for (var k = 0; k < num; k++) points[k] = [];
               colors = [];
               prevOO = null;
            }
            currentO3 = new TV3(atom.x, atom.y, atom.z);
            currentChain = atom.chain;
            currentResi = atom.resi;
            colors.push(atom.color);
         } else { // OP2
            if (!currentO3) {prevOO = null; continue;} // for 5' phosphate (e.g. 3QX3)
            var O = new TV3(atom.x, atom.y, atom.z);
            O.sub(currentO3);
            O.normalize().multiplyScalar(nucleicAcidWidth);  // TODO: refactor
            if (prevOO != undefined && O.dot(prevOO) < 0) {
               O.negate();
            }
            prevOO = O;
            for (var j = 0; j < num; j++) {
               var delta = -1 + 2 / (num - 1) * j;
               points[j].push(new TV3(currentO3.x + prevOO.x * delta, 
                 currentO3.y + prevOO.y * delta, currentO3.z + prevOO.z * delta));
            }
            currentO3 = null;
         }
      }
   }
   if (currentO3) {
      for (var j = 0; j < num; j++) {
         var delta = -1 + 2 / (num - 1) * j;
         points[j].push(new TV3(currentO3.x + prevOO.x * delta, 
           currentO3.y + prevOO.y * delta, currentO3.z + prevOO.z * delta));
      }
   }
   if (fill) this.drawStrip(group, points[0], points[1], colors, div, thickness); 
   for (var j = 0; !thickness && j < num; j++)
      this.drawSmoothCurve(group, points[j], 1 ,colors, div);
};

GLmol.prototype.drawDottedLines = function(group, points, color) {
    var geo = new THREE.Geometry();
    var step = 0.3;

    for (var i = 0, lim = Math.floor(points.length / 2); i < lim; i++) {
        var p1 = points[2 * i], p2 = points[2 * i + 1];
        var delta = p2.clone().sub(p1);
        var dist = delta.length();
        delta.normalize().multiplyScalar(step);
        var jlim =  Math.floor(dist / step);
        for (var j = 0; j < jlim; j++) {
           var p = new TV3(p1.x + delta.x * j, p1.y + delta.y * j, p1.z + delta.z * j);
           geo.vertices.push(p);
        }
        if (jlim % 2 == 1) geo.vertices.push(p2);
    }

    var mat = new THREE.LineBasicMaterial({'color': color.getHex()});
    mat.linewidth = 2;
    var line = new THREE.Line(geo, mat, THREE.LinePieces);
    group.add(line);
};

GLmol.prototype.getAllAtoms = function() {
   var ret = [];
   for (var i in this.atoms) {
      ret.push(this.atoms[i].serial);
   }
   return ret;
};

// Probably I can refactor using higher-order functions.
GLmol.prototype.getHetatms = function(atomlist) {
   var ret = [];
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (atom.hetflag) ret.push(atom.serial);
   }
   return ret;
};

GLmol.prototype.removeSolvents = function(atomlist) {
   var ret = [];
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (atom.resn != 'HOH') ret.push(atom.serial);
   }
   return ret;
};

GLmol.prototype.getProteins = function(atomlist) {
   var ret = [];
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (!atom.hetflag) ret.push(atom.serial);
   }
   return ret;
};

// TODO: Test
GLmol.prototype.excludeAtoms = function(atomlist, deleteList) {
   var ret = [];
   var blackList = new Object();
   for (var _i in deleteList) blackList[deleteList[_i]] = true;

   for (var _i in atomlist) {
      var i = atomlist[_i];

      if (!blackList[i]) ret.push(i);
   }
   return ret;
};

GLmol.prototype.getSidechains = function(atomlist) {
   var ret = [];
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (atom.hetflag) continue;
      if (atom.atom == 'C' || atom.atom == 'O' || (atom.atom == 'N' && atom.resn != "PRO")) continue;
      ret.push(atom.serial);
   }
   return ret;
};

GLmol.prototype.getAtomsWithin = function(atomlist, extent) {
   var ret = [];

   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (atom.x < extent[0][0] || atom.x > extent[1][0]) continue;
      if (atom.y < extent[0][1] || atom.y > extent[1][1]) continue;
      if (atom.z < extent[0][2] || atom.z > extent[1][2]) continue;
      ret.push(atom.serial);      
   }
   return ret;
};

GLmol.prototype.getExtent = function(atomlist) {
   var xmin = ymin = zmin = 9999;
   var xmax = ymax = zmax = -9999;
   var xsum = ysum = zsum = cnt = 0;

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

GLmol.prototype.getResiduesById = function(atomlist, resi) {
   var ret = [];
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (resi.indexOf(atom.resi) != -1) ret.push(atom.serial);
   }
   return ret;
};

GLmol.prototype.getResidueBySS = function(atomlist, ss) {
   var ret = [];
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (ss.indexOf(atom.ss) != -1) ret.push(atom.serial);
   }
   return ret;
};

GLmol.prototype.getChain = function(atomlist, chain) {
   var ret = [], chains = {};
   chain = chain.toString(); // concat if Array
   for (var i = 0, lim = chain.length; i < lim; i++) chains[chain.substr(i, 1)] = true;
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (chains[atom.chain]) ret.push(atom.serial);
   }
   return ret;
};

// for HETATM only
GLmol.prototype.getNonbonded = function(atomlist, chain) {
   var ret = [];
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (atom.hetflag && atom.bonds.length == 0) ret.push(atom.serial);
   }
   return ret;
};

GLmol.prototype.colorByAtom = function(atomlist, colors) {
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      var c = colors[atom.elem];
      if (c == undefined) c = this.ElementColors[atom.elem];
      if (c == undefined) c = this.defaultColor;
      atom.color = c;
   }
};

// MEMO: Color only CA. maybe I should add atom.cartoonColor.
GLmol.prototype.colorByStructure = function(atomlist, helixColor, sheetColor, colorSidechains) {
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (!colorSidechains && (atom.atom != 'CA' || atom.hetflag)) continue;
      if (atom.ss[0] == 's') atom.color = sheetColor;
      else if (atom.ss[0] == 'h') atom.color = helixColor;
   }
};

GLmol.prototype.colorByBFactor = function(atomlist, colorSidechains) {
   var minB = 1000, maxB = -1000;

   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (atom.hetflag) continue;
      if (colorSidechains || atom.atom == 'CA' || atom.atom == 'O3\'') {
         if (minB > atom.b) minB = atom.b;
         if (maxB < atom.b) maxB = atom.b;
      }
   }

   var mid = (maxB + minB) / 2;

   var range = (maxB - minB) / 2;
   if (range < 0.01 && range > -0.01) return;
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (atom.hetflag) continue;
      if (colorSidechains || atom.atom == 'CA' || atom.atom == 'O3\'') {
         var color = new TCo(0);
         if (atom.b < mid)
            color.setHSL(0.667, (mid - atom.b) / range, 1);
         else
            color.setHSL(0, (atom.b - mid) / range, 1);
         atom.color = color.getHex();
      }
   }
};

GLmol.prototype.colorByChain = function(atomlist, colorSidechains) {
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      if (atom.hetflag) continue;
      if (colorSidechains || atom.atom == 'CA' || atom.atom == 'O3\'') {
         var color = new TCo(0);
         color.setHSL((atom.chain.charCodeAt(0) * 5) % 17 / 17.0, 1, 0.9);
         atom.color = color.getHex();
      }
   }
};

GLmol.prototype.colorByResidue = function(atomlist, residueColors) {
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      c = residueColors[atom.resn]
      if (c != undefined) atom.color = c;
   }
};

GLmol.prototype.colorAtoms = function(atomlist, c) {
   for (var i in atomlist) {
      var atom = this.atoms[atomlist[i]]; if (atom == undefined) continue;

      atom.color = c;
   }
};

GLmol.prototype.colorByPolarity = function(atomlist, polar, nonpolar) {
   var polarResidues = ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS'];
   var nonPolarResidues = ['GLY', 'PRO', 'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TYR', 'TRP'];
   var colorMap = {};
   for (var i in polarResidues) colorMap[polarResidues[i]] = polar;
   for (i in nonPolarResidues) colorMap[nonPolarResidues[i]] = nonpolar;
   this.colorByResidue(atomlist, colorMap);   
};

// TODO: Add near(atomlist, neighbor, distanceCutoff)
// TODO: Add expandToResidue(atomlist)

GLmol.prototype.colorChainbow = function(atomlist, colorSidechains) {
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

GLmol.prototype.drawSymmetryMates2 = function(group, asu, matrices) {
   if (matrices == undefined) return;
   asu.matrixAutoUpdate = false;

   var cnt = 1;
   this.protein.appliedMatrix = new THREE.Matrix4();
   for (var i = 0; i < matrices.length; i++) {
      var mat = matrices[i];
      if (mat == undefined || mat.isIdentity()) continue;
      //console.log(mat);
      var symmetryMate = THREE.SceneUtils.cloneObject(asu);
      symmetryMate.matrix = mat;
      group.add(symmetryMate);
      for (var j = 0; j < 16; j++) this.protein.appliedMatrix.elements[j] += mat.elements[j];
      cnt++;
   }
   this.protein.appliedMatrix.multiplyScalar(cnt);
};


GLmol.prototype.drawSymmetryMatesWithTranslation2 = function(group, asu, matrices) {
   if (matrices == undefined) return;
   var p = this.protein;
   asu.matrixAutoUpdate = false;

   for (var i = 0; i < matrices.length; i++) {
      var mat = matrices[i];
      if (mat == undefined) continue;

      for (var a = -1; a <=0; a++) {
         for (var b = -1; b <= 0; b++) {
             for (var c = -1; c <= 0; c++) {
                var translationMat = new THREE.Matrix4().makeTranslation(
                   p.ax * a + p.bx * b + p.cx * c,
                   p.ay * a + p.by * b + p.cy * c,
                   p.az * a + p.bz * b + p.cz * c);
                var symop = mat.clone().multiply(translationMat);
                if (symop.isIdentity()) continue;
                var symmetryMate = THREE.SceneUtils.cloneObject(asu);
                symmetryMate.matrix = symop;
                group.add(symmetryMate);
             }
         }
      }
   }
};

GLmol.prototype.defineRepresentation = function(representationOptions) {
    var all = this.getAllAtoms();
    var hetatm = this.removeSolvents(this.getHetatms(all));
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

    if (mainchainMode == 'ribbon') {
        this.drawCartoon(asu, all, doNotSmoothen);
        this.drawCartoonNucleicAcid(asu, all);
    } else if (mainchainMode == 'thickRibbon') {
        this.drawCartoon(asu, all, doNotSmoothen, this.thickness);
        this.drawCartoonNucleicAcid(asu, all, null, this.thickness);
    } else if (mainchainMode == 'strand') {
        this.drawStrand(asu, all, null, null, null, null, null, doNotSmoothen);
        this.drawStrandNucleicAcid(asu, all);
    } else if (mainchainMode == 'chain') {
        this.drawMainchainCurve(asu, all, this.curveWidth, 'CA', 1);
        this.drawMainchainCurve(asu, all, this.curveWidth, 'O3\'', 1);
    } else if (mainchainMode == 'cylinderHelix') {
        this.drawHelixAsCylinder(asu, all, 1.6);
        this.drawCartoonNucleicAcid(asu, all);
    } else if (mainchainMode == 'tube') {
        this.drawMainchainTube(asu, all, 'CA');
        this.drawMainchainTube(asu, all, 'O3\''); // FIXME: 5' end problem!
    } else if (mainchainMode == 'bonds') {
        this.drawBondsAsLine(asu, all, this.lineWidth);
    }

   if (hetatmMode == 'stick') {
      this.drawBondsAsStick(target, hetatm, this.cylinderRadius, this.cylinderRadius, true);
   } else if (hetatmMode == 'sphere') {
       this.drawAtomsAsSphere(target, hetatm, this.sphereRadius);
   } else if (hetatmMode == 'line') {
       this.drawBondsAsLine(target, hetatm, this.curveWidth);
   } else if (hetatmMode == 'icosahedron') {
       this.drawAtomsAsIcosahedron(target, hetatm, this.sphereRadius);
   } else if (hetatmMode == 'ballAndStick') {
       this.drawBondsAsStick(target, hetatm, this.cylinderRadius / 2.0, this.cylinderRadius, true, false, 0.3);
   } else if (hetatmMode == 'ballAndStick2') {
       this.drawBondsAsStick(target, hetatm, this.cylinderRadius / 2.0, this.cylinderRadius, true, true, 0.3);
   }

   if (sideChains == 'line') {
       this.drawBondsAsLine(this.modelGroup, this.getSidechains(all), this.lineWidth);
   }

   // if (projectionMode == 'perspective') {
   //     this.camera = this.perspectiveCamera;
   // } else if (projectionMode == 'orthoscopic') {
   //     this.camera = this.orthoscopicCamera;
   // }

   this.setBackground(representationOptions.backgroundColor);

};

GLmol.prototype.getView = function() {
   if (!this.modelGroup) return [0, 0, 0, 0, 0, 0, 0, 1];
   var pos = this.modelGroup.position;
   var q = this.rotationGroup.quaternion;
   return [pos.x, pos.y, pos.z, this.rotationGroup.position.z, q.x, q.y, q.z, q.w];
};

GLmol.prototype.setView = function(arg) {
   if (!this.modelGroup || !this.rotationGroup) return;
   this.modelGroup.position.x = arg[0];
   this.modelGroup.position.y = arg[1];
   this.modelGroup.position.z = arg[2];
   this.rotationGroup.position.z = arg[3];
   this.rotationGroup.quaternion.x = arg[4];
   this.rotationGroup.quaternion.y = arg[5];
   this.rotationGroup.quaternion.z = arg[6];
   this.rotationGroup.quaternion.w = arg[7];
   this.show();
};

GLmol.prototype.setBackground = function(hex, a) {
   a = a | 1.0;
   this.bgColor = hex;
   this.renderer.setClearColorHex(hex, a);
   this.scene.fog.color = new TCo(hex);
};

GLmol.prototype.initializeScene = function() {
   // CHECK: Should I explicitly call scene.deallocateObject?
   this.scene = new THREE.Scene();
   this.scene.fog = new THREE.Fog(this.bgColor, 100, 200);

   this.modelGroup = new THREE.Object3D();
   this.rotationGroup = new THREE.Object3D();
   this.rotationGroup.useQuaternion = true;
   this.rotationGroup.quaternion = new THREE.Quaternion(1, 0, 0, 0);
   this.rotationGroup.add(this.modelGroup);

   this.scene.add(this.rotationGroup);
   this.setupLights(this.scene);
};

GLmol.prototype.zoomInto = function(atomlist, keepSlab) {
   var tmp = this.getExtent(atomlist);
   var center = new TV3(tmp[2][0], tmp[2][1], tmp[2][2]);//(tmp[0][0] + tmp[1][0]) / 2, (tmp[0][1] + tmp[1][1]) / 2, (tmp[0][2] + tmp[1][2]) / 2);
   if (this.protein.appliedMatrix) {center = this.protein.appliedMatrix.multiplyVector3(center);}
   this.modelGroup.position = center.multiplyScalar(-1);
   var x = tmp[1][0] - tmp[0][0], y = tmp[1][1] - tmp[0][1], z = tmp[1][2] - tmp[0][2];

   var maxD = Math.sqrt(x * x + y * y + z * z);
   if (maxD < 25) maxD = 25;

   if (!keepSlab) {
      this.slabNear = -maxD / 1.9;
      this.slabFar = maxD / 3;
   }

   this.rotationGroup.position.z = maxD * 0.35 / Math.tan(Math.PI / 180.0 * this.camera.fov / 2) - 150;
   this.rotationGroup.quaternion = new THREE.Quaternion(1, 0, 0, 0);
};

GLmol.prototype.rebuildScene = function(representationOptions) {
   time = new Date();

   var view = this.getView();
   this.initializeScene();
   this.defineRepresentation(representationOptions);
   this.setView(view);

   //console.log("builded scene in " + (+new Date() - time) + "ms");
};

GLmol.prototype.loadMoleculeStr = function(repressZoom, source) {
   var time = new Date();

   this.protein = {sheet: [], helix: [], biomtChains: '', biomtMatrices: [], symMat: [], pdbID: '', title: ''};
   this.atoms = [];

   this.parsePDB2(source);
   console.log("parsed in " + (+new Date() - time) + "ms");
   
   this.rebuildScene(true);
   if (repressZoom == undefined || !repressZoom) this.zoomInto(this.getAllAtoms());

   this.show();
 };

GLmol.prototype.setSlabAndFog = function() {
   var center = this.rotationGroup.position.z - this.camera.position.z;
   if (center < 1) center = 1;
   this.camera.near = center + this.slabNear;
   if (this.camera.near < 1) this.camera.near = 1;
   this.camera.far = center + this.slabFar;
   if (this.camera.near + 1 > this.camera.far) this.camera.far = this.camera.near + 1;
   if (this.camera instanceof THREE.PerspectiveCamera) {
      this.camera.fov = this.fov;
   } else {
      this.camera.right = center * Math.tan(Math.PI / 180 * this.fov);
      this.camera.left = - this.camera.right;
      this.camera.top = this.camera.right / this.ASPECT;
      this.camera.bottom = - this.camera.top;
   }
   this.camera.updateProjectionMatrix();
   this.scene.fog.near = this.camera.near + this.fogStart * (this.camera.far - this.camera.near);
//   if (this.scene.fog.near > center) this.scene.fog.near = center;
   this.scene.fog.far = this.camera.far;
};

GLmol.prototype.enableMouse = function() {
   var me = this, glDOM = $(this.renderer.domElement); 

   // TODO: Better touch panel support. 
   // Contribution is needed as I don't own any iOS or Android device with WebGL support.
   glDOM.bind('mousedown touchstart', function(ev) {
      ev.preventDefault();
      if (!me.scene) return;
      var x = ev.pageX, y = ev.pageY;
      if (ev.originalEvent.targetTouches && ev.originalEvent.targetTouches[0]) {
         x = ev.originalEvent.targetTouches[0].pageX;
         y = ev.originalEvent.targetTouches[0].pageY;
      }
      if (x == undefined) return;
      me.isDragging = true;
      me.mouseButton = ev.which;
      me.mouseStartX = x;
      me.mouseStartY = y;
      me.cq = me.rotationGroup.quaternion;
      me.cz = me.rotationGroup.position.z;
      me.currentModelPos = me.modelGroup.position.clone();
      me.cslabNear = me.slabNear;
      me.cslabFar = me.slabFar;
    });

   glDOM.bind('DOMMouseScroll mousewheel', function(ev) { // Zoom
      ev.preventDefault();
      if (!me.scene) return;
      var scaleFactor = (me.rotationGroup.position.z - me.CAMERA_Z) * 0.85;
      if (ev.originalEvent.detail) { // Webkit
         me.rotationGroup.position.z += scaleFactor * ev.originalEvent.detail / 10;
      } else if (ev.originalEvent.wheelDelta) { // Firefox
         me.rotationGroup.position.z -= scaleFactor * ev.originalEvent.wheelDelta / 400;
      }
      //console.log(ev.originalEvent.wheelDelta, ev.originalEvent.detail, me.rotationGroup.position.z);
      me.show();
   });
   glDOM.bind("contextmenu", function(ev) {ev.preventDefault();});
   $('body').bind('mouseup touchend', function(ev) {
      me.isDragging = false;
   });

   glDOM.bind('mousemove touchmove', function(ev) { // touchmove
      ev.preventDefault();
      if (!me.scene) return;
      if (!me.isDragging) return;
      var mode = 0;
      var modeRadio = $('input[name=' + me.id + '_mouseMode]:checked');
      if (modeRadio.length > 0) mode = parseInt(modeRadio.val());

      var x = ev.pageX, y = ev.pageY;
      if (ev.originalEvent.targetTouches && ev.originalEvent.targetTouches[0]) {
         x = ev.originalEvent.targetTouches[0].pageX;
         y = ev.originalEvent.targetTouches[0].pageY;
      }
      if (x == undefined) return;
      var dx = (x - me.mouseStartX) / me.WIDTH;
      var dy = (y - me.mouseStartY) / me.HEIGHT;
      var r = Math.sqrt(dx * dx + dy * dy);
      if (mode == 3 || (me.mouseButton == 3 && ev.ctrlKey)) { // Slab
          me.slabNear = me.cslabNear + dx * 100;
          me.slabFar = me.cslabFar + dy * 100;
      } else if (mode == 2 || me.mouseButton == 3 || ev.shiftKey) { // Zoom
         var scaleFactor = (me.rotationGroup.position.z - me.CAMERA_Z) * 0.85; 
         if (scaleFactor < 80) scaleFactor = 80;
         me.rotationGroup.position.z = me.cz - dy * scaleFactor;
      } else if (mode == 1 || me.mouseButton == 2 || ev.ctrlKey) { // Translate
         var scaleFactor = (me.rotationGroup.position.z - me.CAMERA_Z) * 0.85;
         if (scaleFactor < 20) scaleFactor = 20;
         var translationByScreen = new TV3(- dx * scaleFactor, - dy * scaleFactor, 0);
         var q = me.rotationGroup.quaternion;
         var qinv = new THREE.Quaternion(q.x, q.y, q.z, q.w).inverse().normalize(); 
         var translation = qinv.multiplyVector3(translationByScreen);
         me.modelGroup.position.x = me.currentModelPos.x + translation.x;
         me.modelGroup.position.y = me.currentModelPos.y + translation.y;
         me.modelGroup.position.z = me.currentModelPos.z + translation.z;
      } else if ((mode == 0 || me.mouseButton == 1) && r != 0) { // Rotate
         var rs = Math.sin(r * Math.PI) / r;
         me.dq.x = Math.cos(r * Math.PI); 
         me.dq.y = 0;
         me.dq.z =  rs * dx; 
         me.dq.w =  rs * dy;
         me.rotationGroup.quaternion = new THREE.Quaternion(1, 0, 0, 0); 
         me.rotationGroup.quaternion.multiply(me.dq);
         me.rotationGroup.quaternion.multiply(me.cq);
      }
      me.show();
   });
};


GLmol.prototype.show = function() {
   if (!this.scene) return;

   var time = new Date();
   this.setSlabAndFog();
   this.renderer.render(this.scene, this.camera);
   //console.log("rendered in " + (+new Date() - time) + "ms");
};

return GLmol;
});