/*!
 * iview is an interactive WebGL visualizer for protein-ligand complex. iview is based on GLmol, three.js and jQuery.
 * http://github.com/HongjianLi/istar
 * Copyright (c) 2012-2014 Chinese University of Hong Kong
 * License: Apache License 2.0
 * Hongjian Li, Kwong-Sak Leung, Takanori Nakane and Man-Hon Wong.
 * iview: an interactive WebGL visualizer for protein-ligand complex.
 * BMC Bioinformatics, 15(1):56, 2014.
 *
 * GLmol
 * https://github.com/biochem-fan/GLmol
 * Copyright (c) 2011-2012 biochem_fan
 * License: dual license of MIT or LGPL3
 *
 * three.js
 * https://github.com/mrdoob/three.js
 * Copyright (c) 2010-2012 three.js Authors. All rights reserved.
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * jQuery
 * http://jquery.org
 * Copyright (c) 2011 John Resig
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/*
Adapted from
https://github.com/HongjianLi/istar/blob/master/public/iview/iview.js
version : 27b5c4905f9d9269ace950933e05c050ca3df9cb
We've made a couple small modifications
 - Added methods for setting the topology / coordinates without a PDB file.
   Instead, we want a closer connection to the parsed formats inside the mdtraj
   Trajectory and Topology objects.
 - Also, I changed the constructor to take a jQuery element instead of the id
   tag.
*/

if (typeof jQuery === 'undefined') { throw new Error('iview requires jQuery') }

var iview = function ($el) {
	this.container = $el;
	this.container.widthInv  = 1 / this.container.width();
	this.container.heightInv = 1 / this.container.height();
	this.container.whratio = this.container.width() / this.container.height();
	this.container.hwratio = this.container.height() / this.container.width();
	this.renderer = new THREE.WebGLRenderer({
		canvas: this.container.get(0),
		antialias: true,
                preserveDrawingBuffer: true,
	});
	this.effects = {
        // 'anaglyph': new THREE.AnaglyphEffect(this.renderer),
        // 'parallax barrier': new THREE.ParallaxBarrierEffect(this.renderer),
        // 'oculus rift': new THREE.OculusRiftEffect(this.renderer),
        // 'stereo': new THREE.StereoEffect(this.renderer),
		'none': this.renderer,
	};

	this.camera_z = -150;
	this.perspectiveCamera = new THREE.PerspectiveCamera(20, this.container.whratio, 1, 800);
	this.perspectiveCamera.position.set(0, 0, this.camera_z);
	this.perspectiveCamera.lookAt(new THREE.Vector3(0, 0, 0));
	this.orthographicCamera = new THREE.OrthographicCamera();
	this.orthographicCamera.position.set(0, 0, this.camera_z);
	this.orthographicCamera.lookAt(new THREE.Vector3(0, 0, 0));
	this.cameras = {
		 perspective: this.perspectiveCamera,
		orthographic: this.orthographicCamera,
	};
	this.camera = this.perspectiveCamera;

	this.slabNear = -50; // relative to the center of rot
	this.slabFar  = +50;

	// Default values
	this.sphereGeometry = new THREE.SphereGeometry(1, 32, 32);
	this.cylinderGeometry = new THREE.CylinderGeometry(1, 1, 1, 32, 1);
	this.sphereRadius = 1.5;
	this.cylinderRadius = 0.4;
	this.linewidth = 2;
	this.curveWidth = 3;
	this.helixSheetWidth = 1.3;
	this.coilWidth = 0.3;
	this.thickness = 0.4;
	this.axisDIV = 5; // 3
	this.strandDIV = 6;
	this.tubeDIV = 8;
	this.options = {
		camera: 'perspective',
		background: 'black',
		colorBy: 'spectrum',
		primaryStructure: 'nothing',
		secondaryStructure: 'cylinder & plate',
		surface: 'nothing',
		wireframe: 'no',
		opacity: '0.8',
		ligands: 'stick',
		waters: 'dot',
		ions: 'dot',
		labels: 'no',
		effect: 'none',
	};
	this.labelGeo = new THREE.Geometry();
	for (var i = 0; i < 6; ++i) {
		this.labelGeo.vertices.push(new THREE.Vector3(0, 0, 0));
	}
	this.labelGeo.faces.push(new THREE.Face3(0, 1, 2));
	this.labelGeo.faces.push(new THREE.Face3(0, 2, 3));
	this.labelGeo.faceVertexUvs[0].push([new THREE.Vector2(0, 0), new THREE.Vector2(1, 1), new THREE.Vector2(0, 1)]);
	this.labelGeo.faceVertexUvs[0].push([new THREE.Vector2(0, 0), new THREE.Vector2(1, 0), new THREE.Vector2(1, 1)]);
	this.labelByResn = function (atom) {
		return atom.chain + ':' + atom.resn + ':' + atom.resi;
	};
	this.labelByName = function (atom) {
		return atom.name;
	};

	var me = this;
	this.container.bind('contextmenu', function (e) {
		e.preventDefault();
	});
	this.container.bind('mouseup touchend', function (e) {
		me.isDragging = false;
	});
	this.container.bind('mousedown touchstart', function (e) {
		e.preventDefault();
		if (!me.scene) return;
		var x = e.pageX, y = e.pageY;
		if (e.originalEvent.targetTouches && e.originalEvent.targetTouches[0]) {
			x = e.originalEvent.targetTouches[0].pageX;
			y = e.originalEvent.targetTouches[0].pageY;
		}
		me.isDragging = true;
		me.mouseButton = e.which;
		me.mouseStartX = x;
		me.mouseStartY = y;
		me.cq = me.rot.quaternion.clone();
		me.cz = me.rot.position.z;
		me.cp = me.mdl.position.clone();
		me.cslabNear = me.slabNear;
		me.cslabFar = me.slabFar;
	});
	this.container.bind('mousemove touchmove', function (e) {
		e.preventDefault();
		if (!me.scene) return;
		if (!me.isDragging) return;
		var x = e.pageX, y = e.pageY;
		if (e.originalEvent.targetTouches && e.originalEvent.targetTouches[0]) {
			x = e.originalEvent.targetTouches[0].pageX;
			y = e.originalEvent.targetTouches[0].pageY;
		}
		var dx = (x - me.mouseStartX) * me.container.widthInv;
		var dy = (y - me.mouseStartY) * me.container.heightInv;
		if (!dx && !dy) return;
		var mode = $('#mode .active').text().trim();
		if (mode === 'slab' || e.ctrlKey && e.shiftKey) { // Slab
			me.slabNear = me.cslabNear + dx * 100;
			me.slabFar  = me.cslabFar  + dy * 100;
		} else if (mode === 'translate' || e.ctrlKey || me.mouseButton == 3) { // Translate
			var scaleFactor = (me.rot.position.z - me.camera_z) * 0.85;
			if (scaleFactor < 20) scaleFactor = 20;
			me.mdl.position.copy(me.cp).add(new THREE.Vector3(-dx * scaleFactor, -dy * scaleFactor, 0).applyQuaternion(me.rot.quaternion.clone().inverse().normalize()));
		} else if (mode === 'zoom' || e.shiftKey || me.mouseButton == 2) { // Zoom
			var scaleFactor = (me.rot.position.z - me.camera_z) * 0.85;
			if (scaleFactor < 80) scaleFactor = 80;
			me.rot.position.z = me.cz - dy * scaleFactor;
		} else { // Rotate
			var r = Math.sqrt(dx * dx + dy * dy);
			var rs = Math.sin(r * Math.PI) / r;
			me.rot.quaternion.set(1, 0, 0, 0).multiply(new THREE.Quaternion(Math.cos(r * Math.PI), 0, rs * dx, rs * dy)).multiply(me.cq);
		}
		me.render();
	});
	this.container.bind('mousewheel', function (e) {
		e.preventDefault();
		if (!me.scene) return;
		me.rot.position.z -= e.originalEvent.wheelDelta * 0.025;
		me.render();
	});
	this.container.bind('DOMMouseScroll', function (e) {
		e.preventDefault();
		if (!me.scene) return;
		me.rot.position.z += e.originalEvent.detail;
		me.render();
	});
};

iview.prototype = {

	constructor: iview,

	vdwRadii: { // Hu, S.Z.; Zhou, Z.H.; Tsai, K.R. Acta Phys.-Chim. Sin., 2003, 19:1073.
		 H: 1.08,
		HE: 1.34,
		LI: 1.75,
		BE: 2.05,
		 B: 1.47,
		 C: 1.49,
		 N: 1.41,
		 O: 1.40,
		 F: 1.39,
		NE: 1.68,
		NA: 1.84,
		MG: 2.05,
		AL: 2.11,
		SI: 2.07,
		 P: 1.92,
		 S: 1.82,
		CL: 1.83,
		AR: 1.93,
		 K: 2.05,
		CA: 2.21,
		SC: 2.16,
		TI: 1.87,
		 V: 1.79,
		CR: 1.89,
		MN: 1.97,
		FE: 1.94,
		CO: 1.92,
		NI: 1.84,
		CU: 1.86,
		ZN: 2.10,
		GA: 2.08,
		GE: 2.15,
		AS: 2.06,
		SE: 1.93,
		BR: 1.98,
		KR: 2.12,
		RB: 2.16,
		SR: 2.24,
		 Y: 2.19,
		ZR: 1.86,
		NB: 2.07,
		MO: 2.09,
		TC: 2.09,
		RU: 2.07,
		RH: 1.95,
		PD: 2.02,
		AG: 2.03,
		CD: 2.30,
		IN: 2.36,
		SN: 2.33,
		SB: 2.25,
		TE: 2.23,
		 I: 2.23,
		XE: 2.21,
		CS: 2.22,
		BA: 2.51,
		LA: 2.40,
		CE: 2.35,
		PR: 2.39,
		ND: 2.29,
		PM: 2.36,
		SM: 2.29,
		EU: 2.33,
		GD: 2.37,
		TB: 2.21,
		DY: 2.29,
		HO: 2.16,
		ER: 2.35,
		TM: 2.27,
		YB: 2.42,
		LU: 2.21,
		HF: 2.12,
		TA: 2.17,
		 W: 2.10,
		RE: 2.17,
		OS: 2.16,
		IR: 2.02,
		PT: 2.09,
		AU: 2.17,
		HG: 2.09,
		TL: 2.35,
		PB: 2.32,
		BI: 2.43,
		PO: 2.29,
		AT: 2.36,
		RN: 2.43,
		FR: 2.56,
		RA: 2.43,
		AC: 2.60,
		TH: 2.37,
		PA: 2.43,
		 U: 2.40,
		NP: 2.21,
		PU: 2.56,
		AM: 2.56,
		CM: 2.56,
		BK: 2.56,
		CF: 2.56,
		ES: 2.56,
		FM: 2.56,
	},

	covalentRadii: { // http://en.wikipedia.org/wiki/Covalent_radius
		 H: 0.31,
		HE: 0.28,
		LI: 1.28,
		BE: 0.96,
		 B: 0.84,
		 C: 0.76,
		 N: 0.71,
		 O: 0.66,
		 F: 0.57,
		NE: 0.58,
		NA: 1.66,
		MG: 1.41,
		AL: 1.21,
		SI: 1.11,
		 P: 1.07,
		 S: 1.05,
		CL: 1.02,
		AR: 1.06,
		 K: 2.03,
		CA: 1.76,
		SC: 1.70,
		TI: 1.60,
		 V: 1.53,
		CR: 1.39,
		MN: 1.39,
		FE: 1.32,
		CO: 1.26,
		NI: 1.24,
		CU: 1.32,
		ZN: 1.22,
		GA: 1.22,
		GE: 1.20,
		AS: 1.19,
		SE: 1.20,
		BR: 1.20,
		KR: 1.16,
		RB: 2.20,
		SR: 1.95,
		 Y: 1.90,
		ZR: 1.75,
		NB: 1.64,
		MO: 1.54,
		TC: 1.47,
		RU: 1.46,
		RH: 1.42,
		PD: 1.39,
		AG: 1.45,
		CD: 1.44,
		IN: 1.42,
		SN: 1.39,
		SB: 1.39,
		TE: 1.38,
		 I: 1.39,
		XE: 1.40,
		CS: 2.44,
		BA: 2.15,
		LA: 2.07,
		CE: 2.04,
		PR: 2.03,
		ND: 2.01,
		PM: 1.99,
		SM: 1.98,
		EU: 1.98,
		GD: 1.96,
		TB: 1.94,
		DY: 1.92,
		HO: 1.92,
		ER: 1.89,
		TM: 1.90,
		YB: 1.87,
		LU: 1.87,
		HF: 1.75,
		TA: 1.70,
		 W: 1.62,
		RE: 1.51,
		OS: 1.44,
		IR: 1.41,
		PT: 1.36,
		AU: 1.36,
		HG: 1.32,
		TL: 1.45,
		PB: 1.46,
		BI: 1.48,
		PO: 1.40,
		AT: 1.50,
		RN: 1.50,
		FR: 2.60,
		RA: 2.21,
		AC: 2.15,
		TH: 2.06,
		PA: 2.00,
		 U: 1.96,
		NP: 1.90,
		PU: 1.87,
		AM: 1.80,
		CM: 1.69,
	},

	atomColors: {
		 H: new THREE.Color(0xFFFFFF),
		HE: new THREE.Color(0xD9FFFF),
		LI: new THREE.Color(0xCC80FF),
		BE: new THREE.Color(0xC2FF00),
		 B: new THREE.Color(0xFFB5B5),
		 C: new THREE.Color(0x909090),
		 N: new THREE.Color(0x3050F8),
		 O: new THREE.Color(0xFF0D0D),
		 F: new THREE.Color(0x90E050),
		NE: new THREE.Color(0xB3E3F5),
		NA: new THREE.Color(0xAB5CF2),
		MG: new THREE.Color(0x8AFF00),
		AL: new THREE.Color(0xBFA6A6),
		SI: new THREE.Color(0xF0C8A0),
		 P: new THREE.Color(0xFF8000),
		 S: new THREE.Color(0xFFFF30),
		CL: new THREE.Color(0x1FF01F),
		AR: new THREE.Color(0x80D1E3),
		 K: new THREE.Color(0x8F40D4),
		CA: new THREE.Color(0x3DFF00),
		SC: new THREE.Color(0xE6E6E6),
		TI: new THREE.Color(0xBFC2C7),
		 V: new THREE.Color(0xA6A6AB),
		CR: new THREE.Color(0x8A99C7),
		MN: new THREE.Color(0x9C7AC7),
		FE: new THREE.Color(0xE06633),
		CO: new THREE.Color(0xF090A0),
		NI: new THREE.Color(0x50D050),
		CU: new THREE.Color(0xC88033),
		ZN: new THREE.Color(0x7D80B0),
		GA: new THREE.Color(0xC28F8F),
		GE: new THREE.Color(0x668F8F),
		AS: new THREE.Color(0xBD80E3),
		SE: new THREE.Color(0xFFA100),
		BR: new THREE.Color(0xA62929),
		KR: new THREE.Color(0x5CB8D1),
		RB: new THREE.Color(0x702EB0),
		SR: new THREE.Color(0x00FF00),
		 Y: new THREE.Color(0x94FFFF),
		ZR: new THREE.Color(0x94E0E0),
		NB: new THREE.Color(0x73C2C9),
		MO: new THREE.Color(0x54B5B5),
		TC: new THREE.Color(0x3B9E9E),
		RU: new THREE.Color(0x248F8F),
		RH: new THREE.Color(0x0A7D8C),
		PD: new THREE.Color(0x006985),
		AG: new THREE.Color(0xC0C0C0),
		CD: new THREE.Color(0xFFD98F),
		IN: new THREE.Color(0xA67573),
		SN: new THREE.Color(0x668080),
		SB: new THREE.Color(0x9E63B5),
		TE: new THREE.Color(0xD47A00),
		 I: new THREE.Color(0x940094),
		XE: new THREE.Color(0x429EB0),
		CS: new THREE.Color(0x57178F),
		BA: new THREE.Color(0x00C900),
		LA: new THREE.Color(0x70D4FF),
		CE: new THREE.Color(0xFFFFC7),
		PR: new THREE.Color(0xD9FFC7),
		ND: new THREE.Color(0xC7FFC7),
		PM: new THREE.Color(0xA3FFC7),
		SM: new THREE.Color(0x8FFFC7),
		EU: new THREE.Color(0x61FFC7),
		GD: new THREE.Color(0x45FFC7),
		TB: new THREE.Color(0x30FFC7),
		DY: new THREE.Color(0x1FFFC7),
		HO: new THREE.Color(0x00FF9C),
		ER: new THREE.Color(0x00E675),
		TM: new THREE.Color(0x00D452),
		YB: new THREE.Color(0x00BF38),
		LU: new THREE.Color(0x00AB24),
		HF: new THREE.Color(0x4DC2FF),
		TA: new THREE.Color(0x4DA6FF),
		 W: new THREE.Color(0x2194D6),
		RE: new THREE.Color(0x267DAB),
		OS: new THREE.Color(0x266696),
		IR: new THREE.Color(0x175487),
		PT: new THREE.Color(0xD0D0E0),
		AU: new THREE.Color(0xFFD123),
		HG: new THREE.Color(0xB8B8D0),
		TL: new THREE.Color(0xA6544D),
		PB: new THREE.Color(0x575961),
		BI: new THREE.Color(0x9E4FB5),
		PO: new THREE.Color(0xAB5C00),
		AT: new THREE.Color(0x754F45),
		RN: new THREE.Color(0x428296),
		FR: new THREE.Color(0x420066),
		RA: new THREE.Color(0x007D00),
		AC: new THREE.Color(0x70ABFA),
		TH: new THREE.Color(0x00BAFF),
		PA: new THREE.Color(0x00A1FF),
		 U: new THREE.Color(0x008FFF),
		NP: new THREE.Color(0x0080FF),
		PU: new THREE.Color(0x006BFF),
		AM: new THREE.Color(0x545CF2),
		CM: new THREE.Color(0x785CE3),
		BK: new THREE.Color(0x8A4FE3),
		CF: new THREE.Color(0xA136D4),
		ES: new THREE.Color(0xB31FD4),
		FM: new THREE.Color(0xB31FBA),
	},

	defaultAtomColor: new THREE.Color(0xCCCCCC),

	stdChainColors: {
		A: new THREE.Color(0xC0D0FF),
		B: new THREE.Color(0xB0FFB0),
		C: new THREE.Color(0xFFC0C8),
		D: new THREE.Color(0xFFFF80),
		E: new THREE.Color(0xFFC0FF),
		F: new THREE.Color(0xB0F0F0),
		G: new THREE.Color(0xFFD070),
		H: new THREE.Color(0xF08080),
		I: new THREE.Color(0xF5DEB3),
		J: new THREE.Color(0x00BFFF),
		K: new THREE.Color(0xCD5C5C),
		L: new THREE.Color(0x66CDAA),
		M: new THREE.Color(0x9ACD32),
		N: new THREE.Color(0xEE82EE),
		O: new THREE.Color(0x00CED1),
		P: new THREE.Color(0x00FF7F),
		Q: new THREE.Color(0x3CB371),
		R: new THREE.Color(0x00008B),
		S: new THREE.Color(0xBDB76B),
		T: new THREE.Color(0x006400),
		U: new THREE.Color(0x800000),
		V: new THREE.Color(0x808000),
		W: new THREE.Color(0x800080),
		X: new THREE.Color(0x008080),
		Y: new THREE.Color(0xB8860B),
		Z: new THREE.Color(0xB22222),
	},

	hetChainColors: {
		A: new THREE.Color(0x90A0CF),
		B: new THREE.Color(0x80CF98),
		C: new THREE.Color(0xCF90B0),
		D: new THREE.Color(0xCFCF70),
		E: new THREE.Color(0xCF90CF),
		F: new THREE.Color(0x80C0C0),
		G: new THREE.Color(0xCFA060),
		H: new THREE.Color(0xC05070),
		I: new THREE.Color(0xC5AE83),
		J: new THREE.Color(0x00A7CF),
		K: new THREE.Color(0xB54C4C),
		L: new THREE.Color(0x56B592),
		M: new THREE.Color(0x8AB52A),
		N: new THREE.Color(0xBE72BE),
		O: new THREE.Color(0x00B6A1),
		P: new THREE.Color(0x00CF6F),
		Q: new THREE.Color(0x349B61),
		R: new THREE.Color(0x0000BB),
		S: new THREE.Color(0xA59F5B),
		T: new THREE.Color(0x009400),
		U: new THREE.Color(0xB00000),
		V: new THREE.Color(0xB0B000),
		W: new THREE.Color(0xB000B0),
		X: new THREE.Color(0x00B0B0),
		Y: new THREE.Color(0xE8B613),
		Z: new THREE.Color(0xC23232),
	},

	backgroundColors: {
		black: new THREE.Color(0x000000),
		 grey: new THREE.Color(0xCCCCCC),
		white: new THREE.Color(0xFFFFFF),
	},

	residueColors: {
		ALA: new THREE.Color(0xC8C8C8),
		ARG: new THREE.Color(0x145AFF),
		ASN: new THREE.Color(0x00DCDC),
		ASP: new THREE.Color(0xE60A0A),
		CYS: new THREE.Color(0xE6E600),
		GLN: new THREE.Color(0x00DCDC),
		GLU: new THREE.Color(0xE60A0A),
		GLY: new THREE.Color(0xEBEBEB),
		HIS: new THREE.Color(0x8282D2),
		ILE: new THREE.Color(0x0F820F),
		LEU: new THREE.Color(0x0F820F),
		LYS: new THREE.Color(0x145AFF),
		MET: new THREE.Color(0xE6E600),
		PHE: new THREE.Color(0x3232AA),
		PRO: new THREE.Color(0xDC9682),
		SER: new THREE.Color(0xFA9600),
		THR: new THREE.Color(0xFA9600),
		TRP: new THREE.Color(0xB45AB4),
		TYR: new THREE.Color(0x3232AA),
		VAL: new THREE.Color(0x0F820F),
		ASX: new THREE.Color(0xFF69B4),
		GLX: new THREE.Color(0xFF69B4),
	},

	defaultResidueColor: new THREE.Color(0xBEA06E),

	polarColor: new THREE.Color(0xCC0000),

	nonpolarColor: new THREE.Color(0x00CCCC),

	polarityColors: {
		ARG: this.polarColor,
		HIS: this.polarColor,
		LYS: this.polarColor,
		ASP: this.polarColor,
		GLU: this.polarColor,
		SER: this.polarColor,
		THR: this.polarColor,
		ASN: this.polarColor,
		GLN: this.polarColor,
		TYR: this.polarColor,
		GLY: this.nonpolarColor,
		PRO: this.nonpolarColor,
		ALA: this.nonpolarColor,
		VAL: this.nonpolarColor,
		LEU: this.nonpolarColor,
		ILE: this.nonpolarColor,
		MET: this.nonpolarColor,
		PHE: this.nonpolarColor,
		CYS: this.nonpolarColor,
		TRP: this.nonpolarColor,
	},

	ssColors: {
		helix: new THREE.Color(0xFF0080),
		sheet: new THREE.Color(0xFFC800),
		 coil: new THREE.Color(0x6080FF),
	},

	defaultBondColor: new THREE.Color(0x2194D6),

	labelVertexShader: '\
uniform float width, height;\n\
varying vec2 vUv;\n\
void main()\n\
{\n\
	mat4 mv = modelViewMatrix;\n\
	mv[0][0] = mv[1][1] = mv[2][2] = 1.0;\n\
	mv[0][1] = mv[0][2] = mv[1][0] = mv[1][2] = mv[2][0] =  mv[2][1] = 0.0;\n\
	mat4 mat = projectionMatrix * mv;\n\
	vUv = uv;\n\
	float aspect = projectionMatrix[1][1] / projectionMatrix[0][0];\n\
	gl_Position = mat * vec4(position, 1.0);\n\
	gl_Position /= gl_Position.w;\n\
	gl_Position += vec4(uv.x * width * 1e-3, uv.y * height * aspect * 1e-3, 0.0, 0.0);\n\
	gl_Position.z = -0.9;\n\
}',

	labelFragmentShader: '\
uniform sampler2D map;\n\
varying vec2 vUv;\n\
void main()\n\
{\n\
	gl_FragColor = texture2D(map, vec2(vUv.x, 1.0 - vUv.y));\n\
	if (gl_FragColor.a < 0.5) discard;\n\
}',

	hasCovalentBond: function (atom0, atom1) {
		var r = this.covalentRadii[atom0.elem] + this.covalentRadii[atom1.elem];
		return atom0.coord.distanceToSquared(atom1.coord) < 1.3 * r * r;
	},

    loadPDB: function (src) {
        var helices = [], sheets = [];
        this.atoms = {};
        var lines = src.split('\n');
        for (var i in lines) {
            var line = lines[i];
            var record = line.substr(0, 6);
            if (record === 'HELIX ') {
                helices.push({
                    chain: line.substr(19, 1),
                    initialResidue: parseInt(line.substr(21, 4)),
                    initialInscode: line.substr(25, 1),
                    terminalResidue: parseInt(line.substr(33, 4)),
                    terminalInscode: line.substr(37, 1),
                });
            } else if (record === 'SHEET ') {
                sheets.push({
                    chain: line.substr(21, 1),
                    initialResidue: parseInt(line.substr(22, 4)),
                    initialInscode: line.substr(26, 1),
                    terminalResidue: parseInt(line.substr(33, 4)),
                    terminalInscode: line.substr(37, 1),
                });
            } else if (record === 'ATOM  ' || record === 'HETATM') {
                var serial = parseInt(line.substr(6, 5));
                this.atoms[serial] = {
                    het: record[0] === 'H',
                    serial: serial,
                    name: line.substr(12, 4).replace(/ /g, ''),
                    alt: line.substr(16, 1),
                    resn: line.substr(17, 3),
                    chain: line.substr(21, 1),
                    resi: parseInt(line.substr(22, 4)),
                    insc: line.substr(26, 1),
                    coord: new THREE.Vector3(parseFloat(line.substr(30, 8)), parseFloat(line.substr(38, 8)), parseFloat(line.substr(46, 8))),
                    b: parseFloat(line.substr(60, 8)),
                    elem: line.substr(76, 2).replace(/ /g, ''),
                    bonds: [],
                    ss: 'coil',
                };
            } else if (record === 'CONECT') {
                var from = parseInt(line.substr(6, 5));
                for (var j = 0; j < 4; ++j) {
                    var to = parseInt(line.substr([11, 16, 21, 26][j], 5));
                    if (isNaN(to)) continue;
                    this.atoms[from].bonds.push(this.atoms[to]);
                }
            } else if (record === 'TER   ') {
                this.lastTerSerial = parseInt(line.substr(6, 5));
            }
        }
        this.peptides = {};
        this.ligands = {};
        this.ions = {};
        this.waters = {};
        this.calphas = {};
        var curChain, curResi, curInsc, curResAtoms = [], me = this;
        var refreshBonds = function (f) {
            var n = curResAtoms.length;
            for (var j = 0; j < n; ++j) {
                var atom0 = curResAtoms[j];
                for (var k = j + 1; k < n; ++k) {
                    var atom1 = curResAtoms[k];
                    if (atom0.alt === atom1.alt && me.hasCovalentBond(atom0, atom1)) {
                        atom0.bonds.push(atom1);
                        atom1.bonds.push(atom0);
                    }
                }
                f && f(atom0);
            }
        };
        var pmin = new THREE.Vector3( 9999, 9999, 9999);
        var pmax = new THREE.Vector3(-9999,-9999,-9999);
        var psum = new THREE.Vector3();
        var cnt = 0;
        for (var i in this.atoms) {
            var atom = this.atoms[i];
            var coord = atom.coord;
            psum.add(coord);
            pmin.min(coord);
            pmax.max(coord);
            ++cnt;
            if (atom.serial <= this.lastTerSerial) {
                this.peptides[atom.serial] = atom;
                if (atom.name === 'CA') this.calphas[atom.serial] = atom;
                for (var j in helices) {
                    var helix = helices[j];
                    if (atom.chain == helix.chain && (atom.resi > helix.initialResidue || (atom.resi == helix.initialResidue && atom.insc >= helix.initialInscode)) && (atom.resi < helix.terminalResidue || (atom.resi == helix.terminalResidue && atom.insc <= helix.terminalInscode))) {
                        atom.ss = 'helix';
                        if (atom.resi == helix.initialResidue && atom.insc == helix.initialInscode) atom.ssbegin = true;
                        if (atom.resi == helix.terminalResidue && atom.insc == helix.terminalInscode) atom.ssend = true;
                    }
                }
                for (var j in sheets) {
                    var sheet = sheets[j];
                    if (atom.chain == sheet.chain && (atom.resi > sheet.initialResidue || (atom.resi == sheet.initialResidue && atom.insc >= sheet.initialInscode)) && (atom.resi < sheet.terminalResidue || (atom.resi == sheet.terminalResidue && atom.insc <= sheet.terminalInscode))) {
                        atom.ss = 'sheet';
                        if (atom.resi == sheet.initialResidue && atom.insc == sheet.initialInscode) atom.ssbegin = true;
                        if (atom.resi == sheet.terminalResidue && atom.insc == sheet.terminalInscode) atom.ssend = true;
                    }
                }
                if (atom.het) continue;
                if (!(curChain == atom.chain && curResi == atom.resi && curInsc == atom.insc)) {
                    refreshBonds(function (atom0) {
                        if (((atom0.name === 'C' && atom.name === 'N') || (atom0.name === 'O3\'' && atom.name === 'P')) && me.hasCovalentBond(atom0, atom)) {
                            atom0.bonds.push(atom);
                            atom.bonds.push(atom0);
                        }
                    });
                    curChain = atom.chain;
                    curResi = atom.resi;
                    curInsc = atom.insc;
                    curResAtoms.length = 0;
                }
                curResAtoms.push(atom);
            } else if ((this.atoms[atom.serial - 1] === undefined || this.atoms[atom.serial - 1].resi !== atom.resi) && (this.atoms[atom.serial + 1] === undefined || this.atoms[atom.serial + 1].resi !== atom.resi)) {
                if (atom.elem === 'O') {
                    this.waters[atom.serial] = atom;
                } else {
                    this.ions[atom.serial] = atom;
                }
            } else {
                this.ligands[atom.serial] = atom;
            }
        }
        refreshBonds();
        this.pmin = pmin;
        this.pmax = pmax;
        this.surfaces = {
            1: undefined,
            2: undefined,
            3: undefined,
            4: undefined,
        };
        this.rebuildScene();
        this.mdl.position.copy(psum).multiplyScalar(-1 / cnt);
        var maxD = pmax.distanceTo(pmin);
        if (maxD < 25) maxD = 25;
        this.slabNear = -maxD * 0.50;
        this.slabFar  =  maxD * 0.25;
        this.rot.position.z = maxD * 0.35 / Math.tan(Math.PI / 180.0 * 10) - 150;
        this.rot.quaternion.set(1, 0, 0, 0);
        this.render();
    },

    /*************************************************************************/

    /**
     * RTM added.
     * Zoom in, rebuildScene, and render.
     */
    zoomInto : function (options) {
		var pmin = new THREE.Vector3( 9999, 9999, 9999);
		var pmax = new THREE.Vector3(-9999,-9999,-9999);
		var psum = new THREE.Vector3();
		var cnt = 0;
		for (var i in this.atoms) {
			var atom = this.atoms[i];
			var coord = atom.coord;
			psum.add(coord);
			pmin.min(coord);
			pmax.max(coord);
			++cnt;
        }
		this.pmin = pmin;
		this.pmax = pmax;
		this.surfaces = {
			1: undefined,
			2: undefined,
			3: undefined,
			4: undefined,
		};
		this.rebuildScene(options);
		this.mdl.position.copy(psum).multiplyScalar(-1 / cnt);
		var maxD = pmax.distanceTo(pmin);
		if (maxD < 25) maxD = 25;
		this.slabNear = -maxD * 0.50;
		this.slabFar  =  maxD * 0.25;
		this.rot.position.z = maxD * 0.35 / Math.tan(Math.PI / 180.0 * 10) - 150;
		this.rot.quaternion.set(1, 0, 0, 0);
		this.render();
    },

    /**
     * RTM added.
     * load a topology. TOOD(rmcgibbo) document the schema
     */
    loadTopology : function (src) {
        this.atoms = src.atoms;
        this.peptides = {};
        this.ligands = {};
        this.ions = {};
        this.waters = {};
        this.calphas = {};
        this.lastTerSerial = src.peptideIndices.length;

        this.top = src;

        var keys = [['peptideIndices', 'peptides'],
                    ['ligandIndices', 'ligands'],
                    ['ionIndices', 'ions',],
                    ['calphaIndices', 'calphas']];
        for (var key_i = 0; key_i < keys.length; key_i++) {
            var key = keys[key_i][0];
            var dest = keys[key_i][1];
            for (var ii = 0; ii < src[key].length; ii++) {
                var i = src[key][ii];
                this[dest][this.atoms[i].serial] = this.atoms[i];
            }
        }

        for (var ii = 0; ii < src.bondIndices.length; ii++) {
            var ai = src.bondIndices[ii][0];
            var aj = src.bondIndices[ii][1];
            this.atoms[ai].bonds.push(this.atoms[aj]);
            this.atoms[aj].bonds.push(this.atoms[ai]);
        }
    },

    /**
     * RTM added.
     * Add attributes to the self.atoms dict. Look at the loadPDB method
     * to see the keys in the atoms dict.
     * `src` should be a dict whose keys are the atom index.
     */
    loadAtomAttributes : function (src) {
        $.extend(true, this.atoms, src);
    },

    /**
     * RTM added
     * Set the xyz coordinates. src should be a list 
     */
    loadCoordinates : function (src) {
        var natoms = src.length;
        for (var i = 0; i < natoms; i++) {
            // convert nanometers to angstroms
            this.atoms[i].coord = (new THREE.Vector3).fromArray(src[i])
                .multiplyScalar(10);
        }
    },

    /*************************************************************************/

	createSphere: function (atom, defaultRadius, forceDefault, scale) {
		var mesh = new THREE.Mesh(this.sphereGeometry, new THREE.MeshLambertMaterial({ color: atom.color }));
		mesh.scale.x = mesh.scale.y = mesh.scale.z = forceDefault ? defaultRadius : (this.vdwRadii[atom.elem] || defaultRadius) * (scale ? scale : 1);
		mesh.position.copy(atom.coord);
		this.mdl.add(mesh);
	},

	createCylinder: function (p0, p1, radius, color) {
		var mesh = new THREE.Mesh(this.cylinderGeometry, new THREE.MeshLambertMaterial({ color: color }));
		mesh.position.copy(p0).add(p1).multiplyScalar(0.5);
		mesh.matrixAutoUpdate = false;
		mesh.lookAt(p0);
		mesh.updateMatrix();
		mesh.matrix.multiply(new THREE.Matrix4().makeScale(radius, radius, p0.distanceTo(p1))).multiply(new THREE.Matrix4().makeRotationX(Math.PI * 0.5));
		this.mdl.add(mesh);
	},

	createRepresentationSub: function (atoms, f0, f01) {
		var ged = new THREE.Geometry();
		for (var i in atoms) {
			var atom0 = atoms[i];
			f0 && f0(atom0);
			for (var j in atom0.bonds) {
				var atom1 = atom0.bonds[j];
				if (atom1.serial < atom0.serial) continue;
				if (atom1.chain === atom0.chain && ((atom1.resi === atom0.resi) || (atom0.name === 'C' && atom1.name === 'N') || (atom0.name === 'O3\'' && atom1.name === 'P'))) {
					f01 && f01(atom0, atom1);
				} else {
					ged.vertices.push(atom0.coord);
					ged.vertices.push(atom1.coord);
				}
			}
		}
		if (ged.vertices.length) {
			ged.computeLineDistances();
			this.mdl.add(new THREE.Line(ged, new THREE.LineDashedMaterial({ linewidth: this.linewidth, color: this.defaultBondColor, dashSize: 0.25, gapSize: 0.125 }), THREE.LinePieces));
		}
	},

	createSphereRepresentation: function (atoms, defaultRadius, forceDefault, scale) {
		var me = this;
		this.createRepresentationSub(atoms, function (atom0) {
			me.createSphere(atom0, defaultRadius, forceDefault, scale);
		});
	},

	createStickRepresentation: function (atoms, atomR, bondR, scale) {
		var me = this;
		this.createRepresentationSub(atoms, function (atom0) {
			me.createSphere(atom0, atomR, !scale, scale);
		}, function (atom0, atom1) {
			if (atom0.color === atom1.color) {
				me.createCylinder(atom0.coord, atom1.coord, bondR, atom0.color);
			} else {
				var mp = atom0.coord.clone().add(atom1.coord).multiplyScalar(0.5);
				me.createCylinder(atom0.coord, mp, bondR, atom0.color);
				me.createCylinder(atom1.coord, mp, bondR, atom1.color);
			}
		});
	},

	createLineRepresentation: function (atoms) {
		var me = this;
		var geo = new THREE.Geometry();
		this.createRepresentationSub(atoms, undefined, function (atom0, atom1) {
			if (atom0.color === atom1.color) {
				geo.vertices.push(atom0.coord);
				geo.vertices.push(atom1.coord);
				geo.colors.push(atom0.color);
				geo.colors.push(atom1.color);
			} else {
				var mp = atom0.coord.clone().add(atom1.coord).multiplyScalar(0.5);
				geo.vertices.push(atom0.coord);
				geo.vertices.push(mp);
				geo.vertices.push(atom1.coord);
				geo.vertices.push(mp);
				geo.colors.push(atom0.color);
				geo.colors.push(atom0.color);
				geo.colors.push(atom1.color);
				geo.colors.push(atom1.color);
			}
		});
		this.mdl.add(new THREE.Line(geo, new THREE.LineBasicMaterial({ linewidth: this.linewidth, vertexColors: true }), THREE.LinePieces));
	},

	subdivide: function (_points, DIV) { // Catmull-Rom subdivision
		var ret = [];
		var points = new Array(); // Smoothing test
		points.push(_points[0]);
		for (var i = 1, lim = _points.length - 1; i < lim; ++i) {
			var p0 = _points[i], p1 = _points[i + 1];
			points.push(p0.smoothen ? p0.clone().add(p1).multiplyScalar(0.5) : p0);
		}
		points.push(_points[_points.length - 1]);
		for (var i = -1, size = points.length, DIVINV = 1 / DIV; i <= size - 3; ++i) {
			var p0 = points[i == -1 ? 0 : i];
			var p1 = points[i + 1], p2 = points[i + 2];
			var p3 = points[i == size - 3 ? size - 1 : i + 3];
			var v0 = p2.clone().sub(p0).multiplyScalar(0.5);
			var v1 = p3.clone().sub(p1).multiplyScalar(0.5);
			for (var j = 0; j < DIV; ++j) {
				var t = DIVINV * j;
				var x = p1.x + t * v0.x
						 + t * t * (-3 * p1.x + 3 * p2.x - 2 * v0.x - v1.x)
						 + t * t * t * (2 * p1.x - 2 * p2.x + v0.x + v1.x);
				var y = p1.y + t * v0.y
						 + t * t * (-3 * p1.y + 3 * p2.y - 2 * v0.y - v1.y)
						 + t * t * t * (2 * p1.y - 2 * p2.y + v0.y + v1.y);
				var z = p1.z + t * v0.z
						 + t * t * (-3 * p1.z + 3 * p2.z - 2 * v0.z - v1.z)
						 + t * t * t * (2 * p1.z - 2 * p2.z + v0.z + v1.z);
				ret.push(new THREE.Vector3(x, y, z));
			}
		}
		ret.push(points[points.length - 1]);
		return ret;
	},

	createCurveSub: function (_points, width, colors, div) {
		if (_points.length == 0) return;
		div = div || 5;
		var points = this.subdivide(_points, div);
		var geo = new THREE.Geometry();
		for (var i = 0, divInv = 1 / div; i < points.length; ++i) {
			geo.vertices.push(points[i]);
			geo.colors.push(new THREE.Color(colors[i == 0 ? 0 : Math.round((i - 1) * divInv)]));
		}
		this.mdl.add(new THREE.Line(geo, new THREE.LineBasicMaterial({ linewidth: width, vertexColors: true }), THREE.LineStrip));
	},

	createCurve: function (atoms, curveWidth, atomName, div) {
		var points = [], colors = [];
		var currentChain, currentResi;
		div = div || 5;
		for (var i in atoms) {
			var atom = atoms[i];
			if (atom.name == atomName && !atom.het) {
				if (currentChain != atom.chain || currentResi + 1 != atom.resi) {
					this.createCurveSub(points, curveWidth, colors, div);
					points = [];
					colors = [];
				}
				points.push(atom.coord);
				colors.push(atom.color);
				currentChain = atom.chain;
				currentResi = atom.resi;
			}
		}
		this.createCurveSub(points, curveWidth, colors, div);
	},

	createStrip: function (p0, p1, colors, div, thickness) {
		if (p0.length < 2) return;
		div = div || this.axisDIV;
		p0 = this.subdivide(p0, div);
		p1 = this.subdivide(p1, div);
		var geo = new THREE.Geometry();
		var vs = geo.vertices, fs = geo.faces;
		var axis, p0v, p1v, a0v, a1v;
		for (var i = 0, lim = p0.length; i < lim; ++i) {
			vs.push(p0v = p0[i]); // 0
			vs.push(p0v); // 1
			vs.push(p1v = p1[i]); // 2
			vs.push(p1v); // 3
			if (i < lim - 1) {
				axis = p1[i].clone().sub(p0[i]).cross(p0[i + 1].clone().sub(p0[i])).normalize().multiplyScalar(thickness);
			}
			vs.push(a0v = p0[i].clone().add(axis)); // 4
			vs.push(a0v); // 5
			vs.push(a1v = p1[i].clone().add(axis)); // 6
			vs.push(a1v); // 7
		}
		var faces = [[0, 2, -6, -8], [-4, -2, 6, 4], [7, 3, -5, -1], [-3, -7, 1, 5]];
		for (var i = 1, lim = p0.length, divInv = 1 / div; i < lim; ++i) {
			var offset = 8 * i, color = new THREE.Color(colors[Math.round((i - 1) * divInv)]);
			for (var j = 0; j < 4; ++j) {
				fs.push(new THREE.Face3(offset + faces[j][0], offset + faces[j][1], offset + faces[j][2], undefined, color));
				fs.push(new THREE.Face3(offset + faces[j][3], offset + faces[j][0], offset + faces[j][2], undefined, color));
			}
		}
		var vsize = vs.length - 8; // Cap
		for (var i = 0; i < 4; ++i) {
			vs.push(vs[i * 2]);
			vs.push(vs[vsize + i * 2]);
		};
		vsize += 8;
		fs.push(new THREE.Face3(vsize, vsize + 2, vsize + 6, undefined, fs[0].color));
		fs.push(new THREE.Face3(vsize + 4, vsize, vsize + 6, undefined, fs[0].color));
		fs.push(new THREE.Face3(vsize + 1, vsize + 5, vsize + 7, undefined, fs[fs.length - 3].color));
		fs.push(new THREE.Face3(vsize + 3, vsize + 1, vsize + 7, undefined, fs[fs.length - 3].color));
		geo.computeFaceNormals();
		geo.computeVertexNormals(false);
		this.mdl.add(new THREE.Mesh(geo, new THREE.MeshLambertMaterial({ vertexColors: THREE.FaceColors, side: THREE.DoubleSide })));
	},

	createStrand: function (atoms, num, div, fill, coilWidth, helixSheetWidth, doNotSmoothen, thickness) {
		num = num || this.strandDIV;
		div = div || this.axisDIV;
		coilWidth = coilWidth || this.coilWidth;
		doNotSmoothen = doNotSmoothen || false;
		helixSheetWidth = helixSheetWidth || this.helixSheetWidth;
		var points = {}; for (var k = 0; k < num; ++k) points[k] = [];
		var colors = [];
		var currentChain, currentResi, currentCA;
		var prevCO = null, ss = null, ssborder = false;
		for (var i in atoms) {
			var atom = atoms[i];
			if ((atom.name === 'O' || atom.name === 'CA') && !atom.het) {
				if (atom.name === 'CA') {
					if (currentChain != atom.chain || currentResi + 1 != atom.resi) {
						for (var j = 0; !thickness && j < num; ++j)
							this.createCurveSub(points[j], 1, colors, div);
						if (fill) this.createStrip(points[0], points[num - 1], colors, div, thickness);
						var points = {}; for (var k = 0; k < num; ++k) points[k] = [];
						colors = [];
						prevCO = null; ss = null; ssborder = false;
					}
					currentCA = atom.coord;
					currentChain = atom.chain;
					currentResi = atom.resi;
					ss = atom.ss; ssborder = atom.ssstart || atom.ssend;
					colors.push(atom.color);
				} else { // O
					var O = atom.coord.clone();
					O.sub(currentCA);
					O.normalize(); // can be omitted for performance
					O.multiplyScalar(ss === 'coil' ? coilWidth : helixSheetWidth);
					if (prevCO != undefined && O.dot(prevCO) < 0) O.negate();
					prevCO = O;
					for (var j = 0, numM1Inv2 = 2 / (num - 1); j < num; ++j) {
						var delta = -1 + numM1Inv2 * j;
						var v = new THREE.Vector3(currentCA.x + prevCO.x * delta, currentCA.y + prevCO.y * delta, currentCA.z + prevCO.z * delta);
						if (!doNotSmoothen && ss === 'sheet') v.smoothen = true;
						points[j].push(v);
					}
				}
			}
		}
		for (var j = 0; !thickness && j < num; ++j)
			this.createCurveSub(points[j], 1, colors, div);
		if (fill) this.createStrip(points[0], points[num - 1], colors, div, thickness);
	},

	createTubeSub: function (_points, colors, radii) {
		if (_points.length < 2) return;
		var circleDiv = this.tubeDIV, axisDiv = this.axisDIV;
		var circleDivInv = 1 / circleDiv, axisDivInv = 1 / axisDiv;
		var geo = new THREE.Geometry();
		var points = this.subdivide(_points, axisDiv);
		var prevAxis1 = new THREE.Vector3(), prevAxis2;
		for (var i = 0, lim = points.length; i < lim; ++i) {
			var r, idx = (i - 1) * axisDivInv;
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
				delta = points[i].clone().sub(points[i + 1]);
				axis1 = new THREE.Vector3(0, -delta.z, delta.y).normalize().multiplyScalar(r);
				axis2 = delta.clone().cross(axis1).normalize().multiplyScalar(r);
				//      var dir = 1, offset = 0;
				if (prevAxis1.dot(axis1) < 0) {
					axis1.negate(); axis2.negate();  //dir = -1;//offset = 2 * Math.PI / axisDiv;
				}
				prevAxis1 = axis1; prevAxis2 = axis2;
			} else {
				axis1 = prevAxis1; axis2 = prevAxis2;
			}
			for (var j = 0; j < circleDiv; ++j) {
				var angle = 2 * Math.PI * circleDivInv * j; //* dir  + offset;
				geo.vertices.push(points[i].clone().add(axis1.clone().multiplyScalar(Math.cos(angle))).add(axis2.clone().multiplyScalar(Math.sin(angle))));
			}
		}
		var offset = 0;
		for (var i = 0, lim = points.length - 1; i < lim; ++i) {
			var c = new THREE.Color(colors[Math.round((i - 1) * axisDivInv)]);
			var reg = 0;
			var r1 = geo.vertices[offset].clone().sub(geo.vertices[offset + circleDiv]).lengthSq();
			var r2 = geo.vertices[offset].clone().sub(geo.vertices[offset + circleDiv + 1]).lengthSq();
			if (r1 > r2) { r1 = r2; reg = 1; };
			for (var j = 0; j < circleDiv; ++j) {
				geo.faces.push(new THREE.Face3(offset + j, offset + (j + reg) % circleDiv + circleDiv, offset + (j + 1) % circleDiv, undefined, c));
				geo.faces.push(new THREE.Face3(offset + (j + 1) % circleDiv, offset + (j + reg) % circleDiv + circleDiv, offset + (j + reg + 1) % circleDiv + circleDiv, undefined, c));
			}
			offset += circleDiv;
		}
		geo.computeFaceNormals();
		geo.computeVertexNormals(false);
		this.mdl.add(new THREE.Mesh(geo, new THREE.MeshLambertMaterial({ vertexColors: THREE.FaceColors, side: THREE.DoubleSide })));
	},

	createTube: function (atoms, atomName, radius) {
		var points = [], colors = [], radii = [];
		var currentChain, currentResi;
		for (var i in atoms) {
			var atom = atoms[i];
			if ((atom.name == atomName) && !atom.het) {
				if (currentChain != atom.chain || currentResi + 1 != atom.resi) {
					this.createTubeSub(points, colors, radii);
					points = []; colors = []; radii = [];
				}
				points.push(atom.coord);
				radii.push(radius || (atom.b > 0 ? atom.b * 0.01 : 0.3));
				colors.push(atom.color);
				currentChain = atom.chain;
				currentResi = atom.resi;
			}
		}
		this.createTubeSub(points, colors, radii);
	},

	createCylinderHelix: function (atoms, radius) {
		var start = null;
		var currentChain, currentResi;
		var others = {}, beta = {};
		for (var i in atoms) {
			var atom = atoms[i];
			if (atom.het) continue;
			if ((atom.ss != 'helix' && atom.ss != 'sheet') || atom.ssend || atom.ssbegin) others[atom.serial] = atom;
			if (atom.ss === 'sheet') beta[atom.serial] = atom;
			if (atom.name != 'CA') continue;
			if (atom.ss === 'helix' && atom.ssend) {
				if (start != null) this.createCylinder(start.coord, atom.coord, radius, atom.color, true);
				start = null;
			}
			currentChain = atom.chain;
			currentResi = atom.resi;
			if (start == null && atom.ss === 'helix' && atom.ssbegin) start = atom;
		}
		if (start != null) this.createCylinder(start.coord, atom.coord, radius, atom.color);
		this.createTube(others, 'CA', 0.3);
		this.createStrand(beta, undefined, undefined, true, 0, this.helixSheetWidth, false, this.thickness * 2);
	},

	createSurfaceRepresentation: function (atoms, type, wireframe, opacity) {
		if (this.surfaces[type] === undefined) {
			var ps = ProteinSurface({
				min: this.pmin,
				max: this.pmax,
				atoms: atoms,
				type: type,
			});
			var verts = ps.verts;
			var faces = ps.faces;
			var geo = new THREE.Geometry();
			geo.vertices = verts.map(function (v) {
				return new THREE.Vector3(v.x, v.y, v.z);
			});
			geo.faces = faces.map(function (f) {
				return new THREE.Face3(f.a, f.b, f.c, undefined, Object.keys(f).map(function (d) {
					return atoms[verts[f[d]].atomid].color;
				}));
			});
			geo.computeFaceNormals();
			geo.computeVertexNormals(false);
			this.surfaces[type] = geo;
		}
		this.mdl.add(new THREE.Mesh(this.surfaces[type], new THREE.MeshLambertMaterial({
			vertexColors: THREE.VertexColors,
			wireframe: wireframe,
			opacity: opacity,
			transparent: true,
		})));
	},

	createLabel: function (text, size, color) {
		var canvas = document.createElement('canvas');
		canvas.style.backgroundColor = 'rgba(0, 0, 0, 0.0)';
		var ctx = canvas.getContext('2d');
		ctx.font = size + 'px Arial';
		canvas.width = ctx.measureText(text).width;
		canvas.height = size;
		ctx.font = size + 'px Arial';
		ctx.fillStyle = color;
		ctx.fillText(text, 0, size);
		var tex = new THREE.Texture(canvas);
		tex.flipY = false;
		tex.needsUpdate = true;
		return new THREE.Mesh(this.labelGeo, new THREE.ShaderMaterial({
			vertexShader: this.labelVertexShader,
			fragmentShader: this.labelFragmentShader,
			uniforms: {
				map: {
					type: 't',
					value: tex,
				},
				width: {
					type: 'f',
					value: tex.image.width,
				},
				height: {
					type: 'f',
					value: tex.image.height,
				},
			},
		}));
	},

	createLabelRepresentation: function (atoms, label) {
		for (var i in atoms) {
			var atom = atoms[i];
			var bb = this.createLabel(label(atom), 20, '#dddddd');
			bb.position.copy(atom.coord);
			this.mdl.add(bb);
		}
	},

	rebuildScene: function (options) {
		$.extend(this.options, options);
		this.scene = new THREE.Scene();
		var directionalLight = new THREE.DirectionalLight(0xFFFFFF, 1.2);
		directionalLight.position.set(0.2, 0.2, -1).normalize();
		var ambientLight = new THREE.AmbientLight(0x202020);
		this.scene.add(directionalLight);
		this.scene.add(ambientLight);
		var mp = this.mdl ? this.mdl.position : new THREE.Vector3();
		var rz = this.rot ? this.rot.position.z : 0;
		var rq = this.rot ? this.rot.quaternion : new THREE.Quaternion();
		this.mdl = new THREE.Object3D();
		this.mdl.position.copy(mp);
		this.rot = new THREE.Object3D();
		this.rot.position.z = rz;
		this.rot.quaternion.copy(rq);
		this.rot.add(this.mdl);
		this.scene.add(this.rot);
		this.camera = this.cameras[this.options.camera];
		var background = this.backgroundColors[this.options.background];
		this.renderer.setClearColor(background);
		this.scene.fog = new THREE.Fog(background, 100, 200);
		switch (this.options.colorBy) {
			case 'spectrum':
				var idx = 0, lastTerSerialInv = 1 / this.lastTerSerial;
				for (var i in this.atoms) {
					var atom = this.atoms[i];
					atom.color = atom.het ? this.atomColors[atom.elem] || this.defaultAtomColor : new THREE.Color().setHSL(2 / 3 * (1 - idx++ * lastTerSerialInv), 1, 0.45);
				}
				break;
			case 'chain':
				for (var i in this.atoms) {
					var atom = this.atoms[i];
					atom.color = (atom.het ? this.hetChainColors : this.stdChainColors)[atom.chain];
				}
				break;
			case 'secondary structure':
				for (var i in this.atoms) {
					var atom = this.atoms[i];
					atom.color = atom.het ? this.atomColors[atom.elem] || this.defaultAtomColor : this.ssColors[atom.ss];
				}
				break;
			case 'B factor':
				if (!this.middB) {
					var minB = 1000, maxB = -1000;
					for (var i in this.atoms) {
						var atom = this.atoms[i];
						if (minB > atom.b) minB = atom.b;
						if (maxB < atom.b) maxB = atom.b;
					}
					this.middB = (maxB + minB) * 0.5;
					this.spanB = (maxB - minB) * 0.5;
					this.spanBinv = 1 / this.spanB;
				}
				for (var i in this.atoms) {
					var atom = this.atoms[i];
					atom.color = atom.b < this.middB ? new THREE.Color().setRGB(1 - (s = (this.middB - atom.b) * this.spanBinv), 1 - s, 1) : new THREE.Color().setRGB(1, 1 - (s = (atom.b - this.middB) * this.spanBinv), 1 - s);
				}
				break;
			case 'residue':
				for (var i in this.atoms) {
					var atom = this.atoms[i];
					atom.color = atom.het ? this.atomColors[atom.elem] || this.defaultAtomColor : this.residueColors[atom.resn] || this.defaultResidueColor;
				}
				break;
			case 'polarity':
				for (var i in this.atoms) {
					var atom = this.atoms[i];
					atom.color = atom.het ? this.atomColors[atom.elem] || this.defaultAtomColor : this.polarityColors[atom.resn] || this.defaultResidueColor;
				}
				break;
			case 'atom':
				for (var i in this.atoms) {
					var atom = this.atoms[i];
					atom.color = this.atomColors[atom.elem] || this.defaultAtomColor;
				}
				break;
		}

		switch (this.options.primaryStructure) {
			case 'lines':
				this.createLineRepresentation(this.peptides);
				break;
			case 'stick':
				this.createStickRepresentation(this.peptides, this.cylinderRadius, this.cylinderRadius);
				break;
			case 'ball & stick':
				this.createStickRepresentation(this.peptides, this.cylinderRadius, this.cylinderRadius * 0.5, 0.3);
				break;
			case 'sphere':
				this.createSphereRepresentation(this.peptides, this.sphereRadius);
				break;
		}

		var doNotSmoothen = false;
		switch (this.options.secondaryStructure) {
			case 'ribbon':
				this.createStrand(this.peptides, 2, undefined, true, undefined, undefined, doNotSmoothen, this.thickness);
				break;
			case 'strand':
				this.createStrand(this.peptides, null, null, null, null, null, doNotSmoothen);
				break;
			case 'cylinder & plate':
				this.createCylinderHelix(this.peptides, 1.6);
				break;
			case 'C alpha trace':
				this.createCurve(this.peptides, this.curveWidth, 'CA', 1);
				break;
			case 'B factor tube':
				this.createTube(this.peptides, 'CA');
				break;
		}

		switch (this.options.wireframe) {
			case 'yes':
				this.options.wireframe = true;
				break;
			case 'no':
				this.options.wireframe = false;
				break;
		}

		this.options.opacity = parseFloat(this.options.opacity);

		switch (this.options.surface) {
			case 'Van der Waals surface':
				this.createSurfaceRepresentation(this.peptides, 1, this.options.wireframe, this.options.opacity);
				break;
			case 'solvent excluded surface':
				this.createSurfaceRepresentation(this.peptides, 2, this.options.wireframe, this.options.opacity);
				break;
			case 'solvent accessible surface':
				this.createSurfaceRepresentation(this.peptides, 3, this.options.wireframe, this.options.opacity);
				break;
			case 'molecular surface':
				this.createSurfaceRepresentation(this.peptides, 4, this.options.wireframe, this.options.opacity);
				break;
		}

		switch (this.options.ligands) {
			case 'line':
				this.createLineRepresentation(this.ligands, this.curveWidth);
				break;
			case 'stick':
				this.createStickRepresentation(this.ligands, this.cylinderRadius, this.cylinderRadius);
				break;
			case 'ball & stick':
				this.createStickRepresentation(this.ligands, this.cylinderRadius, this.cylinderRadius * 0.5, 0.3);
				break;
			case 'sphere':
				this.createSphereRepresentation(this.ligands, this.sphereRadius);
				break;
		}

		switch (this.options.ions) {
			case 'sphere':
				this.createSphereRepresentation(this.ions, this.sphereRadius);
				break;
			case 'dot':
				this.createSphereRepresentation(this.ions, this.sphereRadius, false, 0.3);
				break;
		}

		switch (this.options.waters) {
			case 'sphere':
				this.createSphereRepresentation(this.waters, this.sphereRadius);
				break;
			case 'dot':
				this.createSphereRepresentation(this.waters, this.sphereRadius, false, 0.3);
				break;
		}

		switch (this.options.labels) {
			case 'yes':
				this.createLabelRepresentation(this.calphas, this.labelByResn);
				this.createLabelRepresentation(this.ligands, this.labelByName);
				this.createLabelRepresentation(this.ions, this.labelByName);
				break;
		}

		this.renderer.autoClear = this.options.effect !== 'oculus rift' && this.options.effect !== 'stereo';
		this.effect = this.effects[this.options.effect];
		this.effect.setSize(this.container.width(), this.container.height());
	},

	render: function () {
		var center = this.rot.position.z - this.camera.position.z;
		if (center < 1) center = 1;
		this.camera.near = center + this.slabNear;
		if (this.camera.near < 1) this.camera.near = 1;
		this.camera.far = center + this.slabFar;
		if (this.camera.near + 1 > this.camera.far) this.camera.far = this.camera.near + 1;
		if (this.camera === this.orthographicCamera){
			this.camera.right = center * Math.tan(Math.PI / 180 * 20);
			this.camera.left = -this.camera.right;
			this.camera.top = this.camera.right * this.container.hwratio;
			this.camera.bottom = -this.camera.top;
		}
		this.camera.updateProjectionMatrix();
		this.scene.fog.near = this.camera.near + 0.4 * (this.camera.far - this.camera.near);
		//if (this.scene.fog.near > center) this.scene.fog.near = center;
		this.scene.fog.far = this.camera.far;
		this.effect.render(this.scene, this.camera);
		if (!this.effect.init) {
			this.effect.render(this.scene, this.camera);
			this.effect.init = true;
		}
	},

	exportCanvas: function () {
		this.render();
		window.open(this.renderer.domElement.toDataURL('image/png'));
	},

};
