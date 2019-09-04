/*
	Coordination Complexes for InChI

	(c) 2019 InChI Trust
*/

///<reference path='../../../WebMolKit/src/decl/corrections.d.ts'/>
///<reference path='../../../WebMolKit/src/decl/jquery.d.ts'/>
///<reference path='../../../WebMolKit/src/util/util.ts'/>
///<reference path='../../../WebMolKit/src/util/Geom.ts'/>
///<reference path='../../../WebMolKit/src/ui/Widget.ts'/>
///<reference path='../../../WebMolKit/src/dialog/Dialog.ts'/>
///<reference path='../../../WebMolKit/src/data/Molecule.ts'/>
///<reference path='../../../WebMolKit/src/data/MoleculeStream.ts'/>
///<reference path='../../../WebMolKit/src/data/DotPath.ts'/>
///<reference path='../../../WebMolKit/src/gfx/Rendering.ts'/>
///<reference path='../../../WebMolKit/src/gfx/ArrangeMolecule.ts'/>
///<reference path='../../../WebMolKit/src/gfx/DrawMolecule.ts'/>

///<reference path='../decl/node.d.ts'/>
///<reference path='../decl/electron.d.ts'/>
///<reference path='../data/DotHash.ts'/>

namespace WebMolKit /* BOF */ {

/*
	Displays a single molecule with interactive dot-path display.
*/

interface PathOutline
{
	pblk:DotPathBlock;
	px:number[];
	py:number[];
	//dx:number[];
	//dy:number[];
	ax:number[];
	ay:number[];
	ar:number[];
}

export class ZoomDotMol extends Dialog
{
	private policy:RenderPolicy;
	private effects:RenderEffects;
	private measure:ArrangeMeasurement;

	private dotpath:DotPath;
	private outlines:PathOutline[] = [];

	private divOutline:JQuery;
	private spanOutlines:JQuery[] = []; // one per outline for toggling on/off

	// ------------ public methods ------------

	constructor(private mol:Molecule)
	{
		super();

		this.title = 'Dot Paths';

		this.policy = RenderPolicy.defaultColourOnWhite();
		this.policy.data.pointScale = 35;
		this.effects = new RenderEffects();
		this.measure = new OutlineMeasurement(0, 0, this.policy.data.pointScale);
	}

	// ------------ private methods ------------

	protected populate():void
	{
		MolUtil.expandAbbrevs(this.mol, true);

		let layout = new ArrangeMolecule(this.mol, this.measure, this.policy, this.effects);
		layout.arrange();
		layout.squeezeInto(0, 0, 700, 700);
		let padding = 0.5 * this.policy.data.pointScale, bounds = layout.determineBoundary();
		layout.offsetEverything(-bounds[0] + padding, -bounds[1] + padding);
		let gfx = new MetaVector();
		new DrawMolecule(layout, gfx).draw();
		gfx.width = bounds[2] - bounds[0] + 2 * padding;
		gfx.height = bounds[3] - bounds[1] + 2 * padding;
	
		this.dotpath = new DotPath(this.mol);
		for (let pblk of this.dotpath.paths) this.outlines.push(this.createOutline(pblk, layout));	

		let w = Math.ceil(gfx.width), h = Math.ceil(gfx.height);

		let div = $('<div></div>').appendTo(this.body()).css({'text-align': 'center'});
		this.divOutline = $('<div></div>').appendTo(div);
		this.divOutline.css({'display': 'inline-block', 'position': 'relative', 'width': `${w}px`, 'height': `${h}px`});

		// draw all paths faintly
		let spanBackdrop = $('<div></div>').appendTo(this.divOutline);
		spanBackdrop.css({'position': 'absolute', 'left': '0px', 'top': '0px'});
		$(this.renderOutlines(this.outlines, w, h, 0x808080, MetaVector.NOCOLOUR, false).createSVG()).appendTo(spanBackdrop);

		// draw each path as an on/off toggle
		for (let out of this.outlines)
		{
			let span = $('<div></div>').appendTo(this.divOutline);
			span.css({'position': 'absolute', 'left': '0px', 'top': '0px', 'display': 'none'});
			$(this.renderOutlines([out], w, h, 0x000000, 0xC0C0C0, true).createSVG()).appendTo(span);
			this.spanOutlines.push(span);
		}

		// molecule goes on top
		let spanMol = $('<div></div>').appendTo(this.divOutline);
		spanMol.css({'position': 'absolute', 'left': '0px', 'top': '0px'});
		$(gfx.createSVG()).appendTo(spanMol);

		this.divOutline.mouseenter((event:JQueryMouseEventObject) => this.hoverMouse(event));
		this.divOutline.mouseleave((event:JQueryMouseEventObject) => this.hoverMouse(null));
		this.divOutline.mousemove((event:JQueryMouseEventObject) => this.hoverMouse(event));

		// now render as text

		let table = $('<table></table>').appendTo(div);
		let tr = $('<tr></tr>').appendTo(table);
		tr.append('<th>Row</th>');
		tr.append('<th>Atoms</th>');
		tr.append('<th>Bonds</th>');
		tr.append('<th>Electrons</th>');
		tr.append('<th>Charge</th>');
		let row = 0;
		for (let pblk of this.dotpath.paths)
		{
			row++;
			tr = $('<tr></tr>').appendTo(table);
			let thRow = $('<th></th>').appendTo(tr);
			let tdAtoms = $('<td></td>').appendTo(tr), tdBonds = $('<td></td>').appendTo(tr);
			let tdElectrons = $('<td></td>').appendTo(tr), tdCharge = $('<td></td>').appendTo(tr);

			let chg = 0;
			for (let a of pblk.atoms) chg += this.mol.atomCharge(a);

			thRow.text(row.toString());
			tdAtoms.text(JSON.stringify(pblk.atoms));
			tdBonds.text(JSON.stringify(pblk.bonds));
			tdElectrons.text(pblk.numer + ' / ' + pblk.denom);
			tdCharge.text(chg.toString());
		}

	}

	// creates an outline 
	private createOutline(pblk:DotPathBlock, layout:ArrangeMolecule):PathOutline
	{
		let x:number[] = [], y:number[] = [];
		let scale = this.policy.data.pointScale;
		let ax:number[] = [], ay:number[] = [], ar:number[] = [];

		for (let a of pblk.atoms)
		{
			let pt = layout.getPoint(a - 1);
			let rad = Math.max(0.5 * scale, Math.max(pt.oval.rw, pt.oval.rh));
			const NPT = 36, THPT = TWOPI / NPT;
			for (let n = 0; n < NPT; n++)
			{
				let th = n * THPT;
				x.push(pt.oval.cx + rad * Math.cos(th));
				y.push(pt.oval.cy + rad * Math.sin(th));
			}

			ax.push(pt.oval.cx);
			ay.push(pt.oval.cy);
			ar.push(rad);
		}

		for (let n = 1; n <= this.mol.numBonds; n++)
		{
			let bfr = this.mol.bondFrom(n), bto = this.mol.bondTo(n);
			if (pblk.atoms.indexOf(bfr) < 0 || pblk.atoms.indexOf(bto) < 0) continue;
			let pt1 = layout.getPoint(bfr - 1), pt2 = layout.getPoint(bto - 1);
			let x1 = pt1.oval.cx, y1 = pt1.oval.cy, x2 = pt2.oval.cx, y2 = pt2.oval.cy;
			let dx = x2 - x1, dy = y2 - y1, d = norm_xy(dx, dy), invD = invZ(d);
			let ox = dy * invD * 0.3 * scale, oy = -dx * invD * 0.3 * scale;
			let npWidth = Math.ceil(2 * d / scale) + 1, npHeight = Math.ceil(2 * norm_xy(ox, oy) / scale) + 1;
			for (let n = 0; n <= npWidth; n++)
			{
				x.push(x1 - ox + dx * n / npWidth);
				y.push(y1 - oy + dy * n / npWidth);
				x.push(x1 + ox + dx * n / npWidth);
				y.push(y1 + oy + dy * n / npWidth);
			}
			for (let n = 1; n < npHeight; n++)
			{
				x.push(x1 - ox + 2 * ox * n / npHeight);
				y.push(y1 - oy + 2 * oy * n / npHeight);
				x.push(x2 - ox + 2 * ox * n / npHeight);
				y.push(y2 - oy + 2 * oy * n / npHeight);
			}
		}

		//this.removeDegenerates(x, y);

		let [px, py] = GeomUtil.outlinePolygon(x, y, 0.5 * scale);

		return {'pblk': pblk, 'px': px, 'py': py, 'ax': ax, 'ay': ay, 'ar': ar};
	}

	// draw one or more outlines onto a vector canvas (to be activated as background)
	private renderOutlines(outlines:PathOutline[], w:number, h:number, edgeCol:number, fillCol:number, withAtoms:boolean):MetaVector
	{
		let gfx = new MetaVector();
		gfx.width = w;
		gfx.height = h;
		for (let out of outlines)
		{
			gfx.drawPoly(out.px, out.py, edgeCol, 1, fillCol, false);
			//for (let n = 0; n < out.dx.length; n++) gfx.drawOval(out.dx[n], out.dy[n], 1, 1, MetaVector.NOCOLOUR, 0, 0xFF0000);

			if (withAtoms) for (let n = 0; n < out.ar.length; n++) 
				gfx.drawOval(out.ax[n], out.ay[n], out.ar[n], out.ar[n], MetaVector.NOCOLOUR, 0, 0xE0000000);
		}

		return gfx;
	}

	// toggle on one-or-none of the block outlines
	private hoverMouse(event:JQueryMouseEventObject):void
	{
		let which = -1;
		if (event)
		{
			let x = event.pageX - this.divOutline.offset().left;
			let y = event.pageY - this.divOutline.offset().top;
			for (let n = 0; n < this.outlines.length; n++)
			{
				let out = this.outlines[n];
				if (GeomUtil.pointInPolygon(x, y, out.px, out.py)) {which = n; break;}
			}
		}
		for (let n = 0; n < this.outlines.length; n++) this.spanOutlines[n].css('display', n == which ? 'inline' : 'none');
	}

	// clip out points that are too close to each other for comfort
	private removeDegenerates(x:number[], y:number[]):void
	{
		// it's O(N^2), could be improved a lot
		const sz = x.length;
		const threshSq = sqr(2);
		let mask = Vec.booleanArray(true, sz);
		for (let i = 0; i < sz - 1; i++) if (mask[i])
			for (let j = i + 1; j < sz; j++) if (mask[j])
		{
			let dsq = norm2_xy(x[i] - x[j], y[i] - y[j]);
			if (dsq < threshSq) mask[j] = false;
		}
		//console.log('MASK:'+mask);		
		for (let n = sz - 1; n >= 0; n--) if (!mask[n])
		{
			x.splice(n, 1);
			y.splice(n, 1);
		}
		/*let z = '';		
		for (let n = 0; n < x.length; n++) z += ', {' + x[n].toFixed(2) + 'f,' + y[n].toFixed(2) + 'f}';
		console.log('PP:'+z);*/
	}
}

/* EOF */ }