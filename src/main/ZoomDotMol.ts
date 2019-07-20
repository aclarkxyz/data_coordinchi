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

namespace WebMolKit /* BOF */ {

/*
	Displays a single molecule with interactive dot-path display.
*/

interface PathOutline
{
	pblk:DotPathBlock;
	px:number[];
	py:number[];
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
		$(this.renderOutlines(this.outlines, w, h, 0x808080, MetaVector.NOCOLOUR).createSVG()).appendTo(spanBackdrop);

		// draw each path as an on/off toggle
		for (let out of this.outlines)
		{
			let span = $('<div></div>').appendTo(this.divOutline);
			span.css({'position': 'absolute', 'left': '0px', 'top': '0px', 'display': 'none'});
			$(this.renderOutlines([out], w, h, 0x000000, 0xC0C0C0).createSVG()).appendTo(span);
			this.spanOutlines.push(span);
		}

		// molecule goes on top
		let spanMol = $('<div></div>').appendTo(this.divOutline);
		spanMol.css({'position': 'absolute', 'left': '0px', 'top': '0px'});
		$(gfx.createSVG()).appendTo(spanMol);

		this.divOutline.mouseenter((event:JQueryMouseEventObject) => this.hoverMouse(event));
		this.divOutline.mouseleave((event:JQueryMouseEventObject) => this.hoverMouse(null));
		this.divOutline.mousemove((event:JQueryMouseEventObject) => this.hoverMouse(event));
	}

	// creates an outline 
	private createOutline(pblk:DotPathBlock, layout:ArrangeMolecule):PathOutline
	{
		let x:number[] = [], y:number[] = [];
		for (let a of pblk.atoms)
		{
			let pt = layout.getPoint(a - 1);
			let rad = Math.max(0.5 * this.policy.data.pointScale, Math.max(pt.oval.rw, pt.oval.rh));
			const NPT = 36, THPT = TWOPI / NPT;
			for (let n = 0; n < NPT; n++)
			{
				let th = n * THPT;
				x.push(pt.oval.cx + rad * Math.cos(th));
				y.push(pt.oval.cy + rad * Math.sin(th));
			}
		}
		let [px, py] = GeomUtil.convexHull(x, y, 2);

		return {'pblk': pblk, 'px': px, 'py': py};
	}

	// draw one or more outlines onto a vector canvas
	private renderOutlines(outlines:PathOutline[], w:number, h:number, edgeCol:number, fillCol:number):MetaVector
	{
		let gfx = new MetaVector();
		gfx.width = w;
		gfx.height = h;
		for (let out of outlines) gfx.drawPoly(out.px, out.py, edgeCol, 1, fillCol, false);
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
}

/* EOF */ }