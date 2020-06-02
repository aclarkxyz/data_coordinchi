/*
	Coordination Complexes for InChI

	(c) 2019 InChI Trust
*/

///<reference path='../../../WebMolKit/src/decl/corrections.d.ts'/>
///<reference path='../../../WebMolKit/src/decl/jquery/index.d.ts'/>
///<reference path='../../../WebMolKit/src/util/util.ts'/>
///<reference path='../../../WebMolKit/src/util/Geom.ts'/>
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
	Renders a "card" that displays a structure and various information about it, into a DOM object.
*/

export class MoleculeCard
{
	public dom:JQuery = null; // resulting object for the whole card
	public spanMol:JQuery = null; // sub-area for the molecule rendering

	// ------------ public methods ------------

	constructor(public mol:Molecule, public molExpanded:Molecule, public inchi:string, public dhash:string, public dimsz:number, public policy:RenderPolicy)
	{
	}

	// builds the card; the resulting content is found within the 'dom' member
	public generate():void
	{
		this.dom = $('<div></div>');
		this.dom.css({'display': 'inline-block', 'margin': '0.5em', 'padding': '0.5em'});

		let divMol = $('<div></div>').appendTo(this.dom).css({'text-align': 'center'});
		let measure = new OutlineMeasurement(0, 0, this.policy.data.pointScale);		
		let layout = new ArrangeMolecule(this.mol, measure, this.policy);
		layout.arrange();
		layout.squeezeInto(0, 0, this.dimsz, this.dimsz);
		let gfx = new MetaVector();
		new DrawMolecule(layout, gfx).draw();
		gfx.normalise();
		this.spanMol = $('<span></span>').appendTo(divMol);
		this.spanMol.css({'display': 'inline-block', 'font-size': '0' /* only way to squish the baseline gap */});
		$(gfx.createSVG()).appendTo(this.spanMol);

		let divFormula = $('<div></div>').appendTo(this.dom).css({'text-align': 'center', 'font-size': '70%', 'font-weight': 'bold'});
		divFormula.html(MolUtil.molecularFormula(this.molExpanded, ['<sub>', '</sub>', '<sup>', '</sup>']));
		let chg = 0;
		for (let n = 1; n <= this.molExpanded.numAtoms; n++) chg += this.molExpanded.atomCharge(n);
		divFormula.append(' [a=' + this.molExpanded.numAtoms + ',b=' + this.molExpanded.numBonds + ',c=' + chg + ']');

		if (this.inchi)
		{
			let maxWidth = Math.max(this.dimsz, gfx.boundHighX() - gfx.boundLowX());
			let divInChI = $('<div></div>').appendTo(this.dom).css({'text-align': 'left', 'font-size': '70%', 'max-width': maxWidth + 'px', 'word-wrap': 'break-word'});

			let bits = /^(InChI=1S?\/)([\w\d\.]+)(\/.*)$/.exec(this.inchi);
			if (!bits)
			{
				let span = $('<span></span>').appendTo(divInChI);
				span.css({'color': 'white', 'background-color': '#4E1A09'});
				span.text(this.inchi);
			}
			else if (!this.sameFormula(bits[2], this.mol))
			{
				divInChI.append(escapeHTML(bits[1]));
				let span = $('<span></span>').appendTo(divInChI);
				span.css({'background-color': '#E0E000'});
				span.text(bits[2]);
				divInChI.append(escapeHTML(bits[3]));
			}
			else // all good
			{
				divInChI.text(this.inchi);
			}
		}
		if (this.dhash)
		{
			let maxWidth = Math.max(this.dimsz, gfx.boundHighX() - gfx.boundLowX());
			let divHash = $('<div></div>').appendTo(this.dom);
			divHash.css({'text-align': 'left', 'font-size': '70%', 'max-width': maxWidth + 'px', 'word-wrap': 'break-word', 'margin-top': '0.1em'});
			divHash.css({'border-top': '1px solid #C0C0C0'});
			divHash.text(this.dhash);
		}
	}

	// ------------ private methods ------------

	// returns true if the InChI-derived formula is the same as the molecule's formula
	private sameFormula(formula:string, mol:Molecule):boolean
	{
		let makeMap = (str:string):{[id:string] : number} =>
		{
			let map:{[id:string] : number} = {};
			let mul = 1;
			while (str)
			{			
				let grp = /^[\.\s]+(.*)$/.exec(str);
				if (grp) 
				{
					str = grp[1];
					mul = 1;
					continue;
				}
				grp = /^([A-Z][a-z]?)(\d+)(.*)/.exec(str);
				if (grp)
				{
					let el = grp[1], num = parseInt(grp[2]) * mul;
					map[el] = (map[el] || 0) + num;
					str = grp[3];
					continue;
				}
				grp = /^([A-Z][a-z]?)(.*)/.exec(str);
				if (grp)
				{
					let el = grp[1], num = 1 * mul;
					map[el] = (map[el] || 0) + num;
					str = grp[2];
					continue;
				}
				grp = /^(\d+)(.*)/.exec(str);
				if (grp)
				{
					mul = parseInt(grp[1]);
					str = grp[2];
					continue;
				}
				throw 'Unparseable formula: ' + str;
			}
			return map;
		};

		let map1 = makeMap(formula);
		let map2 = makeMap(MolUtil.molecularFormula(mol));
		if (map1 == null || map2 == null) return false;
		if (Object.keys(map1).length != Object.keys(map2).length) return false;
		if (map1.size != map2.size) return false;
		for (let el in map1) if (map1[el] != map2[el])
		{
			/*console.log('M1:'+formula);
			console.log('    '+JSON.stringify(map1));
			console.log('M2:'+MolUtil.molecularFormula(mol));
			console.log('    '+JSON.stringify(map2));
			throw "fnord";*/
			return false;
		}

		return true;
	}

}

/* EOF */ }