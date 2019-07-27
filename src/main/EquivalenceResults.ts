/*
	Coordination Complexes for InChI

	(c) 2019 InChI Trust
*/

///<reference path='../../../WebMolKit/src/decl/corrections.d.ts'/>
///<reference path='../../../WebMolKit/src/decl/jquery.d.ts'/>
///<reference path='../../../WebMolKit/src/util/util.ts'/>
///<reference path='../../../WebMolKit/src/sketcher/Sketcher.ts'/>
///<reference path='../../../WebMolKit/src/ui/ClipboardProxy.ts'/>
///<reference path='../../../WebMolKit/src/ui/Widget.ts'/>
///<reference path='../../../WebMolKit/src/data/Molecule.ts'/>
///<reference path='../../../WebMolKit/src/data/MoleculeStream.ts'/>
///<reference path='../../../WebMolKit/src/data/MDLWriter.ts'/>
///<reference path='../../../WebMolKit/src/gfx/Rendering.ts'/>
///<reference path='../../../WebMolKit/src/gfx/ArrangeMolecule.ts'/>
///<reference path='../../../WebMolKit/src/gfx/DrawMolecule.ts'/>

///<reference path='../decl/node.d.ts'/>
///<reference path='../decl/electron.d.ts'/>

///<reference path='../data/AnalyseMolecule.ts'/>
///<reference path='ZoomDotMol.ts'/>

namespace WebMolKit /* BOF */ {

/*
	Equivalence results: reads a datasheet that has rows of structures that are supposed to be the same as each other, but not to any
	of the other rows.
*/

interface EquivalenceRow
{
	molList:Molecule[];
	inchiList:string[];
	dhashList:string[];
}

export class EquivalenceResults
{
	private divContent:JQuery;
	private tableRows:JQuery;
	
	private colMol:number[] = [];
	private colInChI:number[] = [];
	private rows:EquivalenceRow[] = [];

	protected policy:RenderPolicy;
	protected effects:RenderEffects;
	protected measure:ArrangeMeasurement;

	// ------------ public methods ------------

	constructor(private ds:DataSheet, private callbackDone:() => void)
	{
	}

	public render(parent:JQuery):void
	{
		this.divContent = $('<div></div>').appendTo(parent);

		this.policy = RenderPolicy.defaultColourOnWhite();
		this.policy.data.pointScale = 15;
		this.effects = new RenderEffects();
		this.measure = new OutlineMeasurement(0, 0, this.policy.data.pointScale);		

		this.startAnalysis();
	}

	// ------------ private methods ------------

	private startAnalysis():void
	{
		for (let n = 0; n < this.ds.numCols; n++) if (this.ds.colType(n) == DataSheet.COLTYPE_MOLECULE)
		{
			this.colMol.push(n);
			this.colInChI.push(this.ds.findColByName(this.ds.colName(n) + 'InChI', DataSheet.COLTYPE_STRING));
		}
		for (let n = 0; n < this.ds.numRows; n++)
		{
			let eqr:EquivalenceRow = {'molList': [], 'inchiList': [], 'dhashList': []};
			for (let i = 0; i < this.colMol.length; i++) if (this.ds.notNull(n, this.colMol[i]))
			{
				let mol = this.ds.getMolecule(n, this.colMol[i]);
				eqr.molList.push(mol);
				eqr.inchiList.push(this.ds.getString(n, this.colInChI[i]));
				eqr.dhashList.push(new DotHash(new DotPath(mol)).calculate());
			}
			this.rows.push(eqr);
		}

		this.tableRows = $('<table></table>').appendTo(this.divContent);

		this.processRow(0);
	}

	private processRow(row:number):void
	{
		if (row >= this.ds.numRows)
		{
			this.callbackDone();
			return;
		}

		let tr = $('<tr></tr>').appendTo(this.tableRows);
		let th = $('<th></th>').appendTo(tr).css({'text-align': 'left', 'vertical-align': 'top'});
		th.text('Row ' + (row + 1));

		let td = $('<td></td>').appendTo(tr).css({'text-align': 'left', 'vertical-align': 'top', 'border': '1px solid #808080', 'background-color': '#F0F0F0'});
		let flex = $('<div></div>').appendTo(td);
		flex.css({'display': 'flex', 'flex-wrap': 'wrap', 'justify-content': 'flex-start', 'align-items': 'flex-start'});

		let eqr = this.rows[row];
		let nmol = eqr.molList.length;
		for (let n = 0; n < nmol; n++)
		{
			let mol = eqr.molList[n], inchi = eqr.inchiList[n], dhash = eqr.dhashList[n]

			// optionally take this opportunity to make sure that permuted versions are the same
			//this.investigateHashes(mol, dhash);

			let [card, spanMol] = this.generateCard(mol, inchi, dhash, 300);
			card.css({'background-color': 'white', 'border': '1px solid black', 'box-shadow': '3px 3px 5px #808080'});
			flex.append(card);
			spanMol.mouseenter(() => spanMol.css({'background-color': '#C0C0C0', 'border-radius': '5px'}));
			spanMol.mouseleave(() => spanMol.css('background-color', 'transparent'));
			spanMol.click(() => new ZoomDotMol(mol).open());
		}

		// compare InChI's within same row: they are expected to be the same
		for (let i = 0; i < nmol - 1; i++) for (let j = i + 1; j < nmol; j++) 
		{
			let inchi1 = eqr.inchiList[i], inchi2 = eqr.inchiList[j];
			let dhash1 = eqr.dhashList[i], dhash2 = eqr.dhashList[j];
			let badInChI = inchi1 && inchi2 && inchi1 != inchi2, badHash = dhash1 != dhash2
			if (badInChI || badHash)
			{
				let mol1 = this.ds.getMolecule(row, this.colMol[i]), mol2 = this.ds.getMolecule(row, this.colMol[j]);
				let card1 = this.generateCard(mol1, inchi1, dhash1, 200)[0], card2 = this.generateCard(mol2, inchi2, dhash2, 200)[0];
				
				let dualCard = $('<div></div>').appendTo(flex);
				dualCard.css({'display': 'inline-block', 'margin': '0.5em'});
				dualCard.css({'background-color': 'white', 'border': '1px solid black', 'box-shadow': '3px 3px 5px #800000'});

				let divHdr = $('<div></div>').appendTo(dualCard).css({'text-align': 'center', 'color': '#FF0000'});
				if (badInChI && badHash) divHdr.text('InChI & dots both different');
				else if (badInChI) divHdr.text('InChI codes differ');
				else if (badHash) divHdr.text('dot-hashes differ');

				let divMols = $('<div></div>').appendTo(dualCard).css({'text-align': 'center'});
				divMols.append(card1);
				divMols.append(card2);
			}
		}

		// compare InChI's between different rows: they are expected to be different
		for (let n = 0; n < this.rows.length; n++) if (n != row)
			for (let i = 0; i < nmol; i++) for (let j = 0; j < this.rows[n].molList.length; j++)
		{
			let other = this.rows[n];

			let inchi1 = eqr.inchiList[i], inchi2 = other.inchiList[j];
			let dhash1 = eqr.dhashList[i], dhash2 = other.dhashList[j];
			let badInChI = inchi1 && inchi2 && inchi1 == inchi2, badHash = dhash1 == dhash2

			if (badInChI || badHash)
			{
				let mol1 = this.ds.getMolecule(row, this.colMol[i]), mol2 = this.ds.getMolecule(n, this.colMol[j]);
				let card1 = this.generateCard(mol1, inchi1, dhash1, 200)[0], card2 = this.generateCard(mol2, inchi2, dhash2, 200)[0];
				
				let dualCard = $('<div></div>').appendTo(flex);
				dualCard.css({'display': 'inline-block', 'margin': '0.5em'});
				dualCard.css({'background-color': 'white', 'border': '1px solid black', 'box-shadow': '3px 3px 5px #800080'});

				let divHdr = $('<div></div>').appendTo(dualCard).css({'text-align': 'center', 'color': '#800080'});
				if (badInChI && badHash) divHdr.text('InChI & dots falsely equivalent');
				else if (badInChI) divHdr.text('InChI codes falsely equivalent');
				else if (badHash) divHdr.text('Dot-hashes falsely equivalent');

				let divMols = $('<div></div>').appendTo(dualCard).css({'text-align': 'center'});
				divMols.append(card1);
				divMols.append(card2);
			}
		}

		setTimeout(() => this.processRow(row + 1), 1);
	}

	// creates a rectangular DOM block that shows a molecule & its InChI, if available
	private generateCard(mol:Molecule, inchi:string, dhash:string, dimsz:number):[JQuery, JQuery]
	{
		let div = $('<div></div>');
		div.css({'display': 'inline-block', 'margin': '0.5em', 'padding': '0.5em'});

		let divMol = $('<div></div>').appendTo(div).css({'text-align': 'center'});
		let layout = new ArrangeMolecule(mol, this.measure, this.policy, this.effects);
		layout.arrange();
		layout.squeezeInto(0, 0, dimsz, dimsz);
		let gfx = new MetaVector();
		new DrawMolecule(layout, gfx).draw();
		gfx.normalise();
		let spanMol = $('<span></span>').appendTo(divMol);
		spanMol.css('display', 'inline-block');
		$(gfx.createSVG()).appendTo(spanMol);

		let divFormula = $('<div></div>').appendTo(div).css({'text-align': 'center', 'font-size': '70%', 'font-weight': 'bold'});
		divFormula.html(MolUtil.molecularFormula(mol, ['<sub>', '</sub>', '<sup>', '</sup>']));

		if (inchi)
		{
			let maxWidth = Math.max(dimsz, gfx.boundHighX() - gfx.boundLowX());
			let divInChI = $('<div></div>').appendTo(div).css({'text-align': 'left', 'font-size': '70%', 'max-width': maxWidth + 'px', 'word-wrap': 'break-word'});

			let bits = /^(InChI=1S?\/)([\w\d\.]+)(\/.*)$/.exec(inchi);
			if (!bits)
			{
				let span = $('<span></span>').appendTo(divInChI);
				span.css({'color': 'white', 'background-color': '#4E1A09'});
				span.text(inchi);
			}
			else if (!this.sameFormula(bits[2], mol))
			{
				divInChI.append(escapeHTML(bits[1]));
				let span = $('<span></span>').appendTo(divInChI);
				span.css({/*'color': 'white', */'background-color': '#E0E000'});
				span.text(bits[2]);
				divInChI.append(escapeHTML(bits[3]));
			}
			else // all good
			{
				divInChI.text(inchi);
			}
		}
		if (dhash)
		{
			let maxWidth = Math.max(dimsz, gfx.boundHighX() - gfx.boundLowX());
			let divHash = $('<div></div>').appendTo(div);
			divHash.css({'text-align': 'left', 'font-size': '70%', 'max-width': maxWidth + 'px', 'word-wrap': 'break-word', 'margin-top': '0.1em'});
			divHash.css({'border-top': '1px solid #C0C0C0'});
			divHash.text(dhash);
		}

		return [div, spanMol];
	}

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

	// try some variations on dotpath-hash generation
	private investigateHashes(mol:Molecule, dhash:string)
	{
		if (mol.numAtoms <= 1) return;
		for (let count = 20, n = 0; count > 0; count--, n++)
		{
			let a1 = (n % mol.numAtoms) + 1, a2 = ((n + 3) % mol.numAtoms) + 1;
			mol.swapAtoms(a1, a2);
			let phash = new DotHash(new DotPath(mol)).calculate();
			if (dhash != phash) 
			{
				console.log('ORIGINAL:' + dhash);
				console.log('PERMUTED:' + phash);
				throw 'Dot hashes differ';
			}
		}
	}
}

/* EOF */ }