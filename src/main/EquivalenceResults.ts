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

///<reference path='../decl/node.d.ts'/>
///<reference path='../decl/electron.d.ts'/>

///<reference path='../data/AnalyseMolecule.ts'/>

namespace WebMolKit /* BOF */ {

/*
	Equivalence results: reads a datasheet that has rows of structures that are supposed to be the same as each other, but not to any
	of the other rows.
*/

export class EquivalenceResults
{
	private divContent:JQuery;
	private tableRows:JQuery;

	/*private colMol:number;
	private colFormula:number;
	private colError:number;
	private colWarning:number;
	private colFixed:number;
	private tableResults:JQuery = null;*/
	
	private colMol:number[] = [];
	private colInChI:number[] = [];

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

		for (let n = 0; n < this.colMol.length; n++) if (this.ds.notNull(row, this.colMol[n]))
		{
			let mol = this.ds.getMolecule(row, this.colMol[n]);
			let inchi = this.colInChI[n] >= 0 ? this.ds.getString(row, this.colInChI[n]) : null;
			let card = this.generateCard(mol, inchi);
			flex.append(card);
		}

		setTimeout(() => this.processRow(row + 1), 1);
	}

	// creates a rectangular DOM block that shows a molecule & its InChI, if available
	private generateCard(mol:Molecule, inchi:string):JQuery
	{
		let div = $('<div></div>');
		div.css({'display': 'inline-block', 'margin': '0.5em', 'padding': '0.5em'});
		div.css({'background-color': 'white', 'border': '1px solid black', 'box-shadow': '3px 3px 5px #808080'});

		let divMol = $('<div></div>').appendTo(div).css({'text-align': 'center'});
		let layout = new ArrangeMolecule(mol, this.measure, this.policy, this.effects);
		layout.arrange();
		layout.squeezeInto(0, 0, 300, 300);
		let gfx = new MetaVector();
		new DrawMolecule(layout, gfx).draw();
		gfx.normalise();
		$(gfx.createSVG()).appendTo(divMol);

		if (inchi)
		{
			let maxWidth = Math.max(300, gfx.boundHighX() - gfx.boundLowX());
			let divInChI = $('<div></div>').appendTo(div).css({'text-align': 'left', 'max-width': maxWidth + 'px', 'word-wrap': 'break-word'});
			divInChI.text(inchi);
		}

		return div;
	}
}

/* EOF */ }