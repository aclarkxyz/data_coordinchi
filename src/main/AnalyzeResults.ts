/*
	Coordination Complexes for InChI

	(c) 2019 InChI Trust
*/

///<reference path='../../../WebMolKit/src/decl/corrections.d.ts'/>
///<reference path='../../../WebMolKit/src/decl/jquery/index.d.ts'/>
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
	Analyze results: takes an input datasheet and determines some validity information.
*/

export class AnalyzeResults
{
	private divContent:JQuery;

	private colMol:number;
	private colFormula:number;
	private colError:number;
	private colWarning:number;
	private colFixed:number;
	private tableResults:JQuery = null;
	
	protected policy:RenderPolicy;
	protected effects:RenderEffects;
	protected measure:ArrangeMeasurement;

	// ------------ public methods ------------

	constructor(private ds:DataSheet, private callbackDone:() => void)
	{
	}

	public render(parent:JQuery):void
	{
		this.divContent = $('<div/>').appendTo(parent);

		this.policy = RenderPolicy.defaultColourOnWhite();
		this.policy.data.pointScale = 15;
		this.effects = new RenderEffects();
		this.measure = new OutlineMeasurement(0, 0, this.policy.data.pointScale);		

		this.startAnalysis();
	}

	// ------------ private methods ------------

	private startAnalysis():void
	{
		// !! TEMPORARY special deal...
		/*let colCorr = this.ds.findColByName('Corrected', DataSheet.COLTYPE_MOLECULE), colBond = this.ds.findColByName('BondAnnot', DataSheet.COLTYPE_STRING);
		if (colCorr >= 0 && colBond >= 0)
		{
			console.log('SPECIAL CORRECTIONS');
			this.ds.ensureColumn('Formula', DataSheet.COLTYPE_STRING, 'Molecular formula implied by structure');
			for (let n = this.ds.numRows - 1; n >= 0; n--) if (this.ds.notNull(n, colCorr))
			{
				this.ds.insertRow(n + 1);
				for (let i = 0; i < this.ds.numCols; i++) this.ds.setObject(n + 1, i, this.ds.getObject(n, i));
				let mol = this.ds.getMolecule(n, colCorr);
				this.ds.setMolecule(n + 1, 'Molecule', mol);
				this.ds.setString(n + 1, colBond, 'corrected');

				let formula = MolUtil.molecularFormula(mol);
				this.ds.setString(n, 'Formula', formula);
				this.ds.setString(n + 1, 'Formula', formula);
			}
			this.ds.deleteColumn(colCorr);
		}*/

		this.colMol = this.ds.firstColOfType(DataSheetColumn.Molecule);
		this.colFormula = this.ds.ensureColumn('Formula', DataSheetColumn.String, 'Molecular formula implied by structure');
		this.colError = this.ds.ensureColumn('Errors', DataSheetColumn.String, 'Fatal flaws with the incoming molecule');
		this.colWarning = this.ds.ensureColumn('Warnings', DataSheetColumn.String, 'Questionable attributes of the incoming molecule');
		this.colFixed = this.ds.ensureColumn('Fixes', DataSheetColumn.String, 'Corrections that could be made unambiguously');

		this.setupHeadings();
		this.processRow(0);
	}

	private setupHeadings():void
	{
		this.tableResults = $('<table/>').appendTo(this.divContent);
		this.tableResults.css({'border': 'none'});

		let tr = $('<tr/>').appendTo(this.tableResults);
		tr.css('background-color', '#C0C0C0');

		for (let n = -1; n < this.ds.numCols; n++)
		{
			let title = n < 0 ? '#' : this.ds.colName(n), ct = n < 0 ? DataSheetColumn.Integer : this.ds.colType(n);
			let th = $('<th/>').appendTo(tr);
			th.css('text-align', ct == DataSheetColumn.String ? 'left' : 'center');
			th.text(title);
		}
	}

	private processRow(row:number):void
	{
		if (row >= this.ds.numRows)
		{
			this.callbackDone();
			return;
		}

		let mol = this.ds.getMolecule(row, this.colMol), formula = this.ds.getString(row, this.colFormula);
		let anal = new AnalyseMolecule(mol, formula);
		try {anal.perform();}
		catch (ex)
		{
			console.log('Failure: ' + ex);
			console.log(ex.stack);
		}

		let error:string[] = [], warning:string[] = [], fixed:string[] = [];
		for (let result of anal.results)
		{
			if (result.type == AnalyseMoleculeType.BadValence)
				error.push('Valence:atom=' + result.atom + '[' + (result.atom == 0 ? '?' : mol.atomElement(result.atom)) + ']:val=' + result.value);
			else if (result.type == AnalyseMoleculeType.OddOxState)
				warning.push('OxState:atom=' + result.atom + '[' + (result.atom == 0 ? '?' : mol.atomElement(result.atom)) + ']:ox=' + result.value);
			else if (result.type == AnalyseMoleculeType.WrongFormula)
				error.push('Formula:' + result.text);
			else if (result.type == AnalyseMoleculeType.NonElement)
				error.push('NonElement:atom=' + result.atom + '[' + result.text + ']');
			else if (result.type == AnalyseMoleculeType.FixCarbonyl)
				fixed.push('FixedCarbonyl:atom=' + result.atom);
			
			else throw '?' + result.type;
		}

		this.ds.setMolecule(row, this.colMol, anal.mol);
		this.ds.setString(row, this.colError, error.join('\n'));
		this.ds.setString(row, this.colWarning, warning.join('\n'));
		this.ds.setString(row, this.colFixed, fixed.join('\n'));

		let tr = $('<tr/>').appendTo(this.tableResults);
		tr.css('background-color', row % 2 ? '#F0F0F0' : '#F8F8F8');

		for (let n = -1; n < this.ds.numCols; n++)
		{
			let title = n < 0 ? '#' : this.ds.colName(n), ct = n < 0 ? DataSheetColumn.Integer : this.ds.colType(n);
			let td = $('<td/>').appendTo(tr);
			let align = ct == DataSheetColumn.String ? 'left' : 'center';
			td.css({'text-align': align, 'vertical-align': 'center', 'white-space': 'pre-wrap'});
			if (n >= 0) td.css({'margin-left': '0.3em', 'margin-right': '0.3em'});

			if (n < 0) td.text((row + 1).toString());
			else if (this.ds.isNull(row, n)) {}
			else if (ct == DataSheetColumn.Molecule)
			{
				if (n == this.colMol)
				{
					let layout = new ArrangeMolecule(mol, this.measure, this.policy, this.effects);
					layout.arrange();
					layout.squeezeInto(0, 0, 300, 300);
					let gfx = new MetaVector();
					new DrawMolecule(layout, gfx).draw();
					gfx.normalise();

					let div = $('<div style="display: inline-block; position: relative;"/>').appendTo(td);
					let domSVG = $(gfx.createSVG()).appendTo(div);
				}
				else td.html('<i>source</i>'); // only want to show the main molecule
			}
			else td.text(this.ds.getObject(row, n).toString());

			// special deal for formula
			if (n == this.colFormula && this.ds.isNull(row, n))
			{
				let divBtn = $('<div/>').appendTo(td);
				divBtn.css({'text-align': 'center'});
				let btnFormula = $('<button class="wmk-button wmk-button-default">Use</button>').appendTo(divBtn);
				btnFormula.click(() => 
				{
					td.text(anal.calcFormula);
					this.ds.setString(row, this.colFormula, anal.calcFormula);
				});
			}
		}

		setTimeout(() => this.processRow(row + 1), 1);
	}
}

/* EOF */ }