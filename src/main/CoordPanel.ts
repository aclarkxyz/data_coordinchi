/*
	Coordination Complexes for InChI

	(c) 2019 InChI Trust
*/

///<reference path='../../../WebMolKit/src/decl/corrections.d.ts'/>
///<reference path='../../../WebMolKit/src/decl/jquery.d.ts'/>
///<reference path='../../../WebMolKit/src/util/util.ts'/>
///<reference path='../../../WebMolKit/src/sketcher/Sketcher.ts'/>
///<reference path='../../../WebMolKit/src/ui/ClipboardProxy.ts'/>
///<reference path='../../../WebMolKit/src/data/Molecule.ts'/>
///<reference path='../../../WebMolKit/src/data/MoleculeStream.ts'/>
///<reference path='../../../WebMolKit/src/data/MDLWriter.ts'/>

///<reference path='../decl/node.d.ts'/>
///<reference path='../decl/electron.d.ts'/>

///<reference path='../data/AnalyseMolecule.ts'/>
///<reference path='MainPanel.ts'/>

namespace WebMolKit /* BOF */ {

/*
	Drawing window: dedicated entirely to the sketching of a molecular structure.
*/

export class CoordPanel extends MainPanel
{
	private proxyClip = new ClipboardProxy();

	private inputFile:JQuery;
	private btnAnalyse:JQuery;
	private divResults:JQuery;

	// current task
	private filename:string = null;
	private ds:DataSheet = null;
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

	constructor(root:JQuery)
	{
		super(root);

		const {clipboard} = require('electron');
		this.proxyClip.getString = ():string => clipboard.readText();
		this.proxyClip.setString = (str:string):void => clipboard.writeText(str);
		this.proxyClip.canAlwaysGet = ():boolean => true;

		document.title = 'Coordination Analysis';

		this.policy = RenderPolicy.defaultColourOnWhite();
		this.policy.data.pointScale = 15;
		this.effects = new RenderEffects();
		this.measure = new OutlineMeasurement(0, 0, this.policy.data.pointScale);		

		this.build();
	}

	public loadFile(filename:string):void
	{
		//const process = require('process');
		//console.log('CWD:'+process.cwd());
		const path = require('path');
		if (!path.isAbsolute(filename))
		{
			//console.log('RESOLVE:'+path.resolve());	
			filename = path.relative(path.resolve() + path.sep, filename);
		}

		this.inputFile.val(filename);
	}

	public menuAction(cmd:string):void
	{
		/*if (cmd == 'new') openNewWindow('DrawPanel');
		else if (cmd == 'open') this.actionFileOpen();
		else if (cmd == 'save') this.actionFileSave();
		else if (cmd == 'saveAs') this.actionFileSaveAs();
		else if (cmd == 'exportSVG') this.actionFileExportSVG();
		else if (cmd == 'undo') this.sketcher.performUndo();
		else if (cmd == 'redo') this.sketcher.performRedo();
		else if (cmd == 'cut') this.actionCopy(true);
		else*/ if (cmd == 'copy') document.execCommand('copy');
		/*else if (cmd == 'copyMDL') this.actionCopyMDL();
		else if (cmd == 'paste') this.actionPaste();
		else if (cmd == 'delete') new MoleculeActivity(this.sketcher, ActivityType.Delete, {}).execute();
		else if (cmd == 'selectAll') new MoleculeActivity(this.sketcher, ActivityType.SelectAll, {}).execute();
		else if (cmd == 'zoomFull') this.sketcher.autoScale();
		else if (cmd == 'zoomIn') this.sketcher.zoom(1.25);
		else if (cmd == 'zoomOut') this.sketcher.zoom(0.8);
		else console.log('MENU:' + cmd);*/
	}

	// ------------ private methods ------------

	private build():void
	{
		let divMain = $('<div></div>').appendTo(this.root);
		divMain.css({'padding': '0.5em'});

		let divInput = $('<div></div>').appendTo(divMain);
		divInput.css({'width': '100%', 'display': 'flex'});
		let spanTitle = $('<span>File:</span>').appendTo(divInput);
		spanTitle.css({'align-self': 'center'});
		this.inputFile = $('<input type="text" size="40"></input>').appendTo(divInput);
		this.inputFile.css({'flex-grow': '1', 'font': 'inherit', 'margin': '0 0.5em 0 0.5em'});
		this.inputFile.keypress((event:JQueryEventObject) => {if (event.keyCode == 13) this.runAnalysis();});
		let btnPick = $('<button class="wmk-button wmk-button-default">Pick</button>').appendTo(divInput);
		btnPick.css({'align-self': 'center'});
		btnPick.click(() => this.pickFilename());

		let paraRun = $('<p></p>').appendTo(divMain);
		paraRun.css({'text-align': 'center'});
		this.btnAnalyse = $('<button class="wmk-button wmk-button-primary">Analyse</button>').appendTo(paraRun);
		this.btnAnalyse.click(() => this.runAnalysis());

		this.divResults = $('<div></div>').appendTo(divMain);

		this.inputFile.focus();
	}

	private pickFilename():void
	{
		const electron = require('electron');
		const dialog = electron.remote.dialog; 
		let params:any =
		{
			'title': 'Open DataSheet',
			'filters':
			[
				{'name': 'Molecular DataSheet', 'extensions': ['ds']}
			]
		};
		dialog.showOpenDialog(params, (filenames:string[]):void =>
		{
			if (filenames.length > 0) this.loadFile(filenames[0]);
		});
	}

	private runAnalysis():void
	{
		this.filename = null;
		this.ds = null;
		this.divResults.empty();

		this.filename = this.inputFile.val();
		if (!this.filename) return;

		const fs = require('fs');
		let strXML = '';
		try {strXML = fs.readFileSync(this.filename).toString();}
		catch (ex) {throw 'Unable to read file: ' + this.filename;}
		this.ds = DataSheetStream.readXML(strXML);		
		if (this.ds == null) throw 'Unable to parse file ' + this.filename;

		this.startAnalysis();
	}

	// ------------ private methods ------------

	private startAnalysis():void
	{
		this.btnAnalyse.prop('disabled', true);

		this.colMol = this.ds.firstColOfType(DataSheet.COLTYPE_MOLECULE);
		this.colFormula = this.ds.ensureColumn('Formula', DataSheet.COLTYPE_STRING, 'Molecular formula implied by structure');
		this.colError = this.ds.ensureColumn('Errors', DataSheet.COLTYPE_STRING, 'Fatal flaws with the incoming molecule');
		this.colWarning = this.ds.ensureColumn('Warnings', DataSheet.COLTYPE_STRING, 'Questionable attributes of the incoming molecule');
		this.colFixed = this.ds.ensureColumn('Fixes', DataSheet.COLTYPE_STRING, 'Corrections that could be made unambiguously');

		this.setupHeadings();
		this.processRow(0);
	}

	private setupHeadings():void
	{
		this.tableResults = $('<table></table>').appendTo(this.divResults);
		this.tableResults.css({'border': 'none'});

		let tr = $('<tr></tr>').appendTo(this.tableResults);
		tr.css('background-color', '#C0C0C0');

		for (let n = -1; n < this.ds.numCols; n++)
		{
			let title = n < 0 ? '#' : this.ds.colName(n), ct = n < 0 ? DataSheet.COLTYPE_INTEGER : this.ds.colType(n);
			let th = $('<th></th>').appendTo(tr);
			th.css('text-align', ct == DataSheet.COLTYPE_STRING ? 'left' : 'center');
			th.text(title);
		}
	}

	private processRow(row:number):void
	{
		if (row >= this.ds.numRows)
		{
			this.finishAnalysis();
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
				error.push('Valence:atom=' + result.atom + '[' + (result.atom == 0 ? '?' : mol.atomElement(result.atom)) + ']:' + result.value);
			else if (result.type == AnalyseMoleculeType.OddOxState)
				warning.push('OxState:atom=' + result.atom + '[' + (result.atom == 0 ? '?' : mol.atomElement(result.atom)) + ']:' + result.value);
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

		let tr = $('<tr></tr>').appendTo(this.tableResults);
		tr.css('background-color', row % 2 ? '#F0F0F0' : '#F8F8F8');

		for (let n = -1; n < this.ds.numCols; n++)
		{
			let title = n < 0 ? '#' : this.ds.colName(n), ct = n < 0 ? DataSheet.COLTYPE_INTEGER : this.ds.colType(n);
			let td = $('<td></td>').appendTo(tr);
			let align = ct == DataSheet.COLTYPE_STRING ? 'left' : 'center';
			td.css({'text-align': align, 'vertical-align': 'center', 'white-space': 'pre-wrap'});
			if (n >= 0) td.css({'margin-left': '0.3em', 'margin-right': '0.3em'});

			if (n < 0) td.text((row + 1).toString());
			else if (this.ds.isNull(row, n)) {}
			else if (ct == DataSheet.COLTYPE_MOLECULE)
			{
				let layout = new ArrangeMolecule(mol, this.measure, this.policy, this.effects);
				layout.arrange();
				layout.squeezeInto(0, 0, 300, 300);
				let gfx = new MetaVector();
				new DrawMolecule(layout, gfx).draw();
				gfx.normalise();

				let div = $('<div style="display: inline-block; position: relative;"></div>').appendTo(td);
				let domSVG = $(gfx.createSVG()).appendTo(div);
			}
			else td.text(this.ds.getObject(row, n).toString());

			// special deal for formula
			if (n == this.colFormula && this.ds.isNull(row, n))
			{
				let divBtn = $('<div></div>').appendTo(td);
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

	private finishAnalysis():void
	{
		this.btnAnalyse.prop('disabled', false);

		let paraSave = $('<p></p>').appendTo(this.divResults);
		paraSave.css({'text-align': 'center'});
		let btnSave = $('<button class="wmk-button wmk-button-primary">Save</button>').appendTo(paraSave);
		btnSave.click(() => this.saveFile());
	}

	// write the current file back to disk
	private saveFile():void
	{
		let strXML = DataSheetStream.writeXML(this.ds);

		const fs = require('fs');
		try {fs.writeFileSync(this.filename, strXML);}
		catch (ex) {throw 'Unable to write file: ' + this.filename;}
	}
}

/* EOF */ }