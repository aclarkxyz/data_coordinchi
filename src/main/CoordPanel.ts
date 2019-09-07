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
///<reference path='../data/CallInChI.ts'/>
///<reference path='WindowPanel.ts'/>
///<reference path='AnalyzeResults.ts'/>
///<reference path='EquivalenceResults.ts'/>
///<reference path='CustomStructures.ts'/>

namespace WebMolKit /* BOF */ {

/*
	Coordination panel: content for the main window.
*/

export class CoordPanel extends WindowPanel
{
	private proxyClip = new ClipboardProxy();
	private callInChI:CallInChI;

	private divHeader:JQuery;
	private divSetup:JQuery;
	private divCustom:JQuery;
	private divSummary:JQuery;
	private divResults:JQuery;
	private inputFile:JQuery;
	private chkFailOnly:JQuery;
	private chkInChIFail:JQuery;
	private inputStartAt:JQuery;
	private btnRun:JQuery;
	private btnCancel:JQuery;
	private btnDraw:JQuery;

	// current task
	private filename:string = null;
	private ds:DataSheet = null;
	private task:EquivalenceResults = null;
	private custom:CustomStructures = null;

	// ------------ public methods ------------

	constructor(root:JQuery)
	{
		super(root);

		const {clipboard} = require('electron');
		this.proxyClip.getString = ():string => clipboard.readText();
		this.proxyClip.setString = (str:string):void => clipboard.writeText(str);
		this.proxyClip.canAlwaysGet = ():boolean => true;

		const remote = require('electron').remote;
		this.callInChI = new CallInChI(remote.getGlobal('INCHI_EXEC'));

		document.title = 'Coordination Analysis';

		this.build();
	}

	public selectFile(filename:string):void
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
		//divMain.css({'padding': '0.5em'});

		this.divHeader = $('<div></div>').appendTo(divMain);
		this.divSetup = $('<div></div>').appendTo(divMain);
		this.divCustom = $('<div></div>').appendTo(divMain);
		this.divSummary = $('<div></div>').appendTo(divMain);
		this.divResults = $('<div></div>').appendTo(divMain);

		this.buildHeader();
		this.buildSetup();

// !!
/*		let mol = Molecule.fromString(
			'SketchEl!(6,6)\n' +
			'Pt=0.0000,0.0000;-1,0,i0\n' +
			'Cl=1.5000,0.0000;0,0,i0\n' +
			'Cl=0.0000,1.5000;0,0,i0\n' +
			'Cl=0.0000,-1.5000;0,0,i0\n' +
			'C=-1.5990,0.7500;0,0,i2\n' +
			'C=-1.5990,-0.7500;0,0,i2\n' +
			'1-2=1,0\n' +
			'1-3=1,0\n' +
			'1-4=1,0\n' +
			'1-5=0,0\n' +
			'5-6=2,0\n' +
			'6-1=0,0\n' +
			'!End');	
		let dh = new DotHash(new DotPath(mol));
		console.log('DOTHASH:'+dh.calculate());	
*/
	}

	private buildHeader():void
	{
		let divTitle = $('<div></div>').appendTo(this.divHeader);
		divTitle.css({'margin': '0.25em', 'text-align': 'center'});
		divTitle.css({'font-size': '2.5em', 'font-family': '"Helvetica Neue", Tahoma, Geneva, sans-serif', 'font-weight': 'bold', 'text-shadow': '2px 2px 1px #808080'});
		divTitle.html('<big>C</big>OORDINATION <big>I</big>N<big>C</big>H<big>I</big>');

		let divInfo = $('<div></div>').appendTo(this.divHeader).css({'margin': '0.5em', 'text-align': 'center'});
		let spanInfo = $('<span></span>').appendTo(divInfo).css({'max-width': '40em', 'display': 'inline-block', 'text-align': 'left'});
		spanInfo.append('Validation tools. Runs through training sets or user-specified molecules and examines the performance of standard InChI. ');
		spanInfo.append('Datasets are focused on exotic coordination bonds which cause trouble for most contemporary cheminformatics algorithms. ');
	}

	private buildSetup():void
	{
		this.divSetup.empty();
		let divMain = $('<div></div>').appendTo(this.divSetup).css('padding', '0.5em');

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

		// !! TODO: area for defining user molecules

		// options

		let divOptions = $('<div></div>').appendTo(divMain).css({'text-align': 'center', 'padding': '0.5em'});
		let spanOptions = $('<div></div>').appendTo(divOptions).css({'text-align': 'left', 'display': 'inline-block'});

		let makeCheck = (txt:string):JQuery =>
		{
			let div = $('<div></div>').appendTo(spanOptions);
			let label = $('<label></label>').appendTo(div);
			let chk = $('<input type="checkbox"></input>').appendTo(label);
			label.append(txt);
			return chk;
		};
		let makeInput = (txt:string, width:number):JQuery =>
		{
			let div = $('<div></div>').appendTo(spanOptions);
			div.append(txt);
			let input = $('<input type="text"></input>').appendTo(div);
			input.css({'font': 'inherit', 'margin-left': '0.5em'});
			input.attr('size', width.toString());
			return input;
		};

		this.chkFailOnly = makeCheck('Show failure cases only');
		this.chkInChIFail = makeCheck('Count standard InChI clashes as failures');
		this.inputStartAt = makeInput('Start at row #', 10);

		// action buttons

		let divRun = $('<div></div>').appendTo(divMain);
		divRun.css({'display': 'flex', 'justify-content': 'center'});
		
		this.btnRun = $('<button class="wmk-button wmk-button-primary">Run</button>').appendTo(divRun).css({'margin': '0.5em'});
		this.btnRun.click(() => this.runAnalysis());
		
		this.btnCancel = $('<button class="wmk-button wmk-button-default">Cancel</button>').appendTo(divRun).css({'margin': '0.5em'});
		this.btnCancel.click(() => this.cancelAnalysis());

		this.btnDraw = $('<button class="wmk-button wmk-button-primary">Draw</button>').appendTo(divRun).css({'margin': '0.5em 0.5em 0.5em 1.5em'});
		this.btnDraw.click(() => this.drawStructure());

		this.btnRun.prop('disabled', false);
		this.btnCancel.prop('disabled', true);
		this.btnDraw.prop('disabled', false);

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
			if (filenames.length > 0) this.selectFile(filenames[0]);
		});
	}

	private runAnalysis():void
	{
		this.divSummary.empty();
		this.divResults.empty();

		let areaSummary = $('<div></div>').appendTo(this.divSummary).css('padding', '0.5em');
		areaSummary.css({'border-top': '1px solid #808080', 'border-bottom': '1px solid #808080'});
		areaSummary.text('Loading...');
		let areaResults = $('<div></div>').appendTo(this.divResults).css('padding', '0.5em');

		this.btnRun.prop('disabled', true);
		this.btnCancel.prop('disabled', false);
		this.btnDraw.prop('disabled', true);

		setTimeout(() =>
		{
			this.loadFile();
			areaSummary.empty();
			if (!this.ds) return;
			let opt:EquivalenceResultsOptions =
			{
				'failOnly': this.chkFailOnly.prop('checked'),
				'inchiFail': this.chkInChIFail.prop('checked'),
				'startAt': parseInt(this.inputStartAt.val())
			};
			this.task = new EquivalenceResults(this.ds, this.callInChI, opt, () => this.finishedResults());
	
			this.task.render(areaSummary, areaResults);
			// TODO: fold AnalyseMolecule features into this 
		}, 1);
	}

	// ask the task to cancel, then make re-running an option
	private cancelAnalysis():void
	{
		if (this.task) this.task.cancel();
		this.task = null;
		this.btnRun.prop('disabled', false);
		this.btnCancel.prop('disabled', true);
		this.btnDraw.prop('disabled', false);
	}

	// give the user a chance to draw a molecule
	private drawStructure():void
	{
		if (!this.custom)
		{
			this.custom = new CustomStructures(this.proxyClip);
			this.custom.render(this.divCustom);
		}
		this.custom.sketchNew();
	}

	// obtains the file contents, and sets this.ds if successful
	private loadFile():void
	{
		this.filename = this.inputFile.val();
		this.ds = null;
		if (!this.filename) return;

		const fs = require('fs');
		let strXML = '';
		try {strXML = fs.readFileSync(this.filename).toString();}
		catch (ex) {throw 'Unable to read file: ' + this.filename;}
		this.ds = DataSheetStream.readXML(strXML);		
		if (this.ds == null) throw 'Unable to parse file ' + this.filename;
	}

	// task is done
	private finishedResults():void
	{
		this.btnRun.prop('disabled', false);
		this.btnCancel.prop('disabled', true);

		let paraSave = $('<p></p>').appendTo(this.divSummary);
		paraSave.css({'text-align': 'center'});
		let btnSave = $('<button class="wmk-button wmk-button-primary">Save</button>').appendTo(paraSave);
		btnSave.click(() => this.saveFile());
	}

	// write the current file back to disk
	private saveFile():void
	{
		const electron = require('electron'), fs = require('fs');
		const dialog = electron.remote.dialog; 
		let params:any =
		{
			'title': 'Save Dataset',
			'defaultPath': this.filename,
			'filters':
			[
				{'name': 'DataSheet XML', 'extensions': ['ds']},
				{'name': 'MDL SDfile', 'extensions': ['sdf']}
			]
		};
		dialog.showSaveDialog(params, (filename:string):void =>
		{
			if (!filename) return;
			
			let content = '';
			if (filename.endsWith('.ds')) content = DataSheetStream.writeXML(this.ds);
			else if (filename.endsWith('.sdf')) content = new MDLSDFWriter(this.ds).write();
			else
			{
				alert('Filename has unknown extension: ' + filename);
				return;
			}
			
			try {fs.writeFileSync(filename, content);}
			catch (ex) {alert('Unable to write file: ' + filename);}
		});
	}
}

/* EOF */ }