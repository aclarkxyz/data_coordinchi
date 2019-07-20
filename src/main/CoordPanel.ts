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
///<reference path='WindowPanel.ts'/>
///<reference path='AnalyzeResults.ts'/>
///<reference path='EquivalenceResults.ts'/>

namespace WebMolKit /* BOF */ {

/*
	Coordination panel: content for the main window.
*/

export class CoordPanel extends WindowPanel
{
	private proxyClip = new ClipboardProxy();

	private inputFile:JQuery;
	private btnAnalyse:JQuery;
	private btnEquivalence:JQuery;
	private divResults:JQuery;
	private divOutcome:JQuery;

	// current task
	private filename:string = null;
	private ds:DataSheet = null;

	// ------------ public methods ------------

	constructor(root:JQuery)
	{
		super(root);

		const {clipboard} = require('electron');
		this.proxyClip.getString = ():string => clipboard.readText();
		this.proxyClip.setString = (str:string):void => clipboard.writeText(str);
		this.proxyClip.canAlwaysGet = ():boolean => true;

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

		let divRun = $('<div></div>').appendTo(divMain);
		divRun.css({'display': 'flex', 'justify-content': 'center'});
		
		this.btnAnalyse = $('<button class="wmk-button wmk-button-primary">Analyse</button>').appendTo(divRun).css({'margin': '0.5em'});
		//addTooltip(this.btnAnalyse, 'Show the datasheet and look for structure validity issues.');
		this.btnAnalyse.click(() => this.runAnalysis());
		
		this.btnEquivalence = $('<button class="wmk-button wmk-button-primary">Equivalence</button>').appendTo(divRun).css({'margin': '0.5em'});
		//addTooltip(this.btnEquivalence, 'Display equivalences between related structures.');
		this.btnEquivalence.click(() => this.runEquivalence());

		this.divResults = $('<div></div>').appendTo(divMain);
		this.divOutcome = $('<div></div>').appendTo(divMain);

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
		this.loadFile();
		if (!this.ds) return;
		this.preProcess();
		new AnalyzeResults(this.ds, () => this.finishedResults()).render(this.divResults);
	}

	private runEquivalence():void
	{
		this.loadFile();
		if (!this.ds) return;
		this.preProcess();
		new EquivalenceResults(this.ds, () => this.finishedResults()).render(this.divResults);
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

	// get things ready for a task
	private preProcess():void
	{
		this.divResults.empty();
		this.divOutcome.empty();

		this.btnAnalyse.prop('disabled', true);
		this.btnEquivalence.prop('disabled', true);
	}

	// task is done
	private finishedResults():void
	{
		this.btnAnalyse.prop('disabled', false);
		this.btnEquivalence.prop('disabled', false);

		let paraSave = $('<p></p>').appendTo(this.divOutcome);
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