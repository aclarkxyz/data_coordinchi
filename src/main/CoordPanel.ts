/*
	Coordination Complexes for InChI

	(c) 2019 InChI Trust
*/

///<reference path='../../../WebMolKit/src/decl/corrections.d.ts'/>
///<reference path='../../../WebMolKit/src/decl/jquery/index.d.ts'/>
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
///<reference path='../data/ExportContent.ts'/>
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
	private modeStereo = false;

	private divHeader:JQuery;
	private divSetup:JQuery;
	private divCustom:JQuery;
	private divSummary:JQuery;
	private divResults:JQuery;
	private inputFile:JQuery;
	private chkStereo:JQuery;
	private chkFailOnly:JQuery;
	private chkInChIFail:JQuery;
	private chkPermute:JQuery;
	private inputStartAt:JQuery;
	private inputEndAt:JQuery;
	private btnRun:JQuery;
	private btnCancel:JQuery;
	private btnDraw:JQuery;
	private btnExport:JQuery;

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
		this.modeStereo = !!remote.getGlobal('PARAM_STEREO');
		this.callInChI = new CallInChI(remote.getGlobal('INCHI_EXEC'), this.modeStereo);


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
		/*else if (cmd == 'copyMDL') this.actionCopyMDL();*/
		else if (cmd == 'paste') this.actionPaste();
		/*else if (cmd == 'delete') new MoleculeActivity(this.sketcher, ActivityType.Delete, {}).execute();
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

		this.divHeader = $('<div/>').appendTo(divMain);
		this.divSetup = $('<div/>').appendTo(divMain);
		this.divCustom = $('<div/>').appendTo(divMain);
		this.divSummary = $('<div/>').appendTo(divMain);
		this.divResults = $('<div/>').appendTo(divMain);

		this.buildHeader();
		this.buildSetup();
	}

	private buildHeader():void
	{
		let divTitle = $('<div/>').appendTo(this.divHeader);
		divTitle.css({'margin': '0.25em', 'text-align': 'center'});
		divTitle.css({'font-size': '2.5em', 'font-family': '"Helvetica Neue", Tahoma, Geneva, sans-serif', 'font-weight': 'bold', 'text-shadow': '2px 2px 1px #808080'});
		divTitle.html('<big>C</big>OORDINATION <big>I</big>N<big>C</big>H<big>I</big>');

		let divInfo = $('<div/>').appendTo(this.divHeader).css({'margin': '0.5em', 'text-align': 'center'});
		let spanInfo = $('<span/>').appendTo(divInfo).css({'max-width': '40em', 'display': 'inline-block', 'text-align': 'left'});
		spanInfo.append('Validation tools. Runs through training sets or user-specified molecules and examines the performance of standard InChI. ');
		spanInfo.append('Datasets are focused on exotic coordination bonds which cause trouble for most contemporary cheminformatics algorithms. ');
	}

	private buildSetup():void
	{
		this.divSetup.empty();
		let divMain = $('<div/>').appendTo(this.divSetup).css('padding', '0.5em');

		let divInput = $('<div/>').appendTo(divMain);
		divInput.css({'width': '100%', 'display': 'flex'});
		let spanTitle = $('<span>File:</span>').appendTo(divInput);
		spanTitle.css({'align-self': 'center'});
		this.inputFile = $('<input type="text" size="40"/>').appendTo(divInput);
		this.inputFile.css({'flex-grow': '1', 'font': 'inherit', 'margin': '0 0.5em 0 0.5em'});
		this.inputFile.keypress((event:JQueryEventObject) => {if (event.keyCode == 13) this.runAnalysis();});
		let btnPick = $('<button class="wmk-button wmk-button-default">Pick</button>').appendTo(divInput);
		btnPick.css({'align-self': 'center'});
		btnPick.click(() => this.pickFilename());

		// options

		let divOptions = $('<div/>').appendTo(divMain).css({'text-align': 'center', 'padding': '0.5em'});
		let spanOptions = $('<div/>').appendTo(divOptions).css({'text-align': 'left', 'display': 'inline-block'});

		let makeCheck = (txt:string, value:boolean):JQuery =>
		{
			let div = $('<div/>').appendTo(spanOptions);
			let label = $('<label/>').appendTo(div);
			let chk = $('<input type="checkbox"/>').appendTo(label);
			label.append(txt);
			chk.prop('checked', value);
			return chk;
		};
		let makeInput = (divParent:JQuery, txt:string, width:number):JQuery =>
		{
			divParent.append(txt);
			let input = $('<input type="text"/>').appendTo(divParent);
			input.css({'font': 'inherit', 'margin': '0 0.5em 0 0.5em'});
			input.attr('size', width.toString());
			return input;
		};

		this.chkStereo = makeCheck('Evaluate stereochemistry', this.modeStereo);
		this.chkFailOnly = makeCheck('Show failure cases only', false);
		this.chkInChIFail = makeCheck('Count standard InChI clashes as failures', false);
		this.chkPermute = makeCheck('Perform internal permutation checks', true);

		let divRange = $('<div/>').appendTo(spanOptions);
		this.inputStartAt = makeInput(divRange, 'Start at row #', 5);
		this.inputEndAt = makeInput(divRange, 'to', 5);

		// action buttons

		let divRun = $('<div/>').appendTo(divMain);
		divRun.css({'display': 'flex', 'justify-content': 'center'});
		
		this.btnRun = $('<button class="wmk-button wmk-button-primary">Run</button>').appendTo(divRun).css({'margin': '0.5em'});
		this.btnRun.click(() => this.runAnalysis());
		
		this.btnCancel = $('<button class="wmk-button wmk-button-default">Cancel</button>').appendTo(divRun).css({'margin': '0.5em'});
		this.btnCancel.click(() => this.cancelAnalysis());

		this.btnDraw = $('<button class="wmk-button wmk-button-primary">Draw</button>').appendTo(divRun).css({'margin': '0.5em 0.5em 0.5em 1.5em'});
		this.btnDraw.click(() => this.drawStructure());

		this.btnExport = $('<button class="wmk-button wmk-button-default">Export</button>').appendTo(divRun).css({'margin': '0.5em 0.5em 0.5em 1.5em'});
		this.btnExport.click(() => this.exportContent());

		this.btnRun.prop('disabled', false);
		this.btnCancel.prop('disabled', true);
		this.btnDraw.prop('disabled', false);
		this.btnExport.prop('disabled', false);

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
		dialog.showOpenDialog(params).then(value =>
		{
			if (value.canceled) return;
			this.selectFile(value.filePaths[0]);
		});
	}

	private runAnalysis():void
	{
		this.divSummary.empty();
		this.divResults.empty();

		let areaSummary = $('<div/>').appendTo(this.divSummary).css('padding', '0.5em');
		areaSummary.css({'border-top': '1px solid #808080', 'border-bottom': '1px solid #808080'});
		areaSummary.text('Loading...');
		let areaResults = $('<div/>').appendTo(this.divResults).css('padding', '0.5em');

		this.btnRun.prop('disabled', true);
		this.btnCancel.prop('disabled', false);
		this.btnDraw.prop('disabled', true);
		this.btnExport.prop('disabled', true);

		setTimeout(() =>
		{
			this.loadFile();
			areaSummary.empty();
			if (!this.ds) return;
			let opt:EquivalenceResultsOptions =
			{
				'stereochemistry': this.chkStereo.prop('checked'),
				'failOnly': this.chkFailOnly.prop('checked'),
				'inchiFail': this.chkInChIFail.prop('checked'),
				'permute': this.chkPermute.prop('checked'),
				'startAt': parseInt(this.inputStartAt.val().toString()),
				'endAt': parseInt(this.inputEndAt.val().toString()),
			};
			this.task = new EquivalenceResults(this.ds, this.callInChI, opt, () => this.finishedResults());
			this.task.render(areaSummary, areaResults);
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
		this.btnExport.prop('disabled', false);
	}

	// give the user a chance to draw a molecule
	private drawStructure():void
	{
		if (!this.custom)
		{
			this.custom = new CustomStructures(this.callInChI, this.proxyClip);
			this.custom.render(this.divCustom);
		}
		this.custom.sketchNew();
	}

	// write the file to SDfile for use in other environments
	private exportContent():void
	{
		this.loadFile();
		if (!this.ds) return;

		this.btnRun.prop('disabled', true);
		this.btnCancel.prop('disabled', false);
		this.btnDraw.prop('disabled', true);
		this.btnExport.prop('disabled', true);

		let sdFN = this.filename;
		let dot = sdFN.lastIndexOf('.');
		if (dot >= 0) sdFN = sdFN.substring(0, dot);
		sdFN += '.sdf';

		let stereochemistry = this.chkStereo.prop('checked');
		new ExportContent(this.ds, sdFN, stereochemistry).perform();

		this.btnRun.prop('disabled', false);
		this.btnCancel.prop('disabled', true);
		this.btnDraw.prop('disabled', false);
		this.btnExport.prop('disabled', false);
	}

	// obtains the file contents, and sets this.ds if successful
	private loadFile():void
	{
		this.filename = this.inputFile.val().toString();
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
		this.btnDraw.prop('disabled', false);
		this.btnExport.prop('disabled', false);

		let paraSave = $('<p/>').appendTo(this.divSummary);
		paraSave.css({'text-align': 'center'});
		
		let btnSaveData = $('<button class="wmk-button wmk-button-primary">Save Data</button>').appendTo(paraSave);
		btnSaveData.click(() => this.saveDataFile());
		paraSave.append(' ');
		let btnSaveView = $('<button class="wmk-button wmk-button-default">Save HTML</button>').appendTo(paraSave);
		btnSaveView.click(() => this.saveHTMLFile());
	}

	// write the current file back to disk
	private saveDataFile():void
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
		dialog.showSaveDialog(params).then(value =>
		{
			if (value.canceled) return;
			
			let content = '';
			if (value.filePath.endsWith('.ds')) content = DataSheetStream.writeXML(this.ds);
			else if (value.filePath.endsWith('.sdf')) content = new MDLSDFWriter(this.ds).write();
			else
			{
				alert('Filename has unknown extension: ' + value.filePath);
				return;
			}
			
			try {fs.writeFileSync(value.filePath, content);}
			catch (ex) {alert('Unable to write file: ' + value.filePath);}
		});
	}

	// save an HTML representation of the results
	private saveHTMLFile():void
	{
		const STYLE = 
			'<style type="text/css">\n' +
			'body\n' +
			'{\n' +
			'  color: #000000;\n' +
			'  font-family: sans-serif;\n' +
			'}\n' +
			'</style>';
		let head = '<title>Coordination InChI Validation</title>\n' + STYLE + '\n';
		let body = this.divResults[0].outerHTML;
		let html = '<html>\n<head>' + head + '</head>\n<body>\n' + body + '\n</body>\n</html>';

		const electron = require('electron'), fs = require('fs');
		const dialog = electron.remote.dialog; 
		let params:any =
		{
			'title': 'Save Results as HTML',
			'filters':
			[
				{'name': 'HTML', 'extensions': ['html']},
			]
		};
		dialog.showSaveDialog(params).then(value =>
		{
			if (value.canceled) return;
			
			try {fs.writeFileSync(value.filePath, html);}
			catch (ex) {alert('Unable to write file: ' + value.filePath);}
		});
	}

	private actionPaste():void
	{
		this.proxyClip.triggerPaste();
	}
}

/* EOF */ }