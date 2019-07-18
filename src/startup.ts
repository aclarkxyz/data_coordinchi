/*
	Coordination Complexes for InChI

	(c) 2019 InChI Trust
*/

///<reference path='../../WebMolKit/src/decl/corrections.d.ts'/>
///<reference path='../../WebMolKit/src/decl/jquery.d.ts'/>
///<reference path='../../WebMolKit/src/util/util.ts'/>

///<reference path='decl/node.d.ts'/>
///<reference path='main/CoordPanel.ts'/>

namespace WebMolKit /* BOF */ {

$ = require('./jquery.js');

/*
	Startup: gets the ball rolling, and provide some high level window handling.
*/

let BASE_APP = ''; // base URL location for the app's program files (could be URL or filename)

export function runCoordInChI(rootID:string):void
{
	let root = $('#' + rootID);

	// node/electron imports; note these are defined inside the function so as not to perturb normal web-access, which does not
	// include these libraries
	const path = require('path');
	const electron = require('electron');
	const process = require('process');

	BASE_APP = path.normalize('file:/' + __dirname);

	let url = window.location.href.substring(0, window.location.href.lastIndexOf('/'));
	initWebMolKit(path.normalize(url + '/res'));

	// unpack web params: if present, they determine where to go from here
 	let params = window.location.search.substring(1).split('&');
	let panelClass:string = null; // default is straight to molecule editing
	let filename:string = null;
	for (let p of params)
	{
		let eq = p.indexOf('=');
		if (eq < 0) continue;
		let key = p.substring(0, eq), val = decodeURIComponent(p.substring(eq + 1));
		if (key == 'panel') panelClass = 'WebMolKit.' + val;
		else if (key == 'fn') filename = val;
	}	

	if (!panelClass)
	{
		let dw = new CoordPanel(root);
		if (filename) dw.selectFile(filename);
	}
	else
	{
		let constructor = eval(panelClass);
		let dw:WindowPanel = new constructor(root);
		if (filename) dw.selectFile(filename);
	}
}

// high level functionality for opening a window, with a given panel as content
export function openNewWindow(panelClass:string, filename?:string):void
{
	const electron = require('electron');
	const WEBPREF = {'nodeIntegration': true};
	let bw = new electron.remote.BrowserWindow({'width': 800, 'height': 700, 'icon': 'app/img/icon.png', 'webPreferences': WEBPREF});
	let url = BASE_APP + '/index.html?panel=' + panelClass;
	if (filename) url += '&fn=' + encodeURIComponent(filename);
	bw.loadURL(url);
	/*bw.on('closed', function() {bw = null;});*/
}

/* EOF */ }