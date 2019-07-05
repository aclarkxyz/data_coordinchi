/*
	Coordination Complexes for InChI

	(c) 2019 InChI Trust
*/

///<reference path='../../../WebMolKit/src/decl/corrections.d.ts'/>
///<reference path='../../../WebMolKit/src/decl/jquery.d.ts'/>
///<reference path='../../../WebMolKit/src/util/util.ts'/>

namespace WebMolKit /* BOF */ {

/*
	Base class for "main windows": an object that takes up the entire browser window document, responds to resizing, etc.
*/

export class WindowPanel
{
	constructor(public root:JQuery)
	{
		$('body').css('overflow', 'scroll');

		root.css('width', '100%');
		root.css('height', document.documentElement.clientHeight + 'px');
		$(window).resize(() => this.onResize()); 

		root.on('menuAction', (event:any, cmd:string) => this.menuAction(cmd));
	}

	// stub: may be called early on to provide a source file upon which to work
	public loadFile(filename:string):void
	{
	}

	// minimum required functionality for resizing windows; override to capture
	protected onResize()
	{
		this.root.css('height', document.documentElement.clientHeight + 'px');
	}

	// stub: override this to receive menu events
	public menuAction(cmd:string):void
	{
	}
} 

/* EOF */ }