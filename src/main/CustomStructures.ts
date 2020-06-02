/*
	Coordination Complexes for InChI

	(c) 2019 InChI Trust
*/

///<reference path='../../../WebMolKit/src/decl/corrections.d.ts'/>
///<reference path='../../../WebMolKit/src/decl/jquery/index.d.ts'/>
///<reference path='../../../WebMolKit/src/util/util.ts'/>
///<reference path='../../../WebMolKit/src/util/Geom.ts'/>
///<reference path='../../../WebMolKit/src/ui/Widget.ts'/>
///<reference path='../../../WebMolKit/src/ui/ClipboardProxy.ts'/>
///<reference path='../../../WebMolKit/src/dialog/EditCompound.ts'/>
///<reference path='../../../WebMolKit/src/data/Molecule.ts'/>
///<reference path='../../../WebMolKit/src/data/MoleculeStream.ts'/>
///<reference path='../../../WebMolKit/src/data/DotPath.ts'/>
///<reference path='../../../WebMolKit/src/gfx/Rendering.ts'/>
///<reference path='../../../WebMolKit/src/gfx/ArrangeMolecule.ts'/>
///<reference path='../../../WebMolKit/src/gfx/DrawMolecule.ts'/>

///<reference path='../decl/node.d.ts'/>
///<reference path='../decl/electron.d.ts'/>
///<reference path='../data/DotHash.ts'/>
///<reference path='../data/CallInChI.ts'/>
///<reference path='ZoomDotMol.ts'/>
///<reference path='MoleculeCard.ts'/>

namespace WebMolKit /* BOF */ {

/*
	Manages a list of user-provided structures, for which analysis is shown.
*/

interface StructureInfo
{
	mol:Molecule;
	card:MoleculeCard
}

export class CustomStructures
{
	private divMain:JQuery;
	private structures:StructureInfo[] = [];
	private policy = RenderPolicy.defaultColourOnWhite();

	// ------------ public methods ------------

	constructor(private callInChI:CallInChI, private proxyClip:ClipboardProxy)
	{
		this.policy.data.pointScale = 15;
	}

	public render(parent:JQuery):void
	{
		this.divMain = $('<div></div>').appendTo(parent);
		this.updateOutline();
	}

	// bring up the sketcher dialog for a new compound
	public sketchNew():void
	{
		let dlg = new EditCompound(new Molecule());
		dlg.title = 'Sketch Structure';
		dlg.defineClipboard(this.proxyClip);
		dlg.onSave(() => {this.appendStructure(dlg.getMolecule()); dlg.close();});
		dlg.open();
	}

	// add new structure to the list
	public appendStructure(mol:Molecule):void
	{
		let molExpanded = mol.clone();
		MolUtil.expandAbbrevs(molExpanded, true);

		let inchi:string = null;
		if (this.callInChI.isAvailable) inchi = this.callInChI.calculate(molExpanded).inchi;

		let dhash = new DotHash(new DotPath(molExpanded)).calculate();

		let card = new MoleculeCard(mol, molExpanded, inchi, dhash, 300, this.policy);
		card.generate();

		card.dom.css({'background-color': 'white', 'border': '1px solid black', 'box-shadow': '3px 3px 5px #808080'});
		this.divMain.append(card.dom);
		card.spanMol.mouseenter(() => card.spanMol.css({'background-color': '#C0C0C0', 'border-radius': '5px'}));
		card.spanMol.mouseleave(() => card.spanMol.css('background-color', 'transparent'));
		card.spanMol.click(() => new ZoomDotMol(mol).open());

		card.dom.contextmenu(() => this.contextMenu(card));

		this.structures.push({'mol': mol, 'card': card});
		this.updateOutline();
	}

	// ------------ private methods ------------

	private updateOutline():void
	{
		if (this.structures.length > 0)
		{
			this.divMain.css({'padding': '0.5em', 'border-top': '1px solid #808080'});
		}
		else
		{
			this.divMain.css({'padding': '0', 'border': '0'});
		}
	}

	private contextMenu(card:MoleculeCard):void
	{
		let electron = require('electron');
		let menu = new electron.remote.Menu();

		menu.append(new electron.remote.MenuItem({'label': 'Copy', 'click': () => this.copyNative(card)}));
		menu.append(new electron.remote.MenuItem({'label': 'Copy Molfile V2000', 'click': () => this.copyMDLV2000(card)}));
		//menu.append(new electron.remote.MenuItem({'label': 'Copy Molfile V3000', 'click': () => this.copyMDLV3000()}));
		menu.append(new electron.remote.MenuItem({'label': 'Remove', 'click': () => this.removeStructure(card)}));

		menu.popup({'window': electron.remote.getCurrentWindow()});
	}

	private copyNative(card:MoleculeCard):void
	{
		let str = MoleculeStream.writeNative(card.mol);
		let clipboard = require('electron').clipboard;
		clipboard.writeText(str);
	}

	private copyMDLV2000(card:MoleculeCard):void
	{
		// TODO: pay more attention to how the special fields are written, and at some point S-groups for abbreviations
		let str = MoleculeStream.writeMDLMOL(card.molExpanded);
		let clipboard = require('electron').clipboard;
		clipboard.writeText(str);
	}

	private removeStructure(card:MoleculeCard):void
	{
		for (let n = 0; n < this.structures.length; n++) if (this.structures[n].card === card)
		{
			this.structures[n].card.dom.remove();
			this.structures.splice(n, 1);
			break;
		}
		this.updateOutline();
	}
}

/* EOF */ }